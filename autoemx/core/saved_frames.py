#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle detection in saved frame images, with other segmentation settings than those of the acquisition.

Used to reanalyse particle size statistics (``autoemx.runners.reanalyze_particle_statistics``), and to segment the
frames searched by EDS again with the settings of the particle size statistics measured after them
(``EMXSp_Composition_Analyzer.collect_particle_stats_after_eds``), so that all the frames of the statistics are
segmented alike. The results are added to the frames under a segmentation id: none is replaced.
"""

import logging
import os
from types import SimpleNamespace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import autoemx.utils.constants as cnst
from autoemx import microscope_drivers as EM_driver
from autoemx.config.ledger_schemas import FrameRecord, FrameSegmentation, SampleLedger, SegmentationConfig
from autoemx.config.runtime_configs import MeasurementConfig, PowderMeasurementConfig
from autoemx.core.em_runtime import image_utilities
from autoemx.core.particle_outlines import load_outlines, save_outlines

__all__ = ["SavedFrames", "load_frame_image", "resegment_frames", "segment_frames"]


class SavedFrames:
    """Stand-in for the microscope controller: the particle finder analyses saved frames, at their position."""

    def __init__(self, sample_id: str, microscope_id: str):
        self.sample_id = sample_id
        self.im_width = self.im_height = 0
        self.pixel_size_um = 1.0
        self.measurement_cfg = MeasurementConfig(type=MeasurementConfig.PARTICLE_STATS_MEAS_TYPE_KEY)
        self.frame_navigator = SimpleNamespace(current_frame_label=None, _current_pos=None, current_pos_is_measured=False)
        try:
            EM_driver.load_microscope_driver(microscope_id)
            self._transform = np.asarray(EM_driver.image_to_stage_coords_transform)
        except Exception:
            self._transform = np.array([1, -1])
            logging.warning(f"Microscope driver of '{microscope_id}' not available: particle positions may be mirrored.")

    @property
    def current_frame_label(self):
        return self.frame_navigator.current_frame_label

    def show_frame(self, label: Optional[str], center_mm, image: np.ndarray) -> None:
        """Set the frame being analysed: its label, position and size."""
        self.frame_navigator.current_frame_label = label
        self.frame_navigator._current_pos = (center_mm.x, center_mm.y) if center_mm is not None else None
        self.frame_navigator.current_pos_is_measured = center_mm is not None
        self.im_height, self.im_width = image.shape[:2]

    def convert_pixel_pos_to_mm(self, pos_pixels):
        return image_utilities.convert_pixel_pos_to_mm(
            pos_pixels, self.im_width, self.im_height, self.pixel_size_um,
            self.frame_navigator._current_pos or (0.0, 0.0), self._transform,
        )

    def save_frame_image(self, *args, **kwargs):
        pass  # the frame images are those saved during the acquisition


def load_frame_image(images_dir: str, frame) -> Optional[np.ndarray]:
    """Greyscale image saved for a frame, or None if it has none."""
    if not frame.image_file:
        return None
    path = os.path.join(images_dir, frame.image_file)
    if not os.path.exists(path):
        return None
    return image_utilities.load_raw_saved_image(path)[..., 0]  # saved as RGB with identical channels


def segment_frames(
    frames: Sequence[FrameRecord],
    saved_frames: SavedFrames,
    powder_meas_cfg: PowderMeasurementConfig,
    images_dir: str,
    segmentation_id: int,
    sample_dir: Optional[str] = None,
) -> Tuple[Dict[str, FrameSegmentation], List[str]]:
    """
    Segment frames again from their saved images with the segmentation settings of ``powder_meas_cfg``, each among
    the diameters searched in it: the result of segmentation ``segmentation_id`` of each frame (by frame id), and the
    frames without saved image. With ``sample_dir``, the outlines of the particles of each frame are saved (see
    particle_outlines).
    """
    from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder

    finder = EM_Particle_Finder(saved_frames, powder_meas_cfg.model_copy(deep=True), results_dir=images_dir,
                                verbose=False, development_mode=False)
    results: Dict[str, FrameSegmentation] = {}
    without_image: List[str] = []
    for frame in frames:
        image = load_frame_image(images_dir, frame)
        if image is None:
            without_image.append(frame.frame_id)
            continue
        saved_frames.show_frame(frame.label, frame.center_mm, image)
        results[frame.frame_id], outlines = finder.segment_frame(image, frame.pixel_size_um, frame.search_window_um,
                                                                 segmentation_id, with_outlines=True)
        if sample_dir is not None:
            save_outlines(sample_dir, segmentation_id, frame.frame_id, outlines)
        logging.info(f"Frame {frame.frame_id}: {results[frame.frame_id].n_particles} particles "
                     f"({powder_meas_cfg.par_segmentation_model}).")
    return results, without_image


def resegment_frames(
    ledger: SampleLedger,
    sample_dir: str,
    frames: Sequence[FrameRecord],
    powder_meas_cfg: PowderMeasurementConfig,
    source: str = "reanalysis",
) -> Tuple[int, List[str]]:
    """
    Segment frames of a sample again, from their saved images, with the segmentation settings of
    ``powder_meas_cfg``: the results are added to the frames (those already segmented with these settings are kept;
    no result is replaced), and the outlines of their particles saved (see particle_outlines). Frames already
    segmented with these settings but without saved outlines (segmented before AutoEMX saved them) are segmented again
    for their outlines, kept only if the same particles are found.

    Returns
    -------
    segmentation_id : int
        Id of the segmentation with these settings (added to the ledger if new).
    without_image : list of str
        Frames without saved image, not segmented.
    """
    segmentation_id = ledger.append_segmentation(
        SegmentationConfig.from_powder_meas_cfg(powder_meas_cfg, 0, source=source))
    to_segment = [f for f in frames if f.segmentation(segmentation_id) is None]
    without_outlines = [f for f in frames if f.segmentation(segmentation_id) is not None and f.image_file
                        and load_outlines(sample_dir, segmentation_id, f.frame_id) is None]
    if not to_segment and not without_outlines:
        return segmentation_id, []
    saved_frames = SavedFrames(ledger.sample_id, ledger.configs.microscope_cfg.ID)
    images_dir = os.path.join(sample_dir, cnst.IMAGES_DIR)
    results, without_image = segment_frames(to_segment, saved_frames, powder_meas_cfg, images_dir, segmentation_id,
                                            sample_dir)
    for frame_id, result in results.items():
        ledger.add_frame_segmentation(frame_id, result)
    if without_outlines:
        _add_missing_outlines(ledger, sample_dir, without_outlines, saved_frames, powder_meas_cfg, images_dir,
                              segmentation_id)
    return segmentation_id, without_image


def _add_missing_outlines(ledger, sample_dir, frames, saved_frames, powder_meas_cfg, images_dir, segmentation_id):
    """Outlines of frames segmented before they were saved: segmented again, kept if the same particles are found
    (so that the particles left out, identified by their number in the frame, stay the same)."""
    from autoemx.core.particle_outlines import outlines_path

    results, _ = segment_frames(frames, saved_frames, powder_meas_cfg, images_dir, segmentation_id, sample_dir)
    mismatched = []
    for frame_id, result in results.items():
        stored = ledger.frame(frame_id).segmentation(segmentation_id)
        same = (len(result.particle_areas_um2) == len(stored.particle_areas_um2)
                and np.allclose(result.particle_areas_um2, stored.particle_areas_um2, rtol=1e-6))
        if not same:
            os.remove(outlines_path(sample_dir, segmentation_id, frame_id))
            mismatched.append(frame_id)
    if mismatched:
        logging.warning(f"⚠️ Particle outlines not added to frames {mismatched}: segmented again, they give other "
                        "particles than recorded (shown as circles).")
