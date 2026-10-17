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
) -> Tuple[Dict[str, FrameSegmentation], List[str]]:
    """
    Segment frames again from their saved images with the segmentation settings of ``powder_meas_cfg``, each among
    the diameters searched in it: the result of segmentation ``segmentation_id`` of each frame (by frame id), and the
    frames without saved image.
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
        results[frame.frame_id] = finder.segment_frame(image, frame.pixel_size_um, frame.search_window_um,
                                                       segmentation_id)
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
    no result is replaced).

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
    if not to_segment:
        return segmentation_id, []
    results, without_image = segment_frames(
        to_segment, SavedFrames(ledger.sample_id, ledger.configs.microscope_cfg.ID), powder_meas_cfg,
        os.path.join(sample_dir, cnst.IMAGES_DIR), segmentation_id,
    )
    for frame_id, result in results.items():
        ledger.add_frame_segmentation(frame_id, result)
    return segmentation_id, without_image
