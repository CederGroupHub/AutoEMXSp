#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: reanalyse a particle size statistics measurement from its saved frame images, with other segmentation
settings (brightness threshold, segmentation model, minimum number of pixels per particle).

The frames of each size band (and the frames searched by EDS on the sample) are segmented again, as during the
measurement, from the images saved in ``SEM images``; their particles replace those of the measurement in the ledger,
and the statistics are rebuilt with the same size bands, frames and scanned area. The statistics are also saved with
the suffix ``_reanalysed``.

Frames recorded without image (frames without particles of measurements made before AutoEMX 1.6.3) are kept as
frames without particles; frames whose image is missing are left out. Particles measured by EDS stay in the ledger.

Created on Fri Oct 16 2026
@author: Andrea
"""

import logging
import os
from types import SimpleNamespace
from typing import List, Optional

import numpy as np
import pandas as pd

import autoemx.utils.constants as cnst
from autoemx import microscope_drivers as EM_driver
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import ParticleSizeBand, SampleLedger
from autoemx.config.runtime_configs import MeasurementConfig, PowderMeasurementConfig
from autoemx.core.composition_analysis.ledger import LedgerModule
from autoemx.core.em_runtime import image_utilities
from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder, min_valid_diameter_um
from autoemx.core.particle_statistics import recalculate_particle_stats, save_particle_statistics
from autoemx.utils import get_sample_dir

__all__ = ["reanalyze_particle_statistics"]


class _SavedFrames:
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
        pass  # the frame images are those saved by the measurement


def _load_frame(images_dir: str, frame) -> Optional[np.ndarray]:
    """Greyscale image saved for a frame, or None if it has none."""
    if not frame.image_file:
        return None
    path = os.path.join(images_dir, frame.image_file)
    if not os.path.exists(path):
        return None
    return image_utilities.load_raw_saved_image(path)[..., 0]  # saved as RGB with identical channels


def _finder(frames: _SavedFrames, powder_meas_cfg: PowderMeasurementConfig, images_dir: str) -> EM_Particle_Finder:
    return EM_Particle_Finder(frames, powder_meas_cfg, results_dir=images_dir, verbose=False, development_mode=False)


def _reanalyse_bands(ledger: SampleLedger, frames: _SavedFrames, cfg: PowderMeasurementConfig, images_dir: str):
    """Segment again the frames of the size bands of the measurement. Returns the new bands and frame counts."""
    finder = _finder(frames, cfg, images_dir)
    finder.attach_stats_ledger(ledger)
    bands, n_empty_kept, missing = [], 0, []
    for band in ledger.particle_stats.bands:
        # The band starts at its smallest particles with enough pixels (raised if min_par_pixels was raised)
        pixel_sizes = [f.pixel_size_um for f in band.frames] or [0.0]
        d_min_um = max(band.d_min_um, min_valid_diameter_um(max(pixel_sizes), cfg.min_par_pixels))
        if d_min_um > band.d_min_um * (1 + 1e-9):
            logging.warning(
                f"⚠️ Size band {band.d_min_um:.4g}–{band.d_max_um:.4g} µm now starts at {d_min_um:.4g} µm "
                f"(min_par_pixels {cfg.min_par_pixels}): smaller particles of this band are not measured."
            )
        new_band = ParticleSizeBand(d_min_um=min(d_min_um, band.d_max_um), d_max_um=band.d_max_um,
                                    frame_width_um=band.frame_width_um, n_frames=0, scanned_area_um2=0.0)
        finder._set_size_band_area_limits(new_band.d_min_um, new_band.d_max_um, band.frame_width_um or 0.0)
        finder._stats_band = new_band
        for frame in band.frames:
            image = _load_frame(images_dir, frame)
            if image is not None:
                frames.show_frame(frame.frame_id, frame.center_mm, image)
                finder._move_and_get_particles_stats_in_frame(frame_image=image, pixel_size=frame.pixel_size_um)
                new_band.frames[-1] = new_band.frames[-1].model_copy(update={"image_file": frame.image_file})
            elif frame.n_particles == 0:
                new_band.add_frame(frame.model_copy(update={"particle_marks_px": []}), [])  # kept: no particle
                n_empty_kept += 1
            else:
                missing.append(frame.frame_id)
        finder._stats_band = None
        bands.append(new_band)
    return bands, n_empty_kept, missing


def _reanalyse_eds_frames(ledger: SampleLedger, frames: _SavedFrames, overrides: dict, images_dir: str):
    """Segment again the frames searched by EDS. Returns the new EDS scans and the frames kept as recorded."""
    eds_cfg = (ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig()).model_copy(update=overrides)
    finder = _finder(frames, eds_cfg, images_dir)
    rebuilt = SampleLedger.model_construct(eds_size_scans=[])
    analyzer = SimpleNamespace(EM_controller=SimpleNamespace(particle_finder=finder))
    kept: List[str] = []
    for scan in ledger.eds_size_scans:
        eds_cfg.min_area_par = float(np.pi * scan.d_min_um**2 / 4)
        eds_cfg.max_area_par = float(np.pi * scan.d_max_um**2 / 4)
        for frame in scan.frames:
            image = _load_frame(images_dir, frame)
            if image is None:
                finder.pending_eds_frames.append(((scan.d_min_um, scan.d_max_um, scan.frame_width_um or 0.0), frame))
                kept.append(frame.frame_id)
                continue
            frames.show_frame(frame.frame_id.removeprefix("EDS-"), frame.center_mm, image)
            finder._get_particles_coordinates_in_frame(frame_image=image, pixel_size=frame.pixel_size_um,
                                                       results_dir=images_dir)
            window, new_frame = finder.pending_eds_frames[-1]
            finder.pending_eds_frames[-1] = (window, new_frame.model_copy(update={
                "frame_id": frame.frame_id, "image_file": frame.image_file,
            }))
    LedgerModule._store_eds_search_frames(analyzer, rebuilt)
    return rebuilt.eds_size_scans, kept


def reanalyze_particle_statistics(
    sample_ID: str,
    results_path: str,
    par_brightness_thresh: Optional[int] = None,
    par_segmentation_model: Optional[str] = None,
    min_par_pixels: Optional[int] = None,
    include_eds_frames: bool = True,
    output_file_suffix: str = "_reanalysed",
    verbose: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Reanalyse the particle size statistics measurement of a sample from its saved frame images.

    Parameters
    ----------
    sample_ID : str
        Sample ID (name of the sample folder).
    results_path : str
        Folder containing the sample folder (also searched in its subfolders).
    par_brightness_thresh, par_segmentation_model, min_par_pixels : optional
        Segmentation settings to use instead of those of the measurement (see ``PowderMeasurementConfig``).
    include_eds_frames : bool, optional
        Whether to also reanalyse the frames searched by EDS on the sample (default: True).
    output_file_suffix : str, optional
        Suffix of the files saved (default: '_reanalysed').
    verbose : bool, optional
        If True (default), prints the statistics.

    Returns
    -------
    pandas.DataFrame or None
        Summary statistics, or None if no particle is left.

    Raises
    ------
    FileNotFoundError
        If the sample folder or its ledger is not found.
    ValueError
        If the sample has no completed particle size statistics measurement.
    """
    sample_dir = get_sample_dir(results_path, sample_ID)
    ledger_path = os.path.join(sample_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
    if not os.path.exists(ledger_path):
        raise FileNotFoundError(f"No ledger found at {ledger_path}")
    ledger = load_sample_ledger(ledger_path)
    measurement = ledger.particle_stats
    if measurement is None or not measurement.bands:
        raise ValueError(
            f"Sample '{sample_ID}' has no completed particle size statistics measurement"
            + (" (run it again to complete it)." if ledger.particle_stats_progress is not None else ".")
        )
    images_dir = os.path.join(sample_dir, cnst.IMAGES_DIR)

    overrides = {key: value for key, value in (
        ("par_brightness_thresh", par_brightness_thresh),
        ("par_segmentation_model", par_segmentation_model),
        ("min_par_pixels", min_par_pixels),
    ) if value is not None}
    cfg = PowderMeasurementConfig.model_validate({**measurement.powder_meas_cfg.model_dump(), **overrides})

    # Particles of the measurement are replaced (those measured by EDS stay)
    measured_ids = {pid for band in measurement.bands for pid in band.particle_ids}
    eds_ids = {s.acquisition_details.particle_id for s in ledger.spectra
               if s.acquisition_details is not None and s.acquisition_details.particle_id is not None}
    ledger.particles = [p for p in ledger.particles if p.id not in measured_ids or p.id in eds_ids]

    frames = _SavedFrames(sample_ID, ledger.configs.microscope_cfg.ID)
    bands, n_empty_kept, missing = _reanalyse_bands(ledger, frames, cfg, images_dir)
    eds_kept: List[str] = []
    if include_eds_frames and ledger.eds_size_scans:
        ledger.eds_size_scans, eds_kept = _reanalyse_eds_frames(ledger, frames, overrides, images_dir)

    if n_empty_kept:
        logging.warning(f"⚠️ {n_empty_kept} frame(s) without saved image kept as frames without particles.")
    if missing:
        logging.warning(f"⚠️ Frames left out, their image is missing: {missing}")
    if eds_kept:
        logging.warning(f"⚠️ {len(eds_kept)} frame(s) searched by EDS without saved image kept as recorded.")
    if measurement.excluded_particle_ids:
        logging.warning("⚠️ Particles left out of the measurement are not left out of the reanalysis (new particles).")

    ledger.particle_stats = measurement.model_copy(update={"bands": bands, "powder_meas_cfg": cfg})
    result = recalculate_particle_stats(
        ledger, frame_ids_to_filter=measurement.excluded_frame_ids,
        include_eds_auto_particles=measurement.powder_meas_cfg.include_eds_auto_particles and include_eds_frames,
    )
    if result is None:
        logging.warning("No particles found in the reanalysed frames.")
        ledger.to_json_file(ledger_path)
        return None
    ledger.particle_stats = result
    ledger.to_json_file(ledger_path)
    if verbose:
        logging.info(f"✅ Reanalysed {sum(b.n_frames for b in bands)} frame(s): {result.size_stats.n_par_analysed} particles.")

    used_ids = set(result.particle_ids)
    return save_particle_statistics(
        [p for p in ledger.particles if p.id in used_ids], sample_dir, sample_ID,
        output_file_suffix=output_file_suffix, verbose=verbose, bands=result.bands,
        excluded_frame_ids=result.excluded_frame_ids, image_format=cfg.size_stats_plot_format,
        confidence_intervals=result.confidence_intervals_95, eds_scans=result.eds_scans,
    )
