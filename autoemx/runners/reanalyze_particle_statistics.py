#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: reanalyse a particle size statistics measurement from its saved frame images, with other segmentation
settings (brightness threshold, segmentation model and its options, minimum number of pixels per particle).

The frames of each size band (and the frames searched by EDS on the sample) are segmented again, as during the
measurement, from the images saved in ``SEM images``: the results are added to the frames under a segmentation id
(those of other segmentations are kept, see ``SampleLedger.segmentations``), and the statistics rebuilt from them with
the same size bands, frames and scanned area are added as a new result, made the active one. The statistics are also
saved with the suffix ``_reanalysed``.

Frames recorded without image (frames without particles of measurements made before AutoEMX 1.6.3) are kept as
frames without particles; frames whose image is missing are left out. Particles measured by EDS stay in the ledger.

Created on Fri Oct 16 2026
@author: Andrea
"""

import logging
import os
from typing import Any, Dict, Optional

import pandas as pd

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import FrameSegmentation
from autoemx.config.runtime_configs import PowderMeasurementConfig
from autoemx.core.particle_statistics import recalculate_particle_stats, save_particle_statistics, size_stats_inputs
from autoemx.core.saved_frames import resegment_frames
from autoemx.utils import get_sample_dir

__all__ = ["reanalyze_particle_statistics"]


def reanalyze_particle_statistics(
    sample_ID: str,
    results_path: str,
    par_brightness_thresh: Optional[int] = None,
    par_segmentation_model: Optional[str] = None,
    min_par_pixels: Optional[int] = None,
    par_segmentation_params: Optional[Dict[str, Any]] = None,
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
    par_segmentation_params : dict, optional
        Options of the segmentation model (e.g. ``{"checkpoint": "facebook/sam2.1-hiera-small"}`` for SAM2). Default:
        those of the measurement if the model is the same, else none (the model defaults).
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
    active = ledger.size_result
    measurement = ledger.measurement(active.measurement_id) if active is not None else None
    if measurement is None:
        measurement = next((m for m in reversed(ledger.size_measurements) if m.status == "completed"), None)
    if measurement is None and not ledger.frames_of("eds_search"):
        raise ValueError(
            f"Sample '{sample_ID}' has no completed particle size statistics measurement, nor frames searched by EDS"
            + (" (run the measurement again to complete it)." if ledger.size_measurements else ".")
        )

    overrides = {key: value for key, value in (
        ("par_brightness_thresh", par_brightness_thresh),
        ("par_segmentation_model", par_segmentation_model),
        ("min_par_pixels", min_par_pixels),
    ) if value is not None}
    base_cfg = (measurement.powder_meas_cfg if measurement is not None
                else ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig())
    if par_segmentation_params is not None:
        overrides["par_segmentation_params"] = dict(par_segmentation_params)
    elif par_segmentation_model is not None and par_segmentation_model != base_cfg.par_segmentation_model:
        overrides["par_segmentation_params"] = {}  # options of another model
    cfg = PowderMeasurementConfig.model_validate({**base_cfg.model_dump(), **overrides})
    # Without measurement (sample measured only by EDS), the frames searched by EDS
    include_eds = measurement is None or (base_cfg.include_eds_auto_particles and include_eds_frames)
    frames = ledger.frames_of("size_stats", measurement.measurement_id) if measurement is not None else []
    if include_eds:
        frames += ledger.frames_of("eds_search")

    # Results added to the frames (those of other segmentations are kept)
    segmentation_id, without_image = resegment_frames(ledger, sample_dir, frames, cfg)
    n_empty_kept, missing = 0, []
    for frame_id in without_image:
        frame = ledger.frame(frame_id)
        if all(s.n_particles == 0 for s in frame.segmentations):
            # Recorded without image because it had no particle (measurements before AutoEMX 1.6.3)
            ledger.add_frame_segmentation(frame_id, FrameSegmentation(segmentation_id=segmentation_id,
                                                                      d_range_um=frame.search_window_um))
            n_empty_kept += 1
        else:
            missing.append(frame_id)
    if n_empty_kept:
        logging.warning(f"⚠️ {n_empty_kept} frame(s) without saved image kept as frames without particles.")
    if missing:
        logging.warning(f"⚠️ Frames left out, their image is missing: {missing}")
    if active is not None and active.excluded_particles:
        logging.warning("⚠️ Particles left out of the statistics are not left out of the reanalysis (new particles).")

    eds_segmentation_id = segmentation_id if include_eds and ledger.frames_of("eds_search") else None
    measurement_id = measurement.measurement_id if measurement is not None else None
    result = recalculate_particle_stats(
        ledger, frame_ids_to_filter=active.excluded_frame_ids if active is not None else None,
        include_eds_auto_particles=include_eds, measurement_id=measurement_id,
        segmentation_id=segmentation_id, eds_segmentation_id=eds_segmentation_id,
    )
    if result is None:
        logging.warning("No particles found in the reanalysed frames.")
        ledger.to_json_file(ledger_path)
        return None
    result = ledger.append_size_result(result)
    ledger.to_json_file(ledger_path)
    if verbose:
        n_frames = sum(b.n_frames for b in result.bands)
        logging.info(f"✅ Reanalysed {n_frames} frame(s) (segmentation {segmentation_id}): "
                     f"{result.size_stats.n_par_analysed} particles.")

    inputs = size_stats_inputs(ledger, measurement_id, segmentation_id, eds_segmentation_id, include_eds)
    return save_particle_statistics(
        inputs.particles, sample_dir, sample_ID, output_file_suffix=output_file_suffix, verbose=verbose,
        bands=inputs.bands, excluded_frame_ids=result.excluded_frame_ids, image_format=cfg.size_stats_plot_format,
        confidence_intervals=result.confidence_intervals_95, eds_scans=inputs.eds_scans,
    )
