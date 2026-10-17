#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: recalculate the particle size statistics of a sample, filtering out some particles or frames.

The particles are read from the frames of the sample ledger (``ledger.json``), under the segmentations of its
active statistics result; the recalculated statistics are added as a new result (made the active one), recording the
particles (as (frame_id, index)) and frames left out. Frames and earlier results stay in the ledger, so recalculating
without exclusions restores the statistics of the measurement. Samples measured before the ledger stored particle
sizes are read from their particle sizes file (``<sample_ID>_Par_sizes.csv``).

The statistics are also saved with the suffix ``_processed`` (particle sizes, summary statistics and histogram).

Created on Mon Oct 13 2026
@author: Andrea
"""

import logging
import os
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import ParticleInfo, SampleLedger
from autoemx.core.particle_statistics import recalculate_particle_stats, save_particle_statistics, size_stats_inputs
from autoemx.utils import get_sample_dir


def _ledger_with_size_stats_particles(sample_dir: str) -> Optional[SampleLedger]:
    """Ledger of the sample, if it holds frames with particle sizes (of a measurement, or searched by EDS)."""
    ledger_path = os.path.join(sample_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
    if not os.path.exists(ledger_path):
        return None
    ledger = load_sample_ledger(ledger_path)
    return ledger if ledger.frames else None


def _particles_from_csv(sample_dir: str, sample_ID: str) -> List[ParticleInfo]:
    """Particles of the particle sizes file of the sample."""
    par_data_path = os.path.join(sample_dir, f"{sample_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv")
    if not os.path.exists(par_data_path):
        raise FileNotFoundError(f"Particle size file not found at {par_data_path}")
    par_data = pd.read_csv(par_data_path, dtype={cnst.FRAME_ID_DF_KEY: str}, keep_default_na=False)
    # Particles counted in frames searched by EDS have no ID: the file is only fully read back through the ledger
    par_data = par_data[par_data[cnst.PAR_ID_DF_KEY].astype(str) != ""]
    return [
        ParticleInfo(
            id=int(row[cnst.PAR_ID_DF_KEY]),
            frame_id=str(row[cnst.FRAME_ID_DF_KEY]) or None,
            area_um=float(row[cnst.PAR_AREA_UM_KEY]),
            eq_diameter_um=(
                float(row[cnst.PAR_EQ_D_KEY])
                if cnst.PAR_EQ_D_KEY in par_data.columns and not pd.isna(row[cnst.PAR_EQ_D_KEY])
                else float(2.0 * np.sqrt(float(row[cnst.PAR_AREA_UM_KEY]) / np.pi))
            ),
        )
        for _, row in par_data.iterrows()
    ]


def process_particle_statistics(
    sample_ID: str,
    results_path: str,
    particle_ids_to_filter: Optional[Iterable[int]] = None,
    frame_ids_to_filter: Optional[Iterable[str]] = None,
    include_eds_auto_particles: Optional[bool] = None,
    particles_to_filter: Optional[Iterable[Tuple[str, int]]] = None,
    result_id: Optional[int] = None,
    output_file_suffix: str = "_processed",
    verbose: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Recalculate the particle size statistics of a sample without some particles or frames.

    Parameters
    ----------
    sample_ID : str
        Sample ID (name of the sample folder).
    results_path : str
        Folder containing the sample folder (also searched in its subfolders).
    particle_ids_to_filter : iterable of int, optional
        IDs of the particles to remove ('Particle #' column), for samples read from their particle sizes file.
    frame_ids_to_filter : iterable of str, optional
        IDs of the frames whose particles are removed ('Frame ID' column). Only exact matches: the subframes of a
        frame (e.g. 'A0_a1' in 'A0') must be listed separately.
    include_eds_auto_particles : bool, optional
        Whether to include what EDS measured (ledger only). Default: the setting of the active result.
    particles_to_filter : iterable of (frame_id, index), optional
        Particles to remove, by frame and index in the frame ('Frame ID' and 'Particle #' columns), ledger only.
    result_id : int, optional
        Statistics result of the ledger to recalculate (its measurement and segmentations). Default: the active one.
    output_file_suffix : str, optional
        Suffix of the files saved (default: '_processed').
    verbose : bool, optional
        If True (default), prints the statistics.

    Returns
    -------
    pandas.DataFrame or None
        Summary statistics, or None if no particle is left.

    Raises
    ------
    FileNotFoundError
        If the sample folder is not found, or the sample has neither particle size statistics in its ledger nor a
        particle sizes file.
    """
    sample_dir = get_sample_dir(results_path, sample_ID)
    if not os.path.exists(sample_dir):
        raise FileNotFoundError(f"Could not find sample at {sample_dir}. Please check 'sample_ID' and 'results_path'.")

    particle_ids_to_filter = [int(pid) for pid in (particle_ids_to_filter or [])]
    frame_ids_to_filter = [str(frame_id) for frame_id in (frame_ids_to_filter or [])]

    ledger = _ledger_with_size_stats_particles(sample_dir)
    if ledger is not None:
        particles_to_filter = [(str(f), int(i)) for f, i in (particles_to_filter or [])]
        if particle_ids_to_filter:
            logging.warning("⚠️ Particles of the ledger frames are identified by (frame ID, index): use "
                            "particles_to_filter (particle_ids_to_filter ignored).")
        inputs = {}
        if result_id is not None:
            base = next((r for r in ledger.size_results if r.result_id == result_id), None)
            if base is None:
                raise ValueError(f"No statistics result {result_id} in the ledger of sample '{sample_ID}'")
            inputs = {"measurement_id": base.measurement_id, "segmentation_id": base.segmentation_id,
                      "eds_segmentation_id": base.eds_segmentation_id}
            if include_eds_auto_particles is None:
                include_eds_auto_particles = base.powder_meas_cfg.include_eds_auto_particles
            if not include_eds_auto_particles:
                inputs["eds_segmentation_id"] = None
        result = recalculate_particle_stats(
            ledger, particles_to_filter, frame_ids_to_filter, include_eds_auto_particles, **inputs
        )
        if result is None:
            logging.warning("No particles left to build statistics: do not filter out all particles and frames.")
            return None
        excluded = set(map(tuple, result.excluded_particles))
        missing_particles = [key for key in particles_to_filter if key not in excluded]
        missing_frames = [fid for fid in frame_ids_to_filter if fid not in result.excluded_frame_ids]
        inputs = size_stats_inputs(ledger, result.measurement_id, result.segmentation_id, result.eds_segmentation_id,
                                   result.powder_meas_cfg.include_eds_auto_particles)
        particles = inputs.particles
    else:
        logging.info("No particle size statistics in the ledger: reading the particle sizes file.")
        particles = _particles_from_csv(sample_dir, sample_ID)
        missing_particles = [pid for pid in particle_ids_to_filter if pid not in {p.id for p in particles}]
        missing_frames = [fid for fid in frame_ids_to_filter if fid not in {p.frame_id for p in particles}]
        particles = [
            p for p in particles if p.id not in particle_ids_to_filter and p.frame_id not in frame_ids_to_filter
        ]
        if not particles:
            logging.warning("No particles left to build statistics: do not filter out all particles and frames.")
            return None

    if missing_particles or missing_frames:
        logging.warning(
            "⚠️ Not found in the particle sizes (not removed): "
            + ", ".join(
                part for part in (
                    f"particle IDs {missing_particles}" if missing_particles else "",
                    f"frame IDs {missing_frames}" if missing_frames else "",
                ) if part
            )
        )
    else:
        logging.info("✅ All requested particle IDs and frame IDs were found and removed.")

    if ledger is None:
        return save_particle_statistics(particles, sample_dir, sample_ID, output_file_suffix=output_file_suffix,
                                        verbose=verbose)
    ledger.append_size_result(result)
    ledger.to_json_file(os.path.join(sample_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT))
    return save_particle_statistics(particles, sample_dir, sample_ID, output_file_suffix=output_file_suffix,
                                    verbose=verbose, bands=inputs.bands, excluded_frame_ids=result.excluded_frame_ids,
                                    image_format=result.powder_meas_cfg.size_stats_plot_format,
                                    confidence_intervals=result.confidence_intervals_95, eds_scans=inputs.eds_scans,
                                    excluded_particles=result.excluded_particles)
