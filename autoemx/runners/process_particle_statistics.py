#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: recalculate the particle size statistics of a sample, filtering out some particles or frames.

The particles are read from the sample ledger (``ledger.json``) when it holds the particles of a particle size
statistics measurement; the recalculated statistics then replace ``SampleLedger.particle_stats``, recording the
particles and frames left out. Since all particles stay in the ledger, recalculating without exclusions restores
the statistics of the measurement. Samples measured before the ledger stored the particles are read from their
particle sizes file (``<sample_ID>_Par_sizes.csv``).

The statistics are also saved with the suffix ``_processed`` (particle sizes, summary statistics and histogram).

Created on Mon Oct 13 2026
@author: Andrea
"""

import logging
import os
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import ParticleInfo, SampleLedger
from autoemx.core.particle_statistics import recalculate_particle_stats, save_particle_statistics
from autoemx.utils import get_sample_dir


def _ledger_with_size_stats_particles(sample_dir: str) -> Optional[SampleLedger]:
    """Ledger of the sample, if it holds particles of a particle size statistics measurement or frames searched by EDS."""
    ledger_path = os.path.join(sample_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
    if not os.path.exists(ledger_path):
        return None
    ledger = load_sample_ledger(ledger_path)
    if not any(p.selection == "size_stats" for p in ledger.particles) and not ledger.eds_size_scans:
        return None
    return ledger


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
        IDs of the particles to remove ('Particle #' column).
    frame_ids_to_filter : iterable of str, optional
        IDs of the frames whose particles are removed ('Frame ID' column). Only exact matches: the subframes of a
        frame (e.g. 'A0_a1' in 'A0') must be listed separately.
    include_eds_auto_particles : bool, optional
        Whether to include the particles measured by EDS with automatic particle selection (ledger only). Default:
        the setting of the measurement.
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
        result = recalculate_particle_stats(
            ledger, particle_ids_to_filter, frame_ids_to_filter, include_eds_auto_particles
        )
        if result is None:
            logging.warning("No particles left to build statistics: do not filter out all particles and frames.")
            return None
        missing_particles = [pid for pid in particle_ids_to_filter if pid not in result.excluded_particle_ids]
        missing_frames = [fid for fid in frame_ids_to_filter if fid not in result.excluded_frame_ids]
        kept_ids = set(result.particle_ids)
        particles = [p for p in ledger.particles if p.id in kept_ids]
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

    if ledger is not None:
        ledger.particle_stats = result
        ledger.to_json_file(os.path.join(sample_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT))

    return save_particle_statistics(particles, sample_dir, sample_ID, output_file_suffix=output_file_suffix,
                                    verbose=verbose, bands=result.bands if ledger is not None else None,
                                    excluded_frame_ids=result.excluded_frame_ids if ledger is not None else (),
                                    image_format=(
                                        result.powder_meas_cfg.size_stats_plot_format if ledger is not None else "png"
                                    ),
                                    confidence_intervals=result.confidence_intervals_95 if ledger is not None else None,
                                    eds_scans=result.eds_scans if ledger is not None else ())
