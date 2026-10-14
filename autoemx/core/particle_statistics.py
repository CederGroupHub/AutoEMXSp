#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle size statistics: summary statistics of the equivalent diameters of particles, and their export.

The statistics use the particles found by a particle size statistics measurement (selection 'size_stats') and,
if ``PowderMeasurementConfig.include_eds_auto_particles`` (default), the particles of the same sample measured by
EDS with automatic particle selection (selection 'auto'), before or after it. The result is stored in the sample
ledger (``SampleLedger.particle_stats``); a new measurement replaces it, and particles of earlier measurements stay
in the ledger but are left out of the new statistics.

Independent of the microscope: used both at the end of a particle statistics measurement
(``EM_Particle_Finder.get_particle_stats``) and to recalculate statistics from saved particle sizes
(``autoemx.runners.process_particle_statistics``).

Files saved in the sample folder (``<suffix>`` is optional):

- ``<sample_ID>_Par_sizes<suffix>.csv``: ID, frame, area and equivalent diameter of each particle;
- ``<sample_ID>_Par_size_stats<suffix>.csv``: summary statistics of the equivalent diameters;
- ``<sample_ID>_Par_size_distribution_hist<suffix>.<format>``: size distribution plot (histogram and cumulative
  distributions of the equivalent diameters).
"""

import os
from collections import Counter
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import autoemx.utils.constants as cnst
from autoemx._logging import get_logger
from autoemx.config.ledger_schemas import (
    ParticleInfo,
    ParticleSizeBand,
    ParticleSizeStats,
    ParticleStatsResult,
    SampleLedger,
)
from autoemx.config.runtime_configs import PowderMeasurementConfig
from autoemx.utils.helper import print_double_separator

logger = get_logger(__name__)

__all__ = [
    "band_particle_weights",
    "bootstrap_confidence_intervals",
    "edge_correction_factor",
    "build_particle_stats_result",
    "compute_particle_size_stats",
    "particle_sizes_dataframe",
    "recalculate_particle_stats",
    "save_particle_size_histogram",
    "save_particle_statistics",
    "select_size_stats_particles",
]

#: Description of the measure, first column of the statistics file.
SIZE_STATS_MEASURE = "Equivalent particle diameter in μm"

#: Statistics given a confidence interval (name: key in :func:`compute_particle_size_stats`).
CI_QUANTITIES = {"D10": "D10", "D50": "median", "D90": "D90", "D43": "D43", "Dv10": "Dv10", "Dv50": "Dv50", "Dv90": "Dv90"}

#: Number of bootstrap resamples of the confidence intervals.
N_BOOTSTRAP_RESAMPLES = 200

#: Minimum number of resampled units (frames, or particles without frame records) in each size band for confidence
#: intervals: with fewer, the bootstrap underestimates the variability.
MIN_BOOTSTRAP_UNITS = 5


def _weighted_percentile(sorted_values: np.ndarray, sorted_weights: np.ndarray, percent: float) -> float:
    """
    Percentile of weighted values (sorted, positive weights), by linear interpolation of the weighted cumulative
    distribution, each value at the middle of its share of the total weight (below the first and above the last
    middle: the smallest and largest value); with equal weights, ``numpy.percentile`` with method 'hazen'.
    """
    cumulative = np.cumsum(sorted_weights)
    positions = (cumulative - sorted_weights / 2) / cumulative[-1]
    return float(np.interp(percent / 100.0, positions, sorted_values))


def compute_particle_size_stats(
    eq_diameters_um: Sequence[float],
    weights: Optional[Sequence[float]] = None,
) -> Dict[str, Union[int, float]]:
    """
    Summary statistics of particle equivalent diameters, optionally weighted.

    Parameters
    ----------
    eq_diameters_um : sequence of float
        Equivalent diameters of the particles, in µm.
    weights : sequence of float, optional
        Non-negative (number) weight of each particle (see :func:`band_particle_weights`); particles of weight 0
        are left out. Default: all particles count equally.

    Returns
    -------
    dict
        Number-weighted: ``n_par_analysed`` (particles counted), ``mean``, ``stdev`` (population standard
        deviation), ``median``, ``max``, ``min`` and the percentiles ``D10``, ``D25``, ``D75``, ``D90``.
        Volume-weighted (each particle also weighted by the volume of a sphere of its diameter, as measured by
        laser diffraction): ``D43`` (mean, D[4,3]) and the percentiles ``Dv10``, ``Dv50``, ``Dv90``. In µm.

    Raises
    ------
    ValueError
        If no diameter (or no particle of positive weight) is given.
    """
    # Number-weighted percentiles: numpy's default with equal weights (as before AutoEMX 1.6.0), else
    # _weighted_percentile; volume-weighted ones always _weighted_percentile, the largest particles weighing most
    d = np.asarray(eq_diameters_um, dtype=float)
    if weights is None:
        w = np.ones(d.size)
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != d.shape or np.any(w < 0) or not np.all(np.isfinite(w)):
            raise ValueError("weights must be finite, non-negative and one per diameter")
        d, w = d[w > 0], w[w > 0]
    if d.size == 0:
        raise ValueError("No particle diameters to compute statistics from")

    w = w / w.max()  # relative weights (band weights are small numbers, particles per µm²)
    if np.allclose(w, 1.0, rtol=0, atol=1e-12):  # all particles count equally: as numpy
        stats = {
            "n_par_analysed": int(d.size),
            "mean": float(np.mean(d)),
            "stdev": float(np.std(d)),
            "median": float(np.median(d)),
            "max": float(np.max(d)),
            "min": float(np.min(d)),
            "D10": float(np.percentile(d, 10)),
            "D25": float(np.percentile(d, 25)),
            "D75": float(np.percentile(d, 75)),
            "D90": float(np.percentile(d, 90)),
        }
        w = np.ones(d.size)
    else:
        order = np.argsort(d, kind="stable")
        d, w = d[order], w[order]
        mean = float(np.sum(w * d) / np.sum(w))
        stats = {
            "n_par_analysed": int(d.size),
            "mean": mean,
            "stdev": float(np.sqrt(np.sum(w * (d - mean) ** 2) / np.sum(w))),
            "median": _weighted_percentile(d, w, 50),
            "max": float(d[-1]),
            "min": float(d[0]),
            "D10": _weighted_percentile(d, w, 10),
            "D25": _weighted_percentile(d, w, 25),
            "D75": _weighted_percentile(d, w, 75),
            "D90": _weighted_percentile(d, w, 90),
        }

    order = np.argsort(d, kind="stable")
    d_sorted, volume_weights = d[order], (w * d**3)[order]
    stats.update({
        "D43": float(np.sum(volume_weights * d_sorted) / np.sum(volume_weights)),
        "Dv10": _weighted_percentile(d_sorted, volume_weights, 10),
        "Dv50": _weighted_percentile(d_sorted, volume_weights, 50),
        "Dv90": _weighted_percentile(d_sorted, volume_weights, 90),
    })
    return stats


def edge_correction_factor(
    bbox_width_px: int,
    bbox_height_px: int,
    frame_width_px: int,
    frame_height_px: int,
    edge_margin_px: int,
) -> float:
    """
    Miles-Lantuéjoul correction of a particle measured in a frame where particles touching the edge are left out.

    A particle is kept only if its bounding box stays more than ``edge_margin_px`` pixels from every edge of the
    frame, which is less likely the larger it is. Its weight is multiplied by the inverse of the fraction of frame
    positions where it is kept, ``W·H / ((W - w - 2m - 1)(H - h - 2m - 1))``, so that large particles are not
    under-counted.

    Raises
    ------
    ValueError
        If the particle cannot fit in the frame without touching its edge.
    """
    n_x = frame_width_px - bbox_width_px - 2 * edge_margin_px - 1
    n_y = frame_height_px - bbox_height_px - 2 * edge_margin_px - 1
    if n_x <= 0 or n_y <= 0:
        raise ValueError("The particle cannot fit in the frame without touching its edge")
    return float(frame_width_px * frame_height_px / (n_x * n_y))


def _band_members(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand],
) -> Tuple[List[List[int]], List[List[int]]]:
    """
    Indices of the particles of each band, and of those among them found in it by the measurement.

    The particles found in a band by the measurement belong to it; the others (e.g. measured by EDS) belong to the
    band of their diameter (the larger one at a boundary), or to none.
    """
    measurement_band = {pid: i for i, band in enumerate(bands) for pid in band.particle_ids}
    eq_diameters_um = particle_sizes_dataframe(particles)[cnst.PAR_EQ_D_KEY].to_numpy() if particles else []
    members: List[List[int]] = [[] for _ in bands]
    found: List[List[int]] = [[] for _ in bands]
    for k, (particle, d) in enumerate(zip(particles, eq_diameters_um)):
        band_idx = measurement_band.get(particle.id)
        if band_idx is not None:
            found[band_idx].append(k)
        else:
            band_idx = next(
                (i for i, band in enumerate(bands) if band.d_min_um * (1 - 1e-9) <= d <= band.d_max_um * (1 + 1e-9)),
                None,
            )
        if band_idx is not None:
            members[band_idx].append(k)
    return members, found


def band_particle_weights(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand],
    excluded_frame_ids: Iterable[str] = (),
) -> Tuple[np.ndarray, List[ParticleSizeBand]]:
    """
    Weights of particles pooled from several size bands, so that each band counts in proportion to its number of
    particles per unit area scanned (not to the number of particles collected in it).

    The particles found in a band by the measurement belong to it; particles measured by EDS belong to the band
    of their diameter. Each band's total weight is the number of particles the measurement found in it (those
    left out excepted, each counted ``edge_correction`` times: see :func:`edge_correction_factor`) per µm²
    scanned (frames left out excepted, if recorded). It is shared among its particles in proportion to their
    edge correction (1 for EDS particles): EDS particles refine the size distribution within a band without
    changing its weight. Particles in no band, or in a band where the measurement found no particle, get
    weight 0.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles included in the statistics.
    bands : sequence of ParticleSizeBand
        Size bands of the measurement.
    excluded_frame_ids : iterable of str, optional
        Frames left out of the statistics: their area is not counted (if the band records its frames).

    Returns
    -------
    weights : numpy.ndarray
        Weight of each particle (all 1 without bands).
    bands : list of ParticleSizeBand
        The bands with their ``number_density_per_mm2`` and ``n_eds_particles``.
    """
    if not bands:
        return np.ones(len(particles)), []
    members, found_in_band = _band_members(particles, bands)
    corrections = np.array([p.edge_correction or 1.0 for p in particles])
    weights = np.zeros(len(particles))
    weighted_bands = []
    for band, band_members, found in zip(bands, members, found_in_band):
        n_found = len(found)
        area_um2 = band.scanned_area_um2 - band.frame_area_um2(excluded_frame_ids)
        density_per_um2 = float(np.sum(corrections[found])) / area_um2 if area_um2 > 0 else 0.0
        for k in band_members:
            weights[k] = density_per_um2 * corrections[k] / np.sum(corrections[band_members])
        weighted_bands.append(band.model_copy(update={
            "number_density_per_mm2": density_per_um2 * 1e6,
            "n_eds_particles": len(band_members) - n_found,
        }))
    return weights, weighted_bands


def bootstrap_confidence_intervals(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand] = (),
    excluded_frame_ids: Iterable[str] = (),
    confidence: float = 0.95,
    n_resamples: int = N_BOOTSTRAP_RESAMPLES,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, Tuple[float, float]]:
    """
    Confidence intervals of the size statistics (:data:`CI_QUANTITIES`), by bootstrap.

    The frames of each size band are resampled with replacement, with their particles and area, so that both the
    size distribution within a band and its number density vary as between measurements of the same sample;
    the statistics are then weighted as in :func:`band_particle_weights`. Particles of a band without frame
    records (and without bands, all particles) are resampled individually. Particles measured by EDS are kept.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles included in the statistics.
    bands : sequence of ParticleSizeBand, optional
        Size bands of the measurement.
    excluded_frame_ids : iterable of str, optional
        Frames left out of the statistics.
    confidence : float, optional
        Confidence level (default 0.95).
    n_resamples : int, optional
        Number of bootstrap resamples.
    rng : numpy.random.Generator, optional
        Random generator (default: a new one).

    Returns
    -------
    dict
        Name (``D10``, ``D50``, ``D90``, ``D43``, ``Dv10``, ``Dv50``, ``Dv90``) to (low, high) bounds in µm; empty
        if a band has fewer than :data:`MIN_BOOTSTRAP_UNITS` frames (or particles, without frame records).
    """
    rng = np.random.default_rng() if rng is None else rng
    d = particle_sizes_dataframe(particles)[cnst.PAR_EQ_D_KEY].to_numpy() if particles else np.array([])
    corrections = np.array([p.edge_correction or 1.0 for p in particles])

    # Resampling units of each band: (area, particle indices), plus the particles kept as they are
    strata = []
    if not bands:
        strata.append(([(1.0, [k]) for k in range(len(particles))], []))
    else:
        excluded = set(excluded_frame_ids)
        members, found_in_band = _band_members(particles, bands)
        for band, band_members, found in zip(bands, members, found_in_band):
            found_set = set(found)
            kept = [k for k in band_members if k not in found_set]
            if band.frames:
                by_frame: Dict[Optional[str], List[int]] = {}
                for k in found:
                    by_frame.setdefault(particles[k].frame_id, []).append(k)
                frames = [f for f in band.frames if f.frame_id not in excluded]
                units = [(f.area_um2, by_frame.pop(f.frame_id, [])) for f in frames]
                for frame_id in excluded:
                    by_frame.pop(frame_id, None)
                units += [(0.0, [k]) for ks in by_frame.values() for k in ks]  # particles of no recorded frame
            else:
                area_um2 = band.scanned_area_um2 - band.frame_area_um2(excluded)
                units = [(area_um2 / len(found), [k]) for k in found] if found else []
            strata.append((units, kept))

    if not strata or any(len(units) < MIN_BOOTSTRAP_UNITS for units, _ in strata):
        return {}
    estimates: Dict[str, List[float]] = {name: [] for name in CI_QUANTITIES}
    for _ in range(n_resamples):
        indices: List[int] = []
        weights: List[float] = []
        for units, kept in strata:
            if not units:
                continue
            picked = [units[j] for j in rng.integers(len(units), size=len(units))]
            area_um2 = sum(area for area, _ in picked)
            resampled = [k for _, ks in picked for k in ks]
            if not resampled or area_um2 <= 0:
                continue
            band_particles = resampled + kept
            density = corrections[resampled].sum() / area_um2
            weights.extend(density * corrections[band_particles] / corrections[band_particles].sum())
            indices.extend(band_particles)
        if len(indices) < 2:
            continue
        stats = compute_particle_size_stats(d[indices], np.array(weights))
        for name, key in CI_QUANTITIES.items():
            estimates[name].append(stats[key])

    tail_percent = (1 - confidence) / 2 * 100
    return {
        name: (float(np.percentile(values, tail_percent)), float(np.percentile(values, 100 - tail_percent)))
        for name, values in estimates.items() if len(values) >= 2
    }


def ci_relative_half_widths(
    size_stats: Dict[str, Union[int, float]],
    confidence_intervals: Dict[str, Tuple[float, float]],
    quantities: Iterable[str],
) -> Dict[str, float]:
    """Half-width of the confidence interval of each of ``quantities`` (e.g. 'D50'), relative to its value."""
    return {
        name: (confidence_intervals[name][1] - confidence_intervals[name][0]) / 2 / size_stats[CI_QUANTITIES[name]]
        if name in confidence_intervals else float("inf")
        for name in quantities
    }


def select_size_stats_particles(
    particles: Sequence[ParticleInfo],
    include_eds_auto_particles: bool = True,
    size_stats_ids: Optional[Iterable[int]] = None,
) -> List[ParticleInfo]:
    """
    Particles included in the size statistics, in ledger order.

    Particles found by a particle size statistics measurement (selection 'size_stats'), and particles measured by
    EDS with automatic particle selection (selection 'auto') if ``include_eds_auto_particles``. Particles chosen
    manually or from a list, without a recorded selection, or without an area are left out.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles of the sample (ledger particles).
    include_eds_auto_particles : bool, optional
        Whether to include the particles measured by EDS with automatic particle selection (default: True).
    size_stats_ids : iterable of int, optional
        IDs of the particles found by the measurement whose statistics are computed (the particles of its size
        bands). Particles of earlier particle size statistics measurements of the sample are then left out.
        Default: all the 'size_stats' particles.
    """
    run_ids = None if size_stats_ids is None else set(size_stats_ids)
    return [
        p for p in particles
        if p.area_um is not None and (
            (p.selection == "size_stats" and (run_ids is None or p.id in run_ids))
            or (include_eds_auto_particles and p.selection == "auto")
        )
    ]


def build_particle_stats_result(
    particles: Sequence[ParticleInfo],
    powder_meas_cfg: PowderMeasurementConfig,
    n_par_target: Optional[int] = None,
    bands: Optional[Sequence[ParticleSizeBand]] = None,
    excluded_particle_ids: Sequence[int] = (),
    excluded_frame_ids: Sequence[str] = (),
) -> ParticleStatsResult:
    """
    Particle size statistics result of the given particles (already selected), to store in the sample ledger.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles included in the statistics (see :func:`select_size_stats_particles`).
    powder_meas_cfg : PowderMeasurementConfig
        Settings of the measurement.
    n_par_target : int, optional
        Target number of particles of the measurement (per size band).
    bands : sequence of ParticleSizeBand, optional
        Size bands searched by the measurement: the statistics are weighted by band (see
        :func:`band_particle_weights`).
    excluded_particle_ids, excluded_frame_ids : sequence, optional
        Particles and frames left out of the statistics (already removed from ``particles``), recorded in the
        result; the area of the frames left out is not counted.

    Raises
    ------
    ValueError
        If no particle has a positive weight (the measurement found no particle in the bands).
    """
    eq_diameters_um = particle_sizes_dataframe(particles)[cnst.PAR_EQ_D_KEY].to_numpy()
    weights, weighted_bands = band_particle_weights(particles, list(bands or []), excluded_frame_ids)
    return ParticleStatsResult(
        n_par_target=n_par_target,
        powder_meas_cfg=powder_meas_cfg.model_copy(deep=True),
        particle_ids=[p.id for p in particles],
        n_par_by_selection=dict(Counter(str(p.selection) for p in particles)),
        size_stats=ParticleSizeStats(**compute_particle_size_stats(eq_diameters_um, weights)),
        unweighted_size_stats=ParticleSizeStats(**compute_particle_size_stats(eq_diameters_um)),
        n_par_without_weight=int(np.sum(weights == 0)),
        confidence_intervals_95=bootstrap_confidence_intervals(particles, list(bands or []), excluded_frame_ids),
        bands=weighted_bands,
        excluded_particle_ids=list(excluded_particle_ids),
        excluded_frame_ids=list(excluded_frame_ids),
    )


def recalculate_particle_stats(
    ledger: SampleLedger,
    particle_ids_to_filter: Optional[Iterable[int]] = None,
    frame_ids_to_filter: Optional[Iterable[str]] = None,
    include_eds_auto_particles: Optional[bool] = None,
) -> Optional[ParticleStatsResult]:
    """
    Particle size statistics of a sample recalculated from its ledger particles, leaving out some particles or frames.

    Uses the particles found by the particle size statistics measurement of the sample (the particles of its size
    bands, or all 'size_stats' particles if it has none) and, if ``include_eds_auto_particles``, the particles
    measured by EDS with automatic particle selection, including those measured after it. The settings, target
    and bands of the measurement are kept.

    Parameters
    ----------
    ledger : SampleLedger
        Ledger of the sample (not modified).
    particle_ids_to_filter : iterable of int, optional
        IDs of the particles to leave out.
    frame_ids_to_filter : iterable of str, optional
        IDs of the frames whose particles are left out (exact matches).
    include_eds_auto_particles : bool, optional
        Default: the setting of the measurement.

    Returns
    -------
    ParticleStatsResult or None
        The statistics, recording the particles and frames left out (those found among the particles), or None if
        no particle is left (or none found by the measurement in its size bands).
    """
    measurement_result = ledger.particle_stats
    powder_meas_cfg = (
        measurement_result.powder_meas_cfg if measurement_result is not None
        else ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig()
    )
    if include_eds_auto_particles is None:
        include_eds_auto_particles = powder_meas_cfg.include_eds_auto_particles
    bands = measurement_result.bands if measurement_result is not None else []
    size_stats_ids = [pid for band in bands for pid in band.particle_ids] if bands else None
    particles = select_size_stats_particles(ledger.particles, include_eds_auto_particles, size_stats_ids)

    particle_ids_to_filter = [int(pid) for pid in (particle_ids_to_filter or [])]
    frame_ids_to_filter = [str(frame_id) for frame_id in (frame_ids_to_filter or [])]
    excluded_particle_ids = [pid for pid in particle_ids_to_filter if pid in {p.id for p in particles}]
    known_frame_ids = {p.frame_id for p in particles} | {f.frame_id for band in bands for f in band.frames}
    excluded_frame_ids = [fid for fid in frame_ids_to_filter if fid in known_frame_ids]
    particles = [p for p in particles if p.id not in particle_ids_to_filter and p.frame_id not in frame_ids_to_filter]
    if not particles:
        return None

    weights, _ = band_particle_weights(particles, bands, excluded_frame_ids)
    if not np.any(weights > 0):
        logger.warning("No particle found by the measurement is left in its size bands: no statistics.")
        return None
    result = build_particle_stats_result(
        particles,
        powder_meas_cfg.model_copy(update={"include_eds_auto_particles": include_eds_auto_particles}),
        n_par_target=measurement_result.n_par_target if measurement_result is not None else None,
        bands=bands,
        excluded_particle_ids=excluded_particle_ids,
        excluded_frame_ids=excluded_frame_ids,
    )
    return result


def particle_sizes_dataframe(particles: Sequence[ParticleInfo]) -> pd.DataFrame:
    """
    Table of particle sizes: ID, frame, area (µm²) and equivalent diameter (µm) of each particle.

    The equivalent diameter is computed from the area when a particle does not record it.

    Raises
    ------
    ValueError
        If a particle has no area.
    """
    if any(p.area_um is None for p in particles):
        raise ValueError("All ParticleInfo records must have area_um set for PSD export")
    return pd.DataFrame({
        cnst.PAR_ID_DF_KEY: np.array([p.id for p in particles], dtype=int),
        cnst.FRAME_ID_DF_KEY: np.array([p.frame_id if p.frame_id is not None else "" for p in particles], dtype=str),
        cnst.PAR_AREA_UM_KEY: np.array([float(p.area_um) for p in particles], dtype=float),
        cnst.PAR_EQ_D_KEY: np.array(
            [
                float(p.eq_diameter_um) if p.eq_diameter_um is not None else float(2.0 * np.sqrt(float(p.area_um) / np.pi))
                for p in particles
            ],
            dtype=float,
        ),
    })


def save_particle_size_histogram(
    eq_diameters_um: Sequence[float],
    output_dir: str,
    sample_id: str,
    output_file_suffix: str = "",
    weights: Optional[Sequence[float]] = None,
    size_stats: Optional[Dict[str, Union[int, float]]] = None,
    image_format: str = "png",
    bins: int = 20,
) -> str:
    """
    Save the particle size distribution plot as ``<sample_id>_Par_size_distribution_hist<suffix>.<image_format>``.

    Histogram of the number fraction of particles (weighted, if ``weights`` are given) by equivalent diameter,
    with logarithmic bins when the diameters span more than a factor 10, and cumulative number and volume
    distributions; the D10, D50 and D90 of ``size_stats`` (see :func:`compute_particle_size_stats`) are marked.

    Returns
    -------
    str
        Path of the saved image.
    """
    d = np.asarray(eq_diameters_um, dtype=float)
    w = np.ones(d.size) if weights is None else np.asarray(weights, dtype=float)
    d, w = d[w > 0], w[w > 0]
    fractions_percent = w / np.sum(w) * 100

    d_min, d_max = float(d.min()), float(d.max())
    is_log_scale = d_max / d_min > 10
    if np.isclose(d_min, d_max):
        edges = np.array([d_min * 0.95, d_max * 1.05])
    elif is_log_scale:
        edges = np.geomspace(d_min, d_max, bins + 1)
    else:
        edges = np.linspace(d_min, d_max, bins + 1)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(d, bins=edges, weights=fractions_percent, color="tab:blue", alpha=0.6, edgecolor="black")
    if is_log_scale:
        ax.set_xscale("log")
    ax.set_xlabel("Equivalent Diameter (μm)")
    ax.set_ylabel("Number fraction (%)")
    ax.set_title("Particle size distribution" + (" (weighted by size band and frame edge)" if weights is not None else ""))

    order = np.argsort(d, kind="stable")
    ax_cumulative = ax.twinx()
    for distribution_weights, style, label in ((w, "-", "Cumulative, number"), (w * d**3, "--", "Cumulative, volume")):
        cumulative = np.cumsum(distribution_weights[order]) / np.sum(distribution_weights) * 100
        ax_cumulative.step(d[order], cumulative, where="post", color="black", linestyle=style, label=label)
    ax_cumulative.set_ylim(0, 100)
    ax_cumulative.set_ylabel("Cumulative (%)")
    if size_stats is not None:
        for height, (key, label) in zip((0.98, 0.91, 0.84), (("D10", "D10"), ("median", "D50"), ("D90", "D90"))):
            ax.axvline(size_stats[key], color="tab:red", linestyle=":")
            ax.text(size_stats[key], height, f" {label}", color="tab:red", transform=ax.get_xaxis_transform(),
                    va="top", ha="left", fontsize=8)
    ax_cumulative.legend(loc="center right", fontsize=8)

    out_path = os.path.join(
        output_dir, f'{sample_id}_{cnst.PARTICLE_STAT_HIST_FILENAME}{output_file_suffix}.{image_format}'
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def save_particle_statistics(
    particles: Sequence[ParticleInfo],
    output_dir: str,
    sample_id: str,
    output_file_suffix: str = "",
    verbose: bool = False,
    bands: Optional[Sequence[ParticleSizeBand]] = None,
    excluded_frame_ids: Iterable[str] = (),
    image_format: str = "png",
    confidence_intervals: Optional[Dict[str, Tuple[float, float]]] = None,
) -> pd.DataFrame:
    """
    Compute the particle size statistics and save the particle sizes, the statistics and the histogram.

    With size bands, the statistics are weighted by band (see :func:`band_particle_weights`): the particle sizes
    file then gives the weight of each particle (relative, 1 on average over the particles counted), and the
    statistics file the number of bands.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles included in the statistics (each with ``area_um``).
    output_dir : str
        Folder where the files are saved (the sample folder).
    sample_id : str
        Sample ID, prefix of the file names.
    output_file_suffix : str, optional
        String added to the file names.
    verbose : bool, optional
        If True, logs the statistics.
    bands : sequence of ParticleSizeBand, optional
        Size bands of the measurement, to weight the statistics.
    excluded_frame_ids : iterable of str, optional
        Frames left out of the statistics, whose area is not counted.
    image_format : str, optional
        File format of the size distribution plot (default: 'png').
    confidence_intervals : dict, optional
        95 % confidence intervals of the statistics (default: computed, see
        :func:`bootstrap_confidence_intervals`), saved as ``<name>_CI95_low`` and ``<name>_CI95_high``.

    Returns
    -------
    pandas.DataFrame
        One row: the measure described, followed by the statistics of :func:`compute_particle_size_stats` (and
        ``n_bands`` with size bands).

    Raises
    ------
    ValueError
        If no particle is given, a particle has no area, or no particle has a positive weight.
    """
    if not particles:
        raise ValueError("No particles available to save statistics from")

    particle_data = particle_sizes_dataframe(particles)
    areas_um = particle_data[cnst.PAR_AREA_UM_KEY].to_numpy()
    eq_diam_um = particle_data[cnst.PAR_EQ_D_KEY].to_numpy()

    if len(areas_um) > 1 and np.isclose(areas_um.min(), np.partition(areas_um, 1)[1]):
        logger.warning(
            '⚠ The 2 smallest particles have identical area.\n'
            '   This may indicate they are only a few pixels in size.\n'
            '   Consider increasing min_par_pixels, the minimum number of pixels\n'
            '   of a particle for its area to be measured.'
        )

    weights = None
    if bands:
        weights, _ = band_particle_weights(particles, list(bands), excluded_frame_ids)
        n_weighted = int(np.sum(weights > 0))
        particle_data[cnst.PAR_WEIGHT_KEY] = weights * n_weighted / np.sum(weights) if n_weighted else weights
        if n_weighted < len(particles):
            logger.warning(
                f"⚠️ {len(particles) - n_weighted} particle(s) measured by EDS are outside the size bands of the "
                "measurement (or in a band where it found no particle): left out of the statistics."
            )

    particle_data.to_csv(
        os.path.join(output_dir, f"{sample_id}_{cnst.PARTICLE_SIZES_FILENAME}{output_file_suffix}.csv"),
        header=True,
        index=False,
    )

    stats = compute_particle_size_stats(eq_diam_um, weights)
    if confidence_intervals is None:
        confidence_intervals = bootstrap_confidence_intervals(particles, list(bands or []), excluded_frame_ids)
    for name, (low, high) in confidence_intervals.items():
        stats[f"{name}_CI95_low"], stats[f"{name}_CI95_high"] = low, high
    if bands:
        stats["n_bands"] = len(bands)
    stats_df = pd.DataFrame({'Measure': SIZE_STATS_MEASURE, **stats}, index=[0])
    stats_df.to_csv(
        os.path.join(output_dir, f"{sample_id}_{cnst.PARTICLE_STATS_FILENAME}{output_file_suffix}.csv"),
        header=True,
        index=False,
    )

    if verbose:
        print_double_separator("📊 Computed statistics:")
        logger.info(stats_df.T.to_string(header=False))

    save_particle_size_histogram(
        eq_diam_um, output_dir, sample_id, output_file_suffix=output_file_suffix, weights=weights,
        size_stats=stats, image_format=image_format,
    )
    return stats_df
