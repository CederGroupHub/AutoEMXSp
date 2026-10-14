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
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple, Union

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
    "SizeObservations",
    "band_particle_weights",
    "bootstrap_confidence_intervals",
    "edge_correction_factor",
    "build_particle_stats_result",
    "compute_particle_size_stats",
    "overlapping_frame_ids",
    "particle_sizes_dataframe",
    "recalculate_particle_stats",
    "save_particle_size_histogram",
    "save_particle_statistics",
    "select_size_stats_particles",
    "size_observations",
    "weighted_size_observations",
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


class SizeObservations(NamedTuple):
    """
    Particle sizes from the scans of a sample (size bands and frames searched by EDS), with what is needed to
    weight them (see :func:`size_observations`).
    """

    #: Equivalent diameter (µm), area (µm²) and frame-edge correction of each observation.
    eq_diameters_um: np.ndarray
    areas_um2: np.ndarray
    corrections: np.ndarray
    #: Scan of each observation: index in ``scans`` (size bands first, then EDS scans), or -1 for a particle
    #: measured by EDS that only shares the weight of the band of its diameter (no EDS scan recorded).
    scan: np.ndarray
    #: Frame of each observation, and index of its particle in the particles given (None for EDS frame counts).
    frame_ids: List[Optional[str]]
    particle_index: List[Optional[int]]
    #: Size bands of the measurement, then EDS scans.
    scans: List[ParticleSizeBand]
    n_bands: int


def size_observations(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand] = (),
    eds_scans: Sequence[ParticleSizeBand] = (),
    excluded_frame_ids: Iterable[str] = (),
) -> SizeObservations:
    """
    Observations of particle sizes: the particles given (those found by the measurement belong to their size band,
    the others share the weight of the band of their diameter) and the particles counted in the frames searched by
    EDS (``eds_scans``, see ``SampleLedger.eds_size_scans``), except in the frames left out: ``excluded_frame_ids``
    and frames overlapping a larger one (see :func:`overlapping_frame_ids`).
    """
    excluded = _excluded_frames(list(bands) + list(eds_scans), excluded_frame_ids)
    measurement_band = {pid: i for i, band in enumerate(bands) for pid in band.particle_ids}
    kept_index = [k for k, p in enumerate(particles) if not (p.id in measurement_band and p.frame_id in excluded)]
    kept = [particles[k] for k in kept_index]
    sizes = particle_sizes_dataframe(kept) if kept else None
    eq_diameters = list(sizes[cnst.PAR_EQ_D_KEY]) if sizes is not None else []
    areas = list(sizes[cnst.PAR_AREA_UM_KEY]) if sizes is not None else []
    corrections = [p.edge_correction or 1.0 for p in kept]
    scan = [measurement_band.get(p.id, -1) for p in kept]
    frame_ids: List[Optional[str]] = [p.frame_id for p in kept]
    particle_index: List[Optional[int]] = list(kept_index)
    for s, eds_scan in enumerate(eds_scans):
        for frame in eds_scan.frames:
            if frame.frame_id in excluded:
                continue
            for area, correction in zip(frame.particle_areas_um2 or [], frame.particle_edge_corrections or []):
                eq_diameters.append(float(2.0 * np.sqrt(area / np.pi)))
                areas.append(float(area))
                corrections.append(float(correction))
                scan.append(len(bands) + s)
                frame_ids.append(frame.frame_id)
                particle_index.append(None)
    return SizeObservations(
        np.asarray(eq_diameters, dtype=float), np.asarray(areas, dtype=float), np.asarray(corrections, dtype=float),
        np.asarray(scan, dtype=int), frame_ids, particle_index, list(bands) + list(eds_scans), len(bands),
    )


def overlapping_frame_ids(scans: Sequence[ParticleSizeBand]) -> List[str]:
    """
    Frames left out of the statistics because they overlap, on the sample, a frame of a scan with an overlapping
    size range (which counts the same particles): of overlapping frames, the largest is kept (the first recorded,
    if as large). Frames of different size bands, whose size ranges do not overlap, are never compared. Frames
    without position are kept.
    """
    frames = [(s, f) for s, scan in enumerate(scans) for f in scan.frames if f.center_mm is not None]
    order = sorted(range(len(frames)), key=lambda i: (-frames[i][1].frame_width_um, i))
    centers = np.array([[f.center_mm.x * 1000, f.center_mm.y * 1000] for _, f in frames]).reshape(-1, 2)
    half_sizes = np.array([[f.frame_width_um / 2, f.area_um2 / f.frame_width_um / 2] for _, f in frames]).reshape(-1, 2)
    d_ranges = np.array([[scans[s].d_min_um, scans[s].d_max_um] for s, _ in frames]).reshape(-1, 2)
    kept: List[int] = []
    dropped: List[str] = []
    for i in order:
        if kept:
            k = np.array(kept)
            same_sizes = (d_ranges[k, 0] < d_ranges[i, 1] * (1 - 1e-9)) & (d_ranges[i, 0] < d_ranges[k, 1] * (1 - 1e-9))
            overlap = np.all(np.abs(centers[k] - centers[i]) < (half_sizes[k] + half_sizes[i]) * (1 - 1e-3), axis=1)
            if np.any(same_sizes & overlap):
                dropped.append(frames[i][1].frame_id)
                continue
        kept.append(i)
    return dropped


def _excluded_frames(scans: Sequence[ParticleSizeBand], excluded_frame_ids: Iterable[str] = ()) -> set:
    """Frames left out of the statistics: those given, and those overlapping a larger one."""
    return set(excluded_frame_ids) | set(overlapping_frame_ids(scans))


def _in_size_range(d: np.ndarray, scan: ParticleSizeBand) -> np.ndarray:
    return (d >= scan.d_min_um * (1 - 1e-9)) & (d <= scan.d_max_um * (1 + 1e-9))


def _band_of_diameter(d: np.ndarray, bands: Sequence[ParticleSizeBand]) -> np.ndarray:
    """Index of the size band of each diameter (the larger one at a boundary), or -1."""
    band = np.full(len(d), -1)
    for b in range(len(bands) - 1, -1, -1):
        band = np.where(_in_size_range(d, bands[b]), b, band)
    return band


def _observation_weights(obs: SizeObservations, indices: np.ndarray, scan_areas_um2: np.ndarray) -> np.ndarray:
    """
    Weights of the observations ``indices`` (possibly repeated, as in a bootstrap resample), given the area
    scanned by each scan.

    Each observation counts its edge correction over the area of the scans able to observe its size: its own size
    band (size bands do not overlap) or the band of its diameter, and every EDS scan whose size range includes it.
    Particles measured by EDS without EDS scan (``scan == -1``) then share the weight of the band of their
    diameter, in proportion to their correction.
    """
    d, c, scan = obs.eq_diameters_um[indices], obs.corrections[indices], obs.scan[indices]
    n_bands = obs.n_bands
    # Size band of each observation: its own, or that of its diameter (the larger one at a boundary)
    band_of_diameter = _band_of_diameter(d, obs.scans[:n_bands])
    band = np.where((scan >= 0) & (scan < n_bands), scan, np.where(scan >= n_bands, band_of_diameter, -1))
    member_band = np.where(scan == -1, band_of_diameter, -1)

    observed_area = np.where(band >= 0, scan_areas_um2[np.maximum(band, 0)], 0.0)
    for s in range(n_bands, len(obs.scans)):
        observed_area = observed_area + np.where((scan == s) | _in_size_range(d, obs.scans[s]), scan_areas_um2[s], 0.0)
    observed_area = np.where(scan == -1, 0.0, observed_area)
    weights = np.where(observed_area > 0, c / np.where(observed_area > 0, observed_area, 1.0), 0.0)

    # Particles measured by EDS without EDS scan: share the weight of their band
    for b in range(n_bands):
        sharing = (member_band == b) | ((scan == b) & np.any(member_band == b))
        if np.any(member_band == b):
            band_weight = weights[scan == b].sum()
            weights[sharing] = band_weight * c[sharing] / c[sharing].sum()
    return weights


def _scan_areas_um2(scans: Sequence[ParticleSizeBand], excluded_frame_ids: Iterable[str] = ()) -> np.ndarray:
    excluded = list(excluded_frame_ids)
    return np.array([scan.scanned_area_um2 - scan.frame_area_um2(excluded) for scan in scans], dtype=float)


def weighted_size_observations(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand] = (),
    excluded_frame_ids: Iterable[str] = (),
    eds_scans: Sequence[ParticleSizeBand] = (),
) -> Tuple[SizeObservations, np.ndarray, List[ParticleSizeBand], List[ParticleSizeBand]]:
    """
    Observations of particle sizes (see :func:`size_observations`) and their weights, so that the statistics
    describe the sample: each counts its edge correction (:func:`edge_correction_factor`) over the area of the
    scans able to observe its size, so that each size counts in proportion to its number of particles per unit
    area of the sample, whatever the number of particles collected in each size band.

    Returns
    -------
    observations : SizeObservations
    weights : numpy.ndarray
        Weight of each observation (all 1 without bands or EDS scans).
    bands, eds_scans : list of ParticleSizeBand
        The bands and EDS scans with their ``number_density_per_mm2`` (particles they counted per mm²) and, for
        bands, ``n_eds_particles`` (particles measured by EDS sharing their weight).
    """
    obs = size_observations(particles, bands, eds_scans, excluded_frame_ids)
    if not obs.scans:
        return obs, np.ones(len(obs.eq_diameters_um)), [], []
    areas = _scan_areas_um2(obs.scans, _excluded_frames(obs.scans, excluded_frame_ids))
    weights = _observation_weights(obs, np.arange(len(obs.eq_diameters_um)), areas)
    member_band = _band_of_diameter(obs.eq_diameters_um, obs.scans[:obs.n_bands])
    weighted_scans = []
    for s, scan in enumerate(obs.scans):
        counted = obs.corrections[obs.scan == s].sum()
        update = {"number_density_per_mm2": float(counted / areas[s] * 1e6) if areas[s] > 0 else 0.0}
        if s < obs.n_bands:
            update["n_eds_particles"] = int(np.sum((obs.scan == -1) & (member_band == s)))
        weighted_scans.append(scan.model_copy(update=update))
    return obs, weights, weighted_scans[:obs.n_bands], weighted_scans[obs.n_bands:]


def band_particle_weights(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand],
    excluded_frame_ids: Iterable[str] = (),
) -> Tuple[np.ndarray, List[ParticleSizeBand]]:
    """
    Weights of the particles given, pooled from several size bands (see :func:`weighted_size_observations`, without
    EDS scans).

    Returns
    -------
    weights : numpy.ndarray
        Weight of each particle (all 1 without bands).
    bands : list of ParticleSizeBand
        The bands with their ``number_density_per_mm2`` and ``n_eds_particles``.
    """
    _, weights, weighted_bands, _ = weighted_size_observations(particles, bands, excluded_frame_ids)
    return weights, weighted_bands


def bootstrap_confidence_intervals(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand] = (),
    excluded_frame_ids: Iterable[str] = (),
    confidence: float = 0.95,
    n_resamples: int = N_BOOTSTRAP_RESAMPLES,
    rng: Optional[np.random.Generator] = None,
    eds_scans: Sequence[ParticleSizeBand] = (),
) -> Dict[str, Tuple[float, float]]:
    """
    Confidence intervals of the size statistics (:data:`CI_QUANTITIES`), by bootstrap.

    The frames of each scan (size band, or frames searched by EDS) are resampled with replacement, with their
    particles and area, so that both the size distribution and the number density of each scan vary as between
    measurements of the same sample; the statistics are then weighted as in :func:`weighted_size_observations`.
    Particles of a scan without frame records (and without scans, all particles) are resampled individually.
    Particles measured by EDS sharing the weight of a band are kept.

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
    eds_scans : sequence of ParticleSizeBand, optional
        Frames searched by EDS included in the statistics.

    Returns
    -------
    dict
        Name (``D10``, ``D50``, ``D90``, ``D43``, ``Dv10``, ``Dv50``, ``Dv90``) to (low, high) bounds in µm; empty
        if a scan has fewer than :data:`MIN_BOOTSTRAP_UNITS` frames (or particles, without frame records).
    """
    rng = np.random.default_rng() if rng is None else rng
    excluded = _excluded_frames(list(bands) + list(eds_scans), excluded_frame_ids)
    obs = size_observations(particles, bands, eds_scans, excluded)
    kept = np.flatnonzero(obs.scan == -1)

    # Resampling units of each scan: (area, observation indices)
    strata = []
    if not obs.scans:
        strata.append([(1.0, [k]) for k in range(len(obs.eq_diameters_um))])
    else:
        for s, scan in enumerate(obs.scans):
            observed = np.flatnonzero(obs.scan == s)
            if scan.frames:
                by_frame: Dict[Optional[str], List[int]] = {}
                for k in observed:
                    by_frame.setdefault(obs.frame_ids[k], []).append(int(k))
                units = [(f.area_um2, by_frame.pop(f.frame_id, [])) for f in scan.frames if f.frame_id not in excluded]
                for frame_id in excluded:
                    by_frame.pop(frame_id, None)
                units += [(0.0, [k]) for ks in by_frame.values() for k in ks]  # particles of no recorded frame
            else:
                area_um2 = scan.scanned_area_um2 - scan.frame_area_um2(excluded)
                units = [(area_um2 / len(observed), [int(k)]) for k in observed] if len(observed) else []
            strata.append(units)
    if not strata or any(len(units) < MIN_BOOTSTRAP_UNITS for units in strata):
        return {}

    estimates: Dict[str, List[float]] = {name: [] for name in CI_QUANTITIES}
    for _ in range(n_resamples):
        indices: List[int] = list(kept)
        scan_areas = np.zeros(len(strata))
        for s, units in enumerate(strata):
            picked = [units[j] for j in rng.integers(len(units), size=len(units))]
            scan_areas[s] = sum(area for area, _ in picked)
            indices.extend(k for _, ks in picked for k in ks)
        indices_array = np.asarray(indices, dtype=int)
        if not obs.scans:
            weights = np.ones(len(indices_array))
        else:
            weights = _observation_weights(obs, indices_array, scan_areas)
        if np.sum(weights > 0) < 2:
            continue
        stats = compute_particle_size_stats(obs.eq_diameters_um[indices_array], weights)
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
    eds_scans: Sequence[ParticleSizeBand] = (),
) -> ParticleStatsResult:
    """
    Particle size statistics result of the given particles (already selected), to store in the sample ledger,
    weighted as in :func:`weighted_size_observations`.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles included in the statistics (see :func:`select_size_stats_particles`).
    powder_meas_cfg : PowderMeasurementConfig
        Settings of the measurement.
    n_par_target : int, optional
        Target number of particles of the measurement (per size band).
    bands : sequence of ParticleSizeBand, optional
        Size bands searched by the measurement.
    excluded_particle_ids, excluded_frame_ids : sequence, optional
        Particles and frames left out of the statistics (already removed from ``particles``), recorded in the
        result; the area of the frames left out is not counted.
    eds_scans : sequence of ParticleSizeBand, optional
        Frames searched by EDS included in the statistics (``SampleLedger.eds_size_scans``).

    Raises
    ------
    ValueError
        If no observation has a positive weight.
    """
    obs, weights, weighted_bands, weighted_eds_scans = weighted_size_observations(
        particles, list(bands or []), excluded_frame_ids, list(eds_scans)
    )
    used = [particles[k] for k in obs.particle_index if k is not None]  # those of frames left out excepted
    n_par_by_selection = dict(Counter(str(p.selection) for p in used))
    n_eds_frame_particles = int(np.sum(obs.scan >= obs.n_bands))
    if n_eds_frame_particles:
        n_par_by_selection["eds_frames"] = n_eds_frame_particles
    return ParticleStatsResult(
        n_par_target=n_par_target,
        powder_meas_cfg=powder_meas_cfg.model_copy(deep=True),
        particle_ids=[p.id for p in used],
        n_par_by_selection=n_par_by_selection,
        size_stats=ParticleSizeStats(**compute_particle_size_stats(obs.eq_diameters_um, weights)),
        unweighted_size_stats=ParticleSizeStats(**compute_particle_size_stats(obs.eq_diameters_um)),
        n_par_without_weight=int(np.sum(weights == 0)),
        confidence_intervals_95=bootstrap_confidence_intervals(
            particles, list(bands or []), excluded_frame_ids, eds_scans=list(eds_scans)
        ),
        bands=weighted_bands,
        eds_scans=weighted_eds_scans,
        excluded_particle_ids=list(excluded_particle_ids),
        excluded_frame_ids=list(excluded_frame_ids),
        overlapping_frame_ids=overlapping_frame_ids(list(bands or []) + list(eds_scans)),
    )


def recalculate_particle_stats(
    ledger: SampleLedger,
    particle_ids_to_filter: Optional[Iterable[int]] = None,
    frame_ids_to_filter: Optional[Iterable[str]] = None,
    include_eds_auto_particles: Optional[bool] = None,
) -> Optional[ParticleStatsResult]:
    """
    Particle size statistics of a sample recalculated from its ledger, leaving out some particles or frames.

    Uses the particles found by the particle size statistics measurement of the sample (the particles of its size
    bands, or all 'size_stats' particles if it has none) and, if ``include_eds_auto_particles``, what EDS measured
    on the sample with automatic particle selection, before or after it: the particles counted in the frames
    searched by EDS (``SampleLedger.eds_size_scans``) or, for samples measured before these were recorded, the
    particles measured by EDS. A sample measured only by EDS gets statistics from its EDS frames. The settings,
    target and bands of the measurement are kept.

    Parameters
    ----------
    ledger : SampleLedger
        Ledger of the sample (not modified).
    particle_ids_to_filter : iterable of int, optional
        IDs of the particles to leave out.
    frame_ids_to_filter : iterable of str, optional
        IDs of the frames whose particles are left out (exact matches; frames searched by EDS are ``EDS-<label>``).
    include_eds_auto_particles : bool, optional
        Default: the setting of the measurement (of the EDS acquisition, for a sample measured only by EDS).

    Returns
    -------
    ParticleStatsResult or None
        The statistics, recording the particles and frames left out (those found), or None if nothing is left.
    """
    measurement_result = ledger.particle_stats
    has_measurement = measurement_result is not None and bool(measurement_result.bands)
    powder_meas_cfg = (
        measurement_result.powder_meas_cfg if measurement_result is not None
        else ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig()
    )
    if include_eds_auto_particles is None:
        include_eds_auto_particles = powder_meas_cfg.include_eds_auto_particles
    bands = measurement_result.bands if measurement_result is not None else []
    eds_scans = list(ledger.eds_size_scans) if include_eds_auto_particles else []
    size_stats_ids = [pid for band in bands for pid in band.particle_ids] if bands else None
    # EDS frames count the particles measured by EDS: these are then left out
    particles = select_size_stats_particles(
        ledger.particles, include_eds_auto_particles and not eds_scans, size_stats_ids
    )
    if not has_measurement and eds_scans:
        particles = []  # only EDS: its frames

    particle_ids_to_filter = [int(pid) for pid in (particle_ids_to_filter or [])]
    frame_ids_to_filter = [str(frame_id) for frame_id in (frame_ids_to_filter or [])]
    excluded_particle_ids = [pid for pid in particle_ids_to_filter if pid in {p.id for p in particles}]
    known_frame_ids = {p.frame_id for p in particles} | {
        f.frame_id for scan in list(bands) + eds_scans for f in scan.frames
    }
    excluded_frame_ids = [fid for fid in frame_ids_to_filter if fid in known_frame_ids]
    particles = [p for p in particles if p.id not in particle_ids_to_filter and p.frame_id not in frame_ids_to_filter]

    obs, weights, _, _ = weighted_size_observations(particles, bands, excluded_frame_ids, eds_scans)
    if not np.any(weights > 0):
        logger.warning("No particle left in the scans of the sample: no statistics.")
        return None
    return build_particle_stats_result(
        particles,
        powder_meas_cfg.model_copy(update={"include_eds_auto_particles": include_eds_auto_particles}),
        n_par_target=measurement_result.n_par_target if measurement_result is not None else None,
        bands=bands,
        excluded_particle_ids=excluded_particle_ids,
        excluded_frame_ids=excluded_frame_ids,
        eds_scans=eds_scans,
    )


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
    eds_scans: Sequence[ParticleSizeBand] = (),
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
    eds_scans : sequence of ParticleSizeBand, optional
        Frames searched by EDS included in the statistics: their particles are listed without particle ID.

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
    obs, weights, _, _ = weighted_size_observations(particles, list(bands or []), excluded_frame_ids, list(eds_scans))
    if len(obs.eq_diameters_um) == 0:
        raise ValueError("No particles available to save statistics from")
    if not (bands or eds_scans):
        weights = None
    areas_um = obs.areas_um2
    eq_diam_um = obs.eq_diameters_um

    if len(areas_um) > 1 and np.isclose(areas_um.min(), np.partition(areas_um, 1)[1]):
        logger.warning(
            '⚠ The 2 smallest particles have identical area.\n'
            '   This may indicate they are only a few pixels in size.\n'
            '   Consider increasing min_par_pixels, the minimum number of pixels\n'
            '   of a particle for its area to be measured.'
        )

    # Particles, and particles counted in frames searched by EDS (without particle ID)
    particle_data = pd.DataFrame({
        cnst.PAR_ID_DF_KEY: pd.array([particles[k].id if k is not None else None for k in obs.particle_index],
                                     dtype="Int64"),
        cnst.FRAME_ID_DF_KEY: np.array([f if f is not None else "" for f in obs.frame_ids], dtype=str),
        cnst.PAR_AREA_UM_KEY: areas_um,
        cnst.PAR_EQ_D_KEY: eq_diam_um,
    })
    if weights is not None:
        n_weighted = int(np.sum(weights > 0))
        particle_data[cnst.PAR_WEIGHT_KEY] = weights * n_weighted / np.sum(weights) if n_weighted else weights
        if n_weighted < len(weights):
            logger.warning(
                f"⚠️ {len(weights) - n_weighted} particle(s) measured by EDS are outside the size ranges scanned "
                "(or in a size band where the measurement found no particle): left out of the statistics."
            )

    particle_data.to_csv(
        os.path.join(output_dir, f"{sample_id}_{cnst.PARTICLE_SIZES_FILENAME}{output_file_suffix}.csv"),
        header=True,
        index=False,
    )

    stats = compute_particle_size_stats(eq_diam_um, weights)
    if confidence_intervals is None:
        confidence_intervals = bootstrap_confidence_intervals(
            particles, list(bands or []), excluded_frame_ids, eds_scans=list(eds_scans)
        )
    for name, (low, high) in confidence_intervals.items():
        stats[f"{name}_CI95_low"], stats[f"{name}_CI95_high"] = low, high
    if bands:
        stats["n_bands"] = len(bands)
    if eds_scans:
        stats["n_eds_scans"] = len(eds_scans)
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
