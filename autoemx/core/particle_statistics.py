#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle size statistics: summary statistics of the equivalent diameters of particles, and their export.

The statistics use the frames of a particle size statistics measurement (``SampleLedger.size_measurements``) and,
if ``PowderMeasurementConfig.include_eds_auto_particles`` (default), the frames searched by EDS on the same sample,
before or after it, each frame under a segmentation (``SampleLedger.frames``, ``SampleLedger.segmentations``). Each
result is stored in the sample ledger (``SampleLedger.size_results``): results of other segmentations, frames or
measurements are kept. For samples measured before frames searched by EDS were recorded, the particles measured by
EDS with automatic particle selection (selection 'auto') count instead.

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
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import autoemx.utils.constants as cnst
from autoemx._logging import get_logger
from autoemx.config.ledger_schemas import (
    BandSummary,
    FrameRecord,
    ParticleInfo,
    ParticleKey,
    ParticleSizeBand,
    ParticleSizeStats,
    ParticleSizeStatsResult,
    ParticleStatsFrame,
    SampleLedger,
)
from autoemx.config.runtime_configs import PowderMeasurementConfig
from autoemx.utils.helper import print_double_separator

logger = get_logger(__name__)

__all__ = [
    "SizeObservations",
    "size_result_label",
    "bootstrap_confidence_intervals",
    "EdsSearchFrame",
    "edge_correction_factor",
    "eds_search_frames",
    "frames_overlap",
    "SizeStatsInputs",
    "build_size_result",
    "default_result_inputs",
    "eds_segmentation_for",
    "size_distribution",
    "working_frame",
    "compute_particle_size_stats",
    "overlapping_frame_ids",
    "particle_sizes_dataframe",
    "recalculate_particle_stats",
    "save_particle_size_histogram",
    "save_particle_statistics",
    "select_size_stats_particles",
    "size_stats_inputs",
    "size_ranges_overlap",
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
    #: Frame of each observation, index of its particle in the particles given (only for those), and its
    #: (frame_id, index) in its frame (only for particles of frames).
    frame_ids: List[Optional[str]]
    particle_index: List[Optional[int]]
    particle_keys: List[Optional[ParticleKey]]
    #: Size bands of the measurement, then EDS scans.
    scans: List[ParticleSizeBand]
    n_bands: int


def size_observations(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand] = (),
    eds_scans: Sequence[ParticleSizeBand] = (),
    excluded_frame_ids: Iterable[str] = (),
    excluded_particles: Iterable[ParticleKey] = (),
) -> SizeObservations:
    """
    Observations of particle sizes: the particles counted in the frames of the size bands and in the frames
    searched by EDS (``eds_scans``), and the particles given (measured by EDS, sharing the weight of the band of
    their diameter), except in the frames left out (``excluded_frame_ids`` and frames overlapping a larger one, see
    :func:`overlapping_frame_ids`) and the particles left out (``excluded_particles``, as (frame_id, index)).
    """
    excluded_frame_ids = set(excluded_frame_ids)
    excluded = _excluded_frames(list(bands) + list(eds_scans), excluded_frame_ids)
    excluded_keys = {(str(f), int(i)) for f, i in excluded_particles}
    kept_index = [k for k, p in enumerate(particles) if p.area_um is not None and p.frame_id not in excluded_frame_ids]
    kept = [particles[k] for k in kept_index]
    eq_diameters = [float(p.eq_diameter_um) if p.eq_diameter_um is not None else float(2.0 * np.sqrt(p.area_um / np.pi))
                    for p in kept]
    areas = [float(p.area_um) for p in kept]
    corrections = [p.edge_correction or 1.0 for p in kept]
    scan = [-1] * len(kept)
    frame_ids: List[Optional[str]] = [p.frame_id for p in kept]
    particle_index: List[Optional[int]] = list(kept_index)
    particle_keys: List[Optional[ParticleKey]] = [None] * len(kept)
    for s, frame_scan in enumerate(list(bands) + list(eds_scans)):
        for frame in frame_scan.frames:
            if frame.frame_id in excluded:
                continue
            for i, (area, correction) in enumerate(zip(frame.particle_areas_um2 or [],
                                                       frame.particle_edge_corrections or [])):
                if (frame.frame_id, i) in excluded_keys:
                    continue
                eq_diameters.append(float(2.0 * np.sqrt(area / np.pi)))
                areas.append(float(area))
                corrections.append(float(correction))
                scan.append(s)
                frame_ids.append(frame.frame_id)
                particle_index.append(None)
                particle_keys.append((frame.frame_id, i))
    return SizeObservations(
        np.asarray(eq_diameters, dtype=float), np.asarray(areas, dtype=float), np.asarray(corrections, dtype=float),
        np.asarray(scan, dtype=int), frame_ids, particle_index, particle_keys, list(bands) + list(eds_scans),
        len(bands),
    )


def frames_overlap(center_um, half_size_um, centers_um, half_sizes_um) -> np.ndarray:
    """
    Whether a frame (centre and half width and height on the sample, in µm) overlaps each of other frames (arrays
    of centres and half sizes), beyond a touching edge.
    """
    centers_um = np.asarray(centers_um, dtype=float).reshape(-1, 2)
    half_sizes_um = np.asarray(half_sizes_um, dtype=float).reshape(-1, 2)
    return np.all(
        np.abs(centers_um - np.asarray(center_um, dtype=float)) < (half_sizes_um + np.asarray(half_size_um)) * (1 - 1e-3),
        axis=1,
    )


def size_ranges_overlap(d_ranges_um, d_range_um) -> np.ndarray:
    """Whether diameter ranges (array of (d_min, d_max)) overlap a diameter range, beyond a common bound."""
    d_ranges_um = np.asarray(d_ranges_um, dtype=float).reshape(-1, 2)
    lo, hi = float(d_range_um[0]), float(d_range_um[1])
    return (d_ranges_um[:, 0] < hi * (1 - 1e-9)) & (lo < d_ranges_um[:, 1] * (1 - 1e-9))


class EdsSearchFrame(NamedTuple):
    """Frame searched by EDS, on the sample: centre and half width and height (µm), and size range of its scan."""

    frame_id: str
    center_um: Tuple[float, float]
    half_size_um: Tuple[float, float]
    d_range_um: Tuple[float, float]
    frame_width_um: float


def eds_search_frames(eds_scans: Sequence[ParticleSizeBand]) -> Tuple[List[EdsSearchFrame], np.ndarray]:
    """
    Frames searched by EDS with a known position, except those left out of the statistics as overlapping another
    (:func:`overlapping_frame_ids`), and the equivalent diameters (µm) of the particles counted in them.
    """
    excluded = set(overlapping_frame_ids(eds_scans))
    frames, diameters = [], []
    for scan in eds_scans:
        for f in scan.frames:
            if f.frame_id in excluded:
                continue
            diameters += [float(2.0 * np.sqrt(a / np.pi)) for a in f.particle_areas_um2 or []]
            if f.center_mm is not None:
                frames.append(EdsSearchFrame(
                    f.frame_id, (f.center_mm.x * 1000, f.center_mm.y * 1000),
                    (f.frame_width_um / 2, f.area_um2 / f.frame_width_um / 2), (scan.d_min_um, scan.d_max_um),
                    f.frame_width_um,
                ))
    return frames, np.asarray(diameters, dtype=float)


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
            same_sizes = size_ranges_overlap(d_ranges[k], d_ranges[i])
            overlap = frames_overlap(centers[i], half_sizes[i], centers[k], half_sizes[k])
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
    excluded_particles: Iterable[ParticleKey] = (),
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
    obs = size_observations(particles, bands, eds_scans, excluded_frame_ids, excluded_particles)
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


def bootstrap_confidence_intervals(
    particles: Sequence[ParticleInfo],
    bands: Sequence[ParticleSizeBand] = (),
    excluded_frame_ids: Iterable[str] = (),
    confidence: float = 0.95,
    n_resamples: int = N_BOOTSTRAP_RESAMPLES,
    rng: Optional[np.random.Generator] = None,
    eds_scans: Sequence[ParticleSizeBand] = (),
    excluded_particles: Iterable[ParticleKey] = (),
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
    excluded_particles : iterable of (frame_id, index), optional
        Particles of frames left out of the statistics.

    Returns
    -------
    dict
        Name (``D10``, ``D50``, ``D90``, ``D43``, ``Dv10``, ``Dv50``, ``Dv90``) to (low, high) bounds in µm; empty
        if a scan has fewer than :data:`MIN_BOOTSTRAP_UNITS` frames (or particles, without frame records).
    """
    rng = np.random.default_rng() if rng is None else rng
    excluded = _excluded_frames(list(bands) + list(eds_scans), excluded_frame_ids)
    obs = size_observations(particles, bands, eds_scans, excluded, excluded_particles)
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
) -> List[ParticleInfo]:
    """
    Particles of the ledger counted in the size statistics besides those of the frames: the particles measured by
    EDS with automatic particle selection (selection 'auto', with an area), if ``include_eds_auto_particles``, for
    samples without frames searched by EDS (measured before these were recorded), in ledger order.
    """
    if not include_eds_auto_particles:
        return []
    return [p for p in particles if p.area_um is not None and p.selection == "auto"]


def build_size_result(
    powder_meas_cfg: PowderMeasurementConfig,
    bands: Sequence[ParticleSizeBand] = (),
    eds_scans: Sequence[ParticleSizeBand] = (),
    particles: Sequence[ParticleInfo] = (),
    n_par_target: Optional[int] = None,
    measurement_id: Optional[int] = None,
    segmentation_id: Optional[int] = None,
    eds_segmentation_id: Optional[int] = None,
    excluded_frame_ids: Sequence[str] = (),
    excluded_particles: Sequence[ParticleKey] = (),
    missing_frame_ids: Sequence[str] = (),
    precision_reached: Optional[bool] = None,
) -> ParticleSizeStatsResult:
    """
    Particle size statistics result of frames under a segmentation (see :func:`size_stats_inputs`), weighted as in
    :func:`weighted_size_observations`, to store in the sample ledger (``SampleLedger.append_size_result``, which
    gives it its id).

    Parameters
    ----------
    powder_meas_cfg : PowderMeasurementConfig
        Settings of the measurement (of the EDS acquisition, without measurement).
    bands, eds_scans : sequence of ParticleSizeBand, optional
        Size bands of the measurement, and frames searched by EDS, with their frames under the segmentation.
    particles : sequence of ParticleInfo, optional
        Particles measured by EDS counted (samples without frames searched by EDS, see
        :func:`select_size_stats_particles`).
    n_par_target, measurement_id, segmentation_id, eds_segmentation_id, precision_reached : optional
        Recorded in the result.
    excluded_frame_ids, excluded_particles : sequence, optional
        Frames and particles (as (frame_id, index)) left out of the statistics, recorded in the result; the area of
        the frames left out is not counted.
    missing_frame_ids : sequence of str, optional
        Frames left out because they have no result of the segmentation, recorded in the result.

    Raises
    ------
    ValueError
        If no observation has a positive weight.
    """
    bands, eds_scans = list(bands), list(eds_scans)
    obs, weights, weighted_bands, weighted_eds_scans = weighted_size_observations(
        particles, bands, excluded_frame_ids, eds_scans, excluded_particles
    )
    n_par_by_selection = {}
    for key, n in (("size_stats", int(np.sum((obs.scan >= 0) & (obs.scan < obs.n_bands)))),
                   ("eds_frames", int(np.sum(obs.scan >= obs.n_bands))),
                   ("auto", int(np.sum(obs.scan == -1)))):
        if n:
            n_par_by_selection[key] = n

    def summaries(scans: Sequence[ParticleSizeBand]) -> List[BandSummary]:
        return [BandSummary(
            d_min_um=scan.d_min_um, d_max_um=scan.d_max_um, frame_width_um=scan.frame_width_um, n_frames=scan.n_frames,
            scanned_area_um2=scan.scanned_area_um2, n_particles=scan.n_particles,
            number_density_per_mm2=scan.number_density_per_mm2, n_eds_particles=scan.n_eds_particles,
            frame_ids=[frame.frame_id for frame in scan.frames],
        ) for scan in scans]

    return ParticleSizeStatsResult(
        result_id=0,
        measurement_id=measurement_id,
        segmentation_id=segmentation_id,
        eds_segmentation_id=eds_segmentation_id,
        n_par_target=n_par_target,
        powder_meas_cfg=powder_meas_cfg.model_copy(deep=True),
        n_par_by_selection=n_par_by_selection,
        size_stats=ParticleSizeStats(**compute_particle_size_stats(obs.eq_diameters_um, weights)),
        unweighted_size_stats=ParticleSizeStats(**compute_particle_size_stats(obs.eq_diameters_um)),
        n_par_without_weight=int(np.sum(weights == 0)),
        confidence_intervals_95=bootstrap_confidence_intervals(
            particles, bands, excluded_frame_ids, eds_scans=eds_scans, excluded_particles=excluded_particles
        ),
        precision_reached=precision_reached,
        bands=summaries(weighted_bands),
        eds_scans=summaries(weighted_eds_scans),
        excluded_frame_ids=list(excluded_frame_ids),
        excluded_particles=[(str(f), int(i)) for f, i in excluded_particles],
        overlapping_frame_ids=overlapping_frame_ids(bands + eds_scans),
        missing_frame_ids=list(missing_frame_ids),
    )


def working_frame(frame: FrameRecord, segmentation_id: int) -> Optional[ParticleStatsFrame]:
    """A frame of the ledger with its particles under a segmentation (None if not segmented with it)."""
    result = frame.segmentation(segmentation_id)
    if result is None:
        return None
    marks = [(int(x), int(y), int(r), i) for i, (x, y, r) in enumerate(result.particle_marks_px)]
    return ParticleStatsFrame(
        frame_id=frame.frame_id, center_mm=frame.center_mm, frame_width_um=frame.frame_width_um,
        pixel_size_um=frame.pixel_size_um, area_um2=frame.area_um2, n_particles=result.n_particles,
        particle_areas_um2=list(result.particle_areas_um2),
        particle_edge_corrections=list(result.particle_edge_corrections),
        image_file=frame.image_file, particle_marks_px=marks,
    )


class SizeStatsInputs(NamedTuple):
    """What the size statistics of a sample are computed from (see :func:`size_stats_inputs`)."""

    bands: List[ParticleSizeBand]
    eds_scans: List[ParticleSizeBand]
    particles: List[ParticleInfo]
    missing_frame_ids: List[str]


def size_stats_inputs(
    ledger: SampleLedger,
    measurement_id: Optional[int] = None,
    segmentation_id: Optional[int] = None,
    eds_segmentation_id: Optional[int] = None,
    include_eds_auto_particles: bool = True,
) -> SizeStatsInputs:
    """
    Size bands of measurement ``measurement_id`` with their frames segmented with ``segmentation_id``, and the frames
    searched by EDS segmented with ``eds_segmentation_id`` (one scan per size range and frame width), from the
    ledger frames; frames without result of the segmentation are left out (``missing_frame_ids``). Without frames
    searched by EDS (samples measured before these were recorded), the particles measured by EDS with automatic
    selection count, if ``include_eds_auto_particles``.
    """
    bands: List[ParticleSizeBand] = []
    missing: List[str] = []
    measurement = ledger.measurement(measurement_id)
    if measurement is not None:
        for band in measurement.bands:
            frames = []
            for frame_id in band.frame_ids:
                frame = ledger.frame(frame_id)
                working = working_frame(frame, segmentation_id) if frame is not None else None
                if working is None:
                    missing.append(frame_id)
                else:
                    frames.append((frame, working))
            # The band starts at its smallest particles with enough pixels under the segmentation
            d_min = max([band.d_min_um] + [f.segmentation(segmentation_id).d_range_um[0] for f, _ in frames])
            working_band = ParticleSizeBand(d_min_um=min(d_min, band.d_max_um), d_max_um=band.d_max_um,
                                            frame_width_um=band.frame_width_um, n_frames=0, scanned_area_um2=0.0)
            for _, working in frames:
                working_band.add_frame(working)
            bands.append(working_band)

    eds_scans: List[ParticleSizeBand] = []
    eds_frames = ledger.frames_of("eds_search")
    if eds_segmentation_id is not None:
        for frame in eds_frames:
            working = working_frame(frame, eds_segmentation_id)
            if working is None:
                missing.append(frame.frame_id)
                continue
            d_lo, d_hi = frame.segmentation(eds_segmentation_id).d_range_um
            scan = next((s for s in eds_scans if np.allclose((s.d_min_um, s.d_max_um, s.frame_width_um),
                                                             (d_lo, d_hi, frame.frame_width_um))), None)
            if scan is None:
                scan = ParticleSizeBand(d_min_um=d_lo, d_max_um=d_hi, frame_width_um=frame.frame_width_um,
                                        n_frames=0, scanned_area_um2=0.0)
                eds_scans.append(scan)
            scan.add_frame(working)
    particles = [] if eds_frames else select_size_stats_particles(ledger.particles, include_eds_auto_particles)
    return SizeStatsInputs(bands, eds_scans, particles, missing)


def eds_segmentation_for(ledger: SampleLedger, segmentation_id: Optional[int] = None) -> Optional[int]:
    """Segmentation of the frames searched by EDS counted with the frames of segmentation ``segmentation_id``: the
    same if they have its results, else the last one they have (None without frames searched by EDS)."""
    eds_frames = ledger.frames_of("eds_search")
    if not eds_frames:
        return None
    if segmentation_id is not None and any(f.segmentation(segmentation_id) for f in eds_frames):
        return segmentation_id
    return max((s.segmentation_id for f in eds_frames for s in f.segmentations), default=None)


def default_result_inputs(ledger: SampleLedger, include_eds_auto_particles: Optional[bool] = None) -> dict:
    """
    What the size statistics of a sample are computed from by default (keyword arguments of
    :func:`recalculate_particle_stats`): those of the active result; else the last completed measurement under the
    segmentation it was acquired with, and the frames searched by EDS (under the same segmentation if they have its
    results, else under the last one they have) if ``include_eds_auto_particles`` (default: the setting of the
    measurement, or of the EDS acquisition).
    """
    active = ledger.size_result
    if active is not None:
        include = (active.powder_meas_cfg.include_eds_auto_particles if include_eds_auto_particles is None
                   else include_eds_auto_particles)
        eds_segmentation = active.eds_segmentation_id
        if eds_segmentation is None and include:  # frames searched by EDS since (e.g. EDS after the measurement)
            eds_segmentation = eds_segmentation_for(ledger, active.segmentation_id)
        return {"measurement_id": active.measurement_id, "segmentation_id": active.segmentation_id,
                "eds_segmentation_id": eds_segmentation,
                "particles_to_filter": active.excluded_particles, "frame_ids_to_filter": active.excluded_frame_ids,
                "include_eds_auto_particles": include}
    measurement = next((m for m in reversed(ledger.size_measurements) if m.status == "completed"), None)
    cfg = (measurement.powder_meas_cfg if measurement is not None
           else ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig())
    include = cfg.include_eds_auto_particles if include_eds_auto_particles is None else include_eds_auto_particles
    segmentation = measurement.segmentation_id if measurement is not None else None
    eds_segmentation = eds_segmentation_for(ledger, segmentation) if include else None
    return {"measurement_id": measurement.measurement_id if measurement is not None else None,
            "segmentation_id": segmentation, "eds_segmentation_id": eds_segmentation,
            "particles_to_filter": [], "frame_ids_to_filter": [], "include_eds_auto_particles": include}


_UNSET = object()


def recalculate_particle_stats(
    ledger: SampleLedger,
    particles_to_filter: Optional[Iterable[ParticleKey]] = None,
    frame_ids_to_filter: Optional[Iterable[str]] = None,
    include_eds_auto_particles: Optional[bool] = None,
    measurement_id: Any = _UNSET,
    segmentation_id: Any = _UNSET,
    eds_segmentation_id: Any = _UNSET,
) -> Optional[ParticleSizeStatsResult]:
    """
    Particle size statistics of a sample recalculated from its ledger frames, leaving out some particles or frames
    (not stored: see ``SampleLedger.append_size_result``).

    By default (see :func:`default_result_inputs`), from what the active result was computed from (the last
    completed measurement and the frames searched by EDS, without result): the measurement, its segmentation and
    the segmentation of the frames searched by EDS can be given instead (``eds_segmentation_id=None``: frames
    searched by EDS not counted). A sample measured only by EDS gets statistics from its EDS frames (or, measured
    before these were recorded, from its particles measured by EDS with automatic selection).

    Parameters
    ----------
    ledger : SampleLedger
        Ledger of the sample (not modified).
    particles_to_filter : iterable of (frame_id, index), optional
        Particles to leave out. Default (None): those left out of the active result, in the frames segmented by the
        same segmentation (particles are numbered per segmentation).
    frame_ids_to_filter : iterable of str, optional
        IDs of the frames whose particles are left out (exact matches). Default (None): those left out of the active
        result.
    include_eds_auto_particles : bool, optional
        Whether to count what EDS measured. Default: the setting of the active result (of the measurement, or of
        the EDS acquisition).
    measurement_id, segmentation_id, eds_segmentation_id : int or None, optional
        What to compute the statistics from (default: see above).

    Returns
    -------
    ParticleSizeStatsResult or None
        The statistics, recording the particles and frames left out (those found), or None if nothing is left.
    """
    defaults = default_result_inputs(ledger, include_eds_auto_particles)
    include = defaults["include_eds_auto_particles"]
    if frame_ids_to_filter is None:
        frame_ids_to_filter = defaults["frame_ids_to_filter"]
    measurement_id = defaults["measurement_id"] if measurement_id is _UNSET else measurement_id
    segmentation_id = defaults["segmentation_id"] if segmentation_id is _UNSET else segmentation_id
    if eds_segmentation_id is _UNSET:
        eds_segmentation_id = defaults["eds_segmentation_id"] if include else None
    if particles_to_filter is None:
        # Particles are numbered per segmentation: those left out of the active result only in frames segmented alike
        def same_segmentation(frame_id: str) -> bool:
            frame = ledger.frame(frame_id)
            if frame is None:
                return False
            if frame.purpose == "eds_search":
                return eds_segmentation_id == defaults["eds_segmentation_id"]
            return segmentation_id == defaults["segmentation_id"]

        particles_to_filter = [key for key in defaults["particles_to_filter"] if same_segmentation(key[0])]
    inputs = size_stats_inputs(ledger, measurement_id, segmentation_id, eds_segmentation_id, include)

    frames = {f.frame_id: f for scan in inputs.bands + inputs.eds_scans for f in scan.frames}
    frame_ids_to_filter = [str(fid) for fid in (frame_ids_to_filter or [])]
    # Frames of the scans, or of particles measured by EDS (samples without frames searched by EDS)
    legacy_frames = {p.frame_id for p in inputs.particles}
    excluded_frame_ids = [fid for fid in frame_ids_to_filter if fid in frames or fid in legacy_frames]
    excluded_particles = [(str(f), int(i)) for f, i in (particles_to_filter or [])
                          if str(f) in frames and 0 <= int(i) < frames[str(f)].n_particles]
    obs, weights, _, _ = weighted_size_observations(
        inputs.particles, inputs.bands, excluded_frame_ids, inputs.eds_scans, excluded_particles
    )
    if not np.any(weights > 0):
        logger.warning("No particle left in the scans of the sample: no statistics.")
        return None
    measurement = ledger.measurement(measurement_id)
    cfg = (measurement.powder_meas_cfg if measurement is not None
           else ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig())
    return build_size_result(
        cfg.model_copy(update={"include_eds_auto_particles": include}), inputs.bands, inputs.eds_scans,
        inputs.particles, n_par_target=measurement.n_par_target if measurement is not None else None,
        measurement_id=measurement_id if measurement is not None else None,
        segmentation_id=segmentation_id if measurement is not None else None,
        eds_segmentation_id=eds_segmentation_id if inputs.eds_scans else None,
        excluded_frame_ids=excluded_frame_ids, excluded_particles=excluded_particles,
        missing_frame_ids=inputs.missing_frame_ids,
        precision_reached=measurement.precision_reached if measurement is not None else None,
    )


def size_distribution(
    ledger: SampleLedger, result: Optional[ParticleSizeStatsResult] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Equivalent diameters (µm) and weights of the particle sizes behind a statistics result of a sample: what
    :func:`save_particle_size_histogram` plots.

    Parameters
    ----------
    ledger : SampleLedger
        Ledger of the sample.
    result : ParticleSizeStatsResult, optional
        Statistics of the sample (default: the active one, else those computed by default).

    Returns
    -------
    eq_diameters_um, weights : numpy.ndarray
        Diameters with a positive weight, and their weights.
    """
    result = result if result is not None else ledger.size_result
    if result is None:
        result = recalculate_particle_stats(ledger)
        if result is None:
            return np.array([]), np.array([])
    include = result.powder_meas_cfg.include_eds_auto_particles
    inputs = size_stats_inputs(ledger, result.measurement_id, result.segmentation_id, result.eds_segmentation_id,
                               include)
    obs, weights, _, _ = weighted_size_observations(
        inputs.particles, inputs.bands, result.excluded_frame_ids, inputs.eds_scans, result.excluded_particles
    )
    keep = weights > 0
    return obs.eq_diameters_um[keep], weights[keep]


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
    excluded_particles: Iterable[ParticleKey] = (),
) -> pd.DataFrame:
    """
    Compute the particle size statistics and save the particle sizes, the statistics and the histogram.

    With size bands, the statistics are weighted by band (see :func:`band_particle_weights`): the particle sizes
    file then gives the weight of each particle (relative, 1 on average over the particles counted), and the
    statistics file the number of bands.

    Parameters
    ----------
    particles : sequence of ParticleInfo
        Particles measured by EDS included in the statistics (each with ``area_um``), for samples without frames
        searched by EDS.
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
        Frames searched by EDS included in the statistics.
    excluded_particles : iterable of (frame_id, index), optional
        Particles of frames left out of the statistics.

    Each particle of a frame is listed with its frame and its index in the frame; particles measured by EDS
    (``particles``) with their ID.

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
    excluded_particles = list(excluded_particles)
    obs, weights, _, _ = weighted_size_observations(
        particles, list(bands or []), excluded_frame_ids, list(eds_scans), excluded_particles
    )
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

    # Particles of frames (index in their frame), and particles measured by EDS (their ID)
    particle_data = pd.DataFrame({
        cnst.PAR_ID_DF_KEY: pd.array([particles[k].id if k is not None else key[1]
                                      for k, key in zip(obs.particle_index, obs.particle_keys)], dtype="Int64"),
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
            particles, list(bands or []), excluded_frame_ids, eds_scans=list(eds_scans),
            excluded_particles=excluded_particles,
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


def size_result_label(ledger, result) -> str:
    """Readable name of a statistics result: id, measurement, segmentations, particles, frames and particles left out,
    e.g. 'Result 2 · measurement 0 · SAM2 · 80 px · 412 particles · active'."""
    parts = [f"Result {result.result_id}"]
    if result.measurement_id is not None:
        parts.append(f"measurement {result.measurement_id}")
    segmentation = ledger.segmentation(result.segmentation_id)
    if segmentation is not None:
        parts.append(segmentation.label)
    if result.eds_segmentation_id is not None:
        eds = ledger.segmentation(result.eds_segmentation_id)
        if result.eds_segmentation_id == result.segmentation_id or eds is None:
            parts.append("+ EDS frames")
        else:
            parts.append(f"{'+ ' if segmentation is not None else ''}EDS frames {eds.label}")
    parts.append(f"{result.size_stats.n_par_analysed} particles")
    n_out = len(result.excluded_frame_ids) + len(result.excluded_particles)
    if n_out:
        parts.append(f"{len(result.excluded_frame_ids)} frame(s), {len(result.excluded_particles)} particle(s) left out")
    if result.result_id == ledger.active_size_result:
        parts.append("active")
    return " · ".join(parts)
