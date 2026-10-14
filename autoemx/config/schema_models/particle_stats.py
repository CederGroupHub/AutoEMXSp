#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Particle size statistics stored in the sample ledger (``SampleLedger.particle_stats``)."""

from datetime import datetime
from typing import Dict, List, Literal, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from autoemx._version import __version__ as AUTOEMX_VERSION
from autoemx.config.runtime_configs import PowderMeasurementConfig

from .acquisition import Coordinate2D

#: How a particle was chosen: found by a particle size statistics measurement ('size_stats'), or measured by EDS
#: with particle selection mode 'auto', 'manual' or 'list' (PowderMeasurementConfig.par_selection_mode).
ParticleSelection = Literal["size_stats", "auto", "manual", "list"]


class ParticleSizeStats(BaseModel):
    """
    Summary statistics of particle equivalent diameters, in µm (``median`` is the D50).

    ``D43`` (volume-weighted mean, D[4,3]) and the volume-weighted percentiles ``Dv10``, ``Dv50``, ``Dv90`` weight
    each particle by the volume of a sphere of its equivalent diameter, for comparison with laser diffraction
    (not computed before AutoEMX 1.6.0).
    """

    n_par_analysed: int
    mean: float
    stdev: float  # population standard deviation
    median: float
    max: float
    min: float
    D10: float
    D25: float
    D75: float
    D90: float
    D43: Optional[float] = None
    Dv10: Optional[float] = None
    Dv50: Optional[float] = None
    Dv90: Optional[float] = None

    model_config = ConfigDict(extra="forbid")

    @field_validator("n_par_analysed")
    @classmethod
    def validate_n_par(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("n_par_analysed must be positive")
        return v

    @field_validator("mean", "stdev", "median", "max", "min", "D10", "D25", "D75", "D90", "D43", "Dv10", "Dv50", "Dv90")
    @classmethod
    def validate_finite(cls, v: Optional[float]) -> Optional[float]:
        if v is not None and not np.isfinite(v):
            raise ValueError("Particle size statistics must be finite")
        return v


class ParticleStatsFrame(BaseModel):
    """
    One frame analysed by a particle size statistics measurement: ``center_mm`` (stage position of its centre, if
    known), frame width and pixel size, ``area_um2`` (total area of the frame) and ``n_particles`` (particles found).

    Frames searched by EDS (``SampleLedger.eds_size_scans``) also give the area (``particle_areas_um2``) and
    frame-edge correction (``particle_edge_corrections``) of each particle valid for the statistics found in them,
    whether measured by EDS or not.

    ``image_file`` is the frame image (in the ``SEM images`` folder of the sample), saved without annotations
    unless ``MeasurementConfig.annotate_particle_images``; ``particle_marks_px`` are the particles to mark on it
    (x, y, radius in pixels, and particle id or None), drawn by ``autoemx.runners.annotate_particle_images``.
    """

    frame_id: str
    center_mm: Optional[Coordinate2D] = None
    frame_width_um: float
    pixel_size_um: float
    area_um2: float
    n_particles: int
    particle_areas_um2: Optional[List[float]] = None
    particle_edge_corrections: Optional[List[float]] = None
    image_file: Optional[str] = None
    particle_marks_px: Optional[List[Tuple[int, int, int, Optional[int]]]] = None

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def validate_particle_lists(self) -> "ParticleStatsFrame":
        lists = (self.particle_areas_um2, self.particle_edge_corrections)
        if any(values is not None for values in lists):
            if any(values is None or len(values) != self.n_particles for values in lists):
                raise ValueError("particle_areas_um2 and particle_edge_corrections must both give n_particles values")
            if any(not np.isfinite(a) or a <= 0 for a in self.particle_areas_um2):
                raise ValueError("particle_areas_um2 must be finite and positive")
            if any(not np.isfinite(c) or c < 1 for c in self.particle_edge_corrections):
                raise ValueError("particle_edge_corrections must be finite and at least 1")
        return self

    @field_validator("frame_width_um", "pixel_size_um", "area_um2")
    @classmethod
    def validate_positive(cls, v: float) -> float:
        if not np.isfinite(v) or v <= 0:
            raise ValueError("Frame width, pixel size and area must be finite and positive")
        return v

    @field_validator("n_particles")
    @classmethod
    def validate_n_particles(cls, v: int) -> int:
        if v < 0:
            raise ValueError("n_particles must be non-negative")
        return v


class ParticleSizeBand(BaseModel):
    """
    One size band of a particle size statistics measurement: frames searched for particles of equivalent diameter
    between ``d_min_um`` and ``d_max_um``.

    Measurements are split into size bands, each searched with a frame width matched to its largest particles and
    going down to the smallest particles covering ``PowderMeasurementConfig.min_par_pixels`` pixels at that frame
    width (see ``EM_Particle_Finder.get_particle_stats``). ``scanned_area_um2`` is the total area
    of the frames analysed in the band (empty frames included), ``frame_width_um`` the frame width set by the
    automatic navigation (None with manual navigation), ``particle_ids`` the particles found in the band (including
    any later left out of the statistics, see ``ParticleStatsResult.excluded_particle_ids``).

    ``frames`` are the frames analysed in the band (not recorded before AutoEMX 1.6.0); ``n_frames`` and
    ``scanned_area_um2`` are their number and total area. The same model describes the frames searched by EDS
    (``SampleLedger.eds_size_scans``), whose particles are listed in each frame instead of ``particle_ids``.

    ``number_density_per_mm2`` (particles found in the band, left-out ones excepted and corrected for those
    touching the frame edge, per mm² scanned) weights the
    band in the statistics, and ``n_eds_particles`` is the number of particles measured by EDS whose diameter falls
    in the band, which share its weight (see :func:`autoemx.core.particle_statistics.band_particle_weights`).
    """

    d_min_um: float
    d_max_um: float
    frame_width_um: Optional[float] = None
    n_frames: int
    scanned_area_um2: float
    particle_ids: List[int] = Field(default_factory=list)
    frames: List[ParticleStatsFrame] = Field(default_factory=list)
    number_density_per_mm2: Optional[float] = None
    n_eds_particles: int = 0

    model_config = ConfigDict(extra="forbid")

    @field_validator("d_min_um", "d_max_um", "scanned_area_um2")
    @classmethod
    def validate_non_negative(cls, v: float) -> float:
        if not np.isfinite(v) or v < 0:
            raise ValueError("Band diameters and scanned area must be finite and non-negative")
        return v

    @field_validator("frame_width_um")
    @classmethod
    def validate_frame_width(cls, v: Optional[float]) -> Optional[float]:
        if v is not None and (not np.isfinite(v) or v <= 0):
            raise ValueError("frame_width_um must be finite and positive when set")
        return v

    @field_validator("n_frames", "n_eds_particles")
    @classmethod
    def validate_counts(cls, v: int) -> int:
        if v < 0:
            raise ValueError("n_frames and n_eds_particles must be non-negative")
        return v

    @field_validator("number_density_per_mm2")
    @classmethod
    def validate_density(cls, v: Optional[float]) -> Optional[float]:
        if v is not None and (not np.isfinite(v) or v < 0):
            raise ValueError("number_density_per_mm2 must be finite and non-negative when set")
        return v

    @model_validator(mode="after")
    def validate_frames(self) -> "ParticleSizeBand":
        if self.frames and self.n_frames != len(self.frames):
            raise ValueError("n_frames must be the number of frames recorded")
        return self

    def add_frame(self, frame: ParticleStatsFrame, particle_ids: List[int]) -> None:
        """Record a frame analysed in the band, with the ids of the particles found in it."""
        self.frames.append(frame)
        self.particle_ids.extend(particle_ids)
        self.n_frames = len(self.frames)
        self.scanned_area_um2 += frame.area_um2

    def frame_area_um2(self, frame_ids) -> float:
        """Total area of the recorded frames of the band with the given ids."""
        frame_ids = set(frame_ids)
        return float(sum(frame.area_um2 for frame in self.frames if frame.frame_id in frame_ids))


class ParticleStatsResult(BaseModel):
    """
    Result of a particle size statistics measurement.

    ``particle_ids`` are the ledger particles included in the statistics: the particles found by the measurement
    (selection 'size_stats') and, if ``powder_meas_cfg.include_eds_auto_particles``, the particles measured by EDS
    with automatic particle selection (selection 'auto'). ``powder_meas_cfg`` holds the settings at the start of
    the measurement. ``bands`` are the size bands searched by the measurement, largest particles first.
    ``eds_scans`` are the frames searched by EDS on the same sample (``SampleLedger.eds_size_scans``) included in
    the statistics: they then replace the particles measured by EDS, which they count.

    ``size_stats`` are number-weighted by size band: each band counts in proportion to the number of particles per
    unit area scanned in it, whatever the number of particles collected in each band, and each particle in
    proportion to its edge correction (``ParticleInfo.edge_correction``), since particles touching the frame edge
    are left out. ``unweighted_size_stats`` treat all particles equally, for comparison.
    ``n_par_without_weight`` particles are left out of ``size_stats``: particles measured by EDS outside every band,
    or in a band where the measurement found no particle.

    ``confidence_intervals_95`` are the 95 % confidence intervals of D10, D50, D90, D43, Dv10, Dv50 and Dv90 (by
    bootstrap of the frames of each band, see
    :func:`autoemx.core.particle_statistics.bootstrap_confidence_intervals`). ``precision_reached`` tells whether
    a measurement stopping on precision (``powder_meas_cfg.size_stats_ci_tolerance``) reached it (None otherwise,
    and for recalculated statistics). ``excluded_particle_ids`` and ``excluded_frame_ids`` are the particles and
    frames left out when the statistics were recalculated (``autoemx.runners.process_particle_statistics``), and
    ``overlapping_frame_ids`` the frames left out because they overlap a larger frame counting the same particle
    sizes (see :func:`autoemx.core.particle_statistics.overlapping_frame_ids`); the bands keep the frames,
    scanned area and particles of the measurement.
    """

    autoemx_version: str = AUTOEMX_VERSION
    datetime: str = Field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    n_par_target: Optional[int] = None
    powder_meas_cfg: PowderMeasurementConfig
    particle_ids: List[int]
    n_par_by_selection: Dict[str, int] = Field(default_factory=dict)
    size_stats: ParticleSizeStats
    unweighted_size_stats: Optional[ParticleSizeStats] = None
    n_par_without_weight: int = 0
    confidence_intervals_95: Dict[str, Tuple[float, float]] = Field(default_factory=dict)
    precision_reached: Optional[bool] = None
    bands: List[ParticleSizeBand] = Field(default_factory=list)
    eds_scans: List[ParticleSizeBand] = Field(default_factory=list)
    excluded_particle_ids: List[int] = Field(default_factory=list)
    excluded_frame_ids: List[str] = Field(default_factory=list)
    overlapping_frame_ids: List[str] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

    @field_validator("particle_ids")
    @classmethod
    def validate_particle_ids(cls, v: List[int]) -> List[int]:
        # Empty for a sample measured only by EDS (its statistics come from the frames searched by EDS)
        if len(v) != len(set(v)):
            raise ValueError("particle_ids must be unique")
        return v

    @model_validator(mode="after")
    def validate_band_particles(self) -> "ParticleStatsResult":
        band_ids = [particle_id for band in self.bands for particle_id in band.particle_ids]
        if len(band_ids) != len(set(band_ids)):
            raise ValueError("A particle can belong to one size band only")
        return self


class ParticleStatsProgress(BaseModel):
    """
    Particle size statistics measurement in progress (``SampleLedger.particle_stats_progress``): its target,
    settings and size bands so far (the last one possibly unfinished), saved regularly so that an interrupted
    measurement can be resumed. Replaced by ``SampleLedger.particle_stats`` when the measurement ends.
    """

    autoemx_version: str = AUTOEMX_VERSION
    datetime: str = Field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    n_par_target: int
    powder_meas_cfg: PowderMeasurementConfig
    bands: List[ParticleSizeBand] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")
