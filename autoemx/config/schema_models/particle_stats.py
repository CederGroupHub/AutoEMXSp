#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle size statistics stored in the sample ledger: the measurements (``SampleLedger.size_measurements``) and the
statistics computed from frames under a segmentation (``SampleLedger.size_results``), with the working records of
size bands and frames used to compute them (built from ``SampleLedger.frames``, not stored).
"""

from datetime import datetime
from typing import Dict, List, Literal, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from autoemx._version import __version__ as AUTOEMX_VERSION
from autoemx.config.runtime_configs import PowderMeasurementConfig

from .acquisition import Coordinate2D

#: How a particle was chosen: found by a particle size statistics measurement ('size_stats', ledgers before frames
#: were recorded), or measured by EDS with particle selection mode 'auto', 'manual' or 'list'
#: (PowderMeasurementConfig.par_selection_mode).
ParticleSelection = Literal["size_stats", "auto", "manual", "list"]

#: A particle found by a segmentation of a frame: (frame_id, index of the particle in the frame result)
ParticleKey = Tuple[str, int]


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
    Working record of a frame under one segmentation (built from ``SampleLedger.frames``): ``center_mm`` (stage
    position of its centre, if known), frame width and pixel size, ``area_um2`` (total area of the frame), the area
    and frame-edge correction of each particle valid for the statistics, ``image_file`` (in the ``SEM images`` folder
    of the sample) and ``particle_marks_px`` (x, y, radius in pixels, and index of the particle in the frame).
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
    Working record of the frames of a size band of a measurement, or of frames searched by EDS with the same size
    range and frame width, under one segmentation: diameters searched, frame width, frames and scanned area, and
    (once weighted) ``number_density_per_mm2`` (particles counted per mm²) and ``n_eds_particles`` (particles measured
    by EDS sharing the weight of the band, ledgers without frames searched by EDS).
    """

    d_min_um: float
    d_max_um: float
    frame_width_um: Optional[float] = None
    n_frames: int
    scanned_area_um2: float
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

    @model_validator(mode="after")
    def validate_frames(self) -> "ParticleSizeBand":
        if self.frames and self.n_frames != len(self.frames):
            raise ValueError("n_frames must be the number of frames recorded")
        return self

    @property
    def n_particles(self) -> int:
        """Particles counted in the frames of the band."""
        return sum(frame.n_particles for frame in self.frames)

    def add_frame(self, frame: ParticleStatsFrame) -> None:
        """Record a frame analysed in the band."""
        self.frames.append(frame)
        self.n_frames = len(self.frames)
        self.scanned_area_um2 += frame.area_um2

    def frame_area_um2(self, frame_ids) -> float:
        """Total area of the recorded frames of the band with the given ids."""
        frame_ids = set(frame_ids)
        return float(sum(frame.area_um2 for frame in self.frames if frame.frame_id in frame_ids))


class SizeBand(BaseModel):
    """Size band of a particle size statistics measurement: diameters searched, frame width, and its frames."""

    d_min_um: float
    d_max_um: float
    frame_width_um: Optional[float] = None
    frame_ids: List[str] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")


class SizeMeasurement(BaseModel):
    """
    A particle size statistics measurement of the sample: its target, settings, size bands (largest particles
    first, with the frames acquired in each) and the segmentation used while acquiring (which decided when each
    band reached its target). ``status`` is 'in_progress' while acquiring (an interrupted measurement with the same
    target and settings is resumed), then 'completed', or 'interrupted' if a new measurement started instead of
    resuming it (its frames are kept). ``count_eds_frames``: it continued an EDS acquisition, counting the frames
    searched by EDS.
    """

    measurement_id: int
    autoemx_version: str = AUTOEMX_VERSION
    datetime: str = Field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    status: Literal["in_progress", "completed", "interrupted"] = "in_progress"
    n_par_target: int
    powder_meas_cfg: PowderMeasurementConfig
    segmentation_id: int
    count_eds_frames: bool = False
    bands: List[SizeBand] = Field(default_factory=list)
    precision_reached: Optional[bool] = None

    model_config = ConfigDict(extra="forbid")

    @property
    def frame_ids(self) -> List[str]:
        return [frame_id for band in self.bands for frame_id in band.frame_ids]


class BandSummary(BaseModel):
    """Size band (or frames searched by EDS) as counted in a statistics result."""

    d_min_um: float
    d_max_um: float
    frame_width_um: Optional[float] = None
    n_frames: int
    scanned_area_um2: float
    n_particles: int
    number_density_per_mm2: Optional[float] = None
    n_eds_particles: int = 0
    frame_ids: List[str] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")


class ParticleSizeStatsResult(BaseModel):
    """
    Particle size statistics of the sample computed from frames under a segmentation: the frames of measurement
    ``measurement_id`` segmented with ``segmentation_id`` and, if ``eds_segmentation_id`` is set, the frames searched
    by EDS segmented with it (usually the same segmentation); frames without result of the segmentation are left out.
    For samples measured before frames searched by EDS were recorded, the particles measured by EDS with automatic
    selection (``n_par_by_selection['auto']``) share the weight of the band of their diameter.

    ``size_stats`` are number-weighted: each particle counts its frame-edge correction over the area of the frames
    able to observe its size, so that each size counts in proportion to its number of particles per unit area of
    the sample; ``unweighted_size_stats`` treat all particles equally. ``confidence_intervals_95`` are bootstrapped
    by frame. Frames and particles left out (``excluded_frame_ids``, ``excluded_particles`` as (frame_id, index))
    are recorded; ``overlapping_frame_ids`` are left out because a larger frame counts the same sizes there.
    ``powder_meas_cfg`` are the settings of the measurement (of the EDS acquisition, without measurement).
    """

    result_id: int
    autoemx_version: str = AUTOEMX_VERSION
    datetime: str = Field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    measurement_id: Optional[int] = None
    segmentation_id: Optional[int] = None
    eds_segmentation_id: Optional[int] = None
    n_par_target: Optional[int] = None
    powder_meas_cfg: PowderMeasurementConfig
    n_par_by_selection: Dict[str, int] = Field(default_factory=dict)
    size_stats: ParticleSizeStats
    unweighted_size_stats: Optional[ParticleSizeStats] = None
    n_par_without_weight: int = 0
    confidence_intervals_95: Dict[str, Tuple[float, float]] = Field(default_factory=dict)
    precision_reached: Optional[bool] = None
    bands: List[BandSummary] = Field(default_factory=list)
    eds_scans: List[BandSummary] = Field(default_factory=list)
    excluded_frame_ids: List[str] = Field(default_factory=list)
    excluded_particles: List[ParticleKey] = Field(default_factory=list)
    overlapping_frame_ids: List[str] = Field(default_factory=list)
    missing_frame_ids: List[str] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

    def inputs_key(self) -> tuple:
        """What the result was computed from: a new result with the same inputs replaces it."""
        return (self.measurement_id, self.segmentation_id, self.eds_segmentation_id,
                self.powder_meas_cfg.include_eds_auto_particles, tuple(sorted(self.excluded_frame_ids)),
                tuple(sorted(map(tuple, self.excluded_particles))))
