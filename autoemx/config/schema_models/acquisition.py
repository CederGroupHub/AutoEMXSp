#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Any, List, Optional, Tuple
from pathlib import Path, PurePosixPath, PureWindowsPath

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .quantification import QuantificationResult


class Coordinate2D(BaseModel):
    """Simple 2D coordinate pair."""

    x: float
    y: float

    model_config = ConfigDict(extra="forbid")

    @field_validator("x", "y")
    @classmethod
    def validate_coordinate(cls, v: float) -> float:
        if not np.isfinite(v):
            raise ValueError("Coordinates must be finite")
        return v


class SpotCoordinates(BaseModel):
    """Acquisition spot coordinates within a frame image.

    ``machine_coordinates`` are normalized frame-relative coordinates (as used by
    the EM API), not absolute stage positions in mm. They are derived from
    ``pixel_coordinates`` via ``frame_pixel_to_rel_coords``.
    """

    machine_coordinates: Optional[Coordinate2D] = None
    pixel_coordinates: Optional[Tuple[int, int]] = None

    model_config = ConfigDict(extra="forbid")

    @field_validator("pixel_coordinates")
    @classmethod
    def validate_pixel_coordinates(cls, v: Optional[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        if v is None:
            return None
        if len(v) != 2:
            raise ValueError("pixel_coordinates must contain exactly two integers")
        return (int(v[0]), int(v[1]))


class AcquisitionDetails(BaseModel):
    """Structured acquisition context for one spectrum."""

    frame_id: Optional[str] = None
    particle_id: Optional[int] = None
    spot_coordinates: Optional[SpotCoordinates] = None

    model_config = ConfigDict(extra="forbid")

    @field_validator("frame_id")
    @classmethod
    def normalize_frame_id(cls, v: Optional[str]) -> Optional[str]:
        if v is None:
            return None
        normalized = v.strip()
        return normalized if normalized else None

    @field_validator("particle_id")
    @classmethod
    def validate_particle_id(cls, v: Optional[int]) -> Optional[int]:
        if v is None:
            return None
        if v < 0:
            raise ValueError("particle_id must be non-negative")
        return v


class ParticleInfo(BaseModel):
    """Geometry and analysis summary for one segmented particle."""

    id: int
    area_um: Optional[float] = None
    eq_diameter_um: Optional[float] = None
    clusters: List[int] = Field(default_factory=list)
    composition: Optional[str] = None

    model_config = ConfigDict(extra="forbid")

    @field_validator("id")
    @classmethod
    def validate_id(cls, v: int) -> int:
        if v < 0:
            raise ValueError("id must be non-negative")
        return v

    @field_validator("area_um", "eq_diameter_um")
    @classmethod
    def validate_positive_finite(cls, v: Optional[float]) -> Optional[float]:
        if v is None:
            return None
        if not np.isfinite(v) or v <= 0:
            raise ValueError("area_um and eq_diameter_um must be finite and positive when set")
        return v


class SpectrumEntry(BaseModel):
    """Minimal per-spectrum record inside a sample ledger."""

    spectrum_id: Optional[str] = None
    total_counts: Optional[int] = None
    live_acquisition_time: float
    acquisition_details: Optional[AcquisitionDetails] = None
    spectrum_relpath: Optional[str] = None
    instrument_background_relpath: Optional[str] = None
    quantification_results: List[QuantificationResult] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def drop_legacy_fields(cls, data: Any) -> Any:
        """Drop removed legacy fields so old ledger JSON remains loadable."""
        if isinstance(data, dict):
            cleaned = dict(data)
            cleaned.pop("live_time", None)
            cleaned.pop("background_vals", None)
            cleaned.pop("spectrum_pointer_exists", None)
            cleaned.pop("original_spectrum_exists", None)
            cleaned.pop("spectrum_source", None)
            cleaned.pop("metadata", None)
            return cleaned
        return data

    @field_validator("live_acquisition_time")
    @classmethod
    def validate_live_acquisition_time(cls, v: float) -> float:
        if not np.isfinite(v) or v <= 0:
            raise ValueError("live_acquisition_time must be finite and positive")
        return v

    @field_validator("total_counts")
    @classmethod
    def validate_total_counts(cls, v: Optional[int]) -> Optional[int]:
        if v is not None and v < 0:
            raise ValueError("total_counts cannot be negative")
        return v

    @field_validator("spectrum_relpath")
    @classmethod
    def validate_spectrum_relpath(cls, v: Optional[str]) -> Optional[str]:
        return cls._normalize_relpath(v, field_name="spectrum_relpath")

    @field_validator("instrument_background_relpath")
    @classmethod
    def validate_instrument_background_relpath(cls, v: Optional[str]) -> Optional[str]:
        return cls._normalize_relpath(v, field_name="instrument_background_relpath")

    @staticmethod
    def _normalize_relpath(v: Optional[str], *, field_name: str) -> Optional[str]:
        if v is None:
            return None
        rel = v.strip()
        if not rel:
            return None

        # Normalize path separators so ledgers created on Windows can be read on POSIX.
        rel_normalized = rel.replace("\\", "/")

        if PureWindowsPath(rel).is_absolute() or PurePosixPath(rel_normalized).is_absolute():
            raise ValueError(f"{field_name} must be relative to SampleLedger.sample_path")

        path = PurePosixPath(rel_normalized)
        if ".." in path.parts:
            raise ValueError(f"{field_name} cannot contain '..'")

        return path.as_posix()

    @model_validator(mode="after")
    def validate_total_and_background(self) -> "SpectrumEntry":
        if self.spectrum_relpath is None:
            raise ValueError("spectrum_relpath is required for SpectrumEntry")

        return self
