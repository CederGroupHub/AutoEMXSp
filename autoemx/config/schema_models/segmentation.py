#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle segmentations of the frames of a sample (``SampleLedger.segmentations``), and the result of one segmentation
of one frame (``FrameRecord.segmentations``).

Like quantifications, each segmentation has an id: segmenting frames again (with another model or threshold) adds
results with a new id to the frames, and never replaces earlier ones, so that the statistics of different
segmentations can be compared and one chosen.
"""

import json
from datetime import datetime
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from autoemx._version import __version__ as AUTOEMX_VERSION
from autoemx.config.runtime_configs import PowderMeasurementConfig


class SegmentationConfig(BaseModel):
    """
    Settings of a particle segmentation of the frames of a sample: segmentation model, brightness threshold (used by
    the threshold model), minimum number of pixels of a particle for its area to be measured, and options of the
    model (``PowderMeasurementConfig.par_segmentation_params``, e.g. the SAM2 checkpoint). ``source`` tells
    whether frames were segmented with it during an acquisition, or again afterwards from their saved images.
    """

    segmentation_id: int
    label: str = ""
    autoemx_version: str = AUTOEMX_VERSION
    datetime: str = Field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    source: Literal["acquisition", "reanalysis"] = "acquisition"
    par_segmentation_model: str
    par_brightness_thresh: int
    min_par_pixels: int
    model_params: Dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")

    @field_validator("segmentation_id", mode="before")
    @classmethod
    def validate_id(cls, v: Any) -> Any:
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            raise ValueError("segmentation_id must be a non-negative integer")
        return v

    @model_validator(mode="after")
    def default_label(self) -> "SegmentationConfig":
        if not self.label:
            self.label = segmentation_label(self.par_segmentation_model, self.par_brightness_thresh,
                                            self.min_par_pixels, self.model_params)
        return self

    def settings(self) -> Dict[str, Any]:
        """The segmentation settings (not the id, version, date, source or label); the brightness threshold only for
        the threshold model, the only one using it."""
        settings = {"par_segmentation_model": self.par_segmentation_model, "min_par_pixels": self.min_par_pixels,
                    "model_params": self.model_params}
        if self.par_segmentation_model == PowderMeasurementConfig.DEFAULT_PAR_SEGMENTATION_MODEL:
            settings["par_brightness_thresh"] = self.par_brightness_thresh
        return settings

    def fingerprint(self) -> str:
        """Identity of the settings: a segmentation with the same settings is reused, not added again."""
        return json.dumps(self.settings(), sort_keys=True, default=str)

    @classmethod
    def from_powder_meas_cfg(cls, cfg: PowderMeasurementConfig, segmentation_id: int,
                             source: str = "acquisition") -> "SegmentationConfig":
        """Segmentation with the settings of a measurement."""
        return cls(segmentation_id=segmentation_id, source=source, par_segmentation_model=cfg.par_segmentation_model,
                   par_brightness_thresh=cfg.par_brightness_thresh, min_par_pixels=cfg.min_par_pixels,
                   model_params=dict(cfg.par_segmentation_params))


def segmentation_label(model: str, threshold: int, min_par_pixels: int,
                       model_params: Optional[Dict[str, Any]] = None) -> str:
    """Readable name of segmentation settings, e.g. 'threshold_bright · 100 · 80 px', or
    'SAM2 · 80 px · checkpoint=sam2.1-hiera-small' (model options, values shortened to their last path part)."""
    parts = [model] + ([str(threshold)] if model == PowderMeasurementConfig.DEFAULT_PAR_SEGMENTATION_MODEL else [])
    parts.append(f"{min_par_pixels} px")
    if model_params:
        parts.append(", ".join(f"{k}={str(v).rsplit('/', 1)[-1]}" for k, v in sorted(model_params.items())))
    return " · ".join(parts)


class FrameSegmentation(BaseModel):
    """
    Particles found in a frame by one segmentation (``segmentation_id``): area, frame-edge correction and mark
    (x, y and radius in pixels, on the frame image) of each particle valid for the statistics: in ``d_range_um``
    (the diameters searched in the frame, from the smallest with at least ``min_par_pixels`` pixels) and not touching
    the frame edge. Particle ``i`` of the frame is identified by ``(frame_id, i)``.
    """

    segmentation_id: int
    d_range_um: Tuple[float, float]
    particle_areas_um2: List[float] = Field(default_factory=list)
    particle_edge_corrections: List[float] = Field(default_factory=list)
    particle_marks_px: List[Tuple[int, int, int]] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def validate_particles(self) -> "FrameSegmentation":
        if len(self.particle_edge_corrections) != len(self.particle_areas_um2):
            raise ValueError("particle_areas_um2 and particle_edge_corrections must have one value per particle")
        if self.particle_marks_px and len(self.particle_marks_px) != len(self.particle_areas_um2):
            raise ValueError("particle_marks_px must have one mark per particle (or none)")
        if any(not np.isfinite(a) or a <= 0 for a in self.particle_areas_um2):
            raise ValueError("particle_areas_um2 must be finite and positive")
        if any(not np.isfinite(c) or c < 1 for c in self.particle_edge_corrections):
            raise ValueError("particle_edge_corrections must be finite and at least 1")
        lo, hi = self.d_range_um
        if not (np.isfinite(lo) and np.isfinite(hi)) or lo < 0 or hi < lo:
            raise ValueError("d_range_um must be finite, non-negative and increasing")
        return self

    @property
    def n_particles(self) -> int:
        return len(self.particle_areas_um2)
