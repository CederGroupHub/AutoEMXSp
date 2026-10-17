#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Frames acquired on a sample (``SampleLedger.frames``): one record per frame, with its segmentation results."""

from datetime import datetime
from typing import List, Literal, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from .acquisition import Coordinate2D
from .segmentation import FrameSegmentation

#: Why a frame was acquired: searched for particles to measure by EDS, or by a particle size statistics measurement
FramePurpose = Literal["eds_search", "size_stats"]


class FrameRecord(BaseModel):
    """
    A frame acquired on the sample: what the acquisition fixed (position of its centre on the stage, if known, frame
    width, pixel size, area and saved image), and the particles found in it by each segmentation (never replaced:
    segmenting it again adds a result).

    ``frame_id`` is unique in the ledger and names the image file (``<sample>_fr<frame_id>_particles``); ``label`` is
    the label of the frame navigation (``A0``, ``A0_a0``), as in ``AcquisitionDetails.frame_id`` of the spectra
    acquired in it. ``search_window_um`` are the particle diameters searched (particle area limits of the
    acquisition); ``measurement_id`` is the particle size statistics measurement that acquired it.
    """

    frame_id: str
    label: str
    purpose: FramePurpose
    measurement_id: Optional[int] = None
    datetime: Optional[str] = Field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))
    center_mm: Optional[Coordinate2D] = None
    frame_width_um: float
    pixel_size_um: float
    area_um2: float
    search_window_um: Tuple[float, float]
    image_file: Optional[str] = None
    segmentations: List[FrameSegmentation] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

    @field_validator("frame_width_um", "pixel_size_um", "area_um2")
    @classmethod
    def validate_positive(cls, v: float) -> float:
        if not np.isfinite(v) or v <= 0:
            raise ValueError("Frame width, pixel size and area must be finite and positive")
        return v

    @model_validator(mode="after")
    def validate_record(self) -> "FrameRecord":
        ids = [s.segmentation_id for s in self.segmentations]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Frame {self.frame_id}: one result per segmentation")
        if (self.purpose == "size_stats") != (self.measurement_id is not None):
            raise ValueError(f"Frame {self.frame_id}: measurement_id is set for (and only for) size_stats frames")
        return self

    def segmentation(self, segmentation_id: Optional[int]) -> Optional[FrameSegmentation]:
        """Result of a segmentation of the frame, or None if not segmented with it."""
        return next((s for s in self.segmentations if s.segmentation_id == segmentation_id), None)

    def add_segmentation(self, result: FrameSegmentation) -> None:
        """Add the result of a segmentation (one per segmentation id)."""
        if self.segmentation(result.segmentation_id) is not None:
            raise ValueError(f"Frame {self.frame_id} is already segmented with segmentation {result.segmentation_id}")
        self.segmentations.append(result)
