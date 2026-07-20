#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Image loading helpers for CI tests (no tifffile dependency)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Tuple

import cv2
import numpy as np
from PIL import Image


def load_grayscale_image(image_path: Path) -> Tuple[np.ndarray, float]:
    """Load an image as grayscale and optional pixel size from TIFF metadata."""
    path = Path(image_path)
    pixel_size_um = 1.0

    with Image.open(path) as img:
        description = img.tag_v2.get(270) if hasattr(img, "tag_v2") else None
        if isinstance(description, bytes):
            description = description.decode("utf-8", errors="ignore")
        if isinstance(description, str) and description.strip():
            try:
                meta = json.loads(description)
                if isinstance(meta, dict) and "pixel_size_um" in meta:
                    pixel_size_um = float(meta["pixel_size_um"])
            except json.JSONDecodeError:
                pass
        image = np.asarray(img)

    if image.ndim == 3:
        if image.shape[2] == 3:
            image = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
        elif image.shape[2] == 4:
            image = cv2.cvtColor(image, cv2.COLOR_RGBA2GRAY)

    return image, pixel_size_um


def make_dev_controller(results_dir: str) -> Any:
    """Build an EM_Controller in development mode for offline image tests."""
    from autoemx.config.runtime_configs import (
        BulkMeasurementConfig,
        MeasurementConfig,
        MicroscopeConfig,
        PowderMeasurementConfig,
        SampleConfig,
        SampleSubstrateConfig,
    )
    from autoemx.core.em_runtime.controller import EM_Controller

    microscope_cfg = MicroscopeConfig(ID="PhenomXL", type="SEM")
    sample_cfg = SampleConfig(elements=[])
    measurement_cfg = MeasurementConfig()
    sample_substrate_cfg = SampleSubstrateConfig()
    bulk_meas_cfg = BulkMeasurementConfig()
    powder_meas_cfg = PowderMeasurementConfig(
        par_segmentation_model=PowderMeasurementConfig.DEFAULT_PAR_SEGMENTATION_MODEL
    )

    return (
        EM_Controller.from_configs(
            microscope_cfg,
            sample_cfg,
            measurement_cfg,
            sample_substrate_cfg,
            powder_meas_cfg,
            bulk_meas_cfg,
            sample_id="",
            development_mode=True,
        ),
        powder_meas_cfg,
    )
