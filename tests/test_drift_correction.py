#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI tests for image drift estimation and correction."""

from pathlib import Path

import cv2
import numpy as np

from autoemx.core.em_runtime.drift_correction import estimate_drift_of, shift_image

from ci_paths import INPUTS_DIR


def _load_grayscale(image_path: Path) -> np.ndarray:
    image = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise FileNotFoundError(f"Could not load image: {image_path}")
    return image


def test_drift_correction_estimates_nonzero_shift():
    original = _load_grayscale(INPUTS_DIR / "pre_drift.png")
    shifted = _load_grayscale(INPUTS_DIR / "post_drift.png")

    estimated_dx, estimated_dy = estimate_drift_of(original, shifted)
    reshifted = shift_image(shifted, -estimated_dx, -estimated_dy)

    assert np.isfinite(estimated_dx)
    assert np.isfinite(estimated_dy)
    assert abs(estimated_dx) + abs(estimated_dy) > 0.1

    before = np.mean(np.abs(original.astype(np.float32) - shifted.astype(np.float32)))
    after = np.mean(np.abs(original.astype(np.float32) - reshifted.astype(np.float32)))
    assert after < before
