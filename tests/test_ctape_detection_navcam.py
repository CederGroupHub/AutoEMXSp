#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI test for offline C-tape detection on a navigation-camera image."""

import cv2

from autoemx.core.em_runtime.sample_finder import EM_Sample_Finder

from ci_paths import INPUTS_DIR


def test_ctape_detection_navcam(outputs_dir):
    image_path = INPUTS_DIR / "navcam.tiff"
    image_np = cv2.imread(str(image_path), cv2.IMREAD_UNCHANGED)
    assert image_np is not None, f"Failed to load {image_path}"

    sample_finder = EM_Sample_Finder(
        microscope_ID="PhenomXL",
        center_pos=(22.5, 37.5),
        sample_half_width_mm=3,
        substrate_width_mm=12,
        development_mode=True,
        results_dir=str(outputs_dir),
        verbose=False,
    )
    ctape_coords = sample_finder.detect_Ctape(image_np)

    assert ctape_coords is not None
    center_pos, sample_hw_mm = ctape_coords
    assert len(center_pos) == 2
    assert sample_hw_mm > 0
