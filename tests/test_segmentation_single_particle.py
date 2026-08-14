#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI test for single-particle segmentation on a static TIFF."""

import numpy as np

from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder

from _image_helpers import load_grayscale_image, make_dev_controller
from ci_paths import INPUTS_DIR


def test_segmentation_single_particle(outputs_dir):
    image, pixel_size_um = load_grayscale_image(INPUTS_DIR / "example_particle_image.tif")
    em, powder_meas_cfg = make_dev_controller(str(outputs_dir))
    particle_finder = EM_Particle_Finder(em, powder_meas_cfg, results_dir=str(outputs_dir))

    result = particle_finder._get_particle_mask(
        par_image=image,
        pixel_size_um=pixel_size_um,
    )

    assert result is not None
    mask, _par_image = result
    assert mask.shape[:2] == image.shape[:2]
    assert mask.any()


def test_particle_mask_uses_recentered_frame(outputs_dir, monkeypatch):
    """After a microscope recenter, mask and returned image must come from the recapture."""
    height = width = 200
    off_center = np.zeros((height, width), dtype=np.uint8)
    off_center[20:90, 20:90] = 200
    centered = np.zeros((height, width), dtype=np.uint8)
    centered[65:135, 65:135] = 200

    em, powder_meas_cfg = make_dev_controller(str(outputs_dir))
    em.pixel_size_um = 0.1
    particle_finder = EM_Particle_Finder(em, powder_meas_cfg, results_dir=str(outputs_dir))
    particle_finder._im_width = width
    particle_finder._im_height = height

    captures = {"n": 0}

    def fake_capture(par_image=None, pixel_size_um=None):
        captures["n"] += 1
        return off_center if captures["n"] == 1 else centered

    monkeypatch.setattr(
        "autoemx.core.em_runtime.particle_finder.EM_driver.is_microscope_connected",
        lambda: True,
    )
    monkeypatch.setattr(particle_finder, "_capture_ref_frame", fake_capture)
    monkeypatch.setattr(em, "move_to_pos", lambda pos: None)
    monkeypatch.setattr(em, "convert_pixel_pos_to_mm", lambda pix: np.array([0.0, 0.0]))

    result = particle_finder._get_particle_mask()
    assert result is not None
    mask, returned_image = result

    assert np.array_equal(returned_image, centered)
    assert mask[100, 100] == 255
    assert mask[40, 40] == 0
    assert captures["n"] == 2
