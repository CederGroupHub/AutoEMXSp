#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI test for single-particle segmentation on a static TIFF."""

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
