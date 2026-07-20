#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI test for multi-particle frame segmentation on a static TIFF."""

from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder

from _image_helpers import load_grayscale_image, make_dev_controller
from ci_paths import INPUTS_DIR


def test_segmentation_particles_in_frame(outputs_dir):
    image, pixel_size_um = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    em, powder_meas_cfg = make_dev_controller(str(outputs_dir))
    em.pixel_size_um = pixel_size_um
    particle_finder = EM_Particle_Finder(em, powder_meas_cfg, results_dir=str(outputs_dir))

    result = particle_finder._get_particles_on_substrate_mask(
        frame_image=image,
        save_image=True,
    )

    assert result is not None
    mask, _mask_path = result
    assert mask.shape[:2] == image.shape[:2]
    assert mask.any()
