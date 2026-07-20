#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI test for X-ray spot selection on a segmented particle image."""

from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder

from _image_helpers import load_grayscale_image, make_dev_controller
from ci_paths import INPUTS_DIR


def test_xsp_measurement_spot_selection(outputs_dir):
    image, pixel_size_um = load_grayscale_image(INPUTS_DIR / "example_particle_image.tif")
    em, powder_meas_cfg = make_dev_controller(str(outputs_dir))
    particle_finder = EM_Particle_Finder(em, powder_meas_cfg, results_dir=str(outputs_dir))

    spots = particle_finder.get_XS_acquisition_spots_coord_list(
        n_tot_sp_collected=0,
        par_image=image,
        pixel_size_um=pixel_size_um,
    )

    assert spots is not None
    assert len(spots) >= 1
