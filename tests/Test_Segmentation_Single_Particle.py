#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle Segmentation Test Script
=================================

This script is designed to test the particle segmentation algorithms used in the 
`autoemxsp` framework. It loads a sample particle image in TIFF format, processes it, 
and runs segmentation using the `EM_Particle_Finder` module. The segmentation results 
are saved to a dedicated output directory for further validation.

The purpose of this script is to provide a reproducible test workflow for verifying 
algorithm performance using controlled input data.

Usage
-----
Simply place your test image (`example_particle_image.tif`) in the script directory 
and run:

    python3 test_particle_segmentation.py

The results will be stored in the `outputs` folder created in the same directory.

Author: Andrea
Created: Mon Oct 13 15:38:36 2025
"""

import os
import json
import cv2
import tifffile

from autoemxsp.core.EM_particle_finder import EM_Particle_Finder
from autoemxsp.core.EM_controller import EM_Controller
from autoemxsp.tools.config_classes import (
    PowderMeasurementConfig,
    MicroscopeConfig,
    MeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
    BulkMeasurementConfig
)
from autoemxsp.tools.utils import print_single_separator

# -------------------------------------------------------------------------
# Configuration Parameters
# -------------------------------------------------------------------------
microscope_ID = 'PhenomXL'
microscope_type = 'SEM'

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, 'outputs')

# Ensure output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# -------------------------------------------------------------------------
# Load Test Image
# -------------------------------------------------------------------------
test_image_path = os.path.join(script_dir, 'inputs', 'example_particle_image.tif')
IMAGE_PIXEL_SIZE_um = 1 #um

with tifffile.TiffFile(test_image_path) as tif:
    # Read image data into a numpy array
    image = tif.asarray()
    # Extract ImageDescription metadata from page 0
    description_str = tif.pages[0].description

# Parse the JSON description (currently not used for pixel size)
description_dict = json.loads(description_str)

# Load pixel size
ps_key = 'pixel_size_um'
if ps_key in description_dict:
    pixel_size_um = description_dict['pixel_size_um']
else:
    print("Pixel size could not be extracted from loaded image.")
    print("ENSURE VALUE OF 'IMAGE_PIXEL_SIZE_um' IS CORRECT IF NECESSARY TO THE SEGMENTATION MODEL")
    pixel_size_um = IMAGE_PIXEL_SIZE_um
    print_single_separator()


# -------------------------------------------------------------------------
# Initialize Microscope and Measurement Configurations
# -------------------------------------------------------------------------
microscope_cfg = MicroscopeConfig(microscope_ID, microscope_type)
sample_cfg = SampleConfig('Segmentation_test', [])
measurement_cfg = MeasurementConfig()
sample_substrate_cfg = SampleSubstrateConfig()
bulk_meas_cfg = BulkMeasurementConfig()

# Select here segmentation model 
powder_meas_cfg = PowderMeasurementConfig(par_segmentation_model = PowderMeasurementConfig.DEFAULT_PAR_SEGMENTATION_MODEL)

# Controller handles EM workflow
EM = EM_Controller(
    microscope_cfg,
    sample_cfg,
    measurement_cfg,
    sample_substrate_cfg,
    powder_meas_cfg,
    bulk_meas_cfg,
    development_mode=True
)

# Particle finder instance
particle_finder = EM_Particle_Finder(EM, powder_meas_cfg, results_dir=output_dir)

# -------------------------------------------------------------------------
# Pre-process Image (convert to grayscale if necessary)
# -------------------------------------------------------------------------
if image.ndim == 3:
    if image.shape[2] == 3:
        # RGB to grayscale
        image = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
    elif image.shape[2] == 4:
        # RGBA to grayscale (drop alpha channel)
        image = cv2.cvtColor(image, cv2.COLOR_RGBA2GRAY)

# -------------------------------------------------------------------------
# Run Particle Segmentation and Save Results
# -------------------------------------------------------------------------
particle_finder._get_particle_mask(
    par_image=image,
    pixel_size_um=pixel_size_um
)

print_single_separator()
print(f"Segmentation completed. Results saved in: {output_dir}")