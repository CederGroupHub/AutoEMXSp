#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 16:23:45 2025

@author: Andrea
"""
from autoemxsp.core.EM_controller import EM_Sample_Finder

import cv2
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
image_path = os.path.join(script_dir, 'inputs', 'navcam.tiff')

output_dir = os.path.join(script_dir, 'outputs')

# Ensure output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load TIFF image as NumPy array
image_np = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)

sample_finder = EM_Sample_Finder(
    microscope_ID='PhenomXL',
    center_pos=(22.5,37.5),
    sample_half_width_mm=3,
    substrate_width_mm=12,
    development_mode = True,
    results_dir=output_dir,
    verbose=True
)
Ctape_coords = sample_finder.detect_Ctape(image_np)