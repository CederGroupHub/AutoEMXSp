#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loads the available segmentation models.

This module is responsible for loading segmentation model modules that follow a specific 
structure and provide a `segment_particles()` function. To add a new segmentation model, 
create a new Python file in the same directory as this script, and ensure the file 
contains a correctly defined `segment_particles()` function as described below.

### Adding New Segmentation Models

See segmentation_model_template.py for a detailed description of the required structure

1. **Create a new module**:
   In this directory, create a new Python file named `your_model_name.py`

2. **Define `segment_particles()`**:
   In your new module, define a function named `segment_particles()` that takes the following parameters:
   - `frame_image` (np.array): A grayscale image containing particles.
   - `powder_meas_config` (PowderMeasurementConfig): A configuration object with relevant parameters for segmentation (e.g., threshold values).
   - `save_image` (bool): An optional flag to save the masked image.
   - `EM` (EM_controller): An optional object for saving the segmented image.

   The function should return a binary mask where pixel values represent detected particles.

Created on Thu Oct  9 09:34:39 2025

@author: Andrea
"""


import os
import importlib

import autoemxsp.tools.constants as cnst

# Dictionary to hold the available models
PAR_SEGMENTATION_MODEL_REGISTRY = {}

def load_models():
    """
    Loads all available segmentation model modules from the specified directory.

    This function will iterate through all `.py` files in the directory, import
    them, and store them in the MODEL_REGISTRY for later use.

    Returns
    -------
    module_names : list
        A list of the names of the modules that have been successfully loaded.
    """
    module_names = []
    # Loop through all Python files in the directory
    for model_file in os.listdir(os.path.dirname(__file__)):
        # Ignore the 'ml_segmentation_template.py' file
        if model_file == 'segmentation_model_template.py':
            continue
        
        if model_file.endswith(".py") and model_file != "__init__.py":
            model_name = model_file[:-3]  # Remove .py extension
            module = importlib.import_module(f'autoemxsp.core.{cnst.PAR_SEGMENTATION_MODELS_DIR}.{model_name}')
            PAR_SEGMENTATION_MODEL_REGISTRY[model_name] = module  # Save the module to the registry
            module_names.append(model_name)
    return module_names

# Automatically load user-added models
AVAILABLE_SEGMENTATION_MODELS = load_models()