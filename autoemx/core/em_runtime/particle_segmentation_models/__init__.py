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

   The function should return a binary mask where pixel values represent detected particles, or a label image
   (one integer per particle).

3. **Optional dependencies** (e.g. torch): import them inside `segment_particles()`, never at the top of the
   module, so that the module always imports (ledgers recorded with the model stay loadable where it is not
   installed). Define `is_available() -> bool` (whether the dependencies are installed) and `INSTALL_HINT` (how to
   install them): only installed models are offered in the GUI, and using a model that is not installed raises an
   ImportError with the hint.

Modules whose name starts with '_' are helpers, not models.

Created on Thu Oct  9 09:34:39 2025

@author: Andrea
"""


import os
import importlib

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
    for model_file in sorted(os.listdir(os.path.dirname(__file__))):
        # Ignore the 'ml_segmentation_template.py' file
        if model_file == 'segmentation_model_template.py':
            continue
        
        if model_file.endswith(".py") and not model_file.startswith("_"):
            model_name = model_file[:-3]  # Remove .py extension
            module = importlib.import_module(f'{__package__}.{model_name}')
            PAR_SEGMENTATION_MODEL_REGISTRY[model_name] = module  # Save the module to the registry
            module_names.append(model_name)
    return module_names

# Automatically load user-added models
AVAILABLE_SEGMENTATION_MODELS = load_models()


def is_model_available(model_name: str) -> bool:
    """Whether the dependencies of a registered model are installed (models without `is_available()` are)."""
    module = PAR_SEGMENTATION_MODEL_REGISTRY[model_name]
    check = getattr(module, "is_available", None)
    if check is None:
        return True
    try:
        return bool(check())
    except Exception:
        return False


def available_models():
    """Registered models whose dependencies are installed."""
    return [name for name in AVAILABLE_SEGMENTATION_MODELS if is_model_available(name)]


def install_hint(model_name: str) -> str:
    """How to install the dependencies of a model ('' if none given)."""
    return getattr(PAR_SEGMENTATION_MODEL_REGISTRY[model_name], "INSTALL_HINT", "")


def require_model(model_name: str):
    """The module of a model, or ImportError (with how to install it) if its dependencies are not installed."""
    if not is_model_available(model_name):
        hint = install_hint(model_name)
        raise ImportError(f"Particle segmentation model '{model_name}' is not installed."
                          + (f" Install it with: {hint}" if hint else ""))
    return PAR_SEGMENTATION_MODEL_REGISTRY[model_name]
