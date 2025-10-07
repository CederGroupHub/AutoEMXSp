#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Electron Microscope Driver Loader

This module provides a function to dynamically load and inject
microscope-specific driver parameters and functions for SEM operation
from the 'EM_driver' directory.

Usage from other modules:
    from autoemxsp import EM_driver
    EM_driver.load_microscope_driver(microscope_ID='PhenomXL')

Author: Andrea Giunto
Created on: Mon Jan 20 15:40:42 2025
"""

import os
import sys
import importlib

def load_microscope_driver(microscope_ID: str) -> None:
    """
    Dynamically load driver parameters and functions for a given microscope.

    This function imports all public attributes from the driver file
    named '{microscope_ID}.py' (located inside the EM_driver directory)
    and injects them into the current module's namespace.

    Args
    ----
        microscope_ID (str): The name of the microscope (must match a .py file in the EM_driver directory).

    Raises
    ------
        ValueError: If the driver file cannot be found or imported.

    Warning
    -------
        This function injects variables and functions into the module namespace.
        Use with care to avoid name collisions.
    """
    # Build the path to the driver file
    drivers_dir = os.path.join(os.path.dirname(__file__))
    driver_file = os.path.join(drivers_dir, f"{microscope_ID}.py")
    if not os.path.isfile(driver_file):
        raise ValueError(
            f"Could not find the microscope driver file at '{driver_file}'.\n"
            f"Please ensure microscope_ID ('{microscope_ID}') matches a .py file in 'EM_driver'."
        )

    # Import the driver module dynamically
    module_name = f"autoemxsp.EM_driver.{microscope_ID}"
    try:
        pkg = __package__ if __package__ else __name__
        mod = importlib.import_module(module_name, package=pkg)
    except ModuleNotFoundError as e:
        raise ValueError(
            f"Could not import driver module '{module_name}'.\n"
            f"Tried to import: {module_name} (relative to package '{pkg}')."
        ) from e

    # Inject all public attributes from the driver module into this module
    thismod = sys.modules[__name__]
    for k in dir(mod):
        if not k.startswith('_'):
            setattr(thismod, k, getattr(mod, k))