#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Electron Microscope Driver Loader

This module provides a function to dynamically load and inject
microscope-specific driver parameters and functions for SEM operation
from the 'microscope_drivers' directory.

Usage from other modules:
    from autoemx import EM_driver
    EM_driver.load_microscope_driver(microscope_ID='PhenomXL')

Author: Andrea Giunto
Created on: Mon Jan 20 15:40:42 2025
"""

import os
import importlib
from types import ModuleType
from typing import Any, Optional

_active_driver_module: Optional[ModuleType] = None
_active_driver_id: Optional[str] = None

def load_microscope_driver(microscope_ID: str) -> None:
    """
    Dynamically load driver parameters and functions for a given microscope.

    This function imports all public attributes from the driver file of the microscope
    (located inside the microscope_drivers directory) and injects them into the current module's namespace.
    The driver file is '{microscope_ID}.py', unless the detector model of the microscope
    (autoemx/calibrations/<microscope_ID>/detector_model.json) names another driver, shared by
    instruments of the same model.

    Args
    ----
        microscope_ID (str): The name of the microscope.

    Raises
    ------
        ValueError: If the driver file cannot be found or imported.

    Warning
    -------
        This function injects variables and functions into the module namespace.
        Use with care to avoid name collisions.
    """
    from autoemx.calibrations import microscope_driver_name

    # Build the path to the driver file
    driver_name = microscope_driver_name(microscope_ID)
    drivers_dir = os.path.join(os.path.dirname(__file__))
    driver_file = os.path.join(drivers_dir, f"{driver_name}.py")
    if not os.path.isfile(driver_file):
        raise ValueError(
            f"Could not find the microscope driver file at '{driver_file}'.\n"
            f"Please ensure the driver of microscope '{microscope_ID}' ('{driver_name}') matches a .py file in "
            "'microscope_drivers' (the driver is the microscope ID, unless set as 'driver' in its detector_model.json)."
        )

    # Import the driver module dynamically
    module_name = f"autoemx.microscope_drivers.{driver_name}"
    pkg = __package__ if __package__ else __name__
    try:
        mod = importlib.import_module(module_name, package=pkg)
    except ModuleNotFoundError as e:
        raise ValueError(
            f"Could not import driver module '{module_name}'.\n"
            f"Tried to import: {module_name} (relative to package '{pkg}')."
        ) from e

    global _active_driver_module, _active_driver_id
    _active_driver_module = mod
    _active_driver_id = microscope_ID


def get_active_driver_module() -> ModuleType:
    """Return currently loaded microscope driver module."""
    if _active_driver_module is None:
        raise RuntimeError("No microscope driver loaded. Call load_microscope_driver first.")
    return _active_driver_module


def connect_to_microscope(warn_if_unavailable: bool = True) -> bool:
    """Connect using the active driver, if supported."""
    mod = get_active_driver_module()
    connect_fn = getattr(mod, "connect_to_microscope", None)
    if connect_fn is None:
        raise AttributeError(
            f"Driver '{_active_driver_id}' does not define connect_to_microscope()."
        )
    return bool(connect_fn(warn_if_unavailable=warn_if_unavailable))


def is_microscope_connected() -> bool:
    """Return True if the active driver reports a live microscope connection."""
    if _active_driver_module is None:
        return False

    is_connected_fn = getattr(_active_driver_module, "is_microscope_connected", None)
    if is_connected_fn is not None:
        try:
            return bool(is_connected_fn())
        except Exception:
            return False

    return bool(getattr(_active_driver_module, "is_at_EM", False))


def __getattr__(name: str) -> Any:
    """Delegate unresolved attributes to the active microscope driver module.

    Raises AttributeError (not RuntimeError) when no driver is loaded, so that hasattr() and getattr() with a
    default keep working, e.g. for introspection tools.
    """
    if _active_driver_module is None or name.startswith("__"):
        raise AttributeError(
            f"module '{__name__}' has no attribute '{name}' (no microscope driver loaded: "
            "call load_microscope_driver first)."
        )
    mod = _active_driver_module
    try:
        return getattr(mod, name)
    except AttributeError as exc:
        raise AttributeError(f"Active driver '{_active_driver_id}' has no attribute '{name}'.") from exc


def __dir__() -> list[str]:
    attrs = set(globals().keys())
    if _active_driver_module is not None:
        attrs.update(attr for attr in dir(_active_driver_module) if not attr.startswith("_"))
    return sorted(attrs)