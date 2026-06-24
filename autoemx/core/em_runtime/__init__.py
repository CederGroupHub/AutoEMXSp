#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Runtime modules for electron microscopy acquisition and control."""

__all__ = [
    "MicroscopeController",
    "FrameNavigator",
    "SpectrumAcquisition",
    "EM_Sample_Finder",
    "EM_Particle_Finder",
    "XSpSpotSelectionContext",
    "XSpSpotSelectorCallback",
    "validate_xsp_spot_pixels",
]


def __getattr__(name):
    if name == "MicroscopeController":
        from autoemx.core.em_runtime.microscope_controller import MicroscopeController
        return MicroscopeController
    if name == "FrameNavigator":
        from autoemx.core.em_runtime.frame_navigator import FrameNavigator
        return FrameNavigator
    if name == "SpectrumAcquisition":
        from autoemx.core.em_runtime.spectrum_acquisition import SpectrumAcquisition
        return SpectrumAcquisition
    if name == "EM_Sample_Finder":
        from autoemx.core.em_runtime.sample_finder import EM_Sample_Finder
        return EM_Sample_Finder
    if name == "EM_Particle_Finder":
        from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder
        return EM_Particle_Finder
    if name in ("XSpSpotSelectionContext", "XSpSpotSelectorCallback", "validate_xsp_spot_pixels"):
        from autoemx.core.em_runtime import xsp_spot_selection
        return getattr(xsp_spot_selection, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
