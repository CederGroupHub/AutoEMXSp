#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting module for EDS spectrum analysis.

This package provides classes and functions for physically-accurate modeling, fitting, and
analysis of X-ray energy-dispersive spectroscopy (EDS) spectra.

Main Components:
    - **fitter.XSp_Fitter**: Main orchestrator for spectrum fitting
    - **detector_response.DetectorResponseFunction**: Detector response and convolution handling
    - **peaks.Peaks_Model**: Peak modeling and parameterization
    - **background.Background_Model**: Background continuum modeling

Typical Usage::

    from autoemx.core.fitter import XSp_Fitter

    fitter = XSp_Fitter(
        spectrum_vals, energy_vals,
        spectrum_lims=spectrum_lims,
        microscope_ID='PhenomXL',
        ...
    )
    fit_result, fitted_lines = fitter.fit_spectrum(plot_result=True)
"""

from .fitter import XSp_Fitter
from .detector_response import DetectorResponseFunction
from .peaks import Peaks_Model
from .background import Background_Model

__all__ = [
    'XSp_Fitter',
    'DetectorResponseFunction',
    'Peaks_Model',
    'Background_Model',
]
