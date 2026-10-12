#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate the background functions (modified Duncumb model) of the detector model.

Refits spectra of standards with the chosen background parameters (P, F, beta) fitted freely, and fits the functions of
the detector model (of the mean atomic number Z) to their values: P = a2*Z**2 + a1*Z + a0,
F = 1 - exp(-(Z + z_offset) / z_scale), beta = sigmoid + Gaussian. See autoemx.runners.calibrate_detector_model.

Typical usage:
    1. Acquire high-count spectra of bulk standards, ideally pure elements (beta is attributed to their Z), covering
       a range of atomic numbers (e.g. with Run_Experimental_Standard_Collection.py).
    2. List the standards below.
    3. Since the parameters are correlated, calibrate one function at a time, updating detector_model.json after each
       step (order of Giunto et al., Supplementary Note "Calibration of P, F and beta background model parameters"):
         a. params = ['P', 'F', 'beta'], fit_functions = ['P']   (P sets the curvature of the continuum)
         b. params = ['beta'],           fit_functions = ['beta']
         c. params = ['F'],              fit_functions = ['F']
       Repeat the sequence if the values change significantly.
    4. Run the script, check the plots, and paste the printed 'background' into
       autoemx/calibrations/<microscope_ID>/detector_model.json.

Created on Mon Oct 12 2026

@author: Andrea
"""

# =============================================================================
# Input
# =============================================================================

# Standards: paths of the sample folders
samples = [
    '/absolute/path/to/Std_measurements/C',
    '/absolute/path/to/Std_measurements/Al',
    '/absolute/path/to/Std_measurements/Cu',
]

# Background parameters fitted freely in the spectra, among 'P', 'F', 'beta' (the others are set by the detector model)
params = ['P', 'F', 'beta']

# Parameters whose functions are fitted to their values (among params)
fit_functions = ['P']

# Coefficients fitted for each parameter (None: all). The others keep their values. E.g. {'P': ['a0']}
# P: 'a2', 'a1', 'a0'; F: 'z_offset', 'z_scale'; beta: 'sigmoid_offset', 'sigmoid_amplitude', 'sigmoid_steepness',
# 'sigmoid_center', 'gaussian_amplitude', 'gaussian_center', 'gaussian_sigma'
coeffs = None

max_spectra_per_sample = 15  # None: all. Each standard was averaged over ~15 spectra of 250k counts in the paper

# =============================================================================
# Run
# =============================================================================

from autoemx.runners.calibrate_detector_model import calibrate_background_params

calibrate_background_params(samples=samples, params=params, fit_functions=fit_functions, coeffs=coeffs,
                            max_spectra_per_sample=max_spectra_per_sample)
