#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate the peak shape (incomplete charge collection) functions of the detector model.

Refits spectra of standards with the chosen peak shape parameters fitted freely, and fits the sigmoids of the detector
model (functions of the X-ray energy) to their values. gamma and f_tail (low-energy tailing) apply to lines below
1.84 keV, R_e and F_loss (incomplete charge collection) to lines between 1.18 and 5 keV.
See autoemx.runners.calibrate_detector_model for details.

Typical usage:
    1. Acquire high-count spectra of bulk standards, ideally pure elements, with lines spanning the energy range
       (e.g. with Run_Experimental_Standard_Collection.py).
    2. List the standards below, with the elements whose peak shapes are fitted freely.
    3. Choose the parameters to fit, and the sigmoid coefficients to fit to their values (more lines than
       coefficients are needed). Since the parameters are correlated, preferably calibrate one parameter at a time,
       updating detector_model.json after each one, and repeat until the values are stable (see the documentation for
       the recommended order).
    4. Run the script, check the plots, and paste the printed 'peak_shape' into
       autoemx/calibrations/<microscope_ID>/detector_model.json.

Created on Mon Oct 12 2026

@author: Andrea
"""

# =============================================================================
# Input
# =============================================================================

# Standards: {path of the sample folder: elements whose peak shapes are fitted freely}
samples = {
    '/absolute/path/to/Std_measurements/C': ['C'],
    '/absolute/path/to/Std_measurements/Al': ['Al'],
}

# Peak shape parameters fitted freely, among 'gamma', 'f_tail', 'R_e', 'F_loss'
params = ['gamma']

# Sigmoid coefficients fitted, among 'offset', 'amplitude', 'center', 'steepness' (the others keep their values)
sigmoid_coeffs = ['offset', 'amplitude', 'center', 'steepness']

max_spectra_per_sample = 10  # None: all. Ten spectra of 250k counts per standard were used in the paper

# =============================================================================
# Run
# =============================================================================

from autoemx.runners.calibrate_detector_model import calibrate_peak_shape_params

calibrate_peak_shape_params(samples=samples, params=params, sigmoid_coeffs=sigmoid_coeffs,
                            max_spectra_per_sample=max_spectra_per_sample)
