#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run quantification and analysis of X-ray spectra for a list of samples.

This script provides automated batch quantification and clustering/statistical
analysis of acquired X-ray spectra for multiple samples. It is robust to missing files or
errors in individual samples, making it suitable for unattended batch processing.

Run `run_quantification_and_analysis.py` directly to process the list of sample IDs with the defined configuration options.

Notes
-----
- Only the `sample_ID` is required if acquisition output is saved in the default directory;
  otherwise, specify `results_path`.
- Clustering/statistical analysis runs by default (runner default: `run_analysis=True`).
- Designed to continue processing even if some samples are missing or have errors.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""
# =============================================================================
# Examples
# =============================================================================
sample_IDs = [
    'Wulfenite_example',
    ]

is_known_precursor_mixture = None # Loads value from config file, if unspecified. Otherwise set to True or False.

import os
results_dir = os.path.dirname(os.path.abspath(__file__)) # Default: save and load results in the same folder as this script. Set to None to use the current working directory, or replace with another path.

# =============================================================================
# Options
# =============================================================================
max_analytical_error = 5 # w% Threhsold value of analytical error above which spectra are filtered out. Only used at the analysis stage, so it does not affect the quantification
min_bckgrnd_cnts = 5 # Minimum value of background counts that a reference peak (used for quantification) has to possess in order for measurement to be valid
    # Spectra not satisfying this are flagged (quant_flag = 8) and not quantified if interrupt_fits_bad_spectra = True. If False, they are still quantified, and filtered out later in the clustering stage
    # If too many spectra end up being flagged, decrease min_bckgrnd_cnts or increase the spectra target total counts
    # If you change min_bckgrnd_cnts, you can requantify the previously-interrupted spectra by setting interrupt_fits_bad_spectra = False


num_CPU_cores = None # Number of cores used during fitting and quantification. If None, selects automatically half the available cores
force_requantification = False # If True, re-quantifies all spectra regardless of existing quantification runs/settings.
interrupt_fits_bad_spectra = True # Interrupts the fit and quantification of spectra when it finds they will lead to large quantification errors. Used to speed up computations. If False, previously interrupted spectra are re-quantified without interruption.
use_project_specific_std_dict = None # If True, loads standards from project folder (i.e. results_dir) during quantification.

output_filename_suffix = '' # Suffix added to Analysis folder and Data.csv file

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze

comp_analyzer = batch_quantify_and_analyze(
    sample_IDs=sample_IDs,
    quantification_method = 'PB',
    min_bckgrnd_cnts = min_bckgrnd_cnts,
    results_path=results_dir,
    output_filename_suffix=output_filename_suffix,
    max_analytical_error=max_analytical_error,
    num_CPU_cores = num_CPU_cores,
    force_requantification=force_requantification,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
    use_project_specific_std_dict = use_project_specific_std_dict,
    is_known_precursor_mixture = is_known_precursor_mixture,
)