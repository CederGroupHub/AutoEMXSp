#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch quantification and analysis of X-ray spectra for a list of samples.

This script provides automated batch quantification and (optionally) clustering/statistical
analysis of acquired X-ray spectra for multiple samples. It is robust to missing files or
errors in individual samples, making it suitable for unattended batch processing.

Run this file directly to process the list of sample IDs with the defined configuration options.

Notes
-----
- Only the `sample_ID` is required if acquisition output is saved in the default directory;
  otherwise, specify `results_path`.
- Designed to continue processing even if some samples are missing or have errors.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""
from autoemxsp.runners import batch_quantify_and_analyze

# =============================================================================
# Examples
# =============================================================================
sample_IDs = [
    'Wulfenite_example',
    'K-412_NISTstd_example'
    ]

# =============================================================================
# Paper data (Download data from github repository:
#       https://github.com/CederGroupHub/AutoEMXSp/tree/main/paper_data)   
# =============================================================================
# sample_IDs = [
    # Minerals
    # 'Alamosite_mineral',
    # 'Albite_mineral',
    # 'Anhydrite_mineral',
    # 'Anorthite_mineral',
    # 'Benitoite_mineral',
    # 'Bornite_mineral',
    # 'Chalcopyrite_mineral',
    # 'CoOlivine_mineral',
    # 'FeOlivine_mineral',
    # 'Fluorphlogopite_mineral',
    # 'Jadeite_mineral',
    # 'K-412_NISTstd_mineral',
    # 'Labradorite_mineral',
    # 'MnOlivine_mineral',
    # 'Nepheline_mineral',
    # 'Orthoclase_mineral',
    # 'Rhodonite_mineral',
    # 'ScPO4_mineral',
    # 'Wulfenite_mineral',
    # 'YIG_mineral',
    # 'YPO4_mineral',
    # ]

results_path = None # Looks in default Results folder if left unspecified
# =============================================================================
# Options
# =============================================================================
max_analytical_error = 5 # w%
min_bckgrnd_cnts = 5

run_clustering_analysis = True

quantify_only_unquantified_spectra = False # Set to True if running on Data.csv file that has already been quantified. Used to quantify discarded unqiantified spectra
interrupt_fits_bad_spectra = True # Interrupts the fit and quantification of spectra when it finds they will lead to large quantification errors. Used to speed up computations

output_filename_suffix = ''

# =============================================================================
# Run
# =============================================================================
comp_analyzer = batch_quantify_and_analyze(
    sample_IDs=sample_IDs,
    quantification_method = 'PB',
    min_bckgrnd_cnts = min_bckgrnd_cnts,
    results_path=results_path,
    output_filename_suffix=output_filename_suffix,
    max_analytical_error=max_analytical_error,
    run_analysis=run_clustering_analysis,
    quantify_only_unquantified_spectra=quantify_only_unquantified_spectra,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
)