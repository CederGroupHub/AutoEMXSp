#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated X-Ray Experimental Standard Acquisition and Analysis

This script configures and runs automated collection and fitting
of EDS/WDS spectra from experimental standards (i.e., samples of known composition)
to generate reference values of peak-to-background ratios.

Requirements:
    - Proper instrument calibration files and instrument driver for the selected microscope

Typical usage:
    - Edit the 'std_list' list to define your standards
        Use preferably bulk standards. Powder standards are also acceptable, especially when making standards for analysing known precursor mixtures.
    - Adjust configuration parameters as needed
    - Run the script to collect experimental standards for one or multiple samples at a time

Created on Fri Aug 20 09:34:34 2025

@author: Andrea
"""
# =============================================================================
# General Configuration
# =============================================================================
spectrum_lims = (14, 1100)  # eV

use_instrument_background = False

min_bckgrnd_cnts = 10

output_filename_suffix = ''

exp_std_dir = None # Defines directory where measurements are saved. If None, uses default path.
exp_std_dir = '/Users/Andrea_1/Desktop/Work/Projects/SEM EDX automation/EDX standards/Auto measurements'  # Default directory for saving experimental standards
# =============================================================================
# Sample Definitions
# =============================================================================

std_ids = ['NdF3']

#TODO ensure that you can modify all the necessary parameters from here instead of having to open the stdanrd ledger.json file in the standard folder
# els_to_use_for_mean_PB_calc: elements to include as standards; default is ["all"] if sample is "bulk", None otherwise. (all elements/lines). Use ["Fe", "O"] for specific elements, []/["none"] for none.

# =============================================================================

update_std_library = True
# If True, writes fitted corrected PB reference values into EDS_Stds_{voltage}keV.json.
# If False, acquisition/fitting still runs and CSV outputs are saved, but the standards
# dictionary JSON is left unchanged.

# Substrate elements (may depend on target_Xsp_counts)
els_substrate = None#['C', 'O', 'Al']  # Elements from substrate that may be present in the spectrum. Loads from legder if None.

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.batch_fit_experimental_stds import batch_fit_experimental_stds

exp_std_maker = batch_fit_experimental_stds(
    std_ids=std_ids,
    spectrum_lims=spectrum_lims,
    use_instrument_background=use_instrument_background,
    min_bckgrnd_cnts=min_bckgrnd_cnts,
    update_std_library = update_std_library,
    els_substrate=els_substrate,
    output_filename_suffix=output_filename_suffix,
    verbose=True,
    stds_path = exp_std_dir
)