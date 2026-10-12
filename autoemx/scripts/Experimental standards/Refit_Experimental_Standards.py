#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refit Experimental Standards

This script refits the X-ray spectra of experimental standards (i.e., samples of known composition) that were
already acquired, e.g. with Run_Experimental_Standard_Collection.py, to recompute their reference values of
peak-to-background (P/B) ratios. Use it after changing fitting settings, the detector model, or the settings of the
standards below, without measuring the standards again.

Typical usage:
    - List the IDs of the standards in 'std_ids': their folders in exp_std_dir (default: ./Std_measurements).
    - Adjust the settings of the standards as needed (None keeps the values saved in each standard's ledger).
    - Set update_std_library = True to write the new reference P/B values into the standards file of the microscope.
    - Run the script.

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
# =============================================================================
# Sample Definitions
# =============================================================================

std_ids = ['Al']

# =============================================================================
# Standard settings. None keeps the value saved in each standard's ledger; values given apply to this run only.
# =============================================================================
els_to_use_for_mean_PB_calc = None  # Elements whose P/B ratios are used for the mean reference values: e.g. ["Fe", "O"]; ["all"] for all, [] or ["none"] for none
min_acceptable_PB_ratio = None  # Minimum P/B ratio for a peak to be accepted as a standard (cnts/(cnts*keV^-1))
quant_flags_accepted = None  # Quantification flags of the spectra used for the standard, e.g. [0]
generate_separate_std_dict = None  # If True, updates a separate standards file in the project folder instead of the main one

# =============================================================================

update_std_library = False
# If True, writes fitted corrected PB reference values into EDS_Stds_{voltage}keV.json, in the local calibration
# folder of the microscope (~/.autoemx/calibrations/<microscope_ID>/).
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
    els_to_use_for_mean_PB_calc=els_to_use_for_mean_PB_calc,
    min_acceptable_PB_ratio=min_acceptable_PB_ratio,
    quant_flags_accepted=quant_flags_accepted,
    generate_separate_std_dict=generate_separate_std_dict,
    output_filename_suffix=output_filename_suffix,
    verbose=True,
    stds_path = exp_std_dir
)