#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch quantification and analysis of X-ray spectra for known precursor mixes.

Used to assess the extent of intermixing of precurors prior to a reaction.
See example at:
    L. N. Walters et al., Synthetic Accessibility and Sodium Ion Conductivity of the Na 8– x A x P 2 O 9 (NAP)
    High-Temperature Sodium Superionic Conductor Framework, Chem. Mater. 37, 6807 (2025).

This script provides automated batch quantification and (optionally) clustering/statistical
analysis of acquired X-ray spectra for multiple powder mixes. It is robust to missing files or
errors in individual samples, making it suitable for unattended batch processing.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

from autoemxsp.runners import batch_quantify_and_analyze

# =============================================================================
# Samples
# =============================================================================  
sample_IDs = ['known_powder_mixture_example']

results_path = None # Looks in default Results folder if left unspecified

# =============================================================================
# Options
# =============================================================================

max_analytical_error = 5 # w%
min_bckgrnd_cnts = 5

run_analysis = True

quantify_only_unquantified_spectra = False
interrupt_fits_bad_spectra = True
is_known_precursor_mixture = True

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
    run_analysis=run_analysis,
    quantify_only_unquantified_spectra=quantify_only_unquantified_spectra,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
    is_known_precursor_mixture = is_known_precursor_mixture
)