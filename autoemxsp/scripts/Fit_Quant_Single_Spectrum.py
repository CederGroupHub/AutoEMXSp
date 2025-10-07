#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting and quantification of a single X-ray spectrum.

For spectrum-level analysis of fitting and quantification performance.

Run this module directly to process a spectrum with the defined configuration options

Notes
-----
- Only the `sample_ID` and 'spectrum_ID' are required if acquisition output is saved in the default Results directory;
  otherwise, specify `results_path`.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

from autoemxsp.runners import fit_and_quantify_spectrum

# =============================================================================
# Sample and spectrum to process
# =============================================================================
sample_ID = 'Wulfenite_example'

spectrum_ID = 4  # Value reported in 'Spectrum #' column in Data.csv

results_path = None # Looks in default Results folder if left unspecified

# =============================================================================
# Options
# =============================================================================
is_particle = True
is_standard = False
quantify_plot = False
plot_signal = True
zoom_plot = False
line_to_plot = 'O_Ka'
fit_tol = 1e-4


max_undetectable_w_fr = 0
use_instrument_background = False
force_single_iteration = False
interrupt_fits_bad_spectra = False

# Params loaded from configuration file when left unspecified
spectrum_lims = None # (80, 1100)
els_substrate = None # ['C', 'O', 'Al']

quantifier = fit_and_quantify_spectrum(
    sample_ID=sample_ID,
    spectrum_ID=spectrum_ID,
    is_standard = is_standard,
    spectrum_lims = spectrum_lims,
    use_instrument_background=use_instrument_background,
    quantify_plot=quantify_plot,
    plot_signal=plot_signal,
    zoom_plot=zoom_plot,
    line_to_plot=line_to_plot,
    els_substrate=els_substrate,
    fit_tol=fit_tol,
    is_particle=is_particle,
    max_undetectable_w_fr=max_undetectable_w_fr,
    force_single_iteration=force_single_iteration,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
    results_path = results_path
)

#%% Optionally print atomic and mass fractions of reference compound
# from autoemxsp.tools.utils import print_element_fractions_table

# This uses the elements from the loaded sample config
# print_element_fractions_table('ZnF2')  # Or use: print_element_fractions_table(''.join(sample_cfg.elements))