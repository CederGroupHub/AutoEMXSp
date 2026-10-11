#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting and quantification of a single X-ray spectrum from the ledger.

Example: Fit and quantify spectrum #1 from the Wulfenite sample.

This example demonstrates the ledger-native single-spectrum runner,
which loads spectra and configurations from the sample's ledger.json
and spectra/ folder, then quantifies a single requested spectrum.

Created on April 24, 2026

@author: Andrea
"""

import os

from autoemx.runners.fit_and_quantify_spectrum_from_ledger import fit_and_quantify_spectrum_from_ledger

# =============================================================================
# Sample and spectrum to process
# =============================================================================
sample_ID = 'Wulfenite_example'

spectrum_ID = 1  # Spectrum # from ledger metadata

# Use examples/Results as the results directory
results_path = os.path.join(os.path.dirname(__file__), 'Results')

# =============================================================================
# Fitting and quantification options
# =============================================================================
is_particle = True
is_standard = False
quantify = True
plot_signal = True
zoom_plot = False
line_to_plot = 'W_Ma'
fit_tol = 1e-4

max_undetectable_w_fr = 0
use_instrument_background = False
force_single_iteration = False
interrupt_fits_bad_spectra = False

# Params loaded from ledger configuration when left unspecified
spectrum_lims = None  # (80, 1100) — uses ledger value if None
els_substrate = None  # ['C', 'O', 'Al'] — uses ledger value if None

# =============================================================================
# Run
# =============================================================================
quantifier = fit_and_quantify_spectrum_from_ledger(
    sample_ID=sample_ID,
    spectrum_ID=spectrum_ID,
    results_path=results_path,
    is_standard=is_standard,
    spectrum_lims=spectrum_lims,
    use_instrument_background=use_instrument_background,
    quantify=quantify,
    plot_signal=plot_signal,
    zoom_plot=zoom_plot,
    line_to_plot=line_to_plot,
    els_substrate=els_substrate,
    fit_tol=fit_tol,
    is_particle=is_particle,
    max_undetectable_w_fr=max_undetectable_w_fr,
    force_single_iteration=force_single_iteration,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
)
