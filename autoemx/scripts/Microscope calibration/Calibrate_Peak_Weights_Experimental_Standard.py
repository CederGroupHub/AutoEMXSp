#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate relative X-ray line weights for one experimental standard.

This script is a thin configuration wrapper around:
`autoemx.runners.calibrate_peak_weights_experimental_standard`.

Created on Thu Jun 4 2026

@author: Andrea

Outputs (in the standard folder):
- *_peak_weights_areas_table_*.csv: line-by-row table with Meas 1..N, Mean, Stdev, Er (%).
- *_peak_weights_ratios_table_*.csv: same layout normalized to an auto-selected reference peak.
- *_peak_weights_summary_*.csv: correlations and sum-locked diagnostics.
"""

from autoemx.runners.calibrate_peak_weights_experimental_standard import (
    calibrate_peak_weights_experimental_standard,
)


# =============================================================================
# Standard and output configuration
# =============================================================================
std_ID = "PrPO4"  # Sample ID of the standard to calibrate, matching the ledger entry
results_path = None  # Uses cwd/Results if None

# Optional subset of spectra. Use None to fit every spectrum in the ledger.
spectrum_IDs_to_fit = None  # Example: [1, 2, 3, 4, 5]


# =============================================================================
# Calibrations of peak weights
# =============================================================================
fitted_extra_peaks = ['Pr' + '_' + line for line in ['Ma1', 'Mg', 'M2N4', 'Mb', 'Mz1', 'Mz2']] # ensure to include original reference peak(s)
# ['W' + '_' + line for line in ['La1', 'Lb4', 'Ln', 'Lb1', 'Lb6', 'L3N2', 'Lb3', 'L3N3', 'Lb2', 'Lb7', 'L3O2', 'Lu', 'Lb10', 'Lb5']]
# ['Mn' + '_' + line for line in ['La1', 'Ll', 'Ln', 'Lb1', 'Lg5', 'Lb3', 'Lb4']]


# =============================================================================
# Fitting options
# =============================================================================
spectrum_lims = None
els_substrate = None
is_particle = False
use_instrument_background = False
fit_tol = 1e-4


# =============================================================================
# Correlation signaling options
# =============================================================================
abs_corr_threshold = 0.90
sum_locked_neg_corr_threshold = -0.90
sum_cv_threshold = 0.12
n_sig_figs = 4


# =============================================================================
# Run
# =============================================================================
calibrate_peak_weights_experimental_standard(
    std_ID=std_ID,
    free_area_el_lines=fitted_extra_peaks,
    spectrum_IDs_to_fit=spectrum_IDs_to_fit,
    results_path=results_path,
    spectrum_lims=spectrum_lims,
    els_substrate=els_substrate,
    is_particle=is_particle,
    use_instrument_background=use_instrument_background,
    fit_tol=fit_tol,
    abs_corr_threshold=abs_corr_threshold,
    sum_locked_neg_corr_threshold=sum_locked_neg_corr_threshold,
    sum_cv_threshold=sum_cv_threshold,
    n_sig_figs=n_sig_figs,
)
