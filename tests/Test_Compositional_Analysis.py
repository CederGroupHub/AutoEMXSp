#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-sample clustering and analysis of X-ray spectra.

This script loads configurations and acquired X-ray spectra for a single sample,
performs clustering/statistical analysis, and prints results.

Run this file directly to analyze the specified sample.

Notes
-----
- Requires `sample_ID` (and optionally `results_path` if not using the default directory).
- Designed to be robust and flexible for both batch and single-sample workflows.


Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

from autoemxsp.runners import analyze_sample

# =============================================================================
# Initializations - Uses default values if variable is set to None
# =============================================================================
els_excluded_clust_plot = None
ref_formulae = None
clustering_features = None
k_finding_method = None

# =============================================================================
# Examples
# =============================================================================
# sample_ID = 'Wulfenite_example'
sample_ID = 'K-412_NISTstd_example'


results_path = None # Looks in default Results folder if left unspecified
# =============================================================================
# Clustering and Plotting options
# =============================================================================
k_forced = None

max_analytical_error_percent = 5 # w%
quant_flags_accepted = [0, -1] #8 #, 4, 5, 6, 7, 8]

plot_custom_plots = False
show_unused_compositions_cluster_plot = True

output_filename_suffix = ''

# =============================================================================
# Run
# =============================================================================
comp_analyzer = analyze_sample(
    sample_ID=sample_ID,
    results_path=results_path,
    ref_formulae=ref_formulae,
    k_forced = k_forced,
    els_excluded_clust_plot=els_excluded_clust_plot,
    k_finding_method = k_finding_method,
    max_analytical_error_percent=max_analytical_error_percent,
    quant_flags_accepted=quant_flags_accepted,
    plot_custom_plots=plot_custom_plots,
    show_unused_compositions_cluster_plot=show_unused_compositions_cluster_plot,
    output_filename_suffix=output_filename_suffix,
)