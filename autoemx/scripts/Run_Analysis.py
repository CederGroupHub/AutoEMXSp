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

# =============================================================================
# Sample Definition
# =============================================================================
sample_ID = 'Wulfenite_example'

import os
results_dir = os.path.dirname(os.path.abspath(__file__)) # Default: save and load results in the same folder as this script. Set to None to use the current working directory, or replace with another path.

# =============================================================================
# Clustering options
# =============================================================================
clustering_features = None # 'w_fr', 'at_fr'. Uses default value if variable is set to None

# Clustering algorithm: 'kmeans' or 'dbscan'. Uses the saved value if set to None.
clustering_method: str | None = None

# DBSCAN parameters (only used when clustering_method == 'dbscan').
# Recognized keys: 'eps' (float), 'min_samples' (int), 'metric' (str).
# Unspecified keys keep their existing/default values. Ignored for k-means.
dbscan_params: dict | None = None  # e.g. {'eps': 0.04, 'min_samples': 4}

# Number of clusters to use, if manually specified (k-means only).
# If None, the number of clusters will be determined automatically.
k_forced: int | None = None  

# Method used to determine the number of clusters (see ClusteringConfig.ALLOWED_K_FINDING_METHODS).
# Only applied if `k_forced` is None. Forces re-computation of the optimal k value. (k-means only)
k_finding_method: str | None = None  

# Behavior:
# - If both `k_finding_method` and `k_forced` are None, clustering configurations
#   are loaded directly from the saved `Comp_analysis_configs.json` file.

# Whether to compute matrix decomposition for intermixed phases. Slow if many candidate phases are provided.
do_matrix_decomposition = True

# =============================================================================
# Spectral Filtering options
# =============================================================================
max_analytical_error_percent = 5 # w%
quant_flags_accepted = [0, -1] #8 #, 4, 5, 6, 7, 8]  # See docs: Quantification flags (quant_flag).

# =============================================================================
# Plotting options
# =============================================================================
ref_formulae = None # List of candidate compositions. If the first entry is "" or None, the rest are appended to the 
                    # list loaded from Comp_analysis_configs.json; otherwise, the provided list replaces it.
                    # Uses values loaded from Comp_analysis_configs.json if ref_formulae = None.
els_excluded_clust_plot = None # List of elements to exclude from the 3D clustering plot. Uses default values if variable is set to None
plot_custom_plots = False # If True, auto-creates sample_folder/custom_plot.py (first run) and uses it for plotting.
show_unused_compositions_cluster_plot = True
output_filename_suffix = ''

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.analyze_sample import analyze_sample

comp_analyzer = analyze_sample(
    sample_ID=sample_ID,
    results_path=results_dir,
    ref_formulae=ref_formulae,
    k_forced = k_forced,
    clustering_method = clustering_method,
    dbscan_params = dbscan_params,
    clustering_features = clustering_features,
    els_excluded_clust_plot=els_excluded_clust_plot,
    k_finding_method = k_finding_method,
    do_matrix_decomposition = do_matrix_decomposition,
    max_analytical_error_percent=max_analytical_error_percent,
    quant_flags_accepted=quant_flags_accepted,
    plot_custom_plots=plot_custom_plots,
    show_unused_compositions_cluster_plot=show_unused_compositions_cluster_plot,
    output_filename_suffix=output_filename_suffix,
)