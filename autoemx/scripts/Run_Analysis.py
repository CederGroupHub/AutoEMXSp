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

Citation
--------
If you use this package, please cite:
A. Giunto et al., Accurate SEM-EDS Quantification, Automation, and Machine Learning Enable
High-Throughput Compositional Characterization of Powders,
Nature Communications 17, 9735 (2026).
https://doi.org/10.1038/s41467-026-76633-x
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
# If 'eps' is not set, it defaults to 0.05 in Euclidean geometry and 0.3 in Aitchison geometry.
# Unspecified keys keep their existing/default values. Ignored for k-means.
dbscan_params: dict | None = None  # e.g. {'eps': 0.04, 'min_samples': 4}

# Clustering geometry: 'euclidean', 'aitchison' or 'auto'. Uses the saved value if set to None.
# 'aitchison' clusters CLR-transformed compositions (log-ratios), which weighs relative
# differences in minor elements more fairly. Applies to both k-means and DBSCAN; with DBSCAN,
# 'eps' is then in CLR units (default 0.3 if unset). Reported centroids remain in fractions.
# 'auto' uses Euclidean when at most 2 elements are present or when a major element is near zero
# (< 1%) in at least 10% of the spectra, and Aitchison otherwise.
clustering_geometry: str | None = None

# Aitchison parameters (only used when clustering_geometry is 'aitchison' or 'auto').
# Recognized keys: 'detection_limit_percent' (float > 0, default 0.5). Before the log-ratio transform,
#   exact zeros, and all values below this limit for trace elements (median < 2x the limit), are set to
#   0.65 x the limit, so noisy trace elements do not dominate the clustering.
#   'auto_near_zero_percent' (default 1.0) and 'auto_max_near_zero_fraction' (default 0.10): thresholds
#   used by 'auto' to decide when major elements are too often near zero for log-ratios.
aitchison_params: dict | None = None  # e.g. {'detection_limit_percent': 0.5}

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
min_bckgrnd_cnts = None # Minimum background counts under the reference peaks (quant_flag 8). None keeps the saved value. Applied with no need to requantify
quant_flags_accepted = [0, -1] #8 #, 4, 5, 6, 7, 8]  # See docs: Quantification flags (quant_flag).

# =============================================================================
# Plotting options
# =============================================================================
ref_formulae = None # List of candidate compositions. If the first entry is "" or None, the rest are appended to the 
                    # list loaded from Comp_analysis_configs.json; otherwise, the provided list replaces it.
                    # Uses values loaded from Comp_analysis_configs.json if ref_formulae = None.
els_excluded_clust_plot = None # List of elements to exclude from the 3D clustering plot. Uses default values if variable is set to None
els_to_plot = None # List of elements forced onto the clustering plot axes (opposite of els_excluded_clust_plot).
                   # With 2 or 3 entries, those become the axes. With 1 entry, remaining axes are filled from
                   # detectable elements. Must yield 2 or 3 axes. Uses defaults if None.
plot_custom_plots = False # If True, auto-creates sample_folder/custom_plot.py (first run) and uses it for plotting.
show_unused_compositions_cluster_plot = True
intermixing_violin_plots = None # Violin plots of the intermixing of the reported mixtures. None keeps the saved value (by default, plotted for known phase mixtures)
output_filename_suffix = ''
activate_new_analysis = "ask" # If these settings create a new analysis, whether it becomes the active one (settings reused later, shown in the GUI): "ask", True or False.
                              # The new analysis is saved either way; change the active one with Set_Active_Analysis.py

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
    clustering_geometry = clustering_geometry,
    aitchison_params = aitchison_params,
    clustering_features = clustering_features,
    els_excluded_clust_plot=els_excluded_clust_plot,
    els_to_plot=els_to_plot,
    k_finding_method = k_finding_method,
    do_matrix_decomposition = do_matrix_decomposition,
    max_analytical_error_percent=max_analytical_error_percent,
    min_bckgrnd_cnts=min_bckgrnd_cnts,
    quant_flags_accepted=quant_flags_accepted,
    plot_custom_plots=plot_custom_plots,
    show_unused_compositions_cluster_plot=show_unused_compositions_cluster_plot,
    intermixing_violin_plots=intermixing_violin_plots,
    output_filename_suffix=output_filename_suffix,
    activate_new_analysis=activate_new_analysis,
)