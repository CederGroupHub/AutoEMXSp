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
# sample_ID = 'known_powder_mixture_example'

# =============================================================================
# Paper data (Download data from github repository:
#       https://github.com/CederGroupHub/AutoEMXSp/tree/main/paper_data)   
# =============================================================================
# --- Minerals
# sample_ID = 'Alamosite_mineral'
# sample_ID = 'Albite_mineral'
# sample_ID = 'Anhydrite_mineral'
# sample_ID = 'Anorthite_mineral'
# sample_ID = 'Benitoite_mineral'
# sample_ID = 'Bornite_mineral'
# sample_ID = 'Chalcopyrite_mineral'
# sample_ID = 'CoOlivine_mineral'
# sample_ID = 'FeOlivine_mineral'
# sample_ID = 'Fluorphlogopite_mineral'
# sample_ID = 'Jadeite_mineral'
# sample_ID = 'K-412_NISTstd_mineral'
# sample_ID = 'Labradorite_mineral'
# sample_ID = 'MnOlivine_mineral'
# sample_ID = 'Nepheline_mineral'
# sample_ID = 'Orthoclase_mineral'
# sample_ID = 'Rhodonite_mineral'
# sample_ID = 'ScPO4_mineral'
# sample_ID = 'Wulfenite_mineral'
# sample_ID = 'YIG_mineral'
# sample_ID = 'YPO4_mineral'

# --- Commercial precursors
# sample_ID = 'Al2O3_precursor'
# sample_ID = 'AlPO4_precursor'
# sample_ID = 'Ba3(PO4)2_precursor'
# sample_ID = 'Co3O4_precursor'
# sample_ID = 'CuO_precursor'
# sample_ID = 'Fe2O3_precursor'
# sample_ID = 'Ga2O3_precursor'
# sample_ID = 'GeO2_precursor'
# sample_ID = 'HfO2_precursor'
# sample_ID = 'In2O3_precursor'
# sample_ID = 'KCl_precursor'
# sample_ID = 'Li2WO4_precursor'
# sample_ID = 'LiCoPO4_precursor'
# sample_ID = 'LiNiCoMnO2_precursor'
# sample_ID = 'MgF2_precursor'
# sample_ID = 'MgO_precursor'
# sample_ID = 'MnO_precursor'
# sample_ID = 'MnO2_precursor'
# sample_ID = 'Mn2O3_precursor'
# sample_ID = 'MoO2_precursor'
# sample_ID = 'Na2MoO4_precursor'
# sample_ID = 'Na4P2O7_precursor'
# sample_ID = 'NaNO3_precursor'
# sample_ID = 'Ni(OH)2_precursor'
# sample_ID = 'NiO_precursor'
# sample_ID = 'PbO_precursor'
# sample_ID = 'Sb2O3_precursor'
# sample_ID = 'SiO2_precursor'
# sample_ID = 'SnO2_precursor'
# sample_ID = 'Ta2O5_precursor'
# sample_ID = 'TiN_precursor'
# sample_ID = 'TiO2_precursor'
# sample_ID = 'WO3_precursor'
# sample_ID = 'ZnF2_precursor'
# sample_ID = 'ZnO_precursor'
# sample_ID = 'ZrO2_precursor'

# sample_ID = 'MnO-Mn2O3-Fe2O3_mix'
# k_finding_method = 'calinski_harabasz'
# els_excluded_clust_plot = []

## Synthetic samples
# sample_ID = 'NaGe2(PO4)3'
# sample_ID = 'NaSn2(PO4)3'
# sample_ID = 'Na0.4Zr1.4Ta0.6(PO4)3'
# sample_ID = 'NaZrTi(PO4)3'
# sample_ID = 'NaTiSn(PO4)3'
# sample_ID = 'MnAgO2'
# sample_ID = 'CaCo(PO3)4'
# sample_ID = 'MgTi2NiO6'
# sample_ID = 'K4MgFe3(PO4)5'
# sample_ID = 'KNaTi2(PO5)2'
# sample_ID = 'Hf2Sb2Pb4O13'
# sample_ID = 'K2TiCr(PO4)3'
# sample_ID = 'MgCuP2O7'
# sample_ID = 'MgTi4(PO4)6'
# sample_ID = 'NaCaMgFe(SiO3)4'
# sample_ID = 'Bi2Fe4O9'
# sample_ID = 'Bi25FeO39'
# sample_ID = 'LaNbO4'

# --- Other
# sample_ID = 'Nepheline+LMNO mixture'
# sample_ID = 'NASICON synthetic mixture'


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