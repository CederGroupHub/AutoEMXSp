#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-sample clustering and analysis of X-ray spectra.

This module loads configurations and acquired X-ray spectra for a single sample,
performs clustering/statistical analysis, and prints results. It is suitable for
both interactive use and integration into larger workflows.

Import this module in your own code and call the
`analyze_sample()` function, passing the sample ID (and optional arguments)
to perform analysis programmatically.

Workflow:
    - Loads sample configuration from `Spectra_collection_info.json`
    - Loads acquired spectral data from `Data.csv`
    - Performs clustering/statistical analysis
    - Prints summary results

Notes
-----
- Requires `sample_ID` (and optionally `results_path` if not using the default directory).
- Designed to be robust and flexible for both batch and single-sample workflows.

Typical usage:
    - Edit the `sample_ID` and options in the script, or
    - Import and call `analyze_sample()` with your own arguments.
    
Parameters
----------
sample_ID : str
    Sample identifier.
results_path : str, optional
    Directory where results are stored. If None, sets to default autoemxsp/Results
output_filename_suffix : str, optional
    Suffix for output files.
ref_formulae : list of str, optional
    Reference formulae for clustering.
els_excluded_clust_plot : list of str, optional
    Elements to exclude from cluster plot.
clustering_features : list of str, optional
    Features to use for clustering.
k_finding_method : str, optional
    Method for determining optimal number of clusters.
k_forced : int, optional
    Forced number of clusters.
max_analytical_error_percent : float, optional
    Maximum analytical error allowed for clustering.
quant_flags_accepted : list of int, optional
    Accepted quantification flags.
plot_custom_plots : bool, optional
    Whether to use custom plots.
show_unused_compositions_cluster_plot : bool, optional
    Whether to show unused compositions in cluster plot.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import os
import time
import logging
from typing import Optional, List

from autoemxsp.tools.utils import (
    print_single_separator,
    print_double_separator,
    get_sample_dir,
    load_configurations_from_json,
    extract_spectral_data,
)
import autoemxsp.tools.constants as cnst
from autoemxsp.tools.config_classes import config_classes_dict
from autoemxsp.core.EMXSp_composition_analyser import EMXSp_Composition_Analyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def analyze_sample(
    sample_ID: str,
    results_path: str = None,
    output_filename_suffix: str = "",
    ref_formulae: Optional[List[str]] = None,
    els_excluded_clust_plot: Optional[List[str]] = None,
    clustering_features: Optional[List[str]] = None,
    k_finding_method: Optional[str] = None,
    k_forced: Optional[int] = None,
    max_analytical_error_percent: float = 5,
    quant_flags_accepted: Optional[List[int]] = None,
    plot_custom_plots: bool = False,
    show_unused_compositions_cluster_plot: bool = True,
) -> None:
    """
    Run clustering and analysis for a single sample.

    Parameters
    ----------
    sample_ID : str
        Sample identifier.
    results_path : str, optional
        Directory where results are loaded and stored. If None, defaults to autoemxsp/Results
    output_filename_suffix : str, optional
        Suffix for output files.
    ref_formulae : list of str, optional
        Reference formulae for clustering.
    els_excluded_clust_plot : list of str, optional
        Elements to exclude from cluster plot.
    clustering_features : list of str, optional
        Features to use for clustering.
    k_finding_method : str, optional
        Method for determining optimal number of clusters.
    k_forced : int, optional
        Forced number of clusters.
    max_analytical_error_percent : float, optional
        Maximum analytical error allowed for clustering.
    quant_flags_accepted : list of int, optional
        Accepted quantification flags.
    plot_custom_plots : bool, optional
        Whether to use custom plots.
    show_unused_compositions_cluster_plot : bool, optional
        Whether to show unused compositions in cluster plot.
        
    Returns
    -------
    comp_analyzer : EMXSp_Composition_Analyzer
        The composition analysis object containing the results and methods for further analysis.
    """
    if results_path is None:
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        results_path = os.path.join(parent_dir, cnst.RESULTS_DIR)
        
    print_double_separator()
    logging.info(f"Sample '{sample_ID}'")
    
    sample_dir = get_sample_dir(results_path, sample_ID)
    spectral_info_f_path = os.path.join(sample_dir, f'{cnst.ACQUISITION_INFO_FILENAME}.json')
    try:
        configs, metadata = load_configurations_from_json(spectral_info_f_path, config_classes_dict)
    except FileNotFoundError:
        logging.error(f"Could not find {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return
    except Exception as e:
        logging.error(f"Error loading {spectral_info_f_path}. Skipping sample '{sample_ID}': {e}")
        return

    sample_processing_time_start = time.time()

    # Retrieve configuration objects for this sample
    try:
        microscope_cfg      = configs[cnst.MICROSCOPE_CFG_KEY]
        sample_cfg          = configs[cnst.SAMPLE_CFG_KEY]
        measurement_cfg     = configs[cnst.MEASUREMENT_CFG_KEY]
        sample_substrate_cfg= configs[cnst.SAMPLESUBSTRATE_CFG_KEY]
        quant_cfg           = configs[cnst.QUANTIFICATION_CFG_KEY]
        clustering_cfg      = configs[cnst.CLUSTERING_CFG_KEY]
        plot_cfg            = configs[cnst.PLOT_CFG_KEY]
        powder_meas_cfg     = configs.get(cnst.POWDER_MEASUREMENT_CFG_KEY, None)  # Optional
        bulk_meas_cfg     = configs.get(cnst.BULK_MEASUREMENT_CFG_KEY, None)  # Optional
    except KeyError as e:
        logging.error(f"Missing configuration '{e.args[0]}' in {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return
    
    # --- Modify Clustering Configuration
    if quant_flags_accepted is not None:
        clustering_cfg.quant_flags_accepted = quant_flags_accepted
    clustering_cfg.max_analytical_error_percent = max_analytical_error_percent
    if ref_formulae is not None:
        clustering_cfg.ref_formulae = ref_formulae
    if clustering_features is not None:
        clustering_cfg.features = clustering_features
    if isinstance(k_forced, int):
        clustering_cfg.k = k_forced
        clustering_cfg.k_finding_method = 'forced'
    elif k_finding_method == 'forced':
        raise ValueError(f"'k_finding_method' must be one of {clustering_cfg.ALLOWED_K_FINDING_METHODS}, but not 'forced', if 'k_forced' is set to None")
    elif k_finding_method is not None:
        clustering_cfg.k_finding_method = k_finding_method

    # --- Modify Plot Configuration
    plot_cfg.show_unused_comps_clust = show_unused_compositions_cluster_plot
    plot_cfg.use_custom_plots = plot_custom_plots
    if els_excluded_clust_plot is not None:
        plot_cfg.els_excluded_clust_plot = els_excluded_clust_plot

    # Load 'Data.csv' into a DataFrame
    data_path = os.path.join(sample_dir, f'{cnst.DATA_FILENAME}.csv')
    try:
        spectra_quant, spectral_data, sp_coords, _ = extract_spectral_data(data_path)
    except Exception as e:
        logging.error(f"Could not load spectral data for '{sample_ID}': {e}")
        return
    
    if spectra_quant is None:
        logging.error(f"No quantification data found in {data_path}")
        return
    
    # --- Run Composition Analysis or Spectral Acquisition
    comp_analyzer = EMXSp_Composition_Analyzer(
        microscope_cfg=microscope_cfg,
        sample_cfg=sample_cfg,
        measurement_cfg=measurement_cfg,
        sample_substrate_cfg=sample_substrate_cfg,
        quant_cfg=quant_cfg,
        clustering_cfg=clustering_cfg,
        powder_meas_cfg=powder_meas_cfg,
        bulk_meas_cfg=bulk_meas_cfg,
        plot_cfg=plot_cfg,
        is_acquisition=False,
        development_mode=False,
        output_filename_suffix=output_filename_suffix,
        verbose=True,
        results_dir=sample_dir
    )

    comp_analyzer.spectra_quant = spectra_quant
    comp_analyzer.sp_coords = sp_coords
    comp_analyzer.spectral_data = spectral_data

    # Perform analysis and print results
    try:
        analysis_successful, _, _ = comp_analyzer.analyse_data(max_analytical_error_percent, k = clustering_cfg.k)
    except Exception as e:
        logging.exception(f'Error during clustering analysis for {sample_ID}: {e}')
        return

    total_process_time = (time.time() - sample_processing_time_start)
    
    if analysis_successful:
        comp_analyzer.print_results()
        print_single_separator()
        logging.info(f"Sample '{sample_ID}' successfully analysed in {total_process_time:.1f} sec.")
    else:
        print_single_separator()
        logging.info(f"Analysis was not successful for '{sample_ID}'.")
    
    return comp_analyzer