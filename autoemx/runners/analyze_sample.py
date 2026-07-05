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
    - Loads sample configuration and spectral data from ``ledger.json`` (primary source)
    - Falls back to ``Data.csv`` only when no ledger exists (one-time migration)
    - Performs clustering/statistical analysis
    - Prints summary results

Notes 
-----
- Requires `sample_ID` (and optionally `results_path` if not using the default directory).
- Designed to be robust and flexible for both batch and single-sample workflows.

Typical usage:
    - Edit the `sample_ID` and options in the script, or
    - Import and call `analyze_sample()` with your own arguments.
    

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import os
import time
import logging
from typing import Optional, List

import autoemx.utils as utils
import autoemx.utils.constants as cnst
from autoemx.utils.plotting_helpers import (
    ensure_custom_plot_file,
    refresh_custom_plot_template_file,
)
from autoemx.config import config_classes_dict, load_sample_ledger
from autoemx.config.ledger_schemas import ClusteringConfig, DBSCANParams
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["analyze_sample", "refresh_custom_plot_template"]


def refresh_custom_plot_template(
    sample_ID: str,
    results_path: Optional[str] = None,
    overwrite: bool = True,
) -> str:
    """(Re)create sample-local custom_plot.py from the packaged template."""
    custom_plot_file, was_written = refresh_custom_plot_template_file(
        sample_ID=sample_ID,
        results_path=results_path,
        overwrite=overwrite,
    )
    if was_written:
        logging.info("Custom plot template written to: %s", custom_plot_file)
    else:
        logging.info("Custom plot template already exists and was not overwritten: %s", custom_plot_file)

    return custom_plot_file


def _ensure_custom_plot_file(sample_dir: str, plot_cfg) -> None:
    """Create a sample-local custom plot template and store its path in plot config."""
    if not plot_cfg.use_custom_plots:
        return

    custom_plot_file, was_written = ensure_custom_plot_file(
        sample_dir=sample_dir,
        custom_plot_file=plot_cfg.custom_plot_file,
    )
    if was_written:
        logging.info("Created custom plot template: %s", custom_plot_file)

    plot_cfg.custom_plot_file = custom_plot_file

def analyze_sample(
    sample_ID: str,
    results_path: Optional[str] = None,
    output_filename_suffix: str = "",
    ref_formulae: Optional[List[str]] = None,
    els_excluded_clust_plot: Optional[List[str]] = None,
    clustering_features: Optional[str] = None,
    clustering_method: Optional[str] = None,
    dbscan_params: Optional[dict] = None,
    k_finding_method: Optional[str] = None,
    k_forced: Optional[int] = None,
    do_matrix_decomposition: bool = True,
    max_analytical_error_percent: float = 5,
    quant_flags_accepted: Optional[List[int]] = None,
    show_plots: bool = True,
    plot_custom_plots: bool = False,
    show_unused_compositions_cluster_plot: bool = True,
) -> Optional[EMXSp_Composition_Analyzer]:
    """
    Run clustering and analysis for a single sample.
    ----------
    sample_ID : str
        Sample identifier.
    results_path : str, optional
        Directory where results are loaded and stored. If None, defaults to autoemx/Results
    output_filename_suffix : str, optional
        Suffix for output files.
    ref_formulae : list of str, optional
        Reference formulae for clustering. If the first entry is "" or None, the rest are appended to the 
        list loaded from Comp_analysis_configs.json; otherwise, the provided list replaces it.
    els_excluded_clust_plot : list of str, optional
        Elements to exclude from cluster plot.
    clustering_features : list of str, optional
        Features to use for clustering.
    clustering_method : str, optional
        Clustering algorithm to use. One of ``"kmeans"`` or ``"dbscan"``. If None,
        the value stored in the sample's clustering config is kept.
    dbscan_params : dict, optional
        Overrides for DBSCAN parameters (only used when ``clustering_method="dbscan"``).
        Recognized keys: ``eps`` (float), ``min_samples`` (int), ``metric`` (str).
        Unspecified keys keep their existing/default values.
    k_finding_method : str, optional
        Method for determining optimal number of clusters. Set to "forced" if a value of 'k' is specified manually.
            Allowed methods are "silhouette", "calinski_harabasz", "elbow".
    k_forced : int, optional
        Forced number of clusters.
    do_matrix_decomposition : bool, optional
        Whether to compute matrix decomposition for intermixed phases. Slow if many candidate phases are provided. Default: True..
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
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)
        
    utils.print_double_separator()
    logging.info(f"Sample '{sample_ID}'")
    
    sample_dir = utils.get_sample_dir(results_path, sample_ID)
    ledger_path = os.path.join(sample_dir, f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
    spectral_info_f_path = ledger_path
    try:
        ledger = load_sample_ledger(ledger_path)
        configs = {
            cnst.MICROSCOPE_CFG_KEY: ledger.configs.microscope_cfg,
            cnst.SAMPLE_CFG_KEY: ledger.configs.sample_cfg,
            cnst.MEASUREMENT_CFG_KEY: ledger.configs.measurement_cfg,
            cnst.SAMPLESUBSTRATE_CFG_KEY: ledger.configs.sample_substrate_cfg,
            cnst.PLOT_CFG_KEY: ledger.configs.plot_cfg,
        }
        if ledger.quantifications:
            active_quant_id = ledger.active_quant
            active_quant_config = next(
                (
                    quant_config
                    for quant_config in ledger.quantifications
                    if quant_config.quantification_id == active_quant_id
                ),
                ledger.quantifications[-1],
            )
            configs[cnst.QUANTIFICATION_CFG_KEY] = config_classes_dict[cnst.QUANTIFICATION_CFG_KEY](
                **active_quant_config.options
            )
            active_clustering_analysis = active_quant_config.get_active_clustering_analysis()
            active_clustering_config = (
                active_clustering_analysis.config if active_clustering_analysis is not None else None
            )
            if active_clustering_config is not None:
                configs[cnst.CLUSTERING_CFG_KEY] = active_clustering_config
        else:
            configs[cnst.QUANTIFICATION_CFG_KEY] = config_classes_dict[cnst.QUANTIFICATION_CFG_KEY]()
        if ledger.configs.measurement_cfg.powder_meas_cfg is not None:
            configs[cnst.POWDER_MEASUREMENT_CFG_KEY] = ledger.configs.measurement_cfg.powder_meas_cfg
        if ledger.configs.measurement_cfg.bulk_meas_cfg is not None:
            configs[cnst.BULK_MEASUREMENT_CFG_KEY] = ledger.configs.measurement_cfg.bulk_meas_cfg
        metadata = {}
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
        clustering_cfg      = configs.get(cnst.CLUSTERING_CFG_KEY)
        plot_cfg            = configs[cnst.PLOT_CFG_KEY]
        powder_meas_cfg     = configs.get(cnst.POWDER_MEASUREMENT_CFG_KEY, None)  # Optional
        bulk_meas_cfg     = configs.get(cnst.BULK_MEASUREMENT_CFG_KEY, None)  # Optional
    except KeyError as e:
        logging.error(f"Missing configuration '{e.args[0]}' in {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return

    if clustering_cfg is None:
        clustering_cfg = ClusteringConfig()
    if quant_cfg is None:
        quant_cfg = config_classes_dict[cnst.QUANTIFICATION_CFG_KEY]()

    # --- Modify Clustering Configuration
    forced_key = "forced"
    allowed_k_finding_methods = ("silhouette", "calinski_harabasz", "elbow", forced_key)
    if quant_flags_accepted is not None:
        clustering_cfg.quant_flags_accepted = quant_flags_accepted
    clustering_cfg.max_analytical_error_percent = max_analytical_error_percent
    if ref_formulae is not None:
        if ref_formulae and (ref_formulae[0] == "" or ref_formulae[0] is None):
            # Append mode: skip the first empty entry
            clustering_cfg.ref_formulae.extend(ref_formulae[1:])
        else:
            # Replace mode
            clustering_cfg.ref_formulae = ref_formulae
    if clustering_features is not None:
        clustering_cfg.features = clustering_features
    if clustering_method is not None:
        clustering_cfg.method = clustering_method
    if dbscan_params is not None:
        # Merge overrides onto existing DBSCAN params and re-validate (model_copy(update=...)
        # bypasses field validators in pydantic v2).
        merged_dbscan = {**clustering_cfg.dbscan.model_dump(), **dbscan_params}
        clustering_cfg.dbscan = DBSCANParams.model_validate(merged_dbscan)
    if isinstance(k_forced, int):
        # Forces the k to be the provided number of clusters
        clustering_cfg.k_forced = k_forced
        clustering_cfg.k_finding_method = forced_key
    elif k_finding_method == forced_key:
        raise ValueError(
            f"'k_finding_method' must be one of {allowed_k_finding_methods}, "
            f"but not {forced_key}, if 'k_forced' is set to None"
        )
    elif k_finding_method is not None:
        # If k_forced is None and a method is specified, force recomputation of k in each run.
        clustering_cfg.k_forced = None
        clustering_cfg.k_finding_method = k_finding_method
    else:
        # If a finding method is not specified and k_forced is None, simply loads the default values from clustering_cfg
        pass
    
    if do_matrix_decomposition is not None:
        clustering_cfg.do_matrix_decomposition = do_matrix_decomposition
    
    # --- Modify Plot Configuration
    plot_cfg.show_plots = show_plots # show plots by default, but can be turned off for batch processing
    plot_cfg.show_unused_comps_clust = show_unused_compositions_cluster_plot
    plot_cfg.use_custom_plots = plot_custom_plots
    if els_excluded_clust_plot is not None:
        plot_cfg.els_excluded_clust_plot = els_excluded_clust_plot
    _ensure_custom_plot_file(sample_dir, plot_cfg)

    # Spectral source resolution is delegated to analyser._load_or_create_ledger():
    # --- Run Composition Analysis or Spectral Acquisition
    comp_analyzer = EMXSp_Composition_Analyzer(
        microscope_cfg=microscope_cfg,
        sample_id=sample_ID,
        sample_cfg=sample_cfg,
        measurement_cfg=measurement_cfg,
        sample_substrate_cfg=sample_substrate_cfg,
        quant_cfg=quant_cfg,
        initial_clustering_cfg=clustering_cfg,
        powder_meas_cfg=powder_meas_cfg,
        bulk_meas_cfg=bulk_meas_cfg,
        plot_cfg=plot_cfg,
        is_acquisition=False,
        development_mode=False,
        output_filename_suffix=output_filename_suffix,
        verbose=True,
        results_dir=sample_dir
    )

    # analyse_data calls _sync_in_memory_spectra_from_ledger / _load_or_create_ledger
    # and always hydrates spectra + quantification records from ledger-managed sources.

    has_ledger = os.path.exists(ledger_path)
    source_label = "Data.csv (first-run migration)" if (not has_ledger) else "ledger.json"
    logging.info(f"Running analysis for '{sample_ID}' (source: {source_label}).") 

    # Perform analysis and print results
    try:
        analysis_successful, _, _ = comp_analyzer.analyse_data(
            max_analytical_error_percent,
            k=comp_analyzer.clustering_cfg.k_forced if comp_analyzer.clustering_cfg.k_finding_method == forced_key else None,
        )
    except Exception as e:
        logging.exception(f'Error during clustering analysis for {sample_ID}: {e}')
        return

    total_process_time = (time.time() - sample_processing_time_start)
    
    if analysis_successful:
        comp_analyzer.print_results()
        utils.print_single_separator()
        logging.info(f"Sample '{sample_ID}' successfully analysed in {total_process_time:.1f} sec.")
    else:
        utils.print_single_separator()
        logging.info(f"Analysis was not successful for '{sample_ID}'.")
    
    return comp_analyzer