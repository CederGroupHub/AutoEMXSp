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
    - Loads sample configuration and spectral data from ``ledger.json``. Results saved in the
      legacy ``Data.csv`` format are converted to a ``ledger.json`` automatically
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
from typing import Optional, List, Union

import autoemx.utils as utils
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.utils.plotting_helpers import (
    ensure_custom_plot_file,
    refresh_custom_plot_template_file,
)
from autoemx.config import load_sample_ledger
from autoemx.config.ledger_io import ledger_runtime_configs
from autoemx.runners.set_active_analysis import set_active_analysis
from autoemx.config.ledger_schemas import (
    AitchisonParams,
    ClusterMergeParams,
    ClusteringConfig,
    DBSCANParams,
    MixtureParams,
)
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

def _ask_yes_no(question: str) -> bool:
    """Ask a yes/no question in the console (Enter = yes). Returns True if no answer can be given."""
    try:
        answer = input(question)
    except (EOFError, OSError):
        logging.info("No console input available: answering yes.")
        return True
    return answer.strip().lower() not in ("n", "no")


def _keep_or_revert_new_analysis(
    ledger_path: str,
    quantification_id: Optional[int],
    previous_index: Optional[int],
    n_previous_analyses: int,
    activate_new_analysis: Union[bool, str],
) -> None:
    """Decide (or ask) whether a clustering analysis created by this run stays active; else restore the previous one."""
    if quantification_id is None:
        return
    ledger = load_sample_ledger(ledger_path)
    quant_config = next((q for q in ledger.quantifications if q.quantification_id == quantification_id), None)
    if quant_config is None or len(quant_config.clustering_analyses) <= n_previous_analyses:
        return  # No new analysis: an existing one with the same settings was reused
    new_clustering_id = quant_config.clustering_analyses[-1].config.clustering_id
    activate = activate_new_analysis
    if activate == "ask":
        activate = _ask_yes_no(
            f"New clustering analysis {new_clustering_id} (quantification {quantification_id}). "
            "Make it the active analysis? [Y/n] "
        )
    if activate:
        return
    if previous_index is None:
        logging.info(f"Clustering analysis {new_clustering_id} stays active: no analysis was active before.")
        return
    previous_clustering_id = quant_config.clustering_analyses[previous_index].config.clustering_id
    set_active_analysis(quantification_id, previous_clustering_id, sample_path=os.path.dirname(ledger_path))
    logging.info(
        f"Clustering analysis {new_clustering_id} saved but not active: the active analysis is still clustering "
        f"{previous_clustering_id}. Use Set_Active_Analysis.py to change it."
    )


def analyze_sample(
    sample_ID: str,
    results_path: Optional[str] = None,
    output_filename_suffix: str = "",
    ref_formulae: Optional[List[str]] = None,
    els_excluded_clust_plot: Optional[List[str]] = None,
    els_to_plot: Optional[List[str]] = None,
    clustering_features: Optional[str] = None,
    clustering_method: Optional[str] = None,
    dbscan_params: Optional[dict] = None,
    clustering_geometry: Optional[str] = None,
    aitchison_params: Optional[dict] = None,
    k_finding_method: Optional[str] = None,
    k_forced: Optional[Union[int, bool]] = None,
    max_k: Optional[int] = None,
    auto_merge_clusters: Optional[bool] = None,
    cluster_merge_params: Optional[dict] = None,
    do_matrix_decomposition: bool = True,
    mixture_params: Optional[dict] = None,
    max_analytical_error_percent: float = dflt.max_analytical_error_percent,
    quant_flags_accepted: Optional[List[int]] = None,
    min_bckgrnd_cnts: Optional[float] = None,
    show_plots: bool = True,
    plot_custom_plots: bool = False,
    show_unused_compositions_cluster_plot: bool = True,
    plot_best_mixture: Optional[bool] = None,
    intermixing_violin_plots: Optional[bool] = None,
    show_legend_clustering: Optional[bool] = None,
    activate_new_analysis: Union[bool, str] = True,
) -> Optional[EMXSp_Composition_Analyzer]:
    """
    Run clustering and analysis for a single sample.

    Parameters
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
    els_to_plot : list of str, optional
        Elements forced onto the clustering plot axes (opposite of ``els_excluded_clust_plot``).
        When provided with 2 or 3 elements, those become the plot axes. A single element is
        kept and filled from remaining detectable elements until 2 or 3 axes are available.
        Must yield exactly 2 or 3 axes; unknown elements raise ``ValueError``.
    clustering_features : list of str, optional
        Features to use for clustering.
    clustering_method : str, optional
        Clustering algorithm to use. One of ``"kmeans"`` or ``"dbscan"``. If None,
        the value stored in the sample's clustering config is kept.
    dbscan_params : dict, optional
        Overrides for DBSCAN parameters (only used when ``clustering_method="dbscan"``).
        Recognized keys: ``eps`` (float), ``min_samples`` (int), ``metric`` (str).
        If ``eps`` is not set, it defaults to 0.05 in Euclidean and 0.3 in Aitchison geometry.
        Unspecified keys keep their existing/default values.
    clustering_geometry : str, optional
        Geometry in which clustering distances are computed. One of ``"euclidean"``,
        ``"aitchison"`` (clusters CLR-transformed compositions; applies to both k-means and
        DBSCAN, and DBSCAN ``eps`` is then in CLR units, default 0.3) or ``"auto"`` (Euclidean when at most
        two elements are present or a major element is near zero in many spectra, Aitchison
        otherwise; see ``AitchisonParams``). Reported centroids and cluster statistics remain
        in fraction space. If None, the saved value is kept.
    aitchison_params : dict, optional
        Overrides for Aitchison parameters (used when ``clustering_geometry`` is ``"aitchison"`` or ``"auto"``).
        Recognized keys: ``detection_limit_percent`` (float > 0, default 0.5). Before the log-ratio
        transform, exact zeros, and all values below the limit for trace elements (median below
        2x the limit), are set to 0.65 x the limit.
        ``auto_near_zero_percent`` (default 1.0) and ``auto_max_near_zero_fraction`` (default 0.10)
        set the near-zero criterion used by ``"auto"``.
    k_finding_method : str, optional
        Method for determining optimal number of clusters. Set to "forced" if a value of 'k' is specified manually.
            Allowed methods are "silhouette", "calinski_harabasz", "elbow".
    k_forced : int or bool, optional
        Controls the number of clusters:
            - ``int``: force clustering to use exactly this number of clusters.
            - ``False``: force recomputation of the optimal number of clusters,
              discarding any previously saved forced k. Uses ``k_finding_method``
              when provided, otherwise the saved (or default) finding method.
            - ``None`` (default): reuse the clustering settings saved in the ledger.
    max_k : int, optional
        Maximum number of clusters tested when k is determined automatically (k-means only).
        If None, the saved value is kept (default 6).
    auto_merge_clusters : bool, optional
        Whether to merge k-means clusters that are pieces of one continuous population, e.g. a mixture
        line cut into several clusters (only when k is chosen automatically; see ``ClusterMergeParams``).
        If None, the saved value is kept (True for new samples, False for samples analysed before this option).
    cluster_merge_params : dict, optional
        Overrides for the cluster-merging criteria (see ``ClusterMergeParams``). Recognized keys:
        ``dip_alpha``, ``max_gap_ratio``, ``dbscan_min_samples``, ``dbscan_eps_factor``.
        Unspecified keys keep their existing/default values.
    do_matrix_decomposition : bool, optional
        Whether to compute matrix decomposition for intermixed phases. Slow if many candidate phases are provided. Default: True..
    mixture_params : dict, optional
        Overrides for the decomposition of clusters into mixtures of candidate phases (see ``MixtureParams``).
        Main keys: ``max_n_phases`` (int >= 2, default 4) and ``max_recon_error`` (float, default 0.4).
        Pairs of candidate phases are always tested; combinations of more phases are only tested when no
        combination with fewer phases reconstructs the cluster with an error below ``max_recon_error``.
        Use ``{'max_n_phases': 2}`` to only consider binary mixtures. Other keys (reconstruction error and
        confidence definitions, single-phase and free-NMF criteria) are described in ``MixtureParams``.
        Unspecified keys keep their existing/default values.
    max_analytical_error_percent : float, optional
        Maximum analytical error allowed for clustering.
    quant_flags_accepted : list of int, optional
        Accepted quantification flags. See the user documentation page
        "Quantification flags (quant_flag)".
    min_bckgrnd_cnts : float, optional
        Minimum background counts under the reference peaks: spectra below are flagged (quant_flag 8), with no need to
        re-quantify (ClusteringConfig.min_bckgrnd_cnts). If None, the saved value is kept.
    plot_custom_plots : bool, optional
        Whether to use custom plots.
    show_unused_compositions_cluster_plot : bool, optional
        Whether to show unused compositions in cluster plot.
    plot_best_mixture : bool, optional
        Whether to save the plot of the top-ranked mixture of each cluster. If None, the saved value is kept.
    intermixing_violin_plots : bool, optional
        Whether to plot the intermixing of the reported mixtures as violin plots. If None, the saved value is kept
        (by default, they are plotted for samples quantified as known phase mixtures).
    show_legend_clustering : bool, optional
        Whether to show the legend in the clustering plot. If None, the saved value is kept.
    activate_new_analysis : bool or "ask", optional
        If the settings create a new clustering analysis, whether it becomes the active one (whose settings are
        reused by later analyses, and shown by default in the GUI). "ask" asks once the analysis is done (if no
        answer can be given, e.g. without a console, it becomes active). The new analysis is saved either way;
        the active one can be changed later with set_active_analysis. Default: True.
        
    Returns
    -------
    comp_analyzer : EMXSp_Composition_Analyzer
        The composition analysis object containing the results and methods for further analysis.
    """
    if results_path is None:
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)
        
    utils.print_double_separator(f"Sample '{sample_ID}'")
    
    sample_dir = utils.get_sample_dir(results_path, sample_ID)
    ledger_path = os.path.join(sample_dir, f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
    spectral_info_f_path = ledger_path
    try:
        ledger = load_sample_ledger(ledger_path)
        ledger_cfgs = ledger_runtime_configs(ledger)
    except Exception as e:
        logging.error(f"Error loading {spectral_info_f_path}. Skipping sample '{sample_ID}': {e}")
        return

    sample_processing_time_start = time.time()

    # Retrieve configuration objects for this sample
    clustering_cfg       = ledger_cfgs.clustering_cfg
    plot_cfg             = ledger_cfgs.plot_cfg


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
    if clustering_geometry is not None:
        clustering_cfg.geometry = ClusteringConfig.validate_geometry(clustering_geometry)
    if aitchison_params is not None:
        merged_aitchison = {**clustering_cfg.aitchison.model_dump(), **aitchison_params}
        clustering_cfg.aitchison = AitchisonParams.model_validate(merged_aitchison)
    if k_forced is False:
        # Force recomputation of the optimal number of clusters, discarding any
        # previously saved forced k. Prefer an explicitly provided finding method,
        # then the saved one, falling back to the default when it was 'forced'.
        clustering_cfg.k_forced = None
        if k_finding_method is not None and k_finding_method != forced_key:
            clustering_cfg.k_finding_method = k_finding_method
        elif clustering_cfg.k_finding_method == forced_key:
            clustering_cfg.k_finding_method = ClusteringConfig.model_fields["k_finding_method"].default
    elif k_forced is True:
        raise ValueError(
            "k_forced=True is not supported. Pass an integer to force a specific "
            "number of clusters, False to force recomputation of the optimal k, or "
            "None to reuse the saved clustering settings."
        )
    elif isinstance(k_forced, int):
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
    
    if max_k is not None:
        clustering_cfg.max_k = ClusteringConfig.validate_max_k(int(max_k))
    if min_bckgrnd_cnts is not None:
        clustering_cfg.min_bckgrnd_cnts = ClusteringConfig.validate_min_bckgrnd_cnts(float(min_bckgrnd_cnts))
    if do_matrix_decomposition is not None:
        clustering_cfg.do_matrix_decomposition = do_matrix_decomposition
    if auto_merge_clusters is not None:
        clustering_cfg.auto_merge_clusters = bool(auto_merge_clusters)
    if cluster_merge_params is not None:
        merged_merge = {**clustering_cfg.cluster_merge.model_dump(), **cluster_merge_params}
        clustering_cfg.cluster_merge = ClusterMergeParams.model_validate(merged_merge)
    if mixture_params is not None:
        merged_mixture = {**clustering_cfg.mixture.model_dump(), **mixture_params}
        clustering_cfg.mixture = MixtureParams.model_validate(merged_mixture)
    
    # --- Modify Plot Configuration
    plot_cfg.show_plots = show_plots # show plots by default, but can be turned off for batch processing
    plot_cfg.show_unused_comps_clust = show_unused_compositions_cluster_plot
    plot_cfg.use_custom_plots = plot_custom_plots
    if els_excluded_clust_plot is not None:
        plot_cfg.els_excluded_clust_plot = els_excluded_clust_plot
    if els_to_plot is not None:
        plot_cfg.els_to_plot = els_to_plot
    if plot_best_mixture is not None:
        plot_cfg.plot_best_mixture = bool(plot_best_mixture)
    if intermixing_violin_plots is not None:
        plot_cfg.intermixing_violin_plots = bool(intermixing_violin_plots)
    if show_legend_clustering is not None:
        plot_cfg.show_legend_clustering = bool(show_legend_clustering)
    _ensure_custom_plot_file(sample_dir, plot_cfg)

    # --- Run Composition Analysis
    comp_analyzer = EMXSp_Composition_Analyzer.from_ledger(
        ledger,
        ledger_cfgs,
        sample_id=sample_ID,
        output_filename_suffix=output_filename_suffix,
    )

    if activate_new_analysis not in (True, False, "ask"):
        raise ValueError(f"activate_new_analysis must be True, False or 'ask', got {activate_new_analysis!r}")
    # Active analysis before this run, to restore it if a new analysis is not made active
    active_quant_config = ledger_cfgs.active_quant_config
    previous_quantification_id = active_quant_config.quantification_id if active_quant_config is not None else None
    previous_index = active_quant_config.active_clustering_analysis_index if active_quant_config is not None else None
    n_previous_analyses = len(active_quant_config.clustering_analyses) if active_quant_config is not None else 0

    has_ledger = os.path.exists(ledger_path)
    source_label = "Data.csv (legacy format)" if (not has_ledger) else "ledger.json"
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
        utils.print_single_separator(f"Sample '{sample_ID}' successfully analysed in {total_process_time:.1f} sec.")
    else:
        utils.print_single_separator(f"Analysis was not successful for '{sample_ID}'.")
        comp_analyzer._warn_peak_overlaps()

    _keep_or_revert_new_analysis(
        ledger_path, previous_quantification_id, previous_index, n_previous_analyses, activate_new_analysis
    )
    return comp_analyzer