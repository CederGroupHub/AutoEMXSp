#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch quantification and analysis of X-ray spectra for a list of samples.

This module provides automated batch quantification and (optionally) clustering/statistical
analysis of acquired X-ray spectra for multiple samples. It is robust to missing files or
errors in individual samples, making it suitable for unattended batch processing.

Import this module in your own code and call the
`batch_quantify_and_analyze()` function, passing your desired sample IDs and
options as arguments. This enables integration into larger workflows or pipelines.

Workflow
--------
- Loads sample configurations and spectral data from ``ledger.json``. Results saved in the
  legacy ``Data.csv`` format are converted to a ``ledger.json`` automatically.
- Performs quantification and optionally clustering/statistical analysis and saves results.

Notes
-----
- Only the ``sample_ID`` is required if acquisition output is saved in the default directory;
  otherwise, specify ``results_path``.
- Designed to continue processing even if some samples are missing or have errors.

Typical usage:
    - Edit the `sample_IDs` list and parameter options in the script, or
    - Import and call `batch_quantify_and_analyze()` with your own arguments.
    
Returns
-------
quant_results : list()
    List of EMXSp_Composition_Analyzer, the composition analysis object containing the results and methods for further analysis.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import os
import time
import logging
import traceback
import warnings
from typing import List, Optional

from autoemx.utils.helper import (
    print_double_separator,
    get_sample_dir,
)
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.config.ledger_io import ledger_runtime_configs, load_sample_ledger
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["batch_quantify_and_analyze"]

def batch_quantify_and_analyze(
    sample_IDs: List[str],
    els_sample: Optional[List[str]] = None,
    quantification_method: Optional[str] = None,
    results_path: Optional[str] = None,
    min_bckgrnd_cnts: Optional[float] = None,
    min_total_counts_fraction: Optional[float] = None,
    identify_missing_elements: Optional[bool] = None,
    output_filename_suffix: str = "",
    use_instrument_background: bool = dflt.use_instrument_background,
    max_analytical_error: float = dflt.max_analytical_error_percent,
    run_analysis: bool = True,
    num_CPU_cores: Optional[int] = None,
    force_requantification: bool = False,
    interrupt_fits_bad_spectra: bool = False,
    use_project_specific_std_dict: Optional[bool] = None,
    is_known_phase_mixture: Optional[bool] = None,
    intermixing_violin_plots: Optional[bool] = None,
    standards_dict: Optional[dict] = None,
    spectrum_lims: Optional[tuple] = None,
    els_substrate: Optional[List[str]] = None,
    fit_tolerance: Optional[float] = None,
    requantify_only_unquantified_spectra: bool = False,
    max_spectra_to_quantify: Optional[int] = None,
    is_known_precursor_mixture: Optional[bool] = None,
) -> List[EMXSp_Composition_Analyzer]:
    """
    Batch quantification and analysis for a list of samples.

    Parameters
    ----------
    sample_IDs : List[str]
        List of sample identifiers.
    els_sample : list(str), optional
        List of elements in the sample. If the first entry is "" or None, the rest of the list is appended to the 
        list loaded from the ledger; otherwise, the provided list replaces it.
    els_substrate : list(str), optional
        List of substrate elements. If None, uses the substrate elements from the
        active quantification configuration in each ledger. If the first entry is
        "" or None, the remaining entries are appended; otherwise, the provided
        list replaces the ledger value.
    quantification_method : str, optional
        Method to use for quantification. Uses quant_cfg.method if unspecified. Currently only supports 'PB'.
    spectrum_lims : tuple, optional
        Lower and upper channel-index limits for spectrum fitting. If None, uses
        ``spectrum_lims`` from the active quantification configuration in each ledger.
    results_path : str, optional
        Base directory where results are stored. Default: autoemx/Results
    min_bckgrnd_cnts : float, optional
        Minimum background counts under the reference peaks: spectra below are flagged (quant_flag 8)
        (ClusteringConfig.min_bckgrnd_cnts). Changing it needs no re-quantification.
        If None, keeps the value of the sample's active analysis. Default: None
    min_total_counts_fraction : float, optional
        Minimum accepted total spectrum counts as a fraction of ``target_acquisition_counts``.
        Spectra below this threshold are flagged ``quant_flag = 2`` ("Total counts too low").
        If None, uses the value already on the quantification config (default ``0.9``).
        Lower this to quantify shorter acquisitions; ``0`` disables the check.
    identify_missing_elements : bool, optional
        If True, each spectrum is checked for elements missing from its fit before quantification
        (e.g. contaminants, or elements of an unknown material), and a report of the elements found is printed
        for each sample and for the batch. If None, uses the value of the quantification config (default True;
        off for quantification runs saved before this option existed).
    output_filename_suffix : str, optional
        Suffix to append to output filenames.
    use_instrument_background : bool, optional
        Whether to use instrument background if present (Default: False).
    max_analytical_error : float, optional
        Maximum allowed analytical error for analysis.
    run_analysis : bool, optional
        Whether to run clustering/statistical analysis after quantification.
    num_CPU_cores : bool | None, optional
        Number of CPU cores to use during fitting and quantification. If None, half of the available cores are used.
    force_requantification : bool, optional
        If True, re-quantifies all spectra regardless of whether an existing
        quantification run with matching settings is already available.
    interrupt_fits_bad_spectra : bool, optional
        Controls early-exit behaviour during iterative spectral fitting.

        If ``True`` (default), the fit is aborted as soon as any of the following
        conditions is detected mid-iteration:

        - Reduced chi-squared exceeds 20 % of total spectrum counts (poor fit, flag 4).
        - Analytical error exceeds 50 w% (flag 5).
        - Excessive X-ray absorption around reference peaks (flag 6).
        - Total counts below ``min_total_counts_fraction`` of target (flag 2; default 90%).
        - Low-energy background counts below threshold (flag 3).

        The aborted spectrum is stored with ``QuantificationDiagnostics.interrupted=True``
        and no composition is saved.  This speeds up batch quantification significantly
        when many spectra are expected to be unreliable.

        If ``False``, early-exit is disabled for the current run.  Any spectrum whose
        previous record has ``interrupted=True`` (from a prior run with
        ``interrupt_fits_bad_spectra=True``) is re-quantified, and its ledger record is
        overwritten with the new result.
    use_project_specific_std_dict : bool, optional
        If True, loads standards from project folder (i.e. results_dir) during quantification.
        Default: None. Loads it from quant_cfg file
    is_known_phase_mixture : bool, optional
        Whether the sample is a mixture of known phases (the candidate phases). If True, spectra are quantified with
        the reference P/B values of these phases, and the analysis characterizes their intermixing.
        Default: None, keeps the value of the sample's active quantification run. See example at:
            L. N. Walters et al., Synthetic Accessibility and Sodium Ion Conductivity of the Na 8– x A x P 2 O 9 (NAP)
            High-Temperature Sodium Superionic Conductor Framework, Chem. Mater. 37, 6807 (2025).
    intermixing_violin_plots : bool, optional
        Whether to plot the intermixing of the reported mixtures as violin plots (PlotConfig.intermixing_violin_plots).
        Default: None, keeps the saved setting (by default, plotted for known phase mixtures).
    standards_dict : dict, optional
        Dictionary of reference PB values from experimental standards. Default : None.
        If None, dictionary of standards is loaded from the XSp_calibs/Your_Microscope_ID directory.
        Provide standards_dict only when providing different standards from those normally used for quantification.
    fit_tolerance : float, optional
        Tolerance for fit convergence. If None, uses ``fit_tolerance`` from the active quantification
        configuration in each ledger (default 1e-4).
    requantify_only_unquantified_spectra : bool, optional
        Reuse the latest matching quantification run, but re-quantify only spectra without a
        composition result (never quantified, or previously skipped/flagged). Default: False.
    max_spectra_to_quantify : int, optional
        Quantify at most this number of spectra per sample, among those that need quantification
        (in acquisition order), e.g. for quick tests. None (default) quantifies all of them.
    is_known_precursor_mixture : bool, optional
        Deprecated, use `is_known_phase_mixture`.
            
    Returns
    -------
    quant_results : list()
        List of EMXSp_Composition_Analyzer, the composition analysis object containing the results and methods for further analysis.
    """
    if is_known_precursor_mixture is not None:
        warnings.warn("'is_known_precursor_mixture' is deprecated, use 'is_known_phase_mixture' instead.", FutureWarning, stacklevel=2)
        if is_known_phase_mixture is None:
            is_known_phase_mixture = is_known_precursor_mixture

    if results_path is None:
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)

    non_empty_sample_IDs = [s for s in sample_IDs if s]  # filters out None, empty strings, etc.

    if els_sample is not None and len([el for el in els_sample if el]) > 0 and len(non_empty_sample_IDs) > 1:
        logging.info(f"⚠️ Warning: More than one sample ID provided. Using the same provided sample elements list for all samples: {els_sample}.\nEnsure this behaviour is intended, or provide sample-specific element lists by leaving els_sample as None and specifying them in the ledgers.")
    if els_substrate is not None and len([el for el in els_substrate if el]) > 0 and len(non_empty_sample_IDs) > 1:
        logging.info(f"⚠️ Warning: More than one sample ID provided. Using the same provided substrate elements list for all samples: {els_substrate}.\nEnsure this behaviour is intended, or provide sample-specific substrate element lists by leaving els_substrate as None and specifying them in the ledgers.")
    
    quant_results = []
    for sample_ID in sample_IDs:
        try:
            sample_dir = get_sample_dir(results_path, sample_ID)
        except Exception as e:
            logging.warning("Failed to get sample directory for %s: %s", sample_ID, e)
            continue
        
        ledger_path = os.path.join(sample_dir, f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
        spectral_info_f_path = ledger_path
        
        print("\n \n")
        print_double_separator()
        print_double_separator()
        logging.info(f"Sample '{sample_ID}'")
        
        try:
            ledger = load_sample_ledger(ledger_path)
            ledger_cfgs = ledger_runtime_configs(ledger)
            active_quant_config = ledger_cfgs.active_quant_config
        except FileNotFoundError:
            logging.warning(f"Could not find {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
            continue
        except Exception as e:
            logging.warning(f"Error loading {spectral_info_f_path}. Skipping sample '{sample_ID}': {e}")
            continue

        sample_processing_time_start = time.time()
    
        # Retrieve configuration objects for this sample
        microscope_cfg       = ledger_cfgs.microscope_cfg
        sample_cfg           = ledger_cfgs.sample_cfg
        measurement_cfg      = ledger_cfgs.measurement_cfg
        sample_substrate_cfg = ledger_cfgs.sample_substrate_cfg
        quant_cfg            = ledger_cfgs.quant_cfg
        clustering_cfg       = ledger_cfgs.clustering_cfg
        plot_cfg             = ledger_cfgs.plot_cfg
        powder_meas_cfg      = ledger_cfgs.powder_meas_cfg
        bulk_meas_cfg        = ledger_cfgs.bulk_meas_cfg
        
        # Sample elements
        if els_sample is not None and (els_sample[0] == "" or els_sample[0] is None):
            sample_cfg.elements = sample_cfg.elements + els_sample[1:]
        elif els_sample is not None and len(els_sample) > 0:
            sample_cfg.elements = els_sample

        # Sample substrate elements. By default, retain the elements recorded for
        # the active quantification rather than the ledger's current base config.
        ledger_substrate_elements = (
            list(active_quant_config.substrate_elements)
            if active_quant_config is not None
            else list(sample_substrate_cfg.elements)
        )
        if els_substrate is None:
            sample_substrate_cfg.elements = ledger_substrate_elements
        elif len(els_substrate) > 0 and (els_substrate[0] == "" or els_substrate[0] is None):
            sample_substrate_cfg.elements = ledger_substrate_elements + els_substrate[1:]
        elif len(els_substrate) > 0:
            sample_substrate_cfg.elements = els_substrate

        
        if min_bckgrnd_cnts is not None:
            clustering_cfg.min_bckgrnd_cnts = min_bckgrnd_cnts
        if min_total_counts_fraction is not None:
            quant_cfg.min_total_counts_fraction = min_total_counts_fraction
        if identify_missing_elements is not None:
            quant_cfg.identify_missing_elements = bool(identify_missing_elements)
        if quantification_method is not None:
            quant_cfg.method = quantification_method
        if spectrum_lims is not None:
            quant_cfg.spectrum_lims = spectrum_lims
        if fit_tolerance is not None:
            quant_cfg.fit_tolerance = float(fit_tolerance)
        if use_project_specific_std_dict is not None:
            quant_cfg.use_project_specific_std_dict = use_project_specific_std_dict
        if use_instrument_background is not None:
            quant_cfg.use_instrument_background = use_instrument_background

        if is_known_phase_mixture is not None:
            quant_cfg.is_known_phase_mixture = bool(is_known_phase_mixture)
        if intermixing_violin_plots is not None:
            plot_cfg.intermixing_violin_plots = bool(intermixing_violin_plots)

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
            standards_dict=standards_dict,
            output_filename_suffix=output_filename_suffix,
            verbose=True,
            results_dir=sample_dir
        )

        logging.info(f"Quantifying spectra for '{sample_ID}'.")
        
        try:
            comp_analyzer.run_quantification(
                force_requantification=force_requantification,
                requantify_only_unquantified_spectra=requantify_only_unquantified_spectra,
                interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
                num_CPU_cores=num_CPU_cores,
                max_spectra_to_quantify=max_spectra_to_quantify,
                warn_peak_overlaps_at_end=not run_analysis,
            )
        except FileNotFoundError as e:
            logging.warning(f"Skipping sample '{sample_ID}': {e}")
            continue
        except Exception:
            tb_str = traceback.format_exc()  # get full traceback as a string
            logging.warning(
                f"Error during spectral quantification for '{sample_ID}'. Skipping sample.\nFull traceback:\n{tb_str}"
            )
            continue

        # Perform analysis and print results
        if run_analysis:
            try:
                comp_analyzer.output_filename_suffix = output_filename_suffix
                analysis_successful, _, _ = comp_analyzer.analyse_data(max_analytical_error)
            except Exception as e:
                logging.exception(f"Error during clustering analysis for '{sample_ID}'. Rerun separately if needed: {e}")
                continue
            if analysis_successful:
                comp_analyzer.print_results()
            else:
                logging.info(f"Analysis was not successful for '{sample_ID}'.")
                comp_analyzer._warn_peak_overlaps()

        total_process_time = (time.time() - sample_processing_time_start) / 60
        print_double_separator()
        logging.info(f"Sample '{sample_ID}' successfully quantified in {total_process_time:.1f} min.")
        try:
            quantified_spectra = len(comp_analyzer._load_or_create_ledger().spectra)
        except Exception:
            quantified_spectra = len(comp_analyzer.spectra_quant_records)
        logging.info(f"{quantified_spectra} spectra have been quantified and saved for '{sample_ID}'.")

        
        quant_results.append(comp_analyzer)

    # Elements found in the batch
    checked = [a for a in quant_results if a._identify_missing_elements()]
    if len(checked) > 1:
        print_double_separator()
        logging.info("Element identification — batch summary:\n" + "\n".join(
            EMXSp_Composition_Analyzer.format_element_identification_summary(a.element_identification_summary(), a.sample_cfg.ID)
            for a in checked
        ))

    return quant_results
