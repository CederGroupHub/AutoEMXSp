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
- Loads sample configurations and spectral data from ``ledger.json`` (primary source).
- Falls back to ``Data.csv`` only when no ledger exists; the first quantification run
  then migrates the CSV data into a ``ledger.json`` for all subsequent runs.
- Performs quantification and optionally clustering/statistical analysis and saves results.

Notes
-----
- Only the ``sample_ID`` is required if acquisition output is saved in the default directory;
  otherwise, specify ``results_path``.
- When a ``ledger.json`` is present it is always used as the authoritative source of
  spectral data and quantification results; ``Data.csv`` is only read as a one-time
  migration path when no ledger has been created yet.
- ``spectra_quant`` (the in-memory composition buffer consumed by ``analyse_data``) is
  always populated by ``run_quantification`` from ledger records — the runner never
  needs to pre-build it.
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
from typing import Any, Dict, List, Optional, cast

from autoemx.utils.helper import (
    print_double_separator,
    get_sample_dir,
)
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.config.ledger_io import load_sample_ledger
import autoemx.config as config_module
from autoemx.config.ledger_schemas import ClusteringConfig # type: ignore
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer

config_classes_dict = cast(Dict[str, Any], getattr(config_module, "config_classes_dict"))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["batch_quantify_and_analyze"]

def batch_quantify_and_analyze(
    sample_IDs: List[str],
    quantification_method: Optional[str] = None,
    results_path: Optional[str] = None,
    min_bckgrnd_cnts: Optional[float] = None,
    output_filename_suffix: str = "",
    use_instrument_background: bool = dflt.use_instrument_background,
    max_analytical_error: float = 5,
    run_analysis: bool = True,
    num_CPU_cores: Optional[int] = None,
    force_requantification: bool = False,
    interrupt_fits_bad_spectra: bool = False,
    use_project_specific_std_dict: Optional[bool] = None,
    is_known_precursor_mixture: Optional[bool] = None,
    standards_dict: Optional[dict] = None,
) -> List[EMXSp_Composition_Analyzer]:
    """
    Batch quantification and analysis for a list of samples.

    Parameters
    ----------
    sample_IDs : List[str]
        List of sample identifiers.
    quantification_method : str, optional
        Method to use for quantification. Uses quant_cfg.method if unspecified. Currently only supports 'PB'.
    results_path : str, optional
        Base directory where results are stored. Default: autoemx/Results
    min_bckgrnd_cnts : float, optional
        Minimum number of background counts underneath reference peaks below which spectra are flagged.
        If None, leaves it unchanged. Default: None
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
    is_known_precursor_mixture : bool, optional
        Whether sample is a mixture of two known powders. Used to characterize extent of intermixing in powders.
        See example at:
            L. N. Walters et al., Synthetic Accessibility and Sodium Ion Conductivity of the Na 8– x A x P 2 O 9 (NAP)
            High-Temperature Sodium Superionic Conductor Framework, Chem. Mater. 37, 6807 (2025).
    standards_dict : dict, optional
        Dictionary of reference PB values from experimental standards. Default : None.
        If None, dictionary of standards is loaded from the XSp_calibs/Your_Microscope_ID directory.
        Provide standards_dict only when providing different standards from those normally used for quantification.
            
    Returns
    -------
    quant_results : list()
        List of EMXSp_Composition_Analyzer, the composition analysis object containing the results and methods for further analysis.
    """
    if results_path is None:
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)
        
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
        except FileNotFoundError:
            logging.warning(f"Could not find {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
            continue
        except Exception as e:
            logging.warning(f"Error loading {spectral_info_f_path}. Skipping sample '{sample_ID}': {e}")
            continue

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
            logging.warning(f"Missing configuration '{e.args[0]}' in {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
            continue

        if clustering_cfg is None:
            clustering_cfg = ClusteringConfig()
        if quant_cfg is None:
            quant_cfg = config_classes_dict[cnst.QUANTIFICATION_CFG_KEY]()
        
        if min_bckgrnd_cnts is not None:
            clustering_cfg.min_bckgrnd_cnts = min_bckgrnd_cnts
        if quantification_method is not None:
            quant_cfg.method = quantification_method
        if use_project_specific_std_dict is not None:
            quant_cfg.use_project_specific_std_dict = use_project_specific_std_dict
        if use_instrument_background is not None:
            quant_cfg.use_instrument_background = use_instrument_background

        # Change is_known_precursor_mixture if provided
        if is_known_precursor_mixture is not None:
            if powder_meas_cfg is not None:
                powder_meas_cfg.is_known_powder_mixture_meas = is_known_precursor_mixture
            else:
                from autoemx.config.runtime_configs import PowderMeasurementConfig
                powder_meas_cfg = PowderMeasurementConfig()
                powder_meas_cfg.is_known_powder_mixture_meas = is_known_precursor_mixture

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
                interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
                num_CPU_cores=num_CPU_cores,
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

        total_process_time = (time.time() - sample_processing_time_start) / 60
        print_double_separator()
        logging.info(f"Sample '{sample_ID}' successfully quantified in {total_process_time:.1f} min.")
        try:
            quantified_spectra = len(comp_analyzer._load_or_create_ledger().spectra)
        except Exception:
            quantified_spectra = len(comp_analyzer.spectra_quant_records)
        logging.info(f"{quantified_spectra} spectra have been quantified and saved for '{sample_ID}'.")

        
        quant_results.append(comp_analyzer)
    
    return quant_results
