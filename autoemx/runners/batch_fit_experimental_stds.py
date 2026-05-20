#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch fitting of experimental standards spectra (PB ratios) from existing ledgers.

This script is analogous to batch_acquire_experimental_stds.py, but instead of running collection and quantification,
it directly fits and updates standards from existing data using run_exp_std_fit.

Usage:
    - Edit the 'std_list' list to define your standards (same as in batch_acquire_experimental_stds.py)
    - Adjust configuration parameters as needed
    - Run the script, or import and call batch_fit_experimental_stds() to fit/update standards for one or multiple samples

Created on May 20, 2026
@author: Andrea
"""

import logging
import os
import time
from typing import List, Dict, Any, Optional, Tuple, cast
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
from autoemx.utils import print_double_separator, get_sample_dir
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    QuantificationOptionsConfig,
    ClusteringConfig,
    PowderMeasurementConfig,
    BulkMeasurementConfig,
    ExpStandardsConfig
)
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflts
import autoemx.calibrations as calibs

import autoemx.config as config_module
config_classes_dict = cast(Dict[str, Any], getattr(config_module, "config_classes_dict"))


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["batch_fit_experimental_stds"]

def batch_fit_experimental_stds(
    std_ids: List[str],
    spectrum_lims: Tuple[float, float] = dflts.spectrum_lims,
    use_instrument_background: bool = dflts.use_instrument_background,
    min_bckgrnd_cnts: float = 5,
    update_std_library = True,
    els_substrate: Optional[List[str]] = None,
    output_filename_suffix: str = '',
    verbose: bool = True,
    stds_path: Optional[str] = None,
) -> List[Optional[EMXSp_Composition_Analyzer]]:
    """
    Batch acquisition (and optional quantification) of X-ray spectra for a list of powder samples.

    Parameters
    ----------
    std_ids : list of str
        List of identifiers for the experimental standards.
    spectrum_lims : tuple of float, optional
        Lower and upper energy limits for spectrum fitting in eV.
        Default is (14, 1100) (QuantificationOptionsConfig.spectrum_lims).
    use_instrument_background : bool, optional
        Whether to use instrument background files during fitting.
        If False, background is computed during fitting.
        Default is False (QuantificationOptionsConfig.use_instrument_background).
    min_bckgrnd_cnts : float, optional
        Minimum background counts required for a spectrum not to be filtered out.
        Default is 5 (ClusteringConfig.min_bckgrnd_cnts).
    update_std_library : bool, optional
        If True, the experimental standard library will be updated with the newly fitted PB ratios.
        Default is True.
    els_substrate : Optional[List[str]], optional
        List of elements in the substrate. If None, defaults to all elements in the sample.
        Default is None.
    output_filename_suffix : str, optional
        Suffix for the output filename. Default is ''.
    verbose : bool, optional
        If True, print verbose output. Default is True.
    stds_path : Optional[str], optional
        Directory containing experimental standard data. If None, uses default directory.
        Default is None.

    Returns
    -------
    results : list(EMXSp_Composition_Analyzer)
        A list of the composition analysis objects (one per sample) containing the results and methods for further analysis.
    """
    if stds_path is None:
        stds_path = os.path.join(os.getcwd(), cnst.STDS_DIR)

    if els_substrate is None:
        els_substrate = dflts.substrate_els

    results: List[Optional[EMXSp_Composition_Analyzer]] = []
    comp_analyzer: Optional[EMXSp_Composition_Analyzer] = None
    
    for std_id in std_ids:
        try:
            std_path = get_sample_dir(stds_path, std_id)
        except Exception as e:
            logging.warning("Failed to get sample directory for %s: %s", std_id, e)
            continue

        ledger_path = os.path.join(std_path, f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
        
        print("\n \n")
        print_double_separator()
        print_double_separator()
        logging.info(f"Standard '{std_id}'")

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
            logging.warning(f"Could not find {ledger_path}. Skipping sample '{std_id}'.")
            continue
        except Exception as e:
            logging.warning(f"Error loading {ledger_path}. Skipping sample '{std_id}': {e}")
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
            logging.warning(f"Missing configuration '{e.args[0]}' in {ledger_path}. Skipping sample '{std_id}'.")
            continue

        if quant_cfg is None:
            quant_cfg = config_classes_dict[cnst.QUANTIFICATION_CFG_KEY]()
        
        if min_bckgrnd_cnts is not None:
            clustering_cfg.min_bckgrnd_cnts = min_bckgrnd_cnts
        if use_project_specific_std_dict is not None:
            quant_cfg.use_project_specific_std_dict = use_project_specific_std_dict
        if use_instrument_background is not None:
            quant_cfg.use_instrument_background = use_instrument_background

        # --- Run Composition Analysis or Spectral Acquisition
        comp_analyzer = EMXSp_Composition_Analyzer(
            microscope_cfg=microscope_cfg,
            sample_id=std_id,
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

        # --- Run Composition Analyzer
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
            exp_stds_cfg=exp_stds_cfg,
            is_acquisition=True,
            development_mode=development_mode,
            output_filename_suffix=output_filename_suffix,
            verbose=verbose,
            results_dir=exp_std_dir
        )
        
        try:
            # Always run fitting (run_fitting=True) since this runner is for batch fitting from existing ledgers
            comp_analyzer.run_exp_std_fit(run_fitting=True, update_std_library=update_std_library)
            results.append(comp_analyzer)
        except Exception as e:
            results.append(None)
            logging.exception(f"Sample '{sample_ID}': acquisition/fitting failed: {e}")
            continue
        
    
    # Put microscope in standby after completion
    if comp_analyzer is not None and not development_mode and len(stds) > 1:
        try:
            comp_analyzer.EM_controller.standby()
        except Exception as e:
            logging.warning(f"Could not put microscope in standby: {e}")

    # Generate a standards-coverage report only when standards library was updated.
    if update_std_library:
        try:
            from autoemx.runners.extract_experimental_standards_details import extract_experimental_standards_details  # type: ignore

            if exp_std_dir is not None:
                standards_dir = Path(exp_std_dir).expanduser().resolve()
            else:
                standards_dir = Path(__file__).resolve().parents[1] / cnst.CALIBS_DIR / microscope_ID

            standards_json_filename = f"{measurement_type}_{cnst.STD_FILENAME}_{int(beam_energy):d}keV.json"
            standards_json_path = standards_dir / standards_json_filename

            report_kwargs = {
                "microscope_ID": microscope_ID,
                "voltage": beam_energy,
                "report_output_dir": str(standards_dir),
                "print_report": True,
                "print_per_peak": False,
            }
            if standards_json_path.exists():
                report_kwargs["standards_json_path"] = str(standards_json_path)

            report_path = extract_experimental_standards_details(**report_kwargs)
            logging.info("Standards details report generated at: %s", report_path)
        except Exception as e:
            logging.warning("Could not generate standards details report after collection: %s", e)
    
    return results
    