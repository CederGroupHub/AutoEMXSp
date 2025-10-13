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

Workflow:
    - Loads sample configurations from `Spectra_collection_info.json`
    - Loads acquired spectral data from `Data.csv`
    - Performs quantification (optionally only on unquantified spectra)
    - Optionally performs clustering/statistical analysis and saves results

Notes
-----
- Only the `sample_ID` is required if acquisition output is saved in the default directory;
  otherwise, specify `results_path`.
- Designed to continue processing even if some samples are missing or have errors.

Typical usage:
    - Edit the `sample_IDs` list and parameter options in the script, or
    - Import and call `batch_quantify_and_analyze()` with your own arguments.
    
Parameters
----------
sample_IDs : List[str]
    List of sample identifiers.
quantification_method : str, optional
    Method to use for quantification. Uses quant_cfg.method if unspecified. Currently only supports 'PB'.
results_path : str, optional
    Base directory where results are stored. Default: autoemxsp/Results
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
quantify_only_unquantified_spectra : bool, optional
    If True, only quantify spectra that lack analytical error.
interrupt_fits_bad_spectra : bool, optional
    If True, interrupt fitting if bad spectra are encountered. Speeds up computations
is_known_precursor_mixture : bool, optional
    Whether sample is a mixture of two known powders. Used to characterize extent of intermixing in powders.
    See example at:
        L. N. Walters et al., Synthetic Accessibility and Sodium Ion Conductivity of the Na 8– x A x P 2 O 9 (NAP)
        High-Temperature Sodium Superionic Conductor Framework, Chem. Mater. 37, 6807 (2025).
    
Returns
-------
quant_results : list()
    List of EMXSp_Composition_Analyzer, the composition analysis object containing the results and methods for further analysis.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import os
import pandas as pd
import warnings
import time
import logging
from typing import List, Optional

from autoemxsp.tools.utils import (
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

def batch_quantify_and_analyze(
    sample_IDs: List[str],
    quantification_method: str = None,
    results_path: str = None,
    min_bckgrnd_cnts: float = None,
    output_filename_suffix: str = "",
    use_instrument_background: bool = False,
    max_analytical_error: float = 5,
    run_analysis: bool = True,
    quantify_only_unquantified_spectra: bool = False,
    interrupt_fits_bad_spectra: bool = False,
    is_known_precursor_mixture: Optional[bool] = None,
    standards_dict: dict = None,
) -> None:
    """
    Batch quantification and analysis for a list of samples.

    Parameters
    ----------
    sample_IDs : List[str]
        List of sample identifiers.
    quantification_method : str, optional
        Method to use for quantification. Uses quant_cfg.method if unspecified. Currently only supports 'PB'.
    results_path : str, optional
        Base directory where results are stored. Default: autoemxsp/Results
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
    quantify_only_unquantified_spectra : bool, optional
        If True, only quantify spectra that lack analytical error.
    interrupt_fits_bad_spectra : bool, optional
        If True, interrupt fitting if bad spectra are encountered. Speeds up computations
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
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        results_path = os.path.join(parent_dir, cnst.RESULTS_DIR)
        
    quant_results = []
    for sample_ID in sample_IDs:
        sample_dir = get_sample_dir(results_path, sample_ID)
        spectral_info_f_path = os.path.join(sample_dir, f"{cnst.ACQUISITION_INFO_FILENAME}.json")
        data_path = os.path.join(sample_dir, f"{cnst.DATA_FILENAME}{cnst.DATA_FILEEXT}")
        
        print_double_separator()
        logging.info(f"Sample '{sample_ID}'")
        
        try:
            configs, metadata = load_configurations_from_json(spectral_info_f_path, config_classes_dict)
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
            clustering_cfg      = configs[cnst.CLUSTERING_CFG_KEY]
            plot_cfg            = configs[cnst.PLOT_CFG_KEY]
            powder_meas_cfg     = configs.get(cnst.POWDER_MEASUREMENT_CFG_KEY, None)  # Optional
            bulk_meas_cfg     = configs.get(cnst.BULK_MEASUREMENT_CFG_KEY, None)  # Optional
        except KeyError as e:
            logging.warning(f"Missing configuration '{e.args[0]}' in {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
            continue
        
        if min_bckgrnd_cnts is not None:
            quant_cfg.min_bckgrnd_cnts = min_bckgrnd_cnts
        if quantification_method is not None:
            quant_cfg.method = quantification_method
            
        # Load 'Data.csv' into a DataFrame
        try:
            spectra_quant, spectral_data, sp_coords, original_df = extract_spectral_data(data_path)
        except Exception as e:
            logging.warning(f"Could not load spectral data for '{sample_ID}': {e}")
            continue

        if use_instrument_background:
            if getattr(spectral_data, 'get', None) and spectral_data.get(cnst.BACKGROUND_DF_KEY, []) == []:
                warnings.warn(
                    "Background column not found in input data. "
                    "Spectral background will be computed instead."
                )
                
        # Change is_known_precursor_mixture if provided
        if is_known_precursor_mixture:
            powder_meas_cfg.is_known_precursor_mixture = is_known_precursor_mixture

        # Decide which spectra to quantify
        quant_cfg.interrupt_fits_bad_spectra = interrupt_fits_bad_spectra
        if quantify_only_unquantified_spectra and cnst.AN_ER_DF_KEY in original_df.columns:
            to_quantify_mask = original_df[cnst.AN_ER_DF_KEY].isna()
            if not to_quantify_mask.any():
                logging.info("All spectra already quantified. Skipping.")
                continue
            indices_to_quantify = list(original_df.index[to_quantify_mask])
            logging.info(f"Quantifying {len(indices_to_quantify)} unquantified spectra out of {len(original_df)}.")
        else:
            indices_to_quantify = list(original_df.index)
            logging.info(f"Quantifying all {len(original_df)} spectra.")
            
        # Subset the data for quantification
        spectral_data_sub = {key: [spectral_data[key][i] for i in indices_to_quantify] for key in cnst.LIST_SPECTRAL_DATA_QUANT_KEYS}
        sp_coords_sub = [sp_coords[i] for i in indices_to_quantify]

        # Handle output filename for analyzer
        if quantify_only_unquantified_spectra:
            temp_suffix = output_filename_suffix + '_temp'
            output_suffix = temp_suffix
        else:
            output_suffix = output_filename_suffix

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
            standards_dict=standards_dict,
            output_filename_suffix=output_suffix,
            verbose=True,
            results_dir=sample_dir
        )

        comp_analyzer.sp_coords = sp_coords_sub
        for key in cnst.LIST_SPECTRAL_DATA_QUANT_KEYS:
            comp_analyzer.spectral_data[key] = spectral_data_sub[key]
        
        try:
            comp_analyzer.run_quantification()
        except Exception as e:
            logging.warning(f"Error during spectral quantification for '{sample_ID}'. Skipping sample: {e}")
            continue

        if quantify_only_unquantified_spectra:
            # Read the temporary file produced by the analyzer
            temp_data_path = os.path.join(sample_dir, f"{cnst.DATA_FILENAME}{output_suffix}{cnst.DATA_FILEEXT}")
            try:
                temp_df = pd.read_csv(temp_data_path)
            except Exception as e:
                logging.error(f"Could not read temporary quantification file '{temp_data_path}' for '{sample_ID}': {e}")
                continue
            temp_df.index = indices_to_quantify

            # Update only the quantified rows in the original DataFrame
            # Use vectorized assignment for performance
            try:
                original_df.loc[indices_to_quantify, temp_df.columns] = temp_df.values
            except Exception as e:
                logging.error(f"Failed to update quantified rows in DataFrame for '{sample_ID}': {e}")
                continue

            # Convert quantification flag column to int if present
            if cnst.QUANT_FLAG_DF_KEY in original_df.columns:
                original_df[cnst.QUANT_FLAG_DF_KEY] = original_df[cnst.QUANT_FLAG_DF_KEY].astype(int)

            # Save the updated Data.csv with counter logic to avoid overwriting
            new_data_path = os.path.join(sample_dir, f'{cnst.DATA_FILENAME}{output_filename_suffix}{cnst.DATA_FILEEXT}')
            cntr = 1
            while os.path.exists(new_data_path):
                cntr += 1
                new_data_path = os.path.join(sample_dir, f'{cnst.DATA_FILENAME}{output_filename_suffix}_{cntr}{cnst.DATA_FILEEXT}')

            try:
                original_df.to_csv(new_data_path, index=False)
            except Exception as e:
                logging.error(f"Failed to save updated {cnst.DATA_FILENAME}{cnst.DATA_FILEEXT} for '{sample_ID}': {e}")
                continue

            # Remove the temporary file
            try:
                os.remove(temp_data_path)
            except Exception as e:
                logging.warning(f"Could not remove temporary file '{temp_data_path}': {e}")

            # Reload data for analysis if needed
            if run_analysis:
                try:
                    spectra_quant, spectral_data, sp_coords, _ = extract_spectral_data(new_data_path)
                    comp_analyzer.spectra_quant = spectra_quant
                    comp_analyzer.sp_coords = sp_coords
                    comp_analyzer.spectral_data = spectral_data
                except Exception as e:
                    logging.warning(f"Could not reload data for analysis for '{sample_ID}': {e}")
                    continue

        # Perform analysis and print results
        if run_analysis:
            try:
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
        logging.info(f"{len(indices_to_quantify)} spectra have been quantified and saved for '{sample_ID}'.")
        
        quant_results.append(comp_analyzer)
    
    return quant_results
