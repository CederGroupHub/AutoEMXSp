#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch fitting of multiple X-ray spectra for fitting parameter extraction.

This module allows running the fitting step for multiple spectra across multiple samples.
It uses the `fit_and_quantify_spectrum` function internally with quantification disabled,
then extracts the values of desired fitting parameters.

Features:
---------
- Accepts a list of sample IDs and spectrum IDs.
- Supports 'all' spectra mode for each sample.

Example:
--------
>>> from autoemxsp.runners import batch_fit_spectra
>>> batch_fit_spectra(
...     sample_IDs=["Sample1", "Sample2"],
...     spectrum_IDs="all",
...     plot_signal=False
... )

Author:
-------
Andrea (adapted for batch processing)
"""

import os
import logging
import pandas as pd
from datetime import datetime
from typing import List, Optional

import autoemxsp.tools.constants as cnst
from autoemxsp.tools.utils import get_sample_dir
from autoemxsp.runners.fit_and_quantify_spectrum import fit_and_quantify_spectrum

# Configure logging (same style as fit_and_quantify_spectrum)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def batch_fit_spectra(sample_IDs,
                      spectrum_IDs,
                      is_standard: bool,
                      fit_params_vals_to_extract: Optional[List[str]] = None,
                      spectrum_lims: tuple = None,
                      samples_path: str = None,
                      output_path: str = 'Fitting output',
                      use_instrument_background: bool = False,
                      quantify_plot: bool = True,
                      plot_signal: bool = True,
                      zoom_plot: bool = False,
                      line_to_plot: str = '',
                      els_substrate: list = None,
                      fit_tol: float = 1e-4,
                      is_particle: bool = True,
                      max_undetectable_w_fr: float = 0,
                      force_single_iteration: bool = False,
                      interrupt_fits_bad_spectra: bool = False,
                      print_results: bool = True,
                      quant_verbose: bool = True,
                      fitting_verbose: bool = True
):
    """
    Run fitting for multiple spectra across multiple samples to extract values of fitting parameters.
        
    Parameters
    ----------
    sample_IDs : list of str
        List of sample identifiers.
    spectrum_IDs : list of int or str
        List of spectrum IDs to process (values reported in 'Spectrum #' column in Data.csv),
        or 'all' to process all spectra in each sample.
    is_standard : bool
        Defines whether measurement is of a standard (i.e., well defined composition) or not
    fit_params_vals_to_extract : list of str, optional
        List of fitting parameter names whose value to extract and save
    samples_path : str, optional
        Base directory where results are stored. Default: autoemxsp/Results
    output_path : str, optional
        Directory where the extracted values of fitted parameters are saved. Default: /Fitting output 
    use_instrument_background : bool, optional
        Whether to use instrument background if present.
    quantify_plot : bool, optional
        Whether to plot quantification results.
    plot_signal : bool, optional
        Whether to plot the signal.
    zoom_plot : bool, optional
        Whether to zoom on a specific line.
    line_to_plot : str, optional
        Line to zoom in plot.
    els_substrate : list, optional
        List of substrate elements.
    fit_tol : float, optional
        Fit tolerance.
    is_particle : bool, optional
        If True, treat sample as particle (powder).
    max_undetectable_w_fr : float, optional
        Maximum allowed weight fraction for undetectable elements (default: 0).
    force_single_iteration : bool, optional
        If True, quantification will be run for a single iteration only (default: False).
    interrupt_fits_bad_spectra : bool, optional
        If True, interrupt fitting if bad spectra are detected (default: False).
    print_results : bool, optional
        If True, prints all fitted parameters and their values (default: True).
    quant_verbose : bool, optional
        If True, prints quantification operations
    fitting_verbose : bool, optional
        If True, prints fitting operations
        
    Returns
    -------
    quantifier : XSp_Quantifier
        The quantifier object containing the results, fit parameters, and methods for further analysis and plotting.
    """
    
    if samples_path is None:
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        samples_path = os.path.join(parent_dir, cnst.RESULTS_DIR)

    logging.info("Starting batch fitting process...")
    logging.info(f"Results path: {samples_path}")
    
    extracted_par_vals = {}
    for sample_ID in sample_IDs:
        logging.info(f"Processing sample '{sample_ID}'...")
        sample_dir = get_sample_dir(samples_path, sample_ID)
        data_filename = cnst.STDS_MEAS_FILENAME if is_standard else cnst.DATA_FILENAME
        data_path = os.path.join(sample_dir, f"{data_filename}.csv")

        if not os.path.exists(data_path):
            logging.warning(f"Data file not found for sample '{sample_ID}'. Skipping.")
            continue

        try:
            df = pd.read_csv(data_path)
        except Exception as e:
            logging.warning(f"Could not read {data_path} for sample '{sample_ID}': {e}")
            continue

        if cnst.SP_ID_DF_KEY not in df.columns:
            logging.warning(f"Column '{cnst.SP_ID_DF_KEY}' not found in {data_path}. Skipping sample '{sample_ID}'.")
            continue

        # Determine spectra to process
        if spectrum_IDs == 'all' or (isinstance(spectrum_IDs, list) and len(spectrum_IDs) == 1 and spectrum_IDs[0] == 'all'):
            spectra_to_process = df[cnst.SP_ID_DF_KEY].unique()
            logging.info(f"Found {len(spectra_to_process)} spectra for sample '{sample_ID}'.")
        else:
            spectra_to_process = spectrum_IDs
            logging.info(f"Processing specified spectra for sample '{sample_ID}': {spectra_to_process}")
        
        sample_fit_results = []
        for sp_id in spectra_to_process:
            logging.info(f"Fitting Sample '{sample_ID}', Spectrum {sp_id} (fit only, no quantification)")
            try:
                quantifier = fit_and_quantify_spectrum(
                    sample_ID=sample_ID,
                    spectrum_ID=sp_id,
                    is_standard = is_standard,
                    results_path=samples_path,
                    spectrum_lims = spectrum_lims,
                    use_instrument_background = use_instrument_background,
                    quantify_plot = False,
                    plot_signal = plot_signal,
                    zoom_plot = zoom_plot,
                    line_to_plot = line_to_plot,
                    els_substrate = els_substrate,
                    fit_tol = fit_tol,
                    is_particle = is_particle,
                    max_undetectable_w_fr = max_undetectable_w_fr,
                    force_single_iteration = force_single_iteration,
                    interrupt_fits_bad_spectra = interrupt_fits_bad_spectra,
                    print_results=print_results,
                    quant_verbose = quant_verbose,
                    fitting_verbose = fitting_verbose
                )
            except Exception as e:
                logging.exception(f"Error fitting spectrum {sp_id} for sample '{sample_ID}': {e}")
                sample_fit_results.append(None)
            else:
                if fit_params_vals_to_extract and quantifier.bad_quant_flag is None:
                    params = quantifier.fit_result.params
                    extracted_vals = {}
                    for param_name in fit_params_vals_to_extract:
                        if param_name in params:
                            extracted_vals[param_name] = params[param_name].value
                        else:
                            extracted_vals[param_name] = pd.NA
                    sample_fit_results.append({'sp_id' : sp_id, **extracted_vals})
                else:
                    sample_fit_results.append(None)
                    
        if fit_params_vals_to_extract:
            # Create DataFrame from current sample_fit_results
            filtered_results = [item for item in sample_fit_results if item is not None] # Remove None entries
            temp_df = pd.DataFrame(filtered_results)
            
            # Calculate mean and std for numeric columns (excluding sp_id)
            mean_vals = temp_df.drop(columns=['sp_id']).mean(numeric_only=True).to_dict()
            std_vals = temp_df.drop(columns=['sp_id']).std(numeric_only=True).to_dict()
            
            # Append mean and std rows directly to sample_fit_results
            sample_fit_results.append({'sp_id': 'mean', **mean_vals})
            sample_fit_results.append({'sp_id': 'std', **std_vals})
            
            # Now create final DataFrame
            results_df = pd.DataFrame(sample_fit_results)
            
            # Save without index
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            now = datetime.now()
            now_formatted = now.strftime("%Y%m%d_%Hh%Mm")
            file_path = os.path.join(output_path, f"{now_formatted}_{sample_ID}_FitParamVals.csv")
            results_df.to_csv(file_path, index=False)
            
            extracted_par_vals[sample_ID] = sample_fit_results 
        else:
            extracted_par_vals[sample_ID] = None
                    
    logging.info("Batch fitting process completed.")
    
    return extracted_par_vals