#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting and quantification of a single X-ray spectrum.

For spectrum-level analysis of fitting and quantification performance.

Import this module in your own code and call the
`fit_and_quantify_spectrum()` function, passing your desired 'sample_ID', 'spectrum_ID' and
options as arguments. This enables integration into larger workflows or pipelines.

Workflow:
    - Loads sample configurations from `Spectra_collection_info.json`
    - Loads acquired spectral data from `Data.csv`
    - Performs quantification (optionally only on unquantified spectra)
    - Optionally performs clustering/statistical analysis and saves results

Notes
-----
- Only the `sample_ID` and 'spectrum_ID' are required if acquisition output is saved in the default Results directory;
  otherwise, specify `results_path`.

Parameters
----------
sample_ID : str
    Sample identifier.
spectrum_ID : int
    Value reported in 'Spectrum #' column in Data.csv.
is_standard : bool
    Defines whether measurement is of a standard (i.e., well defined composition) or not
results_path : str, optional
    Base directory where results are stored. Default: autoemxsp/Results
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
    
Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import os
import warnings
import time
import logging

from autoemxsp.tools.utils import (
    print_double_separator,
    get_sample_dir,
    load_configurations_from_json,
    extract_spectral_data
)
import autoemxsp.tools.constants as cnst
from autoemxsp.tools.config_classes import config_classes_dict, ExpStandardsConfig
from autoemxsp.core.XSp_quantifier import XSp_Quantifier
from autoemxsp.core.EMXSp_composition_analyser import EMXSp_Composition_Analyzer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def fit_and_quantify_spectrum(
    sample_ID: str,
    spectrum_ID: int,
    is_standard: bool = False,
    spectrum_lims: tuple = None,
    results_path: str = None,
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
    standards_dict: dict = None,
    print_results: bool = True,
    quant_verbose: bool = True,
    fitting_verbose: bool = True
):
    """
    Fit and (optionally) quantify a single spectrum.

    Parameters
    ----------
    sample_ID : str
        Sample identifier.
    spectrum_ID : int
        Value reported in 'Spectrum #' column in Data.csv.
    is_standard : bool
        Defines whether measurement is of a standard (i.e., well defined composition) or not
    results_path : str, optional
        Base directory where results are stored. Default: autoemxsp/Results
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
    standards_dict : dict, optional
        Dictionary of reference PB values from experimental standards. Default : None.
        If None, dictionary of standards is loaded from the XSp_calibs/Your_Microscope_ID directory.
        Provide standards_dict only when providing different standards from those normally used for quantification.
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
    if results_path is None:
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        results_path = os.path.join(parent_dir, cnst.RESULTS_DIR)
        
    sample_dir = get_sample_dir(results_path, sample_ID)
    spectral_info_f_path = os.path.join(sample_dir, f"{cnst.ACQUISITION_INFO_FILENAME}.json")
    data_filename = cnst.STDS_MEAS_FILENAME if is_standard else cnst.DATA_FILENAME
    data_path = os.path.join(sample_dir, f"{data_filename}.csv")
    
    print_double_separator()
    logging.info(f"Sample '{sample_ID}', spectrum {spectrum_ID}")
    
    try:
        configs, metadata = load_configurations_from_json(spectral_info_f_path, config_classes_dict)
    except FileNotFoundError:
        logging.warning(f"Could not find {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return
    except Exception as e:
        logging.warning(f"Error loading {spectral_info_f_path}. Skipping sample '{sample_ID}': {e}")
        return

    sample_processing_time_start = time.time()
    
    # Retrieve configuration objects
    try:
        microscope_cfg      = configs[cnst.MICROSCOPE_CFG_KEY]
        sample_cfg          = configs[cnst.SAMPLE_CFG_KEY]
        measurement_cfg     = configs[cnst.MEASUREMENT_CFG_KEY]
        sample_substrate_cfg= configs[cnst.SAMPLESUBSTRATE_CFG_KEY]
        quant_cfg           = configs[cnst.QUANTIFICATION_CFG_KEY]
        clustering_cfg      = configs[cnst.CLUSTERING_CFG_KEY]
        powder_meas_cfg     = configs.get(cnst.POWDER_MEASUREMENT_CFG_KEY, None)  # Optional
        exp_stds_cfg     = configs.get(cnst.EXP_STD_MEASUREMENT_CFG_KEY, None)  # Optional
    except KeyError as e:
        logging.warning(f"Missing configuration '{e.args[0]}' in {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return
    
    # Load experimental standard dictionary if the sample is a known powder precursor mix
    if powder_meas_cfg and powder_meas_cfg.is_known_powder_mixture_meas:
        if not exp_stds_cfg:
            exp_stds_cfg = ExpStandardsConfig()
        comp_analyzer = EMXSp_Composition_Analyzer(microscope_cfg,
                                                   sample_cfg,
                                                   measurement_cfg,
                                                   sample_substrate_cfg,
                                                   quant_cfg,
                                                   clustering_cfg,
                                                   powder_meas_cfg = powder_meas_cfg,
                                                   exp_stds_cfg = exp_stds_cfg,
                                                   standards_dict = standards_dict)
        stds_dict = comp_analyzer.XSp_std_dict
    else:
        stds_dict = None
    
    # Load 'Data.csv' into a DataFrame
    try:
        _, spectral_data, _, original_df = extract_spectral_data(data_path)
    except Exception as e:
        logging.warning(f"Could not load spectral data for '{sample_ID}': {e}")
        return

    # Extract the row corresponding to spectrum_ID
    if cnst.SP_ID_DF_KEY not in original_df.columns:
        logging.error(f"Column '{cnst.SP_ID_DF_KEY}' not in Data.csv for sample '{sample_ID}'.")
        return

    df_match = original_df[original_df[cnst.SP_ID_DF_KEY] == spectrum_ID]
    if df_match.empty:
        logging.warning(f"Spectrum ID {spectrum_ID} not found in Data.csv for sample '{sample_ID}'.")
        return
    
    sp_idx = df_match.index[0]

    # Extract spectrum data from spectral_data (not from df_row)
    try:
        spectrum = spectral_data[cnst.SPECTRUM_DF_KEY][sp_idx]
    except Exception as e:
        logging.warning(f"Spectrum data not found for spectrum ID {spectrum_ID} in sample '{sample_ID}': {e}")
        return

    # Background extraction
    background = None
    bkg_list = spectral_data.get(cnst.BACKGROUND_DF_KEY)
    
    if use_instrument_background:
        if (
            bkg_list is not None
            and isinstance(bkg_list, (list, tuple))
            and len(bkg_list) > sp_idx
            and bkg_list[sp_idx] is not None
        ):
            background = bkg_list[sp_idx]
        else:
            warnings.warn(
                "Instrument background not found or empty for this spectrum. "
                "Spectral background will be computed instead."
            )

    # Collection time
    if 'Live_time' in spectral_data and len(spectral_data['Live_time']) > sp_idx:
        sp_collection_time = spectral_data['Live_time'][sp_idx]
    else:
        sp_collection_time = None

    # Calibration and configuration parameters
    try:
        beam_energy = measurement_cfg.beam_energy_keV
        emergence_angle = measurement_cfg.emergence_angle
        el_to_quantify = sample_cfg.elements
        offset = microscope_cfg.energy_zero
        scale = microscope_cfg.bin_width
    except Exception as e:
        logging.error(f"Error extracting calibration/configuration parameters: {e}")
        return
    
    # Substrate elements
    if els_substrate is None:
        els_substrate = sample_substrate_cfg.elements
    # Spectral limits
    if spectrum_lims is None:
        spectrum_lims = quant_cfg.spectrum_lims
        
    sp_start, sp_end = spectrum_lims
    spectrum_vals = spectrum

    # Quantification (replace with your quantifier class as needed)
    quantifier = XSp_Quantifier(
        spectrum_vals = spectrum_vals,
        spectrum_lims = spectrum_lims,
        microscope_ID = microscope_cfg.ID,
        meas_type = measurement_cfg.type,
        meas_mode = measurement_cfg.mode,
        det_ch_offset=offset,
        det_ch_width=scale,
        beam_e=beam_energy,
        emergence_angle=emergence_angle,
        background_vals=background,
        els_sample=el_to_quantify,
        els_substrate=els_substrate,
        els_w_fr=None,
        is_particle=is_particle,
        sp_collection_time=sp_collection_time,
        max_undetectable_w_fr=max_undetectable_w_fr,
        fit_tol=fit_tol,
        verbose=quant_verbose,
        fitting_verbose=fitting_verbose,
        standards_dict = stds_dict
    )
    
    try:
        if quantify_plot:
            quant_result, _, flag = quantifier.quantify_spectrum(
                force_single_iteration=force_single_iteration,
                interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
                print_result=print_results
                )
        else:
            quantifier.initialize_and_fit_spectrum(
                print_results=print_results
            )
    except Exception as e:
        logging.exception(f"Error during spectral quantification for '{sample_ID}': {e}")
        return

    if plot_signal:
        line_to_zoom = line_to_plot if zoom_plot else ''
        quantifier.plot_quantified_spectrum(
            plot_title=f"{sample_ID}_#{spectrum_ID}",
            peaks_to_zoom=line_to_zoom,
            annotate_peaks='main'
        )

    total_process_time = (time.time() - sample_processing_time_start) / 60
    print_double_separator()
    logging.info(f"Sample '{sample_ID}' successfully quantified in {total_process_time:.1f} min.")
        
    return quantifier