#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting and quantification of a single X-ray spectrum.

For spectrum-level analysis of fitting and quantification performance.

Import this module in your own code and call the
`fit_and_quantify_spectrum()` function, passing your desired 'sample_ID', 'spectrum_ID' and
options as arguments. This enables integration into larger workflows or pipelines.

Workflow:
    - Loads sample configurations and spectral data from `ledger.json` (primary source)
    - Falls back to legacy `Data.csv` + config JSON only when no ledger exists
    - Performs quantification (optionally only on unquantified spectra)
    - Optionally performs clustering/statistical analysis and saves results

Notes
-----
- Only the `sample_ID` and 'spectrum_ID' are required if acquisition output is saved in the default Results directory;
  otherwise, specify `results_path`.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import os
import warnings
import logging
from pathlib import Path
from typing import Optional, List  

import numpy as np

from autoemx.utils import (
    print_double_separator,
    get_sample_dir,
)
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.config import config_classes_dict
from autoemx.config.ledger_schemas import ClusteringConfig
from .fit_and_quantify_spectrum import fit_and_quantify_spectrum
from autoemx.config import load_sample_ledger

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = [
    "fit_and_quantify_spectrum_from_ledger",
    "fit_and_quantify_spectrum_from_datacsv",
    "fit_and_quantify_spectrum_fromDatacsv",
]

def fit_and_quantify_spectrum_from_ledger(
    sample_ID: str,
    spectrum_ID: int,
    els_sample: Optional[List[str]] = None,
    els_w_frs: Optional[List[str]] = None,
    els_substrate: Optional[List[str]] = None,
    is_standard: bool = False,
    spectrum_lims: Optional[tuple] = None,
    results_path: Optional[str] = None,
    use_instrument_background: bool = dflt.use_instrument_background,
    quantify_plot: bool = True,
    plot_signal: bool = True,
    zoom_plot: bool = False,
    line_to_plot: str = '',
    fit_tol: float = 1e-4,
    is_particle: bool = True,
    max_undetectable_w_fr: float = 0,
    force_single_iteration: bool = False,
    interrupt_fits_bad_spectra: bool = False,
    standards_dict: Optional[dict] = None,
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
        Value reported in the ledger-backed Spectrum ID sequence (legacy-compatible
        with the historical 'Spectrum #' identifier from Data.csv).
    els_sample : list(str), optional
        List of elements in the sample. If the first entry is "" or None, the rest of the list is appended to the 
        list loaded from Comp_analysis_configs.json; otherwise, the provided list replaces it.
    els_substrate : list(str), optional
        List of substrate elements. If the first entry is "" or None, the rest of the list is appended to the 
        list loaded from Comp_analysis_configs.json; otherwise, the provided list replaces it.
    is_standard : bool
        Defines whether measurement is from an experimental standard (i.e., sample of known composition)
    results_path : str, optional
        Base directory where results are stored. Default: autoemx/Results
    use_instrument_background : bool, optional
        Whether to use instrument background if present. Default: False
    quantify_plot : bool, optional
        Whether to quantify the spectrum.
    plot_signal : bool, optional
        Whether to plot the fitted spectrum.
    zoom_plot : bool, optional
        Whether to zoom on a specific line.
    line_to_plot : str, optional
        Line to zoom on.
    fit_tol : float, optional
        scipy fit tolerance. Defines conditions of fit convergence
    is_particle : bool, optional
        If True, treats sample as particle (powder). Uses particle geometry fitting parameters
    max_undetectable_w_fr : float, optional
        Maximum allowed weight fraction for undetectable elements (default: 0). Total mass fraction of fitted
        elements is forced to be between [1-max_undetectable_w_fr, 1]
    force_single_iteration : bool, optional
        If True, quantification will be run for a single iteration only (default: False).
    interrupt_fits_bad_spectra : bool, optional
        If True, interrupt fitting if spectrum is detected to lead to poor quantification (default: False).
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
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)
        
    try:
        sample_dir = get_sample_dir(results_path, sample_ID)
    except Exception as e:
        logging.warning("Failed to get sample directory for %s: %s", sample_ID, e)
        return

    ledger_path = os.path.join(sample_dir, f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
    spectral_info_f_path = ledger_path
    
    print_double_separator()
    logging.info(f"Sample '{sample_ID}', spectrum {spectrum_ID}")
    
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
        if ledger.configs.powder_meas_cfg is not None:
            configs[cnst.POWDER_MEASUREMENT_CFG_KEY] = ledger.configs.powder_meas_cfg
        if ledger.configs.bulk_meas_cfg is not None:
            configs[cnst.BULK_MEASUREMENT_CFG_KEY] = ledger.configs.bulk_meas_cfg
        metadata = {}
    except FileNotFoundError:
        logging.warning(f"Could not find {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return
    except Exception as e:
        logging.warning(f"Error loading {spectral_info_f_path}. Skipping sample '{sample_ID}': {e}")
        return

    # Retrieve configuration objects
    try:
        microscope_cfg      = configs[cnst.MICROSCOPE_CFG_KEY]
        sample_cfg          = configs[cnst.SAMPLE_CFG_KEY]
        measurement_cfg     = configs[cnst.MEASUREMENT_CFG_KEY]
        sample_substrate_cfg= configs[cnst.SAMPLESUBSTRATE_CFG_KEY]
        quant_cfg           = configs[cnst.QUANTIFICATION_CFG_KEY]
        clustering_cfg      = configs.get(cnst.CLUSTERING_CFG_KEY)
    except KeyError as e:
        logging.warning(f"Missing configuration '{e.args[0]}' in {spectral_info_f_path}. Skipping sample '{sample_ID}'.")
        return

    if clustering_cfg is None:
        clustering_cfg = ClusteringConfig()
    
    # Resolve requested spectrum directly from ledger entries.
    target_entry = None
    for idx, spectrum_entry in enumerate(ledger.spectra):
        entry_id = spectrum_entry.spectrum_id
        if entry_id is None:
            if idx == int(spectrum_ID):
                target_entry = spectrum_entry
                break
            continue
        try:
            if int(float(entry_id)) == int(spectrum_ID):
                target_entry = spectrum_entry
                break
        except (TypeError, ValueError):
            if str(entry_id) == str(spectrum_ID):
                target_entry = spectrum_entry
                break

    if target_entry is None:
        logging.warning(
            f"Spectrum ID {spectrum_ID} not found in ledger-backed data for sample '{sample_ID}'."
        )
        return

    spectrum_relpath = target_entry.spectrum_relpath
    if not spectrum_relpath:
        logging.warning(
            f"Missing spectrum pointer path for spectrum ID {spectrum_ID} in sample '{sample_ID}'."
        )
        return

    spectrum_pointer = Path(sample_dir) / spectrum_relpath
    try:
        spectrum = np.asarray(ledger._load_counts_from_pointer_file(spectrum_pointer), dtype=float)
    except Exception as e:
        logging.warning(
            f"Could not load spectrum data for spectrum ID {spectrum_ID} in sample '{sample_ID}': {e}"
        )
        return

    background = None
    if use_instrument_background:
        bkg_relpath = target_entry.instrument_background_relpath
        if bkg_relpath:
            bkg_path = Path(sample_dir) / bkg_relpath
            try:
                background = np.asarray(np.load(bkg_path, allow_pickle=False), dtype=float)
            except Exception:
                background = None
        if background is None:
            warnings.warn(
                "Instrument background not found or empty for this spectrum. "
                "Spectral background will be computed instead."
            )

    sp_collection_time = target_entry.live_acquisition_time

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
    
    # Sample elements
    if els_sample is None:
        els_sample = el_to_quantify
    elif (els_sample[0] == "" or els_sample[0] is None):
        els_sample = el_to_quantify + els_sample[1:]
            
    # Substrate elements
    if els_substrate is None:
        els_substrate = sample_substrate_cfg.elements
    elif (els_substrate[0] == "" or els_substrate[0] is None):
        els_substrate = sample_substrate_cfg.elements + els_substrate[1:]
        
    # Spectral limits
    if spectrum_lims is None:
        spectrum_lims = quant_cfg.spectrum_lims
    
    # Ensure spectrum_lims are integers for array slicing
    if spectrum_lims is not None:
        try:
            spectrum_lims = (int(spectrum_lims[0]), int(spectrum_lims[1]))
        except (TypeError, ValueError, IndexError):
            logging.error(f"Invalid spectrum_lims: {spectrum_lims}")
            return
    else:
        logging.error("spectrum_lims could not be determined")
        return

    # Keep static typing explicit for downstream function signature.
    if els_sample is None:
        els_sample = []
    if els_substrate is None:
        els_substrate = []

    if els_w_frs is None and is_standard:
        if hasattr(sample_cfg, 'w_frs') and sample_cfg.w_frs is not None:
            els_w_frs = sample_cfg.w_frs
        else:
            raise ValueError(
                "No 'w_frs' (element weight fractions) found in sample_cfg for standard sample. "
                "You must pass 'els_w_frs' directly, or modify the ledger to include 'w_frs'."
            )
            
    quantifier = fit_and_quantify_spectrum(
        spectrum_vals = spectrum,
        spectrum_lims = spectrum_lims,
        microscope_ID = microscope_cfg.ID,
        meas_type = measurement_cfg.type,
        meas_mode = measurement_cfg.mode,
        det_ch_offset=offset,
        det_ch_width=scale,
        beam_energy = beam_energy,
        emergence_angle = emergence_angle,
        sp_collection_time = sp_collection_time,
        sample_ID = sample_ID,
        els_sample = els_sample,
        els_w_frs = els_w_frs,
        els_substrate = els_substrate,
        background_vals=background,
        fit_tol = fit_tol,
        is_particle = is_particle,
        quantify_plot = quantify_plot,
        max_undetectable_w_fr = max_undetectable_w_fr,
        force_single_iteration = force_single_iteration,
        interrupt_fits_bad_spectra = interrupt_fits_bad_spectra,
        standards_dict = standards_dict,
        plot_signal = plot_signal,
        plot_title = f"{sample_ID}_#{spectrum_ID}",
        zoom_plot = zoom_plot,
        line_to_plot = line_to_plot,
        print_results = print_results,
        quant_verbose = quant_verbose,
        fitting_verbose = fitting_verbose
    )

    return quantifier


def fit_and_quantify_spectrum_fromDatacsv(*args, **kwargs):
    """Deprecated camelCase alias; use fit_and_quantify_spectrum_from_ledger."""
    warnings.warn(
        "fit_and_quantify_spectrum_fromDatacsv is deprecated; "
        "use fit_and_quantify_spectrum_from_ledger instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fit_and_quantify_spectrum_from_ledger(*args, **kwargs)


def fit_and_quantify_spectrum_from_datacsv(*args, **kwargs):
    """Deprecated snake_case alias; use fit_and_quantify_spectrum_from_ledger."""
    warnings.warn(
        "fit_and_quantify_spectrum_from_datacsv is deprecated; "
        "use fit_and_quantify_spectrum_from_ledger instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return fit_and_quantify_spectrum_from_ledger(*args, **kwargs)