#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting and quantification of a single X-ray spectrum.

For spectrum-level analysis of fitting and quantification performance.

Created on Tue Jul 29 13:18:16 2025

@author: Andrea
"""

import time
import logging
import warnings
from typing import Optional, Any

from autoemx.utils import print_double_separator
from autoemx.core.quantifier import XSp_Quantifier
from autoemx.core.quantifier.peak_overlaps import warn_peak_overlaps
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.utils import (
    print_double_separator,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["fit_and_quantify_spectrum"]

def fit_and_quantify_spectrum(
    spectrum_vals,
    spectrum_lims,
    microscope_ID,
    meas_type,
    meas_mode,
    det_ch_offset,
    det_ch_width,
    beam_energy,
    emergence_angle,
    sp_collection_time = None,
    sample_ID = '',
    els_sample: Optional[list] = None,
    els_w_frs: Optional[list] = None,
    els_substrate: Optional[list] = None,
    background_vals=None,
    fit_tol: float = 1e-4,
    is_particle: bool = False,
    quantify: bool = True,
    max_undetectable_w_fr: float = 0,
    force_single_iteration: bool = False,
    interrupt_fits_bad_spectra: bool = False,
    standards_dict = None,
    plot_signal: bool = True,
    plot_title = '',
    zoom_plot: bool = False,
    line_to_plot: str = '',
    print_results: bool = True,
    quant_verbose: bool = True,
    fitting_verbose: bool = True,
    free_area_el_lines=None,
    peaks_to_annotate: str = 'main',
    identify_missing_elements: bool = dflt.identify_missing_elements,
    quantify_plot: Optional[bool] = None,
):
    """
    Fit and (optionally) quantify a single spectrum.

    Parameters
    ----------
    spectrum_vals : numpy.ndarray
        The measured EDS spectrum (counts per channel).
    spectrum_lims : tuple of int
        Tuple specifying the start and end indices for the spectrum region to analyze.
    microscope_ID : str
        Microscope identifier for detector calibration data.
    meas_type : str
        Measurement type (e.g., 'EDS').
    meas_mode : str
        Measurement mode, defining detector calibrations and (optionally) beam current (e.g., 'point').
    det_ch_offset : float
        Detector channel energy offset (keV).
    det_ch_width : float
        Detector channel width (keV).
    beam_energy : float
        Electron beam energy (keV).
    emergence_angle : float
        Detector emergence (take-off) angle (degrees).
    sp_collection_time : float or None
        Live time of spectrum acquisition (in seconds).
    sample_ID : str
        Sample identifier.
    els_sample : list, optional
        List of elements in the sample.
    els_substrate : list, optional
        List of substrate elements.
    background_vals : array-like or None, optional
        Background spectrum to subtract. If None, the background will be modeled during fitting.
    fit_tol : float, optional
        scipy fit tolerance. Defines conditions of fit convergence
    is_particle : bool, optional
        If True, treats sample as particle (powder). Uses particle geometry fitting parameters
    quantify : bool, optional
        If True (default), fit and quantify the spectrum. If False, only fit it (e.g. to measure the P/B ratios of a standard).
    max_undetectable_w_fr : float, optional
        Maximum allowed weight fraction for undetectable elements (default: 0). Total mass fraction of fitted
        elements is forced to be between [1-max_undetectable_w_fr, 1]
    force_single_iteration : bool, optional
        If True, quantification will be run for a single iteration only (default: False).
    interrupt_fits_bad_spectra : bool, optional
        If True, interrupt fitting if spectrum is detected to lead to poor quantification (default: False).
    standards_dict: dict, optional
        Dictionary of standard values to use for quantification.
    plot_signal : bool, optional
        Whether to plot the fitted spectrum.
    plot_title : str
        String printed as plot title.
    zoom_plot : bool, optional
        Whether to zoom on a specific line.
    line_to_plot : str, optional
        Line to zoom on.
    peaks_to_annotate : str, optional
        Which peaks to annotate in `plot_quantified_spectrum` ('all', 'main', 'most).
    identify_missing_elements : bool, optional
        If True, look for elements missing from the fit before fitting (see XSp_Quantifier). Elements confirmed
        by the fit are added (and quantified if a standard is available); weaker detections are reported.
        Default: True (defaults.identify_missing_elements).
    print_results : bool, optional
        If True, prints all fitted parameters and their values (default: True).
    quant_verbose : bool, optional
        If True, prints quantification operations
    fitting_verbose : bool, optional
        If True, prints fitting operations
    quantify_plot : bool, optional
        Deprecated, use `quantify`.
        
    Returns
    -------
    quantifier : XSp_Quantifier
        The quantifier object containing the results, fit parameters, and methods for further analysis and plotting.
    """
    if quantify_plot is not None:
        warnings.warn("'quantify_plot' is deprecated, use 'quantify' instead.", FutureWarning, stacklevel=2)
        quantify = quantify_plot
    sample_processing_time_start = time.time()
    
    # Quantification
    quantifier = XSp_Quantifier(
        spectrum_vals = spectrum_vals,
        spectrum_lims = spectrum_lims,
        microscope_ID = microscope_ID,
        meas_type = meas_type,
        meas_mode = meas_mode,
        det_ch_offset=det_ch_offset,
        det_ch_width=det_ch_width,
        beam_e=beam_energy,
        emergence_angle=emergence_angle,
        background_vals=background_vals,
        els_sample=els_sample,
        els_substrate=els_substrate,
        els_w_fr=els_w_frs,
        is_particle=is_particle,
        sp_collection_time=sp_collection_time,
        max_undetectable_w_fr=max_undetectable_w_fr,
        fit_tol=fit_tol,
        verbose=quant_verbose,
        fitting_verbose=fitting_verbose,
        standards_dict = standards_dict,
        free_area_el_lines=free_area_el_lines,
        identify_missing_elements=identify_missing_elements,
    )

    def _warn_peak_overlaps():
        warn_peak_overlaps(
            meas_type,
            quantifier.els_to_quantify,
            quantifier.els_substrate,
            beam_energy_keV=beam_energy,
            energy_range_keV=(quantifier.energy_vals[0], quantifier.energy_vals[-1]),
            microscope_ID=microscope_ID,
            free_area_el_lines=free_area_el_lines,
        )

    if quantify:
        _warn_peak_overlaps()

    try:
        if quantify:
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

    # Measured P/B ratios are reported when fitting a standard (known composition) without quantifying it
    if els_w_frs is not None and not quantify:
        fitted_peaks_info = getattr(quantifier, 'fitted_peaks_info', None) or {}
        pb_lines = XSp_Quantifier.format_measured_pb(
            (
                # Escape/pileup peaks (keys like 'Fe_Ka1_esc') are skipped
                (*el_line.split('_'), peak_info.get(cnst.PB_RATIO_KEY))
                for el_line, peak_info in fitted_peaks_info.items()
                if len(el_line.split('_')) == 2
            ),
            standard_elements=els_w_frs,
        )
        if not fitted_peaks_info:
            message = "No fitted peaks available — cannot report Measured PB."
        elif pb_lines:
            message = "\n".join(["Measured PB:"] + [f"  {pb_line}" for pb_line in pb_lines])
        else:
            message = "Fit completed (no reportable P/B peaks)."

        print_double_separator(message=message)

    if plot_signal:
        line_to_zoom = line_to_plot if zoom_plot else ''
        quantifier.plot_quantified_spectrum(
            plot_title=plot_title,
            peaks_to_zoom=line_to_zoom,
            annotate_peaks=peaks_to_annotate,
        )

    total_process_time = (time.time() - sample_processing_time_start)
    print_double_separator()
    time_str = f"{total_process_time/60:.1f} min" if total_process_time > 100 else f"{total_process_time:.1f} sec"
    quant_str = 'quantified' if quantify else 'fitted'
    logging.info(f"Sample '{sample_ID}' successfully {quant_str} in {time_str}.")

    if quantify:
        _warn_peak_overlaps()

    return quantifier