#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
X-ray Spectrum Quantification Module

Created on Thu Jun 27 14:23:22 2024

@author: Andrea Giunto

This module provides classes and functions for matrix correction and quantitative analysis
of X-ray spectra, implementing the peak-to-background (P/B) method and related ZAF corrections
for scanning electron microscopy (SEM) energy-dispersive X-ray spectroscopy (EDS).

Class Structure and Interactions
-------------------------------
The main classes are:
- **XSp_Quantifier**
  Performs quantification of X-ray spectra using matrix corrections. Interfaces with spectral fitting
  routines from XSp_Fitter to extract peak/background intensities and applies the full suite
  of correction factors from Quant_Corrections to obtain quantitative elemental concentrations.
  Additionally provides functions to plot and print the results.
  
- **Quant_Corrections**
  Provides methods for calculating matrix correction factors (Z, A, R) for the P/B method.
  This class is instanced by XSp_Quantifier to calculate the matrix correction factors.

Typical Usage
-------------
1. **Initialize the quantifier object:**
    quantifier = XSp_Quantifier(...)  # See class docs for initialization

2. **Quantify a spectrum:**
    quant_result = quantifier.quantify_spectrum()
        quant_result contains quantified atomic and weight fractions, analytical error, and spectral
        fitting metrics (reduced chi-square, and R-squared).

3. **Print quantification results:**
    quantifier.print_quant_result(quant_result)

3. **Plot spectrum:**
    quantifier.plot_quantified_spectrum()

Customization & Calibration
---------------------------
Detector calibration, physical constants, and mass absorption coefficients are handled via the calibs module and supporting utility functions.
Users should calibrate the values in the calibs module corresponding to their microscope and EDS settings prior usage. 

Dependencies
------------
numpy, scipy, sympy, pandas
lmfit, pymatgen.core.Element (from supporting libraries)
XSp_Fitter (required for spectral fitting; see separate module)
calibs, lib modules

How the classes interact
------------------------
XSp_Quantifier is the main user-facing class for quantification. It uses Quant_Corrections to compute all correction factors and matrix effects,
and relies on the XSp_Fitter module to extract peak/background intensities from measured spectra.
Quant_Corrections provides the core physical models for all matrix corrections, and is initialized with the relevant sample and measurement parameters.
The XSp_Fitter module must be installed and available for spectral fitting and deconvolution.

"""
# =============================================================================
# Standard library imports
# =============================================================================
import os
import re
import warnings
import traceback
from typing import Any, Optional, Dict, Tuple, Sequence, List, Union

# =============================================================================
# Third-party library imports
# =============================================================================
import numpy as np
import pandas as pd
import sympy as sp
import matplotlib.pyplot as plt
import lmfit
from pymatgen.core import Element

# =============================================================================
# Local application/library imports
# =============================================================================
import autoemx.calibrations as calibs 
import autoemx.utils.constants as cnst
from autoemx.utils.helper import (
    print_nice_1d_row,
    print_single_separator,
    print_double_separator,
    EDSError,
    weight_to_atomic_fr,
    atomic_to_weight_fr
)
from autoemx.data.Xray_lines import get_el_xray_lines
from autoemx.data.Xray_absorption_coeffs import xray_mass_absorption_coeff
from autoemx.data.mean_ionization_potentials import J_df
from autoemx.core.fitter import (
    XSp_Fitter,
    Background_Model,
    Peaks_Model
)
from autoemx.config.ledger_schemas import (
    FitResult as QuantFitResult,
    FittedPeakResult,
    QuantificationDiagnostics,
    QuantificationResult,
)
from .corrections import Quant_Corrections

from autoemx._logging import get_logger
logger = get_logger(__name__)

#%% XSp_Quantifier class
class XSp_Quantifier:
    """
    Class for quantitative analysis of EDS spectra.

    Attributes
    ----------
    spectrum_vals : numpy.ndarray
        The measured EDS spectrum (counts per channel).
    energy_vals : numpy.ndarray
        Energy scale corresponding to `spectrum_vals`.
    spectrum_lims : tuple of int
        Start and end indices for the spectrum region to analyze.
    fit_background : bool
        Whether to fit a background model. True, if background_vals are not provided.
    background_vals : numpy.ndarray
        Fitted or provided background spectrum counts.
    els_sample : list of str
        All elements present in the sample, including undetectable elements.
    els_to_quantify : list of str
        Elements to be quantified (excluding undetectable).
    els_substrate : list of str
        Elements present in the substrate or sample holder appearing in spectra (excluding undetectable).
    els_w_fr : dict or None, optional
        Dictionary of fixed elemental mass fractions, by element symbol, e.g. {'Si':0.33, 'O':0.67}.
    max_undetectable_w_fr : float
        Maximum allowed total mass fraction for undetectable elements during fitting.
    force_total_w_fr : bool
        Whether total mass fraction is normalized 1 during fitting.
        False, if undetectable elements are present in the sample. True otherwise.
    is_particle : bool
        Whether the sample is a particle (affects fitting and quantification corrections).
    sp_collection_time : float or None
        Live time of spectrum acquisition (in seconds).
    fit_tol : float
        Tolerance for spectrum fitting convergence.
    bad_quant_flag : int or None
        Flag indicating quantification issues, or None if successful.
    microscope_ID : str
        Microscope identifier for calibration.
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
    verbose : bool
        Print information during quantification.
    fitting_verbose : bool
        If True, print detailed information during each fitting step.
            
    Class attributes
    ----------------
    xray_quant_ref_lines : list(str)
        List of X-ray lines used as reference
    """
    #  Reference lines for quantification
    xray_quant_ref_lines = ['Ka1', 'La1', 'Ma1', 'Mz1']
    
    def __init__(
        self,
        spectrum_vals,
        spectrum_lims,
        microscope_ID,
        meas_type,
        meas_mode,
        det_ch_offset,
        det_ch_width,
        beam_e,
        emergence_angle,
        energy_vals=None,
        background_vals=None,
        els_sample=None,
        els_substrate=None,
        els_w_fr=None,
        is_particle=False,
        sp_collection_time=None,
        max_undetectable_w_fr=0.10,
        fit_tol=1e-4,
        standards_dict = None,
        verbose=False,
        fitting_verbose=False,
        free_area_el_lines=None,
    ):
        """
        Initialize an XSp_Quantifier for quantitative EDS spectrum analysis.
    
        Parameters
        ----------
        spectrum_vals : array-like
            The measured EDS spectrum (counts per channel).
        spectrum_lims : tuple of int
            Tuple specifying the start and end indices for the spectrum region to analyze.
        microscope_ID : str
            Microscope identifier for calibration (e.g., 'PhenomXL').
        meas_type : str
            Measurement type (e.g., 'EDS').
        meas_mode : str
            Measurement mode, defining detector calibrations and (optionally) beam current (e.g., 'point').
        det_ch_offset : float
            Detector channel energy offset (keV).
        det_ch_width : float
            Detector channel width (keV).
        beam_e : float
            Electron beam energy (in keV).
        emergence_angle : float
            Detector emergence (take-off) angle in degrees.
        energy_vals : array-like or None, optional
            The energy scale corresponding to `spectrum_vals`. If None, it will be calculated from
            detector calibration parameters (det_ch_offset, det_ch_width).
        background_vals : array-like or None, optional
            Background spectrum to subtract. If None, the background will be modeled during fitting.
        els_sample : list of str or None, optional
            List of element symbols present in the sample (including those to quantify).
        els_substrate : list of str or None, optional
            List of element symbols present in the substrate or sample holder.
        els_w_fr : dict or None, optional
            Dictionary of fixed elemental mass fractions, by element symbol, e.g. {'Si':0.33, 'O':0.67}.
        is_particle : bool, optional
            Set to True if the sample is a particle (affects quantification corrections).
        sp_collection_time : float or None, optional
            Live time of spectrum acquisition (in seconds).
        max_undetectable_w_fr : float, optional
            Maximum total mass fraction allowed for undetectable elements (default: 0.10).
        fit_tol : float, optional
            Tolerance for spectrum fitting convergence.
        standards_dict: dict, optional
            Dictionary of standard values to use for quantification.
        verbose : bool, optional
            If True, print information during quantification.
        fitting_verbose : bool, optional
            If True, print detailed information during each fitting step.
    
        Notes
        -----
        - If `background_vals` is not provided, a background model will be fitted and subtracted.
        - If `energy_vals` is not provided, the energy scale is calculated from detector calibration parameters.
        - If undetectable elements are present in `els_sample`, normalization of mass fractions is relaxed.
        - Calibration data for the specified microscope and EDS mode is loaded at initialization.
        - The class stores all relevant spectrum and calibration data as attributes for downstream quantification routines.
        """
        # Handle mutable default arguments
        if els_sample is None:
            els_sample = []
        if els_substrate is None:
            els_substrate = ['C', 'O', 'Al']
        if els_w_fr is None:
            els_w_fr = {}

        # Load microscope calibrations for this instrument and mode
        calibs.load_microscope_calibrations(microscope_ID, meas_mode, load_detector_channel_params=False)
        
        # EDS and instrument parameters
        self.microscope_ID = microscope_ID
        self.meas_mode = meas_mode
        self.meas_type = meas_type
        
        # Not loaded from calibs because these values are constantly re-calibrated
        # and may be different from current values for previously collected spectra
        self.det_ch_offset = det_ch_offset
        self.det_ch_width = det_ch_width
        
        self.beam_energy = beam_e
        self.emergence_angle = emergence_angle


        # Store original total counts before spectrum limits are applied
        self.tot_sp_counts = sum(spectrum_vals)
        self.sp_collection_time = sp_collection_time

        # Background handling and spectrum slicing
        self.spectrum_lims = spectrum_lims
        sp_start, sp_end = spectrum_lims
        if background_vals is None:
            self.fit_background = True
            self.spectrum_vals = np.array(spectrum_vals)[sp_start:sp_end]
        else:
            self.fit_background = False
            self.background_vals = np.array(background_vals)[sp_start:sp_end]
            self.spectrum_vals = (np.array(spectrum_vals) - np.array(background_vals))[sp_start:sp_end]

        # Energy values
        if energy_vals is not None:
            self.energy_vals = np.array(energy_vals)[sp_start:sp_end]
        elif det_ch_offset is not None and det_ch_width is not None:
            self.energy_vals = np.array([det_ch_offset + det_ch_width * i for i in range(sp_start, sp_end)])
        else:
            raise ValueError("Need to provide an array of energy values, or the detector bin width and energy offset.")
            
        # Sample characteristics
        self.is_particle = is_particle

        # Elements to include in the analysis
        self.els_sample = list(els_sample)  # All elements present in the sample, including undetectable
        self.els_to_quantify = [el for el in els_sample if el not in calibs.undetectable_els]
        self.els_w_fr = {el: w_fr for el, w_fr in els_w_fr.items() if el not in calibs.undetectable_els}
        self.els_substrate = [el for el in els_substrate if el not in calibs.undetectable_els]

        # Fitting parameters
        self.fit_tol = fit_tol
        self.force_total_w_fr = not any(el in calibs.undetectable_els for el in self.els_sample)
        self.max_undetectable_w_fr = max_undetectable_w_fr

        # Standards
        self.standards = standards_dict  # If None, it is loaded during quantification
        self.bad_quant_flag = None # Initialise, for spectra that are not quantified
        self.iterations_run: Optional[int] = None
        self.quant_converged: Optional[bool] = None
        self.ref_lines_for_quant: List[str] = []
        self.missing_reference_peaks: List[str] = []

        self.verbose = verbose
        self.fitting_verbose = fitting_verbose
        self.free_area_el_lines = free_area_el_lines
        
    #%% Fit spectrum
    # =============================================================================
    def _initialize_spectrum_fitter(self) -> None:
        """
        Initialize the XSp_Fitter instance with the current spectrum and settings.
    
        This method sets up the XSp_Fitter as `self.fitter` using the instance's
        spectrum data and configuration attributes. No fitting is performed.
    
        Parameters
        ----------
        None
    
        Returns
        -------
        None
    
        Notes
        -----
        This method only initializes the fitter. To perform fitting, call `self._fit_spectrum()`.
        """
        self.fitter = XSp_Fitter(
            spectrum_vals=self.spectrum_vals,
            energy_vals=self.energy_vals,
            spectrum_lims=self.spectrum_lims,
            microscope_ID = self.microscope_ID,
            meas_mode=self.meas_mode,
            det_ch_offset=self.det_ch_offset,
            det_ch_width=self.det_ch_width,
            beam_e=self.beam_energy,
            emergence_angle=self.emergence_angle,
            fit_background=self.fit_background,
            is_particle=self.is_particle,
            els_to_quantify=self.els_to_quantify,
            els_substrate=self.els_substrate,
            els_w_fr=self.els_w_fr,
            force_fr_total=self.force_total_w_fr,
            tot_sp_counts=self.tot_sp_counts,
            sp_collection_time=self.sp_collection_time,
            xray_quant_ref_lines=self.xray_quant_ref_lines,
            print_evolving_params=False,
            verbose=self.fitting_verbose,
            free_area_el_lines=self.free_area_el_lines,
        )
    
    
    def initialize_and_fit_spectrum(
        self, 
        params: Optional[lmfit.Parameters] = None,
        print_results: Optional[bool] = False
    ) -> None:
        """
        Perform a complete fit of the spectrum provided to the instance of XSp_Quantifier.
        This method is intended for single-iteration fits where sample composition is known and constrained
        (e.g., EDS standards).
    
        The fit results and all relevant fitted parameters are stored as instance attributes.
    
        Parameters
        ----------
        params : Optional[lmfit.Parameters], optional
            Parameters object to pass to the spectrum fitter. If None, default parameters are used.
        print_result : bool, optional
            If True, prints the fitting results.
    
        Returns
        -------
        None
    
        Nested Calls
        ------------
        Calls the following methods to further store values extracted from the fit:
            - self._initialize_spectrum_fitter(): initializes the spectrum fitter
            - self._fit_spectrum(): fits the spectrum and stores the fit results
        """
        # Get initial value of K. Must be done before initialising fitter
        initial_par_vals: Optional[dict] = None
        if self.is_particle and self.fit_background:
            K_val = self.get_starting_K_val()
            if K_val is not None:
                initial_par_vals = {'K': K_val}
    
        # Initialize the fitter (no fitting performed yet)
        self._initialize_spectrum_fitter()

        # Now perform the fit and store results
        self._fit_spectrum(
            params=params,
            initial_par_vals=initial_par_vals,
            f_tol=self.fit_tol,
            print_result=print_results,
            print_extended_result=self.fitting_verbose
        )
        
        bad_fit_flag = self._check_if_unreliable_quant(
            iter_cntr = 1, analytical_er = 0, interrupt_fits_bad_spectra = False
        )
        
        return bad_fit_flag
        

    def _fit_spectrum(
        self,
        params: Optional[lmfit.Parameters] = None,
        initial_par_vals: Optional[Dict[str, float]] = None,
        f_tol: float = 1e-4,
        n_iter: Optional[int] = None,
        print_result: bool = True,
        print_extended_result: bool = False
    ) -> Tuple[lmfit.Parameters, np.ndarray, float, float]:
        """
        Perform a single spectrum fitting iteration using the already-initialized
        XSp_Fitter instance (i.e., self.fitter).
    
        This method fits the spectrum, updates relevant instance attributes with the results,
        and returns key quantitative outputs.
    
        Parameters
        ----------
        params : Optional[lmfit.Parameters], optional
            Parameters object to pass to the spectrum fitter. If None, default parameters are used.
        initial_par_vals : Optional[Dict[str, float]], optional
            Initial parameter values for the fit. If None, default initial values are used.
        f_tol : float, optional
            Function tolerance for the fitting algorithm (default is 1e-4).
        n_iter : int, optional
            Iteration number (for display purposes only).
        print_result : bool, optional
            If True, print the fitting result summary.
        print_extended_result : bool, optional
            If True and print_result = True, print extended fitting results.
    
        Returns
        -------
        None
    
        Side Effects
        ------------
        Updates the following instance attributes:
            - self.fit_result: lmfit.ModelResult object from the spectrum fitting.
            - self.fitted_els: List of elements for which at least one peak was fitted.
            - self.fitted_els_quant: List of elements to quantify that were successfully fitted.
            - self.fitted_xray_lines: Information on fitted X-ray lines.
            - self.fit_components: Dictionary of fitted spectral components and their values.
    
        Nested Calls
        ------------
        Calls the following methods to further store values extracted from the fit:
            - self._store_background_vals() (if background is fitted)
            - self._assemble_peaks_info()
        """
        fit_result, fitted_lines = self.fitter.fit_spectrum(
            parameters=params,
            initial_par_vals=initial_par_vals,
            function_tolerance=f_tol,
            n_iter=n_iter,
            print_result=print_result,
            print_result_extended=print_extended_result
        )
    
        fit_components = fit_result.eval_components(x=self.energy_vals)
        self.fit_components = fit_components
    
        self.fit_result = fit_result
        self.fitted_els = self.fitter.fitted_els
        self.fitted_els_quant = [el for el in self.els_to_quantify if el in self.fitted_els]
        self.fitted_xray_lines = fitted_lines
    
        if self.fit_background:
            self._store_background_vals(fit_result, fit_components)
    
        self._assemble_peaks_info()
    
    
    def get_starting_K_val(self) -> Optional[float]:
        """
        Estimate the initial background scaling factor K by fitting the high-energy portion of the spectrum.
    
        This method fits only the region of the spectrum above a threshold energy to avoid regions heavily affected by absorption.
        It is intended to provide an optimal starting value for K, especially for particle spectra, to prevent the algorithm from
        compensating for background intensity via particle geometry parameters, which are only intended for background shape fitting.
    
        Returns
        -------
        K_val : Optional[float]
            Estimated initial value for the background scaling factor K,
            or None if a suitable value could not be determined.
    
        Notes
        -----
        - Uses XSp_Fitter with particle geometry disabled for this quick fit.
        - Prints diagnostic information if `self.verbose` is True.
        """
        K_val = None
    
        # Determine energy threshold based on beam energy
        if self.beam_energy > 10:  # keV
            en_thresh = 5  # keV
        elif self.beam_energy >= 7.5:  # keV
            en_thresh = 3  # keV
        else:
            if self.verbose:
                logger.warning("⚠️ Initial background scaling factor K could not be estimated.")
                logger.warning(f"⚠️ Current beam energy of {self.beam_energy} keV too low for reliable high-energy background fitting.")
                logger.warning("⚠️ Beam energy needs to be at least 7.5 keV.")
            return None
    
        high_energy_indices = self.energy_vals > en_thresh
    
        if self.fitting_verbose:
            print_double_separator()
            logger.info(f"🔬 Fit of spectrum above {en_thresh} keV to get initial background scaling factor K...")
            logger.info("ℹ️ Turned off particle morphology parameters to avoid affecting value of K.")
    
        if any(high_energy_indices):
            energy_vals = self.energy_vals[high_energy_indices]
            spectrum_vals = self.spectrum_vals[high_energy_indices]
            low_en_spectrum_lim = self.spectrum_lims[0] + np.argmax(high_energy_indices)
    
            # Initialize XSp_Fitter without particle geometry
            fitter = XSp_Fitter(
                spectrum_vals=spectrum_vals,
                energy_vals=energy_vals,
                spectrum_lims=(low_en_spectrum_lim, self.spectrum_lims[-1]),
                microscope_ID=self.microscope_ID,
                meas_mode=self.meas_mode,
                det_ch_offset=self.det_ch_offset,
                det_ch_width=self.det_ch_width,
                beam_e=self.beam_energy,
                emergence_angle=self.emergence_angle,
                fit_background=self.fit_background,
                is_particle=False,
                els_to_quantify=self.els_to_quantify,
                els_substrate=self.els_substrate,
                els_w_fr=self.els_w_fr,
                force_fr_total=False,
                tot_sp_counts=self.tot_sp_counts,
                sp_collection_time=self.sp_collection_time,
                xray_quant_ref_lines=self.xray_quant_ref_lines,
                print_evolving_params=False,
                verbose=False
            )
    
            try:
                fit_result, _ = fitter.fit_spectrum(function_tolerance=1e-5)
                # If the fit is good, extract K
                if fit_result.redchi < self.tot_sp_counts / 1000:
                    K_val = fit_result.params['K'].value
                if self.verbose:
                    if K_val is not None:
                        logger.info(f"✅ Found K = {K_val:.2f}")
                    else:
                        logger.warning("⚠️ Failed to find initial K value")
            except Exception as e:
                if self.verbose:
                    logger.error("❌ An error occurred during the quick background fit for K estimation:")
                    logger.error(f"❌ {type(e).__name__}: {e}")
                K_val = None
        else:
            if self.verbose:
                logger.warning("⚠️ No suitable high-energy data for K estimation.")
    
        return K_val
    
    
    def _store_background_vals(self, fit_result, fit_components):
        """
        Stores the evaluated background values (with and without detector response) for the current fit.
        
        Used to exctract the values of background counts below reference peaks for quantification.
    
        This method:
        - Calculates and stores the background values using the original energy grid.
        - Re-initializes the background model to clear any previous absorption attenuation.
        - Evaluates and stores the background values on a finer energy grid, with detector response disabled.
    
        Parameters
        ----------
        fit_result : lmfit.model.ModelResult
            Result object from the spectrum fit.
        fit_components : dict
            Dictionary of evaluated fit components from the fit.
        """
    
        # Find relevant component names for absorption attenuation and generated background
        abs_att_param_name = [s for s in fit_components.keys() if '_abs_att' in s][0]
        gen_bckgrnd_param_name = [s for s in fit_components.keys() if '_generated_b' in s][0]
        bcksctr_param_name = [s for s in fit_components.keys() if '_backscattering_corr' in s][0]
        sp_param_name = [s for s in fit_components.keys() if '_stopping_p' in s][0]
        det_eff_param_name = '_det_efficiency'
        
        # Store background values evaluated on the original energy grid
        self.background_vals = (
            fit_components[gen_bckgrnd_param_name]
            * fit_components[abs_att_param_name]
            * fit_components[det_eff_param_name]
            * fit_components[bcksctr_param_name]
            * fit_components[sp_param_name]
        )
    
        # Define a finer energy grid (0.5 eV step)
        deltaE_finer = 0.5e-3  # 0.5 eV in keV
        energy_vals_finer = np.arange(self.energy_vals[0], self.energy_vals[-1] + deltaE_finer, deltaE_finer)
        self.energy_vals_finer = energy_vals_finer
    
        # Prepare fit parameters with detector response disabled
        params_wo_det_response = fit_result.params.copy()
        params_wo_det_response['apply_det_response'].value = 0  # Disable detector response convolution
    
        # Evaluate fit components on the finer energy grid without detector response
        Background_Model._clear_cached_abs_att_variables() # Clear cached absorption attenuation values
        fit_components_wo_det_response = fit_result.eval_components(
            params=params_wo_det_response, x=energy_vals_finer
        )
        self.background_vals_wo_det_response = (
            fit_components_wo_det_response[gen_bckgrnd_param_name]
            * fit_components_wo_det_response[abs_att_param_name]
            * fit_components_wo_det_response[det_eff_param_name]
            * fit_components_wo_det_response[bcksctr_param_name]
            * fit_components_wo_det_response[sp_param_name]
        )


    #%% Extract information from fitted peaks, including measured peak-to-background ratios
    # =============================================================================
    def _get_peak_info(self, el_line: str) -> Tuple[float, float, float, float, float, float, float, float]:
        """
        Returns fitted Gaussian parameters and peak/background ratio for a given characteristic X-ray line.
    
        Parameters
        ----------
        el_line : str
            Characteristic X-ray line (e.g., 'Si_Ka', 'Mn_La', 'Fe_Ka_esc').
    
        Returns
        -------
        area : float
            Area parameter of Gaussian peak.
        sigma : float
            Gaussian sigma of the fitted peak.
        center : float
            Fitted peak position.
        th_energy : float
            Theoretical X-ray energy of the line (in keV).
        height : float
            Peak height.
        PB_ratio : float
            Peak-to-background ratio (cnts/ (cnts/eV)).
        peak_int : float
            Integrated peak intensity (sum of counts over fitted peak component).
        bckgrnd_int : float
            Interpolated background intensity at the peak energy (cnts/eV).
        """
        # Parse element and line
        el_line_str_components = el_line.split("_")
        el, line = el_line_str_components[:2]
        # Determine theoretical energy
        if len(el_line_str_components) == 3:
            if 'esc' in el_line_str_components[2]:  # Escape peak
                th_energy = get_el_xray_lines(el)[line]['energy (keV)'] - 1.74
            elif 'pileup' in el_line_str_components[2]:  # Pileup peak
                th_energy = get_el_xray_lines(el)[line]['energy (keV)'] * 2
            else:
                th_energy = get_el_xray_lines(el)[line]['energy (keV)']
        else:
            th_energy = get_el_xray_lines(el)[line]['energy (keV)']
    
        area = self.fit_result.params[el_line + '_area'].value
        is_peak_absent = False
    
        if area != 0:
            sigma = self.fit_result.params[el_line + '_sigma'].value
            center = self.fit_result.params[el_line + '_center'].value
    
            if self.energy_vals[0] < th_energy < self.energy_vals[-1]:
                peak_int = np.sum(self.fit_components[el_line + "_"])
                if self.fit_background:
                    bckgrnd_int = np.interp(th_energy, self.energy_vals_finer, self.background_vals_wo_det_response)
                else:
                    bckgrnd_int = np.interp(th_energy, self.energy_vals, self.background_vals)
                bckgrnd_int /= (self.det_ch_width * 1000)  # Normalize to cnts/eV so that it's tranferable across detectors with different bin_widths
                PB_ratio = peak_int / bckgrnd_int if bckgrnd_int != 0 else 0
                height = np.max(self.fit_components[el_line + "_"])
            else:
                is_peak_absent = True
        else:
            is_peak_absent = True
    
        if is_peak_absent:
            sigma = 0.0
            height = 0.0
            PB_ratio = 0.0
            peak_int = 0.0
            bckgrnd_int = 0.0
            center = th_energy
    
        return area, sigma, center, th_energy, height, PB_ratio, peak_int, bckgrnd_int


    def _assemble_peaks_info(self) -> None:
        """
        Stores information for all fitted peaks using _get_peak_info and
        generates list of quantified elements and corresponding characteristic lines
    
        Returns
        -------
        None
    
        Side Effects
        ------------
        Updates the following instance attributes:
            - self.fitted_peaks_info : dict
                For each el_line, a dict containing 'area', 'sigma', 'center', 'fwhm',
                'peak_intensity', 'background_intensity', 'th_energy', 'height', 'PB_ratio'.
            - self.ref_lines_for_quant : list
                List of X-ray lines used for quantification of each element.
            - self.fitted_els_quant : list
                List of elements that can be quantified in this spectrum.
        """
        fitted_peaks_info = {}
        for el_line in self.fitted_xray_lines:
            area, sigma, center, th_energy, height, PB_ratio, peak_int, bckgrnd_int = self._get_peak_info(el_line)
            fwhm = 2.355 * sigma
            fitted_peaks_info[el_line] = {
                cnst.PEAK_AREA_KEY : area,
                cnst.PEAK_SIGMA_KEY : sigma,
                cnst.PEAK_CENTER_KEY : center,
                cnst.PEAK_FWHM_KEY : fwhm,
                cnst.PEAK_INTENSITY_KEY : peak_int,
                cnst.BACKGROUND_INT_KEY : bckgrnd_int,
                cnst.PEAK_TH_ENERGY_KEY : th_energy,
                cnst.PEAK_HEIGHT_KEY : height,
                cnst.PB_RATIO_KEY : PB_ratio
            }
    
        self.fitted_peaks_info = fitted_peaks_info
    
        # Determine which elements are present and their corresponding reference peak for quantification
        ref_lines_for_quant = []
        els_absent = []
        for el in self.fitted_els_quant:
            el_line_qnt = self._get_el_line_to_quantify(el)
            if el_line_qnt:
                ref_lines_for_quant.append(el_line_qnt)
            else:
                els_absent.append(el)
    
        self.ref_lines_for_quant = ref_lines_for_quant # List of el_lines used for quantification
        self.fitted_els_quant = [el for el in self.fitted_els_quant if el not in els_absent] # Updates list of quantified elements
        self.missing_reference_peaks = els_absent


    #%% Setup compositional quantification
    # =============================================================================
    def _get_el_line_to_quantify(self, el: str) -> Optional[str]:
        """
        Returns the characteristic X-ray line (e.g., 'Si_Ka', 'Mn_La') to use for quantification
        for the given element. Prefers Ka lines, then La, then Ma, based on presence and relative intensity.
    
        Parameters
        ----------
        el : str
            Element symbol (e.g., 'Si', 'Mn').
    
        Returns
        -------
        el_line_qnt : str or None
            Characteristic X-ray line to use for quantification, or None if not available.
    
        Warns
        -----
        UserWarning
            If no Ka, La, or Ma line is detected for the element.
    
        Notes
        -----
        - If multiple lines are present, prefers those above 2 keV and with overvoltage > 1.65.
          For 15 keV, this corresponds to Ka1 lines up to Z=31 (Ga), La1 lines up to Z=78 (Pt),
          and Ma1 lines for Au and heavier elements.
        - Among candidates, selects the line with the highest fitted peak area.
        - Rejects quantification when a higher-energy ideal reference line for the element is
          missing from the fit, indicating the fitted low-energy line is likely a detector artifact.
          The warning is suppressed when the missing line lies outside the fitted spectrum range
          (e.g. after cutting the spectrum below that line energy).
        """
        IDEAL_REF_LINE_OVERVOLTAGE = 1.65
        IDEAL_REF_LINE_ENERGY_THRESHOLD = 2
        # Get list of fitted lines for element el (only Ka, La, Ma lines considered)
        el_lines_list = [
            el_line for el_line in self.fitted_xray_lines
            if any(el + '_' + line == el_line for line in self.xray_quant_ref_lines)
        ]
        n_lines = len(el_lines_list)
    
        if n_lines == 0:
            warnings.warn(
                f'Element {el} was not quantified because it did not possess Ka, La, nor Ma lines'
            )
            el_line_qnt = None
        elif n_lines == 1:
            el_line_qnt = el_lines_list[0]
        else:
            # Get energies of lines
            lines_energies = [
                self.fitted_peaks_info[el_line][cnst.PEAK_TH_ENERGY_KEY] for el_line in el_lines_list
            ]
            # Define ideal reference lines as those above 2 keV and overvoltage > 1.65
            best_lines = [
                el_line for el_line, energy in zip(el_lines_list, lines_energies)
                if energy > IDEAL_REF_LINE_ENERGY_THRESHOLD and self.beam_energy / energy > IDEAL_REF_LINE_OVERVOLTAGE
            ]
    
            # Select ideal el_line for quantification
            if len(best_lines) == 1:
                el_line_qnt = best_lines[0]
            elif len(best_lines) > 1:
                # Select line with largest intensity among ideal lines
                el_line_qnt = max(best_lines, key=lambda el_line: self.fitted_peaks_info[el_line][cnst.PEAK_AREA_KEY])
            else:
                # Select line with largest intensity among all available lines
                el_line_qnt = max(el_lines_list, key=lambda el_line: self.fitted_peaks_info[el_line][cnst.PEAK_AREA_KEY])
    
        if el_line_qnt is not None:
            selected_energy = self.fitted_peaks_info[el_line_qnt][cnst.PEAK_TH_ENERGY_KEY]
            el_xray_lines = get_el_xray_lines(el)
            for ref_line in self.xray_quant_ref_lines:
                if ref_line not in el_xray_lines:
                    continue
                line_energy = el_xray_lines[ref_line]['energy (keV)']
                if (
                    line_energy > IDEAL_REF_LINE_ENERGY_THRESHOLD
                    and self.beam_energy / line_energy > IDEAL_REF_LINE_OVERVOLTAGE
                    and line_energy > selected_energy
                    and self.energy_vals[0] < line_energy < self.energy_vals[-1]
                    and f"{el}_{ref_line}" not in el_lines_list
                ):
                    warnings.warn(
                        f'Element {el} was not quantified because {el}_{ref_line} '
                        f'(ideal reference line at {line_energy:.2f} keV) is missing from the fit, '
                        f'while {el_line_qnt} was fitted — the element is likely absent '
                        f'(possible detector artifact).'
                    )
                    el_line_qnt = None
                    break
    
        return el_line_qnt
    
    
    def _load_EDS_standards(self) -> None:
        """
        Load EDS standards for the current beam energy and EDS mode.
    
        This method loads the standards from the corresponding microscope calibration folder
        and stores the relevant standards for the current EDS mode in `self.standards`.
    
        The resulting format is a dictionary mapping 'element_line' to standard P/B values.
    
        Raises
        ------
        KeyError
            If the current EDS mode (`self.meas_mode`) is not found in the loaded standards.
        """
        standards, _ = calibs.load_standards(self.meas_type, self.beam_energy)
        try:
            self.standards = standards[self.meas_mode]
        except KeyError as e:
            raise KeyError(
                f"EDS mode '{self.meas_mode}' not found in standards for beam energy {self.beam_energy}.\n"
                f"Available modes: {list(standards.keys())}"
            ) from e


    #%% Launch quantification
    # =============================================================================
    def quantify_spectrum(
        self,
        force_single_iteration=False,
        interrupt_fits_bad_spectra=True,
        print_result=True
    ):
        """
        Quantifies a spectrum, using iterative fitting to refine elemental fractions.
        At each iteration, elemental fractions are quantified and enforced in the next iteration.
        Algorithm converges when quantified elemental fractions all change by less than 0.01%.
        
        Instead, if elemental fractions are provided (i.e., fitting of experimental standards),
        background values are provided, or if force_single_iteration = True, only a single
        iteration is performed.
    
        Parameters
        ----------
        force_single_iteration : bool, optional
            If True, only a single fit iteration is performed, even if some elemental fractions are undefined.
        interrupt_fits_bad_spectra : bool, optional
            If True, fitting is interrupted early if the spectrum is deemed unreliable.
        print_result : bool, optional
            If True, prints the quantification results.
    
        Returns
        -------
        quant_result : dict or None
            Dictionary with quantification results, or None if fit failed.
        min_bckgrnd_ref_lines : float
            Minimum background counts across reference peaks, or 0 if fit failed.
        bad_quant_flag : int or None
            Flag indicating quantification issues, or None if successful.
        """
    
        is_fit_valid = True
        bad_quant_flag = None
        initial_weights_dict = {}
        iter_counter = 1 # Iteration counter
        max_iterations = 30  # Maximum allowed iterations
        converged = True

    
        # Get initial value for parameter 'K' before initializing the fitter (Iteration 0)
        initial_param_values = None
        if self.is_particle and self.fit_background:
            k_val = self.get_starting_K_val()
            if k_val is not None:
                initial_param_values = {'K': k_val}
    
        # Initialize the spectrum fitter (no fitting is performed yet)
        self._initialize_spectrum_fitter()
    
        # Check if all elemental fractions are defined
        has_undefined_fractions = not set(self.els_to_quantify).issubset(self.els_w_fr.keys())
    
        if force_single_iteration:
            fit_iteratively = False
            if has_undefined_fractions:
                missing = set(self.els_to_quantify) - set(self.els_w_fr.keys())
                warnings.warn(
                    "Not all elemental weight fractions are defined during fitting.\n"
                    "This may lead to fitting and quantification errors.\n"
                    f"Missing elemental fractions for: {', '.join(missing)}\n"
                    "If measuring standards, ensure all elemental fractions are specified."
                    "Alternatively, set force_single_iteration = False to iteratively fit unknown elemental fractions.",
                    UserWarning
                )
        else:
            # Determine if iterative fitting is needed
            fit_iteratively = has_undefined_fractions and self.fit_background
            if not has_undefined_fractions:
                warnings.warn(
                    "All elemental weight fractions are defined within 'els_w_fr'.\n"
                    "Spectrum will not be quantified iteratively.",
                    UserWarning
                )
    
        # Set fit tolerance
        if fit_iteratively:
            initial_fit_tolerance = 1e-2  # Quick fit: elemental fractions are likely far off during the first iteration, so fitting with high precision is unnecessary        else:
        else:
            initial_fit_tolerance = self.fit_tol # Single-iteration fitting
        
        # Perform initial fit (Iteration 1)
        try:
            fitted_params, weight_fractions, sample_Z = self._fit_quant_spectrum_iter(
                initial_par_vals=initial_param_values,
                f_tol=initial_fit_tolerance,
                n_iter=iter_counter
            )
        except Exception as e:
            is_fit_valid = False
            tb_str = traceback.format_exc()  # get full traceback as a string
            logger.error("❌ Fit and quantification iteration unsuccessful due to the following error:")
            logger.error(f"❌ {type(e).__name__}: {e}")
            logger.debug(tb_str)
    
        # Iteratively fit and quantify to converge to a solution
        if is_fit_valid and fit_iteratively:
            w_fr_change_convergence = 0.0001 # ZAF corrections converge when change in mass fraction is less than 0.01%
            diff_mass_fractions = 1  # To monitor convergence
            converged = False

            # Normalize mass fractions
            prev_weight_fractions = self._normalise_mass_fractions(weight_fractions)
    
            while iter_counter < max_iterations and diff_mass_fractions > w_fr_change_convergence:
                iter_counter += 1
                # Fix elemental fractions to values from previous iteration (normalized)
                for el, w_fr in zip(self.fitted_els_quant, prev_weight_fractions):
                    fitted_params['f_' + el].value = w_fr
                    fitted_params['f_' + el].vary = False
                    fitted_params['f_' + el].expr = None
    
                if iter_counter == 2:
                    # Fix sum fraction parameters not linked to the model anymore
                    sum_params = [p for p in fitted_params if 'sum_' in p]
                    for param in sum_params:
                        fitted_params[param].vary = False
    
                    # For particles, reset geometric factors to avoid local minima
                    if self.is_particle:
                        if initial_param_values is None:
                            initial_param_values = {}
                        initial_param_values['rhoz_par_slope'] = 0
                        initial_param_values['rhoz_par_offset'] = 0
                        initial_param_values['rhoz_lim'] = 0.001
                else:
                    initial_param_values = None
    
                # Adjust peak weights after the first re-iteration
                if iter_counter > 2:
                    fitted_params = self._update_peak_weights(
                        fitted_params, iter_counter, initial_weights_dict
                    )
    
                # Perform spectrum fit for this iteration
                fitted_params, weight_fractions, sample_Z = self._fit_quant_spectrum_iter(
                    params=fitted_params,
                    initial_par_vals=initial_param_values,
                    f_tol=self.fit_tol,
                    n_iter=iter_counter
                )
    
                # After 2nd re-iteration, check for unreliable quantification and possibly interrupt
                if iter_counter > 3 and interrupt_fits_bad_spectra:
                    analytical_error = sum(weight_fractions) - 1
                    bad_quant_flag = self._check_if_unreliable_quant(
                        iter_counter, analytical_error, interrupt_fits_bad_spectra
                    )
                    if bad_quant_flag is not None:
                        is_fit_valid = False
                        break
    
                # Check convergence of mass fractions
                norm_mass_fractions = self._normalise_mass_fractions(weight_fractions)
                diff_mass_fractions = np.max(np.abs(prev_weight_fractions - norm_mass_fractions))
    
                # Update for next iteration
                prev_weight_fractions = norm_mass_fractions

            converged = diff_mass_fractions <= w_fr_change_convergence
    
            if self.verbose:
                logger.info(f"✅ Spectrum fitted with {iter_counter} iterations")
    
        # Assemble and print quantification results
        if is_fit_valid:
            if self.verbose:
                self.fitter.print_result(extended=self.fitting_verbose)
    
            analytical_error = sum(weight_fractions) - 1
    
            # Update quantification flag if necessary
            bad_quant_flag = self._check_if_unreliable_quant(
                iter_counter, analytical_error, interrupt_fits_bad_spectra
            )
    
            # If not converged, set quantification flag
            if bad_quant_flag is None and iter_counter == max_iterations:
                bad_quant_flag = -1
    
            # Assemble results dictionary
            quant_result = self._assemble_quantification_result(weight_fractions, analytical_error)
    
            if print_result:
                if bad_quant_flag is None:
                    bad_quant_flag = 0
                self.print_quant_result(quant_result, sample_Z, bad_quant_flag)
    
            # Get minimum background counts for reference peaks
            min_bckgrnd_ref_lines = self._get_min_bckgrnd_cnts_ref_quant_lines()
    
        else:
            quant_result = None
            min_bckgrnd_ref_lines = 0
            converged = False
        
        self.bad_quant_flag = bad_quant_flag
        self.iterations_run = iter_counter
        self.quant_converged = converged and is_fit_valid
        
        return quant_result, min_bckgrnd_ref_lines, bad_quant_flag


    def _get_min_bckgrnd_cnts_ref_quant_lines(self):
        """
        Returns the minimum background counts measured around the reference peaks used for quantification.
    
        For each reference line, the minimum background value is found within ±1 FWHM of the peak center.
        This helps to detect excessive absorption around reference peaks, which can cause quantification errors.
        """
        min_bckgrnd_ref_lines = float('inf')  # Initialize to a very large value
    
        for el_line in self.ref_lines_for_quant:
            # Get peak center and FWHM for this reference line
            peak_center = self.fitted_peaks_info[el_line][cnst.PEAK_CENTER_KEY]
            peak_fwhm = self.fitted_peaks_info[el_line][cnst.PEAK_FWHM_KEY]
    
            # Find indices within ±1 FWHM of the peak center
            peak_indices = [
                i for i, energy in enumerate(self.energy_vals)
                if (peak_center - peak_fwhm) < energy < (peak_center + peak_fwhm)
            ]
    
            # Find the minimum background value in this region
            if peak_indices:
                min_background = min(self.background_vals[i] for i in peak_indices)
                # Update the minimum across all reference lines
                min_bckgrnd_ref_lines = min(min_bckgrnd_ref_lines, min_background)
    
        return min_bckgrnd_ref_lines


    def _get_min_bckgrnd_cnts_by_ref_quant_line(self) -> Dict[str, float]:
        """Return minimum background counts for each reference line used for quantification."""
        min_background_by_line: Dict[str, float] = {}

        if not hasattr(self, 'background_vals'):
            return min_background_by_line

        for el_line in self.ref_lines_for_quant:
            peak_center = self.fitted_peaks_info[el_line][cnst.PEAK_CENTER_KEY]
            peak_fwhm = self.fitted_peaks_info[el_line][cnst.PEAK_FWHM_KEY]

            peak_indices = [
                i for i, energy in enumerate(self.energy_vals)
                if (peak_center - peak_fwhm) < energy < (peak_center + peak_fwhm)
            ]

            if peak_indices:
                min_background_by_line[el_line] = float(min(self.background_vals[i] for i in peak_indices))

        return min_background_by_line


    def _get_r_squared(self) -> float:
        """Compute R-squared from the current fit."""
        return float(1 - self.fit_result.residual.var() / np.var(self.spectrum_vals))


    def export_fit_result(self) -> QuantFitResult:
        """Build a schema-ready fit summary from the current fitted peak state."""
        if not hasattr(self, 'fitted_peaks_info') or not hasattr(self, 'fit_result'):
            raise ValueError("Fit results are not available. Fit the spectrum before exporting fit data.")

        fitted_peaks: Dict[str, FittedPeakResult] = {}
        for el_line, peak_info in self.fitted_peaks_info.items():
            try:
                element, line = el_line.split('_', 1)
            except ValueError as exc:
                raise ValueError(f"Unexpected fitted peak key format: '{el_line}'") from exc

            fitted_peaks[el_line] = FittedPeakResult(
                element=element,
                line=line,
                area=peak_info.get(cnst.PEAK_AREA_KEY),
                sigma=peak_info.get(cnst.PEAK_SIGMA_KEY),
                center=peak_info.get(cnst.PEAK_CENTER_KEY),
                fwhm=peak_info.get(cnst.PEAK_FWHM_KEY),
                peak_intensity=peak_info.get(cnst.PEAK_INTENSITY_KEY),
                background_intensity=peak_info.get(cnst.BACKGROUND_INT_KEY),
                theoretical_energy=peak_info.get(cnst.PEAK_TH_ENERGY_KEY),
                height=peak_info.get(cnst.PEAK_HEIGHT_KEY),
                pb_ratio=peak_info.get(cnst.PB_RATIO_KEY),
            )

        return QuantFitResult(
            r_squared=self._get_r_squared(),
            reduced_chi_squared=float(self.fit_result.redchi),
            fitted_peaks=fitted_peaks,
        )


    def export_quantification_result(
        self,
        quantification_id: int,
        quant_result: Optional[Dict[str, Any]],
        quant_flag: Optional[int] = None,
        comment: Optional[str] = None,
    ) -> QuantificationResult:
        """Build a persisted quantification result for the current fitted spectrum."""
        min_background_by_line = self._get_min_bckgrnd_cnts_by_ref_quant_line()
        min_background_ref_lines = None
        if min_background_by_line:
            min_background_ref_lines = float(min(min_background_by_line.values()))

        return QuantificationResult(
            quantification_id=quantification_id,
            quant_flag=quant_flag,
            comment=comment,
            composition_atomic_fractions=(
                dict(quant_result[cnst.COMP_AT_FR_KEY]) if quant_result else None
            ),
            composition_weight_fractions=(
                dict(quant_result[cnst.COMP_W_FR_KEY]) if quant_result else None
            ),
            analytical_error=(
                float(quant_result[cnst.AN_ER_KEY]) if quant_result and cnst.AN_ER_KEY in quant_result else None
            ),
            fit_result=self.export_fit_result() if hasattr(self, 'fit_result') else None,
            diagnostics=QuantificationDiagnostics(
                iterations_run=self.iterations_run,
                converged=self.quant_converged,
                interrupted=(quant_result is None),
                min_background_ref_lines=min_background_ref_lines,
                missing_reference_peaks=list(self.missing_reference_peaks),
            ),
        )
    
    
    def _assemble_quantification_result(self, weight_fractions, analytical_er):
        """
        Assemble the quantification results into a dictionary.
    
        Parameters
        ----------
        weight_fractions : list of float
            List of elemental weight fractions.
        analytical_er : float
            Computed analytical error.
    
        Returns
        -------
        quant_result : dict
            Dictionary containing atomic and weight fractions, analytical error, reduced chi-square, and R-squared metrics.
        """
        # Convert weight fractions to atomic fractions
        atomic_fractions = weight_to_atomic_fr(weight_fractions, self.fitted_els_quant, verbose=False)
    
        # Initialize result dictionary with keys for atomic and weight fractions
        quant_result = {
            cnst.COMP_AT_FR_KEY: {},
            cnst.COMP_W_FR_KEY: {}
        }
    
        # Fill in the rounded atomic and weight fractions for each element
        for el, at_fr, w_fr in zip(self.fitted_els_quant, atomic_fractions, weight_fractions):
            quant_result[cnst.COMP_AT_FR_KEY][el] = round(at_fr, 4)
            quant_result[cnst.COMP_W_FR_KEY][el] = round(w_fr, 4)
    
        # Add analytical error, reduced chi-square, and R-squared to results
        quant_result[cnst.AN_ER_KEY] = round(analytical_er, 4)
        quant_result[cnst.REDCHI_SQ_KEY] = round(self.fit_result.redchi, 1)
        quant_result[cnst.R_SQ_KEY] = round(self._get_r_squared(), 6)
    
        return quant_result
    
    
    def _fit_quant_spectrum_iter(
        self,
        params: Optional[lmfit.Parameters] = None,
        initial_par_vals: Optional[Dict[str, float]] = None,
        f_tol: float = 1e-4,
        n_iter: Optional[int] = None
    ) -> Tuple[lmfit.Parameters, np.ndarray, float]:
        """
        Perform a single spectrum fit and quantification iteration for iterative quantification workflows.
    
        This method assumes that the spectrum fitter (`self.fitter`) is already set up. It performs
        a fit (without printing results), copies the fitted parameters, and calculates mass fractions
        and mean atomic number.
    
        Parameters
        ----------
        params : Optional[lmfit.Parameters], optional
            Parameters object to pass to the spectrum fitter. If None, default parameters are used.
        initial_par_vals : Optional[Dict[str, float]], optional
            Initial parameter values for the fit. If None, default initial values are used.
        f_tol : float, optional
            Function tolerance for the fitting algorithm (default is 1e-4).
        n_iter : Optional[int], optional
            Maximum number of fitting iterations. If None, the default is used.
    
        Returns
        -------
        fitted_params : lmfit.Parameters
            Copy of the fitted parameters object from the spectrum fit.
        weight_fractions : np.ndarray
            Quantified mass fractions for each element.
        sample_Z : float
            Mean atomic number for the quantified sample.
    
        Notes
        -----
        - This method is intended for iterative quantification workflows, where the spectrum fitter
          (`self.fitter`) has already been initialized.
        - Calls `self._fit_spectrum()` and `self._quantify_mass_fractions()` internally.
        """
        self._fit_spectrum(
            params=params,
            initial_par_vals=initial_par_vals,
            f_tol=f_tol,
            n_iter=n_iter,
            print_result=False
        )
        fitted_params = self.fit_result.params.copy()
        weight_fractions, sample_Z = self._quantify_mass_fractions()
    
        return fitted_params, weight_fractions, sample_Z
    
    
    #%% Quantification algorithm
    # =============================================================================
    def _initialize_k_ratios(self, k_ratios: np.ndarray) -> np.ndarray:
        """
        Normalizes k-ratios so their sum is 1.
    
        Parameters
        ----------
        k_ratios : np.ndarray
            Array of initial k-ratio values.
    
        Returns
        -------
        norm_conc : np.ndarray
            Normalized k-ratios (sum to 1).
    
        Notes
        -----
        A potential future addition is to assign any missing concentration (if the sum is less than 1,
        due to unquantified elements) to an 'unquantified' element, thereby reflecting the analytical total error.
        """
        tot_conc = np.sum(k_ratios)
        norm_conc = k_ratios / tot_conc
    
        return norm_conc
    
    
    def _get_k_ratios(self):
        """
        Calculates k-ratios for quantification lines using only the average measured standard P/B for each element.
    
        Returns
        -------
        k_ratios : list of float
            List of k-ratios (one per reference quantification line), using only the 'Mean' standard.
        
        Potential improvements
        -----
        Placeholder sections are included for possible corrections such as substrate signal contamination 
        correction. For example, corrections from Essani et al.:
            M. Essani, E. Brackx, and E. Excoffier, A method for the correction of size effects in
            microparticles using a peak-to-background approach in electron-probe microanalysis,
            Spectrochim. Acta - Part BAt. Spectrosc. 169, 105880 (2020).
        """
        if self.standards is None:
            self._load_EDS_standards()
    
        k_ratios = []

        # --- Placeholder: Correction for substrate signal contamination ---
        # # Calculate correction in function of substrate peak intensity 
        # self._calc_sub_bckgrnd_correction()
        # ---------------------------------------------------------------
    
        for el_line in self.ref_lines_for_quant:
            # Retrieve standard measurements for this reference
            if el_line in self.standards:
                std_vals_list = self.standards[el_line]
            else:
                raise EDSError(
                    f"The {el_line} characteristic X-ray is not present in the standards database "
                    f"of the '{self.meas_mode}' EDS mode"
                )
    
            # Get measured PB ratio (=0 if element was not found in spectrum)
            if el_line in self.fitted_peaks_info:
                meas_PB_ratio = self.fitted_peaks_info[el_line][cnst.PB_RATIO_KEY]
            else:
                meas_PB_ratio = 0
                    
            # --- Placeholder: Correction for substrate signal contamination ---
            # if self.is_particle and self.is_substrate_peak_present and meas_PB_ratio > 0:
            #     meas_PB_ratio = self._correct_PB_for_sub_bckgrnd(meas_PB_ratio, el_line)
            # ---------------------------------------------------------------
            
            # Only use the standard where Std == 'Mean'
            mean_std = next((std for std in std_vals_list if std.get(cnst.STD_ID_KEY) == cnst.STD_MEAN_ID_KEY and cnst.COR_PB_DF_KEY in std and std[cnst.COR_PB_DF_KEY] is not None), None)
            if mean_std is not None:
                std_PB_ratio = mean_std[cnst.COR_PB_DF_KEY]
                if std_PB_ratio <= 0:
                    raise EDSError(
                        f"'{cnst.STD_MEAN_ID_KEY}' PB ratio for {el_line} standard is not >0, unphysical."
                    )
                else:
                    k_ratio_val = meas_PB_ratio / std_PB_ratio
                k_ratios.append(k_ratio_val)
            else:
                raise EDSError(
                    f"No valid standard with '{cnst.STD_ID_KEY}' == '{cnst.STD_MEAN_ID_KEY}' found for {el_line} in the standards database."
                )
        return k_ratios


    def _normalise_mass_fractions(self, weight_fractions):
        """
        Normalizes the list of elemental weight fractions according to total mass fraction constraints.
    
        If total mass fraction exceeds 1, or if forced by settings, the fractions are normalized to sum to 1.
        If total mass fraction is less than the minimum allowed (to account for undetectable elements),
        the fractions are scaled to sum to the minimum allowed value.
        Otherwise, the fractions are returned unchanged.
    
        Parameters
        ----------
        weight_fractions : list or array-like of float
            Elemental weight fractions to be normalized.
    
        Returns
        -------
        w_frs : numpy.ndarray
            Normalized weight fractions.
        """
        min_total = 1 - self.max_undetectable_w_fr
        max_total = 1
    
        total = sum(weight_fractions)
    
        if self.force_total_w_fr or total > max_total:
            # Normalize so that the sum of fractions is exactly 1
            w_frs = np.array(weight_fractions) / total
        elif total < min_total:
            # Scale so that the sum matches the minimum allowed total mass fraction
            w_frs = np.array(weight_fractions) / total * min_total
        else:
            # No normalization needed
            w_frs = np.array(weight_fractions)
    
        return w_frs


    def _quantify_mass_fractions(self):
        """
        Calculates and returns the ZAF-corrected elemental mass fractions for the sample using the PB method.
    
        This function uses the measured PB k-ratios (i.e., PB_sample / PB_standard), applying ZAF corrections.
    
        Returns
        -------
        weight_fractions : numpy.ndarray
            Array of ZAF-corrected elemental mass fractions.
        sample_Z : dict
            Dictionary containing sample mean atomic numbers under different averaging conventions.
        """
        # Get theoretical energies for each X-ray peak used in quantification
        el_peak_energies = np.array([
            self.fitted_peaks_info[el_line][cnst.PEAK_TH_ENERGY_KEY]
            for el_line in self.ref_lines_for_quant
        ])
    
        # Calculate k-ratios
        k_ratios = self._get_k_ratios()
    
        # Get ZAF-corrected mass fractions and sample mean atomic numbers
        weight_fractions, Z_sample = self._correct_ZAF(
            k_ratios, el_peak_energies=el_peak_energies
        )
    
        return weight_fractions, Z_sample


    def _correct_ZAF(self, k_ratios, el_peak_energies):
        """
        Applies iterative ZAF corrections to k-ratios to determine converged elemental mass fractions.
    
        The iteration procedure is based on K.F.J. Heinrich (1972), but uses ZAF corrections as in
        J.L. Lábár & S. Török (1992), without parabolic approximation.
    
        Parameters
        ----------
        k_ratios : array-like
            Measured k-ratios for each quantified element/line.
        el_peak_energies : array-like
            Theoretical energies of the X-ray peaks used for quantification.
    
        Returns
        -------
        weight_fractions : numpy.ndarray
            Converged ZAF-corrected elemental mass fractions.
        sample_Z : dict
            Dictionary containing sample mean atomic numbers under different averaging conventions.
            
        References
        ----------
        K. F. J. Heinrich, A Simple Correction Procedure for Quantitative Electron Probe Microanalysis, 1972.
        J. L. Lábár and S. Török, A peak‐to‐background method for electron‐probe x‐ray microanalysis applied to
            individual small particles, X-Ray Spectrom. 21, 183 (1992).
        """
        # Initialize quantification correction class
        correction = Quant_Corrections(
            self.fitted_els_quant,
            self.beam_energy,
            self.emergence_angle,
            self.meas_mode,
            el_peak_energies,
            verbose=self.verbose
        )

    
        # Iterative ZAF correction parameters
        ZAF_cntr = 0 # Iteration counter
        max_iter = 20 # Max number of iterations
        
        # Initialize convergence parameters
        max_diff = 0.5      # Convergence counter: start with a large value to ensure at least one iteration
        converge_tol = 1e-4 # Convergence condition: stop when max difference in elemental fractions is below 0.01%
    
        # Start with initial guess for weight fractions (usually just k-ratios)
        weight_fractions = self._initialize_k_ratios(k_ratios)
        self.are_w_fr_norm = False
    
        if self.verbose:
            print_double_separator()
            logger.info('🔬 Quantification with ZAF correction:')
            print_single_separator()
            print_nice_1d_row('', self.fitted_els_quant)
            print_nice_1d_row('Initial W_fr', k_ratios)
            logger.info(f"ℹ️  Initial analytical error: {(sum(k_ratios) - 1) * 100:.2f}%")
    
        while max_diff > converge_tol and ZAF_cntr < max_iter:
            ZAF_cntr += 1
            if self.verbose:
                print_single_separator()
                logger.debug(f"  ▶️ Step: {ZAF_cntr}")
    
            # Calculate ZAF factors and sample mean Z
            ZAF_pb_factors, sample_Z = correction.get_ZAF_mult_f_pb(weight_fractions)
    
            # Update weight fractions
            new_weight_fractions = k_ratios * ZAF_pb_factors
            if self.verbose:
                print_nice_1d_row('New W_fr', new_weight_fractions)
                logger.info(f"ℹ️  New analytical error: {(sum(new_weight_fractions) - 1) * 100:.2f}%")
    
            max_diff = max(abs(new_weight_fractions - weight_fractions))
            weight_fractions = new_weight_fractions.copy()
    
        if ZAF_cntr == max_iter:
            print_single_separator()
            logger.warning(f'⚠️ ZAF correction did not converge within {max_iter} iterations.')
        elif self.verbose:
            print_single_separator()
            logger.info(f"✅ ZAF correction converged in {ZAF_cntr} steps.")
    
        return weight_fractions, sample_Z


    def _update_peak_weights(self, fitted_params, iter_cntr, initial_weights_dict):
        """
        Adjusts the weights of dependent X-ray peaks in the spectrum fit parameters
        to account for absorption differences between the measured material and pure standards.
    
        This method updates the area expressions for dependent peaks (such as satellite lines), 
        so that their areas are correctly scaled according to the absorption profile of the sample.
        It uses the absorption correction factors evaluated for both the dependent and reference peaks.
        
        This method is not used when the background values are provided by the user, since the
        absorption attenuation profile is unknown.
    
        Parameters
        ----------
        fitted_params : dict
            Dictionary of fit parameters, as used by the fitting engine.
        iter_cntr : int
            Current iteration number in the quantification procedure.
        initial_weights_dict : dict
            Dictionary to cache the initial weights of dependent peaks as calculated for pure standards.
    
        Returns
        -------
        fitted_params : dict
            The updated dictionary of fit parameters with adjusted peak weights.
        """
    
        area_key = Peaks_Model.area_key
        area_weight_pattern = rf"(.*){area_key}\*(\d+\.?\d*)"  # Pattern for ref_line_prefix + area_key * weight
        area_param_name_pattern = rf"(.*){area_key}$"          # Pattern for line_prefix + area_key
    
        # TODO: Add energy shifts for pileup and escape peaks. Right now these are ignored since mostly negligible.
        pile_up_str = self.fitter.pileup_peaks_str
        escape_up_str = self.fitter.escape_peaks_str
        fixed_peaks = [
            'Ti_Ln', 'Ti_Ll', 'Ti_Lb1', 'Fe_Lb1', 'Co_Lb1',
            'Zn_Lb1', 'Cu_Ll', 'Cu_Ln'
        ]  # TODO: Should calibrate these peaks' areas and remove them from this list.
    
        # Create a copy of the parameters to modify for calculations.
        fitted_params_for_calcs = fitted_params.copy()
        fitted_params_for_calcs['apply_det_response'].value = 0  # Remove detector response for absorption calculation.
        Background_Model._clear_cached_abs_att_variables() # Clear cached absorption attenuation values
        
        for param in fitted_params:
            # Adjust area parameter for all dependent peaks.
            if (
                area_key in param
                and fitted_params[param].expr is not None
                and not any(s in param for s in [pile_up_str, escape_up_str, *fixed_peaks])
            ):
                match_param_name = re.match(area_param_name_pattern, fitted_params[param].name)
                match_weight = re.match(area_weight_pattern, fitted_params[param].expr)
                if match_param_name and match_weight:
                    # Get absorption value at the dependent peak
                    line_prefix = match_param_name.group(1)
                    el, line = line_prefix.split('_')[:2]
                    line_en = np.array([get_el_xray_lines(el)[line]['energy (keV)']])
                    line_abs_val = Background_Model._abs_attenuation_phirho(line_en, **fitted_params_for_calcs)
    
                    # Get absorption value at the reference peak
                    ref_line_prefix = match_weight.group(1)
                    ref_line = ref_line_prefix.split('_')[1]
                    ref_line_en = np.array([get_el_xray_lines(el)[ref_line]['energy (keV)']])
                    ref_line_abs_val = Background_Model._abs_attenuation_phirho(ref_line_en, **fitted_params_for_calcs)
    
                    # Get or cache the weight in pure material
                    if iter_cntr == 3:
                        el_fr_param = {f'f_{el}': 1}
                        Background_Model(
                            self.is_particle,
                            beam_energy=self.beam_energy,
                            emergence_angle=self.emergence_angle,
                        )  # Re-initialize absorption attenuation globals
                        line_abs_val_pure = Background_Model._abs_attenuation_phirho(
                            line_en, det_angle=self.emergence_angle, adr_abs=0, **el_fr_param
                        )
                        ref_line_abs_val_pure = Background_Model._abs_attenuation_phirho(
                            ref_line_en, det_angle=self.emergence_angle, adr_abs=0, **el_fr_param
                        )
                        weight_NIST = get_el_xray_lines(el)[line]['weight']  # NIST weight for dependent peak
                        weight_wo_absorption = weight_NIST * ref_line_abs_val_pure / line_abs_val_pure
                        initial_weights_dict[line_prefix] = weight_wo_absorption
                    else:
                        weight_wo_absorption = initial_weights_dict[line_prefix]
    
                    # Adjust weight based on fitted absorption profile
                    updated_weight = (weight_wo_absorption * line_abs_val / ref_line_abs_val)[0]
                    fitted_params[param].expr = ref_line_prefix + area_key + f"*{updated_weight:.6f}"
                else:
                    continue
    
        return fitted_params


    def _check_if_unreliable_quant(self, iter_cntr, analytical_er, interrupt_fits_bad_spectra):
        """
        Checks for conditions indicating unreliable quantification, such as poor fit quality,
        high analytical error, or excessive absorption effects. Returns a flag if quantification should be halted.
    
        Parameters
        ----------
        iter_cntr : int
            Current iteration number in the quantification process.
        analytical_er : float
            Analytical error (fractional).
        interrupt_fits_bad_spectra : bool
            If True, prints a message when halting quantification due to detected spectral issues.
    
        Returns
        -------
        bad_quant_flag : int or None
            Returns:
                1 if reduced chi-squared is too high,
                2 if analytical error is too high,
                3 if absorption increase is excessive,
                None if quantification is considered reliable.
        """
        # Thresholds for unreliable quantification
        redchi_threshold = 0.16    # Threshold of reduced-chi squared value as % of total counts
        redchi_threshold_val = self.tot_sp_counts * redchi_threshold / 100
        an_err_threshold = 0.5    # Analytical error threshold (50 w%)
        abs_increase_threshold = 0.7  # Absorption increase threshold (170% of bulk absorption)
    
        bad_quant_flag = None  # Default: quantification is considered reliable
        abs_att_param_name = '_abs_attenuation_phirho' # TODO does not work for models different than 'phirho' absorption attenuation model
    
        # 1. Check for extremely poor fit (high reduced chi-squared)
        if self.fit_result.redchi > redchi_threshold_val:
            bad_quant_flag = 1
            if interrupt_fits_bad_spectra:
                logger.warning(f"⚠️ Quantification stopped at iteration #{iter_cntr} due to reduced chi-squared being "
                      f"{self.fit_result.redchi:.1f} > {redchi_threshold_val:.1f}")
    
        # 2. Check for excessive analytical error
        elif analytical_er > an_err_threshold:
            bad_quant_flag = 2
            if interrupt_fits_bad_spectra:
                logger.warning(f"⚠️ Quantification stopped at iteration #{iter_cntr} due to analytical error being "
                      f"{analytical_er*100:.1f}% > {an_err_threshold*100:.1f}%")
    
        # 3. For particles, check for excessive absorption around reference peaks
        # TODO does not work for models different than 'phirho' absorption attenuation model
        elif (
            self.is_particle and
            abs_att_param_name in [comp.func.__name__ for comp in self.fit_result.model.components]
        ):
            # Find lowest-energy reference peak
            en_val_ref_peak = min(
                self.fitted_peaks_info[ref_peak][cnst.PEAK_TH_ENERGY_KEY]
                for ref_peak in self.ref_lines_for_quant
            )
            # Use a 1 keV energy window around this peak
            ref_energy_vals = self.fitter.energy_vals[
                (self.fitter.energy_vals >= en_val_ref_peak - 0.5) &
                (self.fitter.energy_vals <= max(en_val_ref_peak + 0.5, 1))
            ]
    
            # Evaluate fitted absorption envelope (without detector response)
            params_wo_det_response = self.fit_result.params.copy()
            params_wo_det_response['apply_det_response'].value = 0
            fit_components_wo_det_response = self.fit_result.eval_components(
                params=params_wo_det_response, x=ref_energy_vals
            )
            fitted_abs_val = sum(fit_components_wo_det_response[abs_att_param_name])
    
            # Evaluate bulk absorption envelope (reset particle parameters)
            params_wo_det_response['rhoz_par_slope'].value = 0
            params_wo_det_response['rhoz_par_offset'].value = 0
            params_wo_det_response['rhoz_lim'].value = 0.001
            fit_components_wo_det_response = self.fit_result.eval_components(
                params=params_wo_det_response, x=ref_energy_vals
            )
            bulk_abs_val = sum(fit_components_wo_det_response[abs_att_param_name])
    
            # Compute increase in absorption and check against threshold
            abs_increase = 1 - fitted_abs_val / bulk_abs_val
            if abs_increase > abs_increase_threshold:
                bad_quant_flag = 3
                if interrupt_fits_bad_spectra:
                    logger.warning(f"⚠️ Quantification stopped at iteration #{iter_cntr} due to absorption around reference peaks being "
                          f"{abs_increase*100:.1f}% > {abs_increase_threshold*100:.1f}%")
    
        return bad_quant_flag
    
    #%% Post-quantification functions
    # =============================================================================
    def print_quant_result(
        self,
        quant_result: dict,
        Z_sample: dict = None,
        quant_flag: int | None = None
    ) -> None:
        """
        Print the quantification results, including fit quality metrics and elemental composition.
    
        Parameters
        ----------
        quant_result : dict
            Dictionary containing quantification results with keys:
                - cnst.REDCHI_SQ_KEY: Reduced Chi-square of the fit.
                - cnst.R_SQ_KEY: R-squared of the fit.
                - cnst.COMP_W_FR_KEY: Weight fractions (as decimals).
                - cnst.COMP_AT_FR_KEY: Atomic fractions (as decimals).
                - cnst.AN_ER_KEY: Analytical error (as decimal).
        Z_sample : dict, optional
            Dictionary containing mean atomic numbers (optional), e.g.:
                - 'Statham2016': Mean atomic number (Z̅) calculated according to Statham (2016).
                - 'mass-averaged': Mean atomic number (Z̅) weighted by composition.
        quant_flag : int | None, optional
            Flag signaling reliability of quantifification. Anything different than 0 signals
            potential unreliability of the quantified composition. See description of quant_flags]
            in config/classes.py
    
        Returns
        -------
        None
    
        References
        ----------
        Statham P, Penman C, Duncumb P. IOP Conf Ser Mater Sci Eng. 2016;109(1):0–10.
        """
        # Print a double separator for visual clarity
        print_double_separator()
        print_double_separator()
        
        logger.info('📊 Fit result:')
        
        logger.info(f"  Reduced Chi-square: {quant_result[cnst.REDCHI_SQ_KEY]:.2f}")
        logger.info(f"  R-squared: {quant_result[cnst.R_SQ_KEY]:.5f}")
        
        print('')
        
        logger.info('📊 Quantification result:')
        
        # Print list of fitted elements
        print_nice_1d_row('', self.fitted_els_quant)
        
        # Print atomic fractions as percentages
        at_fr_percent = [v * 100 for v in quant_result[cnst.COMP_AT_FR_KEY].values()]
        print_nice_1d_row('At_fr (%)', at_fr_percent)
        
        # Print weight fractions as percentages
        w_fr_percent = [v * 100 for v in quant_result[cnst.COMP_W_FR_KEY].values()]
        print_nice_1d_row('W_fr (%)', w_fr_percent)
        
        print('')
    
        # Print analytical error as a percentage (w%)
        an_err_percent = quant_result[cnst.AN_ER_KEY] * 100
        logger.info(f"  ℹ️  Analytical error: {an_err_percent:.2f} w%")
        
        if quant_flag is not None:
            logger.info(f"  ℹ️  Quantification flag: {quant_flag}")
    
        # Print mean atomic numbers (Z̅) only if provided
        if Z_sample is not None:
            print('')
            if 'Statham2016' in Z_sample:
                logger.info(f"  ℹ️  Z̅_Statham2016: {Z_sample['Statham2016']:.2f}")
            if 'mass-averaged' in Z_sample:
                logger.info(f"  ℹ️  Z̅_w (mass-averaged): {Z_sample['mass-averaged']:.2f}")
            
    
    def plot_quantified_spectrum(
        self,
        annotate_peaks: str = 'all',
        plot_bckgrnd_cnts_ref_peaks: bool = True,
        plot_initial_guess: bool = False,
        plot_title: Optional[str] = None,
        peaks_to_zoom: Optional[Union[str, List[str]]] = None
    ) -> None:
        """
        Plot the quantified spectrum.
        The background counts under reference peaks are highlighted for spectra used for quantification.
    
        Parameters
        ----------
        annotate_peaks : str, optional
            Which peaks to annotate. Options: 'all', 'most', 'main', 'none'. Default is 'all'.
        plot_bckgrnd_cnts_ref_peaks : bool, optional
            If True, plot vertical lines that illustrate value of background counts used for quantification.
            This is shown underneath the corresponding reference characteristic peaks for each element.
        plot_initial_guess : bool, optional
            If True, plot the initial guess as well. Default is False.
        plot_title : str, optional
            Title printed at the top of the plot. Default is None.
        peaks_to_zoom : str or list of str, optional
            Peak label (e.g. 'Si_Ka1') or list of labels to zoom in on. If provided, creates a new figure for each.
    
        Returns
        -------
        None.
        """
        # Accept a single string or a list for peaks_to_zoom
        if isinstance(peaks_to_zoom, str):
            if peaks_to_zoom != '':
                peaks_to_zoom = [peaks_to_zoom]
            else:
                peaks_to_zoom = []
        elif peaks_to_zoom is None:
            peaks_to_zoom = []
    
        # Plot data points + fit adding Phenom background
        if not self.fit_background:
            plt.figure()
            plt.plot(self.energy_vals, self.background_vals + self.spectrum_vals, 'o', label='data')
            fitted_points = self.fit_result.eval()
            plt.plot(self.energy_vals, self.background_vals + fitted_points, color='C1', label='spectrum fit')
        else:
            Background_Model(
                is_particle=self.is_particle,
                beam_energy=self.beam_energy,
                emergence_angle=self.emergence_angle
            )  # Re-initialise Background_Model variables
            fig = self.fit_result.plot()
            axes = fig.get_axes()
            axes[0].set_title('Residual plot')
    
        plt.grid(False)
    
        # Set font size
        fontsize = 12
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.titlesize'] = fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize
    
        # Plot background data
        plt.plot(self.energy_vals, self.background_vals, 'r--', linewidth=2, label='bckgrnd fit')
    
        # Annotate plot with Ka and La peak names
        params = self.fit_result.params
        main_peaks = [param for param in params.keys() if any(line + '_center' in param for line in self.xray_quant_ref_lines)]
        all_peaks = [param for param in params.keys() if '_center' in param and any(el in param for el in (self.els_to_quantify + self.els_substrate))]
        most_peaks = [param for param in all_peaks if any([params[param[:-7] + '_area'] >= 1, any(el in param for el in self.els_to_quantify)])]
        if annotate_peaks == 'most':
            peaks_to_plot = most_peaks
        elif annotate_peaks == 'all':
            peaks_to_plot = all_peaks
        elif annotate_peaks == 'main':
            peaks_to_plot = main_peaks
        elif annotate_peaks is None or str(annotate_peaks).lower() == 'none':
            peaks_to_plot = []
        else:
            logger.warning(f"⚠️ Warning: Unrecognized value for annotate_peaks ('{annotate_peaks}'). No peaks will be annotated.")
            peaks_to_plot = []
    
        y_limit = plt.gca().get_ylim()[1]
        for param in peaks_to_plot:
            el_line = param.split('_center')[0]
            center = self.fitted_peaks_info[el_line][cnst.PEAK_CENTER_KEY]
            interval_indices = np.where((self.energy_vals > center - 0.015) & (self.energy_vals < center + 0.015))[0]
            try:
                max_index = interval_indices[np.argmax(self.spectrum_vals[interval_indices])]
            except Exception:
                pass
            else:
                height = self.spectrum_vals[max_index]
                if not self.fit_background:
                    height += self.background_vals[max_index]
                pos_y = height + y_limit / 100
                plt.text(params[param], pos_y, '— ' + el_line, rotation=90, verticalalignment='bottom', horizontalalignment='center')
    
        # Highlight the background counts with vertical lines and a vertical legend handle
        if plot_bckgrnd_cnts_ref_peaks:
            first_line = True
            for el_line in self.ref_lines_for_quant:
                if first_line:
                    first_line = False
                    bckgrnd_cnts_label = 'background counts'
                else:
                    bckgrnd_cnts_label = ''
                    
                peak_center = self.fitted_peaks_info[el_line][cnst.PEAK_TH_ENERGY_KEY]
                if self.fit_background:
                    peak_bck_val = np.interp(peak_center, self.energy_vals_finer, self.background_vals_wo_det_response)
                else:
                    peak_bck_val = np.interp(peak_center, self.energy_vals, self.background_vals)
                plt.vlines(peak_center, ymin=0, ymax=peak_bck_val, color='red', alpha=1, label = bckgrnd_cnts_label)
    
        # Add initial guess
        if plot_initial_guess:
            init_params = self.fit_result.init_params
            plt.plot(self.energy_vals, self.fitter.spectrum_mod.eval(init_params, x=self.energy_vals), label='initial guess', color='black', linestyle=':')
            if self.fit_background:
                plt.plot(self.energy_vals, self.fitter.background_mod.eval(init_params, x=self.energy_vals), color='black', linestyle=':')
    
        plt.xlabel('Energy (keV)')
        plt.ylabel('Counts')
        if plot_title:
            plt.title(plot_title)
        plt.legend()
        plt.show()
    
        # ---- ZOOMED-IN PLOTS: create a new figure for each requested peak ----
        for peak in peaks_to_zoom:
            if peak not in self.fitted_peaks_info:
                logger.warning(f'⚠️ You have attempted to zoom on a peak, using the line {peak}.')
                logger.warning('⚠️ This line is absent from the list of fitted peaks, so the plot was not zoomed.')
                logger.info(f'ℹ️ The available peak lines are {self.fitted_xray_lines}')
            else:
                self.plot_zoomed_peak(peak, plot_title=plot_title)
    
    def plot_zoomed_peak(
        self,
        zoom_peak: str,
        plot_title: Optional[str] = None,
    ) -> None:
        """
        Create a new figure zoomed in on a specific peak.
    
        Parameters
        ----------
        zoom_peak : str
            The el_line string of the peak to zoom on (e.g. 'Si_Ka1').
        plot_title : str, optional
            Title for the zoomed plot.
    
        Returns
        -------
        None.
        """
        if zoom_peak not in self.fitted_peaks_info:
            logger.warning(f"⚠️ Peak '{zoom_peak}' not found in fitted peaks.")
            return
    
        fig_zoom, ax_zoom = plt.subplots()
        ax_zoom.plot(self.energy_vals, self.spectrum_vals, 'o', label='Data points')
        fitted_points = self.fit_result.eval()
        ax_zoom.plot(self.energy_vals, fitted_points, color='C1', label='Fitted model')
        ax_zoom.plot(self.energy_vals, self.background_vals, 'r--', linewidth=2, label='Background')
    
        el_line = zoom_peak
        peak_fwhm = self.fitted_peaks_info[zoom_peak][cnst.PEAK_FWHM_KEY]
        peak_center = self.fitted_peaks_info[el_line][cnst.PEAK_CENTER_KEY]
        peak_PB_ratio = self.fitted_peaks_info[el_line][cnst.PB_RATIO_KEY]
        xl_lim = peak_center - 1.5 * peak_fwhm
        xr_lim = peak_center + 1.5 * peak_fwhm
    
        # Find max data point within the zoom x-range
        in_range = (self.energy_vals >= xl_lim) & (self.energy_vals <= xr_lim)
        if np.any(in_range):
            max_point = np.max(self.spectrum_vals[in_range])
        else:
            max_point = np.max(self.spectrum_vals)  # fallback if no points in range
    
        ax_zoom.set_xlim(xl_lim, xr_lim)
        ax_zoom.set_ylim(0, max_point * 1.1)
        ax_zoom.text(peak_center, 0 + max_point * 0.1, "%s P/B: %.1f" % (el_line, peak_PB_ratio), fontsize=12,
                     color='black', horizontalalignment='center', verticalalignment='center')
    
        ax_zoom.set_xlabel('Energy (keV)')
        ax_zoom.set_ylabel('Counts')
        title_prefix = f"{plot_title} - " if plot_title else ""
        ax_zoom.set_title(title_prefix + f"Zoom on {zoom_peak}")
        plt.show()
        
        
