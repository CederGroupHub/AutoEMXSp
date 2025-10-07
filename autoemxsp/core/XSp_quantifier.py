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
from typing import Optional, Dict, Tuple, Sequence, List, Union

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
import autoemxsp.XSp_calibs as calibs 
import autoemxsp.tools.constants as cnst
from autoemxsp.tools.utils import (
    print_nice_1d_row,
    print_single_separator,
    print_double_separator,
    EDSError,
    weight_to_atomic_fr,
    atomic_to_weight_fr
)
from autoemxsp.lib.Xray_lines import get_el_xray_lines
from autoemxsp.lib.Xray_absorption_coeffs import xray_mass_absorption_coeff
from autoemxsp.lib.mean_ionization_potentials import J_df
from autoemxsp.core.XSp_fitter import (
    XSp_Fitter,
    Background_Model,
    Peaks_Model
)

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
        fitting_verbose=False
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

        self.verbose = verbose
        self.fitting_verbose = fitting_verbose
        
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
                print("Initial background scaling factor K could not be estimated.")
                print(f"Current beam energy of {self.beam_energy} keV too low for reliable high-energy background fitting.")
                print("Beam energy needs to be at least 7.5 keV.")
            return None
    
        high_energy_indices = self.energy_vals > en_thresh
    
        if self.verbose:
            print_double_separator()
            print(f"Fit of spectrum above {en_thresh} keV to get initial background scaling factor K...")
            print("Turned off particle morphology parameters to avoid affecting value of K.")
    
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
                        print(f"Found K = {K_val:.2f}")
                    else:
                        print("Failed to find initial K value")
            except Exception as e:
                if self.verbose:
                    print("An error occurred during the quick background fit for K estimation:")
                    print(f"{type(e).__name__}: {e}")
                K_val = None
        else:
            if self.verbose:
                print("No suitable high-energy data for K estimation.")
    
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
        """
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
                if energy > 2 and self.beam_energy / energy > 1.65
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
        iter_counter = 0 # Iteration counter
        max_iterations = 30  # Maximum allowed iterations

    
        # Get initial value for parameter 'K' before initializing the fitter (Iteration 0)
        initial_param_values = None
        if self.is_particle and self.fit_background:
            k_val = self.get_starting_K_val()
            iter_counter += 1
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
        iter_counter += 1
        try:
            fitted_params, weight_fractions, sample_Z = self._fit_quant_spectrum_iter(
                initial_par_vals=initial_param_values,
                f_tol=initial_fit_tolerance,
                n_iter=iter_counter
            )
            iter_counter += 1
        except Exception as e:
            is_fit_valid = False
            print("Fit and quantification iteration unsuccessful due to the following error:")
            print(f"{type(e).__name__}: {e}")
    
        # Iteratively fit and quantify to converge to a solution
        if is_fit_valid and fit_iteratively:
            w_fr_change_convergence = 0.0001 # ZAF corrections converge when change in mass fraction is less than 0.01%
            diff_mass_fractions = 1  # To monitor convergence

            # Normalize mass fractions
            prev_weight_fractions = self._normalise_mass_fractions(weight_fractions)
    
            while iter_counter < max_iterations and diff_mass_fractions > w_fr_change_convergence:
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
                iter_counter += 1
                prev_weight_fractions = norm_mass_fractions
    
            if self.verbose:
                print(f"Spectrum fitted with {iter_counter} iterations")
    
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
                self.print_quant_result(quant_result, sample_Z)
    
            # Get minimum background counts for reference peaks
            min_bckgrnd_ref_lines = self._get_min_bckgrnd_cnts_ref_quant_lines()
    
        else:
            quant_result = None
            min_bckgrnd_ref_lines = 0
        
        self.bad_quant_flag = bad_quant_flag
        
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
        quant_result[cnst.R_SQ_KEY] = round(
            1 - self.fit_result.residual.var() / np.var(self.spectrum_vals), 6
        )
    
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
            mean_std = next((std for std in std_vals_list if std[cnst.STD_ID_KEY] == cnst.STD_MEAN_ID_KEY), None)
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
                    f"No standard with '{cnst.STD_ID_KEY}' == '{cnst.STD_MEAN_ID_KEY}' found for {el_line} in the standards database."
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
            print('Quantification with ZAF correction:')
            print_single_separator()
            print_nice_1d_row('', self.fitted_els_quant)
            print_nice_1d_row('Initial W_fr', k_ratios)
            print(f"Initial analytical error: {(1 - sum(k_ratios)) * 100:.2f}%")
    
        while max_diff > converge_tol and ZAF_cntr < max_iter:
            ZAF_cntr += 1
            if self.verbose:
                print_single_separator()
                print(f"Step: {ZAF_cntr}")
    
            # Calculate ZAF factors and sample mean Z
            ZAF_pb_factors, sample_Z = correction.get_ZAF_mult_f_pb(weight_fractions)
    
            # Update weight fractions
            new_weight_fractions = k_ratios * ZAF_pb_factors
            if self.verbose:
                print_nice_1d_row('New W_fr', new_weight_fractions)
                print('Analytical error: %.2f w%%' % (sum(new_weight_fractions) * 100 - 100))
    
            max_diff = max(abs(new_weight_fractions - weight_fractions))
            weight_fractions = new_weight_fractions.copy()
    
        if ZAF_cntr == max_iter:
            print_single_separator()
            print(f'ZAF correction did not converge within {max_iter} iterations.')
        elif self.verbose:
            print_single_separator()
            print(f"ZAF correction converged in {ZAF_cntr} steps.")
    
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
                        Background_Model(True)  # Re-initialize absorption attenuation globals
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
        redchi_threshold = 0.2    # Threshold of reduced-chi squared value as % of total counts
        redchi_threshold_val = self.tot_sp_counts * redchi_threshold / 100
        an_err_threshold = 0.5    # Analytical error threshold (50 w%)
        abs_increase_threshold = 0.7  # Absorption increase threshold (170% of bulk absorption)
    
        bad_quant_flag = None  # Default: quantification is considered reliable
        abs_att_param_name = '_abs_attenuation_phirho' # TODO does not work for models different than 'phirho' absorption attenuation model
    
        # 1. Check for extremely poor fit (high reduced chi-squared)
        if self.fit_result.redchi > redchi_threshold_val:
            bad_quant_flag = 1
            if interrupt_fits_bad_spectra:
                print(f"Quantification stopped at iteration #{iter_cntr} due to reduced chi-squared being "
                      f"{self.fit_result.redchi:.1f} > {redchi_threshold_val:.1f}")
    
        # 2. Check for excessive analytical error
        elif analytical_er > an_err_threshold:
            bad_quant_flag = 2
            if interrupt_fits_bad_spectra:
                print(f"Quantification stopped at iteration #{iter_cntr} due to analytical error being "
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
                    print(f"Quantification stopped at iteration #{iter_cntr} due to absorption around reference peaks being "
                          f"{abs_increase*100:.1f}% > {abs_increase_threshold*100:.1f}%")
    
        return bad_quant_flag
    
    #%% Post-quantification functions
    # =============================================================================
    def print_quant_result(
        self,
        quant_result: dict,
        Z_sample: dict = None
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
        
        print('Fit result:')
        
        print(f"Reduced Chi-square: {quant_result[cnst.REDCHI_SQ_KEY]:.2f}")
        print(f"R-squared: {quant_result[cnst.R_SQ_KEY]:.5f}")
        
        print('')
        
        print('Quantification result:\n')
    
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
        print(f"Analytical error: {an_err_percent:.2f} w%\n")
    
        # Print mean atomic numbers (Z̅) only if provided
        if Z_sample is not None:
            if 'Statham2016' in Z_sample:
                print(f"Z̅_Statham2016: {Z_sample['Statham2016']:.2f}")
            if 'mass-averaged' in Z_sample:
                print(f"Z̅_w (mass-averaged): {Z_sample['mass-averaged']:.2f}")
            
    
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
            print(f"Warning: Unrecognized value for annotate_peaks ('{annotate_peaks}'). No peaks will be annotated.")
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
                print(f'You have attempted to zoom on a peak, using the line {peak}.')
                print('This line is absent from the list of fitted peaks, so the plot was not zoomed.')
                print(f'The available peak lines are {self.fitted_xray_lines}')
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
            print(f"Peak '{zoom_peak}' not found in fitted peaks.")
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
        
        
#%% Quantification Corrections class
class Quant_Corrections:
    """
    Implements matrix correction factors for quantitative X-ray microanalysis using the peak-to-background (P/B) method.

    This class provides methods for calculating Z (atomic number), A (absorption), and R (backscattering) correction factors,
    as well as mass absorption coefficients for a given set of elements and measurement conditions. It is designed for
    both standard-based and standardless quantification workflows in electron probe microanalysis (EPMA) or energy-dispersive
    X-ray spectroscopy (EDS).

    References
    ----------
    G. Love, V.D. Scott, "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
        J. Phys. D. Appl. Phys. 11 (1978) 1369–1376. https://doi.org/10.1088/0022-3727/11/10/002
    P.J. Statham, "A ZAF PROCEDURE FOR MICROPROBE ANALYSIS BASED ON MEASUREMENT OF PEAK-TO-BACKGROUND RATIOS",
        in: D.E. Newbury (Ed.), Fourteenth Annu. Conf. Microbeam Anal. Soc., San Francisco Press, 1979: pp. 247–253.
    M. Essani, E. Brackx, E. Excoffier, "A method for the correction of size effects in microparticles using a peak-to-background approach
        in electron-probe microanalysis", Spectrochim. Acta B 169 (2020) 105880. https://doi.org/10.1016/j.sab.2020.105880

    Attributes
    ----------
    elements : list of str
        List of element symbols included in the quantification (excluding undetectable elements).
    energies : np.ndarray
        X-ray energies (keV) for each element/line.
    emergence_angle : float
        Detector emergence angle (degrees).
    beam_energy : float
        Incident electron beam energy (keV).
    meas_mode : str
        Detector mode (for calibration parameters).
    Z_els : np.ndarray
        Atomic numbers for each element.
    W_els : np.ndarray
        Atomic weights for each element.
    els_nu : np.ndarray
        Backscattering coefficients for each element.
    mass_abs_coeffs_lines : list of list of float
        Mass absorption coefficients for each element at each characteristic energy.
    verbose : bool
        If True, enables verbose output for debugging.
    """

    def __init__(
        self,
        elements: Sequence[str],
        beam_energy: float,
        emergence_angle: float,
        meas_mode: str,
        energies: Optional[Union[Sequence[float], np.ndarray]] = None,
        verbose: bool = False
    ) -> None:
        """
        Initialize the Quant_Corrections class for matrix correction calculations.
    
        Parameters
        ----------
        elements : Sequence[str]
            List or sequence of element symbols to include in quantification (e.g., ['Fe', 'Si', 'O']).
        beam_energy : float
            Incident electron beam energy (keV).
        emergence_angle : float
            Detector emergence (take-off) angle (degrees).
        meas_mode : str
            EDS collection mode (used to retrieve calibration parameters).
        energies : Sequence[float] or np.ndarray, optional
            X-ray energies (keV) corresponding to each element/line.
            Generally provided when class is called from XSp_Quantifier.
            If not provided here, energy values must be passed directly to the functions later.
            This is done when measuring experimental standards.
        verbose : bool, optional
            If True, enables verbose output for debugging (default: False).
    
        Notes
        -----
        - Requires microscope calibrations to be loaded through XSp_calibs.load_microscope_calibrations(). This is done automatically
            when this class is called from XSp_Quantifier
        - All numeric arrays are stored as np.ndarray for consistency and performance.
        - Undetectable elements (as defined in `XSp_calibs.undetectable_els`) are automatically excluded from quantification.
        - Mass absorption coefficients are stored as a nested list, where each sub-list contains the coefficients for all
          elements at a given characteristic energy.
        - If `energies` is not provided at initialization, it must be set before using methods that require energy values.
        """
        # Ensure microscope calibrations have been loaded
        if not calibs.microscope_calibrations_loaded:
            raise EDSError("Microscope calibrations have not been loaded."
                           "Ensure the class XSp_Quantifier is initialised before instancing Quant_Corrections."
                           "Alternatively, load calibrations through XSp_calibs.load_microscope_calibrations() first.")
        
        # Filter out undetectable elements and their corresponding energies (if energies provided)
        detectable_mask = [el not in calibs.undetectable_els for el in elements]
        quant_elements = [el for el, keep in zip(elements, detectable_mask) if keep]
    
        if energies is not None:
            quant_energies = [en for en, keep in zip(energies, detectable_mask) if keep]
            self.energies = np.array(quant_energies, dtype=float)
        else:
            self.energies = None

        self.sample_elements = quant_elements
        self.beam_energy = beam_energy
        self.emergence_angle = emergence_angle
        self.meas_mode = meas_mode

        # Atomic numbers and weights for each element
        Z_els = []
        W_els = []
        for el in quant_elements:
            Z_els.append(Element(el).Z)
            W_els.append(Element(el).atomic_mass)
        self.Z_els = np.array(Z_els)
        self.W_els = np.array(W_els)
        
        # ---- Precalculate fixed attributes ----

        # Backscattering coefficients for all elements (vectorized for all quantifiable elements)
        self.els_nu: np.ndarray = self._nu(self.Z_els)
        
        self.mass_abs_coeffs_lines = None # Initialise for computation at first iteration
        
        self.verbose = verbose

    # =============================================================================
    # Main function
    # =============================================================================
    def get_ZAF_mult_f_pb(
        self,
        weight_fractions: np.ndarray,
        el_lines_energies_d: Optional[Dict[str, float]] = None
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Calculate the ZAF multiplicative correction factors for the measured sample P/B ratio.
    
        This method accounts for:
            1. Differences in average Z between the sample and the employed standard,
               which affect continuum intensity.
            2. Second-order corrections for backscattering and absorption due to differential
               mean generation path between characteristic and continuum X-rays.
        Fluorescence and particle-size corrections are ignored.
        
    
        Parameters
        ----------
        weight_fractions : np.ndarray
            Estimated weight fractions of the elements to quantify.
        el_lines_energies_d : dict[str, float], optional
            Dictionary mapping peak labels (e.g., 'Fe_Ka') to energies (keV).
            If None, uses all elements and energies in the class.
    
        Returns
        -------
        ZAF_pb_corrections : np.ndarray
            ZAF correction factors to multiply with the measured sample P/B ratio.
        Z_sample : dict
            Dictionary with various sample mean atomic numbers.
            
        References
        ----------
        [1] Statham P, Penman C, Duncumb P. Improved spectrum simulation for validating SEM-EDS analysis.
            IOP Conf. Ser. Mater. Sci. Eng. 109, 0 (2016).
        [2] Markowicz AA, Van Grieken RE. Composition dependence of bremsstrahlung background in electron-probe
            x-ray microanalysis. Anal. Chem. 56, 2049 (1984).
            
        Potential improvements
        ----------------------
        Include fluorescence corrections for large particles
        Include particle size corrections, from:
            [1] J. L. Lábár and S. Török, A peak‐to‐background method for electron‐probe x‐ray microanalysis
                applied to individual small particles, X-Ray Spectrom. 21, 183 (1992).
            [2] M. Essani, E. Brackx, and E. Excoffier, A method for the correction of size effects in microparticles
                using a peak-to-background approach in electron-probe microanalysis,
                Spectrochim. Acta - Part B At. Spectrosc. 169, 105880 (2020).
        """
        # Convert mass fractions to atomic fractions
        atomic_fractions = weight_to_atomic_fr(weight_fractions, self.sample_elements, verbose=False)
    
        # Normalize mass fractions to avoid divergence in ZAF algorithm
        norm_weight_fractions = atomic_to_weight_fr(atomic_fractions, self.sample_elements)
        
        # Calculate average Z in the sample using different conventions
        Z_sample_w = float(np.sum(norm_weight_fractions * self.Z_els))
        Z_sample_at = float(np.sum(atomic_fractions * self.Z_els))
        Z_sample_Markowicz = float(self._Z_mean_Markowicz1984(atomic_fractions, norm_weight_fractions))
        Z_sample_Statham = float(self._Z_mean_Statham2016(atomic_fractions, norm_weight_fractions))
    
        Z_sample = {
            cnst.Z_MEAN_W_KEY : Z_sample_w,
            cnst.Z_MEAN_AT_KEY : Z_sample_at,
            cnst.Z_MEAN_STATHAM_KEY : Z_sample_Statham,
            cnst.Z_MEAN_MARKOWICZ_KEY : Z_sample_Markowicz
        }
        
        if el_lines_energies_d is not None:
            energies = np.array(list(el_lines_energies_d.values()))
        else:
            energies = None
        
        # Calculate Z, A, and R multiplicative factors for the sample
        Z_vals = self._Z_pb(Z_sample_Statham, norm_weight_fractions, el_lines_energies_d)
        A_vals = self._A_pb(norm_weight_fractions, energies)
        R_vals = self._R_pb(norm_weight_fractions, energies)

        # ZAF multiplicative factor
        ZAF_pb_corrections = Z_vals * A_vals * R_vals
    
        if self.verbose:
            # Print header row with element names
            print_nice_1d_row('', self.sample_elements)
            # Print data rows with appropriate labels
            print_nice_1d_row('At_fr', atomic_fractions)
            print_nice_1d_row('W_fr', weight_fractions)
            print_nice_1d_row('Z_vals', Z_vals)
            print_nice_1d_row('A_vals', A_vals)
            print_nice_1d_row('R_vals', R_vals)
            print_nice_1d_row('Z·A·R', ZAF_pb_corrections)
    
        return ZAF_pb_corrections, Z_sample
    
    
    def _get_energy_vals(self) -> np.ndarray:
        """
        Retrieve the array of X-ray energies used for quantification.
    
        Returns
        -------
        np.ndarray
            Array of X-ray energies (in keV) for each line.
    
        Raises
        ------
        ValueError
            If the energies attribute is not set or is None.
    
        Notes
        -----
        This method ensures that the object has a valid 'energies' attribute before returning it.
        """
        if not hasattr(self, 'energies') or self.energies is None:
            raise ValueError("No energies provided and self.energies is not set.")
        return self.energies
    
    # =============================================================================
    # Atomic number averaging
    # =============================================================================
    def _Z_mean_Markowicz1984(
        self,
        at_frs: Sequence[float],
        w_frs: Sequence[float]
    ) -> float:
        """
        Calculate the average atomic number (Z) in the sample using the Markowicz method, as described in:
        Markowicz AA, Van Grieken RE. "Composition dependence of bremsstrahlung background in electron-probe x-ray microanalysis."
        Anal. Chem. 1984, 56(12), 2049–2051. https://pubs.acs.org/doi/abs/10.1021/ac00276a016
    
        Parameters
        ----------
        at_frs : Sequence[float]
            Atomic fractions of elements in the sample.
        w_frs : Sequence[float]
            Weight fractions of elements in the sample.
    
        Returns
        -------
        Z_mean : float
            The Markowicz mean atomic number for the sample.
        """
        Z_num = 0.0  # Numerator of Markowicz expression
        Z_den = 0.0  # Denominator of Markowicz expression
    
        for el_Z, el_A, w_fr, at_fr in zip(self.Z_els, self.W_els, w_frs, at_frs):
            Z_num += w_fr * el_Z**2 / el_A
            Z_den += w_fr * el_Z / el_A
    
        Z_mean = Z_num / Z_den
        return Z_mean

    
    def _Z_mean_Statham2016(
        self,
        at_frs: Sequence[float],
        w_frs: Sequence[float]
    ) -> float:
        """
        Calculate the average atomic number (Z) in the sample using the Statham method, as described in:
    
        This method implements the mean Z calculation as described in:
            Statham P, Penman C, Duncumb P. "Improved spectrum simulation for validating SEM-EDS analysis."
            IOP Conf Ser Mater Sci Eng. 2016;109(1):0–10.
        
        This formula is practically the same as in, except for the exponent being 0.7 instead of 0.75:
            - J. J. Donovan and N. E. Pingitore, Compositional Averaging of Continuum Intensities in
            Multielement Compounds, Microsc. Microanal. 8, 429 (2002).
            - J. Donovan, A. Ducharme, J. J. Schwab, A. Moy, Z. Gainsforth, B. Wade, and B. McMorran,
            An Improved Average Atomic Number Calculation for Estimating Backscatter and Continuum
            Production in Compounds, Microsc. Microanal. 29, 1436 (2023).
            
        Parameters
        ----------
        at_frs : Sequence[float]
            Atomic fractions of elements in the sample.
        w_frs : Sequence[float]
            Weight fractions of elements in the sample.
    
        Returns
        -------
        Z_mean : float
            The Statham mean atomic number for the sample.
        """
        Z_num = 0.0  # Numerator of Statham expression
        Z_den = 0.0  # Denominator of Statham expression
    
        for el_Z, el_A, w_fr, at_fr in zip(self.Z_els, self.W_els, w_frs, at_frs):
            Z_num += w_fr * el_Z ** 1.75 / el_A
            Z_den += w_fr * el_Z ** 0.75 / el_A
    
        Z_mean = Z_num / Z_den
        return Z_mean
    
    # =============================================================================
    # Continuum intensity atomic number correction Z_c
    # =============================================================================
    def _Z_pb(
        self,
        Z_sample_Statham: float,
        norm_weight_fractions,
        el_lines_energies_d: Optional[Dict[str, float]] = None
    ):
        """
        Calculate generated continuum values for pure elements and for the sample composition,
        and return the Z factor as used in the standard P/B correction.
    
        Parameters
        ----------
        Z_sample_Statham : float
            Average atomic number of the sample (Statham method).
        norm_weight_fractions : array-like
            Normalized mass fractions of each element in the sample.
        el_lines_energies_d : dict[str, float], optional
            Dictionary mapping peak labels (e.g., 'Fe_Ka') to energies (keV).
            If None, uses all elements and energies in the class.
    
        Returns
        -------
        Z_vals : np.ndarray
            Z factor values (sample/standard continuum ratio).
        gen_bckgrnd_vals_sample : np.ndarray
            Generated continuum values for the sample composition.
        gen_bckgrnd_vals_pure_els : np.ndarray
            Generated continuum values for pure elements (standards).
        """
        # Calculate values of generated continuum for pure elements, which the standard PB values refer to
        if el_lines_energies_d is None:
            # Case: applies to all elements in the class
            ens = self._get_energy_vals()
            Z_els = self.Z_els
            W_els = self.W_els
        else:
            # Case: el_lines_energies_d is a dict of {peak_label: energy}
            ens, Z_els, W_els = [], [], []
            for peak_label, en in el_lines_energies_d.items():
                el = peak_label.split('_')[0]
                if el not in self.sample_elements:
                    raise ValueError(f"Element {el} not found in self.sample_elements.")
                index_el = self.sample_elements.index(el)
                ens.append(en)
                Z_els.append(self.Z_els[index_el])
                W_els.append(self.W_els[index_el])
    
        gen_bckgrnd_vals_pure_els = [
            self._gen_bckgrnd_vals(Z_el, 1.0, en, Z_el, W_el)[0]
            for en, Z_el, W_el in zip(ens, Z_els, W_els)
        ]

        # Calculate values of generated continuum for the sample composition, calculated at energies ens
        gen_bckgrnd_vals_sample = self._gen_bckgrnd_vals(
            Z_sample_Statham, norm_weight_fractions, ens, self.Z_els, self.W_els
        )

        # Calculate Z_c
        Z_vals = gen_bckgrnd_vals_sample / gen_bckgrnd_vals_pure_els
    
        return Z_vals
    
    
    def _gen_bckgrnd_vals(
            self,
            Z_sample: float,
            weight_fractions: Union[float, Sequence[float]],
            energies: Union[float, Sequence[float]],
            Z_els: Union[float, Sequence[float]],
            W_els: Union[float, Sequence[float]]
        ) -> np.ndarray:
            """
            Compute the generated continuum background to calculate the Z correction (Z_c) 
            in the P/B method for quantitative electron probe microanalysis.
    
            Z_c accounts for differences continuum intensity arising from differences in mean
            atomic number (Z) between the measured sample and the standard composition.
    
            Parameters
            ----------
            Z_sample : float
                Average atomic number of the sample.
            weight_fractions : float or Sequence[float]
                Mass fractions of each element in the sample.
            energies : float or Sequence[float]
                X-ray energies (keV) for each element/line.
            Z_els : float or Sequence[float]
                Atomic numbers for each element.
            W_els : float or Sequence[float]
                Atomic weights for each element.
    
            Returns
            -------
            np.ndarray
                Generated background values, free of matrix composition effects.
    
            References
            ----------
            Stopping power correction from:
            G. Love, V.D. Scott, Evaluation of a new correction procedure for quantitative electron probe microanalysis,
                J. Phys. D. Appl. Phys. 11 (1978) 1369–1376. https://doi.org/10.1088/0022-3727/11/10/002
            """
            # Ensure all inputs are arrays for vectorized operations
            weight_fractions = np.atleast_1d(weight_fractions).astype(np.float64)
            energies = np.atleast_1d(energies).astype(np.float64)
            Z_els = np.atleast_1d(Z_els).astype(np.float64)
            W_els = np.atleast_1d(W_els).astype(np.float64)
    
            # Initialise background model
            bckgrnd = Background_Model(
                is_particle=False,
                beam_energy=self.beam_energy,
                emergence_angle=self.emergence_angle,
                meas_mode=self.meas_mode
            )
    
            # Get generated background calibrated parameters
            Z = sp.Symbol('Z')
            P_expr, F_expr, beta_expr = sp.sympify(calibs.get_calibrated_background_params(self.meas_mode))
            P_val = float(P_expr.subs(Z, Z_sample).evalf())
            F_val = float(F_expr.subs(Z, Z_sample).evalf())
            beta_val = 0.0
            for el_Z, w_fr in zip(Z_els, weight_fractions):
                beta_component = float(beta_expr.subs(Z, el_Z).evalf())
                beta_val += beta_component * w_fr
    
            # Compute generated background value using Duncumb modification
            mod_Duncumb_gen_bckgrnd = bckgrnd._generated_bckgrnd_DuncumbMod(
                energies, Z=Z_sample, P=P_val, F=F_val, beta=beta_val, apply_det_response=0
            )
            mod_Duncumb_gen_bckgrnd = np.asarray(mod_Duncumb_gen_bckgrnd, dtype=np.float64)
    
            # Stopping power correction (Love & Scott 1978)
            J_els = np.array([J_df.loc[Z_el, J_df.columns[0]] / 1000 for Z_el in Z_els], dtype=np.float64)  # Mean ionization potential J (keV)
            sum_M = np.sum(weight_fractions * Z_els / W_els)
            ln_J = np.sum(weight_fractions * Z_els / W_els * np.log(J_els)) / sum_M
            J_val = np.exp(ln_J)
            U0 = self.beam_energy / energies
            S_vals = (1 + 16.05 * (J_val / energies) ** 0.5 * ((U0 ** 0.5 - 1) / (U0 - 1)) ** 1.07) / sum_M
    
            # Final generated background value, rid of matrix composition effect
            gen_background_vals = mod_Duncumb_gen_bckgrnd / S_vals
    
            return gen_background_vals

    # =============================================================================
    # Absorption attenuation corrections
    # =============================================================================       
    def _get_mass_abs_coeffs_sample(
        self,
        weight_fractions: np.ndarray,
        energies: np.ndarray,
    ) -> np.ndarray:
        """
        Calculate the mass absorption coefficients of the sample,
        using the defined weight fractions and, if provided, the specified energies.
    
        Parameters
        ----------
        weight_fractions : np.ndarray
            Array of mass fractions for each element in the sample.
        energies : np.ndarray
            Array of X-ray energies (keV) at which mass absorption coefficient is computed
    
        Returns
        -------
        np.ndarray
            Mass absorption coefficients for each line energy, weighted by the sample composition.
    
        Notes
        -----
        If self.mass_abs_coeffs_lines is not already set, it will be calculated on the fly
        for the provided energies and sample elements.
        """
        # Compute (first iteration) or retrieve mass absorption coefficients for each line/element
        if getattr(self, 'mass_abs_coeffs_lines', None) is not None:
            mass_abs_coeffs_lines = self.mass_abs_coeffs_lines
        else:
            # First iteration, compute mass absorption coefficients for each element at each energy value.
            # Structure: nested list for computation efficiency.
            # Each sub-list contains the mass absorption coefficients of all elements at each value of energy:
            # e.g. if each energy value corresponds to a characteristic line:
            #   [ [mu_Fe@FeKa, mu_Si@FeKa, ...], [mu_Fe@SiKa, mu_Si@SiKa, ...], ... ]
            # Indices follow the order of: energies, elements.
            mass_abs_coeffs_lines = [
                [xray_mass_absorption_coeff(el, en) for el in self.sample_elements]
                for en in energies
            ]
            self.mass_abs_coeffs_lines: List[List[float]] = mass_abs_coeffs_lines

        mass_abs_coeffs_lines = np.asarray(mass_abs_coeffs_lines)
    
        # Weighted sum to get sample mass absorption coefficients
        mass_abs_coeffs_sample = np.dot(mass_abs_coeffs_lines, weight_fractions)
        return mass_abs_coeffs_sample

        
    def _A_pb(
        self,
        weight_fractions: np.ndarray,
        energies: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Second-order absorption correction for the P/B ratio to account for differences in mean generation depths
        of characteristic x-rays and continuum.
    
        This correction factor (A_c) should be multiplied by the measured P/B to obtain the P/B without matrix effects from absorption.
        Because the depth of generation of the continuum is larger than that of characteristic x-rays,
        the continuum is absorbed more. Thus, A_pb < 1.
    
        Reference
        ---------
        P.J. Statham, "A ZAF PROCEDURE FOR MICROPROBE ANALYSIS BASED ON MEASUREMENT OF PEAK-TO-BACKGROUND RATIOS",
        in: D.E. Newbury (Ed.), Fourteenth Annu. Conf. Microbeam Anal. Soc., San Francisco Press, San Francisco, 1979: pp. 247–253.
        https://archive.org/details/1979-mas-proc-san-antonio/page/246/mode/2up
    
        Parameters
        ----------
        mass_abs_coeffs_sample : np.ndarray
            Mass absorption coefficients for the sample (for each line).
        energies : np.ndarray, optional
            X-ray energies (keV) at which A_pb is computed, corresponding to the quantified line energies.
            If not provided, self.energies are used instead.
    
        Returns
        -------
        np.ndarray
            Absorption correction factors (A_c), to be multiplied with measured P/B ratios.
            
        Raises
        ------
        ValueError
            If neither `energies` nor `self.energies` are provided.
        """
        if energies is None:
            energies = self._get_energy_vals()
        else:
            energies = np.asarray(energies)
            
        mass_abs_coeffs_sample = self._get_mass_abs_coeffs_sample(weight_fractions, energies)
        
        # Convert emergence angle to radians for np.sin
        emergence_angle_rad = np.deg2rad(self.emergence_angle)
        chi = mass_abs_coeffs_sample / np.sin(emergence_angle_rad)
        gamma = (self.beam_energy ** 1.65 - np.asarray(energies) ** 1.65)
        x = chi * gamma
    
        # Absorption fraction for characteristic X-ray
        f_char = 1 / (1 + 3.0e-6 * x + 4.5e-13 * x ** 2)
        # Absorption fraction for continuum (higher than characteristic X-rays)
        f_cont = 1 / (1 + 3.34e-6 * x + 5.59e-13 * x ** 2)
    
        # Multiplicative factor for PB ratio
        A_c = f_char / f_cont
    
        return A_c


    # =============================================================================
    # Backscattering electron corrections
    # =============================================================================    
    def _R_pb(self, weight_fractions, energies: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Second-order backscattering correction for P/B ratio to account for differences in mean generation depths of
        characteristic x-rays and continuum.
    
        This correction factor (R_pb) should be multiplied by the measured P/B to obtain the P/B without matrix effects.
        Because the depth of generation of the continuum is larger than that of characteristic x-rays,
        the continuum loses more intensity due to backscattering. Thus, R_pb < 1.
    
        Reference
        ---------
        P.J. Statham, "A ZAF PROCEDURE FOR MICROPROBE ANALYSIS BASED ON MEASUREMENT OF PEAK-TO-BACKGROUND RATIOS",
        in: D.E. Newbury (Ed.), Fourteenth Annu. Conf. Microbeam Anal. Soc., San Francisco Press, San Francisco, 1979: pp. 247–253.
        https://archive.org/details/1979-mas-proc-san-antonio/page/246/mode/2up
    
        Parameters
        ----------
        weight_fractions : array-like
            Mass fractions of each element in the sample.
        energies : np.ndarray, optional
            X-ray energies (keV) at which A_pb is computed, corresponding to the quantified line energies.
            If not provided, self.energies are used instead.
    
        Returns
        -------
        np.ndarray
            Backscattering correction factors (R_pb), to be multiplied with measured P/B ratios.
            
        Raises
        ------
        ValueError
            If neither `energies` nor `self.energies` are provided.
        """
        if energies is None:
            energies = self._get_energy_vals()
        else:
            energies = np.asarray(energies)
        
        # Backscattering correction for characteristic X-ray
        R_P_vals = self._R_p(weight_fractions, energies=energies)
        # Backscattering correction for continuum
        R_B_vals = self._R_b(weight_fractions, R_P_vals, energies=energies)
        # Multiplicative factor for PB ratio
        R_vals = R_P_vals / R_B_vals
        return R_vals
    
    
    def _R_b(self, weight_fractions: np.ndarray, R_P_vals: np.ndarray, energies: np.ndarray) -> np.ndarray:
        """
        Statham's formula for second-order backscattering correction for the P/B ratio to account for
        differences in mean generation depths of characteristic x-rays and continuum.
    
        The parameter 'nu' is averaged by weighting on the mass fractions, according to Love (1978).
    
        References
        ----------
        M. Essani, E. Brackx, E. Excoffier,
        "A method for the correction of size effects in microparticles using a peak-to-background approach
        in electron-probe microanalysis", Spectrochim. Acta B 169 (2020) 105880.
        https://doi.org/10.1016/j.sab.2020.105880
    
        G. Love, V.D. Scott, "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
        J. Phys. D. Appl. Phys. 11 (1978) 1369–1376.
    
        Parameters
        ----------
        weight_fractions : np.ndarray
            Mass fractions of each element in the sample.
        R_P_vals : np.ndarray
            Backscattering correction factors for characteristic X-ray lines.
        energies : np.ndarray
            X-ray energies (keV) at which A_pb is computed, corresponding to the quantified line energies.
    
        Returns
        -------
        np.ndarray
            Backscattering correction factors for continuum (R_B).
        """
        # Weighted average of nu for the sample (Love, 1978)
        nu_sample = np.sum(weight_fractions * self.els_nu)
        
        # Statham/Essani formula for continuum backscattering correction
        factor_Statham = (2 / (1 + nu_sample)) ** 0.63 * (0.79 + 0.44 * energies / self.beam_energy)
        R_B_vals = 1 - (1 - R_P_vals) * factor_Statham
        
        return R_B_vals

    
    def _R_p(self, weight_fractions: np.ndarray, energies: np.ndarray) -> np.ndarray:
        """
        Backscattering correction factor for characteristic X-rays.
    
        The parameter 'nu' is averaged by weighting on the mass fractions, according to Love (1978).
    
        References
        ----------
        M. Essani, E. Brackx, E. Excoffier,
        "A method for the correction of size effects in microparticles using a peak-to-background approach
        in electron-probe microanalysis", Spectrochim. Acta B 169 (2020) 105880.
        https://doi.org/10.1016/j.sab.2020.105880
    
        G. Love, V.D. Scott, "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
        J. Phys. D. Appl. Phys. 11 (1978) 1369–1376.
    
        Parameters
        ----------
        weight_fractions : np.ndarray
            Mass fractions of each element in the sample.
        energies : np.ndarray
            X-ray energies (keV) at which A_pb is computed, corresponding to the quantified line energies.
    
        Returns
        -------
        np.ndarray
            Backscattering correction factors for characteristic X-rays (R_p).
        """
        # Weighted average of nu for the sample (Love, 1978)
        nu_sample = np.sum(weight_fractions * self.els_nu)
    
        I_vals, G_vals = self._return_IG(energies)
    
        # Compute the correction factor
        R_p_vals = 1 - nu_sample * (I_vals + nu_sample * G_vals) ** 1.67
    
        return R_p_vals
    
    
    def _nu(self, Z_vals: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate the backscattering coefficient (nu) for a given atomic number or array of atomic numbers.
    
        For a compound, nu should be averaged over the constituent elements, weighted by their mass fractions.
    
        Reference
        ---------
        G. Love, V.D. Scott, "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
        J. Phys. D. Appl. Phys. 11 (1978) 1369–1376. https://doi.org/10.1088/0022-3727/11/10/002
    
        Parameters
        ----------
        Z_vals : float or np.ndarray
            Atomic number(s) for which to calculate the backscattering coefficient.
    
        Returns
        -------
        float or np.ndarray
            Backscattering coefficient(s) (nu) for the given atomic number(s).
        """
        Z = np.asarray(Z_vals)
        nu20 = (-52.3791 + 150.48371 * Z - 1.67373 * Z ** 2 + 0.00716 * Z ** 3) * 1e-4
        G_nu20 = (-1112.8 + 30.289 * Z - 0.15498 * Z ** 2) * 1e-4
        nu_vals = nu20 * (1 + G_nu20 * np.log(self.beam_energy / 20))
        return nu_vals

    
    def _return_IG(self, energies: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate the I and G functions of overvoltage needed for backscattering correction.
    
        Reference
        ---------
        G. Love, V.D. Scott, "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
        J. Phys. D. Appl. Phys. 11 (1978) 1369–1376.
        https://doi.org/10.1088/0022-3727/11/10/002
    
        Parameters
        ----------
        energies : np.ndarray
            X-ray energies (keV) at which A_pb is computed, corresponding to the quantified line energies.
    
        Returns
        -------
        I_vals : np.ndarray
            I function values for each energy line.
        G_vals : np.ndarray
            G function values for each energy line.
        """
        U0 = self.beam_energy / energies
        log_U0 = np.log(U0)
        I_vals = (
            0.33148 * log_U0
            + 0.05596 * log_U0 ** 2
            - 0.06339 * log_U0 ** 3
            + 0.00947 * log_U0 ** 4
        )
        G_vals = (
            1 / U0
            * (
                2.87898 * log_U0
                - 1.51307 * log_U0 ** 2
                + 0.81312 * log_U0 ** 3
                - 0.08241 * log_U0 ** 4
            )
        )
        return I_vals, G_vals