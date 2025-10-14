#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
X-ray Spectrum Fitting Module

Created on Thu Jun 27 10:17:26 2024

@author: Andrea Giunto

This module provides classes and functions for physically-accurate modeling, fitting, and analysis of X-ray energy-dispersive spectroscopy (EDS) spectra.

Overview
--------
The module is designed to enable robust, quantitative analysis of EDS spectra, including background continuum, detector response, and detailed peak shapes (including escape and pileup peaks, ICC effects, and low-energy tails). It is suitable for both bulk and particle samples, and supports flexible calibration and constraint schemes.

Class Structure and Interactions
-------------------------------
The main classes are:
- **XSp_Fitter**  
  The main workflow class. Given a measured spectrum, calibration information, and a list of elements, it builds the complete lmfit model (background + peaks), sets up all parameters and constraints, and runs the fit. It provides methods for plotting, reporting, and extracting fit results. This class orchestrates the use of all other classes and should be the main entry point for typical users.

- **Peaks_Model**  
  Manages the construction and parameterization of all spectral peaks (characteristic X-ray lines, escape peaks, pileup peaks, etc.). It supports advanced peak shapes (e.g., skewed/tail Gaussians, ICC convolution), constraints between related peaks, and caching for efficient repeated use. It is typically instantiated by the spectrum fitter and used to build the composite peak model.

- **Background_Model**  
  Handles the computation and parameterization of the spectral background, including physical effects such as X-ray generation, absorption, detector efficiency, and backscattering. Used to build the background component of the overall spectral model.

- **DetectorResponseFunction**  
  Provides static and class methods for handling the detector's instrumental response, including convolution matrices for energy resolution and incomplete charge collection (ICC). This class is initialized with calibration data and is used by both background and peak models to accurately simulate detector effects.

Typical Usage
-------------
1. **Initialize the fitter:**
   ```python
   fitter = XSp_Fitter(
       spectrum_vals, energy_vals, els_to_quantify=['Fe', 'Ni'], microscope_ID='PhenomXL'
   )
   
2. **Fit the spectrum:**
   ```python
   fit_result, fitted_lines = fitter.fit_spectrum(plot_result=True, print_result=True)
   )
    
3. **Inspect and use results:**
   Use fit_result for detailed analysis.
   Plot or print results using fitter.plot_result() and fitter.print_result().
   Access fitted parameters, background components, and diagnostic information.
   
Customization & Calibration
---------------------------
Detector calibration, physical constants, and peak shape calibration are handled via the calibs module and are loaded automatically based on the specified microscope and EDS mode.
Advanced users may customize which peaks are freely calibrated, which are constrained, and how background/peak models are parameterized by modifying the relevant class parameters or by subclassing.

Dependencies
------------
numpy, scipy, lmfit, matplotlib, and supporting modules for calibration and physical constants.

**How the classes interact:**
------------------------
- `XSp_Fitter` is the main user-facing class. It creates and coordinates instances of `Background_Model` and `Peaks_Model`, and uses `DetectorResponseFunction` to ensure all detector effects are handled consistently.
- `DetectorResponseFunction` is a utility class used by both `Background_Model` and `Peaks_Model` to convolve model components with the detector response.
- `Peaks_Model` and `Background_Model` each build their respective parts of the overall spectrum model, which are then combined by the fitter for the full fit.

**In short:**  
---------
Instantiate `XSp_Fitter` with your data and settings, then call `fit_spectrum()`. The module will handle background, detector response, and peak modeling for you, providing a comprehensive, physically-based EDS spectrum fit.
"""

# =============================================================================
# Standard library imports
# =============================================================================
import os
import re
import time
import json
import warnings
from itertools import combinations

# =============================================================================
# Third-party library imports
# =============================================================================
import numpy as np
from pathlib import Path
import sympy as sp
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.signal import find_peaks, peak_prominences
from scipy.integrate import quad, trapezoid
from scipy.optimize import root_scalar
from lmfit import Model, Parameters, Parameter
from lmfit.models import GaussianModel
from pymatgen.core import Element

# =============================================================================
# Local module imports
# =============================================================================
from autoemxsp.tools.utils import (
    RefLineError, print_single_separator, print_double_separator,
    weight_to_atomic_fr, load_msa
)
import autoemxsp.tools.constants as cnst
import autoemxsp.XSp_calibs as calibs 
from autoemxsp.lib.Xray_lines import get_el_xray_lines
from autoemxsp.lib.Xray_absorption_coeffs import xray_mass_absorption_coeff
from autoemxsp.lib.mean_ionization_potentials import J_df

parent_dir = str(Path(__file__).resolve().parent.parent)

#%% XSp_Fitter class
class XSp_Fitter:
    """
    Fitter for EDS spectra.

    Handles EDS spectral fitting, including background modeling, element quantification,
    and correction for experimental conditions.

    Attributes
    ----------
    spectrum_vals : array-like
        Measured spectrum intensity values.
    energy_vals : array-like
        Corresponding energy values (in keV).
    els_to_quantify : list of str
        Elements to quantify in the sample.
    els_w_fr : dict
        Elements with fixed mass fractions.
    els_substrate : list of str
        Elements present in the substrate.
    fit_background : bool
        Whether to fit a background continuum.
    xray_quant_ref_lines : tuple of str
        X-ray lines used for quantification.
    is_particle : bool
        If True, fit considers absorption and mass effect of particles.
    microscope_ID : str
        Identifier for microscope calibration and detector efficiency.
    meas_mode : str
        EDS acquisition mode.
    spectrum_lims : tuple
        Spectrum limits.
    force_fr_total : bool
        Normalize total fitted elemental fraction to 1 if True.
    beam_energy : float
        Beam energy in keV.
    emergence_angle : float
        X-ray emergence angle in degrees.
    tot_sp_counts : int or None
        Total spectrum counts.
    sp_collection_time : float or None
        Spectrum collection time in seconds.
    print_evolving_params : bool
        Print evolving fit parameters (for debugging).
    verbose : bool
        Verbose output.
    """
    
    # Suffixes for escape and pile-up peaks
    escape_peaks_str = '_escSiKa'
    pileup_peaks_str = '_pileup'

    def __init__(
        self,
        spectrum_vals,
        energy_vals,
        spectrum_lims,
        microscope_ID,
        meas_mode,
        det_ch_offset,
        det_ch_width,
        beam_e,
        emergence_angle,
        fit_background=True,
        is_particle=False,
        els_to_quantify=None,
        els_substrate=None,
        els_w_fr=None,
        force_fr_total=True,
        tot_sp_counts=None,
        sp_collection_time=None,
        xray_quant_ref_lines=None,
        print_evolving_params=False,
        verbose=False
    ):
        """
        Initialize the EDS spectrum fitter.

        Parameters
        ----------
        spectrum_vals : array-like
            Measured spectrum intensity values.
        energy_vals : array-like
            Corresponding energy values (in keV).
        spectrum_lims : tuple of int
            Tuple specifying the start and end indices for the spectrum region to analyze.
        microscope_ID : str, optional
            Identifier for microscope calibration and detector efficiency.
        meas_mode : str, optional
            EDS acquisition mode.
        det_ch_offset : float
            Detector channel energy offset (keV).
        det_ch_width : float
            Detector channel width (keV).
        beam_e : float, optional
            Beam energy in keV.
        emergence_angle : float, optional
            X-ray emergence angle in degrees.
        fit_background : bool, optional
            Whether to fit a background continuum (default: True).
            If False, the spectrum_vals have to be provided stripped of the background.
        is_particle : bool, optional
            If True, fit considers absorption and mass effect of particles.
        els_to_quantify : list of str or None, optional
            Elements to quantify in the sample.
        els_substrate : list of str or None, optional
            Elements present in the substrate affecting spectrum (default: ['C', 'O', 'Al']).
        els_w_fr : dict or None, optional
            Elements with fixed mass fractions, to be used when fitting samples with known
            composition (e.g., standards).
        force_fr_total : bool, optional
            Normalize total fitted elemental fraction to 1 (default: True).
            Set to False if undetectable elements are present in the sample.
        tot_sp_counts : int or None, optional
            Total spectrum counts.
        sp_collection_time : float or None, optional
            Spectrum collection time in seconds.
        xray_quant_ref_lines : tuple of str or None, optional
            X-ray lines used for quantification (default: ('Ka1', 'La1', 'Ma1')).
        print_evolving_params : bool, optional
            Print evolving fit parameters (for debugging).
        verbose : bool, optional
            Verbose output (default: False).
        """
        # Handle mutable default arguments
        if els_to_quantify is None:
            els_to_quantify = []
        if els_substrate is None:
            els_substrate = ['C', 'O', 'Al']
        if els_w_fr is None:
            els_w_fr = {}
        if xray_quant_ref_lines is None:
            xray_quant_ref_lines = ('Ka1', 'La1', 'Ma1')

        # Input spectral data
        self.spectrum_vals = spectrum_vals
        self.energy_vals = energy_vals
        self.tot_sp_counts = tot_sp_counts
        self.sp_collection_time = sp_collection_time

        # Load microscope calibration parameters
        self.microscope_ID = microscope_ID
        self.meas_mode = meas_mode
        calibs.load_microscope_calibrations(microscope_ID, meas_mode, load_detector_channel_params=False)

        # Remove duplicates and undetectable elements from quantification and substrate lists
        self.els_to_quantify = [el for el in dict.fromkeys(els_to_quantify) if el not in calibs.undetectable_els]
        self.els_substrate = [el for el in dict.fromkeys(els_substrate)
                              if el not in calibs.undetectable_els and el not in self.els_to_quantify]

        # Elements with fixed mass fraction (e.g., for standards)
        self.els_w_fr = {el: w_fr for el, w_fr in els_w_fr.items() if el not in calibs.undetectable_els}
        self.force_fr_total = force_fr_total

        # List of all elements present in the spectrum (sample + substrate)
        self.els_to_fit_list = self.els_to_quantify + self.els_substrate
        self.num_els = len(self.els_to_fit_list)

        # EDS acquisition and geometry parameters
        self.emergence_angle = emergence_angle
        self.beam_energy = beam_e

        # X-ray lines used as references for dependent peaks
        self.xray_quant_ref_lines = xray_quant_ref_lines

        # If True, account for absorption and mass effects for particles
        self.is_particle = is_particle

        # Prepare list of X-ray lines to be fitted
        self._define_xray_lines()

        self.fit_background = fit_background

        # Reset detector response function for each spectrum
        DetectorResponseFunction.det_res_conv_matrix = None
        DetectorResponseFunction.icc_conv_matrix = None
        DetectorResponseFunction.setup_detector_response_vars(
            det_ch_offset, det_ch_width, spectrum_lims, microscope_ID, verbose=verbose
        )

        self.print_evolving_params = print_evolving_params
        self.verbose = verbose
        
        
        
    def _define_xray_lines(self):
        """
        Defines the list of elemental X-ray lines to be fitted.
    
        For each element in self.els_to_fit_list, collects all X-ray lines with sufficient overvoltage
        (beam_energy / xray_energy > threshold) and within the spectrum energy range. 
        Also adds escape and pile-up peaks for prominent lines.
    
        Sets:
            self.el_lines_list : list of str
                All X-ray line identifiers to fit (including escape and pile-up peaks).
            self.el_lines_weight_refs_dict : dict
                Maps each X-ray line to its reference line for weight calculation.
        """
        # Minimum overvoltage for peak inclusion, considering spectrum cutoff
        min_overvoltage = max(self.beam_energy / self.energy_vals[-1] * 0.99, 1.2)
        min_energy = self.energy_vals[0]
        max_energy = self.energy_vals[-1]
    
        # Energy threshold for peak inclusion (accounts for detector response broadening)
        peak_en_threshold = min_energy - 3 * DetectorResponseFunction._det_sigma(min_energy)
    
        el_lines_list = []
        el_lines_weight_refs_dict = {}
    
        for el in self.els_to_fit_list:
            el_xRays_dict = get_el_xray_lines(el)  # Get dict of X-ray lines for this element
    
            for xray_line, xray_info in el_xRays_dict.items():
                line_en = xray_info['energy (keV)']
    
                # Only include lines with sufficient overvoltage and within energy bounds
                if self.beam_energy / line_en > min_overvoltage and line_en > peak_en_threshold:
                    xRay_line_str = f"{el}_{xray_line}"
                    el_lines_list.append(xRay_line_str)
                    el_ref_line = self._get_reference_xray_line(el, xray_line, el_xRays_dict)
                    el_lines_weight_refs_dict[xRay_line_str] = el_ref_line
    
                    # Add escape peaks for strong lines above Si K edge
                    # TODO: Replace fixed threshold (0.3) with calibrated value based on counts
                    if xray_info['weight'] > 0.3 and line_en > 1.74:
                        escape_peak_str = xRay_line_str + self.escape_peaks_str
                        el_lines_list.append(escape_peak_str)
                        el_lines_weight_refs_dict[escape_peak_str] = el_ref_line
    
                    # Add pile-up peaks for strong lines within spectrum range
                    # TODO: Replace fixed threshold (0.3) with calibrated value based on counts
                    if xray_info['weight'] > 0.3 and 2 * line_en < max_energy:
                        pileup_peak_str = xRay_line_str + self.pileup_peaks_str
                        el_lines_list.append(pileup_peak_str)
                        el_lines_weight_refs_dict[pileup_peak_str] = el_ref_line
    
        # Reference lines (e.g., Ka1, La1, Ma1) are added first to enforce weight dependencies
        ref_lines = [el_line for el_line in el_lines_list if any(ref in el_line for ref in self.xray_quant_ref_lines)]
        other_lines = list(set(el_lines_list) - set(ref_lines))
    
        # Store the results
        self.el_lines_list = ref_lines + other_lines
        self.el_lines_weight_refs_dict = el_lines_weight_refs_dict
        


    def _get_reference_xray_line(self, el, line, el_xRays_dict):
        """
        Determines the appropriate reference X-ray line for a given characteristic line.
        
        The intensity of the dependent line will be expressed as that of the reference line
        multiplied by a fixed weight, according to NIST conventions.
    
        Parameters
        ----------
        el : str
            Element symbol.
        line : str
            X-ray line identifier (e.g., 'Ka1', 'La1').
        el_xRays_dict : dict
            Dictionary of X-ray lines for the element.
    
        Returns
        -------
        el_ref_line : str
            Reference line string (e.g., 'Fe_Ka1').
    
        Raises
        ------
        RefLineError
            If a suitable reference line cannot be found or if multiple ambiguous references exist.
        """
        # For N lines, use M as the reference group
        if line[0] == 'N':
            ref_line_start = 'M'
        else:
            ref_line_start = line[0]
    
        # Find reference lines in the quantification reference list
        ref_line_l = [ref_line for ref_line in self.xray_quant_ref_lines if ref_line_start == ref_line[0]]
    
        el_line = f"{el}_{line}"
        if len(ref_line_l) == 0:
            raise RefLineError(f"K, L or M references not found for {el_line} line.")
        elif len(ref_line_l) > 1 and ref_line_start in ['K', 'L']:
            raise RefLineError(f"Multiple reference lines found for {el_line}. Only one should be present.")
        elif ref_line_start == 'M':
            # For M lines, select Ma for Z > 58, otherwise Mz (per NIST conventions)
            if len(ref_line_l) > 2:
                raise RefLineError(f"Multiple reference lines found for {el_line}. Only one should be present.")
            else:
                if Element(el).Z > 58:
                    ref_line = [ref_line for ref_line in ref_line_l if ref_line.startswith('Ma')][0]
                else:
                    ref_line = [ref_line for ref_line in ref_line_l if ref_line.startswith('Mz')][0]
        elif ref_line_start == 'K':
            ref_line = ref_line_l[0]
        elif ref_line_start == 'L':
            if 'La1' in el_xRays_dict:
                ref_line = ref_line_l[0]  # Use La1 if available
            else:
                # For 11 < Z < 20, use Ll since La1 is not present according to NIST
                ref_line = 'Ll'
    
        el_ref_line = f"{el}_{ref_line}"
        return el_ref_line



    def _get_fraction_pars(self, elements):
        """
        Create an lmfit Parameters object for elemental mass fractions, which influence background fitting.
    
        Elemental mass fractions affect:
            - The backscattering correction factor (via average atomic number, Z)
            - The mass absorption coefficient (μ), which controls absorption attenuation
            - The generated background model (via parameters P, F, and beta)
    
        The behavior depends on `self.els_w_fr`, which specifies elements with fixed mass fractions:
            - If `els_w_fr` is defined, only those elements are constrained to their given values.
            - All other elements in `elements` are treated as trace components and fitted freely.
            - The total sum of fractions may or may not be constrained to 1, depending on initialization.
    
        Returns
        -------
        fr_params : lmfit.Parameters
            Parameter set representing elemental mass fractions, with values either fixed (from els_w_fr)
            or fitted (for trace elements). Used to compute average Z and μ for background modeling.
        """
        fr_params = Parameters()
    
        if self.els_w_fr:
            # Elements without assigned mass fraction (to be fitted as trace elements)
            trace_els = [el for el in elements if el not in self.els_w_fr]
            
            # Include both trace elements and those with fixed fractions
            if len(trace_els) > 0:
                elements = list(set(trace_els) | set(self.els_w_fr))
            else:
                elements = list(self.els_w_fr)
    
        num_els_to_fit = len(elements)
    
        total_fr = 0
        last_el_fr_expr = '1'  # Used to constrain total fraction to 1 for the last element
    
        for i, el in enumerate(elements):
            par_name = 'f_' + el
            if self.els_w_fr:
                if el in trace_els:
                    if self.verbose:
                        print(f"No fraction was assigned to element {el}. Assuming trace element.")
                    # Distribute remaining mass fraction equally among trace elements
                    val = (1 - sum(self.els_w_fr[el] for el in elements if el in self.els_w_fr)) / len(trace_els)
                else:
                    val = self.els_w_fr[el]
                fr_params.add(par_name, value=val, vary=False)
            else:
                # Optionally constrain the total fraction to 1 by expressing the last element as the remainder
                if self.force_fr_total and i == num_els_to_fit - 1:
                    fr_params.add('f_' + elements[-1], expr=last_el_fr_expr, min=0, max=1)
                else:
                    w_fr = 1 / num_els_to_fit  # Even initial guess
    
                    # Update expression for the last element's fraction
                    last_el_fr_expr += '-' + par_name
                    total_fr += w_fr
    
                    # To enforce sum of fractions <= 1, use cumulative sum parameters (lmfit convention)
                    if i == 0:
                        fr_params.add(par_name, value=w_fr, min=0, max=1)
                        sum_par_name_prev = par_name
                    else:
                        sum_par_name = 'sum' + ''.join(f'_{elements[j]}' for j in range(i + 1))
                        fr_params.add(sum_par_name, value=total_fr, min=0, max=1)
                        fr_params.add(par_name, expr=sum_par_name + '-' + sum_par_name_prev, vary=True, min=0, max=1)
                        sum_par_name_prev = sum_par_name
    
        return fr_params
    
    
    def _initialise_Background_Model(self):
        """
        Instantiate and initialize the background model for the current spectrum.
        
        This function re-initializes stored variables in Background_Model prior to new calculations.
        It is called before every new iteration when quantifying a spectrum iteratively.
    
        Returns
        -------
        bckgrnd_model_and_pars : Background_Model
            Instance of the background model, initialized with current spectrum and fitting parameters.
        """
        # Reinitialize global variables used in background computation
        bckgrnd_model_and_pars = Background_Model(
            self.is_particle,
            self.sp_collection_time,
            self.tot_sp_counts,
            self.beam_energy,
            self.emergence_angle,
            self.els_w_fr,
            self.meas_mode,
            self.energy_vals
        )
        return bckgrnd_model_and_pars
    
    
    def _get_background_mod_pars(self, fitted_elements, fr_pars):
        """
        Generate the background lmfit model and its parameters by calling Background_Model class functions.
    
        The model includes:
            - Generated continuum
            - Backscattered electron correction
            - X-ray absorption attenuation
            - Detector efficiency
            - Detector strobe (zero) peak
    
        Parameters
        ----------
        fitted_elements : list of str
            Elements being fitted in the spectrum.
        fr_pars : lmfit.Parameters
            Parameters object representing elemental fractions.
    
        Returns
        -------
        background_mod : lmfit.Model
            The complete background model.
        background_pars : dict
            Dictionary of background model parameters.
        """
        # Initialize background model for the current spectrum
        bckgrnd_model_and_pars = self._initialise_Background_Model()
        background_mod, background_pars = bckgrnd_model_and_pars.get_full_background_mod_pars(fr_pars)
    
        self.background_mod = background_mod
    
        return background_mod, background_pars
    
    
    def _make_spectrum_mod_pars(self, print_initial_pars=False):
        """
        Generate the peaks lmfit models and parameters by calling Peaks_Model class functions.
    
        Parameters
        ----------
        print_initial_pars : bool, optional
            If True, prints a table of the initial fit parameters and their constraints.
    
        Returns
        -------
        spectrum_mod : lmfit.Model
            Model to be used for fitting the EDS spectrum.
        spectrum_pars : lmfit.Parameters
            Parameters for fitting the EDS spectrum.
        """
        params = Parameters()
    
        # Initialize X-ray peaks model and parameters
        peaks_mod_pars = Peaks_Model(
            spectrum_vals = self.spectrum_vals,
            energy_vals = self.energy_vals,
            microscope_ID = self.microscope_ID,
            meas_mode = self.meas_mode,
            fitting_model = None,
            fitting_pars = params,
            xray_weight_refs_dict=self.el_lines_weight_refs_dict,
            is_particle=self.is_particle,
        )
        
        fitted_peaks = []
        for el_line in self.el_lines_list:
            is_peak_present = peaks_mod_pars._add_peak_model_and_pars(el_line)
            if is_peak_present:
                fitted_peaks.append(el_line)
    
        # Fix centers and sigma of overlapping reference peaks
        peaks_mod_pars._fix_overlapping_ref_peaks()
    
        # Identify elements with peaks present in the spectrum
        fitted_elements = [el for el in self.els_to_fit_list if any(el + '_' in peak for peak in fitted_peaks)]
        fitted_els_to_quantify = [el for el in fitted_elements if el in self.els_to_quantify]
        
        if self.verbose:
            if len(fitted_elements) == 0:
                warnings.warn("No peak from the provided elements was found in the spectrum.", UserWarning)
            elif len(fitted_els_to_quantify) == 0:
                warnings.warn("No peak from the provided elements to quantify was found in the spectrum.", UserWarning)
    
        # Retrieve spectrum model and parameters
        spectrum_mod, spectrum_pars = peaks_mod_pars.get_peaks_mod_pars()
    
        if self.fit_background:
            # Choose elements for which to define fraction parameters (in order of preference)
            if len(fitted_els_to_quantify) > 0:
                fr_pars_els = fitted_els_to_quantify
            elif len(self.els_to_quantify) > 0:
                fr_pars_els = self.els_to_quantify
            elif len(self.els_substrate) > 0:
                fr_pars_els = self.els_substrate
            else:
                raise ValueError(
                    f"No valid element to fit was given to fit the spectrum background. "
                    f"Please provide at least one element (in sample or substrate) that is not {calibs.undetectable_els}, "
                    f"or change the list 'undetectable_els' at calibs.__init__.py"
                )
    
            # Add elemental fraction parameters
            fr_pars = self._get_fraction_pars(fr_pars_els)
            spectrum_pars.update(fr_pars)
    
            # Add background model and parameters
            background_mod, background_pars = self._get_background_mod_pars(fitted_els_to_quantify, fr_pars)
            if spectrum_mod is None:
                spectrum_mod = background_mod
            else:
                spectrum_mod += background_mod
            spectrum_pars.update(background_pars)
        else:
            self.background_mod = None
    
        # Store attributes for later use
        self.spectrum_mod = spectrum_mod
        self.spectrum_pars = spectrum_pars
        self.fitted_els = fitted_elements  # Used by Quantifier class
    
        # Optionally print the initial parameters table for debugging
        if print_initial_pars:
            spectrum_pars.pretty_print()
        
    
    def _iteration_callback(self, params, iter, resid, *args, **kws):
        """
        Callback function to monitor fit iterations during optimization.
    
        Parameters
        ----------
        params : lmfit.Parameters
            Current parameter values.
        iter : int
            Current iteration number.
        resid : np.ndarray
            Residual array at this iteration.
        *args, **kws :
            Additional arguments (ignored).
        """
        self.iteration_counter += 1
    
        # Print progress every 20 iterations if verbose mode is enabled
        if self.verbose and self.iteration_counter % 20 == 0:
            reduced_chi_square = np.sum(resid**2) / (len(self.energy_vals) - len(self.spectrum_pars))
            print(f"Iter. #: {self.iteration_counter}. Residual sum of squares: {reduced_chi_square:.5e}")
    
        # Print evolving parameter values for debugging, if enabled
        if self.print_evolving_params:
            print_single_separator()
            print(f"Params changed in iteration #{self.iteration_counter}")
            if self.iteration_counter == 1:
                # Store initial parameter values for comparison
                self.param_values = {param: params[param].value for param in params}
            else:
                for param in params:
                    par_value = params[param].value
                    # Print only parameters that have changed and are being varied
                    if par_value != self.param_values[param] and params[param].vary:
                        print(f"{param}: {par_value}")
                        self.param_values[param] = par_value
    
                        # Check for NaN values in background or spectrum if background parameter changes
                        if param == 'rhoz_par_slope':
                            bckngrd_contains_nan = np.any(np.isnan(self.background_mod.eval(params=params, x=self.energy_vals)))
                            print("Background contains nan vals: ", bckngrd_contains_nan)
                            sp_contains_nan = np.any(np.isnan(self.spectrum_mod.eval(params=params, x=self.energy_vals)))
                            print("Spectrum contains nan vals: ", sp_contains_nan)
                            # Uncomment for plotting if needed during development
                            # plt.figure()
                            # plt.plot(self.energy_vals, self.background_mod.eval(params=params, x=self.energy_vals))
                        
                        
        
    def fit_spectrum(self, parameters=None, initial_par_vals=None, function_tolerance=1e-3,
                     plot_result=False, print_result=False, print_result_extended=False, n_iter=None):
        """
        Fit the EDS spectrum using lmfit.
    
        Parameters
        ----------
        parameters : lmfit.Parameters, optional
            Parameters object to use for fitting. If None, parameters are generated internally.
        initial_par_vals : dict, optional
            Dictionary of initial parameter values to override defaults.
        function_tolerance : float, optional
            ftol used in scipy.optimize.leastsq (default: 1e-3).
        plot_result : bool, optional
            Whether to plot the fitted spectrum (total fit and background).
        print_result : bool, optional
            Whether to print the quality of the fit.
        print_result_extended : bool, optional
            Whether to print extended fit results.
        n_iter : int, optional
            Iteration number (for display purposes).
    
        Returns
        -------
        fit_result : lmfit.ModelResult
            Contains all information about the result of the fit.
        fitted_lines : list of str
            List of fitted X-ray lines.
        """
        # Initialize iteration counter for callback tracking
        self.iteration_counter = 0
    
        # Build or assign spectrum model and parameters
        if parameters is None:
            self._make_spectrum_mod_pars()
        else:
            self.spectrum_pars = parameters
            if self.fit_background:
                # Re-initialize background model if fitting iteratively
                self._initialise_Background_Model()
    
        params = self.spectrum_pars
    
        # Display fitting progress if verbose
        if self.verbose:
            print_double_separator()
            print_double_separator()
            if n_iter:
                print(f"Iteration #{n_iter}")
            print_single_separator()
            print('Fitting spectrum...')
            start_time = time.time()
    
        # Set user-specified initial parameter values
        if initial_par_vals:
            for par, val in initial_par_vals.items():
                self.spectrum_pars[par].value = val
    
        # Setting full_output=False prevents calculation of uncertainties for faster fitting, but throws error in original versiopn of lmfit
        try:
            fit_result = self.spectrum_mod.fit(
                self.spectrum_vals,
                params,
                x=self.energy_vals,
                iter_cb=self._iteration_callback,
                verbose=True,
                fit_kws={'ftol': function_tolerance, 'full_output': False}
            )
        except:
            # main lmfit branch does not accept 'full_output': False
            fit_result = self.spectrum_mod.fit(
                self.spectrum_vals,
                params,
                x=self.energy_vals,
                iter_cb=self._iteration_callback,
                verbose=True,
                fit_kws={'ftol': function_tolerance}
            )
    
        if self.verbose:
            fitting_time = time.time() - start_time
            print(f'Fit completed in {fitting_time:.1f} s with {self.iteration_counter} steps')
    
        # Identify fitted X-ray lines
        used_params = list(fit_result.params.keys())
        fitted_lines = ['_'.join(param.split('_')[:-1]) for param in used_params if "area" in param]
    
        self.fit_result = fit_result
    
        # Plot fit result if requested and running interactively
        if plot_result and self.verbose:
            # self.verbose is always False when called by SEMEDS_analyser (prevents plotting in batch)
            self.plot_result()
    
        # Print fit results if requested or in verbose mode
        if print_result or self.verbose:
            self.print_result(print_only_independent_params=False, extended=print_result_extended)
    
        return fit_result, fitted_lines
    
     
    def plot_result(self):
        """
        Plot the fitted EDS spectrum and its individual background components.
    
        Displays the total fit, background components, and a residual plot.
        """
        fig = self.fit_result.plot(xlabel='Energy (keV)', ylabel='Counts')
    
        if self.fit_background:
            # Retrieve individual background components from the fit
            components = self.fit_result.eval_components(x=self.energy_vals)
            abs_att_param_name = [s for s in components.keys() if '_abs_att' in s][0]
            gen_bckgrnd_param_name = [s for s in components.keys() if '_generated_bckgrnd' in s][0]
            det_eff_par_name = '_det_efficiency'
            bcksctr_corr_par_name = '_backscattering_correction'
            stop_power_par_name = '_stopping_power'
            det_zero_peak_par_name = 'det_zero_peak_'
    
            # Extract and compute background components
            gen_background = components[gen_bckgrnd_param_name]
            abs_attenuation = components[abs_att_param_name]
            det_efficiency = components[det_eff_par_name]
            bcksctr_corr = components[bcksctr_corr_par_name]
            stopping_power = components[stop_power_par_name]
            det_zero_peak = components[det_zero_peak_par_name]
            total_background_component = gen_background * abs_attenuation * det_efficiency * bcksctr_corr
    
            # Plot each background component
            plt.plot(self.energy_vals, gen_background, 'y--', label='Generated Background')
            plt.plot(self.energy_vals, abs_attenuation * 100, 'r--', label='Absorption (x100)')
            plt.plot(self.energy_vals, det_efficiency * 100, 'b--', label='Detector efficiency (x100)')
            plt.plot(self.energy_vals, stopping_power * 100, 'c--', label='Stopping power (x100)')
            plt.plot(self.energy_vals, bcksctr_corr * 100, 'm--', label='Backscatter corr. (x100)')
            plt.plot(self.energy_vals, det_zero_peak, 'y--', label='Detector zero peak')
            plt.plot(self.energy_vals, total_background_component, 'k--', label='Tot. background')
    
            axes = fig.get_axes()
            axes[0].set_title('Residual plot')
    
        plt.legend()
        plt.show()

            
        
    def print_result(self, extended=False, print_only_independent_params=False):
        """
        Print a summary of the fit results.
        
        Parameters
        ----------
        extended : bool, optional
            If True, prints the full fit report or parameter table. Otherwise, prints only summary statistics
            (Reduced chi-squared and R-squared).
        print_only_independent_params : bool, optional
            If True and extended is True, prints only parameter names and their values.
        """
        print_double_separator()
        if extended:
            if print_only_independent_params:
                print('Parameter           Value')
                for name, param in self.fit_result.params.items():
                    print(f'{name:20s} {param.value:11.10f}')
            else:
                # Print the full fit report from lmfit
                print(self.fit_result.fit_report())
        else:
            reduced_chi_square = self.fit_result.redchi
            r_squared = 1 - self.fit_result.residual.var() / np.var(self.spectrum_vals)
            print('Fit results:')
            print(f"Reduced Chi-square: {reduced_chi_square:.2f}")
            print(f"R-squared: {r_squared:.5f}")


#%% Peaks_Model class
class Peaks_Model:
    """
    Model for X-ray spectral peaks in EDS spectra.

    Handles the construction, parameterization, and constraints of peak models for EDS spectral fitting,
    including area weighting, escape/pileup peaks, and peak shape calibration.

    Attributes
    ----------
    spectrum_vals : array-like or None
        Measured spectrum values (e.g., counts).
    energy_vals : array-like or None
        Corresponding energy values (e.g., keV).
    fitting_model : lmfit.Model or None
        Composite model used for fitting.
    fitting_params : lmfit.Parameters or None
        Parameters for the current model.
    xray_weight_refs_dict : dict
        Maps each secondary line to its reference line for area weighting.
    xray_weight_refs_lines : list
        Unique reference lines extracted from xray_weight_refs_dict.
    microscope_ID : str
        Identifier for the microscope/calibration to use.
    meas_mode : str
        EDS mode, e.g., 'point', 'map', etc.
    is_particle : bool
        Whether the sample is a particle (affects some constraints).
    free_area_el_lines : list
        Elements/lines whose area is fitted freely (for peak intensity weight calibration).
    free_peak_shapes_els : list
        Elements whose peak shapes are calibrated (for shape calibration).
    fixed_peaks_dict : dict
        Tracks dependencies for overlapping peaks.
    
    Class attributes
    ----------------
    icc_freq_spectra : dict (class variable)
        Cache for ICC convolution spectra, shared across all instances.
    center_key, sigma_key, area_key, center_offset_key, sigma_broadening_key,
    gamma_key, tail_fraction_key, F_loss_key, R_e_key, pileup_peaks_str,
    pileup_peaks_str: str
        Standardized parameter names for model building.
        
    Notes
    -----
    - `icc_freq_spectra` is a class variable, shared across all instances and used for static method access.
    - Requires XSp_Fitter to be initialized first for correct DetectorResponseFunction loading.
    """

    # Class-level cache for ICC convolution spectra.
    icc_freq_spectra = {}

    # Standardized parameter keys for all peaks
    center_key = 'center'
    sigma_key = 'sigma'
    area_key = 'area'
    center_offset_key = 'cen_offset'
    sigma_broadening_key = 'sigma_broad'
    gamma_key = 'gamma'
    tail_fraction_key = 'f_tail'
    F_loss_key = 'F_loss'
    R_e_key = 'R_e'
    escape_peaks_str = XSp_Fitter.escape_peaks_str
    pileup_peaks_str = XSp_Fitter.pileup_peaks_str

    def __init__(
        self,
        spectrum_vals,
        energy_vals,
        microscope_ID,
        meas_mode,
        fitting_model,
        fitting_pars,
        xray_weight_refs_dict=None,
        is_particle=False,
        free_area_el_lines=None,
        free_peak_shapes_els=None
    ):
        """
        Initialize a Peaks_Model instance for X-ray spectral peak modeling.
    
        Parameters
        ----------
        spectrum_vals : array-like
            Measured spectrum values (e.g., counts).
        energy_vals : array-like
            Corresponding energy values (e.g., keV) for the spectrum.
        microscope_ID : str
            Identifier for the microscope/calibration to use.
        meas_mode : str
            EDS mode, e.g., 'point', 'map', etc.
        fitting_model : lmfit.Model
            Composite model for all peaks (should be built prior to fitting).
        fitting_pars : lmfit.Parameters
            Parameters for the current model (should be built prior to fitting).
        xray_weight_refs_dict : dict or None, optional
            Dictionary mapping each secondary line to its reference line for area weighting.
            Example: {'Fe_Kb1': 'Fe_Ka1'}.
        is_particle : bool, optional
            Whether the sample is a particle (affects some constraints).
        free_area_el_lines : list or None, optional
            Elements/lines whose area is fitted freely (for peak intensity weight calibration).
        free_peak_shapes_els : list or None, optional
            Elements whose peak shapes are calibrated (for peak shape calibration).

        Notes
        -----
        - If `xray_weight_refs_dict` is not provided, it defaults to an empty dictionary.
        - If `free_area_el_lines` is not provided, it defaults to ['Ge_Lb1'].
        - If `free_peak_shapes_els` is not provided, it defaults to an empty list.
        - `icc_freq_spectra` is reset for each new instance.
        """
        # Store spectrum and energy values
        self.spectrum_vals = spectrum_vals
        self.energy_vals = energy_vals

        # Model and its parameters (set during fitting/building)
        self.fitting_model = fitting_model
        self.fitting_params = fitting_pars

        # Sample/experiment settings
        self.is_particle = is_particle
        self.meas_mode = meas_mode
        calibs.load_microscope_calibrations(microscope_ID, meas_mode, load_detector_channel_params=False)

        # X-ray line references for area weighting
        if xray_weight_refs_dict is None:
            xray_weight_refs_dict = {}
        self.xray_weight_refs_dict = xray_weight_refs_dict
        self.xray_weight_refs_lines = list(set(self.xray_weight_refs_dict.values()))

        # Elements/lines whose area is fitted freely (for weight calibration)
        if free_area_el_lines is None:
            free_area_el_lines = ['Ge_Lb1']
        self.free_area_el_lines = free_area_el_lines

        # Elements whose peak shapes are calibrated (for shape calibration)
        if free_peak_shapes_els is None:
            free_peak_shapes_els = []
        self.free_peak_shapes_els = free_peak_shapes_els

        # Reset icc spectra cache
        self.clear_cached_icc_spectra()

        # Bookkeeping for overlapping/fixed peaks (populated elsewhere)
        self.fixed_peaks_dict = {}
    
    @staticmethod
    def clear_cached_icc_spectra():
        """
        Clear the class-level cache for ICC convolution spectra.
    
        This method resets the `icc_freq_spectra` dictionary shared by all instances
        of Peaks_Model, ensuring that all cached ICC spectra are removed prior new
        calculations that require new computations of icc_freq_spectra.
        """
        Peaks_Model.icc_freq_spectra = {}
    
    
    def get_peaks_mod_pars(self):
        """
        Return the current composite peak model and its associated fit parameters.
    
        Returns
        -------
        model : lmfit.Model or None
            The composite model representing all fitted peaks, or None if not yet defined.
        pars : lmfit.Parameters
            The set of parameters for the current model.
        """
        return self.fitting_model, self.fitting_params
    
    
    # =============================================================================
    # Peak models
    # ============================================================================= 
    @staticmethod
    def _gaussian(x, area, center, sigma):
        """
        Return a Gaussian function.
    
        Parameters
        ----------
        x : array-like
            Independent variable (energy).
        area : float
            Total area under the peak.
        center : float
            Peak center (mean).
        sigma : float
            Standard deviation (width).
    
        Returns
        -------
        model : array-like
            Gaussian function evaluated at x.
        """
        return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x-center)**2 / (2*sigma**2))
    
    
    @staticmethod
    def _gaussian_tail(x, area, center, sigma, gamma):
        """
        Return a skewed tail function for low-Z X-ray peaks.
    
        Reference
        ---------
        J. Osán et al., "Evaluation of energy‐dispersive x‐ray spectra of low‐Z elements from electron‐probe microanalysis of individual particles,"
        X-Ray Spectrom. 30 (2001) 419–426. https://doi.org/10.1002/xrs.523
    
        Parameters
        ----------
        x : array-like
            Independent variable (energy).
        area : float
            Total area under the peak.
        center : float
            Peak center (mean).
        sigma : float
            Standard deviation (width).
        gamma : float
            Tail parameter.
    
        Returns
        -------
        model : array-like
            Skewed tail function evaluated at x.
        """
        return area / (2 * gamma * sigma * np.exp(-1 / (2 * gamma**2))) \
            * np.exp((x - center) / (gamma * sigma)) \
            * erfc((x - center) / (np.sqrt(2) * sigma) + 1 / (np.sqrt(2) * gamma))
    
    
    @staticmethod
    def _gaussian_with_tail(x, area, center, sigma, gamma, f_tail):
        """
        Return a Gaussian with a low-energy tail.
    
        Reference
        ---------
        J. Osán et al., X-Ray Spectrom. 30 (2001) 419–426. https://doi.org/10.1002/xrs.523
    
        Parameters
        ----------
        x : array-like
            Independent variable (energy).
        area : float
            Total area under the peak.
        center : float
            Peak center (mean).
        sigma : float
            Standard deviation (width).
        gamma : float
            Tail parameter.
        f_tail : float
            Fractional intensity of the tail.
    
        Returns
        -------
        model : array-like
            Gaussian with tail evaluated at x.
        """
        return Peaks_Model._gaussian(x, area, center, sigma) + f_tail * Peaks_Model._gaussian_tail(x, area, center, sigma, gamma)
    
    
    @staticmethod
    def _gaussian_shelf(x, area, center, sigma):
        """
        Return a shelf function for X-ray peaks.
        
        Note: This function is not currently used in the model pipeline.
    
        Reference
        ---------
        J. Osán et al., X-Ray Spectrom. 30 (2001) 419–426. https://doi.org/10.1002/xrs.523
    
        Parameters
        ----------
        x : array-like
            Independent variable (energy).
        area : float
            Total area under the peak.
        center : float
            Peak center (mean).
        sigma : float
            Standard deviation (width).
    
        Returns
        -------
        model : array-like
            Shelf function evaluated at x.
        """
        return area / (2 * center) * erfc((x - center) / (np.sqrt(2) * sigma))
    
    
    @staticmethod
    def _gaussian_with_tail_and_shelf(x, area, center, sigma, gamma, f_tail, f_shelf):
        """
        Return a Gaussian with both a low-energy tail and a shelf component.
        
        Note: This function is not currently used in the model pipeline.
    
        Reference
        ---------
        J. Osán et al., X-Ray Spectrom. 30 (2001) 419–426. https://doi.org/10.1002/xrs.523
    
        Parameters
        ----------
        x : array-like
            Independent variable (energy).
        area : float
            Total area under the peak.
        center : float
            Peak center (mean).
        sigma : float
            Standard deviation (width).
        gamma : float
            Tail parameter.
        f_tail : float
            Fractional intensity of the tail.
        f_shelf : float
            Fractional intensity of the shelf.
    
        Returns
        -------
        model : array-like
            Gaussian with tail and shelf evaluated at x.
        """
        return (Peaks_Model._gaussian(x, area, center, sigma)
                + f_tail * Peaks_Model._gaussian_tail(x, area, center, sigma, gamma)
                + f_shelf * Peaks_Model._gaussian_shelf(x, area, center, sigma))
    
    
    @staticmethod
    def _gaussian_with_tail_and_icc(x, area, center, sigma, R_e, F_loss, gamma, f_tail):
        """
        Compute a Gaussian peak with a low-energy tail, convolved with the incomplete charge collection (ICC) model.
    
        ICC model as described in:
        Redus, R. H., & Huber, A. C. (2015). Response Function of Silicon Drift Detectors for Low Energy X-rays.
        In Advances in X-ray Analysis (AXA) (pp. 274–282). International Centre for Diffraction Data (ICDD).
    
        Parameters
        ----------
        x : array-like
            Energy axis.
        area : float
            Total area under the peak.
        center : float
            Peak center (mean, in keV).
        sigma : float
            Standard deviation (width).
        R_e : float
            Effective recombination parameter for ICC.
        F_loss : float
            Fractional charge loss parameter for ICC.
        gamma : float
            Tail parameter for the skewed Gaussian.
        f_tail : float
            Fractional intensity of the tail.
    
        Returns
        -------
        g_with_icc : array-like
            Gaussian with tail, convolved with the ICC distribution.
        """
        # Use high precision for keys to ensure accurate caching of ICC distributions and avoid redundant expensive computations
        key_center = f"{center:.6f}"      # Precision of 0.1 eV
        key_F_loss = f"{F_loss:.6f}"
        key_R_e = f"{R_e * 1e7:.5f}"      # Precision of 1 nm
        key_f_tail = f"{f_tail:.6f}"
        key_gamma = f"{gamma:.6f}"
        key = f"{key_center}_{key_F_loss}_{key_R_e}_{key_f_tail}_{key_gamma}"
    
        # Use the class variable for ICC cache
        icc_cache = Peaks_Model.icc_freq_spectra
        icc_n_vals_distr = icc_cache.get(key)
        if icc_n_vals_distr is None:
            icc_n_vals_distr = DetectorResponseFunction.get_icc_spectrum(x, center, R_e, F_loss)
            icc_cache[key] = icc_n_vals_distr
    
        # Compute Gaussian with tail
        g_vals = Peaks_Model._gaussian_with_tail(x, area, center, sigma, gamma, f_tail)
    
        # Convolve skewed Gaussian with ICC distribution
        g_with_icc = np.convolve(np.array(g_vals), np.array(icc_n_vals_distr), mode='same')
    
        return g_with_icc
    
    
    @staticmethod
    def _gaussian_with_icc(x, area, center, sigma, R_e, F_loss):
        """
        Compute a Gaussian peak convolved with the incomplete charge collection (ICC) model.
    
        Based on:
        Redus, R. H., & Huber, A. C. (2015). Response Function of Silicon Drift Detectors for Low Energy X-rays.
        In Advances in X-ray Analysis (AXA) (pp. 274–282). International Centre for Diffraction Data (ICDD).
    
        Parameters
        ----------
        x : array-like
            Energy axis.
        area : float
            Total area under the peak.
        center : float
            Peak center (mean, in keV).
        sigma : float
            Standard deviation (width).
        R_e : float
            Effective recombination parameter for ICC.
        F_loss : float
            Fractional charge loss parameter for ICC.
    
        Returns
        -------
        g_with_icc : array-like
            Gaussian convolved with the ICC distribution.
        """
        # Use high precision for keys to ensure accurate caching of ICC distributions and avoid redundant expensive computations
        key_center = f"{center:.6f}"      # Precision of 0.1 eV
        key_F_loss = f"{F_loss:.6f}"
        key_R_e = f"{R_e * 1e7:.5f}"      # Precision of 1 nm
        key = f"{key_center}_{key_F_loss}_{key_R_e}"
    
        # Use the class variable for ICC cache
        icc_cache = Peaks_Model.icc_freq_spectra
        icc_n_vals_distr = icc_cache.get(key)
        if icc_n_vals_distr is None:
            icc_n_vals_distr = DetectorResponseFunction.get_icc_spectrum(x, center, R_e, F_loss)
            icc_cache[key] = icc_n_vals_distr
    
        # Compute Gaussian
        g_vals = Peaks_Model._gaussian(x, area, center, sigma)
    
        # Convolve Gaussian with ICC distribution
        g_with_icc = np.convolve(np.array(g_vals), np.array(icc_n_vals_distr), mode='same')
    
        return g_with_icc
    
    # =============================================================================
    # Functions for spectral data analysis: initial parameter estimation and peak identification
    # =============================================================================
    def _identify_peak(self, line_energy, sigma_peak, el_line):
        """
        Identify and characterize the presence of a peak in the spectrum for a given characteristic X-ray line.
        
        Determines whether the specified peak is present, if it is small (to be modeled without a tail),
        or if it is likely overlapping with a larger peak (in which case a fixed initial area is assigned).
        
        Parameters
        ----------
        line_energy : float
            Nominal energy of the characteristic X-ray line (in keV).
        sigma_peak : float
            Estimated standard deviation (width) of the peak.
        el_line : str
            Element and line identifier (e.g., 'Fe_Ka1').
        
        Returns
        -------
        peak_indices : np.ndarray
            Indices of points within ±3σ of the peak energy.
        peak_height : list of float
            Heights of detected peaks.
        is_small_peak : bool or None
            True if the peak is present but small, False if not small, None if not found.
        is_peak_overlapping : bool
            True if the peak is likely overlapping with a larger peak, False otherwise.
        
        Potential improvements:
        ----
        - For improved reliability, peak identification should ideally be performed for the entire group of
            characteristic X-ray lines of an element. This approach helps correctly identify peaks in cases of potential overlap.
            For example, La1 line should not be included if Ka1 line is not present in the spectrum for that element
            (considering appropriate intensity ratios).
        - Peak detection shoulkd account for shift of low-energy peaks (<0.5 keV) due to detector incomplete charge collection
        """
        channel_width_keV = self.energy_vals[1] - self.energy_vals[0]

        # Indices of points within ±3σ of the peak energy
        peak_indices = np.where(
            (self.energy_vals > line_energy - sigma_peak * 3) &
            (self.energy_vals < line_energy + sigma_peak * 3)
        )[0]
    
        if len(peak_indices) == 0:
            return peak_indices, [], None, False
    
        peak_sp_vals = self.spectrum_vals[peak_indices]
    
        # Smooth the curve to reduce noise
        filter_len = 2
        peak_sp_vals_blurred = np.convolve(peak_sp_vals, np.ones(filter_len) / filter_len, mode='valid')
    
        # Minimum prominence for peak detection (10% of local baseline, or 4 counts)
        min_prominence = np.mean(np.array([*peak_sp_vals_blurred[:2], *peak_sp_vals_blurred[-2:]])) / 10
        min_prominence = max(min_prominence, 4)
    
        # Estimate number of points constituting the peak (FWHM increases at higher energy)
        peak_len = int(sigma_peak / channel_width_keV * 6) # Uses 6 sigmas as peak width to be conservative
    
        # Find peaks in the blurred spectrum
        peaks, _ = find_peaks(
            peak_sp_vals_blurred,
            prominence=min_prominence,
            wlen=peak_len,
            distance=5
        )
        peak_pos = peak_indices[peaks]
        peak_height = list(self.spectrum_vals[peak_pos])
    
        # Assess peak prominence to determine if the peak is small (less than 50 counts prominence)
        peak_prom, _, _ = peak_prominences(peak_sp_vals_blurred, peaks, wlen=peak_len)
        if peak_prom.size != 0:
            is_small_peak = peak_prom[0] < 50
        else:
            is_small_peak = None
        
        if peak_height:
            is_peak_overlapping = False
        else:
            # Fit a line to the region to check for overlapping
            peak_energy_vals = self.energy_vals[peak_indices]
            slope = np.polyfit(peak_energy_vals, peak_sp_vals, deg=1)[0]
            # Make overlap threshold scale with channel width (100 for 10 eV, proportional otherwise)
            overlap_threshold = 100 * (channel_width_keV / 0.01)
            if abs(slope) < overlap_threshold:
                is_peak_overlapping = False
            else:
                is_peak_overlapping = True
    
        return peak_indices, peak_height, is_small_peak, is_peak_overlapping
    
    
    def _estimate_peak_area(self, line_energy, sigma_peak, el_line):
        """
        Estimate the area of a characteristic X-ray peak in the spectrum.
    
        Parameters
        ----------
        line_energy : float
            Nominal energy of the characteristic X-ray line (in keV).
        sigma_peak : float
            Estimated standard deviation (width) of the peak.
        el_line : str
            Element and line identifier (e.g., 'Fe_Ka1').
    
        Returns
        -------
        peak_area : float
            Estimated area of the peak.
        is_small_peak : bool or None
            True if the peak is present but small, False if not small, None if not found.
        is_peak_overlapping : bool
            True if the peak is likely overlapping with a larger peak, False otherwise.
        """
        peak_indices, peak_height, is_small_peak, is_peak_overlapping = self._identify_peak(
            line_energy, sigma_peak, el_line
        )
    
        if len(peak_indices) == 0:
            peak_area = 0
            return peak_area, is_small_peak, is_peak_overlapping
    
        if peak_height:
            peak_area = peak_height[0] * (sigma_peak * np.sqrt(2 * np.pi)) * 0.7
            # factor that reduces overestimation of peak area, often due to central peak point being higher than the Gaussian curve
        else:
            if is_peak_overlapping:
                peak_area = 10
            else:
                peak_area = 0
    
        return peak_area, is_small_peak, is_peak_overlapping
    
    # =============================================================================
    # Create lmfit peak model and parameters
    # =============================================================================
    def _get_gaussian_center_param(self, line_prefix, line_energy):
        """
        Create a lmfit Parameter for the center (energy) of a Gaussian X-ray peak.
    
        For low-energy peaks (<0.6 keV), the center is allowed to vary more due to incomplete charge collection (ICC),
        which can result in greater peak shifts. For higher energy peaks, the allowed variation is much smaller.
    
        Parameters
        ----------
        line_prefix : str
            Prefix for the parameter name (typically includes element and line identifier).
        line_energy : float
            Nominal energy of the X-ray line (in keV).
    
        Returns
        -------
        center_par : lmfit.Parameter
            Parameter object for the Gaussian peak center, with appropriate bounds.
        """
        if line_energy < 0.6:
            # Low-energy peaks: allow wider movement due to ICC effects
            center_par = Parameter(
                line_prefix + self.center_key,
                value=line_energy,
                vary=True,
                min=line_energy * 0.93,
                max=line_energy * 1.03
            )
        else:
            # High-energy peaks: only minor shifts expected
            center_par = Parameter(
                line_prefix + self.center_key,
                value=line_energy,
                vary=True,
                min=line_energy * 0.99,
                max=line_energy * 1.01
            )
        return center_par


    def _get_gaussian_sigma_param(self, line_prefix, sigma_init):
        """
        Create an lmfit Parameter for the sigma (width) of a Gaussian X-ray peak.
    
        The parameter is allowed to vary between 93% and 110% of the initial sigma value.
    
        Parameters
        ----------
        line_prefix : str
            Prefix for the parameter name (typically includes element and line identifier).
        sigma_init : float
            Initial value for the peak width (sigma).
    
        Returns
        -------
        sigma_par : lmfit.Parameter
            Parameter object for the Gaussian peak width, with appropriate bounds.
        """
        return Parameter(
            line_prefix + self.sigma_key,
            value=sigma_init,
            vary=True,
            min=sigma_init * 0.93,
            max=sigma_init * 1.1
        )


    def _add_peak_model_and_pars(self, el_line):
        """
        Add a peak model and its parameters for a specific X-ray line to the composite model.
        
        Handles normal, escape, and pileup peaks, with proper parameter dependencies and constraints.
        Updates self.fitting_model and self.fitting_params in-place.
        
        Parameters
        ----------
        el_line : str
            X-ray line identifier, e.g. 'Fe_Ka1', 'Si_Ka1_escape', etc.
        
        Returns
        -------
        is_peak_present : bool
            True if the peak is present and added to the model, False otherwise.
        
        Notes
        -----
        This function is designed for maintainability and clarity, with explicit
        comments at each decision point. It is intended for use in EDS spectrum
        modeling where peak shapes, dependencies, and constraints must be carefully
        managed.
        
        References
        ----------
        - NIST X-ray Transition Energies: https://physics.nist.gov/PhysRefData/XrayTrans/Html/
        - Osan, J. et al., X-ray Spectrom. 2001, 30, 8–17.
        - Redus, R. H. et al., X-ray Spectrom. 2015, 44, 283–291.
        """
        
        # Get current value of model and parameters
        model = self.fitting_model
        params = self.fitting_params
        
         # Extract element and X-ray line (e.g., 'Fe', 'Ka1')
        el, line = el_line.split('_', 1)
        line_prefix = el_line + '_'  # String added to parameter names for uniqueness
        
        # --- Identify if escape or pileup peak, and its reference line ---
        is_escape_peak = self.escape_peaks_str in line
        is_pileup_peak = self.pileup_peaks_str in line
    
        escape_ref_line = None
        is_escape_ref_peak = False
        if is_escape_peak:
            # Reference peak for escape: remove escape suffix
            escape_ref_el_line = el_line[:-len(self.escape_peaks_str)]
            escape_ref_line = escape_ref_el_line.split('_')[1]  # e.g. 'Ka1'
            is_escape_ref_peak = escape_ref_el_line in self.xray_weight_refs_lines
    
        pileup_ref_line = None
        is_pileup_ref_peak = False
        if is_pileup_peak:
            # Reference peak for pileup: remove pileup suffix
            pileup_ref_el_line = el_line[:-len(self.pileup_peaks_str)]
            pileup_ref_line = pileup_ref_el_line.split('_')[1]
            is_pileup_ref_peak = pileup_ref_el_line in self.xray_weight_refs_lines
            
        # --- Get theoretical X-ray energy and weight ---
        lines = get_el_xray_lines(el)  # Used for line energy and weight
        if is_escape_peak:
            # Escape peak energy: reference line minus Si Ka escape energy
            line_energy = lines[escape_ref_line]['energy (keV)'] - 1.740
        elif is_pileup_peak:
            # Pileup peak energy: twice the reference line energy
            line_energy = lines[pileup_ref_line]['energy (keV)'] * 2
        else:
            # Normal peak energy
            line_energy = lines[line]['energy (keV)']
        
        # --- Set up initial parameters for this peak ---
        center_par = self._get_gaussian_center_param(line_prefix, line_energy)
        sigma_init = DetectorResponseFunction._det_sigma(line_energy)
        sigma_par = self._get_gaussian_sigma_param(line_prefix, sigma_init)
        
        # Estimate initial peak area and properties
        # If the peak is small, it will be treated as a standard gaussian instead of adding the exponential tail
        initial_area, is_small_peak, is_peak_overlapping = self._estimate_peak_area(line_energy, sigma_init, el_line)
        
        # --- Handle dependent peaks (area is weighted on a reference peak) ---
        if el_line not in self.xray_weight_refs_lines:
            # Get reference line for area weighting
            ref_line = self.xray_weight_refs_dict[el_line]
            ref_line_prefix = ref_line + '_'
            ref_peak_area_param = params.get(ref_line_prefix + self.area_key)
            
            # If reference peak is not present, skip this peak
            if ref_peak_area_param is None or ref_peak_area_param.value == 0:
                peak_m = None
            else:
                if is_escape_peak:
                    # Escape peaks: modeled with Gaussian; tail is negligible for small peaks
                    peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
                    if is_escape_ref_peak:
                        # Fit area for escape reference peaks, constrain by escape probability
                        max_area_escape = calibs.escape_peak_probability[self.meas_mode] * ref_peak_area_param.value
                        params.add(line_prefix + self.area_key, value=initial_area, min=0, max=max_area_escape)
                    else:
                        # Dependent escape peaks: area is weighted relative to reference escape peak
                        weight = lines[escape_ref_line]['weight']
                        escape_ref_line_prefix = ref_line + self.escape_peaks_str + '_'
                        if params.get(escape_ref_line_prefix + self.area_key) is not None:
                            # Reference escape peak exists
                            params.add(line_prefix + self.area_key, expr=escape_ref_line_prefix + self.area_key + f'*{weight}')
                        else:
                            # If reference escape peak is not in spectrum, fit area with constraint
                            max_area_escape = weight * calibs.escape_peak_probability[self.meas_mode] * ref_peak_area_param.value
                            params.add(line_prefix + self.area_key, value=initial_area, min=0, max=max_area_escape)
                elif is_pileup_peak:
                    # Pileup peaks: modeled with Gaussian for simplicity
                    peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
                    if is_pileup_ref_peak:
                        # Fit area for reference pileup peaks, constrain by pileup probability
                        max_area_pileup = calibs.pileup_peak_probability[self.meas_mode] * ref_peak_area_param.value
                        vary_area = max_area_pileup > 0.1
                        params.add(line_prefix + self.area_key, value=initial_area, vary=vary_area, min=0, max=max_area_pileup)
                    else:
                        # Dependent pileup peaks: area is weighted relative to reference pileup peak
                        weight = lines[pileup_ref_line]['weight']
                        pileup_ref_line_prefix = ref_line + self.pileup_peaks_str + '_'
                        if pileup_ref_line_prefix + self.area_key in params.keys():
                            params.add(line_prefix + self.area_key, expr=pileup_ref_line_prefix + self.area_key + f'*{weight}')
                        else:
                            # Reference pileup peak not present
                            peak_m = None
                else:
                    # Regular dependent peak (not escape or pileup)
                    # During peak shape calibration, fix shape parameters to match reference peak
                    calibrate_peak_shape_params = el in self.free_peak_shapes_els
                    peak_m, _ = self._get_peak_model_and_update_pars(
                        params, line_energy, line_prefix,
                        is_calibration=calibrate_peak_shape_params,
                        ref_line_prefix=ref_line_prefix
                    )
                    if el_line in self.free_area_el_lines:
                        # Allow area to vary for calibration of peak weights
                        params.add(line_prefix + self.area_key, value=initial_area, min=0, vary=True)
                    else:
                        # Area is weighted relative to reference peak
                        weight = lines[line]['weight']
                        params.add(line_prefix + self.area_key, expr=ref_line_prefix + self.area_key + f'*{weight}')
    
                if peak_m is not None:
                    # Fix all dimensions of these dependent peaks
                    params.add(line_prefix + self.center_key, expr=f"{line_energy} - {ref_line_prefix}{self.center_offset_key}")
                    params.add(line_prefix + self.sigma_key, expr=f"{sigma_init} * {ref_line_prefix}{self.sigma_broadening_key}")
        
        # --- Handle reference peaks (area fitted directly) ---
        else: # Reference lines with weight = 1, as per tabualted NIST values
            if initial_area == 0:
                # Peak is absent: add area parameter so that dependent peaks can also be exluded from the fit, but do not add model
                params.add(line_prefix + self.area_key, value=initial_area, min=0, vary=False)
                peak_m = None
            else:
                # Peak detected or overlapping with another peak
                vary_area = True
                if line == 'Ll':
                    # Special handling for Ll lines missing their reference La1 lines (elements 12 < Z < 20)
                    ref_line_prefix = el + '_Ka1_'
                    ref_peak_area_param = params.get(ref_line_prefix + self.area_key)
                    if ref_peak_area_param is None or ref_peak_area_param.value == 0:
                        initial_area = 0
                        max_area = 1  # Placeholder
                        vary_area = False
                    elif not self.is_particle:
                        # Area calibrated in bulk standards
                        max_area = calibs.weight_Ll_ref_Ka1[self.meas_mode] * ref_peak_area_param.value
                        is_small_peak = True  # Fix sigma for small peaks
                    else:
                        # Allow larger area for Ll peaks in particles due to different attenuation envelope, potewntially with lower absorption than in the bulk
                        #TODO: In iterative fitting, this weight may be updated based on absorption instead of fixing it to a factor of 5 and letting the area vary
                        # This can be done with the adjust_peak_weights function in the EDS_Spectrum_Quantifier class, but it's not currently implemented for pileup and escape peaks
                        max_area = 5 * calibs.weight_Ll_ref_Ka1[self.meas_mode] * ref_peak_area_param.value
                else:
                    # All other reference peaks
                    max_area = np.inf
            
                if is_small_peak:
                    # Fix sigma for small peaks to avoid overfitting
                    params.add(line_prefix + self.sigma_key, value=sigma_init, vary=False)
                else:
                    params.add(sigma_par)
                
                # Add center and area parameters
                params.add(center_par)
                params.add(line_prefix + self.area_key, value=initial_area, min=0, max=max_area, vary=vary_area)
    
                # Add offset and broadening parameters for dependent peaks to use
                params.add(line_prefix + self.center_offset_key, expr=f"{line_energy} - {line_prefix}{self.center_key}")
                params.add(line_prefix + self.sigma_broadening_key, expr=f"{line_prefix}{self.sigma_key} / {sigma_init}")
    
                # Use calibrated peak shape models, or enable their calibration
                calibrate_peak_shape_params = el in self.free_peak_shapes_els
                peak_m, _ = self._get_peak_model_and_update_pars(
                    params, line_energy, line_prefix,
                    is_calibration=calibrate_peak_shape_params
                )

        # --- Add peak model to composite model if present ---
        if peak_m is not None:
            if model is None:
                model = peak_m  # Initialize model
            else:
                model += peak_m  # Add peak to composite model
    
        # --- Determine whether the peak is present or not ---
        if peak_m is None or (is_small_peak is None and not is_peak_overlapping):
            is_peak_present = False
        else:
            is_peak_present = True
    
        # --- Save updated model and parameters ---
        self.fitting_model = model
        self.fitting_params = params
    
        return is_peak_present
    
    
    def _get_peak_model_and_update_pars(self, params, line_energy, line_prefix, is_calibration=False, ref_line_prefix=None):
        """
        Select and configure the appropriate peak model and update shape parameters
        for a characteristic X-ray line, depending on the line energy and calibration mode.
    
        The peak shape model is selected based on the X-ray energy:
            - For energies < 1.18 keV: skewed Gaussian (Osan et al., 2001)
            - For 1.18 keV ≤ energy ≤ 1.839 keV: skewed Gaussian with ICC convolution
              (Osan et al., 2001 + Redus et al., 2015)
            - For 1.839 keV < energy ≤ 5 keV: Gaussian with ICC convolution (Redus et al., 2015)
            - For energies > 5 keV: plain Gaussian (ICC effects negligible)
    
        Parameters
        ----------
        params : lmfit.Parameters
            Shared parameter set to be updated in-place with peak shape parameters.
        line_energy : float
            Energy of the X-ray line (in keV).
        line_prefix : str
            Prefix for parameter names (typically includes element and line identifier).
        is_calibration : bool, optional
            If True, parameters are set to be varied for calibration; if False, use fixed/calibrated values.
        ref_line_prefix : str or None, optional
            If provided, dependent peak shape parameters are set as expressions tied to the reference line.
    
        Returns
        -------
        peak_m : lmfit.Model
            The peak model for the given energy range.
        params : lmfit.Parameters
            The updated parameter set (same object as input, returned for clarity).
    
        References
        ----------
        Osan, J. et al., X-ray Spectrom. 2001, 30, 8–17.
        Redus, R. H. et al., X-ray Spectrom. 2015, 44, 283–291.
        """
    
        # Get calibrated peak shape parameters for this energy and EDS mode
        gamma, f_tail, R_e, F_loss = calibs.get_calibrated_peak_shape_params(line_energy, self.meas_mode)
    
        if not is_calibration:
            # Use previously calibrated (fixed) peak shape parameters
            if line_energy < 1.18:
                # Skewed Gaussian
                peak_m = Model(Peaks_Model._gaussian_with_tail, prefix=line_prefix)
                params.add(line_prefix + self.gamma_key, value=gamma, vary=False)
                params.add(line_prefix + self.tail_fraction_key, value=f_tail, vary=False)
    
            elif line_energy <= 1.839:
                # Skewed Gaussian with ICC convolution
                peak_m = Model(Peaks_Model._gaussian_with_tail_and_icc, prefix=line_prefix)
                params.add(line_prefix + self.gamma_key, value=gamma, vary=False)
                params.add(line_prefix + self.tail_fraction_key, value=f_tail, vary=False)
                params.add(line_prefix + self.F_loss_key, value=F_loss, vary=False)
                params.add(line_prefix + self.R_e_key, value=R_e, vary=False)
    
            elif line_energy <= 5:
                # Gaussian with ICC convolution
                peak_m = Model(Peaks_Model._gaussian_with_icc, prefix=line_prefix)
                params.add(line_prefix + self.F_loss_key, value=F_loss, vary=False)
                params.add(line_prefix + self.R_e_key, value=R_e, vary=False)
    
            else:
                # Plain Gaussian: ICC effects negligible at high energy
                peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
    
        else:
            # Calibration mode: fit peak shape parameters (for reference or dependent peaks)
            if not ref_line_prefix:
                # Reference peaks: fit shape parameters within physical bounds
                if line_energy < 1.18:
                    peak_m = Model(Peaks_Model._gaussian_with_tail, prefix=line_prefix)
                    params.add(line_prefix + self.gamma_key, value=gamma, vary=True, min=1, max=6)
                    params.add(line_prefix + self.tail_fraction_key, value=f_tail, vary=True, min=0.0001, max=0.15)
    
                elif line_energy <= 1.839:
                    peak_m = Model(Peaks_Model._gaussian_with_tail_and_icc, prefix=line_prefix)
                    params.add(line_prefix + self.gamma_key, value=gamma, vary=True, min=1, max=4.5)
                    params.add(line_prefix + self.tail_fraction_key, value=f_tail, vary=True, min=0.001, max=0.15)
                    params.add(line_prefix + self.F_loss_key, value=F_loss, vary=True, min=0.01, max=0.5)
                    params.add(line_prefix + self.R_e_key, value=R_e, vary=True, min=1e-6, max=7e-5)
    
                elif line_energy <= 5:
                    peak_m = Model(Peaks_Model._gaussian_with_icc, prefix=line_prefix)
                    params.add(line_prefix + self.F_loss_key, value=F_loss, vary=True, min=0.01, max=0.5)
                    params.add(line_prefix + self.R_e_key, value=R_e, vary=True, min=1e-7, max=7e-5)
    
                else:
                    peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
            else:
                # Dependent peaks: tie shape parameters to reference peak using expressions
                is_f_tail_param = params.get(ref_line_prefix + self.tail_fraction_key) is not None
                is_F_loss_param = params.get(ref_line_prefix + self.F_loss_key) is not None
    
                if is_f_tail_param and is_F_loss_param:
                    peak_m = Model(Peaks_Model._gaussian_with_tail_and_icc, prefix=line_prefix)
                    params.add(line_prefix + self.gamma_key, expr=ref_line_prefix + self.gamma_key)
                    params.add(line_prefix + self.tail_fraction_key, expr=ref_line_prefix + self.tail_fraction_key)
                    params.add(line_prefix + self.F_loss_key, expr=ref_line_prefix + self.F_loss_key)
                    params.add(line_prefix + self.R_e_key, expr=ref_line_prefix + self.R_e_key)
                elif is_f_tail_param:
                    peak_m = Model(Peaks_Model._gaussian_with_tail, prefix=line_prefix)
                    params.add(line_prefix + self.gamma_key, expr=ref_line_prefix + self.gamma_key)
                    params.add(line_prefix + self.tail_fraction_key, expr=ref_line_prefix + self.tail_fraction_key)
                elif is_F_loss_param:
                    peak_m = Model(Peaks_Model._gaussian_with_icc, prefix=line_prefix)
                    params.add(line_prefix + self.F_loss_key, expr=ref_line_prefix + self.F_loss_key)
                    params.add(line_prefix + self.R_e_key, expr=ref_line_prefix + self.R_e_key)
                else:
                    peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
    
        return peak_m, params
    
    # =============================================================================
    # Fix parameters of overlapping peaks
    # =============================================================================
    def _fix_overlapping_ref_peaks(self):
        """
        Identify and constrain overlapping reference peaks by tying their center offset and sigma broadening parameters.
    
        This improves fitting stability in cases where reference peaks are close in energy
        (i.e., their centers are within ~3σ), by allowing them to shift together and
        maintain their relative distances, rather than fixing their positions outright.
        This is especially helpful if detector calibration is slightly off.
    
        Updates self.fitting_params in-place and tracks fixed dependencies in self.fixed_peaks_dict.
    
        Notes
        -----
        - Only reference peaks (those with variable center parameters) are considered.
        - Peaks are considered overlapping if their theoretical centers are within ~3σ.
        - For each overlapping pair, dependencies are set so that center offset and sigma broadening
          are shared, minimizing overfitting and improving robustness.
        """
    
        params = self.fitting_params
    
        # Find all center parameters
        center_params = [pname for pname in params if self.center_key in pname]
    
        # Reference peaks: only peaks whose center can vary
        free_peaks = {}
        for pname in center_params:
            if params[pname].vary:
                # Remove '_center' (and preceding underscore) to get peak prefix
                peak_prefix = pname[:-(len(self.center_key) + 1)]
                free_peaks[peak_prefix] = params[pname].value
    
        # Identify overlapping peaks: centers closer than 3σ
        peaks_to_fix = set()
        for peak1, peak2 in combinations(free_peaks, 2):
            center1 = free_peaks[peak1]
            center2 = free_peaks[peak2]
            sigma1 = DetectorResponseFunction._det_sigma(center1)
            # If centers are closer than 3σ, consider them overlapping
            if abs(center1 - center2) < sigma1 * 3:
                peaks_to_fix.add((peak1, peak2))
    
        # Track which peaks have already been fixed (dependent)
        if not hasattr(self, 'fixed_peaks_dict'):
            self.fixed_peaks_dict = {}
    
        # For each overlapping pair, tie center offset and sigma broadening
        for peak1, peak2 in peaks_to_fix:
            fixed_peaks = list(self.fixed_peaks_dict.keys())
            # Neither peak fixed: tie peak2 to peak1
            if peak1 not in fixed_peaks and peak2 not in fixed_peaks:
                params = self._fix_center_sigma_peak(params, peak1, peak2)
            # Both already fixed: nothing to do
            elif peak1 in fixed_peaks and peak2 in fixed_peaks:
                continue
            # One fixed, one not: tie the unfixed to the fixed (or its reference)
            else:
                if peak1 in fixed_peaks:
                    fixed_peak = peak1
                    dep_peak = peak2
                else:
                    fixed_peak = peak2
                    dep_peak = peak1
                # If fixed_peak is independent, use it as reference
                if self.fixed_peaks_dict[fixed_peak] == '':
                    params = self._fix_center_sigma_peak(params, fixed_peak, dep_peak)
                else:
                    # Use the ultimate reference
                    params = self._fix_center_sigma_peak(params, self.fixed_peaks_dict[fixed_peak], dep_peak)
    
        self.fitting_params = params
    
    
    def _fix_center_sigma_peak(self, params, ref_peak, dep_peak):
        """
        Tie the center and sigma of a dependent peak to those of a reference (independent) peak.
    
        - The dependent peak's center is constrained via an expression involving the reference peak's center offset.
        - Both peaks' sigmas are fixed (not varied) to avoid overfitting in the overlap region.
        - The starting area of both peaks is halved to compensate for initial overestimation due to overlap.
        - Records the dependency relationship in self.fixed_peaks_dict for bookkeeping.
    
        Parameters
        ----------
        params : lmfit.Parameters
            The parameters object to update.
        ref_peak : str
            The prefix for the reference (independent) peak.
        dep_peak : str
            The prefix for the dependent peak.
    
        Returns
        -------
        params : lmfit.Parameters
            The updated parameters object (modified in-place).
        """
    
        # Get the current center value of the dependent peak
        dep_peak_center = params[f"{dep_peak}_{self.center_key}"].value
    
        # Constrain dependent peak's center using the reference peak's center offset
        ref_peak_offset = params[f"{ref_peak}_{self.center_offset_key}"].name
        center_expr = f"{dep_peak_center} - {ref_peak_offset}"
        params[f"{dep_peak}_{self.center_key}"].expr = center_expr
    
        # NOTE: You may consider fixing sigma broadening instead of sigma directly.
        # The code below fixes sigma for both peaks (not varied in fit).
        params[f"{ref_peak}_{self.sigma_key}"].vary = False
        params[f"{dep_peak}_{self.sigma_key}"].vary = False
    
        # Optionally, halve the initial area for both peaks to reduce overestimation due to overlap
        params[f"{ref_peak}_{self.area_key}"].value /= 2
        params[f"{dep_peak}_{self.area_key}"].value /= 2
    
        # Bookkeeping: store the dependency relationship
        self.fixed_peaks_dict[ref_peak] = ''        # Reference peak is independent
        self.fixed_peaks_dict[dep_peak] = ref_peak  # Dependent peak is tied to reference
    
        return params
        
    
    

        
    
#%% Background_Model class
class Background_Model:
    """
    Model for calculating and fitting X-ray background in EDS spectra.

    Class Attributes
    ----------------
    cls_beam_e : float or None
        Electron beam energy in keV, accessible by all instances and fitting routines.
    den_int : any
        Cached denominator integral for background calculation.
    num_int : any
        Cached numerator integral for background calculation.
    prev_x : any
        Cached previous energy values.
    prev_rhoz_par_offset : any
        Cached previous absolute rho-z offset.
    prev_rhoz_par_slope : any
        Cached previous rho-z offset slope.
    prev_rhoz_limit : any
        Cached previous rho-z z limit.
    prev_w_frs : any
        Cached previous weight fractions.
    rhoz_values : any
        rhoz_values computed for a given rhoz_limit.

    Instance Attributes
    -------------------
    is_particle : bool
        If True, indicates the background is for a particle (affects fitting).
    sp_collection_time : float or None
        Spectrum collection time in seconds.
    tot_sp_counts : int or None
        Total counts in the spectrum.
    emergence_angle : float
        Detector emergence angle in degrees.
    energy_vals : array-like or None
        Array of energy values for the spectrum.
    meas_mode : str
        EDS mode, e.g., 'point' or 'map'.
    els_w_fr : dict
        Elemental weight fractions.
    """

    # Class variables for shared/cached state
    # Defined as a class variable so it can be accessed and updated by static methods,
    # including those used by lmfit models, which do not access instance variables.
    cls_beam_e: float = None
    den_int = None
    num_int = None
    prev_x = None
    prev_rhoz_par_offset = None
    prev_rhoz_par_slope = None
    prev_rhoz_limit = None
    prev_w_frs = None
    rhoz_values = None

    def __init__(
        self,
        is_particle: bool,
        sp_collection_time: float = None,
        tot_sp_counts: int = None,
        beam_energy: float = 15,
        emergence_angle: float = 28.5,
        els_w_fr: dict = None,
        meas_mode: str = 'point',
        energy_vals=None
    ):
        """
        Initialize a Background_Model instance.

        Parameters
        ----------
        is_particle : bool
            If True, indicates the background is for a particle (affects fitting).
        sp_collection_time : float, optional
            Spectrum collection time in seconds.
        tot_sp_counts : int, optional
            Total counts in the spectrum.
        beam_energy : float, optional
            Electron beam energy in keV (default is 15).
        emergence_angle : float, optional
            Detector emergence angle in degrees (default is 28.5).
        els_w_fr : dict, optional
            Elemental weight fractions.
        meas_mode : str, optional
            EDS mode, e.g., 'point' or 'map' (default is 'point').
        energy_vals : array-like, optional
            Array of energy values for the spectrum.
        """
        # Set class attribute for beam energy
        if beam_energy is not None:
            type(self).cls_beam_e = beam_energy

        self.is_particle = is_particle
        self.sp_collection_time = sp_collection_time
        self.tot_sp_counts = tot_sp_counts
        self.emergence_angle = emergence_angle
        self.energy_vals = energy_vals
        self.meas_mode = meas_mode
        self.els_w_fr = els_w_fr if els_w_fr is not None else {}

        # Reset the class-level caches to ensure they are recalculated in new fits
        self._clear_cached_abs_att_variables()

    @staticmethod
    def _clear_cached_abs_att_variables():
        """
        Reset class-level caches and variables to ensure they are recalculated in new fits.
    
        This method sets all relevant class attributes to None. It should be called before starting
        a new fit or when the cached values are no longer valid.
        """
        Background_Model.den_int = None
        Background_Model.num_int = None
        Background_Model.prev_x = None
        Background_Model.prev_rhoz_par_offset = None
        Background_Model.prev_rhoz_par_slope = None
        Background_Model.prev_rhoz_limit = None
        Background_Model.prev_w_frs = None
        Background_Model.rhoz_values = None
    
    
    @staticmethod
    def _get_els_frs(**el_fr_params):
        """
        Parses element weight fraction parameters and returns element symbols, weight fractions, and atomic fractions.
    
        Parameters
        ----------
        el_fr_params : dict
            Elemental weight fractions in the form {'f_Si': weight_fraction, ...},
            where weight_fraction can be a float or an lmfit.Parameter.
    
        Returns
        -------
        els : list of str
            List of element symbols (e.g., ["Si", "O"]).
        w_frs : list of float
            Corresponding weight fractions for each element.
        at_frs : list of float
            Corresponding atomic fractions for each element, calculated from weight fractions.
        """
        par_pattern = r'f_[A-Z][a-z]*'  # regex to filter the useful parameters
        els = []
        w_frs = []
        for el_fr_param, w_fr in el_fr_params.items():
            if re.match(par_pattern, el_fr_param):
                el = el_fr_param.split("_")[1]
                els.append(el)
                # Accept float or lmfit.Parameter
                val = w_fr.value if hasattr(w_fr, 'value') else w_fr
                w_frs.append(val)

        if len(els) > 0:
            at_frs = weight_to_atomic_fr(w_frs, els, verbose=False)  # Avoid warning if not normalized
        else:
            at_frs = []
    
        return els, w_frs, at_frs
        
    # =============================================================================
    # Atomic number averaging
    # =============================================================================
    @staticmethod
    def get_average_Z(els_symbols, method="Statham"):
        """
        Returns a symbolic formula string for the average atomic number (Z) of a sample,
        using one of several literature methods.
    
        Parameters
        ----------
        els_symbols : list of str
            List of element symbols (e.g., ['Si', 'O']).
        method : str, optional
            Which method to use for averaging Z. Options are:
            - "mass_weighted": mass-fraction-weighted average Z (sum(Z_i * f_i))
            - "Markowicz": Markowicz & Van Grieken (1984) average Z formula
            - "Statham": Statham et al. (2016) average Z formula
            Default is "mass_weighted".
    
        Returns
        -------
        formula_str : str
            Symbolic formula for the chosen average Z, using mass fractions f_{el}.
    
        References
        ----------
        Markowicz AA, Van Grieken RE. Anal Chem 1984 Oct 1;56(12):2049–51.
        Statham P, Penman C, Duncumb P. IOP Conf Ser Mater Sci Eng. 2016;109(1):0–10.
    
        Notes
        -----
        - This function returns a string formula; you must substitute actual fractions for f_{el} in use.
        - Only the Statham method is currently used in the main implementation.
        """
        els = [Element(el) for el in els_symbols]
        if len(els) == 1:
            el = els[0]
            return f"{el.Z} * f_{el.symbol}"
    
        method = method.lower()
        if method == "mass_weighted":
            return " + ".join(f"{el.Z} * f_{el.symbol}" for el in els)
        elif method == "markowicz":
            num = " + ".join(f"{el.Z**2/el.atomic_mass:.3f} * f_{el.symbol}" for el in els)
            den = " + ".join(f"{el.Z/el.atomic_mass:.3f} * f_{el.symbol}" for el in els)
            return f"({num}) / ({den})"
        elif method == "statham":
            num = " + ".join(f"{el.Z**1.75/el.atomic_mass:.3f} * f_{el.symbol}" for el in els)
            den = " + ".join(f"{el.Z**0.75/el.atomic_mass:.3f} * f_{el.symbol}" for el in els)
            return f"({num}) / ({den})"
        else:
            raise ValueError(f"Unknown method '{method}' for averaging of compound atomic number. Choose from 'mass_weighted', 'Markowicz', or 'Statham'.")
        
        
    # =============================================================================
    # Stopping power
    # =============================================================================
    @staticmethod
    def _stopping_power(x, adr_sp=1, **el_fr_params):
        """
        Computes the stopping power correction for a multi-element sample.
    
        Parameters
        ----------
        x : array-like
            Array of incident electron energies (eV or keV, as appropriate).
        adr_sp : int, optional
            If 1, applies the detector response correction function. Default is 1.
        **el_fr_params : dict
            Elemental weight fractions in the form f_Symbol=weight_fraction, e.g., f_Si=0.5.
    
        Returns
        -------
        S_vals : numpy.ndarray
            Stopping power correction values, same shape as `x`.
    
        References
        ----------
        - G. Love, V.D. Scott, J. Phys. D. Appl. Phys. 11 (1978) 1369–1376.
    
        Notes
        -----
        The stopping power is calculated using the method of Love & Scott (1978),
        with the average ionization potential computed over all elements in the sample.
    
        CURRENTLY NOT USED: The function currently returns an array of ones (i.e., no stopping power correction is applied).
        """
        M_vals = []
        lnJ_vals = []
        
        els, w_frs, _ = Background_Model._get_els_frs(**el_fr_params)
        
        # Collect mass-weighted Z and log(J) values
        for el_symbol, w_fr in zip(els, w_frs):
            el = Element(el_symbol)
            Z_el = el.Z
            W_el = el.atomic_mass
            J_el = J_df.loc[Z_el, J_df.columns[0]] / 1000  # J in keV
            M = w_fr * Z_el / W_el
            M_vals.append(M)
            lnJ_vals.append(M * np.log(J_el))
        
        sum_M = sum(M_vals)
        ln_J = sum(lnJ_vals) / sum_M  # Average log ionization potential
        J_val = np.exp(ln_J)          # Average ionization potential
        U0 = Background_Model.cls_beam_e / x
    
        # Love & Scott (1978) stopping power formula
        S_vals = 1 / ((1 + 16.05 * (J_val / x) ** 0.5 * ((U0 ** 0.5 - 1) / (U0 - 1)) ** 1.07) / sum_M)
        
        # Optional: apply detector response function
        if adr_sp == 1:
            S_vals = DetectorResponseFunction._apply_det_response_fncts(S_vals)
        
        # CURRENTLY NOT USED: Return array of ones
        S_vals = np.ones(len(x))
        
        return S_vals
    
    
    def get_stopping_power_mod_pars(self):
        """
        Returns an lmfit Model and its Parameters for the electron stopping power correction.
    
        Returns
        -------
        stopping_p_correction_m : lmfit.Model
            Model describing the stopping power correction.
        params_stopping_p : lmfit.Parameters
            Parameters for the model.
    
        Notes
        -----
        - The model is based on Background_Model._stopping_power.
        - The adr_sp parameter controls convolution with the detector response function.
        """
        stopping_p_correction_m = Model(Background_Model._stopping_power)
        params_stopping_p = stopping_p_correction_m.make_params(
            adr_sp=dict(expr='apply_det_response')
        )
    
        return stopping_p_correction_m, params_stopping_p


    # =============================================================================
    # X-ray absorption attenuation
    # =============================================================================
    @staticmethod
    def _mass_abs_coeff(x, **el_fr_params):
        """
        Computes the total mass absorption coefficient (μ/ρ) for a compound or mixture,
        as the sum of elemental mass absorption coefficients weighted by their mass fractions.
    
        Mass absorption coefficients are retrieved from the Henke database:
        https://henke.lbl.gov/optical_constants/asf.html
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV).
        **el_fr_params : dict
            Elemental weight fractions, e.g., f_Si=0.5.
    
        Returns
        -------
        mass_abs_coeff : np.ndarray
            Total mass absorption coefficient (cm^2/g), same shape as x.
        """
        mass_abs_coeff = np.zeros(len(x))
        els, w_frs, _ = Background_Model._get_els_frs(**el_fr_params)
        for el, w_fr in zip(els, w_frs):
            mass_abs_coeff += xray_mass_absorption_coeff(el, x) * w_fr  # units: cm^2/g
        return mass_abs_coeff


    @staticmethod
    def _abs_attenuation_Philibert(
        x, det_angle=28.5, abs_par=1.2e-6, abs_path_len_scale=1, adr_abs=1, **el_fr_params
    ):
        """
        Computes the absorption attenuation correction using a modified Philibert equation.
    
        Reference
        ---------
        K.F.J. Heinrich, H. Yakowitz,
        "Absorption of Primary X Rays in Electron Probe Microanalysis",
        Anal. Chem. 47 (1975) 2408–2411. https://doi.org/10.1021/ac60364a018
    
        Notes
        -----
        - Strictly, this equation is correct only when fluorescence from the continuum occurs,
          which is not necessarily the case in particles. The primary excitation fraction should be used,
          not the total.
        - Also, E_q is not the energy, but rather the critical excitation potential; so this equation
          does not strictly apply to continuum. However, after Trincavelli (1998), it is acceptable to use
          E for the continuum, as the minimum energy to excite the continuum is E itself; there is no
          critical energy as in characteristic X-rays. These differences are discussed by Small (1987)
          and Trincavelli (1998).
        - abs_par is an empirical parameter, typically ~1.2e-6.
        - abs_path_len_scale can be used to scale the absorption path length.
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV).
        det_angle : float, optional
            Detector takeoff angle in degrees. Default is 28.5.
        abs_par : float, optional
            Empirical absorption parameter. Default is 1.2e-6.
        abs_path_len_scale : float, optional
            Scaling factor for absorption path length. Default is 1.
        adr_abs : int, optional
            If 1, convolves signal with detector response function. Default is 1.
        **el_fr_params : dict
            Elemental weight fractions, e.g., f_Si=0.5.
    
        Returns
        -------
        model : np.ndarray
            Absorption correction factor, same shape as `x`.
        """
        E0 = Background_Model.cls_beam_e
        chi = (
            Background_Model._mass_abs_coeff(x, **el_fr_params)
            / np.sin(np.deg2rad(det_angle))
            * abs_path_len_scale
        )
        gamma = E0**1.65 - x**1.65
        model = (1 + abs_par * gamma * chi) ** -2
        
        if adr_abs == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    
    
    @staticmethod
    def _A_pb(x, mass_abs_coeffs_sample, det_angle, beam_energy):
        """
        Second-order absorption correction for P/B ratio to account for differences in mean generation depths of
        characteristic X-rays and continuum.
    
        Reference
        ---------
        P.J. Statham,
        "A ZAF PROCEDURE FOR MICROPROBE ANALYSIS BASED ON MEASUREMENT OF PEAK-TO-BACKGROUND RATIOS",
        in: D.E. Newbury (Ed.), Fourteenth Annu. Conf. Microbeam Anal. Soc.,
        San Francisco Press, San Francisco, 1979: pp. 247–253.
        https://archive.org/details/1979-mas-proc-san-antonio/page/246/mode/2up
    
        This correction factor needs to be multiplied by the absorption attenuation of characteristic X-rays to obtain
        the attenuation for continuum X-rays.
        Because the depth of generation of the continuum is larger than that of characteristic X-rays, the continuum
        gets absorbed more, and thus A_pb < 1.
    
        Note
        ----
        This function is **not used** in the current implementation because it was found to be inadequate at low energies,
        especially for high-Z elements.
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV).
        mass_abs_coeffs_sample : array-like
            Sample mass absorption coefficients (cm^2/g), same shape as `x`.
        det_angle : float
            Detector takeoff angle (degrees).
        beam_energy : float
            Incident electron beam energy (keV).
    
        Returns
        -------
        A_vals : np.ndarray
            Absorption correction factor for P/B ratio (always < 1).
        """
        chi = mass_abs_coeffs_sample / np.sin(np.deg2rad(det_angle))
        gamma = beam_energy**1.65 - x**1.65
    
        # Absorption correction for characteristic X-ray
        A_P = 1 + 3.34e-6 * (gamma * chi) + 5.59e-13 * (gamma * chi)**2
    
        # Absorption correction for continuum
        A_B = 1 + 3.0e-6 * (gamma * chi) + 4.5e-13 * (gamma * chi)**2
    
        # Multiplicative factor for PB ratio. Always < 1, because the absorption losses in the background are higher compared to the continuum
        A_vals = A_B / A_P
    
        return A_vals
    
    
    @staticmethod
    def _abs_attenuation_phirho(
        x, det_angle=28.5, rhoz_par_offset=0, rhoz_par_slope=0, rhoz_lim=0.001, adr_abs=1, **el_fr_params
    ):
        """
        Absorption correction based on the ionization depth distribution function (ϕ(ρz)) 
        following Packwood & Brown (1981), with coefficients modified by Riveros et al. (1992)
        and simplified equations from Del Giorgio et al. (1990). The absorption coefficient 
        formula is from Statham (1976).
    
        References
        ----------
        - R.H. Packwood, J.D. Brown, X-Ray Spectrom. 10 (1981) 138–146. https://doi.org/10.1002/xrs.1300100311
        - J.A. Riveros et al., in: 1992: pp. 99–105. https://doi.org/10.1007/978-3-7091-6679-6_7
        - M. Del Giorgio et al., X-Ray Spectrom. 19 (1990) 261–267. https://doi.org/10.1002/xrs.1300190603
        - Statham PJ, X-Ray Spectrom 5 (1976) 154–68. https://doi.org/10.1002/xrs.1300050310
    
        Parameters
        ----------
        x : array-like
            Emission energies (keV).
        det_angle : float, optional
            Detector takeoff angle (degrees). Default is 28.5.
        rhoz_par_offset, rhoz_par_slope : float, optional
            Parameters for absorption path correction.
        rhoz_lim : float, optional
            Upper limit for ρz integration. Default is 0.001.
        adr_abs : int, optional
            If 1, convolves signal with detector response function. Default is 1.
        **el_fr_params : dict
            Elemental weight fractions, e.g. f_Si=0.5.
    
        Returns
        -------
        abs_model : np.ndarray
            Absorption correction factor (shape matches x).
    
        Notes
        -----
        For development and diagnostic purposes, this function can optionally plot the
        φ(ρz) curves and associated absorption integrands, to visually inspect the behavior
        of the correction as a function of energy and composition.
        """
        # Get list of elements and corresponding mass fractions 
        els, w_frs, at_frs = Background_Model._get_els_frs(**el_fr_params)
        
        E0 = Background_Model.cls_beam_e
        U0 = E0 / x  # Incident overvoltage
        
        # Get sample mass absorption coefficient
        mu = Background_Model._mass_abs_coeff(x, **el_fr_params)
        
        # Calculate backscattering coefficient for sample, averaged on mass fractions
        nu = Background_Model._nu_sample(E0, els, w_frs)
        
        # Calculate phi0, gamma0 according to Del Giorgio
        phi0 = 1 + (nu * U0 * np.log(U0)) / (U0 - 1)
        gamma0 = (1 + nu) * (U0 * np.log(U0)) / (U0 - 1)
        
        # Calculate alpha, beta as atomic-fraction averages
        alpha = 0
        beta = 0
        for el, at_fr in zip(els, at_frs):
            a_el, b_el = Background_Model._phi_alpha_beta_coeffs(x, E0, el, gamma0)
            alpha += a_el * at_fr
            beta += b_el * at_fr 
    
        recalc_den = False
        recalc_num = False
    
        # Check if any variable that requires recalculation of den and num has changed
        if (
            x is not Background_Model.prev_x
            or Background_Model.prev_rhoz_limit != rhoz_lim
            or Background_Model.prev_w_frs != w_frs
        ):
            Background_Model.prev_x = x
            Background_Model.prev_rhoz_limit = rhoz_lim
            Background_Model.prev_w_frs = w_frs
            Background_Model.rhoz_values = np.linspace(0, rhoz_lim, 10**3)
            recalc_den = True
            recalc_num = True
        elif (
            Background_Model.prev_rhoz_par_offset != rhoz_par_offset
            or Background_Model.prev_rhoz_par_slope != rhoz_par_slope
        ):
            Background_Model.prev_rhoz_par_offset = rhoz_par_offset
            Background_Model.prev_rhoz_par_slope = rhoz_par_slope
            recalc_num = True
    
        if recalc_num:
            Background_Model.num_int = []
            for a, b, p0, g0, m in zip(alpha, beta, phi0, gamma0, mu):
                num_int_val = Background_Model._get_abs_att_num(
                    a, b, p0, g0, m, det_angle, rhoz_par_offset, rhoz_par_slope
                )
                Background_Model.num_int.append(num_int_val)
            Background_Model.num_int = np.array(Background_Model.num_int)
    
        if recalc_den:
            Background_Model.den_int = []
            for a, b, p0, g0 in zip(alpha, beta, phi0, gamma0):
                den_int_val = Background_Model._get_abs_att_den(a, b, p0, g0, det_angle)
                Background_Model.den_int.append(den_int_val)
            Background_Model.den_int = np.array(Background_Model.den_int)
    
        abs_model = Background_Model.num_int / Background_Model.den_int
    
        # Clip values to 1 to avoid artificial gain in intensity.
        # Not needed if the max and min of offset and slope are well fixed during fitting
        abs_model[abs_model > 1] = 1
    
        # === Development/diagnostic plotting of φ(ρz) curves ===
        # # To visually inspect the model, plot for a few representative energies.
        # # Remove or comment out in production.
        # ens_to_plot = np.linspace(np.min(x), np.max(x), 3)
        # for en_to_plot in ens_to_plot:
        #     alpha_val = np.interp(en_to_plot, x, alpha)
        #     beta_val = np.interp(en_to_plot, x, beta)
        #     phi0_val = np.interp(en_to_plot, x, phi0)
        #     gamma0_val = np.interp(en_to_plot, x, gamma0)
        #     mu_val = np.interp(en_to_plot, x, mu)
        #     rhoz_grid = Background_Model.rhoz_values
    
        #     num_vals = (
        #         Background_Model._phi_rhoz(rhoz_grid, alpha_val, beta_val, phi0_val, gamma0_val)
        #         * np.exp(-mu_val * (rhoz_grid + rhoz_par_offset + rhoz_grid * rhoz_par_slope) / np.sin(np.deg2rad(det_angle)))
        #     )
        #     den_vals = Background_Model._phi_rhoz(rhoz_grid, alpha_val, beta_val, phi0_val, gamma0_val)
        #     plt.figure(figsize=(9, 5))
        #     plt.plot(rhoz_grid, den_vals, '--', label='Denominator Integrand (no absorption)', color='tab:red')
        #     plt.plot(rhoz_grid, num_vals, '-', label='Numerator Integrand (with absorption)', color='tab:blue')
        #     plt.title(f'ϕ(ρz) Integrands at E = {en_to_plot:.2f} keV')
        #     plt.xlabel('ρz')
        #     plt.ylabel('Integrand Value')
        #     plt.legend()
        #     plt.grid(True, linestyle=':')
        #     plt.tight_layout()
        #     plt.show()
    
        if adr_abs == 1:
            abs_model = DetectorResponseFunction._apply_det_response_fncts(abs_model)
    
        return abs_model
    
    
    @staticmethod
    def _get_abs_att_num(alpha, beta, phi0, gamma0, mu, det_angle, rhoz_par_offset, rhoz_par_slope):
        """
        Computes the numerator integral for the absorption-attenuation correction
        using the φ(ρz) (phi-rhoz) double Gaussian model, including attenuation.
    
        Parameters
        ----------
        alpha, beta, phi0, gamma0 : float
            Parameters for the φ(ρz) double Gaussian model.
        mu : float
            Linear attenuation coefficient (units consistent with ρz).
        det_angle : float
            Detector takeoff angle in degrees.
        rhoz_par_offset, rhoz_par_slope : float
            Parameters for the absorption path correction.
    
        Returns
        -------
        num_int : float
            Value of the numerator integral over ρz, including attenuation correction.
    
        Notes
        -----
        The numerator is calculated as:
            ∫ φ(ρz) * exp(-μ * (ρz + offset + ρz * slope) / sin(det_angle)) d(ρz)
        over the grid defined by `Background_Model.rhoz_values`.
        """
        rhoz_grid = np.array(Background_Model.rhoz_values)
        num_vals = (
            Background_Model._phi_rhoz(rhoz_grid, alpha, beta, phi0, gamma0)
            * np.exp(-mu * (rhoz_grid + rhoz_par_offset + rhoz_grid * rhoz_par_slope) / np.sin(np.deg2rad(det_angle)))
        )
        num_int = trapezoid(num_vals, rhoz_grid)
        return num_int
        
    
    @staticmethod
    def _get_abs_att_den(alpha, beta, phi0, gamma0, det_angle):
        """
        Computes the denominator integral for the absorption-attenuation correction
        using the φ(ρz) (phi-rhoz) Gaussian model.

        This integral represents the total generated characteristic X-rays as a function
        of depth, prior to absorption correction.

        Parameters
        ----------
        alpha, beta, phi0, gamma0 : float
            Parameters for the φ(ρz) double Gaussian model.
        det_angle : float
            Detector takeoff angle in degrees (not used in denominator calculation, 
            but included for API consistency).

        Returns
        -------
        den_int : float
            Value of the denominator integral over ρz.

        Notes
        -----
        The denominator is calculated as:
            ∫ φ(ρz) d(ρz)
        over the grid defined by `Background_Model.rhoz_values`.

        This is used in matrix correction procedures and absorption calculations.
        """
        rhoz_grid = np.array(Background_Model.rhoz_values)
        den_vals = Background_Model._phi_rhoz(rhoz_grid, alpha, beta, phi0, gamma0)
        den_int = trapezoid(den_vals, rhoz_grid)
        return den_int
    
    
    @staticmethod
    def _plot_phirho(
        en, alpha, beta, phi0, gamma0, mu,
        rhoz_par_offset, rhoz_par_slope, det_angle, rhoz_lim
    ):
        """
        Visualizes the φ(ρz) (phi-rhoz) curves and their absorption-corrected integrands.

        This method is intended **solely for development and visualisation** of φ(ρz) and 
        absorption integrands, to facilitate comparison with literature (e.g., Bastin et al., 1998) 
        and to aid the tuning and understanding of the matrix correction model. 
        It is not used in production calculations.

        Parameters
        ----------
        en : float
            Incident electron energy (keV).
        alpha, beta, phi0, gamma0 : float
            Parameters for the φ(ρz) double Gaussian model.
        mu : float
            Linear attenuation coefficient.
        rhoz_par_offset, rhoz_par_slope : float
            Parameters for the absorption path correction.
        det_angle : float
            Detector takeoff angle (degrees).
        rhoz_lim : float
            Upper limit of ρz for integration and plotting.

        Returns
        -------
        ratio_abs_change : float
            Ratio of normalized absorption integrals (cut/full), useful for diagnostics.

        Notes
        -----
        The integrands and their areas are compared to published results:
        See: G.F. Bastin, J.M. Dijkstra, H.J.M. Heijligers, 
            "PROZA96: an improved matrix correction program for electron probe microanalysis, 
            based on a double Gaussian ϕ(ρz) approach", X-Ray Spectrom. 27 (1998) 3–10.
            https://doi.org/10.1002/(SICI)1097-4539(199801/02)27:1<3::AID-XRS227>3.0.CO;2-L

        This function is for developer use only and should not be called in production code.
        """

        # Define integrand functions
        num_integrand = lambda rhoz: (
            Background_Model._phi_rhoz(rhoz, alpha, beta, phi0, gamma0)
            * np.exp(-mu * (rhoz + rhoz_par_offset + rhoz * rhoz_par_slope) 
            / np.sin(np.deg2rad(det_angle)))
        )
        den_integrand = lambda rhoz: Background_Model._phi_rhoz(rhoz, alpha, beta, phi0, gamma0)

        # Prepare ρz grid and integration limits
        rhoz_grid = np.array(Background_Model.rhoz_values)
        i_lim = np.argmin(np.abs(rhoz_grid - rhoz_lim))
        rhoz_cut = rhoz_grid[:i_lim]

        # Compute integrand values
        num_values = [num_integrand(rhoz) for rhoz in rhoz_grid]
        den_values = [den_integrand(rhoz) for rhoz in rhoz_grid]
        num_values_cut = [num_integrand(rhoz) for rhoz in rhoz_cut]
        den_values_cut = [den_integrand(rhoz) for rhoz in rhoz_cut]

        # Integrate using trapezoidal rule
        num_area = trapezoid(num_values, rhoz_grid) * 1000
        den_area = trapezoid(den_values, rhoz_grid) * 1000
        num_area_cut = trapezoid(num_values_cut, rhoz_cut)
        den_area_cut = trapezoid(den_values_cut, rhoz_cut)

        # Compute diagnostic ratio
        ratio_abs_change = (num_area_cut / den_area_cut) / (num_area / den_area)

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(rhoz_grid, num_values, label='Numerator Integrand (with absorption)', color='b')
        plt.plot(rhoz_grid, den_values, label='Denominator Integrand (no absorption)', color='r')
        text_x = rhoz_grid[int(len(rhoz_grid) * 0.7)]
        text_y = max(den_values) * 0.8
        plt.text(
            text_x, text_y,
            f"Num Integral (Full): {num_area:.3f}\n"
            f"Den Integral (Full): {den_area:.3f}\n",
            fontsize=10, bbox=dict(facecolor='white', alpha=0.7)
        )
        plt.title(f'Numerator and Denominator Integrands vs ρz at {en:.3f} keV')
        plt.xlabel('ρz')
        plt.ylabel('Integrand Value')
        plt.legend()
        plt.grid(True)
        plt.show()

        return ratio_abs_change
    
    @staticmethod
    def _phi_rhoz(rhoz, alpha, beta, phi0, gamma0):
        """
        Computes the ionization depth distribution function ϕ(ρz) as given by:
        R.H. Packwood, J.D. Brown, "A Gaussian expression to describe ϕ(ρz) curves for quantitative electron probe microanalysis",
        X-Ray Spectrom. 10 (1981) 138–146. https://doi.org/10.1002/xrs.1300100311
    
        Parameters
        ----------
        rhoz : float or np.ndarray
            Depth variable (ρz), typically in g/cm².
        alpha : float
            Gaussian width parameter.
        beta : float
            Exponential tail parameter.
        phi0 : float
            Surface value of ϕ(ρz).
        gamma0 : float
            Maximum value of ϕ(ρz).
    
        Returns
        -------
        phi_rhoz : float or np.ndarray
            Value(s) of the ϕ(ρz) distribution at the given depth(s).
        """
        phi_rhoz = gamma0 * np.exp(-alpha**2 * rhoz**2) * (1 - ((gamma0 - phi0) / gamma0) * np.exp(-beta * rhoz))
        return phi_rhoz
    
    
    @staticmethod
    def _phi_alpha_beta_coeffs(x, E0, el, gamma0):
        """
        Computes the alpha and beta coefficients for the φ(ρz) (ionization depth distribution) model.
    
        - Alpha is calculated after Packwood & Brown (1981):
          R.H. Packwood, J.D. Brown, X-Ray Spectrom. 10 (1981) 138–146. https://doi.org/10.1002/xrs.1300100311
        - Beta is calculated after Riveros (1993) and Tirira et al. (1987):
          J.H. Tirira Saa et al., X-Ray Spectrom. 16 (1987) 243–248. https://doi.org/10.1002/xrs.1300160604
          J. Riveros, G. Castellano, X-Ray Spectrom. 22 (1993) 3–10. https://doi.org/10.1002/xrs.1300220103
    
        Parameters
        ----------
        x : float or np.ndarray
            Emission energy (keV).
        E0 : float
            Incident electron energy (keV).
        el : str
            Element symbol (e.g., "Si").
        gamma0 : float
            Maximum value of φ(ρz) (not used in this calculation, but included for API consistency).
    
        Returns
        -------
        alpha : float or np.ndarray
            Gaussian width parameter for φ(ρz).
        beta : float or np.ndarray
            Exponential tail parameter for φ(ρz).
        """
        Z = Element(el).Z
        A = Element(el).atomic_mass
        J = J_df.loc[Z, J_df.columns[0]] / 1000  # keV
    
        # Alpha (Packwood, Tirira)
        alpha = 2.14e5 * Z**1.16 / (A * E0**1.25) * (np.log(1.166 * E0 / J) / (E0 - x))**0.5
    
        # Beta (Riveros 1993, reports the value from Tirira, but it's different)
        beta = 1.1e5 * Z**1.5 / ((E0 - x) * A)
        
        # Beta (Tirira 1987, there must be a mistake in their formula, because values are higer by 1 order of magnitude, depite being lower in their table)
        # beta = -2.73 * 10**5 * Z**1.5 * np.log(0.0180) / ((E0 - x) * A) #
    
        return alpha, beta
    
    
    def get_abs_attenuation_mod_pars(self, model='phirho'):
        """
        Returns an lmfit Model and its Parameters for background signal attenuation due to absorption in the sample.
    
        Parameters
        ----------
        model : str, optional
            Which absorption model to use. Options:
            - 'phirho' (default): depth-distribution-based attenuation
                  (Packwood & Brown, 1981; Riveros et al., 1992; Del Giorgio et al., 1990)
            - 'Philibert': modified Philibert equation
                  (Heinrich & Yakowitz, 1975)
    
        Returns
        -------
        absorption_attenuation_m : lmfit.Model
            Model function describing the attenuation of the background signal.
        params_abs_att : lmfit.Parameters
            Parameters for the model.
    
        References
        ----------
        - Packwood, R.H., & Brown, J.D. (1981). A Gaussian expression to describe ϕ(ρz) curves for quantitative electron probe microanalysis. X-Ray Spectrom. 10, 138–146. https://doi.org/10.1002/xrs.1300100311
        - Riveros, J.A., Castellano, G.E., & Trincavelli, J.C. (1992). Comparison of Φ (ρz) Curve Models in EPMA. In: Microbeam Analysis, pp. 99–105. https://doi.org/10.1007/978-3-7091-6679-6_7
        - Del Giorgio, M., Trincavelli, J., & Riveros, J.A. (1990). Spectral distribution of backscattered electrons: Application to electron probe microanalysis. X-Ray Spectrom. 19, 261–267. https://doi.org/10.1002/xrs.1300190603
        - Heinrich, K.F.J. & Yakowitz, H. (1975). Absorption of Primary X Rays in Electron Probe Microanalysis. Anal. Chem. 47, 2408–2411. https://doi.org/10.1021/ac60364a018
    
        Raises
        ------
        ValueError
            If an unknown model is requested.
        """
        
        model= model.lower()
        if model == 'phirho':
            absorption_attenuation_m = Model(
                Background_Model._abs_attenuation_phirho, independent_vars=['x']
            )
            params_abs_att = absorption_attenuation_m.make_params(
                det_angle={'value': self.emergence_angle, 'vary': False, 'min': 10, 'max': 90},
                adr_abs=dict(expr='apply_det_response'),
                rhoz_lim={'value': 0.001, 'vary': False, 'min': 0, 'max': 0.001},
                rhoz_par_slope={'value': 0, 'vary': self.is_particle, 'min': -1, 'max': 5},
                # a slope larger than 5 can lead to np.inf values and therefore errors in the fit
                rhoz_par_offset={'value': 0, 'vary': self.is_particle, 'min': -0.0005, 'max': 0.0005},
            )
        elif model == 'philibert':
            absorption_attenuation_m = Model(
                Background_Model._abs_attenuation_Philibert, independent_vars=['x']
            )
            params_abs_att = absorption_attenuation_m.make_params(
                det_angle={'value': self.emergence_angle, 'vary': False, 'min': 10, 'max': 90},
                adr_abs=dict(expr='apply_det_response'),
                abs_par={'value': 1.2e-6, 'vary': False, 'min': 0.5e-6, 'max': 2e-6},
                abs_path_len_scale={'value': 1, 'vary': self.is_particle, 'min': 0.01, 'max': 100},
            )
        else:
            raise ValueError(
                f"Unknown model '{model}' for absorption attenuation. "
                "Choose 'phirho' or 'Philibert'."
            )
    
        return absorption_attenuation_m, params_abs_att


    # =============================================================================
    # Electron backscattering correction
    # =============================================================================
    @staticmethod
    def _backscattering_correction(x, adr_bcksctr=1, **el_fr_params):
        """
        Computes the backscattering correction factor for electron probe microanalysis,
        following the method of Essani et al. (2020):
    
            M. Essani, E. Brackx, E. Excoffier,
            "A method for the correction of size effects in microparticles using a peak-to-background approach in electron-probe microanalysis",
            Spectrochim. Acta Part B At. Spectrosc. 169 (2020) 105880.
            https://doi.org/10.1016/j.sab.2020.105880
    
        Note
        ----
        This correction yields unphysical results below 100 eV (values >1).
    
        Parameters
        ----------
        x : array-like
            Array of incident or emission energies (keV).
        adr_bcksctr : int, optional
            If 1, convolves signal with detector response function. Default is 1.
        **el_fr_params : dict
            Elemental weight fractions, e.g., f_Si=0.5.
    
        Returns
        -------
        model : np.ndarray
            Backscattering correction factor, same shape as `x`.
        """
        E0 = Background_Model.cls_beam_e
    
        # Get list of elements and corresponding mass fractions
        els, w_frs, _ = Background_Model._get_els_frs(**el_fr_params)
    
        # Calculate backscattering coefficient for sample
        nu_val = Background_Model._nu_sample(E0, els, w_frs)
    
        # Calculate backscattering correction factor using Statham's formula
        R_c_val = Background_Model._R_c(x, nu_val, E0)
    
        # Backscattering loss correction model (Essani et al., 2020)
        model = 1 - (1 - R_c_val) * (2 / (1 + nu_val))**0.63 * (0.79 + 0.44 * x / E0)
    
        # Optionally apply detector response function
        if adr_bcksctr == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    
    
    @staticmethod
    def _nu_sample(E0, els, w_frs):
        """
        Computes the sample backscattering coefficient (ν_sample) as the mass-fraction-weighted
        sum of elemental backscattering coefficients, following Love & Scott (1978).
    
        Parameters
        ----------
        E0 : float
            Incident electron energy (keV).
        els : list of str
            List of element symbols (e.g., ['Si', 'O']).
        w_frs : list of float
            Corresponding mass fractions for each element.
    
        Returns
        -------
        nu_val : float
            Mass-fraction-averaged backscattering coefficient for the sample.
        """
        nu_val = 0.0
        for el, w_fr in zip(els, w_frs):
            Z = Element(el).Z
            nu_el = Background_Model._nu_el(Z, E0)
            nu_val += nu_el * w_fr
        return nu_val
    
    
    @staticmethod
    def _nu_el(Z, E0):
        """
        Computes the elemental backscattering coefficient ν for a given atomic number Z and incident energy E0,
        according to Love & Scott (1978):
    
            G. Love, V.D. Scott,
            "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
            J. Phys. D. Appl. Phys. 11 (1978) 1369–1376.
            https://doi.org/10.1088/0022-3727/11/10/002
    
        Parameters
        ----------
        Z : int
            Atomic number of the element.
        E0 : float
            Incident electron energy (keV).
    
        Returns
        -------
        nu_val : float
            Backscattering coefficient for the element.
        """
        nu20 = (-52.3791 + 150.48371 * Z - 1.67373 * Z**2 + 0.00716 * Z**3) * 1e-4
        G_nu20 = (-1112.8 + 30.289 * Z - 0.15498 * Z**2) * 1e-4
        nu_val = nu20 * (1 + G_nu20 * np.log(E0 / 20))
        return nu_val   

    
    @staticmethod
    def _R_c(x, nu_val, E0):
        """
        Computes the backscattering correction factor R_c as described in:
    
            M. Essani, E. Brackx, E. Excoffier,
            "A method for the correction of size effects in microparticles using a peak-to-background approach in electron-probe microanalysis",
            Spectrochim. Acta Part B At. Spectrosc. 169 (2020) 105880.
            https://doi.org/10.1016/j.sab.2020.105880
    
        Parameters
        ----------
        x : array-like
            Emission or incident energies (keV).
        nu_val : float
            Sample backscattering coefficient (dimensionless).
        E0 : float
            Incident electron energy (keV).
    
        Returns
        -------
        R_c : np.ndarray
            Backscattering correction factor, same shape as `x`.
        """
        U0 = E0 / x
        I_val, G_val = Background_Model._return_IG(U0)
        R_c = 1 - nu_val * (I_val + nu_val * G_val) ** 1.67
        return R_c
    
    
    @staticmethod
    def _return_IG(U0):
        """
        Computes the I(U0) and G(U0) functions of overvoltage ratio needed for backscattering correction,
        as described by:
    
            G. Love, V.D. Scott,
            "Evaluation of a new correction procedure for quantitative electron probe microanalysis",
            J. Phys. D. Appl. Phys. 11 (1978) 1369–1376.
            https://doi.org/10.1088/0022-3727/11/10/002
    
        Parameters
        ----------
        U0 : float or np.ndarray
            Overvoltage ratio (E0 / x).
    
        Returns
        -------
        I_val : float or np.ndarray
            Value of the I(U0) function.
        G_val : float or np.ndarray
            Value of the G(U0) function.
        """
        log_U0 = np.log(U0)
        I_val = (
            0.33148 * log_U0
            + 0.05596 * log_U0**2
            - 0.06339 * log_U0**3
            + 0.00947 * log_U0**4
        )
        G_val = (
            1 / U0
            * (
                2.87898 * log_U0
                - 1.51307 * log_U0**2
                + 0.81312 * log_U0**3
                - 0.08241 * log_U0**4
            )
        )
        return I_val, G_val
    
    
    def get_backscattering_correction_mod_pars(self):
        """
        Returns an lmfit Model and its Parameters for the correction of background signal loss 
        due to electron backscattering.
    
        Returns
        -------
        bs_correction_m : lmfit.Model
            Model describing the correction for loss of generated background signal due to 
            electron backscattering.
        params_bs_cor : lmfit.Parameters
            Parameters for the model.
    
        Notes
        -----
        - The model is based on the implementation in Background_Model._backscattering_correction.
        - The adr_bcksctr parameter controls convolution with the detector response function.
        """
        bs_correction_m = Model(Background_Model._backscattering_correction)
        params_bs_cor = bs_correction_m.make_params(
            adr_bcksctr=dict(expr='apply_det_response')
        )
    
        return bs_correction_m, params_bs_cor

    
    # =============================================================================
    # Generated background model
    # =============================================================================
    @staticmethod
    def _generated_bckgrnd_Castellano2004(
        x, Z, K=0.035, a1=-73.9, a2=-1.2446, a3=36.502, a4=148.5, a5=0.1293, a6=-0.006624, a7=0.0002906, apply_det_response=1
    ):
        """
        Analytical model for the generated bremsstrahlung background spectrum after:
    
            Castellano, G., Osán, J., & Trincavelli, J. (2004).
            "Analytical model for the bremsstrahlung spectrum in the 0.25–20 keV photon energy range."
            Spectrochimica Acta Part B: Atomic Spectroscopy, 59(3), 313–319.
            https://doi.org/10.1016/j.sab.2003.11.008
            
        Z should be weighted on mass fractions as recommended by:
            Trincavelli, J., & Castellano, G. (2008).
            "The prediction of thick target electron bremsstrahlung spectra in the 0.25-50 keV energy range."
            Spectrochimica Acta Part B: Atomic Spectroscopy, 63(1), 1–8.
            https://doi.org/10.1016/j.sab.2007.11.009
            
        Note
        ----
        - This function is NOT USED in the current implementation.
        - The model gives unphysical results for large Z, as the generated background shape becomes concave at low energy
          instead of increasing as expected.
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV).
        Z : int or float
            Atomic number of the element.
        K, a1, a2, a3, a4, a5, a6, a7 : float, optional
            Model parameters (see Castellano et al., 2004).
        apply_det_response : int, optional
            If 1, convolves signal with detector response function. Default is 1.
    
        Returns
        -------
        model : np.ndarray
            Generated background spectrum (arbitrary units), same shape as `x`.
        """
        E0 = Background_Model.cls_beam_e
    
        model = (
            K * np.sqrt(Z) * ((E0 - x) / x)
            * (a1 + a2 * x + a3 * np.log(Z) + a4 * (E0 ** a5) / Z)
            * (1 + (a6 + a7 * E0) * (Z / x))
        )
    
        # Apply convolution when fitting, but not when calculating background for PB Z correction
        if apply_det_response == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        # Clip to avoid negative values (which should not be present)
        model = np.where(model < 0, 0.01, model)  # Avoid division by zero in PB method
    
        return model
    
    
    @staticmethod
    def _generated_bckgrnd_Trincavelli1998(
        x, Z, K=0.035, a1=-54.86, a2=-1.072, a3=0.2835, a4=30.4, a5=875, a6=0.08, apply_det_response=1
    ):
        """
        Analytical model for the generated bremsstrahlung background spectrum after:
    
            Trincavelli, J., Castellano, G., & Riveros, J. A. (1998).
            "Model for the bremsstrahlung spectrum in EPMA. Application to standardless quantification."
            X-Ray Spectrometry, 27(2), 81–86.
            https://doi.org/10.1002/(SICI)1097-4539(199803/04)27:2<81::AID-XRS253>3.0.CO;2-R
        
        Z should be weighted on mass fractions as recommended by:
            Trincavelli, J., & Castellano, G. (2008).
            "The prediction of thick target electron bremsstrahlung spectra in the 0.25-50 keV energy range."
            Spectrochimica Acta Part B: Atomic Spectroscopy, 63(1), 1–8.
            https://doi.org/10.1016/j.sab.2007.11.009
        
        Notes
        -----
        - This function is NOT USED in the current implementation because it is imprecise at low energy values.
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV).
        Z : int or float
            Atomic number of the element.
        K, a1, a2, a3, a4, a5, a6 : float, optional
            Model parameters (see Trincavelli et al., 1998).
        apply_det_response : int, optional
            If 1, convolves signal with detector response function. Default is 1.
    
        Returns
        -------
        model : np.ndarray
            Generated background spectrum (arbitrary units), same shape as `x`.
        """
        E0 = Background_Model.cls_beam_e
    
        model = (
            K * np.sqrt(Z) * (E0 - x) / x
            * (a1 + a2 * x + a3 * E0 + a4 * np.log(Z) + a5 / (E0 ** a6 * Z ** 2))
        )
    
        # Apply convolution when fitting, but not when calculating background for PB Z correction
        if apply_det_response == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        # Clip to avoid negative values (should not be present)
        model = np.where(model < 0, 0.01, model)  # Avoid division by zero in PB method
    
        return model
    
    
    @staticmethod
    def _generated_bckgrnd_DuncumbMod(
        x, Z, K=0.8, F=1, P=1, beta=0, apply_det_response=1
    ):
        """
        Analytical model for the generated bremsstrahlung background spectrum after:
    
            Duncumb, P., Barkshire, I. R., & Statham, P. J. (2001).
            "Improved X-ray Spectrum Simulation for Electron Microprobe Analysis."
            Microscopy and Microanalysis, 7(4), 341–355.
            https://doi.org/10.1007/S10005-001-0010-6
    
        Z should be mass-fraction weighted as per Statham, P. J. (2016).
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV).
        Z : int or float
            Atomic number (should be mass-fraction weighted for compounds/mixtures).
        K : float, optional
            Model scaling parameter. Default is 0.8.
        F : float, optional
            Model parameter (see Duncumb et al., 2001). Default is 1.
        P : float, optional
            Model exponent parameter. Default is 1.
        beta : float, optional
            Model parameter (controls low-energy shape). Default is 0.
        apply_det_response : int, optional
            If 1, convolves signal with detector response function. Default is 1.
    
        Returns
        -------
        model : np.ndarray
            Generated background spectrum (arbitrary units), same shape as `x`.
        """
        E0 = Background_Model.cls_beam_e
    
        model = K * Z * F * ((E0 - x) / x) ** P * (x / (x + beta))
    
        # Apply convolution when fitting, but not when calculating background for PB Z correction
        if apply_det_response == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    
    
    @staticmethod
    def get_beta_expr(beta_expr, els_symbols):
        """
        Returns a symbolic formula string for the beta parameter from the modified Duncumb's continuum generation model.
        beta depends on atomic number Z, mass-fraction-averaged over a compound or mixture.
    
        Parameters
        ----------
        beta_expr : sympy expression
            A sympy expression involving the symbol 'Z', e.g. 0.5 * Z**1.2.
        els_symbols : list of str
            List of element symbols (e.g., ['Si', 'O']).
    
        Returns
        -------
        beta_expr_full : str
            Symbolic formula for the compound parameter, mass-fraction-averaged over elements.
    
        Example
        -------
        If beta_expr = 0.5 * Z**1.2 and els_symbols = ['Si', 'O'], returns:
            "14.054 * f_Si + 4.594 * f_O"
        """
        Z = sp.Symbol('Z')
        els = [Element(el) for el in els_symbols]
    
        if len(els) == 1:
            el = els[0]
            beta = beta_expr.subs(Z, el.Z).evalf()
            beta_expr_full = f"{beta:.3f} * f_{el.symbol}"
        else:
            beta_vals = [beta_expr.subs(Z, el.Z).evalf() for el in els]
            beta_expr_full = " + ".join([f"{coeff:.3f} * f_{el}" for el, coeff in zip(els_symbols, beta_vals)])
    
        return beta_expr_full
    
    
    def get_generated_background_mod_pars(self, fr_pars, is_calibration=False, model='DuncumbMod'):
        """
        Returns an lmfit Model and its Parameters for the continuum X-ray background generated within the sample.
    
        Parameters
        ----------
        fr_pars : dict
            Dictionary of element mass fraction parameters (e.g., {'f_Si': 0.5, 'f_O': 0.5}).
        is_calibration : bool, optional
            If True, use calibration mode with free parameters. Default is False.
        model : str, optional
            Background model to use. Options:
                - 'Castellano2004'
                - 'Trincavelli1998'
                - 'DuncumbMod' (default)
                - 'Duncumb2001'
    
        Returns
        -------
        bckgrnd_m : lmfit.Model
            Model describing the continuum X-ray signal generated within the sample.
        params_bckgrnd : lmfit.Parameters
            Parameters for the model.
    
        Notes
        -----
        - For Castellano2004 and Trincavelli1998, the mass-weighted average Z is used.
        - For DuncumbMod and Duncumb2001, Statham's average Z is used.
        - Duncumb2001 is identical to DuncumbMod but with beta fixed to 0 (not varying).
        - Parameter logic adapts to calibration and particle mode.
        - See references for model details.
    
        References
        ----------
        Castellano, G., Osán, J., & Trincavelli, J. (2004). Analytical model for the bremsstrahlung spectrum in the 0.25–20 keV photon energy range. Spectrochimica Acta Part B: Atomic Spectroscopy, 59(3), 313–319. https://doi.org/10.1016/j.sab.2003.11.008
        Trincavelli, J., Castellano, G., & Riveros, J. A. (1998). Model for the bremsstrahlung spectrum in EPMA. Application to standardless quantification. X-Ray Spectrometry, 27(2), 81–86. https://doi.org/10.1002/(SICI)1097-4539(199803/04)27:2<81::AID-XRS253>3.0.CO;2-R
        Statham, P., Penman, C., & Duncumb, P. (2016). Improved spectrum simulation for validating SEM-EDS analysis. IOP Conf Ser Mater Sci Eng, 109(1):0–10.
        Duncumb, P., Barkshire, I. R., & Statham, P. J. (2001). Improved X-ray Spectrum Simulation for Electron Microprobe Analysis. Microscopy and Microanalysis, 7(4), 341–355. https://doi.org/10.1007/S10005-001-0010-6
    
        """
        els, _, _ = Background_Model._get_els_frs(**fr_pars)
    
        # Model and Z averaging model selection
        model = model.lower()
        if model == 'castellano2004':
            bckgrnd_m = Model(Background_Model._generated_bckgrnd_Castellano2004)
            Z_par_expr = Background_Model.get_average_Z(els, method='mass_weighted')
            params_bckgrnd = bckgrnd_m.make_params(
                Z=dict(expr=Z_par_expr),
                apply_det_response=dict(value=1, vary=False),
                K=dict(value=0.044, vary=False, min=0.04, max=0.10),
                a1=dict(value=-73.9, vary=True, min=-500, max=0),
                a2=dict(value=-1.2446, vary=True, min=-5, max=10),
                a3=dict(value=36.502, vary=False, min=0, max=100),
                a4=dict(value=148.5, vary=True, min=100, max=200),
                a5=dict(value=0.1293, vary=False, min=0, max=1),
                a6=dict(value=-0.006624, vary=False, min=-0.01, max=0.1),
                a7=dict(value=0.0002906, vary=False, min=0, max=0.001),
            )
            return bckgrnd_m, params_bckgrnd
    
        elif model == 'trincavelli1998':
            bckgrnd_m = Model(Background_Model._generated_bckgrnd_Trincavelli1998)
            Z_par_expr = Background_Model.get_average_Z(els, method='mass_weighted')
            params_bckgrnd = bckgrnd_m.make_params(
                apply_det_response=dict(value=1, vary=False),
                Z=dict(expr=Z_par_expr),
                K=dict(value=0.045, vary=True, min=0.001, max=0.2),
                a1=dict(value=-54.86, vary=False, min=-100, max=0),
                a2=dict(value=-1.072, vary=False, min=-5, max=5),
                a3=dict(value=0.2835, vary=False, min=0, max=100),
                a4=dict(value=30.4, vary=False, min=100, max=200),
                a5=dict(value=875, vary=False, min=0, max=1),
                a6=dict(value=0.08, vary=False),
            )
            return bckgrnd_m, params_bckgrnd
    
        elif model in ['duncumbmod', 'duncumb2001']:
            bckgrnd_m = Model(
                Background_Model._generated_bckgrnd_DuncumbMod if model == 'duncumbmod'
                else Background_Model._generated_bckgrnd_Duncumb
            )
            Z_par_expr = Background_Model.get_average_Z(els, method='Statham')
    
            # Determine scaling factor K
            if self.sp_collection_time is not None and self.sp_collection_time > 0:
                K_val = self.sp_collection_time * calibs.gen_background_time_scaling_factor[self.meas_mode]
            else:
                K_val = self.tot_sp_counts / 1e5  # Calibrated fallback in case collection time is not recorded (should not be the case)
    
            # Retrieve calibrated model parameter expressions
            P_expr, F_expr, beta_expr_Z = calibs.get_calibrated_background_params(self.meas_mode)
    
            # Set up parameters depending on calibration/particle mode
            if is_calibration:
                P_par = dict(value=1.16, vary=True, min=1, max=1.3)
                F_par = dict(value=1, vary=True, min=0.1, max=1.3)
                beta_param = {'value': 0.2, 'vary': True, 'min': 0, 'max': 0.5} if model == 'duncumbmod' else {'value': 0, 'vary': False}
            elif self.is_particle:
                P_par = dict(expr=P_expr)
                F_par = dict(expr=F_expr)
                beta_expr = Background_Model.get_beta_expr(beta_expr_Z, els)
                beta_param = {'expr': beta_expr} if model == 'duncumbmod' else {'value': 0, 'vary': False}
            else:
                P_par = dict(expr=P_expr)
                F_par = dict(expr=F_expr)
                beta_param = {'value': 0.2, 'vary': True, 'min': 0, 'max': 0.5} if model == 'duncumbmod' else {'value': 0, 'vary': False}
    
            params_bckgrnd = bckgrnd_m.make_params(
                apply_det_response=dict(value=1, vary=False),
                Z=dict(expr=Z_par_expr),
                K=dict(value=K_val, vary=True, min=0.001, max=np.inf),
                P=P_par,
                F=F_par,
                beta=beta_param
            )
            return bckgrnd_m, params_bckgrnd
    
        else:
            raise ValueError(
                f"Unknown background model '{model}'. "
                "Choose from 'Castellano2004', 'Trincavelli1998', 'DuncumbMod', or 'Duncumb2001'."
            )
    
    
    # =============================================================================
    # Detector efficiency and zero strobe peak
    # =============================================================================
    @staticmethod
    def _det_efficiency(x, adr_det_eff=1):
        """
        Returns the detector efficiency as a function of energy by interpolating
        the detector efficiency calibration curve.
    
        Parameters
        ----------
        x : array-like
            Photon energies (keV) at which to evaluate the detector efficiency.
        adr_det_eff : int, optional
            If 1, convolves signal with detector response function. Default is 1.
    
        Returns
        -------
        model : np.ndarray
            Detector efficiency at each energy in `x`.
        """
        model = np.interp(
            x,
            DetectorResponseFunction.det_eff_energy_vals,
            DetectorResponseFunction.det_eff_vals
        )
    
        if adr_det_eff == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    

    def get_detector_efficiency_mod_pars(self):
        """
        Returns an lmfit Model and its Parameters for the EDS detector efficiency correction.
    
        Returns
        -------
        detector_efficiency_m : lmfit.Model
            Model describing the loss of signal due to the EDS detector efficiency.
            This model has no adjustable parameters, as we use the efficiency provided by the EDS manufacturer.
        detector_efficiency_pars : lmfit.Parameters
            Parameters for the model.
    
        Notes
        -----
        - The model is based on Background_Model._det_efficiency.
        - The adr_det_eff parameter controls convolution with the detector response function.
        """
        detector_efficiency_m = Model(Background_Model._det_efficiency)
        detector_efficiency_pars = detector_efficiency_m.make_params(
            adr_det_eff=dict(expr='apply_det_response')
        )
    
        return detector_efficiency_m, detector_efficiency_pars
    
    
    def get_det_zero_peak_model_pars(self, amplitude_val):
        """
        Returns an lmfit GaussianModel and its Parameters for modeling the detector zero (strobe) peak.
    
        Parameters
        ----------
        amplitude_val : float
            Initial amplitude estimate for the zero peak.
    
        Returns
        -------
        det_zero_peak_model : lmfit.models.GaussianModel
            Gaussian model for the detector zero (strobe) peak.
        params_det_zero_peak : lmfit.Parameters
            Parameters for the model.
    
        Notes
        -----
        - The amplitude limits are set wider for particles than for bulk samples.
        - The sigma (width) of the zero peak is taken from calibration for the current EDS mode.
        """
        if self.is_particle:
            min_ampl = amplitude_val / 3
            max_ampl = amplitude_val * 5
        else:
            min_ampl = amplitude_val / 2
            max_ampl = amplitude_val * 3
    
        det_zero_peak_model = GaussianModel(prefix='det_zero_peak_')
        zero_strobe_peak_sigma = calibs.zero_strobe_peak_sigma[self.meas_mode]
        params_det_zero_peak = det_zero_peak_model.make_params(
            amplitude=dict(value=amplitude_val, vary=True, min=min_ampl, max=max_ampl),
            center=dict(value=0, vary=False),
            sigma=dict(value=zero_strobe_peak_sigma, vary=False, min=0.05, max=1)
        )
    
        return det_zero_peak_model, params_det_zero_peak
    
    
    # =============================================================================
    # Full background model construction
    # =============================================================================
    def get_full_background_mod_pars(self, fr_pars):
        """
        Constructs the full background model and its parameters for spectral fitting.
    
        Parameters
        ----------
        fr_pars : dict
            Dictionary of element mass fraction parameters (e.g., {'f_Si': 0.5, 'f_O': 0.5}).
    
        Returns
        -------
        background_mod : lmfit.Model
            The full composite model for the background spectrum, including all physical corrections and
            the detector zero (strobe) peak.
        background_pars : lmfit.Parameters
            Combined parameters for the full background model.
    
        Notes
        -----
        - The model includes generated background, absorption attenuation, detector efficiency,
          backscattering correction, stopping power correction, and the detector strobe peak.
        - All sub-models and their parameters are constructed and combined automatically.
        """
    
        # Generated background, model and parameters
        gen_bckgrnd_mod, gen_bckgrnd_pars = self.get_generated_background_mod_pars(fr_pars, model='DuncumbMod')
    
        # Attenuation due to absorption, model and parameters
        abs_att_mod, abs_att_pars = self.get_abs_attenuation_mod_pars(model='phirho')
    
        # Backscattering correction, model and parameters
        bs_cor_mod, bs_cor_pars = self.get_backscattering_correction_mod_pars()
    
        # Stopping power correction, model and parameters
        stopping_p_mod, stopping_p_pars = self.get_stopping_power_mod_pars()
        
        # EDS detector efficiency, model and parameters
        det_eff_mod, det_eff_pars = self.get_detector_efficiency_mod_pars()
    
        # Add detector strobe (zero) peak to the background model
        if self.sp_collection_time is not None and self.sp_collection_time > 0:
            amplitude_val = self.sp_collection_time * calibs.strobe_peak_int_factor[self.meas_mode]
        else:
            amplitude_val = self.tot_sp_counts / (10**4) # In case self.sp_collection_time = None
        det_zero_peak_mod, det_zero_peak_par = self.get_det_zero_peak_model_pars(amplitude_val)
    
        # Full background model construction
        background_mod = (
            gen_bckgrnd_mod
            * abs_att_mod
            * det_eff_mod
            * bs_cor_mod
            * stopping_p_mod
            + det_zero_peak_mod
        )
    
        # Combine all parameters
        background_pars = Parameters()
        background_pars.update(gen_bckgrnd_pars)
        background_pars.update(abs_att_pars)
        background_pars.update(det_eff_pars)
        background_pars.update(bs_cor_pars)
        background_pars.update(stopping_p_pars)
        background_pars.update(det_zero_peak_par)
    
        return background_mod, background_pars
    

#%% DetectorResponseFunction class
class DetectorResponseFunction():
    det_res_conv_matrix = None
    icc_conv_matrix = None
    energy_vals_padding = 30 # Padding added to energy_vals to ensure correct functioning of convolution operation

    @classmethod
    def setup_detector_response_vars(cls, det_ch_offset, det_ch_width, spectrum_lims, microscope_ID, verbose=True):
        """
        Initialize detector response variables for the EDS system.
    
        Loads detector efficiency and convolution matrices for the specified microscope.
        If convolution matrices for the given channel settings do not exist, they are calculated and saved.
    
        Parameters
        ----------
        det_ch_offset : float
            Energy offset for detector channels (in keV).
        det_ch_width : float
            Channel width (in keV).
        spectrum_lims : tuple of int
            (low, high) indices for the usable spectrum region.
        microscope_ID : str
            Identifier for the microscope/calibration directory.
        verbose : bool, optional
            If True, print status messages.
    
        Sets Class Attributes
        ---------------------
        det_eff_energy_vals : np.ndarray
            Energy values for detector efficiency.
        det_eff_vals : np.ndarray
            Detector efficiency values.
        det_res_conv_matrix : np.ndarray
            Detector response convolution matrix (cropped to spectrum_lims).
        icc_conv_matrix : np.ndarray
            ICC convolution matrix (cropped to spectrum_lims).
    
        Notes
        -----
        - Convolution matrices are detector-dependent and cached in JSON for efficiency.
        - If multiple EDS detectors are used, this code should be adapted.
        """
    
        # --- Load EDS detector efficiency spectrum ---
        detector_efficiency_path = os.path.join(
            parent_dir, cnst.XRAY_SPECTRA_CALIBS_DIR, cnst.MICROSCOPES_CALIBS_DIR, microscope_ID, cnst.DETECTOR_EFFICIENCY_FILENAME
        )
        det_eff_energy_vals, det_eff_vals, metadata = load_msa(detector_efficiency_path)
        if metadata['XUNITS'] == 'eV':
            det_eff_energy_vals /= 1000  # Convert to keV
    
        cls.det_eff_energy_vals = det_eff_energy_vals
        cls.det_eff_vals = det_eff_vals
    
        # --- Load or calculate convolution matrices ---
        
        conv_matrices_file_path = os.path.join(
            parent_dir, cnst.XRAY_SPECTRA_CALIBS_DIR, cnst.MICROSCOPES_CALIBS_DIR, microscope_ID, cnst.DETECTOR_CONV_MATRICES_FILENAME
        )
        if os.path.exists(conv_matrices_file_path):
            with open(conv_matrices_file_path, 'r') as file:
                conv_matrices_dict = json.load(file)
        else:
            conv_matrices_dict = {}
    
        # Key for current detector channel settings
        conv_mat_key = f"O{det_ch_offset},W{det_ch_width}"
        conv_matrices = conv_matrices_dict.get(conv_mat_key)
    
        if conv_matrices is None:
            # Generate full energy vector for all detector channels
            full_en_vector = [
                det_ch_offset + j * det_ch_width for j in range(calibs.detector_ch_n)
            ]
    
            # Calculate and save convolution matrices
            det_res_conv_matrix = cls._calc_det_res_conv_matrix(full_en_vector)
            icc_conv_matrix = cls._calc_icc_conv_matrix(full_en_vector)
    
            # Store as lists for JSON serialization
            conv_matrices_dict[conv_mat_key] = (
                det_res_conv_matrix.tolist(), icc_conv_matrix.tolist()
            )
            with open(conv_matrices_file_path, 'w') as file:
                json.dump(conv_matrices_dict, file)
        else:
            if verbose:
                print_single_separator()
                print("Detector response convolution matrices loaded")
            det_res_conv_matrix, icc_conv_matrix = conv_matrices
    
        # --- Crop matrices to match spectrum limits ---
        low_l = spectrum_lims[0] + 1
        high_l = spectrum_lims[1] + cls.energy_vals_padding // 2 + cls.energy_vals_padding - 1
        cls.det_res_conv_matrix = np.array(det_res_conv_matrix)[low_l:high_l, low_l:high_l]
        cls.icc_conv_matrix = np.array(icc_conv_matrix)[low_l:high_l, low_l:high_l]

    # =============================================================================
    # Convolution of signal with detector response function
    # =============================================================================
    @classmethod
    def _apply_padding_with_fit(cls, signal):
        """
        Pad a signal at both ends by linear extrapolation, using a fit to the first and last few points.
    
        This method is used to extend the signal array, avoiding edge effects in convolution.
        Padding values are clipped to be non-negative.
    
        Parameters
        ----------
        signal : np.ndarray or list
            The input 1D signal to pad.
    
        Returns
        -------
        padded_signal : np.ndarray
            The input signal with linear-extrapolated padding added at both ends.
    
        Notes
        -----
        - The number of padding points is determined by cls.energy_vals_padding.
        - Linear fits are performed on the first and last 4 points of the signal.
        - Extrapolated values are clipped at zero to avoid negative padding.
        """
        n_pts_fitted = 4  # Number of points used for linear fit at each end
    
        # --- Padding at the beginning (head) ---
        x_head = signal[:n_pts_fitted]
        x_indices_head = np.arange(len(x_head))
        # Linear fit to the first few points
        slope_head, intercept_head = np.polyfit(x_indices_head, x_head, 1)
        # Indices for extrapolation (negative indices for padding before signal)
        extrapolation_indices_head = np.arange(-cls.energy_vals_padding // 2 + 1, 0)
        # Linear extrapolation for padding
        extrapolated_values_head = slope_head * extrapolation_indices_head + intercept_head
        # Ensure no negative values in padding
        extrapolated_values_head = np.clip(extrapolated_values_head, 0, None)
    
        # --- Padding at the end (tail) ---
        x_tail = signal[-n_pts_fitted:]
        x_indices_tail = np.arange(len(x_tail))
        # Linear fit to the last few points
        slope_tail, intercept_tail = np.polyfit(x_indices_tail, x_tail, 1)
        # Indices for extrapolation (beyond signal end)
        extrapolation_indices_tail = np.arange(len(x_tail) + 1, len(x_tail) + cls.energy_vals_padding)
        # Linear extrapolation for padding
        extrapolated_values_tail = slope_tail * extrapolation_indices_tail + intercept_tail
        # Ensure no negative values in padding
        extrapolated_values_tail = np.clip(extrapolated_values_tail, 0, None)
    
        # --- Combine padding and original signal ---
        padded_signal = np.concatenate([extrapolated_values_head, signal, extrapolated_values_tail])
    
        return padded_signal
    
    
    @classmethod
    def _apply_det_response_fncts(cls, signal):
        """
        Apply both detector resolution and ICC convolution functions to a signal,
        with edge padding by linear extrapolation.
    
        The signal is padded at both ends using a linear fit, then convolved sequentially
        with the detector resolution and ICC convolution matrices. The padding is removed
        from the final result.
    
        Parameters
        ----------
        signal : np.ndarray or list
            The input 1D signal to be convolved.
    
        Returns
        -------
        processed_model : np.ndarray
            The signal after both convolutions, trimmed to exclude the padding.
    
        Notes
        -----
        - Padding is performed by fitting and extrapolating the first and last few points.
        - The convolution matrices must be initialized in the class before use.
        - This approach is more accurate than simple replicative padding, but slightly slower.
        """
        # Pad the signal at both ends using linear fit/extrapolation
        padded_signal = cls._apply_padding_with_fit(signal)
    
        # First, convolve with the detector resolution matrix
        model = np.sum(cls.det_res_conv_matrix * padded_signal, axis=1)
    
        # Then, convolve with the ICC convolution matrix
        model = np.sum(cls.icc_conv_matrix * model, axis=1)
    
        # Remove padding from the result to return only the original signal region
        processed_model = model[cls.energy_vals_padding // 2 - 1 : -cls.energy_vals_padding + 1]
    
        return processed_model
    
    # =============================================================================
    # Detector resolution convolution
    # =============================================================================
    @staticmethod
    def _det_sigma(E):
        """
        Calculate the detector Gaussian sigma for a given X-ray energy.
        
        Requires calibration file to be correctly loaded via setup_detector_response_vars.
    
        Based on:
        N.W.M. Ritchie, "Spectrum Simulation in DTSA-II", Microsc. Microanal. 15 (2009) 454–468.
        https://doi.org/10.1017/S1431927609990407
    
        Parameters
        ----------
        E : float or array-like
            X-ray energy in keV.
    
        Returns
        -------
        sigma : float or np.ndarray
            Standard deviation (sigma) of the detector response at energy E.
    
        Notes
        -----
        - Parameters (conv_eff, elec_noise, F) are calibrated on several elements in bulk standards.
        - See Calculation_pars_peak_fwhm.py for details.
        """
        # Calculate detector sigma using calibrated parameters
        sigma = calibs.conv_eff * np.sqrt(calibs.elec_noise**2 + E * calibs.F / calibs.conv_eff)
        return sigma
        
    
    @classmethod
    def _calc_det_res_conv_matrix(cls, energy_vals, verbose = True):
        """
        Calculate the detector resolution convolution matrix.
    
        Each row of the matrix represents the probability distribution (Gaussian)
        for an input energy, accounting for the detector's energy resolution.
        Padding is added to minimize edge effects during convolution.
    
        Parameters
        ----------
        energy_vals : array-like
            Array of energy values (in keV) for which to compute the convolution matrix.
        verbose : bool, optional
            If True, print status messages.
    
        Returns
        -------
        det_res_conv_matrix : np.ndarray
            The detector resolution convolution matrix (energy_vals x energy_vals).
    
        Notes
        -----
        - Padding is applied to account for convolution spillover at the spectrum edges.
        - The Gaussian sigma is energy-dependent and calculated with _det_sigma().
        - Integration is performed for each matrix element to ensure normalization.
        """
    
        if verbose:
            start_time = time.time()
            print("Calculating convolution matrix for detector resolution")
    
        deltaE = energy_vals[5] - energy_vals[4]
        n_intervals = cls.energy_vals_padding
    
        # Extend the energy axis with padding on both sides to avoid edge effects
        left_pad = [energy_vals[0] - deltaE * i for i in range(n_intervals // 2, 0, -1)]
        right_pad = [energy_vals[-1] + deltaE * i for i in range(1, n_intervals)]
        energy_vals_extended = left_pad + list(energy_vals) + right_pad
    
        def gaussian(E, E0, sigma):
            """Normalized Gaussian function."""
            return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 / sigma**2 * (E - E0)**2)
    
        conv_matrix = []
        for i, en in enumerate(energy_vals_extended):
            # Initialize row for this energy; padding prevents signal loss at edges
            g_vals = [0.0 for _ in range(len(energy_vals_extended))]
            sigma = cls._det_sigma(en)
            for j in range(-n_intervals, n_intervals + 1):
                cen_E = en + j * deltaE
                idx = i + j
                if 0 <= idx < len(g_vals):
                    try:
                        # Integrate the Gaussian over the width of the energy bin
                        int_E, _ = quad(lambda E: gaussian(E, en, sigma), cen_E - deltaE / 2, cen_E + deltaE / 2)
                        g_vals[idx] = int_E
                    except Exception:
                        pass  # Ignore integration errors (should be rare)
            conv_matrix.append(g_vals)
    
        det_res_conv_matrix = np.array(conv_matrix).T  # Transpose for correct orientation
        
        if verbose:
            process_time = time.time() - start_time
            print(f"Calculation executed in {process_time:.1f} s")
    
        return det_res_conv_matrix
    
    
    # =============================================================================
    # Incomplete Charge Collection Spectrum Calculations
    # =============================================================================
    @staticmethod
    def get_icc_spectrum(energy_vals, line_en, R_e=50e-7, F_loss=0.27):
        """
        Generate the ICC (Incomplete Charge Collection) smearing function for a given line energy.
    
        Parameters
        ----------
        energy_vals : array-like
            Energy axis (in keV) of detector channels.
        line_en : float
            X-ray line energy (in keV).
        R_e : float, optional
            Effective recombination parameter (cm).
        F_loss : float, optional
            Fractional charge loss parameter.
    
        Returns
        -------
        icc_n_vals_distr : list
            ICC distribution, mapped to detector channel energies.
        """
        icc_e_vals, icc_n_vals = DetectorResponseFunction._icc_fnct(line_en, R_e, F_loss)
        icc_n_vals_distr = DetectorResponseFunction._distribute_icc_over_EDS_channels(
            energy_vals, icc_e_vals, icc_n_vals
        )
        return icc_n_vals_distr
    
    
    @staticmethod
    def _icc_fnct(line_en, R_e=50e-7, F_loss=0.27):
        """
        Calculate the ICC function n(E) for a given X-ray line energy.
        
        ICC model as described in:
        Redus, R. H., & Huber, A. C. (2015). Response Function of Silicon Drift Detectors for Low Energy X-rays.
        In Advances in X-ray Analysis (AXA) (pp. 274–282). International Centre for Diffraction Data (ICDD).
    
        Parameters
        ----------
        line_en : float
            X-ray line energy (in keV).
        R_e : float, optional
            Effective recombination parameter (cm).
        F_loss : float, optional
            Fractional charge loss parameter.
    
        Returns
        -------
        e_vals : list of float
            Energy values (in keV) for the ICC function.
        n_vals : list of float
            ICC function values at those energies.
        """
        # Absorption coefficient of Si at energy line_en
        alpha = xray_mass_absorption_coeff(element='Si', energies=line_en) * calibs.Si_density  # cm^-1
    
        V_tot = 4 * np.pi / 3 * (R_e) ** 3
    
        def V_1(z):
            return np.pi / 3 * (R_e - z) ** 2 * (2 * R_e + z)
    
        def Q(z):
            return (V_tot - F_loss * V_1(z)) / V_tot
    
        def dQ_dz(z):
            return -np.pi * F_loss / V_tot * (z ** 2 - R_e ** 2)
    
        def N(z):
            return 1 - np.exp(-alpha * z)
    
        def dN_dz(z):
            return alpha * np.exp(-alpha * z)
    
        def get_z(Q_val):
            Q_val_rnd = np.clip(Q_val, Q_min, 1)
            solution = root_scalar(lambda z: Q(z) - Q_val_rnd, method='brentq', bracket=[0, R_e])
            return solution.root
    
        def n(x):
            Q_val = x / line_en
            z_val = get_z(Q_val)
            n_val = dN_dz(z_val) * dQ_dz(z_val) ** -1 / line_en
            return n_val
    
        def get_n_at_line_en(integral_rest, last_E, last_n):
            def n_fnct(n_):
                n_ = np.float64(n_)
                res = np.trapz([last_n, n_], [last_E, line_en])
                return res
            guess = (1 - integral_rest) / (line_en - last_E)
            guess_2 = guess * 2
            solution = root_scalar(lambda n_: n_fnct(n_) - (1 - integral_rest), x0=guess, x1=guess_2, method='secant')
            return solution.root
    
        # Calculate left boundary of ICC smearing function
        Q_min = Q(0)
        E_min = line_en * Q_min
    
        e_vals = list(np.linspace(E_min, line_en, 1000))
        e_vals.pop()  # Remove last energy value corresponding to line_en
        n_vals = [n(en) for en in e_vals]
        signal_integral = trapezoid(n_vals, e_vals)
        n_val_at_E = get_n_at_line_en(signal_integral, e_vals[-1], n_vals[-1])
    
        # Update lists with values at line_en
        e_vals.append(line_en)
        n_vals.append(n_val_at_E)
    
        return e_vals, n_vals


    @staticmethod
    def _distribute_icc_over_EDS_channels(eds_en_vals, icc_en_vals, icc_n_vals):
        """
        Distribute the ICC function over EDS detector channels.
    
        Parameters
        ----------
        eds_en_vals : array-like
            Detector channel energy values (in keV).
        icc_en_vals : list
            ICC function energy values (in keV).
        icc_n_vals : list
            ICC function values.
    
        Returns
        -------
        eds_icc_n_vals : list
            ICC values distributed over the detector channels.
        """
        ch_width = eds_en_vals[1] - eds_en_vals[0]
        icc_en_spacing = icc_en_vals[1] - icc_en_vals[0]
    
        # Determine which channels are affected by ICC
        indices_affected = [
            i for i, en in enumerate(eds_en_vals)
            if icc_en_vals[0] - ch_width / 2 < en <= icc_en_vals[-1] + ch_width / 2
        ]
        # Calculate number of points to add on the right side of the list to make it symmetrical. Needed to avoid shifts during convolution
        n_pts_to_center_data = len(indices_affected) - 1
        n_pts_added = 20  # Pad array of energy values for full overlap during convolution
        
        # Check if reference energy is outside the range of energy values. This can happen with characteristic X-rays outside the energy range
        if len(indices_affected) == 0:
            return [1] # ICC convolution is not applied if peak is outside energy range
    
        indices_affected = (
            list(range(indices_affected[0] - n_pts_added, indices_affected[0])) +
            indices_affected +
            list(range(indices_affected[-1] + 1, indices_affected[-1] + n_pts_added + n_pts_to_center_data + 1))
        )
    
        # Remove negative indices, maintaining symmetry
        if indices_affected[0] < 0:
            index_zero = indices_affected.index(0)
            indices_affected = indices_affected[index_zero: len(indices_affected) - index_zero]
    
        # Remove indices beyond dimension of eds_en_vals if needed
        len_eds_en_vals = len(eds_en_vals)
        if indices_affected[-1] >= len_eds_en_vals:
            index_last = indices_affected.index(len_eds_en_vals)
            n_pts_to_remove = len(indices_affected) - index_last + 1
            indices_affected = indices_affected[n_pts_to_remove: len(indices_affected) - n_pts_to_remove]
    
        # Distribute ICC function over the detector channels
        eds_icc_en_vals = [eds_en_vals[i] for i in indices_affected]
        eds_icc_n_vals = []
        for index, en in enumerate(eds_icc_en_vals):
            interval_boundary_left = en - ch_width / 2
            interval_boundary_right = en + ch_width / 2
            indices_to_int = [
                i for i, e in enumerate(icc_en_vals) if interval_boundary_left < e <= interval_boundary_right
            ]
            e_vals_to_int = [icc_en_vals[i] for i in indices_to_int]
            n_vals_to_int = [icc_n_vals[i] for i in indices_to_int]
            if indices_to_int:
                if len(indices_to_int) > 1:
                    # Integrate ICC function over interval corresponding to energy value en
                    eds_icc_n_val = trapezoid(n_vals_to_int, e_vals_to_int)
                else: # Case of only 1 point within the detector channel
                    # There is no full interval of en_spacing width within the current detector channel
                    # The portion of interval within this channel is added on the next steps
                    eds_icc_n_val = 0
                # Add portion of interval shared with left of en, unless at boundary
                if interval_boundary_left < e_vals_to_int[0] and e_vals_to_int[0] > 0 and indices_to_int[0] != 0:
                    extra_i_left = indices_to_int[0] - 1
                    left_int = trapezoid([icc_n_vals[extra_i_left], n_vals_to_int[0]],
                                         [icc_en_vals[extra_i_left], e_vals_to_int[0]])
                    eds_icc_n_val += left_int * (e_vals_to_int[0] - interval_boundary_left) / icc_en_spacing
                # Add portion of interval shared with right of en, unless at boundary
                if interval_boundary_right > e_vals_to_int[-1] and indices_to_int[-1] != len(icc_n_vals) - 1:
                    extra_i_right = indices_to_int[-1] + 1
                    right_int = trapezoid([n_vals_to_int[-1], icc_n_vals[extra_i_right]],
                                          [e_vals_to_int[-1], icc_en_vals[extra_i_right]])
                    eds_icc_n_val += right_int * (interval_boundary_right - e_vals_to_int[-1]) / icc_en_spacing
            else:
                eds_icc_n_val = 0
            eds_icc_n_vals.append(eds_icc_n_val)
    
        return eds_icc_n_vals


    @classmethod
    def _calc_icc_conv_matrix(cls, energy_vals, verbose = True):
        """
        Calculate the ICC convolution matrix for all detector channels.
    
        Parameters
        ----------
        energy_vals : array-like
            Array of energy values (in keV) for which to compute the convolution matrix.
        verbose : bool, optional
            If True, print status messages.
            
        Returns
        -------
        icc_conv_matrix : np.ndarray
            The ICC convolution matrix (energy_vals x energy_vals).
        """
        if verbose:
            start_time = time.time()
            print("Calculating convolution matrix for incomplete charge collection")
    
        deltaE = energy_vals[5] - energy_vals[4]
        n_intervals = cls.energy_vals_padding
    
        # Extend the energy axis with padding on both sides to avoid edge effects
        left_pad = [energy_vals[0] - deltaE * i for i in range(n_intervals // 2, 0, -1)]
        right_pad = [energy_vals[-1] + deltaE * i for i in range(1, n_intervals)]
        energy_vals_extended = left_pad + list(energy_vals) + right_pad
    
        conv_matrix = []
        len_row = len(energy_vals_extended)
        for i, en in enumerate(energy_vals_extended):
            if verbose:
                print(f'{i}\tEnergy: {en * 1000:.1f} eV')
            icc_n_vals = np.zeros([len_row])
            if en > 0:
                icc_spec = DetectorResponseFunction.get_icc_spectrum(
                    energy_vals_extended, en, calibs.R_e_background, calibs.F_loss_background
                )
                if len(icc_spec) == 0:
                    icc_spec = [0]
                icc_n_vals[i] = 1
                icc_n_vals = np.convolve(icc_n_vals, icc_spec, mode='same')
            conv_matrix.append(icc_n_vals)
    
        icc_conv_matrix = np.array(conv_matrix).T
        
        if verbose:
            process_time = time.time() - start_time
            print(f"Calculation executed in {process_time:.1f} s")
    
        return icc_conv_matrix