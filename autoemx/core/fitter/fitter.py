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
import sys
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
from pymatgen.core import Element

# =============================================================================
# lmfit import and patching
# =============================================================================
# lmfit does not support full_output=False to prevent calculation of uncertainties
# To make fits considerably faster, we patch lmfit to prevent uncertainty calculation

from lmfit.minimizer import Minimizer

def patch_lmfit_fast_mode(verbose = False):
    """Disable all covariance/uncertainty computations globally in lmfit."""
    if getattr(Minimizer, "_fastmode_patched", False):
        return

    patched_something = False
    
    warnings.filterwarnings("ignore", category=UserWarning, module="uncertainties")
    
    if hasattr(Minimizer, "_calculate_uncertainties_correlations"):
        def dummy_uncertainties(self):
            if hasattr(self, "result"):
                res = self.result
                res.errorbars = False
                res.uvars = None
                res.covar = None
                for p in res.params.values():
                    p.stderr = None
                    p.correl = None
            return None
        Minimizer._calculate_uncertainties_correlations = dummy_uncertainties
        patched_something = True

    elif hasattr(Minimizer, "_calculate_uncertainties"):
        def dummy_uncertainties(self):
            if hasattr(self, "result"):
                res = self.result
                res.errorbars = False
                res.uvars = None
                res.covar = None
                for p in res.params.values():
                    p.stderr = None
                    p.correl = None
            return None
        Minimizer._calculate_uncertainties = dummy_uncertainties
        patched_something = True
    else:
        warnings.warn("⚠️ lmfit fast mode patch could not find uncertainty calculation method.")

    if hasattr(Minimizer, "_int2ext_cov_x"):
        Minimizer._int2ext_cov_x = lambda self, cov_int, fvars: cov_int
    else:
        warnings.warn("⚠️ Covariance transform method not found in Minimizer.")

    Minimizer._fastmode_patched = True

    if patched_something and verbose:
        logger.info("✅ lmfit patched for speed: uncertainties/covariance will NOT be calculated")

patch_lmfit_fast_mode()

from lmfit import Model, Parameters, Parameter
from lmfit.models import GaussianModel

# =============================================================================
# Package imports
# =============================================================================
from autoemx.utils.helper import (
    RefLineError, print_single_separator, print_double_separator,
    weight_to_atomic_fr, load_msa
)
import autoemx.utils.constants as cnst
import autoemx.calibrations as calibs 
from autoemx.data.Xray_lines import get_el_xray_lines
from autoemx.data.Xray_absorption_coeffs import xray_mass_absorption_coeff
from autoemx.data.mean_ionization_potentials import J_df

# Import fitting submodule components
from .detector_response import DetectorResponseFunction
from .peaks import Peaks_Model
from .background import Background_Model

from autoemx._logging import get_logger
logger = get_logger(__name__)

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
        verbose=False,
        free_area_el_lines=None,
    ):
        """Initialize the EDS spectrum fitter."""
        # Handle mutable default arguments
        if els_to_quantify is None:
            els_to_quantify = []
        if els_substrate is None:
            els_substrate = ['C', 'O', 'Al']
        if els_w_fr is None:
            els_w_fr = {}
        self.free_area_el_lines = free_area_el_lines
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
            det_ch_offset, det_ch_width, spectrum_lims, microscope_ID, verbose=True
        )

        self.print_evolving_params = print_evolving_params
        self.verbose = verbose
        
        
        
    def _define_xray_lines(self):
        """Defines the list of elemental X-ray lines to be fitted."""
        min_overvoltage = max(self.beam_energy / self.energy_vals[-1] * 0.99, 1.2)
        min_energy = self.energy_vals[0]
        max_energy = self.energy_vals[-1]
    
        peak_en_threshold = min_energy - 3 * DetectorResponseFunction._det_sigma(min_energy)
    
        el_lines_list = []
        el_lines_weight_refs_dict = {}
    
        for el in self.els_to_fit_list:
            el_xRays_dict = get_el_xray_lines(el)

            for xray_line, xray_info in el_xRays_dict.items():
                line_en = xray_info['energy (keV)']
    
                if self.beam_energy / line_en > min_overvoltage and line_en > peak_en_threshold:
                    xRay_line_str = f"{el}_{xray_line}"
                    el_lines_list.append(xRay_line_str)
                    el_ref_line = self._get_reference_xray_line(el, xray_line, el_xRays_dict)
                    el_lines_weight_refs_dict[xRay_line_str] = el_ref_line
    
                    if xray_info['weight'] > 0.3 and line_en > 1.74:
                        escape_peak_str = xRay_line_str + self.escape_peaks_str
                        el_lines_list.append(escape_peak_str)
                        el_lines_weight_refs_dict[escape_peak_str] = el_ref_line
    
                    if xray_info['weight'] > 0.3 and 2 * line_en < max_energy:
                        pileup_peak_str = xRay_line_str + self.pileup_peaks_str
                        el_lines_list.append(pileup_peak_str)
                        el_lines_weight_refs_dict[pileup_peak_str] = el_ref_line
    
        ref_lines = [el_line for el_line in el_lines_list if any(ref in el_line for ref in self.xray_quant_ref_lines)]
        other_lines = list(set(el_lines_list) - set(ref_lines))
    
        self.el_lines_list = ref_lines + other_lines
        self.el_lines_weight_refs_dict = el_lines_weight_refs_dict
        


    def _get_reference_xray_line(self, el, line, el_xRays_dict):
        """Determines the appropriate reference X-ray line for a given characteristic line."""
        if line[0] == 'N':
            ref_line_start = 'M'
        else:
            ref_line_start = line[0]
    
        ref_line_l = [ref_line for ref_line in self.xray_quant_ref_lines if ref_line_start == ref_line[0]]
    
        el_line = f"{el}_{line}"
        if len(ref_line_l) == 0:
            raise RefLineError(f"K, L or M references not found for {el_line} line.")
        elif len(ref_line_l) > 1 and ref_line_start in ['K', 'L']:
            raise RefLineError(f"Multiple reference lines found for {el_line}. Only one should be present.")
        elif ref_line_start == 'M':
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
                ref_line = ref_line_l[0]
            else:
                ref_line = 'Ll'
    
        el_ref_line = f"{el}_{ref_line}"
        return el_ref_line



    def _get_fraction_pars(self, elements):
        """Create an lmfit Parameters object for elemental mass fractions."""
        fr_params = Parameters()
    
        if self.els_w_fr:
            trace_els = [el for el in elements if el not in self.els_w_fr]
            
            if len(trace_els) > 0:
                elements = list(set(trace_els) | set(self.els_w_fr))
            else:
                elements = list(self.els_w_fr)
    
        num_els_to_fit = len(elements)
    
        total_fr = 0
        last_el_fr_expr = '1'
    
        for i, el in enumerate(elements):
            par_name = 'f_' + el
            if self.els_w_fr:
                if el in trace_els:
                    if self.verbose:
                        logger.warning(f"⚠️ No fraction was assigned to element {el}. Assuming trace element.")
                    val = (1 - sum(self.els_w_fr[el] for el in elements if el in self.els_w_fr)) / len(trace_els)
                else:
                    val = self.els_w_fr[el]
                fr_params.add(par_name, value=val, vary=False)
            else:
                if self.force_fr_total and i == num_els_to_fit - 1:
                    fr_params.add('f_' + elements[-1], expr=last_el_fr_expr, min=0, max=1)
                else:
                    w_fr = 1 / num_els_to_fit

                    last_el_fr_expr += '-' + par_name
                    total_fr += w_fr

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
        """Instantiate and initialize the background model for the current spectrum."""
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
        """Generate the background lmfit model and its parameters."""
        bckgrnd_model_and_pars = self._initialise_Background_Model()
        background_mod, background_pars = bckgrnd_model_and_pars.get_full_background_mod_pars(fr_pars)
    
        self.background_mod = background_mod
    
        return background_mod, background_pars
    
    
    def _make_spectrum_mod_pars(self, print_initial_pars=False):
        """Generate the peaks lmfit models and parameters."""
        params = Parameters()
    
        peaks_mod_pars = Peaks_Model(
            spectrum_vals = self.spectrum_vals,
            energy_vals = self.energy_vals,
            microscope_ID = self.microscope_ID,
            meas_mode = self.meas_mode,
            fitting_model = None,
            fitting_pars = params,
            xray_weight_refs_dict=self.el_lines_weight_refs_dict,
            is_particle=self.is_particle,
            free_area_el_lines=self.free_area_el_lines,
        )
        
        fitted_peaks = []
        for el_line in self.el_lines_list:
            is_peak_present = peaks_mod_pars._add_peak_model_and_pars(el_line)
            if is_peak_present:
                fitted_peaks.append(el_line)
    
        peaks_mod_pars._fix_overlapping_ref_peaks()
    
        fitted_elements = [el for el in self.els_to_fit_list if any(el + '_' in peak for peak in fitted_peaks)]
        fitted_els_to_quantify = [el for el in fitted_elements if el in self.els_to_quantify]
        
        if self.verbose:
            if len(fitted_elements) == 0:
                warnings.warn("No peak from the provided elements was found in the spectrum.", UserWarning)
            elif len(fitted_els_to_quantify) == 0:
                warnings.warn("No peak from the provided elements to quantify was found in the spectrum.", UserWarning)
    
        spectrum_mod, spectrum_pars = peaks_mod_pars.get_peaks_mod_pars()
    
        if self.fit_background:
            if len(fitted_els_to_quantify) > 0:
                fr_pars_els = fitted_els_to_quantify
            elif len(self.els_to_quantify) > 0:
                fr_pars_els = self.els_to_quantify
            elif len(self.els_substrate) > 0:
                fr_pars_els = self.els_substrate
            else:
                raise ValueError(
                    f"No valid element to fit was given. "
                    f"Please provide at least one element (in sample or substrate) that is not {calibs.undetectable_els}, "
                    f"or change the list 'undetectable_els' at calibs.__init__.py"
                )
    
            fr_pars = self._get_fraction_pars(fr_pars_els)
            spectrum_pars.update(fr_pars)
    
            background_mod, background_pars = self._get_background_mod_pars(fitted_els_to_quantify, fr_pars)
            if spectrum_mod is None:
                spectrum_mod = background_mod
            else:
                spectrum_mod += background_mod
            spectrum_pars.update(background_pars)
        else:
            self.background_mod = None
    
        self.spectrum_mod = spectrum_mod
        self.spectrum_pars = spectrum_pars
        self.fitted_els = fitted_elements
    
        if print_initial_pars:
            spectrum_pars.pretty_print()
        
    
    def _iteration_callback(self, params, iter, resid, *args, **kws):
        """Callback function to monitor fit iterations during optimization."""
        self.iteration_counter += 1
    
        if self.verbose and self.iteration_counter % 20 == 0:
            reduced_chi_square = np.sum(resid**2) / (len(self.energy_vals) - len(self.spectrum_pars))
            logger.debug(f"  Iter. #: {self.iteration_counter}. Residual sum of squares: {reduced_chi_square:.5e}")
    
        if self.print_evolving_params:
            print_single_separator()
            logger.debug(f"  Params changed in iteration #{self.iteration_counter}")
            if self.iteration_counter == 1:
                self.param_values = {param: params[param].value for param in params}
            else:
                for param in params:
                    par_value = params[param].value
                    if par_value != self.param_values[param] and params[param].vary:
                        logger.debug(f"{param}: {par_value}")
                        self.param_values[param] = par_value

                        if param == 'rhoz_par_slope':
                            bckngrd_contains_nan = np.any(np.isnan(self.background_mod.eval(params=params, x=self.energy_vals)))
                            logger.debug("Background contains nan vals: %s", bckngrd_contains_nan)
                            sp_contains_nan = np.any(np.isnan(self.spectrum_mod.eval(params=params, x=self.energy_vals)))
                            logger.debug("Spectrum contains nan vals: %s", sp_contains_nan)
                        
                        
        
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
        self.iteration_counter = 0
    
        if parameters is None:
            self._make_spectrum_mod_pars()
        else:
            self.spectrum_pars = parameters
            if self.fit_background:
                self._initialise_Background_Model()
    
        params = self.spectrum_pars
    
        if self.verbose:
            print_double_separator()
            print_double_separator()
            if n_iter:
                logger.debug(f"  ▶️ Iteration #{n_iter}")
            print_single_separator()
            logger.info('🔬 Fitting spectrum...')
            start_time = time.time()
    
        if initial_par_vals:
            for par, val in initial_par_vals.items():
                self.spectrum_pars[par].value = val
                
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
            logger.info(f'✅ Fit completed in {fitting_time:.1f} s with {self.iteration_counter} steps')
    
        used_params = list(fit_result.params.keys())
        fitted_lines = ['_'.join(param.split('_')[:-1]) for param in used_params if "area" in param]
    
        self.fit_result = fit_result
    
        if plot_result and self.verbose:
            self.plot_result()

        if print_result or self.verbose:
            self.print_result(print_only_independent_params=False, extended=print_result_extended)
    
        return fit_result, fitted_lines
    
     
    def plot_result(self):
        """Plot the fitted EDS spectrum and its individual background components."""
        fig = self.fit_result.plot(xlabel='Energy (keV)', ylabel='Counts')
    
        if self.fit_background:
            components = self.fit_result.eval_components(x=self.energy_vals)
            abs_att_param_name = [s for s in components.keys() if '_abs_att' in s][0]
            gen_bckgrnd_param_name = [s for s in components.keys() if '_generated_bckgrnd' in s][0]
            det_eff_par_name = '_det_efficiency'
            bcksctr_corr_par_name = '_backscattering_correction'
            stop_power_par_name = '_stopping_power'
            det_zero_peak_par_name = 'det_zero_peak_'
    
            gen_background = components[gen_bckgrnd_param_name]
            abs_attenuation = components[abs_att_param_name]
            det_efficiency = components[det_eff_par_name]
            bcksctr_corr = components[bcksctr_corr_par_name]
            stopping_power = components[stop_power_par_name]
            det_zero_peak = components[det_zero_peak_par_name]
            total_background_component = gen_background * abs_attenuation * det_efficiency * bcksctr_corr
    
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
        """Print a summary of the fit results."""
        print_double_separator()
        if extended:
            if print_only_independent_params:
                logger.info('Parameter           Value')
                for name, param in self.fit_result.params.items():
                    logger.info(f'{name:20s} {param.value:11.10f}')
            else:
                logger.info(self.fit_result.fit_report())
        else:
            reduced_chi_square = self.fit_result.redchi
            r_squared = 1 - self.fit_result.residual.var() / np.var(self.spectrum_vals)
            logger.info('📊 Fit results:')
            logger.info(f"  Reduced Chi-square: {reduced_chi_square:.2f}")
            logger.info(f"  R-squared: {r_squared:.5f}")
