#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peaks model module for EDS spectrum fitting.

Handles construction and parameterization of spectral peaks including characteristic X-ray lines,
escape peaks, and pileup peaks.

Main Class:
    - Peaks_Model: Manages peak model construction and parameterization
"""

import re
import numpy as np
from itertools import combinations
from scipy.special import erfc
from scipy.signal import find_peaks, peak_prominences
from lmfit import Model, Parameters, Parameter
from pymatgen.core import Element

import autoemx.config.defaults as dflt
import autoemx.calibrations as calibs
from autoemx.data.Xray_lines import get_el_xray_lines
from .detector_response import DetectorResponseFunction


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
    - Requires DetectorResponseFunction to be initialized first for correct peak sigma calculations.
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
    escape_peaks_str = '_escSiKa'
    pileup_peaks_str = '_pileup'

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
        """Return a Gaussian function."""
        return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x-center)**2 / (2*sigma**2))
    
    
    @staticmethod
    def _gaussian_tail(x, area, center, sigma, gamma):
        """Return a skewed tail function for low-Z X-ray peaks."""
        return area / (2 * gamma * sigma * np.exp(-1 / (2 * gamma**2))) \
            * np.exp((x - center) / (gamma * sigma)) \
            * erfc((x - center) / (np.sqrt(2) * sigma) + 1 / (np.sqrt(2) * gamma))
    
    
    @staticmethod
    def _gaussian_with_tail(x, area, center, sigma, gamma, f_tail):
        """Return a Gaussian with a low-energy tail."""
        return Peaks_Model._gaussian(x, area, center, sigma) + f_tail * Peaks_Model._gaussian_tail(x, area, center, sigma, gamma)
    
    
    @staticmethod
    def _gaussian_shelf(x, area, center, sigma):
        """Return a shelf function for X-ray peaks."""
        return area / (2 * center) * erfc((x - center) / (np.sqrt(2) * sigma))
    
    
    @staticmethod
    def _gaussian_with_tail_and_shelf(x, area, center, sigma, gamma, f_tail, f_shelf):
        """Return a Gaussian with both a low-energy tail and a shelf component."""
        return (Peaks_Model._gaussian(x, area, center, sigma)
                + f_tail * Peaks_Model._gaussian_tail(x, area, center, sigma, gamma)
                + f_shelf * Peaks_Model._gaussian_shelf(x, area, center, sigma))
    
    
    @staticmethod
    def _gaussian_with_tail_and_icc(x, area, center, sigma, R_e, F_loss, gamma, f_tail):
        """Compute a Gaussian peak with a low-energy tail, convolved with the incomplete charge collection (ICC) model."""
        key_center = f"{center:.6f}"
        key_F_loss = f"{F_loss:.6f}"
        key_R_e = f"{R_e * 1e7:.5f}"
        key_f_tail = f"{f_tail:.6f}"
        key_gamma = f"{gamma:.6f}"
        key = f"{key_center}_{key_F_loss}_{key_R_e}_{key_f_tail}_{key_gamma}"
    
        icc_cache = Peaks_Model.icc_freq_spectra
        icc_n_vals_distr = icc_cache.get(key)
        if icc_n_vals_distr is None:
            icc_n_vals_distr = DetectorResponseFunction.get_icc_spectrum(x, center, R_e, F_loss)
            icc_cache[key] = icc_n_vals_distr
    
        g_vals = Peaks_Model._gaussian_with_tail(x, area, center, sigma, gamma, f_tail)
        g_with_icc = np.convolve(np.array(g_vals), np.array(icc_n_vals_distr), mode='same')
    
        return g_with_icc
    
    
    @staticmethod
    def _gaussian_with_icc(x, area, center, sigma, R_e, F_loss):
        """Compute a Gaussian peak convolved with the incomplete charge collection (ICC) model."""
        key_center = f"{center:.6f}"
        key_F_loss = f"{F_loss:.6f}"
        key_R_e = f"{R_e * 1e7:.5f}"
        key = f"{key_center}_{key_F_loss}_{key_R_e}"
    
        icc_cache = Peaks_Model.icc_freq_spectra
        icc_n_vals_distr = icc_cache.get(key)
        if icc_n_vals_distr is None:
            icc_n_vals_distr = DetectorResponseFunction.get_icc_spectrum(x, center, R_e, F_loss)
            icc_cache[key] = icc_n_vals_distr
    
        g_vals = Peaks_Model._gaussian(x, area, center, sigma)
        g_with_icc = np.convolve(np.array(g_vals), np.array(icc_n_vals_distr), mode='same')
    
        return g_with_icc
    
    # =============================================================================
    # Functions for spectral data analysis
    # =============================================================================
    def _identify_peak(self, line_energy, sigma_peak, el_line):
        """Identify and characterize the presence of a peak in the spectrum."""
        channel_width_keV = self.energy_vals[1] - self.energy_vals[0]

        peak_indices = np.where(
            (self.energy_vals > line_energy - sigma_peak * 3) &
            (self.energy_vals < line_energy + sigma_peak * 3)
        )[0]
    
        if len(peak_indices) == 0:
            return peak_indices, [], None, False
    
        peak_sp_vals = self.spectrum_vals[peak_indices]
    
        filter_len = 2
        peak_sp_vals_blurred = np.convolve(peak_sp_vals, np.ones(filter_len) / filter_len, mode='valid')
    
        min_prominence = np.mean(np.array([*peak_sp_vals_blurred[:2], *peak_sp_vals_blurred[-2:]])) / 10
        min_prominence = max(min_prominence, 4)
    
        peak_len = int(sigma_peak / channel_width_keV * 6)
    
        peaks, _ = find_peaks(
            peak_sp_vals_blurred,
            prominence=min_prominence,
            wlen=peak_len,
            distance=5
        )
        peak_pos = peak_indices[peaks]
        peak_height = list(self.spectrum_vals[peak_pos])
    
        peak_prom, _, _ = peak_prominences(peak_sp_vals_blurred, peaks, wlen=peak_len)
        if peak_prom.size != 0:
            is_small_peak = peak_prom[0] < 50
        else:
            is_small_peak = None
        
        if peak_height:
            is_peak_overlapping = False
        else:
            peak_energy_vals = self.energy_vals[peak_indices]
            slope = np.polyfit(peak_energy_vals, peak_sp_vals, deg=1)[0]
            overlap_threshold = 100 * (channel_width_keV / 0.01)
            if abs(slope) < overlap_threshold:
                is_peak_overlapping = False
            else:
                is_peak_overlapping = True
    
        return peak_indices, peak_height, is_small_peak, is_peak_overlapping
    
    
    def _estimate_peak_area(self, line_energy, sigma_peak, el_line):
        """Estimate the area of a characteristic X-ray peak in the spectrum."""
        peak_indices, peak_height, is_small_peak, is_peak_overlapping = self._identify_peak(
            line_energy, sigma_peak, el_line
        )
    
        if len(peak_indices) == 0:
            peak_area = 0
            return peak_area, is_small_peak, is_peak_overlapping
    
        if peak_height:
            peak_area = peak_height[0] * (sigma_peak * np.sqrt(2 * np.pi)) * 0.7
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
        """Create a lmfit Parameter for the center (energy) of a Gaussian X-ray peak."""
        if line_energy < 0.6:
            center_par = Parameter(
                line_prefix + self.center_key,
                value=line_energy,
                vary=True,
                min=line_energy * 0.93,
                max=line_energy * 1.03
            )
        else:
            center_par = Parameter(
                line_prefix + self.center_key,
                value=line_energy,
                vary=True,
                min=line_energy * 0.99,
                max=line_energy * 1.01
            )
        return center_par


    def _get_gaussian_sigma_param(self, line_prefix, sigma_init):
        """Create an lmfit Parameter for the sigma (width) of a Gaussian X-ray peak."""
        return Parameter(
            line_prefix + self.sigma_key,
            value=sigma_init,
            vary=True,
            min=sigma_init * 0.93,
            max=sigma_init * 1.1
        )


    def _add_peak_model_and_pars(self, el_line):
        """Add a peak model and its parameters for a specific X-ray line to the composite model."""
        
        model = self.fitting_model
        params = self.fitting_params
        
        el, line = el_line.split('_', 1)
        line_prefix = el_line + '_'
        
        is_escape_peak = self.escape_peaks_str in line
        is_pileup_peak = self.pileup_peaks_str in line
    
        escape_ref_line = None
        is_escape_ref_peak = False
        if is_escape_peak:
            escape_ref_el_line = el_line[:-len(self.escape_peaks_str)]
            escape_ref_line = escape_ref_el_line.split('_')[1]
            is_escape_ref_peak = escape_ref_el_line in self.xray_weight_refs_lines
    
        pileup_ref_line = None
        is_pileup_ref_peak = False
        if is_pileup_peak:
            pileup_ref_el_line = el_line[:-len(self.pileup_peaks_str)]
            pileup_ref_line = pileup_ref_el_line.split('_')[1]
            is_pileup_ref_peak = pileup_ref_el_line in self.xray_weight_refs_lines
            
        lines = get_el_xray_lines(el)
        if is_escape_peak:
            if escape_ref_line not in lines:
                raise ValueError(
                    f"escape_ref_line '{escape_ref_line}' not found in X-ray lines for element '{el}'"
                )
            else:
                line_energy = lines[escape_ref_line]['energy (keV)'] - 1.740
        elif is_pileup_peak:
            if pileup_ref_line not in lines:
                raise ValueError(
                    f"pileup_ref_line '{pileup_ref_line}' not found in X-ray lines for element '{el}'"
                )
            else:
                line_energy = lines[pileup_ref_line]['energy (keV)'] * 2
        else:
            line_energy = lines[line]['energy (keV)']
        
        center_par = self._get_gaussian_center_param(line_prefix, line_energy)
        sigma_init = DetectorResponseFunction._det_sigma(line_energy)
        sigma_par = self._get_gaussian_sigma_param(line_prefix, sigma_init)

        initial_area, is_small_peak, is_peak_overlapping = self._estimate_peak_area(line_energy, sigma_init, el_line)
        
        if el_line not in self.xray_weight_refs_lines:
            ref_line = self.xray_weight_refs_dict[el_line]
            ref_line_prefix = ref_line + '_'
            ref_peak_area_param = params.get(ref_line_prefix + self.area_key)
            
            if ref_peak_area_param is None or ref_peak_area_param.value == 0:
                peak_m = None
            else:
                if is_escape_peak:
                    peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
                    try:
                        escape_peak_probability = calibs.escape_peak_probability[self.meas_mode] # type: ignore
                    except Exception:
                        escape_peak_probability = dflt.escape_peak_probability
                    if is_escape_ref_peak:
                        max_area_escape = escape_peak_probability * ref_peak_area_param.value
                        params.add(line_prefix + self.area_key, value=initial_area, min=0, max=max_area_escape)
                    else:
                        weight = lines[escape_ref_line]['weight']
                        escape_ref_line_prefix = ref_line + self.escape_peaks_str + '_'
                        if params.get(escape_ref_line_prefix + self.area_key) is not None:
                            params.add(line_prefix + self.area_key, expr=escape_ref_line_prefix + self.area_key + f'*{weight}')
                        else:
                            max_area_escape = weight * escape_peak_probability * ref_peak_area_param.value
                            params.add(line_prefix + self.area_key, value=initial_area, min=0, max=max_area_escape)
                elif is_pileup_peak:
                    peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
                    try:
                        pileup_peak_probability = calibs.pileup_peak_probability[self.meas_mode] # type: ignore
                    except Exception:
                        pileup_peak_probability = dflt.pileup_peak_probability
                    if is_pileup_ref_peak:
                        max_area_pileup = pileup_peak_probability * ref_peak_area_param.value
                        vary_area = max_area_pileup > 0.1
                        params.add(line_prefix + self.area_key, value=initial_area, vary=vary_area, min=0, max=max_area_pileup)
                    else:
                        weight = lines[pileup_ref_line]['weight']
                        pileup_ref_line_prefix = ref_line + self.pileup_peaks_str + '_'
                        if pileup_ref_line_prefix + self.area_key in params.keys():
                            params.add(line_prefix + self.area_key, expr=pileup_ref_line_prefix + self.area_key + f'*{weight}')
                        else:
                            peak_m = None
                else:
                    calibrate_peak_shape_params = el in self.free_peak_shapes_els
                    peak_m, _ = self._get_peak_model_and_update_pars(
                        params, line_energy, line_prefix,
                        is_calibration=calibrate_peak_shape_params,
                        ref_line_prefix=ref_line_prefix
                    )
                    if el_line in self.free_area_el_lines:
                        params.add(line_prefix + self.area_key, value=initial_area, min=0, vary=True)
                    else:
                        weight = lines[line]['weight']
                        params.add(line_prefix + self.area_key, expr=ref_line_prefix + self.area_key + f'*{weight}')
    
                if peak_m is not None:
                    params.add(line_prefix + self.center_key, expr=f"{line_energy} - {ref_line_prefix}{self.center_offset_key}")
                    params.add(line_prefix + self.sigma_key, expr=f"{sigma_init} * {ref_line_prefix}{self.sigma_broadening_key}")
        
        else:
            if initial_area == 0:
                params.add(line_prefix + self.area_key, value=initial_area, min=0, vary=False)
                peak_m = None
            else:
                vary_area = True
                if line == 'Ll':
                    ref_line_prefix = el + '_Ka1_'
                    ref_peak_area_param = params.get(ref_line_prefix + self.area_key)
                    try:
                        weight_Ll_ref_Ka1 = calibs.weight_Ll_ref_Ka1[self.meas_mode] # type: ignore
                    except Exception:
                        weight_Ll_ref_Ka1 = dflt.weight_Ll_ref_Ka1
                    if ref_peak_area_param is None or ref_peak_area_param.value == 0:
                        initial_area = 0
                        max_area = 1
                        vary_area = False
                    elif not self.is_particle:
                        max_area = weight_Ll_ref_Ka1 * ref_peak_area_param.value
                        is_small_peak = True
                    else:
                        max_area = 5 * weight_Ll_ref_Ka1 * ref_peak_area_param.value
                else:
                    max_area = np.inf
            
                if is_small_peak:
                    params.add(line_prefix + self.sigma_key, value=sigma_init, vary=False)
                else:
                    params.add(sigma_par)
                
                params.add(center_par)
                params.add(line_prefix + self.area_key, value=initial_area, min=0, max=max_area, vary=vary_area)
    
                params.add(line_prefix + self.center_offset_key, expr=f"{line_energy} - {line_prefix}{self.center_key}")
                params.add(line_prefix + self.sigma_broadening_key, expr=f"{line_prefix}{self.sigma_key} / {sigma_init}")
    
                calibrate_peak_shape_params = el in self.free_peak_shapes_els
                peak_m, _ = self._get_peak_model_and_update_pars(
                    params, line_energy, line_prefix,
                    is_calibration=calibrate_peak_shape_params
                )

        if peak_m is not None:
            if model is None:
                model = peak_m
            else:
                model += peak_m
    
        if peak_m is None or (is_small_peak is None and not is_peak_overlapping):
            is_peak_present = False
        else:
            is_peak_present = True
    
        self.fitting_model = model
        self.fitting_params = params
    
        return is_peak_present
    
    
    def _get_peak_model_and_update_pars(self, params, line_energy, line_prefix, is_calibration=False, ref_line_prefix=None):
        """Select and configure the appropriate peak model and update shape parameters."""
    
        gamma, f_tail, R_e, F_loss = calibs.get_calibrated_peak_shape_params(line_energy, self.meas_mode)
    
        if not is_calibration:
            if line_energy < 1.18:
                peak_m = Model(Peaks_Model._gaussian_with_tail, prefix=line_prefix)
                params.add(line_prefix + self.gamma_key, value=gamma, vary=False)
                params.add(line_prefix + self.tail_fraction_key, value=f_tail, vary=False)
    
            elif line_energy <= 1.839:
                peak_m = Model(Peaks_Model._gaussian_with_tail_and_icc, prefix=line_prefix)
                params.add(line_prefix + self.gamma_key, value=gamma, vary=False)
                params.add(line_prefix + self.tail_fraction_key, value=f_tail, vary=False)
                params.add(line_prefix + self.F_loss_key, value=F_loss, vary=False)
                params.add(line_prefix + self.R_e_key, value=R_e, vary=False)
    
            elif line_energy <= 5:
                peak_m = Model(Peaks_Model._gaussian_with_icc, prefix=line_prefix)
                params.add(line_prefix + self.F_loss_key, value=F_loss, vary=False)
                params.add(line_prefix + self.R_e_key, value=R_e, vary=False)
    
            else:
                peak_m = Model(Peaks_Model._gaussian, prefix=line_prefix)
    
        else:
            if not ref_line_prefix:
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
        """Identify and constrain overlapping reference peaks."""
    
        params = self.fitting_params
    
        center_params = [pname for pname in params if self.center_key in pname]
    
        free_peaks = {}
        for pname in center_params:
            if params[pname].vary:
                peak_prefix = pname[:-(len(self.center_key) + 1)]
                free_peaks[peak_prefix] = params[pname].value
    
        peaks_to_fix = set()
        for peak1, peak2 in combinations(free_peaks, 2):
            center1 = free_peaks[peak1]
            center2 = free_peaks[peak2]
            sigma1 = DetectorResponseFunction._det_sigma(center1)
            if abs(center1 - center2) < sigma1 * 3:
                peaks_to_fix.add((peak1, peak2))
    
        if not hasattr(self, 'fixed_peaks_dict'):
            self.fixed_peaks_dict = {}
    
        for peak1, peak2 in peaks_to_fix:
            fixed_peaks = list(self.fixed_peaks_dict.keys())
            if peak1 not in fixed_peaks and peak2 not in fixed_peaks:
                params = self._fix_center_sigma_peak(params, peak1, peak2)
            elif peak1 in fixed_peaks and peak2 in fixed_peaks:
                continue
            else:
                if peak1 in fixed_peaks:
                    fixed_peak = peak1
                    dep_peak = peak2
                else:
                    fixed_peak = peak2
                    dep_peak = peak1
                if self.fixed_peaks_dict[fixed_peak] == '':
                    params = self._fix_center_sigma_peak(params, fixed_peak, dep_peak)
                else:
                    params = self._fix_center_sigma_peak(params, self.fixed_peaks_dict[fixed_peak], dep_peak)
    
        self.fitting_params = params
    
    
    def _fix_center_sigma_peak(self, params, ref_peak, dep_peak):
        """Tie the center and sigma of a dependent peak to those of a reference (independent) peak."""
    
        dep_peak_center = params[f"{dep_peak}_{self.center_key}"].value
    
        ref_peak_offset = params[f"{ref_peak}_{self.center_offset_key}"].name
        center_expr = f"{dep_peak_center} - {ref_peak_offset}"
        params[f"{dep_peak}_{self.center_key}"].expr = center_expr
    
        params[f"{ref_peak}_{self.sigma_key}"].vary = False
        params[f"{dep_peak}_{self.sigma_key}"].vary = False
    
        params[f"{ref_peak}_{self.area_key}"].value /= 2
        params[f"{dep_peak}_{self.area_key}"].value /= 2
    
        self.fixed_peaks_dict[ref_peak] = ''
        self.fixed_peaks_dict[dep_peak] = ref_peak
    
        return params
