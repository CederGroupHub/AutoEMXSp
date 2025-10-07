#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microscope XS Calibration Parameters

This module defines all calibration constants and functions required for
 X-ray spectroscopy (XS) quantification for a specific instrument setup.
 
It centralizes all physical and empirical calibration data needed for
accurate spectrum modeling, including detector characteristics, 
energy calibration, peak shape, background continuum, and system-specific 
factors. The parameters and functions here are referenced throughout the 
quantification pipeline to ensure consistency and traceability.
 
Currently, this module defines calibration parameters specifics to
scanning electron microscope (SEM) energy-dispersive X-ray spectroscopy (EDS).

References:
    - Goldstein et al., "Scanning Electron Microscopy and X-ray Microanalysis", 4th Ed.
    - N.W.M. Ritchie, Spectrum Simulation in DTSA-II, Microsc. Microanal. 15 (2009) 454–468. https://doi.org/10.1017/S1431927609990407

Author: Andrea Giunto
Created on: Mon Jan 20 15:40:42 2025

Contents
--------
- Detector physical and operational parameters (e.g., density, channel count, take-off angle)
- Lists of elements undetectable by the system and available EDS acquisition modes
- System-specific calibration constants for:
    * Escape peak and pileup probabilities
    * Peak shape and tailing (incomplete charge collection, ICC)
    * Gaussian noise and energy calibration
    * Channel-to-energy calibration per EDS mode
    * Background continuum modeling (Duncumb model and modifications)
- Functions to retrieve calibrated parameters for peak shape and background,
  as a function of X-ray energy and EDS mode

Typical Usage
-------------
This module is imported by higher-level spectrum modeling and fitting code.
For example:

    >>> from calibrations import get_calibrated_peak_shape_params, detector_channel_params
    >>> gamma, f_tail, R_e, F_loss = get_calibrated_peak_shape_params(6.4, 'point')
    >>> offset = detector_channel_params['point']['offset']

Notes
-----
- **IMPORTANT:** This file must be **copied and customized** for every new detector or microscope system.
  Update all calibration values, constants, and empirical fits to match the specific hardware and calibration data
  for each instrument. Failure to do so will result in inaccurate quantification.
- All calibration parameters should be reviewed and updated regularly as the
  microscope system is serviced or recalibrated.
- If additional EDS modes or new detector configurations are added, update the
  relevant dictionaries and functions accordingly.
- All formulas and empirical fits are traceable to published literature or 
  in-house calibration data (see References).
"""

import numpy as np
import sympy as sp

# =============================================================================
# EDS Detector Physical Characteristics
# =============================================================================

Si_density = 2.33  # g/cm^3
"""Silicon density at detector operating temperature (235 K)."""

detector_ch_n = 2048
"""Number of energy channels in the EDS detector."""

emergence_angle = 28.5  # degrees
"""X-ray emergence (take-off) angle, in degrees."""

undetectable_els = ['H', 'He', 'Li']
"""
List of elements that are undetectable or poorly detected by this EDS system.
Note: Li and Be detectability depends on detector window and configuration.
"""

# =============================================================================
# Detector-Specific Calibration Parameters
# =============================================================================
R_e_background = 50e-7
F_loss_background = 0.27
"""
Incomplete charge collection (ICC) parameters for background fitting.
Model parameters from:
    Redus, R. H., & Huber, A. C. (2015). Response Function of Silicon Drift Detectors for Low Energy X-rays.
    In Advances in X-ray Analysis (AXA) (pp. 274–282). International Centre for Diffraction Data (ICDD).

Calibration Guidance
--------------------
- These parameters can be manually adjusted by observing the background intensity "jump" at ~1.84 keV,
  which corresponds to the Si absorption edge in the detector response.
- For calibration, use high-count (high-statistics) spectra to minimize the effect of noise.
- The exact values are not critical for spectrum fitting unless they are set to extreme or unphysical values.
  Reasonable estimates suffice for most applications.
"""

# =============================================================================
# Detector Noise and Peak Width Parameters
# =============================================================================
"""
Parameters controlling the detector noise and the width (sigma) of spectral peaks.

These values are empirically calibrated by measuring the widths (sigmas) of
characteristic X-ray peaks at various energies using high-statistics spectra
from bulk standards. The observed sigmas are then fitted to the following model:

    sigma = sqrt(elec_noise**2 + conv_eff * E * F)

where:
    - sigma: Gaussian standard deviation of the peak (in keV)
    - elec_noise: Electronic noise (in eV)
    - conv_eff: Conversion efficiency (in keV)
    - F: Fano factor (dimensionless)
    - E: X-ray energy (in keV)

Calibration Procedure
---------------------
1. Acquire high-count spectra for a series of well-characterized elements with
   characteristic X-ray lines spanning the relevant energy range.
2. For each element, measure the sigma (standard deviation) of the main peak.
    This can be obtained from peak sigma parameters output by lmfit model.
3. Fit the measured sigmas as a function of energy E using the formula above.
   This can be done using the script or notebook: `calculation_pars_peak_fwhm.py`.
4. Update the values below with the fitted parameters.

Notes
-----
- These parameters are system- and detector-specific, and should be recalibrated
  after hardware changes or detector servicing.
- The values are not highly sensitive for routine fitting, but accurate calibration
  improves peak shape modeling and quantitative results.
- See: Ritchie, N.W.M., "Spectrum Simulation in DTSA-II", Microsc. Microanal. 15 (2009) 454–468.
"""

F = 0.12 # From Ritchie, 2009
"""(Fixed) Fano factor (dimensionless), characterizing statistical fluctuations in charge generation."""

elec_noise = 4.627
"""(Fitted) Electronic noise (eV), characterizing the detector's intrinsic noise."""

conv_eff = 0.003724  # keV
"""(Fitted) Conversion efficiency (keV), relating to charge collection efficiency."""


# =============================================================================
# Calibration Parameters Dependent on Detector and EDS Mode
# =============================================================================

available_meas_modes = ('point', 'map')
"""
Tuple of valid EDS acquisition modes for this system.
All mode-dependent calibration parameters must have entries for each mode.
"""


# =============================================================================
# Fine-Tuning and Empirical System Calibrations
# =============================================================================
"""
Empirical parameters for fine-tuning the EDS spectral model.

Calibration Procedure (Common to All Parameters Below)
-----------------------------------------------------
- These values are empirical and not highly sensitive for most quantification tasks; reasonable estimates usually suffice.
- Calibration can be performed manually by observing the fitted values in a few high-quality spectra of EDS standards.
- To do this, allow the parameter(s) to vary freely in your lmfit model (set `vary=True`), fit the spectrum, and observe the optimized value reported by lmfit.
- If desired, repeat for several standard spectra to confirm consistency.
- Adjust the values below to match typical fitted results, or use the average.
- For most applications, precise calibration is not critical unless the parameter is set to an extreme or unphysical value.
"""

weight_Ll_ref_Ka1 = {
    'point': 0.05
}
"""
Estimated upper bound for Ll/Ka1 intensity ratio (for 12 < Z < 20, no La line).
Prevents overestimation of Ll peak, e.g. for Al and Si.
Actual value may fit lower during analysis; precise calibration is not critical.
"""

escape_peak_probability = {
    'point': 0.03  # 3%
}
"""
Estimated upper bound for escape peak probability (fraction).
Actual value may fit lower during analysis; precise calibration is not critical.
"""

pileup_peak_probability = {
    'point': 0.003  # 0.3%
}
"""
Estimated upper bound for pileup peak probability (fraction).
Actual value may fit lower during analysis; precise calibration is not critical.
"""

gen_background_time_scaling_factor = {
    'point': 0.0853
}
"""
Scaling factor (K) in modified Duncumb continuum model, yielding F=1 for bulk standards of high atomic number.
Used to guess initial value of background scaling factor K.
K is allowed to vary within a small interval during fitting.
"""

strobe_peak_int_factor = {
    'point': 1
}
"""Multiplicative factor to convert total spectrum collection time to initial guessed value of strobe peak intensity.
Strobe peak intensity is allowed to vary within a small interval during fitting.
"""

zero_strobe_peak_sigma = {
    'point': 0.16
}
"""Standard deviation (sigma) of Gaussian zero strobe peak.
Fixed during fitting.
"""

# =============================================================================
# Calibration Functions for Peak Shape and Background
# =============================================================================

def get_calibrated_peak_shape_params(E, meas_mode):
    """
    Returns peak shape parameters for incomplete charge collection (ICC) low-energy tailing.
    
    Calibration Procedure
    ---------------------
    The energy dependence and values of these parameters should be calibrated empirically
    using high-count spectra from bulk standards, ideally pure elements. The recommended
    workflow is:
    
    1. For each element of interest, fit spectra individually by adding the element symbol
       to `free_peak_shapes_els` in the `Peaks_Model` class (inside `EDS_spectrum_fitter.py`).
       This allows the ICC peak shape parameters to be freely optimized for that element.
    2. Perform the fitting using your main fitting code or the script `FitSpectra.py`, which
       can automate fitting across multiple spectra and elements and output the optimized
       parameters for each case.
    3. Extract the fitted ICC parameters (gamma, f_tail, R_e, F_loss) from the lmfit output
       for each element/energy.
    4. Fit the energy dependence of each parameter (as a function of X-ray energy E) using
       suitable empirical functions, such as the sigmoid forms implemented below.
    5. Update the fitted coefficients in this function with the new calibration results.
    
    Parameters
    ----------
    E : float
        X-ray energy in keV.
    meas_mode : str
        EDS mode, must be in available_meas_modes.
    
    Returns
    -------
    tuple
        (gamma, f_tail, R_e, F_loss)
        - gamma: Tailing shape parameter (dimensionless)
        - f_tail: Fraction of tailing (dimensionless)
        - R_e: ICC parameter for tailing (dimensionless)
        - F_loss: ICC parameter for tailing (dimensionless)
    
    Raises
    ------
    ValueError
        If meas_mode is not recognized.
    
    Notes
    -----
    - These parameters are system- and detector-specific, and should be recalibrated
      after hardware changes or detector servicing.
    - Accurate calibration of ICC parameters improves modeling of low-energy tailing
      and the quantitative accuracy of EDS fits, especially for light elements.
    - See also: Redus, R. H., & Huber, A. C. (2015). Response Function of Silicon Drift Detectors
      for Low Energy X-rays. In Advances in X-ray Analysis (AXA), ICDD.
    """
    def sigmoid(x, L, x0, k):
        return L / (1 + np.exp(-k * (x - x0)))

    if meas_mode == 'point':
        # Fitted gamma
        gamma = 1 + sigmoid(E, 3.5, 0.88, 9.63)
        # Fitted f_tail
        f_tail = 0.0094 + sigmoid(E, 0.1022, 0.926, -11.42)
        # Fitted R_e
        R_e = (80.67 + sigmoid(E, 101.43, 1.9051, 8.6)) * 1e-7
        # Fitted F_loss
        F_loss = 0.1 + sigmoid(E, 0.1912, 1.7586, -40.3)
        return gamma, f_tail, R_e, F_loss

    raise ValueError(f"Missing peak shape parameters for meas_mode '{meas_mode}'. "
                     "Add to function get_calibrated_peak_shape_params.")


def get_calibrated_background_params(meas_mode):
    """
    Returns parameters for the modified Duncumb background model.
    
    Calibration Procedure
    ---------------------
    The background model parameters (P, F, and beta) should be empirically calibrated
    using high-count spectra from bulk standards. The recommended procedure is:
    
    1. In the `Background_Model` class (inside `EDS_spectrum_fitter.py`),
       set the `is_calibration` parameter to `True` in the function 
       `get_generated_background_mod_pars`. This enables background parameter calibration mode.
    2. For each parameter (P, F, beta), select whether it should be varied freely or held fixed
       by manually setting `vary=True/False` in the relevant parameter definition.
       This can be done to calibrate one parameter at a time or all together.
    3. Fit high-statistics spectra from well-characterized standards, ideally covering a range of
       atomic numbers Z, to optimize the background parameters.
       Perform the fitting using the script `FitSpectra.py`, which can automate fitting across
       multiple spectra and elements and output the optimized parameters for each case.
    4. Extract the fitted values for P, F, and beta (as a function of Z) from the lmfit output.
    5. Fit the Z-dependence of each parameter using suitable empirical functions, such as the
       expressions implemented below (polynomial, exponential, sigmoid, Gaussian, etc.).
    6. Update the coefficients in this function with the new calibration results.
    
    Parameters
    ----------
    meas_mode : str
        EDS mode, must be in available_meas_modes.
    
    Returns
    -------
    tuple
        (gen_bckgrnd_Duncumb_P, gen_bckgrnd_Duncumb_F, beta_expr)
        - gen_bckgrnd_Duncumb_P: str, expression for background P parameter as a function of Z
        - gen_bckgrnd_Duncumb_F: str, expression for background F parameter as a function of Z
        - beta_expr: sympy expression for beta parameter as a function of Z
    
    Raises
    ------
    ValueError
        If meas_mode is not recognized.
    
    Notes
    -----
    - These parameters are system- and detector-specific, and should be recalibrated after
      hardware changes or detector servicing.
    - Accurate calibration of the background model is essential for quantification, especially for
      light elements and low-count regions, since it is used for composition-dependent
      background corrections in the PB method.
    - The functional forms and coefficients should be updated as needed to best fit the
      empirical data for your instrument.
    """
    
    if meas_mode == 'point':
        Z = sp.symbols('Z')
        # Duncumb's generated background P parameter
        gen_bckgrnd_Duncumb_P = str(0.00002041 * Z**2 - 0.004076 * Z + 1.252)
        # Duncumb's generated background F parameter
        gen_bckgrnd_Duncumb_F = str((1 - sp.exp(-(Z + 0.7453) / 6.1543)))
        # Beta parameter (sigmoid + Gaussian)
        sigmoid_part = 0.0419 + 0.2105 / (1 + sp.exp(-0.4788 * (Z - 45.9098)))
        gaussian_part = 0.2807 * sp.exp(-((Z - 28.0063) ** 2) / (2 * 4.6587 ** 2))
        beta_expr = sigmoid_part + gaussian_part
        return gen_bckgrnd_Duncumb_P, gen_bckgrnd_Duncumb_F, beta_expr

    raise ValueError(f"Missing background parameters for meas_mode '{meas_mode}'. "
                     "Add to function get_calibrated_background_params.")