#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Background model module for EDS spectrum fitting.

Handles calculation and fitting of X-ray background continuum, including absorption
attenuation, backscattering correction, and generated background models.

Main Class:
    - Background_Model: Comprehensive background modeling with multiple physics corrections
"""

import re
import numpy as np
from scipy.integrate import trapezoid
from lmfit import Model, Parameters
from lmfit.models import GaussianModel
from pymatgen.core import Element
import sympy as sp

import autoemx.calibrations as calibs
from autoemx.data import (
    xray_mass_absorption_coeff,
    J_df,
)
from autoemx.utils.helper import weight_to_atomic_fr
from .detector_response import DetectorResponseFunction


class Background_Model:
    """
    Model for calculating and fitting X-ray background in EDS spectra.

    Handles background continuum generation, absorption attenuation corrections,
    detector efficiency, backscattering effects, and detector response.

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
        """Initialize a Background_Model instance."""
        if beam_energy is not None:
            type(self).cls_beam_e = beam_energy

        self.is_particle = is_particle
        self.sp_collection_time = sp_collection_time
        self.tot_sp_counts = tot_sp_counts
        self.emergence_angle = emergence_angle
        self.energy_vals = energy_vals
        self.meas_mode = meas_mode
        self.els_w_fr = els_w_fr if els_w_fr is not None else {}

        self._clear_cached_abs_att_variables()

    @staticmethod
    def _clear_cached_abs_att_variables():
        """Reset class-level caches to ensure they are recalculated in new fits."""
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
        """Parse element weight fraction parameters and return elements, fractions, and atomic fractions."""
        par_pattern = r'f_[A-Z][a-z]*'
        els = []
        w_frs = []
        for el_fr_param, w_fr in el_fr_params.items():
            if re.match(par_pattern, el_fr_param):
                el = el_fr_param.split("_")[1]
                els.append(el)
                val = w_fr.value if hasattr(w_fr, 'value') else w_fr
                w_frs.append(val)

        if len(els) > 0:
            at_frs = weight_to_atomic_fr(w_frs, els, verbose=False)
        else:
            at_frs = []
    
        return els, w_frs, at_frs
        
    # =============================================================================
    # Atomic number averaging
    # =============================================================================
    @staticmethod
    def get_average_Z(els_symbols, method="Statham"):
        """Returns a symbolic formula string for the average atomic number (Z) of a sample."""
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
            raise ValueError(f"Unknown method '{method}' for averaging of compound atomic number.")
        
        
    # =============================================================================
    # Stopping power
    # =============================================================================
    @staticmethod
    def _stopping_power(x, adr_sp=1, **el_fr_params):
        """Computes the stopping power correction for a multi-element sample."""
        M_vals = []
        lnJ_vals = []
        
        els, w_frs, _ = Background_Model._get_els_frs(**el_fr_params)
        
        for el_symbol, w_fr in zip(els, w_frs):
            el = Element(el_symbol)
            Z_el = el.Z
            W_el = el.atomic_mass
            J_el = J_df.loc[Z_el, J_df.columns[0]] / 1000
            M = w_fr * Z_el / W_el
            M_vals.append(M)
            lnJ_vals.append(M * np.log(J_el))
        
        sum_M = sum(M_vals)
        ln_J = sum(lnJ_vals) / sum_M
        J_val = np.exp(ln_J)
        U0 = Background_Model.cls_beam_e / x
    
        S_vals = 1 / ((1 + 16.05 * (J_val / x) ** 0.5 * ((U0 ** 0.5 - 1) / (U0 - 1)) ** 1.07) / sum_M)
        
        if adr_sp == 1:
            S_vals = DetectorResponseFunction._apply_det_response_fncts(S_vals)
        
        S_vals = np.ones(len(x))
        
        return S_vals
    
    
    def get_stopping_power_mod_pars(self):
        """Returns an lmfit Model and Parameters for electron stopping power correction."""
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
        """Computes the total mass absorption coefficient (μ/ρ) for a compound."""
        mass_abs_coeff = np.zeros(len(x))
        els, w_frs, _ = Background_Model._get_els_frs(**el_fr_params)
        for el, w_fr in zip(els, w_frs):
            mass_abs_coeff += xray_mass_absorption_coeff(el, x) * w_fr
        return mass_abs_coeff


    @staticmethod
    def _abs_attenuation_Philibert(
        x, det_angle=28.5, abs_par=1.2e-6, abs_path_len_scale=1, adr_abs=1, **el_fr_params
    ):
        """Computes the absorption attenuation correction using a modified Philibert equation."""
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
        """Second-order absorption correction for P/B ratio."""
        chi = mass_abs_coeffs_sample / np.sin(np.deg2rad(det_angle))
        gamma = beam_energy**1.65 - x**1.65
    
        A_P = 1 + 3.34e-6 * (gamma * chi) + 5.59e-13 * (gamma * chi)**2
        A_B = 1 + 3.0e-6 * (gamma * chi) + 4.5e-13 * (gamma * chi)**2
        A_vals = A_B / A_P
    
        return A_vals
    
    
    @staticmethod
    def _abs_attenuation_phirho(
        x, det_angle=28.5, rhoz_par_offset=0, rhoz_par_slope=0, rhoz_lim=0.001, adr_abs=1, **el_fr_params
    ):
        """Absorption correction based on ionization depth distribution (phi-rho model)."""
        els, w_frs, at_frs = Background_Model._get_els_frs(**el_fr_params)
        
        E0 = Background_Model.cls_beam_e
        U0 = E0 / x
        mu = Background_Model._mass_abs_coeff(x, **el_fr_params)
        nu = Background_Model._nu_sample(E0, els, w_frs)
        
        phi0 = 1 + (nu * U0 * np.log(U0)) / (U0 - 1)
        gamma0 = (1 + nu) * (U0 * np.log(U0)) / (U0 - 1)
        
        alpha = 0
        beta = 0
        for el, at_fr in zip(els, at_frs):
            a_el, b_el = Background_Model._phi_alpha_beta_coeffs(x, E0, el, gamma0)
            alpha += a_el * at_fr
            beta += b_el * at_fr 
    
        recalc_den = False
        recalc_num = False
    
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
        abs_model[abs_model > 1] = 1
    
        if adr_abs == 1:
            abs_model = DetectorResponseFunction._apply_det_response_fncts(abs_model)
    
        return abs_model
    
    
    @staticmethod
    def _get_abs_att_num(alpha, beta, phi0, gamma0, mu, det_angle, rhoz_par_offset, rhoz_par_slope):
        """Computes the numerator integral for absorption-attenuation correction."""
        rhoz_grid = np.array(Background_Model.rhoz_values)
        num_vals = (
            Background_Model._phi_rhoz(rhoz_grid, alpha, beta, phi0, gamma0)
            * np.exp(-mu * (rhoz_grid + rhoz_par_offset + rhoz_grid * rhoz_par_slope) / np.sin(np.deg2rad(det_angle)))
        )
        num_int = trapezoid(num_vals, rhoz_grid)
        return num_int
        
    
    @staticmethod
    def _get_abs_att_den(alpha, beta, phi0, gamma0, det_angle):
        """Computes the denominator integral for absorption-attenuation correction."""
        rhoz_grid = np.array(Background_Model.rhoz_values)
        den_vals = Background_Model._phi_rhoz(rhoz_grid, alpha, beta, phi0, gamma0)
        den_int = trapezoid(den_vals, rhoz_grid)
        return den_int
    
    
    @staticmethod
    def _plot_phirho(
        en, alpha, beta, phi0, gamma0, mu,
        rhoz_par_offset, rhoz_par_slope, det_angle, rhoz_lim
    ):
        """Visualizes the φ(ρz) curves and their absorption-corrected integrands (development use only)."""
        import matplotlib.pyplot as plt
        
        num_integrand = lambda rhoz: (
            Background_Model._phi_rhoz(rhoz, alpha, beta, phi0, gamma0)
            * np.exp(-mu * (rhoz + rhoz_par_offset + rhoz * rhoz_par_slope) 
            / np.sin(np.deg2rad(det_angle)))
        )
        den_integrand = lambda rhoz: Background_Model._phi_rhoz(rhoz, alpha, beta, phi0, gamma0)

        rhoz_grid = np.array(Background_Model.rhoz_values)
        i_lim = np.argmin(np.abs(rhoz_grid - rhoz_lim))
        rhoz_cut = rhoz_grid[:i_lim]

        num_values = [num_integrand(rhoz) for rhoz in rhoz_grid]
        den_values = [den_integrand(rhoz) for rhoz in rhoz_grid]
        num_values_cut = [num_integrand(rhoz) for rhoz in rhoz_cut]
        den_values_cut = [den_integrand(rhoz) for rhoz in rhoz_cut]

        num_area = trapezoid(num_values, rhoz_grid) * 1000
        den_area = trapezoid(den_values, rhoz_grid) * 1000
        num_area_cut = trapezoid(num_values_cut, rhoz_cut)
        den_area_cut = trapezoid(den_values_cut, rhoz_cut)

        ratio_abs_change = (num_area_cut / den_area_cut) / (num_area / den_area)

        plt.figure(figsize=(10, 6))
        plt.plot(rhoz_grid, num_values, label='Numerator (with absorption)', color='b')
        plt.plot(rhoz_grid, den_values, label='Denominator (no absorption)', color='r')
        text_x = rhoz_grid[int(len(rhoz_grid) * 0.7)]
        text_y = max(den_values) * 0.8
        plt.text(
            text_x, text_y,
            f"Num Integral: {num_area:.3f}\nDen Integral: {den_area:.3f}",
            fontsize=10, bbox=dict(facecolor='white', alpha=0.7)
        )
        plt.title(f'Integrands vs ρz at {en:.3f} keV')
        plt.xlabel('ρz')
        plt.ylabel('Integrand Value')
        plt.legend()
        plt.grid(True)
        plt.show()

        return ratio_abs_change
    
    @staticmethod
    def _phi_rhoz(rhoz, alpha, beta, phi0, gamma0):
        """Computes the ionization depth distribution function ϕ(ρz)."""
        phi_rhoz = gamma0 * np.exp(-alpha**2 * rhoz**2) * (1 - ((gamma0 - phi0) / gamma0) * np.exp(-beta * rhoz))
        return phi_rhoz
    
    
    @staticmethod
    def _phi_alpha_beta_coeffs(x, E0, el, gamma0):
        """Computes the alpha and beta coefficients for the φ(ρz) model."""
        Z = Element(el).Z
        A = float(Element(el).atomic_mass)
        J = J_df.loc[Z, J_df.columns[0]] / 1000
    
        alpha = 2.14e5 * Z**1.16 / (A * E0**1.25) * (np.log(1.166 * E0 / J) / (E0 - x))**0.5
        beta = 1.1e5 * Z**1.5 / ((E0 - x) * A)
        
        return alpha, beta
    
    
    def get_abs_attenuation_mod_pars(self, model='phirho'):
        """Returns an lmfit Model and Parameters for background attenuation due to absorption."""
        
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
                f"Unknown model '{model}'. Choose 'phirho' or 'Philibert'."
            )
    
        return absorption_attenuation_m, params_abs_att


    # =============================================================================
    # Electron backscattering correction
    # =============================================================================
    @staticmethod
    def _backscattering_correction(x, adr_bcksctr=1, **el_fr_params):
        """Computes the backscattering correction factor."""
        E0 = Background_Model.cls_beam_e
        els, w_frs, _ = Background_Model._get_els_frs(**el_fr_params)
        nu_val = Background_Model._nu_sample(E0, els, w_frs)
        R_c_val = Background_Model._R_c(x, nu_val, E0)
        model = 1 - (1 - R_c_val) * (2 / (1 + nu_val))**0.63 * (0.79 + 0.44 * x / E0)
    
        if adr_bcksctr == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    
    
    @staticmethod
    def _nu_sample(E0, els, w_frs):
        """Computes the sample backscattering coefficient."""
        nu_val = 0.0
        for el, w_fr in zip(els, w_frs):
            Z = Element(el).Z
            nu_el = Background_Model._nu_el(Z, E0)
            nu_val += nu_el * w_fr
        return nu_val
    
    
    @staticmethod
    def _nu_el(Z, E0):
        """Computes the elemental backscattering coefficient."""
        nu20 = (-52.3791 + 150.48371 * Z - 1.67373 * Z**2 + 0.00716 * Z**3) * 1e-4
        G_nu20 = (-1112.8 + 30.289 * Z - 0.15498 * Z**2) * 1e-4
        nu_val = nu20 * (1 + G_nu20 * np.log(E0 / 20))
        return nu_val   

    
    @staticmethod
    def _R_c(x, nu_val, E0):
        """Computes the backscattering correction factor."""
        U0 = E0 / x
        I_val, G_val = Background_Model._return_IG(U0)
        R_c = 1 - nu_val * (I_val + nu_val * G_val) ** 1.67
        return R_c
    
    
    @staticmethod
    def _return_IG(U0):
        """Computes the I(U0) and G(U0) functions of overvoltage ratio."""
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
        """Returns an lmfit Model and Parameters for backscattering correction."""
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
        """Analytical model for generated bremsstrahlung (Castellano et al., 2004)."""
        E0 = Background_Model.cls_beam_e
    
        model = (
            K * np.sqrt(Z) * ((E0 - x) / x)
            * (a1 + a2 * x + a3 * np.log(Z) + a4 * (E0 ** a5) / Z)
            * (1 + (a6 + a7 * E0) * (Z / x))
        )
    
        if apply_det_response == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        model = np.where(model < 0, 0.01, model)
    
        return model
    
    
    @staticmethod
    def _generated_bckgrnd_Trincavelli1998(
        x, Z, K=0.035, a1=-54.86, a2=-1.072, a3=0.2835, a4=30.4, a5=875, a6=0.08, apply_det_response=1
    ):
        """Analytical model for generated bremsstrahlung (Trincavelli et al., 1998)."""
        E0 = Background_Model.cls_beam_e
    
        model = (
            K * np.sqrt(Z) * (E0 - x) / x
            * (a1 + a2 * x + a3 * E0 + a4 * np.log(Z) + a5 / (E0 ** a6 * Z ** 2))
        )
    
        if apply_det_response == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        model = np.where(model < 0, 0.01, model)
    
        return model
    
    
    @staticmethod
    def _generated_bckgrnd_DuncumbMod(
        x, Z, K=0.8, F=1, P=1, beta=0, apply_det_response=1
    ):
        """Analytical model for generated bremsstrahlung (Duncumb et al., 2001)."""
        E0 = Background_Model.cls_beam_e
    
        model = K * Z * F * ((E0 - x) / x) ** P * (x / (x + beta))
    
        if apply_det_response == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    
    
    @staticmethod
    def get_beta_expr(beta_expr, els_symbols):
        """Returns a symbolic formula string for beta parameter from Duncumb model."""
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
        """Returns an lmfit Model and Parameters for continuum X-ray background."""
    
        els, _, _ = Background_Model._get_els_frs(**fr_pars)
    
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
            bckgrnd_m = Model(Background_Model._generated_bckgrnd_DuncumbMod)
            Z_par_expr = Background_Model.get_average_Z(els, method='Statham')
    
            if self.sp_collection_time is not None and self.sp_collection_time > 0:
                K_val = self.sp_collection_time * calibs.gen_background_time_scaling_factor[self.meas_mode]
            else:
                K_val = self.tot_sp_counts / 1e5
    
            P_expr, F_expr, beta_expr_Z = calibs.get_calibrated_background_params(self.meas_mode)
    
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
            raise ValueError(f"Unknown background model '{model}'.")
    
    
    # =============================================================================
    # Detector efficiency and zero strobe peak
    # =============================================================================
    @staticmethod
    def _det_efficiency(x, adr_det_eff=1):
        """Returns the detector efficiency as a function of energy."""
        model = np.interp(
            x,
            DetectorResponseFunction.det_eff_energy_vals,
            DetectorResponseFunction.det_eff_vals
        )
    
        if adr_det_eff == 1:
            model = DetectorResponseFunction._apply_det_response_fncts(model)
    
        return model
    

    def get_detector_efficiency_mod_pars(self):
        """Returns an lmfit Model and Parameters for EDS detector efficiency correction."""
        detector_efficiency_m = Model(Background_Model._det_efficiency)
        detector_efficiency_pars = detector_efficiency_m.make_params(
            adr_det_eff=dict(expr='apply_det_response')
        )
    
        return detector_efficiency_m, detector_efficiency_pars
    
    
    def get_det_zero_peak_model_pars(self, amplitude_val):
        """Returns an lmfit GaussianModel and Parameters for detector zero (strobe) peak."""
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
        """Constructs the full background model and its parameters for spectral fitting."""
    
        gen_bckgrnd_mod, gen_bckgrnd_pars = self.get_generated_background_mod_pars(fr_pars, model='DuncumbMod')
        abs_att_mod, abs_att_pars = self.get_abs_attenuation_mod_pars(model='phirho')
        bs_cor_mod, bs_cor_pars = self.get_backscattering_correction_mod_pars()
        stopping_p_mod, stopping_p_pars = self.get_stopping_power_mod_pars()
        det_eff_mod, det_eff_pars = self.get_detector_efficiency_mod_pars()
    
        if self.sp_collection_time is not None and self.sp_collection_time > 0:
            amplitude_val = self.sp_collection_time * calibs.strobe_peak_int_factor[self.meas_mode]
        else:
            amplitude_val = self.tot_sp_counts / (10**4)
        det_zero_peak_mod, det_zero_peak_par = self.get_det_zero_peak_model_pars(amplitude_val)
    
        background_mod = (
            gen_bckgrnd_mod
            * abs_att_mod
            * det_eff_mod
            * bs_cor_mod
            * stopping_p_mod
            + det_zero_peak_mod
        )
    
        background_pars = Parameters()
        background_pars.update(gen_bckgrnd_pars)
        background_pars.update(abs_att_pars)
        background_pars.update(det_eff_pars)
        background_pars.update(bs_cor_pars)
        background_pars.update(stopping_p_pars)
        background_pars.update(det_zero_peak_par)
    
        return background_mod, background_pars
