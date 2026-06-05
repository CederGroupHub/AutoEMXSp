#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quantification correction models for EDS peak-to-background analysis."""

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import sympy as sp
from pymatgen.core import Element

import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.utils.helper import (
    EDSError,
    atomic_to_weight_fr,
    print_double_separator,
    print_nice_1d_row,
    print_single_separator,
    weight_to_atomic_fr,
)
from autoemx.data import J_df, xray_mass_absorption_coeff
from autoemx.core.fitter import Background_Model


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
