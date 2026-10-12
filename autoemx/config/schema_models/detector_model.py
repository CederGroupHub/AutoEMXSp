#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Detector model of a microscope: the physical and empirical parameters used to model its X-ray spectra.

The values of each microscope are saved in ``autoemx/calibrations/<microscope_ID>/detector_model.json``, next to its
detector efficiency (detector_efficiency.msa). The functional forms (sigmoids for the peak shape, modified Duncumb
model for the background) are defined here; the file holds their coefficients.

To add a microscope, add a folder with these two files, named with its microscope ID, and its driver in
autoemx/microscope_drivers. The driver is the file named as the microscope ID, unless ``driver`` names another one:
instruments of the same model, each with its own calibrated detector model (e.g. 'PhenomXL_LabB'), share the driver
of that model (driver 'PhenomXL').

This model replaces the module XS_calibrations.py of each microscope (AutoEMX <= 1.4.1).

It centralizes all physical and empirical calibration data needed for accurate spectrum modeling, including detector
characteristics, peak shape, background continuum, and system-specific factors. The channel-to-energy calibration
(offset, bin width) is not part of it: it is saved separately in Detector_channel_params_calibs, since it is
recalibrated regularly (Run_SDD_Calibration.py).

Notes
-----
- The values must be customized for every new detector or microscope system: update all calibration values and
  empirical fits to match the specific hardware and calibration data of each instrument. Failure to do so will result
  in inaccurate quantification.
- All calibration parameters should be reviewed and updated regularly, and after the microscope system is serviced.
- If additional EDS modes are added, every mode-dependent parameter needs an entry for the new mode.
- All formulas and empirical fits are traceable to published literature or in-house calibration data:
    - Goldstein et al., "Scanning Electron Microscopy and X-ray Microanalysis", 4th Ed.
    - N.W.M. Ritchie, Spectrum Simulation in DTSA-II, Microsc. Microanal. 15 (2009) 454–468.
      https://doi.org/10.1017/S1431927609990407
    - R. H. Redus, A. C. Huber (2015). Response Function of Silicon Drift Detectors for Low Energy X-rays. In Advances
      in X-ray Analysis (AXA), pp. 274–282. International Centre for Diffraction Data (ICDD).

Fine-tuning parameters (weight_Ll_ref_Ka1, escape_peak_probability, pileup_peak_probability,
gen_background_time_scaling_factor, strobe_peak_int_factor, zero_strobe_peak_sigma)
-------------------------------------------------------------------------------------
- These values are empirical and not highly sensitive for most quantification tasks; reasonable estimates suffice.
- Calibration can be performed manually by observing the fitted values in a few high-quality spectra of EDS standards.
- To do this, allow the parameter(s) to vary freely in your lmfit model (set `vary=True`), fit the spectrum, and
  observe the optimized value reported by lmfit. If desired, repeat for several standard spectra to confirm
  consistency, and use the typical fitted results, or their average.
- Precise calibration is not critical unless the parameter is set to an extreme or unphysical value.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from pydantic import BaseModel, ConfigDict, Field, model_validator

DETECTOR_MODEL_SCHEMA_VERSION = 1


class Sigmoid(BaseModel):
    """scale * (offset + amplitude / (1 + exp(-steepness * (E - center))))"""

    model_config = ConfigDict(extra="forbid")

    offset: float
    amplitude: float
    center: float
    steepness: float
    scale: float = 1.0

    def __call__(self, E):
        return (self.offset + self.amplitude / (1 + np.exp(-self.steepness * (E - self.center)))) * self.scale


class PeakShapeParams(BaseModel):
    """
    Peak shape parameters for incomplete charge collection (ICC) low-energy tailing, as functions of the X-ray energy
    E (keV): gamma (tailing shape), f_tail (fraction of tailing), R_e and F_loss (ICC parameters), all dimensionless.

    Calibration Procedure
    ---------------------
    The energy dependence and values of these parameters should be calibrated empirically using high-count spectra
    from bulk standards, ideally pure elements, with the script Calibrate_Peak_Shape_Params.py: it refits the
    spectra with the chosen parameters (gamma, f_tail, R_e, F_loss) fitted freely for the chosen elements, fits the
    chosen sigmoid coefficients to their values as a function of energy, and prints the coefficients to update.

    Accurate calibration of ICC parameters improves modeling of low-energy tailing and the quantitative accuracy of
    EDS fits, especially for light elements. These parameters are system- and detector-specific, and should be
    recalibrated after hardware changes or detector servicing (Redus & Huber, 2015).
    """

    model_config = ConfigDict(extra="forbid")

    gamma: Sigmoid
    f_tail: Sigmoid
    R_e: Sigmoid
    F_loss: Sigmoid


class DuncumbF(BaseModel):
    """1 - exp(-(Z + z_offset) / z_scale)"""

    model_config = ConfigDict(extra="forbid")

    z_offset: float
    z_scale: float


class BetaParams(BaseModel):
    """
    sigmoid_offset + sigmoid_amplitude / (1 + exp(-sigmoid_steepness * (Z - sigmoid_center)))
    + gaussian_amplitude * exp(-(Z - gaussian_center)**2 / (2 * gaussian_sigma**2))
    """

    model_config = ConfigDict(extra="forbid")

    sigmoid_offset: float
    sigmoid_amplitude: float
    sigmoid_steepness: float
    sigmoid_center: float
    gaussian_amplitude: float
    gaussian_center: float
    gaussian_sigma: float


class BackgroundParams(BaseModel):
    """
    Parameters of the modified Duncumb background model, as functions of the atomic number Z: P (quadratic polynomial,
    coefficients from the highest power), F, and beta.

    Calibration Procedure
    ---------------------
    The background model parameters (P, F, and beta) should be empirically calibrated using high-count spectra from
    bulk standards, ideally pure elements covering a range of atomic numbers Z, with the script
    Calibrate_Background_Params.py: it refits the spectra with the chosen parameters (P, F, beta) fitted freely, fits
    the chosen coefficients of their functions to their values as a function of Z, and prints the coefficients to
    update. One parameter can be calibrated at a time, or all together.

    Accurate calibration of the background model is essential for quantification, especially for light elements and
    low-count regions, since it is used for composition-dependent background corrections in the PB method. These
    parameters are system- and detector-specific, and should be recalibrated after hardware changes or detector
    servicing.
    """

    model_config = ConfigDict(extra="forbid")

    duncumb_P: List[float] = Field(min_length=3, max_length=3)
    duncumb_F: DuncumbF
    beta: BetaParams


class DetectorModel(BaseModel):
    """Detector model of a microscope (see the module docstring)."""

    model_config = ConfigDict(extra="forbid")

    schema_version: int = DETECTOR_MODEL_SCHEMA_VERSION
    driver: Optional[str] = Field(
        default=None,
        description="Microscope driver (autoemx/microscope_drivers/<driver>.py). Default: the microscope ID.",
    )

    # --- EDS detector physical characteristics
    Si_density: float = Field(description="Silicon density at detector operating temperature (235 K), in g/cm^3.")
    detector_ch_n: int = Field(description="Number of energy channels in the EDS detector.")
    emergence_angle: float = Field(description="X-ray emergence (take-off) angle, in degrees.")
    undetectable_els: List[str] = Field(
        description="Elements that are undetectable or poorly detected by this EDS system. Li and Be detectability "
                    "depends on detector window and configuration."
    )
    available_meas_modes: Tuple[str, ...] = Field(
        description="Valid EDS acquisition modes for this system. All mode-dependent parameters must have entries for "
                    "each mode used."
    )

    # --- ICC parameters for background fitting (Redus & Huber, 2015)
    R_e_background: float = Field(
        description="Incomplete charge collection (ICC) parameter for background fitting. Can be manually adjusted by "
                    "observing the background intensity 'jump' at ~1.84 keV, the Si absorption edge in the detector "
                    "response, in high-count spectra. The exact value is not critical unless extreme or unphysical."
    )
    F_loss_background: float = Field(description="ICC parameter for background fitting (see R_e_background).")

    # --- Detector noise and peak width: sigma = conv_eff * sqrt(elec_noise**2 + E * F / conv_eff), E in keV
    F: float = Field(
        description="(Fixed) Fano factor (dimensionless), characterizing statistical fluctuations in charge "
                    "generation (Ritchie, 2009)."
    )
    elec_noise: float = Field(
        description="(Fitted) Electronic noise (electrons, rms), characterizing the detector's intrinsic noise. "
                    "elec_noise and conv_eff are calibrated by fitting "
                    "sigma = conv_eff * sqrt(elec_noise**2 + E * F / conv_eff) to the sigmas of characteristic peaks "
                    "at various energies E (keV), fitted in high-count spectra of bulk standards, with "
                    "Calibrate_Peak_Sigma_Params.py. Recalibrate after hardware changes or detector servicing."
    )
    conv_eff: float = Field(
        description="(Fitted) Conversion efficiency: energy per electron-hole pair (keV), see elec_noise."
    )

    # --- Fine-tuning parameters, by EDS mode (see the module docstring)
    weight_Ll_ref_Ka1: Dict[str, float] = Field(
        description="Estimated upper bound for the Ll/Ka1 intensity ratio (for 12 < Z < 20, no La line). Prevents "
                    "overestimation of the Ll peak, e.g. for Al and Si. The fitted value may be lower."
    )
    escape_peak_probability: Dict[str, float] = Field(
        description="Estimated upper bound for the escape peak probability (fraction). The fitted value may be lower."
    )
    pileup_peak_probability: Dict[str, float] = Field(
        description="Estimated upper bound for the pileup peak probability (fraction). The fitted value may be lower."
    )
    gen_background_time_scaling_factor: Dict[str, float] = Field(
        description="Scaling factor K in the modified Duncumb continuum model, yielding F=1 for bulk standards of high "
                    "atomic number. Used to guess the initial value of the background scaling factor K, which is "
                    "allowed to vary within a small interval during fitting."
    )
    strobe_peak_int_factor: Dict[str, float] = Field(
        description="Multiplicative factor converting the total spectrum collection time to the initial guess of the "
                    "strobe peak intensity, which is allowed to vary within a small interval during fitting."
    )
    zero_strobe_peak_sigma: Dict[str, float] = Field(
        description="Standard deviation (sigma) of the Gaussian zero strobe peak. Fixed during fitting."
    )

    # --- Calibrated functions, by EDS mode
    peak_shape: Dict[str, PeakShapeParams]
    background: Dict[str, BackgroundParams]

    @model_validator(mode="after")
    def _check_schema_version(self) -> "DetectorModel":
        if self.schema_version > DETECTOR_MODEL_SCHEMA_VERSION:
            raise ValueError(
                f"Detector model schema version {self.schema_version} is newer than the one supported by this "
                f"AutoEMX version ({DETECTOR_MODEL_SCHEMA_VERSION}): update AutoEMX."
            )
        return self

    @classmethod
    def from_json_file(cls, file_path: str | Path) -> "DetectorModel":
        with open(file_path, "r", encoding="utf-8") as handle:
            return cls.model_validate(json.load(handle))

    def to_json_file(self, file_path: str | Path, indent: int = 2) -> None:
        with open(file_path, "w", encoding="utf-8") as handle:
            json.dump(self.model_dump(mode="json"), handle, indent=indent)
            handle.write("\n")

    def fingerprint(self) -> str:
        """Short hash of all values of the detector model (not the driver), e.g. for caches of quantities computed from it."""
        payload = json.dumps(self.model_dump(mode="json", exclude={"driver"}), sort_keys=True)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]

    def get_calibrated_peak_shape_params(self, E, meas_mode: str):
        """
        Peak shape parameters for incomplete charge collection (ICC) low-energy tailing.

        Parameters
        ----------
        E : float or ndarray
            X-ray energy in keV.
        meas_mode : str
            EDS mode, must be in available_meas_modes.

        Returns
        -------
        tuple
            (gamma, f_tail, R_e, F_loss), see PeakShapeParams.
        """
        params = self.peak_shape.get(meas_mode)
        if params is None:
            raise ValueError(f"Missing peak shape parameters for meas_mode '{meas_mode}'. "
                             "Add them to 'peak_shape' in the detector model of the microscope.")
        return params.gamma(E), params.f_tail(E), params.R_e(E), params.F_loss(E)

    def get_calibrated_background_params(self, meas_mode: str):
        """
        Parameters of the modified Duncumb background model.

        Parameters
        ----------
        meas_mode : str
            EDS mode, must be in available_meas_modes.

        Returns
        -------
        tuple
            (gen_bckgrnd_Duncumb_P, gen_bckgrnd_Duncumb_F, beta_expr)
            - gen_bckgrnd_Duncumb_P: str, expression for the background P parameter as a function of Z
            - gen_bckgrnd_Duncumb_F: str, expression for the background F parameter as a function of Z
            - beta_expr: sympy expression for the beta parameter as a function of Z
        """
        params = self.background.get(meas_mode)
        if params is None:
            raise ValueError(f"Missing background parameters for meas_mode '{meas_mode}'. "
                             "Add them to 'background' in the detector model of the microscope.")
        Z = sp.symbols('Z')
        a2, a1, a0 = params.duncumb_P
        gen_bckgrnd_Duncumb_P = str(a2 * Z**2 + a1 * Z + a0)
        gen_bckgrnd_Duncumb_F = str((1 - sp.exp(-(Z + params.duncumb_F.z_offset) / params.duncumb_F.z_scale)))
        beta = params.beta
        sigmoid_part = beta.sigmoid_offset + beta.sigmoid_amplitude / (
            1 + sp.exp(-beta.sigmoid_steepness * (Z - beta.sigmoid_center)))
        gaussian_part = beta.gaussian_amplitude * sp.exp(
            -((Z - beta.gaussian_center) ** 2) / (2 * beta.gaussian_sigma ** 2))
        beta_expr = sigmoid_part + gaussian_part
        return gen_bckgrnd_Duncumb_P, gen_bckgrnd_Duncumb_F, beta_expr
