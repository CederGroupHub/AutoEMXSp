#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate the peak width (sigma) parameters of the detector model.

Fits the conversion efficiency (conv_eff) and electronic noise (elec_noise) of the detector to the Gaussian sigmas
of X-ray peaks, measured over a range of energies:

    sigma = conv_eff * sqrt(elec_noise**2 + E * F / conv_eff)

with the Fano factor F fixed (Ritchie, 2009). See autoemx.runners.calibrate_peak_sigma_params for details.

Typical usage:
    1. Acquire high-count spectra of bulk standards, ideally pure elements, with lines spanning the energy range
       (e.g. with Run_Experimental_Standard_Collection.py), and fit them.
    2. List the sample folders below, with the lines to use (None: Ka1, La1 and Ma1 lines).
    3. Run the script, and update 'conv_eff' and 'elec_noise' in autoemx/calibrations/<microscope_ID>/detector_model.json.
    4. If some sigmas are at the bounds of the fit (reported), refit the spectra and calibrate again.

Created on Mon Oct 12 2026

@author: Andrea
"""

# =============================================================================
# Input
# =============================================================================

# Fitted samples: {path of the sample folder: el_lines to use, or None for its Ka1, La1 and Ma1 lines}
samples = {
    '/absolute/path/to/Stds/Al': None,
    '/absolute/path/to/Stds/Cu': ['Cu_Ka1', 'Cu_La1'],
}

# Sigmas measured otherwise (optional): {label: (energy in keV, sigma in keV)}
# Example: {'Mn_Ka1': (5.8987, 0.0535)}
measured_sigmas = None

# Microscope whose detector model is compared. None: that of the samples.
microscope_ID = None

# =============================================================================
# Run
# =============================================================================

from autoemx.runners.calibrate_peak_sigma_params import calibrate_peak_sigma_params

calibrate_peak_sigma_params(samples=samples, measured_sigmas=measured_sigmas, microscope_ID=microscope_ID)
