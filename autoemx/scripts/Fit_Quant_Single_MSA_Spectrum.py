#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting and quantification of a single X-ray spectrum, extracted from a .msa file (typycally exported by EDS software).

For spectrum-level analysis of fitting and quantification performance.

Run this module directly to process a spectrum with the defined configuration options

Created on Thu Jan 15 15:49:51 2026

@author: Andrea

Citation
--------
If you use this package, please cite:
A. Giunto et al., Accurate SEM-EDS Quantification, Automation, and Machine Learning Enable
High-Throughput Compositional Characterization of Powders,
Nature Communications 17, 9735 (2026).
https://doi.org/10.1038/s41467-026-76633-x
"""
import os

from autoemx.utils import load_msa
from autoemx.config import defaults
from autoemx.runners.fit_and_quantify_spectrum import fit_and_quantify_spectrum

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

#%% Sample definition
spectrum_path = os.path.join('input', 'Example_spectrum.msa')

els_sample = ['Bi','Fe','O']
els_substrate = ['C']

is_particle = True

#%% Fit and Quant definitions
quantify = True

spectrum_lims = (14, 1100)
fit_tol = 1e-3
max_undetectable_w_fr = 0

interrupt_fits_bad_spectra = False
force_single_iteration = False

print_results = True
plot_signal = True

#%% Instrument definition
microscope_ID = defaults.microscope_ID
measurement_type = defaults.measurement_type
measurement_mode = defaults.measurement_mode

#%% Extract metadata
_, spectrum_vals, metadata = load_msa(spectrum_path)
# print(metadata)
en_axis_units = metadata['XUNITS']
if 'keV' in en_axis_units:
    en_scaling_factor = 1
elif 'eV' in en_axis_units:
    en_scaling_factor = 1000
else:
    raise ValueError(f"Energy axis unit '{en_axis_units}' unrecognized. Please correct")

offset = float(metadata['OFFSET'])/en_scaling_factor
scale = float(metadata['XPERCHAN'])/en_scaling_factor

beam_en_key = [k for k in metadata.keys() if 'BEAMKV' in k][0]
beam_energy = float(metadata[beam_en_key])

emergence_angle_key = [k for k in metadata.keys() if 'ELEVANGLE' in k][0]
emergence_angle = float(metadata[emergence_angle_key])

livetime_key = [k for k in metadata.keys() if 'LIVETIME' in k][0]
sp_collection_time = float(metadata[livetime_key])

#%% Quantification
quantifier = fit_and_quantify_spectrum(
    spectrum_vals = spectrum_vals,
    microscope_ID = microscope_ID,
    meas_type = measurement_type,
    meas_mode = measurement_mode,
    det_ch_offset = offset,
    det_ch_width = scale,
    beam_energy = beam_energy,
    emergence_angle = emergence_angle,
    sp_collection_time = sp_collection_time,
    background_vals=None,
    els_sample = els_sample,
    els_substrate = els_substrate,
    spectrum_lims = spectrum_lims,
    quantify = quantify,
    plot_signal = plot_signal,
    fit_tol = fit_tol,
    is_particle = is_particle,
    max_undetectable_w_fr = max_undetectable_w_fr,
    force_single_iteration = force_single_iteration,
    interrupt_fits_bad_spectra = interrupt_fits_bad_spectra,
    print_results = print_results,
    quant_verbose = True,
    fitting_verbose = True
    )
