#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated X-Ray Spectral Acquisition and Analysis of Powder Mixtures

This script configures and runs automated collection and (optionally) quantification
of EDS/WDS spectra for powder samples using an electron microscope (EM) with
a specified substrate and calibration.

Typical usage:
    - Edit the 'samples' list to define your standards or unknowns
    - Adjust configuration parameters as needed
    - Run the script to perform spectrum collection and (optionally)
        quantification for one or multiple samples at a time

Note:
    This script is exactly the same as Run_Acquisition_Quant_Analysis.py, except for is_known_powder_mixture_meas = True

Created on Fri Jul 26 09:34:34 2024

@author: Andrea
"""

from autoemxsp.runners import batch_acquire_and_analyze

# =============================================================================
# General Configuration
# =============================================================================

microscope_ID = 'PhenomXL'
microscope_type = 'SEM'
measurement_type = 'EDS'
measurement_mode = 'point'
quantification_method = 'PB'
beam_energy = 15  # keV
spectrum_lims = (14, 1100)  # eV

working_distance = 5 #mm

use_instrument_background = False
interrupt_fits_bad_spectra = True

max_analytical_error_percent = 5
min_bckgrnd_cnts = 5
quant_flags_accepted = [0, -1]

max_n_clusters = 6
show_unused_comps_clust = True

output_filename_suffix = ''

# =============================================================================
# Sample Definitions
# =============================================================================
sample_type = 'powder' # powder, bulk
sample_halfwidth = 3  # mm
sample_substrate_type = 'Ctape' # Ctape, None
sample_substrate_shape = 'circle' # square, circle
sample_substrate_width_mm = 12 # Al stub diameter, in mm

samples = [
    {'ID': 'In2O3_Al2O3', 'els': ['Ca', 'Al', 'Si', 'O'], 'pos': (-37.5, -37.5), 'cnd': ['In2O3', 'Al2O3']},
]

# =============================================================================
# Acquisition and Quantification Options
# =============================================================================

is_manual_navigation = False
is_auto_substrate_detection = True

auto_adjust_brightness_contrast = True
contrast = 4.3877  # Used if auto_adjust_brightness_contrast = False
brightness = 0.4504  # Used if auto_adjust_brightness_contrast = False

quantify_spectra = False

min_n_spectra = 100
max_n_spectra = 200

target_Xsp_counts = 50000
max_XSp_acquisition_time = target_Xsp_counts / 10000 * 5

# Substrate elements (may depend on target_Xsp_counts)
els_substrate = ['C', 'O', 'Al']  # N and F also detectable with >100k counts

# =============================================================================
# Powder sample options
# =============================================================================
powder_meas_cfg_kwargs = dict(
    is_manual_particle_selection = False,
    is_known_powder_mixture_meas = True,
    max_n_par_per_frame=30,
    max_spectra_per_par=5,
    max_area_par=10000.0,
    min_area_par=10.0,
    par_mask_margin=1.0,
    xsp_spots_distance_um=1.0,
    par_brightness_thresh=100,
    par_xy_spots_thresh=100,
    par_feature_selection = 'random',
    par_spot_spacing = 'random'
)

# =============================================================================
# Bulk sample options
# =============================================================================
bulk_meas_cfg_kwargs = dict(
    grid_spot_spacing_um = 100.0, # µm
    min_xsp_spots_distance_um = 5.0, # µm
    randomize_frames = False,
    exclude_sample_margin = False
)

# =============================================================================
# Run
# =============================================================================
comp_analyzer = batch_acquire_and_analyze(
    samples=samples,
    microscope_ID=microscope_ID,
    microscope_type=microscope_type,
    measurement_type=measurement_type,
    measurement_mode=measurement_mode,
    quantification_method = quantification_method,
    sample_type=sample_type,
    sample_halfwidth=sample_halfwidth,
    sample_substrate_type=sample_substrate_type,
    sample_substrate_shape=sample_substrate_shape,
    sample_substrate_width_mm=sample_substrate_width_mm,
    working_distance = working_distance,
    beam_energy=beam_energy,
    spectrum_lims=spectrum_lims,
    use_instrument_background=use_instrument_background,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
    max_analytical_error_percent=max_analytical_error_percent,
    min_bckgrnd_cnts=min_bckgrnd_cnts,
    quant_flags_accepted=quant_flags_accepted,
    max_n_clusters=max_n_clusters,
    show_unused_comps_clust=show_unused_comps_clust,
    is_manual_navigation=is_manual_navigation,
    is_auto_substrate_detection=is_auto_substrate_detection,
    auto_adjust_brightness_contrast=auto_adjust_brightness_contrast,
    contrast=contrast,
    brightness=brightness,
    quantify_spectra=quantify_spectra,
    min_n_spectra=min_n_spectra,
    max_n_spectra=max_n_spectra,
    target_Xsp_counts=target_Xsp_counts,
    max_XSp_acquisition_time=max_XSp_acquisition_time,
    els_substrate=els_substrate,
    powder_meas_cfg_kwargs=powder_meas_cfg_kwargs,
    bulk_meas_cfg_kwargs=bulk_meas_cfg_kwargs,
    output_filename_suffix=output_filename_suffix,
    development_mode=False,
    verbose=True,
)