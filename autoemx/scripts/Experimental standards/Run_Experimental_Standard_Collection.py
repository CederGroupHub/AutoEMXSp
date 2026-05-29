#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated X-Ray Experimental Standard Acquisition and Analysis

This script configures and runs automated collection and fitting
of EDS/WDS spectra from experimental standards (i.e., samples of known composition)
to generate reference values of peak-to-background ratios.

Requirements:
    - Proper instrument calibration files and instrument driver for the selected microscope

Typical usage:
    - Edit the 'std_list' list to define your standards
        Use preferably bulk standards. Powder standards are also acceptable, especially when making standards for analysing known precursor mixtures.
    - Adjust configuration parameters as needed
    - Run the script to collect experimental standards for one or multiple samples at a time

Created on Fri Aug 20 09:34:34 2025

@author: Andrea
"""
# =============================================================================
# General Configuration
# =============================================================================
microscope_ID = 'PhenomXL'
microscope_type = 'SEM'
measurement_type = 'EDS'
spectrum_lims = (14, 1100)  # eV

use_instrument_background = False

min_bckgrnd_cnts = 10

output_filename_suffix = ''

exp_std_dir = None # Defines directory where measurements are saved. If None, uses default path.
# =============================================================================
# Sample Definitions
# =============================================================================

std_list = [
    {'ID': 'Al','formula': 'Al', 'pos': (26.263,-21.261), 'sample_type': 'bulk'},
    # {'ID': 'Al2O3_prec_BM','formula': 'Al2O3', 'pos': (-38.829, 40.011), 'sample_type': 'powder', 'is_manual_meas' : False},
    {'ID': 'Pt','formula': 'Pt0.99 O0.1', 'pos': (), 'sample_type': 'bulk', 'els_to_use_for_mean_PB_calc' : ['all'],'is_manual_meas' : True},
]
# els_to_use_for_mean_PB_calc: elements to include as standards; default is ["all"] if sample is "bulk", None otherwise. (all elements/lines). Use ["Fe", "O"] for specific elements, []/["none"] for none.

sample_substrate_type = 'Ctape'
# =============================================================================
# Acquisition Options and Sample description
working_distance = 8.5 #mm
is_auto_substrate_detection = False

is_manual_meas = False # If True, all samples are measured in manual navigation mode. Set to False to use automatic navigation for all samples. This can be overridden on a per-sample basis by setting 'is_manual_meas' in the sample dict.

fit_during_collection= True
update_std_library = True
# If True, writes fitted corrected PB reference values into EDS_Stds_{voltage}keV.json.
# If False, acquisition/fitting still runs and CSV outputs are saved, but the standards
# dictionary JSON is left unchanged.

sample_substrate_shape = 'circle'
sample_halfwidth = 3  # mm

measurement_mode = 'point'
beam_energy = 15  # keV

auto_adjust_brightness_contrast = True
contrast = None # 4.3877  # Used if auto_adjust_brightness_contrast = False
brightness = None # 0.4504  # Used if auto_adjust_brightness_contrast = False
saved_images_extension = 'png' # Default lightweight output. Set 'tif' for higher-resolution/larger files.
save_raw_images = False # Default saves only the annotated image. Set True to also save raw images.

n_target_spectra = 15
max_n_spectra = 200

target_Xsp_counts = 50000
max_XSp_acquisition_time = target_Xsp_counts / 10000 * 5

# Substrate elements (may depend on target_Xsp_counts)
els_substrate = ['C', 'O', 'Al']  # Contaminants that may be present in the spectrum

# =============================================================================
# Powder options
# =============================================================================
powder_meas_cfg_kwargs = dict(
    is_manual_particle_selection = False,
    is_known_powder_mixture_meas = False,
    img_shift_tracking = True,
    max_n_par_per_frame=30,
    max_spectra_per_par=3,
    max_area_par=10000.0,
    min_area_par=10.0,
    par_mask_margin=1.0,
    xsp_spots_distance_um=1.0,
    par_brightness_thresh=100,
    par_xy_spots_thresh=100,
    par_feature_selection = 'peaks',
    par_spot_spacing = 'random'
)

# =============================================================================
# Bulk options
# =============================================================================
bulk_meas_cfg_kwargs = dict(
    grid_spot_spacing_um = 100.0, # µm
    min_xsp_spots_distance_um = 5.0, # µm
    image_frame_width_um = None, # µm
    randomize_frames = False,
    exclude_sample_margin = False
)

# =============================================================================
# Options for experimental standard collection
# =============================================================================
exp_stds_meas_cfg_kwargs = dict(
    min_acceptable_PB_ratio = 10,
    quant_flags_accepted = [0],
    els_to_use_for_mean_PB_calc = ["all"] if not powder_meas_cfg_kwargs["is_known_powder_mixture_meas"] else ["none"],
    generate_separate_std_dict = powder_meas_cfg_kwargs["is_known_powder_mixture_meas"]
)

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.batch_acquire_experimental_stds import batch_acquire_experimental_stds

exp_std_maker = batch_acquire_experimental_stds(
    stds=std_list,
    microscope_ID=microscope_ID,
    microscope_type=microscope_type,
    measurement_type=measurement_type,
    measurement_mode=measurement_mode,
    sample_halfwidth=sample_halfwidth,
    sample_substrate_type=sample_substrate_type,
    sample_substrate_shape=sample_substrate_shape,
    working_distance = working_distance,
    beam_energy=beam_energy,
    spectrum_lims=spectrum_lims,
    use_instrument_background=use_instrument_background,
    min_bckgrnd_cnts=min_bckgrnd_cnts,
    is_manual_meas=is_manual_meas,
    fit_during_collection= fit_during_collection,
    update_std_library = update_std_library,
    is_auto_substrate_detection=is_auto_substrate_detection,
    auto_adjust_brightness_contrast=auto_adjust_brightness_contrast,
    contrast=contrast,
    brightness=brightness,
    saved_images_extension=saved_images_extension,
    save_raw_images=save_raw_images,
    min_n_spectra=n_target_spectra,
    max_n_spectra=max_n_spectra,
    target_Xsp_counts=target_Xsp_counts,
    max_XSp_acquisition_time=max_XSp_acquisition_time,
    els_substrate=els_substrate,
    powder_meas_cfg_kwargs=powder_meas_cfg_kwargs,
    bulk_meas_cfg_kwargs=bulk_meas_cfg_kwargs,
    exp_stds_meas_cfg_kwargs=exp_stds_meas_cfg_kwargs,
    output_filename_suffix=output_filename_suffix,
    development_mode=False,
    verbose=True,
    exp_std_dir = exp_std_dir
)