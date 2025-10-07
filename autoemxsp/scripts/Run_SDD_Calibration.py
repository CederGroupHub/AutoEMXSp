#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated SDD Detector Calibration

This script automates the collection and fitting of X-ray spectra from experimental standards
to generate an SDD calibration file for AutoEMXsp.

Requirements:
- Proper instrument calibration files and instrument driver for the selected microscope

Typical Usage:
    - Edit the 'std_list' list to define your standards.
        - ID: A preferred ID is fine.
        - formula: Composition of standard, for accurate spectral fitting.
        - ref_el: Element in the composition whose peak should be taken as reference to calibrate the SDD.
        - ref_peak: Characteristic X-ray to take as reference to calibrate the SDD (e.g., Ka1, La1, Ma1, Mz1).
        - pos: Position of standard in the microscope stage.
        - sample_type: 'bulk' or 'powder'. Check autoemxsp.tools.config_classes.SampleConfig for updated support.
        - is_manual_meas: Set to True to manually select spots to measure.

    - Suggestions:
        - Preferably use bulk standards. Powder standards are also acceptable if bulk standards are not available.
        - Choose peaks above 2 keV, well-distanced between themselves.

    - Run the script to collect and fit spectra, with automated generation of SDD calibration file.

Created on Fri Aug 20 09:34:34 2025

@author: Andrea
"""
import os
from datetime import datetime
import pandas as pd

from autoemxsp.runners import batch_acquire_experimental_stds
from autoemxsp.runners import batch_fit_spectra
from autoemxsp.lib.Xray_lines import get_el_xray_lines
import autoemxsp.XSp_calibs as calibs
import autoemxsp.tools.constants as cnst


# =============================================================================
# General Configuration
# =============================================================================
microscope_ID = 'PhenomXL'
microscope_type = 'SEM'
measurement_type = 'EDS'
spectrum_lims = (14, 1100)  # eV

use_instrument_background = False

min_bckgrnd_cnts = 3

now = datetime.now()
now_formatted = now.strftime("%Y%m%d_%Hh%Mm")
output_filename_suffix = f'_{now_formatted}_50kcnts'

# =============================================================================
# Sample Definitions - Add two pure elements to measure
# =============================================================================
Cu_center = (37.863,38.195)
std_list = [
    {'ID': 'Cu','formula': 'Cu', 'ref_el' : 'Cu', 'ref_peak': 'Ka1', 'pos': Cu_center, 'sample_type': 'bulk', 'is_manual_meas' : False},
   {'ID': 'Al','formula': 'Al', 'ref_el' : 'Al', 'ref_peak': 'Ka1', 'pos': tuple(a + b for a, b in zip(Cu_center, (5, 0))), 'sample_type': 'bulk', 'is_manual_meas' : False}, # CAlibration standard center
]

# =============================================================================
# Acquisition Options and Sample description
working_distance = 5.5 #mm
is_auto_substrate_detection = False

fit_during_collection= False
update_std_library = False

sample_substrate_type = 'None'
sample_substrate_shape = 'circle'
sample_halfwidth = 1  # mm

measurement_mode = 'point'
beam_energy = 15  # keV

auto_adjust_brightness_contrast = True
contrast = None # 4.3877  # Used if auto_adjust_brightness_contrast = False
brightness = None # 0.4504  # Used if auto_adjust_brightness_contrast = False

n_target_spectra = 5
max_n_spectra = 10

target_Xsp_counts = 50000
max_XSp_acquisition_time = target_Xsp_counts / 10000 * 5

# Substrate elements (may depend on target_Xsp_counts)
els_substrate = ['C', 'O', 'Al']  # Contaminants that may be present in the spectrum

# =============================================================================
# Powder options
# =============================================================================
powder_meas_cfg_kwargs = dict(
    is_manual_particle_selection = False,
    is_known_powder_mixture_meas = True,
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
    grid_spot_spacing_um = 10.0, # µm
    min_xsp_spots_distance_um = 2.5, # µm
    randomize_frames = False,
    exclude_sample_margin = False
)

# =============================================================================
# Options for experimental standard collection
# =============================================================================
exp_stds_meas_cfg_kwargs = dict(
    min_acceptable_PB_ratio = 10,
    quant_flags_accepted = [0],
    use_for_mean_PB_calc = False
)

# =============================================================================
# Run
# =============================================================================
# Load microscope calibrations for this instrument and mode
calibs.load_microscope_calibrations(microscope_ID, measurement_mode, load_detector_channel_params=True)
eds_calibration_path = os.path.join(calibs.calibration_files_dir, cnst.SDD_CALIBS_MEAS_DIR, now_formatted)

# --- Acquire spectra
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
    fit_during_collection= fit_during_collection,
    update_std_library = update_std_library,
    is_auto_substrate_detection=is_auto_substrate_detection,
    auto_adjust_brightness_contrast=auto_adjust_brightness_contrast,
    contrast=contrast,
    brightness=brightness,
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
    exp_std_dir = eds_calibration_path
)

# --- Fit spectra and extract values of fitting parameters
sample_IDs = list(std_d['ID'] for std_d in std_list)
spectrum_IDs = 'all'
fit_params_vals_to_extract = list(f"{std_d['ref_el']}_{std_d['ref_peak']}_center" for std_d in std_list)
extracted_par_vals = batch_fit_spectra(sample_IDs,
                spectrum_IDs,
                is_standard = True,
                fit_params_vals_to_extract = fit_params_vals_to_extract,
                spectrum_lims = None,
                output_path = os.path.join(eds_calibration_path, 'Fitting output'),
                samples_path = eds_calibration_path,
                use_instrument_background = False,
                quantify_plot = False,
                plot_signal = False,
                zoom_plot = False,
                line_to_plot = '',
                els_substrate = els_substrate,
                fit_tol = 1e-4,
                is_particle = True,
                max_undetectable_w_fr = 0,
                force_single_iteration = False,
                interrupt_fits_bad_spectra = False,
                print_results = False,
                quant_verbose = True,
                fitting_verbose = False
)

# --- Calculate new SDD calibration values
meas_modes_calibs = calibs.detector_channel_params
current_energy_zero = meas_modes_calibs[measurement_mode][cnst.OFFSET_KEY]
current_bin_width = meas_modes_calibs[measurement_mode][cnst.SCALE_KEY]

# Extract mean measured energies for standards
measured_means = {}
for std in std_list:
    el = std['ref_el']
    ref_peak = std['ref_peak']
    param_name = f"{el}_{ref_peak}_center"
    sample_df = pd.DataFrame(extracted_par_vals[std['ID']])  # DataFrame for this sample

    meas_mean = sample_df.loc[sample_df['sp_id'] == 'mean', param_name].values[0]
    measured_means[param_name] = meas_mean

# Assign to calibration variables
x_measured_en = measured_means["Cu_Ka1_center"]
y_measured_en = measured_means["Al_Ka1_center"]

# Theoretical energies
x_th_en = get_el_xray_lines("Cu")["Ka1"]["energy (keV)"]
y_th_en = get_el_xray_lines("Al")["Ka1"]["energy (keV)"]

# Calculate new calibration
i_x = (x_measured_en - current_energy_zero) / current_bin_width
i_y = (y_measured_en - current_energy_zero) / current_bin_width

Dx = x_th_en - x_measured_en
Dy = y_th_en - y_measured_en

new_scale = (Dx - Dy + x_measured_en - y_measured_en) / (i_x - i_y)
new_offset = Dy + y_measured_en - i_y * new_scale

print(f"Current scale: {current_bin_width:.6f}")
print(f"Current offset: {current_energy_zero:.6f}")
print(f"New scale: {new_scale:.6f}")
print(f"New offset: {new_offset:.6f}")

# Add calibration file
calibs.update_detector_channel_params(measurement_mode, new_offset, new_scale)

    
    