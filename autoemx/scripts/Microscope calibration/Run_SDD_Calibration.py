#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated SDD Detector Calibration

This script automates the collection and fitting of X-ray spectra from experimental standards
to generate an SDD calibration file for AutoEMX.

Requirements:
    - Proper instrument calibration files and instrument driver for the selected microscope

Typical Usage:
    - Edit the 'std_list' list to define your standards.
        - ID: A preferred ID is fine.
        - formula: Composition of standard, for accurate spectral fitting.
        - ref_el: Element in the composition whose peak should be taken as reference to calibrate the SDD.
        - ref_peak: Characteristic X-ray to take as reference to calibrate the SDD (e.g., Ka1, La1, Ma1, Mz1).
        - pos: Position of standard in the microscope stage.
        - sample_type: 'bulk' or 'powder'. Check autoemx.config.runtime_configs.SampleConfig for updated support.
        - is_manual_meas: Set to True to manually select spots to measure.

    - Suggestions:
        - Preferably use bulk standards. Powder standards are also acceptable if bulk standards are not available.
        - Choose peaks above 2 keV, well-distanced between themselves.
        - Exactly two standards are required (two-point energy calibration).

    - Set reuse_session to None to acquire new spectra, or to an existing session folder
      name under SDD calibrations/ to skip acquisition and reuse prior ledgers.

    - Run the script to collect and fit spectra, with automated generation of SDD calibration file.

Potential improvements:
    - At the moment the script may run spectral fitting twice when fit_during_collection is True
      and centers are still re-read from the ledger, which is inefficient but acceptable given
      that this code is not run very often.

Created on Fri Aug 20 09:34:34 2025

@author: Andrea
"""
import os
from datetime import datetime

import numpy as np

from autoemx.runners.batch_acquire_experimental_stds import batch_acquire_experimental_stds
from autoemx.runners.batch_fit_spectra import batch_fit_spectra
from autoemx.data.Xray_lines import get_el_xray_lines
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.runtime_configs import SampleConfig
import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.utils import print_double_separator


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

# Set to an existing session folder name under SDD calibrations/ to skip acquisition
# (e.g. '20260521_09h22m'). None acquires new spectra into a timestamped session folder.
reuse_session = None

# =============================================================================
# Sample Definitions - Add two pure elements to measure
# =============================================================================
Cu_center = (37.863, 38.195)
std_list = [
    {'ID': 'Cu', 'formula': 'Cu', 'ref_el': 'Cu', 'ref_peak': 'Ka1', 'pos': Cu_center, 'sample_type': 'bulk', 'is_manual_meas': False},
    {'ID': 'Al', 'formula': 'Al', 'ref_el': 'Al', 'ref_peak': 'Ka1', 'pos': tuple(a + b for a, b in zip(Cu_center, (5, 0))), 'sample_type': 'bulk', 'is_manual_meas': False},
]

if len(std_list) != 2:
    raise ValueError("SDD energy calibration requires exactly two standards in std_list.")

# =============================================================================
# Acquisition Options and Sample description
working_distance = 5.5  # mm
is_auto_substrate_detection = False

is_manual_meas = False  # Default navigation mode; per-sample 'is_manual_meas' overrides this.

fit_during_collection = False
update_std_library = False

sample_substrate_type = 'None'
sample_substrate_shape = 'circle'
sample_halfwidth = 1  # mm

measurement_mode = 'point'
beam_energy = 15  # keV

auto_adjust_brightness_contrast = True
contrast = None  # 4.3877  # Used if auto_adjust_brightness_contrast = False
brightness = None  # 0.4504  # Used if auto_adjust_brightness_contrast = False
saved_images_extension = 'png'  # Default lightweight output. Set 'tif' for higher-resolution/larger files.
save_raw_images = False  # Default saves only the annotated image. Set True to also save raw images.

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
    par_selection_mode='auto',
    is_known_powder_mixture_meas=True,
    img_shift_tracking=True,
    max_n_par_per_frame=30,
    max_spectra_per_par=3,
    max_area_par=10000.0,
    min_area_par=10.0,
    par_mask_margin=1.0,
    xsp_spots_distance_um=1.0,
    par_brightness_thresh=100,
    par_xy_spots_thresh=100,
    par_feature_selection='peaks',
    par_spot_spacing='random',
)

# =============================================================================
# Bulk options
# =============================================================================
bulk_meas_cfg_kwargs = dict(
    grid_spot_spacing_um=10.0,  # µm
    min_xsp_spots_distance_um=2.5,  # µm
    image_frame_width_um=None,  # µm
    randomize_frames=False,
    exclude_sample_margin=False,
)

# =============================================================================
# Options for experimental standard collection
# =============================================================================
exp_stds_meas_cfg_kwargs = dict(
    min_acceptable_PB_ratio=10,
    quant_flags_accepted=[0],
    els_to_use_for_mean_PB_calc=["none"],
    generate_separate_std_dict=False,
)

# =============================================================================
# Run
# =============================================================================
# Load microscope calibrations for this instrument and mode
calibs.load_microscope_calibrations(microscope_ID, measurement_mode, load_detector_channel_params=True)

session_name = reuse_session if reuse_session else now_formatted
eds_calibration_path = os.path.join(
    calibs.calibration_files_dir,
    cnst.SDD_CALIBS_MEAS_DIR,
    session_name,
)

if reuse_session:
    print_double_separator()
    print(f"Reusing existing SDD calibration session: {eds_calibration_path}")
    std_paths = [os.path.join(eds_calibration_path, std['ID']) for std in std_list]
    missing = [path for path in std_paths if not os.path.isdir(path)]
    if missing:
        raise FileNotFoundError(
            "reuse_session is set but the following standard folders are missing:\n"
            + "\n".join(missing)
        )
else:
    # --- Acquire and save spectra
    analyzers = batch_acquire_experimental_stds(
        stds=std_list,
        microscope_ID=microscope_ID,
        microscope_type=microscope_type,
        measurement_type=measurement_type,
        measurement_mode=measurement_mode,
        sample_halfwidth=sample_halfwidth,
        sample_substrate_type=sample_substrate_type,
        sample_substrate_shape=sample_substrate_shape,
        working_distance=working_distance,
        beam_energy=beam_energy,
        spectrum_lims=spectrum_lims,
        use_instrument_background=use_instrument_background,
        min_bckgrnd_cnts=min_bckgrnd_cnts,
        is_manual_meas=is_manual_meas,
        fit_during_collection=fit_during_collection,
        update_std_library=update_std_library,
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
        exp_std_dir=eds_calibration_path,
    )
    if any(an is None for an in analyzers):
        failed = [std['ID'] for std, an in zip(std_list, analyzers) if an is None]
        raise RuntimeError(f"Acquisition failed for standard(s): {failed}")
    std_paths = [an.sample_result_dir for an in analyzers]


def ledger_has_fit_results(ledger, ref_peak):
    for spectrum in getattr(ledger, 'spectra', []):
        for quant_result in getattr(spectrum, 'quantification_results', []):
            fit_result = getattr(quant_result, 'fit_result', None)
            if fit_result and ref_peak in getattr(fit_result, 'fitted_peaks', {}):
                return True
    return False


extracted_par_vals = {}
fit_params_vals_to_extract = [f"{std_d['ref_el']}_{std_d['ref_peak']}_center" for std_d in std_list]

for std, sample_path in zip(std_list, std_paths):
    sample_ID = std['ID']
    ref_peak = f"{std['ref_el']}_{std['ref_peak']}"
    is_particle = std['sample_type'] in SampleConfig.POWDER_SAMPLES_TYPES
    ledger_path = os.path.join(sample_path, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
    if not os.path.exists(ledger_path):
        raise FileNotFoundError(f"Ledger not found for standard '{sample_ID}': {ledger_path}")
    ledger = load_sample_ledger(ledger_path)
    if not ledger_has_fit_results(ledger, ref_peak):
        batch_fit_spectra(
            [sample_ID],
            'all',
            is_standard=True,
            fit_params_vals_to_extract=fit_params_vals_to_extract,
            spectrum_lims=spectrum_lims,
            samples_path=eds_calibration_path,
            use_instrument_background=use_instrument_background,
            plot_signal=False,
            zoom_plot=False,
            line_to_plot='',
            els_substrate=els_substrate,
            fit_tol=1e-4,
            is_particle=is_particle,
            max_undetectable_w_fr=0,
            force_single_iteration=False,
            interrupt_fits_bad_spectra=False,
            print_results=False,
            quant_verbose=True,
            fitting_verbose=False,
        )
        ledger = load_sample_ledger(ledger_path)

    # Collect fit parameter values for this sample, only if quant_flag == 0
    centers = []
    for spectrum in getattr(ledger, 'spectra', []):
        for quant_result in getattr(spectrum, 'quantification_results', []):
            if getattr(quant_result, 'quant_flag', None) == 0:
                fit_result = getattr(quant_result, 'fit_result', None)
                if fit_result and ref_peak in getattr(fit_result, 'fitted_peaks', {}):
                    center = fit_result.fitted_peaks[ref_peak].center
                    if center is not None:
                        centers.append(center)
    extracted_par_vals[sample_ID] = centers


# --- Calculate new SDD calibration values
meas_modes_calibs = calibs.detector_channel_params
current_energy_zero = meas_modes_calibs[measurement_mode][cnst.OFFSET_KEY]
current_bin_width = meas_modes_calibs[measurement_mode][cnst.SCALE_KEY]

# Extract mean measured energies for standards
print_double_separator()
measured_means = {}
for std in std_list:
    el = std['ref_el']
    ref_peak = std['ref_peak']
    param_name = f"{el}_{ref_peak}_center"
    centers = extracted_par_vals[std['ID']]
    valid_centers = [c for c in centers if c is not None and not np.isnan(c)]
    if len(valid_centers) > 0:
        meas_mean = float(np.mean(valid_centers))
        print(f"Center of {std['ID']} was calibrated over {len(valid_centers)} spectra")
    else:
        raise RuntimeError(f"No valid fitted center values found for standard '{std['ID']}' ({param_name})")
    measured_means[param_name] = meas_mean

# Two-point calibration from std_list order (first = x, second = y)
std_x, std_y = std_list
x_param = f"{std_x['ref_el']}_{std_x['ref_peak']}_center"
y_param = f"{std_y['ref_el']}_{std_y['ref_peak']}_center"
x_measured_en = measured_means[x_param]
y_measured_en = measured_means[y_param]

# Theoretical energies
x_th_en = get_el_xray_lines(std_x['ref_el'])[std_x['ref_peak']]["energy (keV)"]
y_th_en = get_el_xray_lines(std_y['ref_el'])[std_y['ref_peak']]["energy (keV)"]

# Calculate new calibration
i_x = (x_measured_en - current_energy_zero) / current_bin_width
i_y = (y_measured_en - current_energy_zero) / current_bin_width

Dx = x_th_en - x_measured_en
Dy = y_th_en - y_measured_en

new_scale = (Dx - Dy + x_measured_en - y_measured_en) / (i_x - i_y)
new_offset = Dy + y_measured_en - i_y * new_scale

print_double_separator()
print(f"Current scale: {current_bin_width:.6f}")
print(f"Current offset: {current_energy_zero:.6f}")
print(f"New scale: {new_scale:.6f}")
print(f"New offset: {new_offset:.6f}")

# Add calibration file
calibs.update_detector_channel_params(measurement_mode, new_offset, new_scale)
