#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration Loader for SEM EDS Analysis

This module loads the calibrations of a microscope: its detector model, from the package folder of the microscope
model (autoemx/calibrations/<microscope_ID>/), whose values are set as attributes of this module; and the detector
channel parameters and experimental standards of the instrument, from its local calibration folder,
~/.autoemx/calibrations/<microscope_ID>/ (see autoemx._paths).

Usage from other modules:
    import autoemx.calibrations as calibs
    calibs.load_microscope_calibrations(microscope_ID ='PhenomXL', meas_mode = 'point', load_detector_channel_params = True)

Author: Andrea Giunto
Created on: Mon Jan 20 15:40:42 2025
"""

import os
import sys
import copy
import json
from typing import Any, Optional
from datetime import datetime


import autoemx.utils.constants as cnst
import autoemx._paths as paths
from autoemx.config.schema_models import EDSStandardsFile
from autoemx.config.schema_models.detector_model import DetectorModel
from autoemx.utils.helper import print_single_separator
from autoemx.utils.legacy import is_legacy_standards_payload, standards_payload_to_model

from autoemx._logging import get_logger
logger = get_logger(__name__)


def ensure_standards_model(
    payload: Any,
    meas_type: str,
    beam_energy: int,
) -> EDSStandardsFile:
    """Normalize arbitrary standards payloads into the typed schema model.

    This keeps backward-compatibility conversion logic in one place for callers
    that may still hold legacy dict payloads in memory.
    """
    if isinstance(payload, EDSStandardsFile):
        return payload
    return standards_payload_to_model(
        payload=payload,
        measurement_type=meas_type,
        beam_energy_keV=int(beam_energy),
    )

microscope_calibrations_loaded = False
loaded_microscope_ID: Optional[str] = None
detector_model: Optional[DetectorModel] = None


def _loaded_microscope_ID() -> str:
    """ID of the microscope whose calibrations are loaded (needed to locate its local calibration folder)."""
    if loaded_microscope_ID is None:
        raise RuntimeError("Microscope calibrations not loaded: call load_microscope_calibrations() first.")
    return loaded_microscope_ID


def load_detector_model(microscope_ID: str) -> DetectorModel:
    """
    Detector model of the microscope, from detector_model.json in the package folder of the microscope model
    (autoemx/calibrations/<microscope_ID>/).

    Raises
    ------
    ValueError
        If the microscope has no calibration folder in the package, or no detector model.
    """
    package_dir = paths.package_calibration_dir(microscope_ID)
    if not package_dir.is_dir():
        raise ValueError(
            f"Could not find calibrations for microscope_ID '{microscope_ID}' at '{package_dir}'.\n"
            f"Please ensure microscope_ID matches a folder in '{paths.PACKAGE_CALIBRATIONS_DIR}'."
        )
    model_path = package_dir / cnst.DETECTOR_MODEL_FILENAME
    if not model_path.exists():
        raise ValueError(f"Could not find the detector model of microscope '{microscope_ID}' at '{model_path}'.")
    try:
        return DetectorModel.from_json_file(model_path)
    except Exception as e:
        raise ValueError(f"Could not load the detector model of microscope '{microscope_ID}' at '{model_path}':\n{e}") from e


def microscope_driver_name(microscope_ID: str) -> str:
    """Driver of the microscope: 'driver' in its detector model if set, otherwise the microscope ID."""
    model_path = paths.package_calibration_dir(microscope_ID) / cnst.DETECTOR_MODEL_FILENAME
    if model_path.exists():
        driver = load_detector_model(microscope_ID).driver
        if driver:
            return driver
    return microscope_ID


def load_microscope_calibrations(
    microscope_ID: str,
    meas_mode: Optional[str] = None,
    load_detector_channel_params: bool = False
) -> None:
    """
    Load the calibrations of a given microscope and EDS mode.

    The detector model of the microscope (see load_detector_model) is set as ``detector_model``, and each of its
    values and calibrated functions is set as an attribute of this module (e.g. ``calibs.elec_noise``,
    ``calibs.get_calibrated_peak_shape_params``).

    Optionally, it also loads the latest detector channel parameters. This is only needed during
    spectra acquisition, not when fitting or quantifying. In the latter cases, the detector channel
    parameters should be loaded from the measurement files directly.

    Parameters
    ----------
    microscope_ID : str
        The name of the microscope (must match a folder in autoemx/calibrations, with its detector model).
    meas_mode : str, Optional
        If provided, it checks whether it is included in the available_meas_modes of the detector model.
    load_detector_channel_params : bool, optional
        If True, loads the latest detector channel parameters as well (default: False).

    Raises
    ------
    ValueError
        If the calibrations or detector model cannot be found, or if the meas_mode is invalid.
    FileNotFoundError
        If no detector channel calibration file is found.
    """
    global microscope_calib_dir, loaded_microscope_ID, detector_model
    microscope_calib_dir = str(paths.package_calibration_dir(microscope_ID))
    model = load_detector_model(microscope_ID)

    # Validate the meas_mode
    if meas_mode is not None and meas_mode not in model.available_meas_modes:
        raise ValueError(
            f"Entered meas_mode '{meas_mode}' is not valid for microscope '{microscope_ID}'.\n"
            f"Available modes: {model.available_meas_modes}\n"
            "Change 'available_meas_modes' in the detector model of the microscope if this is undesired behavior."
        )

    # Values and calibrated functions of the detector model, as attributes of this module
    thismod = sys.modules[__name__]
    for name in DetectorModel.model_fields:
        if name not in ("schema_version", "peak_shape", "background"):
            setattr(thismod, name, getattr(model, name))
    thismod.get_calibrated_peak_shape_params = model.get_calibrated_peak_shape_params
    thismod.get_calibrated_background_params = model.get_calibrated_background_params
    detector_model = model

    global microscope_calibrations_loaded
    microscope_calibrations_loaded = True
    loaded_microscope_ID = microscope_ID

    # Optionally load detector channel parameters and check for the required meas_mode
    if load_detector_channel_params:
        load_latest_detector_channel_params(meas_mode)


def load_latest_detector_channel_params(meas_mode):
    """
    Load the latest detector channel parameters for a given measurement mode.

    This function retrieves the most recent detector channel calibration parameters
    (such as beam current, energy scale, and offset) for the specified measurement mode.
    It ensures that all required calibration keys are present before returning.

    The parameters are loaded via the `get_latest_detector_channel_params()` function,
    which must populate the global `detector_channel_params` dictionary.

    Parameters
    ----------
    meas_mode : str
        The measurement mode for which detector channel parameters are requested.
        Must be one of the modes available in `detector_channel_params`.

    Raises
    ------
    RuntimeError
        If `detector_channel_params` is not loaded into the global namespace.
    ValueError
        If the requested `meas_mode` is not present in `detector_channel_params`.
    KeyError
        If one or more required calibration keys are missing for the given `meas_mode`.
    """
    required_XS_calib_keys = [cnst.BEAM_CURRENT_KEY, cnst.SCALE_KEY, cnst.OFFSET_KEY]
    get_latest_detector_channel_params()
    if 'detector_channel_params' not in globals():
        raise RuntimeError("Failed to load detector_channel_params.")
    if meas_mode not in detector_channel_params:
        raise ValueError(
            f"meas_mode '{meas_mode}' not found in loaded detector_channel_params.\n"
            f"Available modes: {list(detector_channel_params.keys())}"
        )
    # Check that all required keys are present
    missing = [k for k in required_XS_calib_keys if k not in detector_channel_params[meas_mode]]
    if missing:
        raise KeyError(
            f"detector_channel_params for meas_mode '{meas_mode}' is missing required key(s): {missing}\n"
            f"Present keys: {list(detector_channel_params[meas_mode].keys())}"
        )    
        
        
def get_latest_detector_channel_params(verbose: bool = True) -> None:
    """
    Load dictionary of detector channel calibration parameters for each meas_mode.

    Each entry contains:
        - offset: float, energy offset (keV) for channel 0
        - scale: float, energy bin width (keV/channel)
        - spot_size: float, related to beam current

    These parameters should be recalibrated regularly.

    The most recent file in the local calibration folder of the microscope is used
    (~/.autoemx/calibrations/<microscope_ID>/Detector_channel_params_calibs).

    Raises
    ------
        FileNotFoundError: If no detector channel calibration file is found.
    """
    global calibration_files_dir
    calibration_files_dir = str(paths.user_calibration_dir(_loaded_microscope_ID()) / cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR)
    calib_file_path = paths.latest_detector_channel_params_file(_loaded_microscope_ID())
    if calib_file_path is None:
        raise FileNotFoundError(
            f"No detector channel parameter calibration file found in '{calibration_files_dir}'."
        )

    global detector_channel_params
    with open(calib_file_path, 'r') as file:
        detector_channel_params = json.load(file)

    if verbose:
        print_single_separator(f"ℹ️ Using detector calibration file '{calib_file_path}'")
        

def update_detector_channel_params(meas_mode, new_offset, new_scale, verbose: bool = True):
    """
    Update and save detector channel calibration parameters for a given measurement mode.

    This function retrieves the latest detector channel parameters for the specified
    measurement mode, updates the offset and scale values, and saves the updated
    parameters to a timestamped JSON file in the local calibration folder of the microscope.

    Parameters
    ----------
    meas_mode : str
        The measurement mode for which parameters will be updated.
        Must be present in the loaded `detector_channel_params`.
    new_offset : float
        The new detector channel offset value to set.
    new_scale : float
        The new detector channel scale value to set.
    verbose : bool, optional
        If True (default), prints/logs the location of the saved calibration file.

    Raises
    ------
    RuntimeError
        If detector channel parameters cannot be loaded.
    ValueError
        If the specified `meas_mode` does not exist in the parameters.
    KeyError
        If required calibration keys are missing from the parameters.
    """
    # Ensure latest parameters are loaded and valid
    load_latest_detector_channel_params(meas_mode)

    # Create an updated copy of parameters
    new_detector_channel_params = copy.deepcopy(detector_channel_params)
    new_detector_channel_params[meas_mode][cnst.SCALE_KEY] = round(new_scale, 6)
    new_detector_channel_params[meas_mode][cnst.OFFSET_KEY] = round(new_offset, 6)

    # Timestamp for filename
    now_str = datetime.now().strftime("%Y%m%d_%Hh%Mm")
    output_file_path = os.path.join(
        calibration_files_dir,
        f"{now_str}_{cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_FILENAME}.json"
    )

    # Save updated parameters
    with open(output_file_path, "w") as f:
        json.dump(new_detector_channel_params, f, indent=4)

    if verbose:
        print_single_separator(f"✅ Calibration saved to: {output_file_path}")

        

def load_standards_model(
    meas_type: str,
    beam_energy: int,
    std_f_dir: Optional[str] = None,
) -> tuple[EDSStandardsFile, str]:
    """
    Load standards data for a specified technique and beam energy.
    
    Called when performing quantifications.

    Parameters
    ----------
    meas_type : str
        Corresponds to MeasurementConfig.type
    beam_energy : int
        The beam energy (in keV) for which to load the standards.
    std_f_dir : str
        Directory of reference standards file. Default : None; uses the default directory
        Typically only specified when using a special reference file, such as the case of
        measuring the extent of precursor intermixing.

    Returns
    -------
    standards_file : EDSStandardsFile
        Typed standards model.
    standards_dir : str
        Path to the standards file that was loaded.
    
    Raises
    ------
    FileNotFoundError
        If the standards file does not exist.
    ValueError
        If the standards file cannot be parsed as JSON.
        
    Notes
    -----
    Expects a file named '{meas_type}_Stds_{beam_energy}keV.json' in std_f_dir or, by default, in the local
    calibration folder of the microscope (~/.autoemx/calibrations/<microscope_ID>/), where the file of the package
    is copied the first time it is needed.
    beam_energy will be converted to an int for the file name's purpose.
    """
    std_dict_filename = paths.standards_filename(meas_type, beam_energy)
    global standards_dir # So that it can be loaded from EMXSp_comp_analyser

    if std_f_dir is not None:
        standards_dir = os.path.join(std_f_dir, std_dict_filename)
        if not os.path.exists(std_f_dir): # Check if path exists
            logger.warning(f"⚠️ Warning: The provided path for reference standards {std_f_dir} does not exist.")
            std_f_dir = None
        elif not os.path.exists(standards_dir): # Check if reference file exists
            logger.warning(f"⚠️ Warning: The reference standards file {standards_dir} does not exist.")
            std_f_dir = None
        if std_f_dir is None:
            logger.info(f"ℹ️ Using standard reference file at {paths.user_calibration_dir(_loaded_microscope_ID())}")

    if std_f_dir is None:
        standards_dir = str(paths.user_standards_file(_loaded_microscope_ID(), meas_type, beam_energy))
    else:
        standards_dir = os.path.join(std_f_dir, std_dict_filename)
    try:
        with open(standards_dir, "r", encoding="utf-8") as standards_handle:
            standards_payload = json.load(standards_handle)

        payload_is_legacy = is_legacy_standards_payload(standards_payload)
        standards_file = standards_payload_to_model(
            payload=standards_payload,
            measurement_type=meas_type,
            beam_energy_keV=int(beam_energy),
        )

        if payload_is_legacy:
            backup_path = f"{standards_dir}_OLD"
            if os.path.exists(backup_path):
                idx = 1
                while os.path.exists(f"{backup_path}_{idx}"):
                    idx += 1
                backup_path = f"{backup_path}_{idx}"
            os.rename(standards_dir, backup_path)
            standards_file.to_json_file(standards_dir, indent=2)
            logger.info(
                f"ℹ️ Migrated legacy standards file to schema format. "
                f"Legacy backup saved at '{backup_path}'."
            )
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Could not find the standards file for beam energy {beam_energy} keV.\n"
            f"Tried to open: {standards_dir}"
        ) from e
    except json.JSONDecodeError as e:
        raise ValueError(
            f"Could not parse the standards JSON file for beam energy {beam_energy} keV.\n"
            f"File path: {standards_dir}\n"
            f"Error: {e}"
        ) from e
    except Exception as e:
        raise ValueError(
            f"Could not validate the standards JSON file for beam energy {beam_energy} keV.\n"
            f"File path: {standards_dir}\n"
            f"Error: {e}"
        ) from e
    return standards_file, standards_dir


def load_standards(meas_type: str, beam_energy: int, std_f_dir : Optional[str] = None) -> dict:
    """Backward-compatible standards loader returning legacy dict payloads.

    Prefer ``load_standards_model`` for model-first callers.
    """
    standards_file, standards_dir = load_standards_model(
        meas_type=meas_type,
        beam_energy=beam_energy,
        std_f_dir=std_f_dir,
    )
    standards = standards_file.to_standards_dict()
    return standards, standards_dir


