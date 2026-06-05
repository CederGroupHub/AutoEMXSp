#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration Loader for SEM EDS Analysis

This module provides a function to dynamically load and inject
microscope-specific calibration parameters and functions into the caller's namespace.

Usage from other modules:
    import XSp_calibs as calibs 
    calibs.load_microscope_calibrations(microscope_ID ='PhenomXL', meas_mode = 'point', load_detector_channel_params = True)

Author: Andrea Giunto
Created on: Mon Jan 20 15:40:42 2025
"""

import os
import sys
import importlib
import json
from typing import Any, Optional
from datetime import datetime


import autoemx.utils.constants as cnst
from autoemx.config.schema_models import EDSStandardsFile
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
def load_microscope_calibrations(
    microscope_ID: str,
    meas_mode: Optional[str] = None,
    load_detector_channel_params: bool = False
) -> None:
    """
    Dynamically load calibration parameters for a given microscope and EDS mode.

    This function imports all public attributes from the calibration module
    corresponding to the specified microscope_ID and injects them into the
    current module's namespace.
    
    Optionally, it also loads the latest detector channel parameters. This is only needed during
    spectra acquisition, not when fitting or quantifying. In the latter cases, the detector channel
    parameters should be loaded from the measurement files diurectly.


    Parameters
    ----------
    microscope_ID : str
        The name of the microscope (must match a folder in the calibration directory).
    meas_mode : str, Optional
        If provided, it checks whether it is included in the available_meas_modes.
        The EDS mode to use (must be listed in 'available_meas_modes' in the calibration module).
    load_detector_channel_params : bool, optional
        If True, loads the latest detector channel parameters as well (default: False).

    Raises
    ------
    ValueError
        If the calibration directory or module cannot be found, or if the meas_mode is invalid.
    AttributeError
        If the calibration module does not define 'available_meas_modes'.
    FileNotFoundError
        If no detector channel calibration file is found.

    Warning
    -------
    This function injects variables and functions into the module namespace.
    Use with care to avoid name collisions.
    """
    # Build the path and check for microscope calibration directory
    global microscope_calib_dir
    microscope_calib_dir = os.path.join(
        os.path.dirname(__file__),
        microscope_ID
    )
    if not os.path.isdir(microscope_calib_dir):
        raise ValueError(
            f"Could not find the microscope calibration folder at '{microscope_calib_dir}'.\n"
            f"Please ensure microscope_ID ('{microscope_ID}') matches a folder in '{os.path.dirname(__file__)}'."
        )

    # Build the module name for importlib
    module_name = f".{microscope_ID}.XS_calibrations"
    try:
        mod = importlib.import_module(module_name, package=__name__)
    except ModuleNotFoundError as e:
        raise ValueError(
            f"Could not find the calibration module for microscope_ID '{microscope_ID}'.\n"
            f"Tried to import: {module_name} (relative to package '{__name__}')."
        ) from e

    # Inject all public attributes from the calibration module into this module
    thismod = sys.modules[__name__]
    for k in dir(mod):
        if not k.startswith('_'):
            setattr(thismod, k, getattr(mod, k))
    
    # Validate the meas_mode
    if meas_mode is not None:
        if not hasattr(mod, "available_meas_modes"):
            raise AttributeError(
                f"The calibration module for '{microscope_ID}' does not define 'available_meas_modes'."
            )
        if meas_mode not in mod.available_meas_modes:
            raise ValueError(
                f"Entered meas_mode '{meas_mode}' is not valid for microscope '{microscope_ID}'.\n"
                f"Available modes: {mod.available_meas_modes}\n"
                "Change values in the microscope calibration file if this is undesired behavior."
            )
    
    global microscope_calibrations_loaded
    microscope_calibrations_loaded = True
    
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

    Raises
    ------
        FileNotFoundError: If no detector channel calibration file is found.
    """
    import os
    import json
    
    global calibration_files_dir
    calibration_files_dir = os.path.join(microscope_calib_dir, cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR)
    calibration_files = [f for f in os.listdir(calibration_files_dir) if f'{cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_FILENAME}.json' in f]
    if not calibration_files:
        raise FileNotFoundError(
            f"No detector channel parameter calibration file found in '{calibration_files_dir}'."
        )
    calib_file = sorted(calibration_files)[-1]
    calib_file_dir = os.path.join(calibration_files_dir, calib_file)
    
    global detector_channel_params
    with open(calib_file_dir, 'r') as file:
        detector_channel_params = json.load(file)
        
    if verbose:
        print_single_separator()
        logger.info(f"ℹ️ Using detector calibration file '{calib_file}'")
        

def update_detector_channel_params(meas_mode, new_offset, new_scale, verbose: bool = True):
    """
    Update and save detector channel calibration parameters for a given measurement mode.

    This function retrieves the latest detector channel parameters for the specified
    measurement mode, updates the offset and scale values, and saves the updated
    parameters to a timestamped JSON file in the calibration directory.

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
    new_detector_channel_params = detector_channel_params.copy()
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
        print_single_separator()
        logger.info(f"✅ Calibration saved to: {output_file_path}")

        

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
    Expects a file named '{meas_type}_{beam_energy}keV.json' in the directory specified
    by `microscope_calib_dir`.
    beam_energy will be converted to an int for the file name's purpose.
    """
    std_dict_filename = f'{meas_type}_{cnst.STD_FILENAME}_{int(beam_energy):d}keV.json'
    global standards_dir # So that it can be loaded from EMXSp_comp_analyser

    if std_f_dir is not None:
        standards_dir = os.path.join(std_f_dir, std_dict_filename)
        if not os.path.exists(std_f_dir): # Check if path exists
            logger.warning(f"⚠️ Warning: The provided path for reference standards {std_f_dir} does not exist.")
            logger.info(f"ℹ️ Using standard reference file at {microscope_calib_dir}")
            std_f_dir = microscope_calib_dir
        elif not os.path.exists(standards_dir): # Check if reference file exists
            logger.warning(f"⚠️ Warning: The reference standards file {standards_dir} does not exist.")
            logger.info(f"ℹ️ Using standard reference file at {microscope_calib_dir}")
            std_f_dir = microscope_calib_dir
    else:
        std_f_dir = microscope_calib_dir
    
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
            backup_path = f"{standards_dir})OLD"
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


