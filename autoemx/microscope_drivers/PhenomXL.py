#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Electron Microscope Driver
==========================

Hardware driver template for controlling an Electron Microscope.

This reference implementation is based on **PyPhenom 2.0**, but structured for 
adaptation to other microscope APIs by replacing `PyPhenom` calls with equivalents.

Purpose
-------
- Control microscope hardware for `autoemx` workflows.
- Provide functions for stage movement, SEM imaging, EDS acquisition,
  auto-focus/contrast, beam parameter control, and navigation camera operations.
- Serve as a documented reference for developers adapting to other APIs.

IMPORTANT
---------
If replacing with another API:
    - Change name of this module so that it corresponds with the employed MicroscopeConfig.ID (autoemx/config/runtime_configs.py),
      or set it as "driver" in the detector_model.json of the microscopes sharing it (autoemx/calibrations/<ID>/)
    - Maintain function signatures and returned value formats for compatibility.
    - Units in parameters/returns must be consistent across implementations.
    - Error handling should remain explicit and informative (raise `EMError`).


This code has been tested with the following microscope:
- **Phenom User Interface Version**: 2.0.2-rel
- **Phenom Type**: PhenomXL
- **Phenom Software Version**: 6.8.4.rel.e2fd11083c.28328

Created: Tue Jul 29 13:18:16 2025  
Author: Andrea

PyPhenom version: 
"""

import numpy as np
import time
import os
import cv2
import warnings
import traceback
from typing import Optional, Tuple, List
from autoemx.utils import EMError

from autoemx._logging import get_logger
logger = get_logger(__name__)

# =============================================================================
# Stage Physical Limits (example values — adapt to your microscope)
# =============================================================================
# These define the movement range of the microscope stage in its own reference system.
# They are also used when transforming coordinates from image pixels to absolute stage positions.

stage_x_left: float   = -40  # mm, leftmost limit of the stage
stage_x_right: float  =  40  # mm, rightmost limit of the stage
stage_y_top: float    =  40  # mm, topmost limit of the stage
stage_y_bottom: float = -40  # mm, bottommost limit of the stage

# Transformation array to convert image pixel shifts to stage coordinate shifts.
# Encodes axis directionality based on stage limits.
image_to_stage_coords_transform: np.ndarray = np.array([
    np.sign(stage_x_right - stage_x_left),   # X axis: +1 if right > left, else -1
    np.sign(stage_y_bottom - stage_y_top)    # Y axis: +1 if bottom > top, else -1
])

# =============================================================================
# Navigation Camera Physical Parameters
# =============================================================================
navcam_im_w_mm: float = 98
# Width of the navigation camera image, in millimeters.
# This is used to convert between navcam pixel dimensions and real-world stage position.

# Image offsets — apply if misalignment exists between navigation camera and SEM.
# Otherwise set to 0.
navcam_x_offset: int = 2  # X offset (pixels) between navigation camera and SEM alignment.
navcam_y_offset: int = 2  # Y offset (pixels) between navigation camera and SEM alignment.

# =============================================================================
# SEM Acquisition Parameters
# =============================================================================
typical_wd: float = 5
# Typical working distance used in SEM, in millimeters.

im_width: int = 1920
im_height: int = 1200
# SEM image default width and height, in pixels.

# =============================================================================
# EDS Parameters
# =============================================================================
RAW_SPECTRUM_EXT: str = '.msa'  # Format used by the EDS software to export spectra

# =============================================================================
# Microscope Connection
# =============================================================================
# Delay API import/connection until microscope operations are actually requested.
ppi = None
phenom = None
acqScanParams = None
is_at_EM: bool = False


def is_microscope_connected() -> bool:
    """Return True when driver state indicates a live microscope connection."""
    return phenom is not None and ppi is not None and acqScanParams is not None


def connect_to_microscope(warn_if_unavailable: bool = True) -> bool:
    """Connect to the microscope API on demand.

    Returns
    -------
    bool
        True when connected, False when unavailable.
    """
    global ppi, phenom, acqScanParams, is_at_EM

    if is_microscope_connected():
        is_at_EM = True
        return True

    try:
        import PyPhenom as ppi_module

        ppi = ppi_module
        phenom = ppi.Phenom()
        acqScanParams = ppi.ScanParams()
        acqScanParams.size = ppi.Size(im_width, im_height)
        acqScanParams.detector = ppi.DetectorMode.All  # Preferably use backscattered detector
        acqScanParams.nFrames = 1  # Limit frames to 1 to speed acquisition
        acqScanParams.hdr = False
        acqScanParams.scale = 1.0
        is_at_EM = is_microscope_connected()
        logger.info("Microscope driver connected successfully.")
        return True
    except Exception as e:
        ppi = None
        phenom = None
        acqScanParams = None
        is_at_EM = False
        if warn_if_unavailable:
            tb = traceback.format_exc()
            warnings.warn(f"Microscope driver not available: {e}.\nTraceback:\n{tb}")
        return False

# =============================================================================
# Coordinate Conversion Functions
# =============================================================================
'''
Phenom Coordinate System
----------------
The coordinates are expressed in a normalized, aspect-ratio-correct system centered at the image center:

    - The origin (0, 0) is at the image center.
    - The x-axis is horizontal, increasing to the right, ranging from -0.5 (left) to +0.5 (right).
    - The y-axis is vertical, increasing downward, and scaled by the aspect ratio (height/width):
        * Top edge:    y = -0.5 × (height / width)
        * Bottom edge: y = +0.5 × (height / width)
    
    |        (-0.5, -0.5*height/width)         (0.5, -0.5*height/width)
    |                       +-------------------------+
    |                       |                         |
    |                       |                         |
    |                       |           +(0,0)        |-----> +x
    |                       |                         |
    |                       |                         |
    v  +y                   +-------------------------+
            (-0.5,  0.5*height/width)         (0.5, 0.5*height/width)

This ensures the coordinate system is always centered and aspect-ratio-correct, regardless of image size.
'''
    
def frame_rel_to_pixel_coords(
    pts_rel_coords: np.ndarray,
    img_width: int,
    img_height: int
) -> np.ndarray:
    """
    Convert relative coordinates (normalized by the EM API) into pixel coordinates.

    See inline comments for a detailed description of the coordinate system.

    Parameters
    ----------
    pts_rel_coords : np.ndarray
        Shape (N, 2) array with (x_rel, y_rel) coordinates.
        - x_rel in [-0.5, 0.5]
        - y_rel in [-0.5 * (H / W), 0.5 * (H / W)]
    img_width : int
        Image width in pixels.
    img_height : int
        Image height in pixels.

    Returns
    -------
    np.ndarray
        Shape (N, 2) array of (x_px, y_px) pixel coordinates.
    """
    # Coordinate space is normalized, aspect-ratio corrected, and centered at the image centre.
    # Origin (0,0) = image centre.
    # x-axis runs horizontally, increasing to the right.
    # y-axis runs vertically, increasing downward, scaled by aspect ratio (H/W).
    pts_rel_coords = np.atleast_2d(np.asarray(pts_rel_coords, dtype=float))
    aspect_ratio = img_height / img_width
    x_px = (pts_rel_coords[:, 0] + 0.5) * img_width
    y_px = (pts_rel_coords[:, 1] / aspect_ratio + 0.5) * img_height
    return np.column_stack((x_px, y_px))


def frame_pixel_to_rel_coords(
    pts_pixel_coords: np.ndarray,
    img_width: int,
    img_height: int
) -> np.ndarray:
    """
    Convert pixel coordinates into relative normalized coordinates compatible with EM API.

    See inline comments for a detailed description of the coordinate system.

    Parameters
    ----------
    pts_pixel_coords : np.ndarray
        Shape (N, 2) array with (x_px, y_px) pixel coordinates.
    img_width : int
        Image width in pixels.
    img_height : int
        Image height in pixels.

    Returns
    -------
    np.ndarray
        Shape (N, 2) array of normalized (x_rel, y_rel) coordinates.
    """
    pts_pixel_coords = np.atleast_2d(np.asarray(pts_pixel_coords, dtype=float))
    aspect_ratio = img_height / img_width
    x_rel = pts_pixel_coords[:, 0] / img_width - 0.5
    y_rel = (pts_pixel_coords[:, 1] / img_height - 0.5) * aspect_ratio
    return np.column_stack((x_rel, y_rel))

# =============================================================================
# SEM Operational Controls
# =============================================================================
def standby() -> None:
    """Put the SEM into standby mode."""
    return phenom.Standby()

def set_electron_detector_mode(detector_name: str) -> None:
    """
    Set the electron detector mode.

    Notes
    -----
    * Accepted values of detector names are set at `autoemx/config/runtime_configs.py`
      within `MicroscopeConfig.detector_type`. Ensure `ALLOWED_DETECTOR_TYPES`
      is updated when new modes are added.
    * Current version passes the value "BSD", but PyPhenom requires "All".
      This mapping is performed internally.
    """
    # Map "BSD" to "All" for PyPhenom compatibility.
    if detector_name == 'BSD':
        detector_name = 'All'
    viewingMode = phenom.GetSemViewingMode()
    viewingMode.scanParams.detector = getattr(ppi.DetectorMode, detector_name)
    phenom.SetSemViewingMode(viewingMode)

def _get_instrument_mode():
    """Return the current instrument mode (PyPhenom enum)."""
    return phenom.GetInstrumentMode()

def activate() -> None:
    """
    Wake up the instrument from standby mode.
    Required before sending commands if the SEM is inactive.
    """
    if _get_instrument_mode() == ppi.InstrumentMode(2):  # Checks if SEM is in Standby
        phenom.Activate()

def to_SEM(timeout: int = 120) -> None:
    """
    Switch the instrument to SEM mode.

    Parameters
    ----------
    timeout : int
        Maximum time in seconds to wait for SEM mode activation.

    Raises
    ------
    TimeoutError
        If SEM mode is not reached within `timeout` seconds.
    """
    start = time.time()
    while True:
        if time.time() - start > timeout:
            raise TimeoutError(f"SEM did not return to LiveSem within {timeout}s.")
        try:
            phenom.MoveToSem()
            break
        except Exception:
            wait_time = 5
            logger.warning(f"⚠️ Phenom busy. Retrying in {wait_time}s...")
            time.sleep(wait_time)
    time.sleep(1)  # Workaround for API timing

# =============================================================================
# SEM Beam Controls
# =============================================================================
def set_high_tension(voltage_keV: float) -> None:
    """
    Set the SEM acceleration voltage.

    Parameters
    ----------
    voltage_keV : float
        Positive value in kiloelectronvolts (keV).
    """
    # PyPhenom expects volts and an inverted sign convention.
    voltage_V = - voltage_keV * 1000
    phenom.SetSemHighTension(voltage_V)


def set_beam_current(current: float) -> None:
    """
    Set the SEM beam current.

    Parameters
    ----------
    current : float
        Beam current value. Units depend on microscope API requirements.

    Notes
    -----
    Beam current is passed via:
        `autoemx/config/runtime_configs.py`, within `MeasurementConfig.beam_current`,  
        or predefined in:
        `~/.autoemx/calibrations/YOUR_MICROSCOPE/Detector_channel_params_calibs/DATE_detector_channel_params_calibs.json`  
        using key `autoemx.utils.constants.BEAM_CURRENT_KEY`.
    """
    phenom.SetSemSpotSize(current)

# =============================================================================
# EDS / Spectroscopy Functions
# =============================================================================
def get_EDS_analyser_object():
    """
    Retrieve an EDS analyser object for spectral acquisition.

    Returns
    -------
    object
        EDS Job Analyzer instance.

    Notes
    -----
    This analyzer should be passed to `acquire_XS_spectral_data()` for X-ray acquisition.
    """
    # Load API pulse processor settings.
    settings = phenom.LoadPulseProcessorSettings()
    analyzer = ppi.Application.ElementIdentification.EdsJobAnalyzer(phenom)
    analyzer.preset = settings.spot
    return analyzer


def acquire_XS_spectral_data(
    analyzer,
    x: float,
    y: float,
    max_collection_time: float,
    target_counts: int,
    elements: Optional[List[str]] = None,
    msa_file_export_path = None
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[float], Optional[float]]:
    """
    Acquire an EDS spectrum (and optionally background) at given coordinates.

    Parameters
    ----------
    analyzer : object
        EDS Analyzer object from `get_EDS_analyser_object()`.
    x, y : float
        Coordinates in API units (typically mm).
    max_collection_time : float
        Maximum acquisition time in seconds.
    target_counts : int
        Target X-ray counts for acquisition.
    elements : list of str, optional
        Element symbols for quantification. If provided, a background spectrum
        will be quantified using these elements.
    msa_file_export_path: str | None, optional
        If a path is provided, it exports a .msa spectral file with all the metadata. Default: None

    Returns
    -------
    tuple
        (spectrum_data, background_data, real_time, live_time):
            spectrum_data : np.ndarray or None
                Spectrum data array.
            background_data : np.ndarray or None
                Quantified background intensities.
            real_time : float or None
                Actual time measurement took (seconds).
            live_time : float or None
                Detector live time during measurement (seconds).

    Raises
    ------
    EMError
        If acquisition or quantification fails.
    """
    phenom_spectrum = _acquire_spot_spectrum(
        analyzer,
        x,
        y,
        max_collection_time,
        target_counts,
    )

    background_data = _try_extract_instrument_background(phenom_spectrum, elements)

    spectrum_data = phenom_spectrum.spectrum.data
    real_time = phenom_spectrum.metadata.realTime
    live_time = phenom_spectrum.metadata.liveTime

    _export_msa_if_requested(phenom_spectrum, msa_file_export_path)

    return spectrum_data, background_data, real_time, live_time


def _acquire_spot_spectrum(
    analyzer,
    x: float,
    y: float,
    max_collection_time: float,
    target_counts: int,
):
    """Acquire one spot spectrum and return the proprietary spectrum object."""
    try:
        spot_data = analyzer.AddSpot(
            ppi.Position(x, y),
            maxTime=max_collection_time,
            maxCounts=target_counts,
        )
    except Exception as er:
        raise EMError(
            f"The following error was encountered during X-Ray Spectrum acquisition:\n{er}\nSpectrum not recorded."
        )

    analyzer.Wait()
    return spot_data.spotSpectrum


def _try_extract_instrument_background(phenom_spectrum, elements: Optional[List[str]]) -> Optional[np.ndarray]:
    """Attempt proprietary background estimation; warn and return None on failure."""
    if not elements:
        return None

    try:
        ppi_elements = [getattr(ppi.Spectroscopy.Element, el) for el in elements]
    except Exception as er:
        warnings.warn(
            "Could not map requested elements for instrument background estimation; "
            f"falling back to automatic background subtraction. Details: {er}"
        )
        return None

    try:
        phenom_background = ppi.Spectroscopy.Quantify(phenom_spectrum, ppi_elements).background
        return phenom_background.data
    except Exception as er:
        warnings.warn(
            "Instrument background estimation failed in microscope software; "
            f"falling back to automatic background subtraction. Details: {er}"
        )
        return None


def _export_msa_if_requested(phenom_spectrum, msa_file_export_path: Optional[str]) -> None:
    """Write the proprietary .msa file when a valid export path is provided."""
    if msa_file_export_path is None:
        return

    base_dir = os.path.dirname(msa_file_export_path)
    if base_dir:
        os.makedirs(base_dir, exist_ok=True)

    if os.path.exists(base_dir):
        ppi.Spectroscopy.WriteMsaFile(phenom_spectrum, msa_file_export_path)
    else:
        warnings.warn(
            "MSA file could not be exported because the provided path does not exist: "
            f"{base_dir}"
        )
    
# =============================================================================
# Focus and Image Adjustment
# =============================================================================
def auto_focus() -> float:
    """
    Automatically optimize SEM focus.

    Returns
    -------
    float
        Current working distance in millimeters after focus adjustment.
    """
    # Small pause helps autofocus reliability.
    time.sleep(0.2)
    phenom.SemAutoFocus()
    # Return updated working distance.
    return get_current_wd()


def auto_contrast_brightness() -> None:
    """
    Automatically optimize SEM contrast and brightness.
    """
    phenom.SemAutoContrastBrightness()


def adjust_focus(new_wd: float) -> None:
    """
    Set a new working distance.

    Parameters
    ----------
    new_wd : float
        Desired working distance in millimeters.
    """
    phenom.SetSemWD(new_wd * 0.001)

# =============================================================================
# Stage Motion
# =============================================================================
def move_to(x: float, y: float) -> None:
    """
    Move stage to an absolute position.

    Parameters
    ----------
    x : float
        X coordinate in millimeters.
    y : float
        Y coordinate in millimeters.
    """
    phenom.MoveTo(x * 0.001, y * 0.001)
    

def move_by(x: float, y: float) -> None:
    """
    Move stage relative to current position.

    Parameters
    ----------
    x : float
        X shift in millimeters.
    y : float
        Y shift in millimeters.
    """
    phenom.MoveBy(x * 0.001, y * 0.001)


def get_position() -> Tuple[float, float]:
    """
    Get the current absolute stage position.

    Optional driver function: during manual navigation, it is used to record where the user
    moved the stage. Without it, the sample center position is assumed.

    Returns
    -------
    tuple of float
        (x, y) stage coordinates in millimeters.
    """
    pos = phenom.GetStageModeAndPosition().position
    return pos.x * 1000, pos.y * 1000

# =============================================================================
# Imaging Parameters
# =============================================================================
def get_frame_width() -> float:
    """
    Get current SEM frame width.

    Returns
    -------
    float
        Frame width in millimeters.
    """
    return phenom.GetHFW() * 1000


def get_range_frame_width() -> Tuple[float, float]:
    """
    Get allowed SEM frame width range.

    Returns
    -------
    tuple(float, float)
        Minimum and maximum frame width in millimeters.
    """
    range_fw = phenom.GetHFWRange()
    return range_fw.begin * 1000, range_fw.end * 1000


def set_frame_width(frame_width: float) -> None:
    """
    Set SEM horizontal field width.

    Parameters
    ----------
    frame_width : float
        Desired frame width in millimeters.
    """
    phenom.SetHFW(frame_width * 0.001)

# =============================================================================
# Image Acquisition
# =============================================================================
def get_image_data(
    width: int = 1920,
    height: int = 1200,
    frame_avg: int = 1
) -> np.ndarray:
    """
    Acquire SEM image data.

    Parameters
    ----------
    width : int, default=1920
        Image width in pixels.
    height : int, default=1200
        Image height in pixels.
    frame_avg : int, default=1
        Number of frames to average for acquisition.

    Returns
    -------
    np.ndarray
        SEM image data as NumPy array (8-bit if thresholds are used).
    """
    # Note: For threshold-based processing in `autoemx`, image must be 8-bit.
    acq = phenom.SemAcquireImage(width, height, frame_avg)
    return np.asarray(acq.image)

# =============================================================================
# Brightness / Contrast Adjustments
# =============================================================================
def set_brightness(val: float) -> None:
    """Set SEM brightness level."""
    phenom.SetSemBrightness(val)


def set_contrast(val: float) -> None:
    """Set SEM contrast level."""
    phenom.SetSemContrast(val)

# =============================================================================
# Working Distance
# =============================================================================
def get_current_wd() -> float:
    """
    Get current SEM working distance.

    Returns
    -------
    float
        Working distance in millimeters.
    """
    wd_m = phenom.GetSemWD()  # API returns meters
    return wd_m * 1000

# =============================================================================
# Navigation Camera Functions
# =============================================================================
def to_nav() -> bool:
    """
    Switch instrument to navigation camera (optical) mode.

    Returns
    -------
    bool
        True if successfully switched to navigation mode, False otherwise.
    """
    # Wake up SEM if necessary.
    activate()
    try:
        phenom.MoveToNavCam()
    except Exception as e:
        logger.error("❌ Error: Failed to switch to navigation mode: %s", e)
        return False
    return True


def get_navigation_camera_image() -> Optional[np.ndarray]:
    """
    Acquire image from navigation camera.

    Returns
    -------
    np.ndarray or None
        RGB navigation camera image, or None if acquisition fails.
    """
    successful = to_nav()
    if not successful:
        return None

    # Adjust NavCam brightness and contrast for optimal C-tape detection.
    phenom.SetNavCamBrightness(0.34)
    phenom.SetNavCamContrast(0.27)

    # Define NavCam acquisition parameters.
    acqCamParams = ppi.CamParams()
    acqCamParams.size = ppi.Size(912, 912)
    acqCamParams.nFrames = 1

    # Acquire NavCam image.
    acqNavCam = phenom.NavCamAcquireImage(acqCamParams)

    # Temporary save to convert into NumPy format via OpenCV.
    temp_f = 'NavCam_temp.tiff'
    ppi.Save(acqNavCam.image, temp_f)
    navcam_im = cv2.imread(temp_f, cv2.IMREAD_COLOR)
    os.remove(temp_f)

    return navcam_im