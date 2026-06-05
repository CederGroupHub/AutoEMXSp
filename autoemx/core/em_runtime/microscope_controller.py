#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microscope hardware control module for automated electron microscopy operations.

This module provides low-level control of microscope hardware including beam settings,
focus adjustment, brightness/contrast, and stage movement.

Classes
-------
MicroscopeController
    Handles microscope hardware initialization and adjustment.

Created on 2026
@author: Andrea
"""
import time
import warnings

from autoemx.utils.helper import EMError, print_single_separator
from autoemx import microscope_drivers as EM_driver

from autoemx._logging import get_logger
logger = get_logger(__name__)
import numpy as np


class MicroscopeController:
    """
    Low-level microscope hardware controller.
    
    Handles initialization of SEM, beam settings, focus, brightness/contrast adjustment,
    and stage movement.
    
    Parameters
    ----------
    EM_driver_obj : module
        The EM_driver module for hardware communication.
    microscope_cfg : MicroscopeConfig
        Microscope configuration.
    measurement_cfg : MeasurementConfig
        Measurement/acquisition configuration.
    init_wd : float
        Initial working distance in mm.
    im_width : int
        Image width in pixels.
    im_height : int
        Image height in pixels.
    verbose : bool
        If True, print progress information.
    """
    
    def __init__(self, EM_driver_obj, microscope_cfg, measurement_cfg, 
                 init_wd, im_width, im_height, verbose=True):
        self.EM_driver = EM_driver_obj
        self.microscope_cfg = microscope_cfg
        self.measurement_cfg = measurement_cfg
        self.init_wd = init_wd
        self.im_width = im_width
        self.im_height = im_height
        self.verbose = verbose
        
        # Working distance bounds
        self._min_wd = init_wd - measurement_cfg.working_distance_tolerance
        self._max_wd = init_wd + measurement_cfg.working_distance_tolerance
        
        # Refresh time for BCF adjustments
        self.refresh_time = 1  # seconds
        
        # Track last adjustment time
        self._last_EM_adjustment_time = 0.0
    
    
    def initialise_SEM(self) -> None:
        """
        Activate and configure the Scanning Electron Microscope (SEM).
        
        This method performs the following steps:
            1. Wakes up the SEM if necessary.
            2. Switches to SEM mode.
            3. Sets the electron detector mode to backscattered electrons (BSD).
            4. Sets the beam voltage and current according to measurement configuration.
            5. Sets working distance.
            6. Adjusts focus, brightness, and contrast.
        
        Raises
        ------
        EMError
            If an error occurs during SEM activation or configuration.
        """
        if self.verbose:
            print_single_separator()
            logger.info("▶️ Activating SEM, and setting up...")
        
        try:
            # Wake up SEM if necessary
            self.EM_driver.activate()
            
            # Switch to SEM mode
            self.EM_driver.to_SEM()
            
            # Set detector type to BSD (Backscattered electron detector)
            self.EM_driver.set_electron_detector_mode(self.microscope_cfg.detector_type)
            
            # Set beam voltage (high tension) for EDS collection
            if self.measurement_cfg.beam_energy_keV:
                self.EM_driver.set_high_tension(self.measurement_cfg.beam_energy_keV)
            else:
                warnings.warn(
                    "No acceleration voltage provided via measurement_cfg.beam_energy_keV. "
                    "Using current microscope configurations",
                    UserWarning
                )
            
            # Set beam current for EDS collection
            if self.measurement_cfg.beam_current:
                self.EM_driver.set_beam_current(self.measurement_cfg.beam_current)
            else:
                warnings.warn(
                    "No beam current provided via measurement_cfg.beam_current. "
                    "Using current microscope configurations",
                    UserWarning
                )
            
            # Set working distance (needed for reliable autofocus)
            self.EM_driver.adjust_focus(self.init_wd)
            
            # Adjust focus, brightness, and contrast
            if self.verbose:
                print_single_separator()
                logger.info("🔬 Adjusting contrast, brightness, and focus.")
            
            self.adjust_BCF()
            
            if self.verbose:
                logger.info("✅ SEM initialisation completed.")
        
        except KeyError:
            raise
        except Exception as e:
            raise EMError(f"Error during SEM activation: {e}") from e
    
    
    def adjust_BCF(self) -> float:
        """
        Adjust brightness, contrast, and focus (BCF).
        
        If automatic brightness and contrast adjustment is enabled, calls the automatic 
        adjustment method. Otherwise, sets brightness and contrast to fixed values 
        then performs autofocus.
        
        Returns
        -------
        float
            Timestamp of adjustment completion.
            
        Raises
        ------
        EMError
            If an error occurs during adjustment.
        """
        try:
            if self.microscope_cfg.is_auto_BC:
                adj_time = self._autoadjust_BCF()
            else:
                self._set_frame_BC()
                adj_time = self._auto_focus()
                time.sleep(0.5)
            self._last_EM_adjustment_time = adj_time
            return adj_time
        except Exception as e:
            raise EMError(f"Failed to adjust brightness, contrast, and focus: {e}") from e
    
    
    def _set_frame_BC(self) -> None:
        """
        Set frame brightness and contrast to fixed values.
        
        Uses values from self.microscope_cfg.brightness and self.microscope_cfg.contrast.
        
        Raises
        ------
        EMError
            If the EM driver fails to set brightness or contrast.
        """
        try:
            self.EM_driver.set_brightness(self.microscope_cfg.brightness)
            self.EM_driver.set_contrast(self.microscope_cfg.contrast)
        except Exception as e:
            raise EMError(f"Failed to set brightness/contrast: {e}") from e
    
    
    def _autoadjust_BCF(self) -> float:
        """
        Automatically adjust brightness, contrast, and focus.
        
        Returns
        -------
        float
            Timestamp of adjustment completion.
            
        Raises
        ------
        EMError
            If an error occurs during adjustment.
        """
        try:
            self.EM_driver.auto_contrast_brightness()
            self._auto_focus()
            return time.time()
        except Exception as e:
            raise EMError(f"Failed to auto-adjust brightness/contrast/focus: {e}") from e
    
    
    def _auto_focus(self) -> float:
        """
        Automatically adjust focus within allowed bounds.
        
        Calls the EM driver's autofocus method and ensures the resulting working 
        distance is within allowed limits.
        
        Returns
        -------
        float
            Timestamp of adjustment completion.
            
        Raises
        ------
        EMError
            If an error occurs during autofocus.
        """
        try:
            wd = self.EM_driver.auto_focus()
            
            # If WD is out of allowed bounds, clip and readjust
            if not (self._min_wd < wd < self._max_wd):
                logger.warning(f"⚠️ Working distance of {wd:.1f} mm obtained through autofocus was out of accepted limits.")
                wd = float(np.clip(wd, self._min_wd, self._max_wd))
                logger.info(f"✅ WD was set to {wd:.1f} mm")
                self.EM_driver.adjust_focus(wd)
            
            return time.time()
        except Exception as e:
            raise EMError(f"Failed to auto-focus EM: {e}") from e
    
    
    def set_frame_width(self, frame_width: float) -> None:
        """
        Set the frame width at the microscope.
        
        Parameters
        ----------
        frame_width : float
            The desired frame width in millimeters.
            
        Raises
        ------
        EMError
            If an error occurs during frame width setting.
        """
        try:
            self.EM_driver.set_frame_width(frame_width)
        except Exception as e:
            raise EMError(f"Failed to set frame width: {e}") from e
    
    
    def move_to_pos(self, pos: tuple) -> None:
        """
        Move the EM stage to the specified (x, y) position.
        
        Parameters
        ----------
        pos : tuple of float
            Target (x, y) coordinates.
            
        Raises
        ------
        EMError
            If an error occurs during stage movement.
        """
        try:
            x, y = pos
            self.EM_driver.move_to(x, y)
        except Exception as e:
            raise EMError(f"Failed to move to desired position: {e}") from e
    
    
    @staticmethod
    def standby() -> None:
        """Put microscope in standby mode."""
        EM_driver.standby()
