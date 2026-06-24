#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
X-ray spectrum acquisition module for automated electron microscopy.

This module provides functionality for acquiring X-ray spectra (EDS/WDS) at
specific sample locations, with support for manual and automated spot selection.

Classes
-------
SpectrumAcquisition
    Handles X-ray spectrum acquisition and spot coordinate determination.

Created on 2026
@author: Andrea
"""
from typing import Optional, List, Tuple

from autoemx.utils.helper import EMError, Prompt_User

from autoemx._logging import get_logger
logger = get_logger(__name__)


class SpectrumAcquisition:
    """
    Handles X-ray spectrum acquisition and spot selection.
    
    Manages initialization of spectrum analyzers and acquisition of X-ray spectra
    at specified positions, supporting both manual and automated spot selection modes.
    
    Parameters
    ----------
    EM_controller : EM_Controller
        Reference to the parent EM_Controller.
    sample_cfg : SampleConfig
        Sample configuration.
    measurement_cfg : MeasurementConfig
        Measurement configuration.
    EM_driver_obj : module
        EM driver for hardware access.
    im_width : int
        Image width in pixels.
    im_height : int
        Image height in pixels.
    results_dir : str or None
        Directory for saving results.
    verbose : bool
        If True, print progress information.
    """
    
    def __init__(self, EM_controller, sample_cfg, measurement_cfg, EM_driver_obj,
                 im_width, im_height, results_dir=None, verbose=True):
        self.EM_controller = EM_controller
        self.sample_cfg = sample_cfg
        self.measurement_cfg = measurement_cfg
        self.EM_driver = EM_driver_obj
        self.im_width = im_width
        self.im_height = im_height
        self.results_dir = results_dir
        self.verbose = verbose
        
        self.analyzer = None
        self._xsp_callback_particle_active = False
        self._xsp_callback_particle_cntr: Optional[int] = None
    
    
    def _uses_xsp_spot_callback(self) -> bool:
        return self.EM_controller.powder_meas_cfg.par_spot_selection_mode == 'callback'
    
    
    def initialise_XS_analyzer(self, beam_voltage: float = None) -> None:
        """
        Initialize the EDS (Energy Dispersive X-ray Spectroscopy) analyzer.
        
        Optionally sets the electron beam voltage. If unspecified, uses the 
        voltage set during microscope initialization.
        
        Parameters
        ----------
        beam_voltage : float, optional
            Desired EM beam voltage in kilovolts (kV). If provided, sets the
            high tension to this value.
        """
        # Create EDS analyzer object
        self.analyzer = self.EM_driver.get_EDS_analyser_object()
        
        if beam_voltage is not None:
            # Set beam voltage (high tension) for EDS collection
            self.EM_driver.set_high_tension(beam_voltage)
    
    
    def get_XSp_coords(
        self, n_tot_sp_collected: int, frame_navigator, microscope_ctrl
    ) -> Tuple[bool, Optional[List[Tuple[float, float]]], Optional[int]]:
        """
        Determine X-ray spectrum acquisition coordinates for the next spectrum.
        
        Depending on the navigation mode, either prompts the user to select a spot
        manually, or automatically selects the next particle spot for 'powder' samples.
        
        Parameters
        ----------
        n_tot_sp_collected : int
            Current total number of spectra collected (used for labeling).
        frame_navigator : FrameNavigator
            Reference to frame navigator.
        microscope_ctrl : MicroscopeController
            Reference to microscope controller.
        
        Returns
        -------
        success : bool
            True if a spot was selected/found, False if user stopped or no more particles.
        spots_xy_list : list of tuple[float, float] or None
            List of (x, y) coordinates for acquisition spots, or None if unsuccessful.
        particle_cntr : int or None
            The particle counter/index, or None if not applicable.
        """
        
        self.EM_controller.ref_image = None

        if self.measurement_cfg.is_manual_navigation:
            return self._get_manual_coords(n_tot_sp_collected, frame_navigator, microscope_ctrl)
        
        elif self.sample_cfg.is_grid_acquisition:
            return self._get_grid_coords(frame_navigator, microscope_ctrl)
        
        elif self.sample_cfg.is_particle_acquisition:
            return self._get_particle_coords(n_tot_sp_collected, frame_navigator)
        
        logger.error(f"❌ Acquisition mode not implemented for sample type: {self.sample_cfg.type}")
        return False, None, None
    
    
    def _get_manual_coords(
        self, n_tot_sp_collected: int, frame_navigator, microscope_ctrl
    ) -> Tuple[bool, Optional[List[Tuple[float, float]]], Optional[int]]:
        """Get coordinates for manual spot selection."""
        prompt = Prompt_User(
            title="Select X-Ray Acquisition Spot",
            message=f"Center image on point where X-ray spectrum #{n_tot_sp_collected} is acquired."
        )
        prompt.run()
        
        if prompt.execution_stopped:
            logger.warning("⚠️ Execution stopped by the user.")
            return False, None, None
        
        if prompt.ok_pressed:
            spots_xy_list = [(int(self.im_width / 2), int(self.im_height / 2))]
            frame_navigator.current_frame_label = frame_navigator._frame_cntr
            frame_navigator._frame_cntr += 1
            return True, spots_xy_list, None
        
        return False, None, None
    
    
    def _get_grid_coords(
        self, frame_navigator, microscope_ctrl
    ) -> Tuple[bool, Optional[List[Tuple[float, float]]], Optional[int]]:
        """Get coordinates for grid-based acquisition."""
        # Center stage onto next acquisition spot
        movement_success = frame_navigator.go_to_next_frame(microscope_ctrl)
        if not movement_success:
            # Try to recalculate shifted grid
            recalc_success = frame_navigator._calc_bulk_grid_acquisition_spots()
            if recalc_success:
                movement_success = frame_navigator.go_to_next_frame(microscope_ctrl)
                if not movement_success:
                    logger.error("❌ Error moving to next frame")
                    return False, None, None
            else:
                return False, None, None
        
        spots_xy_list = [(int(self.im_width / 2), int(self.im_height / 2))]
        return True, spots_xy_list, None
    
    
    def _get_particle_coords(
        self, n_tot_sp_collected: int, frame_navigator
    ) -> Tuple[bool, Optional[List[Tuple[float, float]]], Optional[int]]:
        """Get coordinates for particle-based acquisition."""
        particle_finder = frame_navigator.particle_finder

        if self._uses_xsp_spot_callback():
            if not self._xsp_callback_particle_active:
                was_particle_found = particle_finder.go_to_next_particle()
                if not was_particle_found:
                    if self.verbose:
                        logger.warning('⚠️ No more particles could be found on the sample.')
                    return False, None, None
                self._xsp_callback_particle_active = True
                self._xsp_callback_particle_cntr = particle_finder.tot_par_cntr
            particle_cntr = self._xsp_callback_particle_cntr
        else:
            was_particle_found = particle_finder.go_to_next_particle()
            if not was_particle_found:
                if self.verbose:
                    logger.warning('⚠️ No more particles could be found on the sample.')
                return False, None, None
            particle_cntr = particle_finder.tot_par_cntr

        try:
            spots_xy_list = particle_finder.get_XS_acquisition_spots_coord_list(
                n_tot_sp_collected
            )
        except Exception as e:
            logger.error(f"❌ Error getting acquisition spot coordinates: {e}")
            if self._uses_xsp_spot_callback():
                self._xsp_callback_particle_active = False
                self._xsp_callback_particle_cntr = None
            return False, None, None

        if self._uses_xsp_spot_callback() and not spots_xy_list:
            self._xsp_callback_particle_active = False
            self._xsp_callback_particle_cntr = None

        self.EM_controller.ref_image = particle_finder.ref_image
        
        return True, spots_xy_list, particle_cntr
    
    
    def acquire_XS_spot_spectrum(
        self,
        x: float,
        y: float,
        max_acquisition_time: float,
        target_acquisition_counts: int,
        elements: Optional[List[str]] = None,
        msa_file_path: Optional[str] = None,
    ) -> Tuple:
        """
        Acquire an X-ray spectrum (EDS/WDS) at the specified position.
        
        Parameters
        ----------
        x, y : float
            X, Y pixel coordinates for spectrum acquisition.
        max_acquisition_time : float
            Maximum allowed acquisition time in seconds.
        target_acquisition_counts : int
            Target total X-ray counts for the spectrum.
        elements : list[str], optional
            Element symbols used by the microscope proprietary quantification to
            estimate instrument background.
        
        Returns
        -------
        spectrum_data : Any
            The acquired spectrum data.
        background_data : Any
            The measured background data.
        
        Raises
        ------
        EMError
            If the spectrum acquisition fails.
        """
        try:
            xs_coords = self.EM_driver.frame_pixel_to_rel_coords(
                (int(x), int(y)), self.im_width, self.im_height
            )
            x_stage, y_stage = xs_coords[0]
            spectrum_data, background_data, _, _ = (
                self.EM_driver.acquire_XS_spectral_data(
                    self.analyzer,
                    x_stage,
                    y_stage,
                    max_acquisition_time,
                    target_acquisition_counts,
                    elements=elements,
                    msa_file_export_path=msa_file_path
                )
            )
            return spectrum_data, background_data
        except Exception as e:
            raise EMError(f"Failed to acquire X-ray spectrum at pixel ({x}, {y}): {e}") from e
