#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated Electron Microscopy (EM) Particle Analysis Toolkit

This module provides the EM_Particle_Finder class for particle detection, selection of X-ray spectra (EDS, WDS) acquisition spots,
and particle size/statistics analysis.
It is designed to interface with an EM driver and an EM_Controller object to streamline the workflow from image acquisition to data export.

Currently only supports SEM.

Main Class
----------
EM_Particle_Finder
    Automates particle detection, mask generation, and X-ray spectra spot selection on detected particles.
    Supports both fully automated and manual (user-guided) collection modes.
    Provides methods for collecting particle statistics and managing frame navigation.

Example Usage
-------------
# EDS spot acquisition workflow:
>>> EM_controller = EM_Controller(...)
>>> EM_controller.initialise_SEM()
>>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, ...)
>>> while particle_finder.go_to_next_particle():
...     xy_spot_list = particle_finder.get_XS_acquisition_spots_coord_list()
...     for (x,y) in xy_spot_list:
...         EM_controller.acquire_XS_spectrum()


# Particle analysis workflow:
>>> EM_controller = EM_Controller(...)
>>> EM_controller.initialise_SEM()
>>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, results_dir="./results")
>>> particle_finder.get_particle_stats(n_par_target=500)

Notes
-----
- Requires an initialised EM_Controller object, unless used in development mode
- Requires a working EM_driver and appropriate hardware configuration, with the following API functions:
    Microscope Status & Image Acquisition
    -------------------------------------
    - is_at_EM
        Boolean attribute; True if running at the actual electron microscope.
    
    - get_image_data(width, height, channel)
        Acquire an image from the microscope with specified dimensions and channels.
    
    - get_frame_width()
        Get the current field of view/frame width (in mm).
    
    Stage & Navigation Control
    --------------------------
    - move_to(x, y)
        Move the microscope stage to the specified (x, y) position.

Created on Wed Jul 31 09:28:07 2024

@author: Andrea
"""
# Standard library imports
import os
import time
import random
import warnings

# Third-party imports
import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from autoemx.core.em_runtime.xsp_spot_selection import (
    XSpSpotSelectionContext,
    XSpSpotSelectorCallback,
    validate_xsp_spot_pixels,
)

# Local project imports
from autoemx.utils.helper import (
    Prompt_User,
    draw_scalebar,
    print_double_separator,
)
from autoemx.config import (
    PowderMeasurementConfig,
)
from autoemx.config.ledger_schemas import Coordinate2D, ParticleInfo
import autoemx.utils.constants as cnst
from autoemx import microscope_drivers as EM_driver
import autoemx.core.em_runtime.particle_segmentation_models as par_seg_models

from autoemx._logging import get_logger
logger = get_logger(__name__)


def split_diameter_range_um(
    d_min_um: float,
    d_max_um: float,
    max_decade: float = 10.0,
) -> List[Tuple[float, float]]:
    """Split a diameter range into bands of at most one order of magnitude.

    Bands are returned largest-first. Example: ``0.5..50`` → ``[(5, 50), (0.5, 5)]``.
    """
    if not np.isfinite(d_min_um) or not np.isfinite(d_max_um):
        raise ValueError("Diameter limits must be finite")
    if d_min_um <= 0 or d_max_um <= 0:
        raise ValueError("Diameter limits must be positive")
    if d_max_um < d_min_um:
        raise ValueError("d_max_um must be greater than or equal to d_min_um")
    if max_decade <= 1.0:
        raise ValueError("max_decade must be greater than 1")

    if np.isclose(d_min_um, d_max_um):
        return [(float(d_min_um), float(d_max_um))]

    bands: List[Tuple[float, float]] = []
    hi = float(d_max_um)
    while hi > d_min_um * (1.0 + 1e-12):
        lo = max(float(d_min_um), hi / max_decade)
        bands.append((lo, hi))
        if lo <= d_min_um + 1e-12:
            break
        hi = lo
    return bands


def diameter_um_to_area_um2(diameter_um: float) -> float:
    """Convert equivalent circular diameter (µm) to area (µm²)."""
    return float(np.pi * (diameter_um / 2.0) ** 2)


def frame_width_um_from_max_area(max_area_par: float) -> float:
    """Search-frame width (µm) from maximum accepted particle area (µm²)."""
    max_par_radius = np.sqrt(max_area_par / np.pi)
    return float(min(20.0 * max_par_radius, 500.0))

#%% Electron Microscope Particle Finder class    
class EM_Particle_Finder:
    """
    Automated particle analysis and X-ray spectra (EDS, WDS) acquisition in an electron microscope (EM).

    This class provides methods for particle detection, selection of X-ray spectra (EDS, WDS) acquisition spots,
    and particle size/statistics analysis.
    It is designed to interface with an EM driver and an EM_Controller object to streamline the workflow from
    image acquisition to data export.


    Main Methods
    ------------
    initialise_SEM()
        Wakes up the SEM microscope, sets measurement parameters, and evaluates locations to scan for particles.
    go_to_next_particle()
        Locates and moves to the next detected particle.
    get_XS_acquisition_spots_coord_list()
        Determines X-ray spectra (EDS, WDS) acquisition spots on the currently detected particle.
    get_particle_stats()
        Scans the full sample, collecting particle size distribution statistics.

    Examples
    --------
    # EDS spot acquisition workflow:
    >>> EM_controller = EM_Controller(...)
    >>> EM_controller.initialise_SEM()
    >>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, ...)
    >>> while particle_finder.go_to_next_particle():
    ...     xy_spot_list = particle_finder.get_XS_acquisition_spots_coord_list()
    ...     for (x,y) in xy_spot_list:
    ...         EM_controller.acquire_XS_spectrum()
    
    
    # Particle analysis workflow:
    >>> EM_controller = EM_Controller(...)
    >>> EM_controller.initialise_SEM()
    >>> particle_finder = EM_Particle_Finder(EM_controller, PowderMeasurementConfig, results_dir="./results")
    >>> particle_finder.get_particle_stats(n_par_target=500)

    Attributes
    ----------
    EM : EM_Controller
        Reference to the parent EM_Controller instance (must be initialised before use).
    powder_meas_cfg : PowderMeasurementConfig
        Configuration object for powder measurement (see dataclass for details).
    par_selection_mode : str
        Particle navigation mode from ``powder_meas_cfg`` ('auto' or 'manual').
    results_dir : str or None
        Directory for saving result images and data.
    verbose : bool
        If True, print progress and information to the console.
    development_mode : bool
        If True, enables extra visualizations and debug image saving.

    Internal Attributes
    -------------------
    _sample_id : str
        Identifier for the sample. Inherited from EM_Controller.
    _im_width : int
        Image width in pixels. Inherited from EM_Controller.
    _im_height : int
        Image height in pixels. Inherited from EM_Controller.
    tot_par_cntr : int
        Counter for the total number of particles analyzed.
    analyzed_pars : list of ParticleInfo
        Particle geometry records collected during particle size statistics.
        When a ledger is attached, this list is the same object as ``ledger.particles``.

    Notes
    -----
    - This class assumes that an EM_Controller object is already initialised and passed as `EM_controller`.
    - All configuration validation is performed in the respective dataclasses.
    - Attributes with a leading underscore are for internal use and should not be accessed directly by users.
    """
    def __init__(
        self,
        EM_controller,
        powder_meas_cfg: PowderMeasurementConfig,
        xsp_spot_selector: Optional[XSpSpotSelectorCallback] = None,
        results_dir: Optional[str] = None,
        verbose: bool = True,
        development_mode: bool = True
    ):
        """
        Initialize an EM_Particle_Finder object for automated/manual electron microscopy particle analysis.

        This class requires that an EM_Controller object has already been fully initialised and passed as `EM_controller`.

        Parameters
        ----------
        EM_controller : EM_Controller
            Reference to the parent EM_Controller instance (must be initialised before use).
        powder_meas_cfg : PowderMeasurementConfig
            Configuration object for powder measurement (see dataclass for details).
        xsp_spot_selector : callable, optional
            Required when ``powder_meas_cfg.par_spot_selection_mode='callback'``.
        results_dir : str, optional
            Directory to save result images and data (default: None).
        verbose : bool, optional
            If True, print progress and information to the console (default: True).
        development_mode : bool, optional
            If True, enables extra visualizations and debug image saving (default: True).
            
        Raises
        ------
        None

        Notes
        -----
        - All configuration validation is performed in the respective dataclasses.
        - This initializer assumes that configuration dataclasses are valid and complete.
        - Additional internal attributes are initialized for particle tracking and statistics.
        """
        # Electron microscope controller object
        self.EM = EM_controller
        self.powder_meas_cfg = powder_meas_cfg
        
        ### Inherit attributes
        if EM_controller is not None:
            # EM_controller is set to None when simply processing data
            self._sample_id = EM_controller.sample_id
            self._im_width = EM_controller.im_width
            self._im_height = EM_controller.im_height
        
        self.xsp_spot_selector = xsp_spot_selector
        if self.powder_meas_cfg.par_spot_selection_mode == 'callback' and self.xsp_spot_selector is None:
            raise ValueError(
                "powder_meas_cfg.par_spot_selection_mode='callback' requires xsp_spot_selector."
            )
        # NOTE: self.tot_par_cntr is initialized below, so _select_par_prompt_title will be set on first use
        if self.is_manual_particle_selection:
            self._select_par_prompt_title = f"Select position for particle #{self.tot_par_cntr}"
            self._select_par_prompt_message = (
                "Center the image around the particle you want\n"
                "to analyse, then click OK."
            )

        # --- General options
        self.results_dir = results_dir
        self.verbose = verbose
        self.development_mode = development_mode  # Set whether in code development mode. Triggers visualisations and image saving

        # --- Variable initializations
        self.tot_par_cntr = 0  # Keeps track of total number of particles analysed
        self.ref_image = None
        self.current_par_area_um2: Optional[float] = None
        self.current_par_coordinates_mm: Optional[tuple] = None
        self._fr_par_cntr = 0
        self._num_par_in_frame = 0
        self.analyzed_pars: List[ParticleInfo] = []
        self._ledger = None  # Optional SampleLedger for particle-stats / list selection
        self._ledger_path: Optional[str] = None
        self._list_particle_ids: Optional[List[int]] = None
        self._list_particle_coords_mm: Optional[List[Tuple[float, float]]] = None
        self._list_target_index: int = 0
        self._list_mode_absolute_ids: bool = False
        self._particles_by_id: Dict[int, ParticleInfo] = {}

    @property
    def is_manual_particle_selection(self) -> bool:
        """Whether particle navigation uses manual centering prompts."""
        return self.powder_meas_cfg.par_selection_mode == 'manual'

    @property
    def is_list_particle_selection(self) -> bool:
        """Whether particle navigation follows a caller-supplied id/coordinate list."""
        return self.powder_meas_cfg.par_selection_mode == 'list'

    @property
    def uses_absolute_particle_ids(self) -> bool:
        """True when ``tot_par_cntr`` is a ledger particle id (list-id mode)."""
        return bool(self._list_mode_absolute_ids)

    @property
    def par_selection_mode(self) -> str:
        return self.powder_meas_cfg.par_selection_mode

    @staticmethod
    def _normalize_stage_coordinates_mm(
        coordinates: Sequence[float],
    ) -> Tuple[float, float]:
        if coordinates is None or len(coordinates) != 2:
            raise ValueError("Stage coordinates must be a sequence of two finite floats (x, y) in mm")
        x, y = float(coordinates[0]), float(coordinates[1])
        if not (np.isfinite(x) and np.isfinite(y)):
            raise ValueError("Stage coordinates must be finite")
        return (x, y)

    def set_particle_targets(
        self,
        particle_ids: Optional[Iterable[int]] = None,
        particle_coords_mm: Optional[Iterable[Sequence[float]]] = None,
        ledger=None,
    ) -> None:
        """Configure list-mode targets (ledger particle ids XOR stage coordinates).

        Parameters
        ----------
        particle_ids :
            Existing ``ParticleInfo.id`` values. Requires ``ledger``; each id must
            exist and have absolute stage ``coordinates``.
        particle_coords_mm :
            Absolute stage positions ``(x, y)`` in mm for new particles.
        ledger :
            Sample ledger used to resolve and validate ``particle_ids``.
        """
        has_ids = particle_ids is not None
        has_coords = particle_coords_mm is not None
        if has_ids == has_coords:
            raise ValueError(
                "set_particle_targets requires exactly one of particle_ids or "
                "particle_coords_mm"
            )

        self._list_target_index = 0
        self._list_particle_ids = None
        self._list_particle_coords_mm = None
        self._list_mode_absolute_ids = False
        self._particles_by_id = {}

        if has_ids:
            if ledger is None:
                raise ValueError("particle_ids list mode requires a SampleLedger")
            ids = [int(pid) for pid in particle_ids]
            if not ids:
                raise ValueError("particle_ids must be a non-empty list")
            by_id = {int(p.id): p for p in ledger.particles}
            missing = [pid for pid in ids if pid not in by_id]
            if missing:
                raise ValueError(
                    f"Particle id(s) not found in ledger: {missing}"
                )
            missing_coords = [
                pid for pid in ids if by_id[pid].coordinates is None
            ]
            if missing_coords:
                raise ValueError(
                    "Particle id(s) in ledger are missing absolute stage "
                    f"coordinates: {missing_coords}"
                )
            self._list_particle_ids = ids
            self._list_mode_absolute_ids = True
            self._particles_by_id = by_id
            self._ledger = ledger
        else:
            coords = [
                self._normalize_stage_coordinates_mm(c) for c in particle_coords_mm
            ]
            if not coords:
                raise ValueError("particle_coords_mm must be a non-empty list")
            self._list_particle_coords_mm = coords
            self._list_mode_absolute_ids = False
            if ledger is not None:
                self._ledger = ledger

    def _default_particle_frame_width_mm(self) -> float:
        """Zoom FOV (mm) when particle size is unknown (coordinate list mode)."""
        max_radius_um = float(np.sqrt(self.powder_meas_cfg.max_area_par / np.pi))
        fw_um = max_radius_um * 2.0 * 1.8 + 10.0
        return float(fw_um / 1000.0)

    def _frame_width_mm_for_particle(self, particle: ParticleInfo) -> float:
        """Zoom FOV (mm) from an existing ParticleInfo size when available."""
        if particle.eq_diameter_um is not None and particle.eq_diameter_um > 0:
            eq_d_um = float(particle.eq_diameter_um)
        elif particle.area_um is not None and particle.area_um > 0:
            eq_d_um = float(2.0 * np.sqrt(float(particle.area_um) / np.pi))
        else:
            return self._default_particle_frame_width_mm()
        fw_um = eq_d_um * 1.8 + 10.0
        return float(fw_um / 1000.0)

    def _go_to_listed_particle(
        self,
        particle_id: Optional[int] = None,
        coordinates: Optional[Sequence[float]] = None,
    ) -> bool:
        """Move to a listed particle id or absolute stage coordinate."""
        self._check_EM_controller_initialization()

        if particle_id is not None and coordinates is not None:
            raise ValueError("Pass particle_id or coordinates, not both")

        if particle_id is None and coordinates is None:
            if self._list_particle_ids is not None:
                if self._list_target_index >= len(self._list_particle_ids):
                    return False
                particle_id = self._list_particle_ids[self._list_target_index]
                self._list_target_index += 1
            elif self._list_particle_coords_mm is not None:
                if self._list_target_index >= len(self._list_particle_coords_mm):
                    return False
                coordinates = self._list_particle_coords_mm[self._list_target_index]
                self._list_target_index += 1
            else:
                raise RuntimeError(
                    "list particle selection requires set_particle_targets(...) "
                    "before go_to_next_particle()"
                )

        if particle_id is not None:
            pid = int(particle_id)
            particle = self._particles_by_id.get(pid)
            if particle is None and self._ledger is not None:
                particle = next(
                    (p for p in self._ledger.particles if int(p.id) == pid),
                    None,
                )
            if particle is None:
                raise ValueError(f"Particle id {pid} not found in ledger")
            if particle.coordinates is None:
                raise ValueError(
                    f"Particle id {pid} has no absolute stage coordinates in the ledger"
                )
            coords = (
                float(particle.coordinates.x),
                float(particle.coordinates.y),
            )
            frame_width_mm = self._frame_width_mm_for_particle(particle)
            self.tot_par_cntr = pid
            self._list_mode_absolute_ids = True
            self.current_par_area_um2 = (
                float(particle.area_um) if particle.area_um is not None else None
            )
            self.current_par_coordinates_mm = coords
            # Keep ledger ParticleInfo in sync with the navigation target.
            particle.coordinates = Coordinate2D(x=coords[0], y=coords[1])
        else:
            coords = self._normalize_stage_coordinates_mm(coordinates)
            frame_width_mm = self._default_particle_frame_width_mm()
            self.tot_par_cntr += 1
            self._list_mode_absolute_ids = False
            self.current_par_coordinates_mm = coords

        self.EM.move_to_pos(coords)
        self.EM.set_frame_width(frame_width_mm)
        self.EM.pixel_size_um = frame_width_mm / self._im_width * 10**3
        self.EM.adjust_BCF()
        return True

    
    def _check_EM_controller_initialization(self) -> None:
        """
        Check whether the associated EM_Controller instance is initialized.
    
        Raises
        ------
        RuntimeError
            If the EM_Controller is not initialized.
    
        Notes
        -----
        This method should be called before performing any operation that
        requires an initialized microscope.
        """
        if not self.EM.is_initialized:
            raise RuntimeError(
                "EM_Controller is not initialized. Please call initialise_SEM() "
                "(or initialise_STEM(), if supported) before using this method."
            )
            
            
    #%% Particle Navigation
    # =============================================================================
    def go_to_next_particle(
        self,
        particle_id: Optional[int] = None,
        coordinates: Optional[Sequence[float]] = None,
    ):
        '''
        Moves the microscope to the next particle, centering and zooming on it.
        It also re-adjusts brightness, contrast and focus.

        Modes
        -----
        - ``auto``: find the next particle in the search-frame grid.
        - ``manual``: prompt the user to center a particle.
        - ``list``: consume the next target from ``set_particle_targets`` (ledger
          particle id or absolute stage coordinates in mm). Explicit ``particle_id``
          / ``coordinates`` arguments also select a specific target directly.

        Parameters
        ----------
        particle_id : int, optional
            Existing ledger particle id to revisit. Updates that ParticleInfo.
        coordinates : sequence of float, optional
            Absolute stage ``(x, y)`` in mm. Creates a new particle on acquisition.

        Returns
        -------
        bool
            True if it could successfully move to a particle.
            False if no more particles are present in the sample or if execution is stopped by the user.
        '''
        if particle_id is not None and coordinates is not None:
            raise ValueError("Pass particle_id or coordinates, not both")

        # Explicit target or configured list mode.
        if (
            self.is_list_particle_selection
            or particle_id is not None
            or coordinates is not None
        ):
            self.current_par_area_um2 = None
            self.current_par_coordinates_mm = None
            return self._go_to_listed_particle(
                particle_id=particle_id,
                coordinates=coordinates,
            )

        self.current_par_area_um2 = None
        self.current_par_coordinates_mm = None
        if not self.is_manual_particle_selection:
            self._check_EM_controller_initialization()
            
            # Automated search of particles
            # Check if all (or the max allowed amount) particles have been analysed in the current frames.
            # If yes, re-calculate particle positions in the next frame
            if self._fr_par_cntr >= self._num_par_in_frame or self._fr_par_cntr >= self.powder_meas_cfg.max_n_par_per_frame:
                # Loop until a frame with particles is found
                were_particles_found = False
                while not were_particles_found:
                    were_particles_found = self._get_particles_coordinates_in_frame()
                    # Check if all frames in sample have been analysed
                    if were_particles_found is None:
                        return False
    
            # Move center of image to the centroid of the particle
            par_pos = self.par_pos_abs_mm[self._fr_par_cntr]
            self.EM.move_to_pos(par_pos)
            self.current_par_coordinates_mm = (float(par_pos[0]), float(par_pos[1]))
    
            # Set frame width to zoom around particle
            frame_width_mm = self.par_fw_mm[self._fr_par_cntr]
            self.EM.set_frame_width(frame_width_mm)
    
            # Adjust EM settings (focus, contrast, brightness)
            self.EM.adjust_BCF()
    
            # Update counter of particles in current frame
            self._fr_par_cntr += 1
    
            # Update counter of total number of particles analysed
            self.tot_par_cntr += 1
    
        else:  # Manually look for particle
            prompt = Prompt_User(self._select_par_prompt_title, self._select_par_prompt_message)
            prompt.run()
    
            if prompt.execution_stopped:  # Check if execution was stopped after the loop
                logger.warning("⚠️ Execution stopped by the user.")
                return False
    
            if prompt.ok_pressed:
                frame_width_mm = EM_driver.get_frame_width()
                self.EM.pixel_size_um = frame_width_mm / self._im_width * 10**3  # um
                self.tot_par_cntr += 1
                current_pos = getattr(getattr(self.EM, "frame_navigator", None), "_current_pos", None)
                if current_pos is not None:
                    self.current_par_coordinates_mm = (float(current_pos[0]), float(current_pos[1]))
    
        return True
    
    
    #%% Particle Segmentation Operations
    # =============================================================================
    def _get_particles_coordinates_in_frame(self, frame_image=None, pixel_size=None, results_dir=None):
        '''
        Detects and extracts center coordinates and widths of suitable particles in the current frame.
    
        If running at the microscope, moves to the next frame and collects the image. If an image and pixel size 
        are provided, it uses these instead. The function applies a mask to segment particles, filters them 
        by area, and determines their positions and frame sizes. It also saves a visualization image with 
        circles drawn around detected particles for development purposes.
    
        Parameters
        ----------
        frame_image : ndarray, optional
            Grayscale image of the current frame. If not provided and running at the EM, the image is acquired.
        pixel_size : float, optional
            Pixel size in micrometers. Required if not running at the EM.
        results_dir : str, optional
            Directory to save results. Required for saving images in development mode.
    
        Returns
        -------
        bool or None
            True if at least one suitable particle was found in the frame.
            False if no particles were found.
            None if there are no more frames available.
        
        Notes
        -----
        - The function saves a visualization image with detected particles for development.
        - All commented-out imshow and development lines are preserved for debugging.
        - Particle positions and frame widths are stored as class attributes for later use.
        
        Potential Improvements
        ---------------------
        - Exclude particles that are close to larger particles in the direction of the EDS (Energy Dispersive Spectroscopy) detector.
          Large particles in the path can absorb or scatter X-rays emitted from smaller particles, degrading the quality and accuracy 
          of the spectral signal for those particles.
        
          Suggested implementation:
          - For each detected particle, determine if there are larger particles located "upstream" (i.e., between the particle and 
            the EDS detector direction).
          - To do this, add a control '_is_particle_shadowed()' together with '_is_particle_area_ok()' to select valid particles.
        '''
        if EM_driver.is_microscope_connected():
            self._check_EM_controller_initialization()
            move_to_frame_success = self.EM.go_to_next_frame()
            if move_to_frame_success:
                # Collect image
                frame_image = EM_driver.get_image_data(self._im_width, self._im_height, 1)
            else:
                # No more frames are available
                return None
        elif frame_image is not None and pixel_size is not None:
            self._im_height, self._im_width = frame_image.shape
            self.EM.pixel_size_um = pixel_size
            if results_dir:
                self.results_dir = results_dir
        else:
            raise ValueError('Function "_get_particles_coordinates_in_frame()" must be run at the microscope, or it needs to be passed both image and its pixel size')
        
        
        ### Select particles on substrate
        # Apply a Gaussian blur to remove single bright pixels
        frame_image = cv2.GaussianBlur(frame_image, (5, 5), 0)
        
        # Get mask of particles on substrate
        par_mask, _ = self._get_particles_on_substrate_mask(frame_image)
    
        # Find connected components
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)
        
        ### Filter out particles that are too small or too big, and store their centroid + size
        par_pos_pixels = []
        par_fw_pixels = []
        par_radius_pixels = []
        fw_margin_um = 5 # Pixel size in frame is a few um, so the particle will be shifted
        # fw_scale_factor = 1.8 # Multiplicative factor to obtain margins from particle
        for i in range(1, num_labels):  # Skip the background component (index 0)
            if self._is_particle_area_ok(stats[i, cv2.CC_STAT_AREA]):
                
                # Append particle centroid
                par_pos_pixels.append(centroids[i])
                
                # Calculate what would be the frame width in pixels in order to contain the particle fully in its width
                fw_width = stats[i, cv2.CC_STAT_WIDTH]
                # Same as above, but to contain the particle in its height
                fw_height = stats[i, cv2.CC_STAT_HEIGHT] / self._im_height * self._im_width
                
                # Select the largest frame width to make sure it contains the particle fully
                par_fw = max([fw_width, fw_height])
                # Append to list of frame_widths
                par_fw_pixels.append(par_fw)
                
                # Save particle radius, for proper circle size in saved image
                par_radius_pixels.append(max(fw_width, stats[i, cv2.CC_STAT_HEIGHT]) / 2)
                
        ### Visualize selected particles. Used for code development
        #     else:
        #         # Blacken pixels corresponding to excluded components
        #         par_mask[labels == i] = 0
        # cv2.imshow('Filtered mask', par_mask)
        
        
        # Save frame image annotating it with the identified particles
        filename = f"{self._sample_id}_fr{self.EM.current_frame_label}_particles"
        im_annotations = [{cnst.ANNOTATION_CIRCLE_KEY: (int(rad), center.astype(int), 2)} for rad, center in zip(par_radius_pixels, par_pos_pixels)]
        self.EM.save_frame_image(filename, im_annotations = im_annotations)
        
        
        # Return false if no particles were detected in the frame
        num_par = len(par_pos_pixels)
        if self.verbose:
            if num_par == 1:
                par_string = 'particle was'
            else:
                par_string = 'particles were'
            logger.info(f"✅ {num_par} {par_string} found in current frame")
        if num_par == 0:
            return False
        
        ### Convert dimensions from pixels to mm
        # Calculate absolute position of particle within the EM, in mm
        par_pos_abs_mm = self.EM.convert_pixel_pos_to_mm(par_pos_pixels)
        # Convert the frame width to mm
        par_fw_mm = (np.array(par_fw_pixels) * self.EM.pixel_size_um + 2 * fw_margin_um) * 10**-3
        
        # Store particle information
        self.par_pos_abs_mm = par_pos_abs_mm
        self.par_fw_mm = par_fw_mm
        self._num_par_in_frame = num_par
        
        # Initialise counter to keep track at how many particles within the frame have been analysed
        self._fr_par_cntr = 0
        
        # Returns True if at least 1 particle was found in the frame, otherwise returns False
        return True


    def _get_particles_on_substrate_mask(self, frame_image, save_image=False, persist_mask=True):
        '''
        Generates a binary mask of detected particles on the substrate from the input frame image.

        This function applies a brightness threshold to the input image to segment particles, 
        finds contours, and fills inner contours to ensure particles are fully masked. 
        Optionally, the mask image can be saved to disk. The function returns the mask and 
        the path where the mask image would be saved.

        Parameters
        ----------
        frame_image : ndarray
            The grayscale input image of the current frame.
        save_image : bool, optional
            If True, saves the binary mask image to disk (default: False).
        persist_mask : bool, optional
            If False, skip writing the mask image even in development mode. Used when a
            caller (e.g. per-particle masking) will save its own mask under a different name.

        Returns
        -------
        par_mask : ndarray
            Either a binary mask of detected particles, or a labels array, where the pixels
            of each particle are identified by a different integer (same as labels returned by
                                                                    cv2.ConnectedComponents)
        mask_img_path : str
            File path for where the mask image is (or would be) saved.
            
        Note
        ----
        The commented `cv2.imshow` line can be enabled for debugging visualization.

        '''
        if self.powder_meas_cfg.par_segmentation_model not in self.powder_meas_cfg.AVAILABLE_PAR_SEGMENTATION_MODELS:
            self.powder_meas_cfg.par_segmentation_model = "threshold_bright"
            warnings.warn(
                f"Chosen particle segmentation model {self.powder_meas_cfg.par_segmentation_model} not available."
                "Defaulting to 'threshold_bright'",
                UserWarning
            )
            
        if self.powder_meas_cfg.par_segmentation_model == "threshold_bright":
            # Apply the threshold to get a binary image
            _, par_mask = cv2.threshold(frame_image, self.powder_meas_cfg.par_brightness_thresh, 255, cv2.THRESH_BINARY)
            
            # Find all contours in the image and fill them
            contours, hierarchy = cv2.findContours(par_mask, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
            for i, contour in enumerate(contours):
                if hierarchy[0][i][3] != -1:  # If contour is inside another contour
                    cv2.drawContours(par_mask, [contour], 0, 255, -1)
                    
        elif self.powder_meas_cfg.par_segmentation_model in self.powder_meas_cfg.AVAILABLE_PAR_SEGMENTATION_MODELS:
            model_module = par_seg_models.PAR_SEGMENTATION_MODEL_REGISTRY[self.powder_meas_cfg.par_segmentation_model]
            par_mask = model_module.segment_particles(frame_image, self.powder_meas_cfg, save_image, self.EM)

        else:
            raise ValueError(f"Unknown error with current particle segmentation model {self.powder_meas_cfg.par_segmentation_model}")
        
        # cv2.imshow('Segmented Particles Mask', par_mask)
        
        mask_img_path = os.path.join(self.results_dir, self._sample_id + f'_fr{self.EM.current_frame_label}' + '_mask.png')
        # Mask is always saved when collecting particles. Avoids double saving
        save_image = save_image and not self.EM.measurement_cfg.type == self.EM.measurement_cfg.PARTICLE_STATS_MEAS_TYPE_KEY
        if persist_mask and (self.development_mode or save_image):
            draw_scalebar(par_mask, self.EM.pixel_size_um)
            cv2.imwrite(mask_img_path, par_mask)
        
        return par_mask, mask_img_path
    
    
    def _get_connected_components_with_stats(self, par_mask: np.ndarray):
        """
        Compute connected components with statistics.
    
        This function accepts either:
          1. A binary mask (0/255 or boolean), in which case it directly calls
             cv2.connectedComponentsWithStats.
          2. A pre-labeled image (integer labels, like output of cv2.connectedComponents),
             in which case stats and centroids are recomputed manually.
    
        Parameters
        ----------
        par_mask : np.ndarray
            Input binary mask (0/255 or bool) or label image (int32).
    
        Returns
        -------
        Same as cv2.connectedComponentsWithStats
        
        num_labels : int
            Number of connected components (including background).
        labels : np.ndarray
            Labeled image of the same size as input.
        stats : np.ndarray
            Statistics for each label. Shape: (num_labels, 5).
            Columns indexed by:
                cv2.CC_STAT_LEFT   (x)
                cv2.CC_STAT_TOP    (y)
                cv2.CC_STAT_WIDTH  (width)
                cv2.CC_STAT_HEIGHT (height)
                cv2.CC_STAT_AREA   (area)
        centroids : np.ndarray
            Centroids of each component. Shape: (num_labels, 2).
        """
    
        # --- Case 1: Binary image ---
        if par_mask.dtype == np.bool_ or np.array_equal(np.unique(par_mask), [0, 255]):
            num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(
                par_mask.astype(np.uint8), connectivity=8, ltype=cv2.CV_32S
            )
    
        else:
            # --- Case 2: Already a label image ---
            labels = par_mask.astype(np.int32, copy=False)
            num_labels = labels.max() + 1
    
            stats = np.zeros((num_labels, 5), dtype=np.int32)
            centroids = np.zeros((num_labels, 2), dtype=np.float64)
    
            for label in range(num_labels):
                mask = labels == label
                if not np.any(mask):
                    continue
    
                ys, xs = np.where(mask)
    
                x_min, x_max = xs.min(), xs.max()
                y_min, y_max = ys.min(), ys.max()
                w = x_max - x_min + 1
                h = y_max - y_min + 1
                area = mask.sum()
    
                stats[label, cv2.CC_STAT_LEFT]   = x_min
                stats[label, cv2.CC_STAT_TOP]    = y_min
                stats[label, cv2.CC_STAT_WIDTH]  = w
                stats[label, cv2.CC_STAT_HEIGHT] = h
                stats[label, cv2.CC_STAT_AREA]   = area
    
                centroids[label] = [xs.mean(), ys.mean()]
    
        return num_labels, labels, stats, centroids

    
    def is_particle_at_frame_edge(self, stats, i):
        '''
        Determines whether a particle is located at or near the edge of the image frame.
    
        A margin is used to account for detection tolerances. If any part of the particle's bounding box 
        is within `pixel_margin` pixels of the image border, it is considered to be at the edge.
    
        Parameters
        ----------
        stats : ndarray
            Connected components statistics array (as returned by OpenCV), where each row corresponds 
            to a detected particle and columns to bounding box info (LEFT, TOP, WIDTH, HEIGHT, etc.).
        i : int
            Index of the particle in the stats array.
    
        Returns
        -------
        bool
            True if the particle is at or near the edge of the image frame, False otherwise.
        '''
        # Apply margin to account for detection tolerances
        pixel_margin = 3
    
        # Check if the particle's bounding box touches or is near any image edge
        is_particle_at_frame_edge = any([
            stats[i, cv2.CC_STAT_LEFT] <= pixel_margin,
            stats[i, cv2.CC_STAT_TOP] <= pixel_margin,
            stats[i, cv2.CC_STAT_LEFT] + stats[i, cv2.CC_STAT_WIDTH] >= self._im_width - pixel_margin,
            stats[i, cv2.CC_STAT_TOP] + stats[i, cv2.CC_STAT_HEIGHT] >= self._im_height - pixel_margin
        ])
    
        return is_particle_at_frame_edge


    def _is_particle_area_ok(self, area):
        '''
        Checks whether the given particle area is within acceptable limits.
    
        The function converts the minimum and maximum allowed particle areas from µm² to pixel units 
        using the current pixel size. It then determines if the provided area is within this range.
        If in manual collection mode and the area is out of bounds, it prints a message to the user.
    
        Parameters
        ----------
        area : float
            The area of the particle in pixels.
    
        Returns
        -------
        bool
            True if the particle area is within the allowed range, False otherwise.
        '''
        # Calculate acceptable particle area dimensions with current pixel size
        min_area_pixels = self.powder_meas_cfg.min_area_par / self.EM.pixel_size_um**2  # e.g., 4um^2 for a 2um x 2um particle
        max_area_pixels = self.powder_meas_cfg.max_area_par / self.EM.pixel_size_um**2  # e.g., 100um^2 for a 10um x 10um particle
    
        is_par_large_enough = area >= min_area_pixels
        if not is_par_large_enough and self.is_manual_particle_selection:
            logger.warning(f"⚠️ Selected particle is too small ({area*self.EM.pixel_size_um**2:.2f} um^2), select another one")
    
        is_par_small_enough = area <= max_area_pixels
        if not is_par_small_enough and self.is_manual_particle_selection:
            logger.warning(f"⚠️ Selected particle is too large ({area*self.EM.pixel_size_um**2:.2f} um^2), select another one")
    
        is_par_area_ok = is_par_large_enough and is_par_small_enough
    
        return is_par_area_ok


    def _capture_ref_frame(self, par_image=None, pixel_size_um=None):
        '''
        Capture (or accept) the grayscale frame used for spot selection.

        At the microscope, a fresh frame is acquired via the driver. Offline, the
        provided ``par_image`` is used and image dimensions / pixel size are
        updated from it. This is the single place that defines how the reference
        frame is obtained, shared by automated and callback spot selection.

        Parameters
        ----------
        par_image : ndarray, optional
            Grayscale frame to use when not running at the microscope.
        pixel_size_um : float, optional
            Pixel size in micrometers. Required together with ``par_image`` offline.

        Returns
        -------
        ndarray
            The grayscale frame image.
        '''
        if EM_driver.is_microscope_connected():
            self._check_EM_controller_initialization()
            return EM_driver.get_image_data(self._im_width, self._im_height, 1)
        if par_image is not None and pixel_size_um is not None:
            self._im_height, self._im_width = par_image.shape
            self.EM.pixel_size_um = pixel_size_um
            return par_image
        raise ValueError('This function must be run at the microscope, or it needs to be passed both image and its pixel size')


    def _segment_center_particle(
        self,
        par_image: np.ndarray,
        *,
        centering: bool = False,
    ) -> Optional[Tuple[np.ndarray, int, np.ndarray, np.ndarray]]:
        """
        Segment the particle at the image center with the configured segmentation model.

        Returns
        -------
        tuple or None
            ``(binary_mask, par_label, stats, centroids)`` for the selected particle.
            Area is computed from this full (non-eroded) mask before any spot-selection
            erosion. Returns None if no suitable particle is found.
        """
        par_mask, _ = self._get_particles_on_substrate_mask(par_image, persist_mask=False)
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)

        if num_labels == 1:
            return None

        center_x = int(self._im_width / 2)
        center_y = int(self._im_height / 2)
        par_label = int(labels[center_y, center_x])

        center_on_particle = par_label > 0
        center_area_ok = center_on_particle and self._is_particle_area_ok(
            stats[par_label, cv2.CC_STAT_AREA]
        )

        # Prefer the particle under the image center. Offline, keep it even when it
        # is outside configured area limits so size/spot selection still target the
        # particle being analysed. At the microscope, oversized/undersized center
        # particles trigger a re-center onto the nearest area-ok neighbour.
        if center_area_ok or (center_on_particle and not EM_driver.is_microscope_connected()):
            binary_mask = np.where(labels == par_label, 255, 0).astype(np.uint8)
            return binary_mask, par_label, stats, centroids

        if centering:
            return None

        distances = np.linalg.norm(centroids[1:] - np.array([center_x, center_y]), axis=1)
        sorted_indices = np.argsort(distances) + 1  # Skip background (index 0)
        for label in sorted_indices:
            label = int(label)
            if not EM_driver.is_microscope_connected():
                binary_mask = np.where(labels == label, 255, 0).astype(np.uint8)
                return binary_mask, label, stats, centroids
            if self._is_particle_area_ok(stats[label, cv2.CC_STAT_AREA]):
                new_center = self.EM.convert_pixel_pos_to_mm(centroids[label])
                self.EM.move_to_pos(new_center)
                recentered_image = self._capture_ref_frame()
                return self._segment_center_particle(recentered_image, centering=True)
        return None

    def _store_current_particle_area(self, area_pixels: float) -> float:
        """Store particle area in µm² from a non-eroded pixel area."""
        area_um = float(area_pixels) * float(self.EM.pixel_size_um) ** 2
        self.current_par_area_um2 = area_um
        return area_um

    def _get_particle_mask(self, par_image=None, pixel_size_um=None, results_dir=None, centering=False):
        '''
        Returns a binary mask of the particle at the center of the image, if the particle is large enough.
        If no particle is present at the center, the function searches for and centers on the next closest 
        sufficiently large particle. This approach makes the function robust against drift or movement of 
        the particle during analysis. If no valid particle is found, the function returns None.
    
        Parameters
        ----------
        par_image : ndarray, optional
            Grayscale image of the current frame. If not provided and running at the EM, the image is acquired.
        pixel_size_um : float, optional
            Pixel size in micrometers. Required if not running at the EM.
        results_dir : str, optional
            Directory to save results. Required for saving the mask in development mode.
        centering : bool, optional
            If True, indicates this is a second attempt at centering on a particle (prevents infinite recursion).
    
        Returns
        -------
        tuple or None
            (par_mask, par_image): Binary mask of the selected particle and the corresponding image.
            None: If no suitable particle is found.
        
        Notes
        -----
        - Uses ``powder_meas_cfg.par_segmentation_model`` (default ``threshold_bright``).
        - Stores ``current_par_area_um2`` from the full mask before any edge erosion.
        - The function is robust to small misalignments: if the center particle is not valid, it finds the next closest.
        - In development mode, it saves the resulting mask with a scalebar overlay.
        - The commented `cv2.imshow` lines can be enabled for debugging visualization.
        '''
        if not EM_driver.is_microscope_connected() and results_dir:
            self.results_dir = results_dir
        par_image = self._capture_ref_frame(par_image, pixel_size_um)

        segmented = self._segment_center_particle(par_image, centering=centering)
        if segmented is None:
            if not centering:
                self.current_par_area_um2 = None
            return None

        par_mask, par_label, stats, _centroids = segmented
        # Area from the full (non-eroded) mask — spot selection may erode later.
        self._store_current_particle_area(stats[par_label, cv2.CC_STAT_AREA])

        if self.development_mode and self.results_dir:
            mask_to_save = draw_scalebar(par_mask.copy(), self.EM.pixel_size_um)
            cv2.imwrite(os.path.join(
                self.results_dir,
                self._sample_id + f'_par{self.tot_par_cntr}_fr{self.EM.current_frame_label}_mask.png'
            ), mask_to_save)

        return (par_mask, par_image)


    def _erode_particle_mask(self, par_mask, margin, erode_inner=False):
        """
        Erode a binary particle mask by a specified margin.
    
        This function erodes the mask either only at the outer boundary or at both the outer boundary and inner holes,
        depending on the `erode_inner` flag.
    
        Parameters
        ----------
        par_mask : np.ndarray
            Binary mask of the particle (dtype uint8, values 0 and 255).
        margin : int
            The erosion margin (in pixels). The structuring element will be of size (2 * margin + 1).
        erode_inner : bool, optional
            If False (default), only the outer boundary of the mask is eroded,
            leaving holes inside the particle unaffected.
            If True, both the outer boundary and any internal holes are eroded.
    
        Returns
        -------
        final_mask : np.ndarray
            The eroded binary mask (same dtype and shape as `par_mask`).
    
        Notes
        -----
        - This method assumes the mask is a binary 8-bit image (0 for background, 255 for foreground).
        - The erosion uses an elliptical structuring element.
    
        Examples
        --------
        >>> eroded = self._erode_particle_mask(par_mask, margin=3, erode_inner=False)
        """
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * margin + 1, 2 * margin + 1))
        if erode_inner:
            # Erode everywhere (outer and inner edges)
            final_mask = cv2.erode(par_mask, kernel)
        else:
            # Erode only the outer contour, not holes
            contours, hierarchy = cv2.findContours(par_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            outer_contour_mask = np.zeros_like(par_mask)
            cv2.drawContours(outer_contour_mask, contours, -1, 255, -1)  # Fill outer boundaries
            eroded_outer_contour_mask = cv2.erode(outer_contour_mask, kernel)
            final_mask = cv2.bitwise_and(par_mask, eroded_outer_contour_mask)
        return final_mask
    
    
    def _find_particle_bright_regions(self, final_mask, par_image):
        """
        Detect bright regions within a particle mask in the input image, in order to acquire X-ray spectra
        from thickest regions of particle and maximize X-ray intensity.
    
        This function applies a mask to the input image, blurs it, normalizes the intensity,
        thresholds to find bright spots, and returns the thresholded image and the minimum area (in pixels)
        for spot detection.
    
        Parameters
        ----------
        final_mask : np.ndarray
            Binary mask (uint8, values 0 and 255) specifying the region of interest (the particle).
        par_image : np.ndarray
            Grayscale image (uint8) of the particle region.
    
        Returns
        -------
        thresholded_image : np.ndarray
            Binary image (uint8, values 0 and 255) showing detected bright spots within the particle mask.
        min_area_pixels : int
            Minimum spot area in pixels, calculated as corresponding to 0.1 μm².
    
        Notes
        -----
        - The function applies a Gaussian blur before normalization and thresholding.
        - The threshold value is taken from `self.powder_meas_cfg.par_xy_spots_thresh`.
        - The minimum area is computed from `self.EM.pixel_size_um`.
        - If `self.development_mode` is True and `self.results_dir` is set, a debug image is saved with a scalebar.
    
        Examples
        --------
        >>> spots_mask, min_area = self._find_particle_bright_regions(final_mask, par_image)
        """
        masked_image = cv2.bitwise_and(par_image, par_image, mask=final_mask)
        # cv2.imshow('CV Image', masked_image)
        if np.all(masked_image == 0):
            return np.zeros_like(masked_image), 0
        masked_im_blurred = cv2.GaussianBlur(masked_image, (5, 5), 0)
        # cv2.imshow('Blurred Image', masked_im_blurred)
        
        # Normalise particle intensity to select thickest regions, regardless of intensity of neighbouring features
        norm_masked_image = (masked_im_blurred / np.max(masked_im_blurred) * 255).astype(np.uint8)
        _, thresholded_image = cv2.threshold(norm_masked_image, self.powder_meas_cfg.par_xy_spots_thresh, 255, cv2.THRESH_BINARY)
        # cv2.imshow('CV Image', thresholded_image)
        # cv2.imshow('CV Image 2', cv2.bitwise_and(par_image, thresholded_image))
        min_area_pixels = int(0.1 / self.EM.pixel_size_um ** 2)
        if self.development_mode and self.results_dir:
            eroded_par_mask = draw_scalebar(thresholded_image, self.EM.pixel_size_um)
            cv2.imwrite(os.path.join(self.results_dir, self._sample_id + f'_par{self.tot_par_cntr}_fr{self.EM.current_frame_label}_maskXY.png'), eroded_par_mask)
        return thresholded_image, min_area_pixels


    def prepare_mask_for_visualization(self, mask: np.ndarray) -> np.ndarray:
        """
        Prepare a segmentation mask for visualization.
    
        Behavior:
        - Binary masks with values in {0, 1} are scaled to [0, 255].
        - Binary masks with values in {0, 255} are returned unchanged.
        - Integer label masks have brightness reversed so higher labels are brighter,
          with a minimum brightness of 30 for positive values. Background (0) remains black.
        - Floating-point masks are normalized to [0, 255] with the same rules for
          positive values and background.
        """
        mask = np.asarray(mask)
        unique_vals = np.unique(mask)
    
        # --- Binary mask [0, 1] ---
        if np.array_equal(unique_vals, [0, 1]):
            return (mask * 255).astype(np.uint8)
    
        # --- Binary mask [0, 255] ---
        elif np.array_equal(unique_vals, [0, 255]):
            return mask.astype(np.uint8)
    
        # --- Integer label masks ---
        elif np.issubdtype(mask.dtype, np.integer):
            max_val = mask.max()
            if max_val > 0:
                scaled = mask.astype(np.float32)
                pos_mask = scaled > 0
    
                # Reverse intensity so higher labels → brighter
                scaled[pos_mask] = (max_val - scaled[pos_mask]) * (255.0 / max_val)
    
                # Minimum brightness for positive labels
                scaled[pos_mask] = np.clip(scaled[pos_mask], 30, 255)
                scaled[~pos_mask] = 0
                return scaled.astype(np.uint8)
            else:
                return mask.astype(np.uint8)
    
        # --- Floating-point masks ---
        else:
            min_val = mask.min()
            max_val = mask.max()
            if max_val > min_val:
                scaled = (mask - min_val) / (max_val - min_val) * 255.0
                scaled = scaled.astype(np.float32)
    
                pos_mask = mask > 0
                scaled[pos_mask] = np.clip(scaled[pos_mask], 30, 255)
                scaled[~pos_mask] = 0
    
                return scaled.astype(np.uint8)
            else:
                return np.zeros_like(mask, dtype=np.uint8)
    
    
    #%% Selection of spots for X-Ray spectra acquisition
    # =============================================================================
    def get_XS_acquisition_spots_coord_list(
            self, n_tot_sp_collected, 
            par_image=None, pixel_size_um=None, results_dir=None):
        '''
        Returns a list of pixel coordinates (in the current image) for X-ray spectrum spot collection on a particle.
    
        The function finds a suitable particle mask, erodes the mask to avoid particle edges, finds bright regions (or peak spots),
        and selects up to `powder_meas_cfg.max_spectra_per_par` spots per particle with a minimum distance between them (determined through 'powder_meas_cfg.par_mask_margin').
        The selection strategy for features and spacing can be controlled via powder_meas_cfg.par_feature_selection and powder_meas_cfg.par_spot_spacing.
    
        Parameters
        ----------
        n_tot_sp_collected : int
            Counter for the total number of spectra collected (used for labeling).
        par_image : ndarray, optional
            Grayscale image of the current frame. If not provided and running at the SEM, the image is acquired.
        pixel_size_um : float, optional
            Pixel size in micrometers. Required if not running at the SEM.
        results_dir : str, optional
            Directory to save result images.
    
        Returns
        -------
        selected_points : list[tuple[int, int]]
            List of selected (x, y) spot coordinates in image pixel units.
    
        Potential Improvements
        ---------------------
        - Select points to ensure all phases are individuated, biasing spot selection to ensure representation
            of phases with different bright/dark contrast. Currently, only the brightest spots are picked,
            which may miss some phases.
        - Add detection of both peaks and valleys to ensure both bright and dark spots are tested.
        - Exclude points that are close to larger particles along the X-ray emission path to the EDS detector,
          as these may absorb emitted X-rays and degrade spectral signal.
        '''
        # In callback mode, spots are supplied by xsp_spot_selector. Skip the
        # particle mask computation entirely and reuse a single reference frame
        # per particle (see _get_callback_xsp_spots).
        if self.powder_meas_cfg.par_spot_selection_mode == 'callback':
            return self._get_callback_xsp_spots(
                n_tot_sp_collected,
                par_image=par_image,
                pixel_size_um=pixel_size_um,
            )

        # --- 1. Acquire or prepare the particle mask and image ---
        if EM_driver.is_microscope_connected():
            par_mask_return = self._get_particle_mask()
        elif par_image is not None and pixel_size_um is not None:
            self._im_height, self._im_width = par_image.shape
            self.EM.pixel_size_um = pixel_size_um
            par_mask_return = self._get_particle_mask(par_image, pixel_size_um)
            if results_dir:
                self.results_dir = results_dir
        else:
            raise ValueError('This function must be run at the microscope, or it needs to be passed both image and its pixel size')
    
        # Check if a particle was detected. If not, return empty list
        if par_mask_return is None:
            self.ref_image = None
            return []
        else:
            par_mask, par_image = par_mask_return
    
        # --- 2. Erode the particle mask to avoid edge effects ---
        margin = max(10, int(self.powder_meas_cfg.par_mask_margin / self.EM.pixel_size_um))
        final_mask = self._erode_particle_mask(par_mask, margin)
    
        # --- 3. Find bright points in image, which indicate highest regions on particle---
        thresholded_image, min_area_pixels = self._find_particle_bright_regions(final_mask, par_image)
    
        # --- 4. Collect candidate points from thresholded image ---
        all_points = self._collect_candidate_points(
            thresholded_image, par_image, self.powder_meas_cfg.par_feature_selection, min_area_pixels
        )
    
        if len(all_points) == 0:
            self.ref_image = None
            return []
    
        # --- 5. Select points with minimum distance and maximum count ---
        min_distance_xsp_spots = (self.powder_meas_cfg.xsp_spots_distance_um / self.EM.pixel_size_um)
        if self.powder_meas_cfg.par_spot_spacing == 'random':
            selected_points = EM_Particle_Finder._select_XSspots_randomly(
                all_points, max_points=self.powder_meas_cfg.max_spectra_per_par, min_distance=min_distance_xsp_spots
            )
        elif self.powder_meas_cfg.par_spot_spacing == 'maximized':
            selected_points = EM_Particle_Finder._select_evenly_spaced_XSspots(
                all_points, max_points=self.powder_meas_cfg.max_spectra_per_par, min_distance=min_distance_xsp_spots
            )
    
        # --- 6. Annotate image and save ---
        if self.development_mode and self.results_dir:
            color_image = cv2.cvtColor(par_image, cv2.COLOR_GRAY2BGR)
            for center in selected_points:
                # Add circle indicating where X-ray spectrum was collected
                cv2.circle(color_image, center, 10, (0, 0, 255), -1)
                label_pos = (center[0] - 30, center[1] - 15)
                cv2.putText(color_image, str(n_tot_sp_collected), label_pos, cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2, cv2.LINE_AA)
                n_tot_sp_collected += 1
            color_image = draw_scalebar(color_image, self.EM.pixel_size_um)
            # cv2.imshow('Selected XS spots', color_image)
            cv2.imwrite(os.path.join(self.results_dir, self._sample_id + f'_par{self.tot_par_cntr}_fr{self.EM.current_frame_label}_xyspots.png'), color_image)

        self.ref_image = par_image
        return [tuple(map(int, pt)) for pt in selected_points]

    def _get_callback_xsp_spots(
        self,
        n_tot_sp_collected: int,
        par_image=None,
        pixel_size_um=None,
    ) -> List[tuple]:
        """Invoke xsp_spot_selector and validate returned pixel coordinates.

        The reference frame is captured once per particle (when ``self.ref_image``
        is None) and reused for every spot batch on that particle, so the
        drift-correction reference stays fixed across that particle's spectra. No
        particle mask is computed in this mode.
        """
        if self.ref_image is None:
            self.ref_image = self._capture_ref_frame(par_image, pixel_size_um)
            # Measure particle size independently of the callback. Uses the configured
            # segmentation model on the raw (non-eroded) reference frame.
            segmented = self._segment_center_particle(self.ref_image)
            if segmented is not None:
                _mask, par_label, stats, _centroids = segmented
                self._store_current_particle_area(stats[par_label, cv2.CC_STAT_AREA])
            else:
                self.current_par_area_um2 = None

        frame_label = getattr(self.EM, 'current_frame_label', '')
        if frame_label is None:
            frame_label = ''

        context = XSpSpotSelectionContext(
            particle_id=self.tot_par_cntr,
            n_tot_sp_collected=n_tot_sp_collected,
            ref_image=self.ref_image,
            pixel_size_um=float(self.EM.pixel_size_um),
            frame_label=str(frame_label),
            im_width=int(self._im_width),
            im_height=int(self._im_height),
            particle_finder=self,
        )
        try:
            raw_spots = self.xsp_spot_selector(context)
        except Exception as exc:
            logger.error(f"❌ xsp_spot_selector callback failed: {exc}")
            raise

        if not raw_spots:
            return []

        selected_points = validate_xsp_spot_pixels(
            raw_spots,
            im_width=int(self._im_width),
            im_height=int(self._im_height),
        )
        return selected_points

        
    def _collect_candidate_points(self, thresholded_image, par_image, feature_selection, min_area_pixels):
        """
        Collect candidate EDS/WDS spot coordinates from thresholded regions.
    
        This function identifies connected components (regions) in the thresholded image,
        filters them by minimum area, and collects candidate points for X-ray spot acquisition.
        The method of selection depends on the `feature_selection` argument:
        - 'random': all pixel coordinates within the component are collected.
        - 'peaks': only the brightest pixel within each component is selected.
    
        Parameters
        ----------
        thresholded_image : np.ndarray
            Binary image (uint8, values 0 and 255) indicating candidate regions for spot selection.
        par_image : np.ndarray
            Grayscale image (uint8) of the particle, used for peak detection.
        feature_selection : str
            'random' to select all pixels in the region, 'peaks' to select local maxima (bright spots).
        min_area_pixels : int
            Minimum area (in pixels) for a component to be considered.
    
        Returns
        -------
        all_points : list of tuple
            List of (x, y) coordinates of candidate points (in pixel units).
    
        Notes
        -----
        - For 'random', all pixels in each sufficiently large component are returned.
        - For 'peaks', local maxima are found by repeatedly masking out regions around each peak.
        - The margin for masking out peaks is determined by `self.powder_meas_cfg.par_mask_margin` and `self.EM.pixel_size_um`.
        - The maximum number of iterations for peak finding is `100 * self.powder_meas_cfg.max_spectra_per_par`.
        - Only components with area >= `min_area_pixels` are considered.
        """
        all_points = []
        num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(thresholded_image, connectivity=8)
        margin = max(10, int(self.powder_meas_cfg.par_mask_margin / self.EM.pixel_size_um))
        for i in range(1, num_labels):  # Skip background
            if stats[i, cv2.CC_STAT_AREA] >= min_area_pixels:
                if feature_selection == 'random':
                    components_coords = np.where(labels == i)
                    all_points += list(zip(components_coords[1], components_coords[0]))
                elif feature_selection == 'peaks':
                    component_mask = (labels == i).astype(np.uint8) * 255
                    component_image = cv2.bitwise_and(par_image, par_image, mask=component_mask)
                    component_image = cv2.GaussianBlur(component_image, (5, 5), 0)
                    max_iter = 100 * self.powder_meas_cfg.max_spectra_per_par
                    iter_cntr = 0
                    while iter_cntr < max_iter:
                        iter_cntr += 1
                        min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(component_image)
                        if max_val < self.powder_meas_cfg.par_xy_spots_thresh:
                            break
                        _, thresh = cv2.threshold(component_image, int(max_val * 0.97), 255, cv2.THRESH_BINARY)
                        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                        for contour in contours:
                            x, y, w, h = cv2.boundingRect(contour)
                            center_x = x + w // 2
                            center_y = y + h // 2
                            center_pnt = (center_x, center_y)
                            all_points.append(center_pnt)
                            cv2.circle(component_image, center_pnt, margin, (0, 0, 0), -1)
        return all_points
    
    
    @staticmethod
    def _select_XSspots_randomly(points, max_points, min_distance):
        '''
        Randomly selects up to `max_points` from a list of acquisition spots, ensuring that each selected point
        is at least `min_distance` away from all others already selected.
    
        Parameters
        ----------
        points : list of tuple or ndarray
            List of (x, y) coordinates to choose from.
        max_points : int
            Maximum number of points to select.
        min_distance : float
            Minimum allowed distance between any two selected points.
    
        Returns
        -------
        selected_points : list
            List of selected (x, y) points.
        '''
        # Shuffle points to unbias the choice of position
        random.shuffle(points)
    
        # Check if some selected points must be cut out
        if len(points) <= max_points:
            return points
    
        selected_points = [points[0]]
        for point in points:
            # Check if enough points have been selected
            if len(selected_points) >= max_points:
                break
    
            # Check if point is distant enough from all selected points
            if point not in selected_points and all(np.linalg.norm(np.array(point) - np.array(p)) > min_distance for p in selected_points):
                selected_points.append(point)
    
        return selected_points
    
    
    @staticmethod
    def _select_evenly_spaced_XSspots(points, max_points, min_distance):
        '''
        Selects up to `max_points` from a list of acquisition spots, maximizing the minimum distance between any two selected points.
        Points are chosen iteratively: each new point is the one farthest (in minimum distance) from those already selected.
    
        Parameters
        ----------
        points : list of tuple or ndarray
            List of (x, y) coordinates to choose from.
        max_points : int
            Maximum number of points to select.
        min_distance : float
            Minimum allowed distance between any two selected points (not strictly enforced, but likely achieved).
    
        Returns
        -------
        selected_points : list
            List of selected (x, y) points.
        '''
        # Shuffle points to unbias the choice of position
        random.shuffle(points)
    
        # Check if some selected points must be cut out
        if len(points) <= max_points:
            return points
    
        selected_points = [points[0]]
        while len(selected_points) < max_points:
            max_min_distance = -1                   
            best_point = None
            for point in points:
                if point not in selected_points:
                    min_dist_to_selected = min(np.linalg.norm(np.array(point) - np.array(p)) for p in selected_points)
                    if min_dist_to_selected > max_min_distance:
                        max_min_distance = min_dist_to_selected
                        best_point = point
            if best_point is not None:
                selected_points.append(best_point)
            else:
                break
        return selected_points
    
    #%% Particle Statistics
    # =============================================================================
    def attach_stats_ledger(self, ledger, ledger_path: Optional[str] = None) -> None:
        """Attach a SampleLedger whose ``particles`` list stores PSD ParticleInfo records."""
        self._ledger = ledger
        self._ledger_path = ledger_path
        if ledger is not None:
            self.analyzed_pars = ledger.particles

    def _persist_stats_ledger(self) -> None:
        """Write the attached particle-stats ledger to disk when available."""
        if self._ledger is None:
            return
        ledger_path = self._ledger_path
        if ledger_path is None and self.results_dir is not None:
            output_dir = (
                os.path.dirname(self.results_dir)
                if os.path.basename(self.results_dir) == cnst.IMAGES_DIR
                else self.results_dir
            )
            ledger_path = os.path.join(output_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
        if ledger_path is not None:
            self._ledger.to_json_file(ledger_path)

    def _apply_size_band(self, d_min_um: float, d_max_um: float) -> float:
        """Update area filters and search FOV for one diameter band; return FOV in mm."""
        min_area = diameter_um_to_area_um2(d_min_um)
        max_area = diameter_um_to_area_um2(d_max_um)
        fw_um = frame_width_um_from_max_area(max_area)
        self.powder_meas_cfg.min_area_par = min_area
        self.powder_meas_cfg.max_area_par = max_area
        self.powder_meas_cfg.par_search_frame_width_um = fw_um
        return fw_um / 1000.0

    def _collect_particle_stats_until(self, n_par_target_total: int) -> None:
        """Scan frames until ``analyzed_pars`` reaches ``n_par_target_total`` or frames end."""
        while len(self.analyzed_pars) < n_par_target_total:
            previous_n_par = len(self.analyzed_pars)
            par_were_found = self._move_and_get_particles_stats_in_frame()
            if par_were_found is None:
                logger.warning(
                    f"⚠️ Could not find {n_par_target_total} particles. "
                    f"Completed statistics using {len(self.analyzed_pars)} particles."
                )
                break
            if par_were_found is False:
                if self.verbose:
                    logger.warning("⚠️ No particle was found in this frame")
            else:
                if self.verbose:
                    tot_n_par_found = len(self.analyzed_pars)
                    n_par_found_frame = tot_n_par_found - previous_n_par
                    plural_frame = n_par_found_frame != 1
                    plural_total = tot_n_par_found != 1
                    logger.info(
                        f"{n_par_found_frame} particle{'s' if plural_frame else ''} "
                        f"{'were' if plural_frame else 'was'} found in this frame."
                    )
                    logger.info(
                        f"A total of {tot_n_par_found} particle{'s' if plural_total else ''} "
                        f"{'have' if plural_total else 'has'} now been analyzed. "
                        f"{n_par_target_total - tot_n_par_found} more to go."
                    )

    def get_particle_stats(self, n_par_target):
        """
        Analyze frames to collect and save particle size statistics.

        Large diameter ranges are automatically split into decade bands (largest
        first). Each band uses an FOV matched to its maximum particle size.
        Top-level frames use normal labels (``A0``); finer bands use hierarchical
        ids (``A0_a0``) for sub-frames inside visited parent FOVs.

        Parameters
        ----------
        n_par_target : int
            Target number of particles to analyze **per size band**.

        Returns
        -------
        pandas.DataFrame or None
            Summary statistics DataFrame, or None if no particles were found.
        """
        d_min_um = 2.0 * np.sqrt(self.powder_meas_cfg.min_area_par / np.pi)
        d_max_um = 2.0 * np.sqrt(self.powder_meas_cfg.max_area_par / np.pi)
        bands = split_diameter_range_um(d_min_um, d_max_um)

        is_manual = bool(
            getattr(getattr(self.EM, "measurement_cfg", None), "is_manual_navigation", False)
        )
        navigator = getattr(self.EM, "frame_navigator", None)

        if is_manual or navigator is None or len(bands) == 1:
            if len(bands) == 1 and navigator is not None and not is_manual:
                fw_mm = self._apply_size_band(bands[0][0], bands[0][1])
                if (
                    navigator.grid_search_fw_mm is None
                    or not np.isclose(float(navigator.grid_search_fw_mm), fw_mm, rtol=1e-6, atol=1e-9)
                ):
                    navigator.configure_top_level_particle_frames(fw_mm)
                else:
                    navigator.reset_frame_cursor()
            self._collect_particle_stats_until(int(n_par_target))
        else:
            parent_frames = None
            for band_idx, (d_lo, d_hi) in enumerate(bands):
                fw_mm = self._apply_size_band(d_lo, d_hi)
                if self.verbose:
                    print_double_separator()
                    logger.info(
                        f"Size band {band_idx + 1}/{len(bands)}: "
                        f"{d_lo:.4g}–{d_hi:.4g} µm "
                        f"(FOV {fw_mm * 1000:.3g} µm)"
                    )
                if band_idx == 0:
                    navigator.configure_top_level_particle_frames(fw_mm)
                else:
                    if not parent_frames:
                        logger.warning(
                            "⚠️ No parent frames from the previous size band; "
                            "skipping remaining finer bands."
                        )
                        break
                    navigator.configure_particle_subframes(parent_frames, fw_mm)

                band_target = len(self.analyzed_pars) + int(n_par_target)
                self._collect_particle_stats_until(band_target)
                parent_frames = navigator.get_visited_band_frames()
                self._persist_stats_ledger()

        n_par_analysed = len(self.analyzed_pars)
        if n_par_analysed == 0:
            logger.error(
                '❌ Could not find any particle. Please check your sample, '
                'or change the contrast/brightness values.'
            )
            return None

        self._persist_stats_ledger()
        return self.save_particle_statistics()

    def save_particle_statistics(self, output_file_suffix=''):
        """
        Process ParticleInfo records, compute summary statistics,
        export results, and produce a particle size histogram.

        Parameters
        ----------
        output_file_suffix : str, optional
            String added to output file name

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the calculated particle size statistics.
        """
        if os.path.basename(self.results_dir) == cnst.IMAGES_DIR:
            output_dir = os.path.dirname(self.results_dir)
        else:
            output_dir = self.results_dir

        if not self.analyzed_pars:
            raise ValueError("No particles available to save statistics from")

        par_IDs = np.array([p.id for p in self.analyzed_pars], dtype=int)
        frame_labels = np.array(
            [p.frame_id if p.frame_id is not None else "" for p in self.analyzed_pars],
            dtype=str,
        )
        areas_um = np.array(
            [float(p.area_um) for p in self.analyzed_pars if p.area_um is not None],
            dtype=float,
        )
        if len(areas_um) != len(self.analyzed_pars):
            raise ValueError("All ParticleInfo records must have area_um set for PSD export")

        eq_diam_um = np.array(
            [
                float(p.eq_diameter_um)
                if p.eq_diameter_um is not None
                else float(2.0 * np.sqrt(float(p.area_um) / np.pi))
                for p in self.analyzed_pars
            ],
            dtype=float,
        )

        n_par_analysed = len(areas_um)

        if n_par_analysed > 1 and np.isclose(areas_um.min(), np.partition(areas_um, 1)[1]):
            logger.warning(
                '⚠ The 2 smallest particles have identical area.\n'
                '   This may indicate they are only 1–2 pixels in size.\n'
                '   Frame width is set based on maximum acceptable particle size.\n'
                '   Consider reducing the maximum particle size so that the minimum\n'
                '   acceptable size is within ~1 order of magnitude.'
            )

        particle_data = pd.DataFrame({
            cnst.PAR_ID_DF_KEY: par_IDs,
            cnst.FRAME_ID_DF_KEY: frame_labels,
            cnst.PAR_AREA_UM_KEY: areas_um,
            cnst.PAR_EQ_D_KEY: eq_diam_um,
        })
        particle_data.to_csv(
            os.path.join(output_dir, f"{self._sample_id}_{cnst.PARTICLE_SIZES_FILENAME}{output_file_suffix}.csv"),
            header=True,
            index=False,
        )

        par_size_distr = {
            'Measure': 'Equivalent particle diameter in μm',
            'n_par_analysed': n_par_analysed,
            'mean': np.mean(eq_diam_um),
            'stdev': np.std(eq_diam_um),
            'median': np.median(eq_diam_um),
            'max': np.max(eq_diam_um),
            'min': np.min(eq_diam_um),
            'D10': np.percentile(eq_diam_um, 10),
            'D25': np.percentile(eq_diam_um, 25),
            'D75': np.percentile(eq_diam_um, 75),
            'D90': np.percentile(eq_diam_um, 90),
        }

        stats_df = pd.DataFrame(par_size_distr, index=[0])
        stats_df.to_csv(
            os.path.join(output_dir, f"{self._sample_id}_{cnst.PARTICLE_STATS_FILENAME}{output_file_suffix}.csv"),
            header=True,
            index=False,
        )

        if self.verbose:
            print_double_separator()
            logger.info("📊 Computed statistics:")
            logger.info(stats_df.T.to_string(header=False))

        self._save_particle_size_histogram(
            eq_diam_um,
            results_dir=output_dir,
            _sample_id=self._sample_id,
            verbose=self.verbose,
            output_file_suffix=output_file_suffix,
        )

        return stats_df


    def _save_particle_size_histogram(self, diameters_um, results_dir=None, _sample_id=None, verbose=False, output_file_suffix = '', bins=20):
        """
        Save a histogram plot of particle equivalent diameters.
    
        This function generates and saves a histogram of the particle size distribution.
        The plot is saved as a PNG file in the specified results directory, with the file name based on the sample ID.
        Optionally, the plot can be displayed interactively.
    
        Parameters
        ----------
        diameters_um : array-like
            Array of equivalent particle diameters in micrometers (μm).
        results_dir : str, optional
            Directory where the histogram PNG file will be saved. If None, uses `self.results_dir`.
        _sample_id : str, optional
            Identifier for the sample, used in the output file name. If None, uses `self._sample_id`.
        verbose : bool, optional
            If True, displays the plot interactively. Default is False.
        output_file_suffix : str, optional
            String added to output file name
        bins : int, optional
            Number of bins to use in the histogram. Default is 20.
    
        Returns
        -------
        None
    
        Notes
        -----
        - The histogram is always saved as a PNG file with the name '{_sample_id}_Par_size_distribution_hist.png'.
        - The function uses matplotlib for plotting.
    
        Potential Improvements / TODO
        -----------------------------
        - Allow user to specify output file format (e.g., SVG, PDF).
        - Add option to overlay statistics (mean, median, percentiles) on the plot.
        - Enable saving both linear and logarithmic scale histograms.
        - Allow customization of plot colors and style.
        - Return the matplotlib Figure object for further manipulation if desired.
        - Add error handling for file I/O issues.
        """
        if results_dir is None:
            results_dir = self.results_dir
        if _sample_id is None:
            _sample_id = self._sample_id
    
        plt.figure()
        plt.hist(diameters_um, bins=bins, edgecolor='black')
        plt.xlabel('Equivalent Diameter (μm)')
        plt.ylabel('Counts')
        plt.title('Particle size distribution')
        out_path = os.path.join(results_dir, f'{_sample_id}_{cnst.PARTICLE_STAT_HIST_FILENAME}{output_file_suffix}.png')
        plt.savefig(out_path)
        plt.close()
        
    
    def _move_and_get_particles_stats_in_frame(self, frame_image=None, pixel_size=None, results_dir=None):
        """
        Move to the next frame and analyze it to detect and characterize particles on the substrate.
    
        This function either acquires a new frame from the electron microscope (if running at the EM)
        or processes a provided image. It detects particles, filters them by size and edge proximity,
        calculates their areas in pixels and micrometers squared, and optionally saves annotated images
        with detected particles and their indices.
    
        Parameters
        ----------
        frame_image : np.ndarray, optional
            Grayscale image of the current frame. If not provided and running at the EM, the image is acquired.
        pixel_size : float, optional
            Pixel size in micrometers. Required if not running at the microscope.
        results_dir : str, optional
            Directory to save result images and masks.
    
        Returns
        -------
        bool or None
            Returns True if at least one valid particle is found in the frame,
            False if no valid particles are detected,
            or None if there are no more frames available (when running at the EM).
    
        Notes
        -----
        - If running at the EM, a new frame is acquired and saved if development_mode=True.
        - If not at the EM, both `frame_image` and `pixel_size` must be provided.
        - The function applies a Gaussian blur to suppress noise before particle detection.
        - Particles are filtered by area and by proximity to the frame edge.
        - The area of each accepted particle is appended to `self.analyzed_pars`.
        - An annotated image with detected particles and their indices is saved for visualization.
        - If no valid particles are found, the mask image is deleted and the function returns False.
        - If at least one particle is found, the microscope focus/contrast/brightness is refreshed if needed.
        """
    
        if EM_driver.is_microscope_connected():
            self._check_EM_controller_initialization()
            move_to_frame_success = self.EM.go_to_next_frame()
            if move_to_frame_success:
                # Collect image
                frame_image = EM_driver.get_image_data(self._im_width, self._im_height, 1)
                if self.development_mode:
                    cv2.imwrite(os.path.join(self.results_dir, self._sample_id + f'_fr_{self.EM.current_frame_label}.png'), frame_image)
            else:
                # No more frames are available
                return None
        elif frame_image is not None and pixel_size is not None:
            self._im_height, self._im_width = frame_image.shape
            self.EM.pixel_size_um = pixel_size
            if results_dir:
                self.results_dir = results_dir
        else:
            raise ValueError('This function must be run at the microscope, or it needs to be passed both image and its pixel size')
        
        # Select particles on substrate
        blurred_image = cv2.GaussianBlur(frame_image, (5, 5), 0)
        
        # Get mask of particles on substrate
        save_mask_img = True
        par_mask, mask_img_path = self._get_particles_on_substrate_mask(blurred_image, save_image=save_mask_img)
    
        # Find connected components
        num_labels, labels, stats, centroids = self._get_connected_components_with_stats(par_mask)
        
        # Store particle area / geometry as ParticleInfo
        par_centroids = []  # Only used for saving image
        par_areas = []
        par_cntr = 0
        frame_label = str(self.EM.current_frame_label) if self.EM.current_frame_label is not None else None
        for i in range(1, num_labels):  # Skip the background component (index 0)
            par_area_pixels = stats[i, cv2.CC_STAT_AREA]
            if self._is_particle_area_ok(par_area_pixels) and not self.is_particle_at_frame_edge(stats, i):
                # If particle is within size limits and is not at the edge, consider it in the statistics
                par_area_um = float(par_area_pixels * self.EM.pixel_size_um**2)
                eq_diameter_um = float(2.0 * np.sqrt(par_area_um / np.pi))
                abs_mm = self.EM.convert_pixel_pos_to_mm(centroids[i])
                global_par_id = len(self.analyzed_pars)
                self.analyzed_pars.append(
                    ParticleInfo(
                        id=global_par_id,
                        area_um=par_area_um,
                        eq_diameter_um=eq_diameter_um,
                        coordinates=Coordinate2D(x=float(abs_mm[0]), y=float(abs_mm[1])),
                        frame_id=frame_label,
                    )
                )
                # Store stats locally to draw circles
                par_areas.append(par_area_pixels)
                par_centroids.append(centroids[i])
                par_cntr += 1
                
        # Save image to visualize selected particles
        text_margin = 30  # pixel margin to determine where to annotate image with particle numbers
        text_pos_scale = 0.9
        im_annotations = []
        if par_cntr > 0:
            first_par_n = len(self.analyzed_pars) - par_cntr
            for i, (center, area) in enumerate(zip(par_centroids, par_areas)):
                ann_dict = {}
                radius_pixel = int(np.sqrt(area / np.pi) * 1.1)  # Equivalent radius for particle (as a circle)
                
                ann_dict[cnst.ANNOTATION_CIRCLE_KEY] = (radius_pixel, center.astype(int), 2)
                
                x_pos_text = int(center[0] + radius_pixel * text_pos_scale)  # Number on the top-right of the circle
                y_pos_text = int(center[1] - radius_pixel * text_pos_scale)
                if x_pos_text > (self._im_width - text_margin) or y_pos_text < text_margin: 
                    # Move number to center
                    x_pos_text = int(center[0])
                    y_pos_text = int(center[1])
                ann_dict[cnst.ANNOTATION_TEXT_KEY] = (str(first_par_n + i), (x_pos_text, y_pos_text))
                im_annotations.append(ann_dict)
                
            filename = f"{self._sample_id}_fr{self.EM.current_frame_label}"
            # Save annotated particle image
            self.EM.save_frame_image(filename + '_particles', im_annotations = im_annotations, frame_image = frame_image)
            
            # Save mask image
            par_mask_to_save = self.prepare_mask_for_visualization(par_mask)
            self.EM.save_frame_image(filename + '_par_mask', im_annotations = im_annotations, frame_image = par_mask_to_save)

        
        # Return True if at least 1 particle was found
        if par_cntr == 0:
            if save_mask_img and os.path.exists(mask_img_path):
                os.remove(mask_img_path)
            return False
        else:
            if not self.EM.measurement_cfg.is_manual_navigation and (time.time() - self.EM.microscope_ctrl._last_EM_adjustment_time > self.EM.microscope_ctrl.refresh_time):
                # Adjust EM focus, contrast and brightness, but only if a particle is actually present
                self.EM.adjust_BCF()
            return True
        