#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated Electron Microscopy (EM) Controller and Sample Image Analyzer

This module provides classes and functions for automated analysis and acquisition in scanning electron microscopy (SEM),
including particle detection, stage navigation, X-ray spectra (EDS, WDS) acquisition, and sample management.
It is designed to interface with an EM driver and streamline the workflow from image acquisition to data export.

Main Classes
------------
EM_Controller
    Automates particle detection, mask generation, and X-ray spectra acquisition on detected particles.
    Supports both fully automated and manual (user-guided) collection modes.
    Provides methods for collecting particle statistics and managing frame navigation.

EM_Sample_Finder
    Provides utilities for locating and managing samples within the EM.
    Currently includes only the detection of the center of a Carbon tape.

Example Usage
-------------
# EDS spot acquisition workflow:
>>> particle_finder = EM_Particle_Finder(...)
>>> particle_finder.initialise_SEM()
>>> while True:
...     particle_finder.go_to_next_particle()
...     xy_spot_list = particle_finder.get_XS_acquisition_spots_coord_list()
...     # Collect X-ray spectra at each (x, y) as needed


# Particle analysis workflow:
>>> from em_analysis import EM_Particle_Finder, EM_Sample_Finder
>>> microscope_cfg = {...}
>>> measurement_cfg = {...}
>>> holder_cfg = {...}
>>> particle_finder = EM_Particle_Finder(
...     sample_ID="Sample1",
...     microscope_d=microscope_cfg,
...     measurement_d=measurement_cfg,
...     sample_holder_d=holder_cfg,
...     results_dir="./results"
... )
>>> particle_finder.initialise_SEM()
>>> particle_finder.get_particle_stats(n_par_target=500)

# C tape detection workflow:
>>> from em_analysis import EM_Sample_Finder
>>> finder = EM_Sample_Finder('MySEM', (0, 0), 3, 12, './results', verbose=True)
>>> Ctape_coords = finder.detect_Ctape()
>>> if Ctape_coords:
...     print("Center (mm):", Ctape_coords[0], "Half-width (mm):", Ctape_coords[1])

Notes
-----
- Requires a working EM_driver and appropriate hardware configuration, with the following API functions:

    Microscope & Driver Setup
    -------------------------
    - load_microscope_driver(microscope_ID)
        Load and initialize the instrument-specific driver.
    
    - activate()
        Wake up and activate the electron microscope.
    
    - to_SEM()
        Switch the microscope to SEM (scanning electron microscopy) mode.
    
    - set_electron_detector_mode(mode)
        Set the electron detector mode (e.g., 'All', 'BSD').
    
    - set_high_tension(voltage)
        Set the accelerating voltage (high tension) for the electron beam.
    
    - set_beam_current(current)
        Set the electron beam current.
    
    - get_range_frame_width()
        Query the allowed minimum and maximum frame widths (in mm).
    
    - set_frame_width(width)
        Set the field of view/frame width (in mm).
    
    - adjust_focus(wd)
        Set the working distance/focus to a specific value (in mm).
    
    Image & Navigation
    ------------------
    - move_to(x, y)
        Move the microscope stage to the specified (x, y) position.
    
    - get_frame_width()
        Get the current frame width (in mm).
    
    - get_navigation_camera_image()
        Acquire an image from the navigation camera.
    
    - auto_focus()
        Run the microscope's autofocus routine and return the resulting working distance.
    
    - auto_contrast_brightness()
        Automatically adjust the microscope's brightness and contrast.
    
    - set_brightness(value)
        Set the frame brightness to a specific value.
    
    - set_contrastss(value)
        Set the frame contrast to a specific value.
        (NOTE: This may be a typo; check if the correct function is set_contrast.)
    
    Spectroscopy & Analysis
    -----------------------
    - get_EDS_analyser_object()
        Create and return an EDS (energy-dispersive X-ray spectroscopy) analyzer object.
    
    - acquire_XS_spectral_data(analyzer, x, y, max_time, target_counts)
        Acquire an X-ray spectrum at the specified position with given analyzer and parameters.
    
    Calibration & Attributes (used as properties/constants)
    -------------------------------------------------------
    - image_to_stage_coords_transform
        Transformation matrix for converting image to stage coordinates.
    
    - navcam_im_w_mm
        Navigation camera image width in mm.
    
    - navcam_x_offset
        Navigation camera X offset.
    
    - navcam_y_offset
        Navigation camera Y offset.
    
    - microscope_calib_dir
        Directory containing microscope calibration files.
    
    - is_at_EM
        Boolean indicating if the code is running at the microscope.
    
    - typical_wd
        Typical working distance (mm).
    
    - im_width
        Image width in pixels.
    
    - im_height
        Image height in pixels.


Created on Wed Jul 31 09:28:07 2024

@author: Andrea
"""
# Standard library imports
import os
import time
import warnings
import json

# Third-party imports
import cv2
import numpy as np
from PIL import Image

# Typing (Python 3.5+)
from typing import Any, List, Optional, Tuple, Union

# Local project imports
import autoemxsp.tools.constants as cnst
from autoemxsp.tools.utils import (
    AlphabetMapper,
    Prompt_User,
    EMError,
    print_single_separator,
    draw_scalebar
)
from autoemxsp.tools.config_classes import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    PowderMeasurementConfig,
    BulkMeasurementConfig
)
from autoemxsp import EM_driver
from autoemxsp.core.EM_particle_finder import EM_Particle_Finder


#%% Electron Microscope Particle Finder class    
class EM_Controller:
    an_circle_key = 'circle'
    an_text_key = 'text'
    """
    Electron Microscope (EM) controller for automated particle analysis and X-ray spectra (EDS, WDS) acquisition.

    This class provides methods for automated image acquisition, stage navigation, and X-ray spectral acquisition.
    It is designed to interface with a microscope driver (`EM_driver`) for stage movement and image capture.

    Configuration is provided via structured dataclasses:
        - MicroscopeConfig
        - SampleConfig
        - MeasurementConfig
        - SampleSubstrateConfig

    Main Methods
    ------------
    initialise_SEM()
        Wakes up the SEM microscope, sets measurement parameters, and evaluates locations to scan for particles.
    initialise_XS_analyzer()
        Gets EDS/WDS analyzer object from EM_driver, and optionally updates electron beam energy.
        Currently only EDS is supported.
    acquire_XS_spot_spectrum(x,y) 
        Acquire and return X-Ray spectral data in position (x,y) in current stage position.
        (x,y) are relative coordinates. See function description for more information.
        Currently only EDS spectra are supported.
        
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
    sample_cfg : SampleConfig
        Sample configuration (see dataclass for details).
    microscope_cfg : MicroscopeConfig
        Microscope configuration (see dataclass for details).
    measurement_cfg : MeasurementConfig
        Measurement/acquisition configuration (see dataclass for details).
    sample_substrate_cfg : SampleSubstrateConfig
        Sample substrate configuration (see dataclass for details).
    powder_meas_cfg : PowderMeasurementConfig
        Configuration for powder measurement.
    bulk_meas_cfg : BulkMeasurementConfig
        Configuration for measurements of bulk or bulk-like samples.
    init_wd : float
        Initial working distance in mm (used to limit autofocus range, set at EM driver).
    im_width : int
        Image width in pixels (set at EM driver).
    im_height : int
        Image height in pixels (set at EM driver).
    grid_search_fw_mm : float
        Frame width (in mm) used during initial grid search.
    pixel_size_um : float or None
        Pixel size in micrometers of current image.
    current_frame_label : str
        Label of currently analyzed frame. Initialized as an empty string for testing
    refresh_time : float
        Time in seconds after which autofocus and brightness are refreshed.
    frame_labels : list
        List of frame labels or identifiers.
    results_dir : str or None
        Directory for saving result images and data.
    verbose : bool
        If True, print progress and information to the console.
    development_mode (bool): Whether class is being used for testing image processing functions, without real-time acquisition.
        Enables use of class outside microscope enviroment

    Internal Attributes
    -------------------
    _center_pos : tuple of float
        Sample center position (x, y) on the microscope stage (from sample_cfg).
    _sample_hw_mm : float
        Half-width of the sample region in mm (from sample_substrate_cfg).
    _min_wd : float
        Minimum allowed working distance in mm.
    _max_wd : float
        Maximum allowed working distance in mm.
    _frame_cntr : int
        Counter for the number of frames analyzed.
    _current_pos : tuple of float
        Current (x, y) position in the sample coordinate system.
    _last_EM_adjustment_time : float
        Timestamp of the last EM adjustment (autofocus/brightness)

    Notes
    -----
    - This class assumes that the global `EM_driver` object is available and properly configured.
    - Configuration is provided via dataclasses. All validation of configuration fields is performed in the dataclasses.
    - Attributes with a leading underscore are for internal use and should not be accessed directly by users.
    """
    def __init__(
        self,
        microscope_cfg: MicroscopeConfig,
        sample_cfg: SampleConfig,
        measurement_cfg: MeasurementConfig,
        sample_substrate_cfg: SampleSubstrateConfig,
        powder_meas_cfg: PowderMeasurementConfig,
        bulk_meas_cfg: BulkMeasurementConfig,
        init_fw: float = 0.5,
        results_dir: Optional[str] = None,
        verbose: bool = True,
        development_mode: Optional[bool] = False, 
    ):
        """
        Initialize an EM_Particle_Finder object for automated/manual electron microscopy particle analysis.
    
        Sets up all parameters for automated particle analysis and X-ray spectra acquisition,
        including sample, microscope, measurement, and sample substrate settings. Also initializes
        options for image acquisition, particle selection, and EDX spot collection, as well as
        internal counters and storage for results.
    
        Parameters
        ----------
        microscope_cfg : MicroscopeConfig, optional
            Microscope configuration dataclass instance.
        sample_cfg : SampleConfig, optional
            Sample configuration dataclass instance.
        measurement_cfg : MeasurementConfig, optional
            Measurement/acquisition configuration dataclass instance.
        sample_substrate_cfg : SampleSubstrateConfig, optional
            Sample substrate configuration dataclass instance.
        powder_meas_cfg : PowderMeasurementConfig
            Configuration for powder measurement.
        bulk_meas_cfg : BulkMeasurementConfig
            Configuration for measurements of bulk or bulk-like samples.
        init_fw : float, optional
            Initial frame width in mm for grid search (default: 0.5).
        results_dir : str, optional
            Directory to save result images and data (default: None).
        verbose : bool, optional
            If True, print progress and information to the console (default: True).
        development_mode (bool): Whether class is being used for testing image processing functions, without real-time acquisition.
            Enables use of class outside microscope enviroment
    
        Raises
        ------
        RuntimeError
            If the microscope driver cannot be loaded.
    
        Notes
        -----
        - All configuration validation is performed in the respective dataclasses.
        - This initializer assumes that configuration dataclasses are valid and complete.
        - Additional internal attributes are initialized for frame tracking and image navigation.
        """

        # --- Sample and system characteristics
        self.sample_cfg = sample_cfg
        self._center_pos = sample_cfg.center_pos
        self.microscope_cfg = microscope_cfg
        self.measurement_cfg = measurement_cfg
        self.sample_substrate_cfg = sample_substrate_cfg
        self._sample_hw_mm = sample_cfg.half_width_mm
        self.powder_meas_cfg = powder_meas_cfg
        self.bulk_meas_cfg = bulk_meas_cfg
        
        # --- Load microscope driver for instrument microscope_ID
        try:
            EM_driver.load_microscope_driver(microscope_cfg.ID)
        except Exception as e:
            raise RuntimeError(f"Failed to load microscope driver: {e}")
        if not development_mode:
            if not EM_driver.is_at_EM:
                raise EMError("Instrument driver could not be loaded")

        # Min and max working distance. To avoid gross failures of autofocus algorithm
        if isinstance(measurement_cfg.working_distance, float):
            self.init_wd = measurement_cfg.working_distance
        else:
            self.init_wd = EM_driver.typical_wd
        self._min_wd = self.init_wd - 1  # in mm
        self._max_wd = self.init_wd + 1  # in mm

        # Image width and height in pixels
        self.im_width = EM_driver.im_width
        self.im_height = EM_driver.im_height

        # Frame width employed during initial grid search, in mm
        self.grid_search_fw_mm = init_fw  # fallback for offline/test mode

        # Time after which auto focus and brightness are refreshed
        self.refresh_time = 1  # 180 secs (i.e., 3 min)

        # --- General options
        self.development_mode = development_mode
        self.results_dir = results_dir
        self.verbose = verbose

        # --- Variable initializations
        self.is_initialized = False
        self.pixel_size_um: Optional[float] = None
        self.frame_labels: List[Any] = []
        self.current_frame_label: str = ''
        self._frame_cntr = 0  # Keeps track of number of frames analysed
        self._current_pos = self._center_pos
        self._bulk_offset_cntr = 0
        self._last_EM_adjustment_time: float = 0.0
    
    #%% Microscope initialization, functions called only once
    # =============================================================================
    def initialise_SEM(
        self, 
    ) -> None:
        """
        Activate and configure the Scanning Electron Microscope (SEM) for particle finding.
    
        This method performs the following steps:
            1. Wakes up the SEM if necessary.
            2. Switches to SEM mode.
            3. Sets the electron detector mode to backscattered electrons (BSD).
            4. Sets the beam voltage and current according to the current measurement mode calibration.
            5. Initializes the EDS analyzer object.
            6. Moves to sample center, sets frame width, and working distance.
            7. Adjusts focus, brightness, and contrast.
            8. Calculates frame centers for searching.
            9. Initializes particle/frame counters.
    
        Raises
        ------
        KeyError
            If the current measurement mode is not found in the calibration dictionary.
        EMError
            If an error occurs during SEM activation or configuration.
    
        Notes
        -----
        This method assumes that EM_driver is an initialized and accessible object that provides
        the necessary SEM control methods.
    
        Potential improvement for SEM initialization:
            - Fix initial brightness and contrast so that C tape is fully below threshold.
            - Go to first frame where it finds particles.
            - Go to first (or biggest) particle and adjust brightness and contrast.
            - Go to second, and cross-check values.
            - If equal, use these brightness and contrast values.
            - If encountering multiple contamination particles, repeat this.
            - Revert to these values every time large scale particle segmentation is needed.
            - Alternatively, brightness and contrast could be passed manually.
        """
        if self.verbose:
            print_single_separator()
            print("Activating SEM, and setting up...")
    
        try:
            # Wake up SEM if necessary
            EM_driver.activate()
    
            # Switch to SEM mode
            EM_driver.to_SEM()
    
            # Set detector type to BSD (Backscattered electron detector)
            EM_driver.set_electron_detector_mode(self.microscope_cfg.detector_type)
    
            # Set beam voltage (high tension) for EDS collection
            if self.measurement_cfg.beam_energy_keV:
                EM_driver.set_high_tension(self.measurement_cfg.beam_energy_keV)
            else:
                warnings.warn("No acceleration voltage was provided via measurement_cfg.beam_energy_keV. Using current microscope configurations",
                              UserWarning)
            
            # Set beam current for EDS collection
            if self.measurement_cfg.beam_current:
                EM_driver.set_beam_current(self.measurement_cfg.beam_current)
            else:
                warnings.warn("No beam current was provided via measurement_cfg.beam_current. Using current microscope configurations",
                              UserWarning)
    
            # Move to sample center
            self.move_to_pos(self._center_pos)
    
            # Set working distance (needed for reliable autofocus)
            EM_driver.adjust_focus(self.init_wd)  # in mm
    
            # Adjust focus, brightness, and contrast
            if self.verbose:
                print_single_separator()
                print("Adjusting contrast, brightness, and focus.")
    
            # If forcing brightness and contrast does not work, algorithm is run with auto brightness and contrast adjustments
            self.adjust_BCF()
            
            # Update for external access
            self.is_initialized = True
    
            if self.verbose:
                print("SEM initialisation completed.")
    
        except KeyError:
            raise
        except Exception as e:
            # Only wrap truly unexpected errors.
            raise EMError(f"Error during SEM activation: {e}") from e
            
    
    def initialise_sample_navigator(self,
        EM_controller,
        exclude_sample_margin: bool = True
    ) -> None:
        """
        Initialize the sample navigator for automated spectra collection.
    
        Only 'powder' and 'bulk' samples are currently supported for automated particle detection and navigation.
        For other sample types, use manual navigation, setting measurement_cfg.is_manual_navigation = True.
        
        Parameters
        ----------
        EM_controller: EM_Controller object
            Required for analysis of samples of type = 'powder'
        exclude_sample_margin : bool, optional
            Whether to use a margin when calculating frame centers. Reduces explored sample area,
            but avoids sample edges, which may be contaminated, or affected by larger substrate signal.
            (default: True).
        
        Raises
        ------
        NotImplementedError
            If the sample type is not 'powder' and manual navigation is not enabled.
    
        """
        if self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE:
            # Set frame width, and update current pixel size
            if getattr(EM_driver, "is_at_EM", True):
                min_fw, max_fw = EM_driver.get_range_frame_width()
                self.grid_search_fw_mm = np.clip(self.powder_meas_cfg.par_search_frame_width_um /1000, min_fw, max_fw)
            self.set_frame_width(self.grid_search_fw_mm)
            
            # Calculate frame centers for particle search
            im_h_to_w_ratio = self.im_height / self.im_width
            self._calc_frame_centers(
                horizontal_spacing_mm=self.grid_search_fw_mm,
                im_h_to_w_ratio = im_h_to_w_ratio,
                center_pos = self._center_pos,
                randomize_frames = True,
                exclude_sample_margin = True
            )
            
            # Initialise particle finder for powder samples
            self.particle_finder = EM_Particle_Finder(
                EM_controller,
                powder_meas_cfg=self.powder_meas_cfg,
                is_manual_particle_selection=self.powder_meas_cfg.is_manual_particle_selection,
                results_dir=self.results_dir,
                verbose=self.verbose,
                development_mode=self.development_mode
            )
        elif self.sample_cfg.type == cnst.S_BULK_SAMPLE_TYPE:
            if getattr(EM_driver, "is_at_EM", True):
                min_fw, max_fw = EM_driver.get_range_frame_width()
                self.grid_search_fw_mm = np.clip(self.bulk_meas_cfg.image_frame_width_um / 1000, min_fw, max_fw)
            self.set_frame_width(self.grid_search_fw_mm)
            # Construct grid of acquisition spots
            self._calc_bulk_grid_acquisition_spots()
            

        elif not self.measurement_cfg.is_manual_navigation:
            raise NotImplementedError(
                "Sample type '{}' is not supported for automated composition analysis. "
                "Use measurement_cfg.is_manual_navigation = True.".format(self.sample_cfg.type)
            )
        # Save image of initial location to show
        initial_image = EM_driver.get_image_data(self.im_width, self.im_height, 1)
        draw_scalebar(initial_image, self.pixel_size_um)
        cv2.imwrite(os.path.join(self.results_dir, cnst.INITIAL_SEM_IM_FILENAME + '.png'), initial_image)
        
    
    
    def _calc_bulk_grid_acquisition_spots(self) -> bool:
        """
        Calculate and apply an offset to the center position for bulk grid acquisition,
        then construct a square grid of acquisition spots.
        
        The offset increases in discrete steps, cycling through three directions:
        (x, 0), (x, x), and (0, x), where x is the offset distance.
        If the offset distance exceeds the grid spot spacing, no further acquisition is performed.
    
        Returns:
            bool: True if the grid was constructed, False if offset exceeds grid spacing.
        """
        # Calculate offset distance in micrometers
        offset_dist_um = (
            self.bulk_meas_cfg.min_xsp_spots_distance_um * np.ceil(self._bulk_offset_cntr / 3)
        )
    
        # Determine offset direction (cycles through 3 directions)
        offset_dir_id = self._bulk_offset_cntr % 3
    
        if offset_dist_um > self.bulk_meas_cfg.grid_spot_spacing_um:
            # Offset is larger than the grid spot spacing; do not proceed
            return False
    
        # Define offset coordinates based on direction
        offset_dist_mm = offset_dist_um / 1000
        if offset_dir_id == 0:
            offset_coords = (offset_dist_mm, 0)
        elif offset_dir_id == 1:
            offset_coords = (offset_dist_mm, offset_dist_mm)
        else:  # offset_dir_id == 2
            offset_coords = (0, offset_dist_mm)
    
        # Apply offset to center position
        center_pos = tuple(np.array(self._center_pos) + np.array(offset_coords))
    
        # Construct grid of acquisition spots (always square grid)
        self._calc_frame_centers(
            horizontal_spacing_mm=self.bulk_meas_cfg.grid_spot_spacing_um / 1000,
            im_h_to_w_ratio=1,
            center_pos=center_pos,
            randomize_frames=self.bulk_meas_cfg.randomize_frames,
            exclude_sample_margin=self.bulk_meas_cfg.exclude_sample_margin,
        )
    
        # Increment offset counter for next call
        self._bulk_offset_cntr += 1
        return True
        
    
    def _calc_frame_centers(self, horizontal_spacing_mm, im_h_to_w_ratio, center_pos, randomize_frames, exclude_sample_margin):
        '''
        Generates and labels a set of evenly spaced scanning locations (frames) within a sample area,
        either circular or rectangular. Optionally avoids the rough sample border (margin), and can 
        randomize the order of the frames to reduce spatial bias.
    
        In general, this function determines where the microscope should position itself to systematically 
        analyze a sample. It calculates a grid of (x, y) positions ("frame centers") covering either a 
        circular or rectangular region, depending on the sample_substrate_cfg.shape ensuring each position is
        within a safe distance from the edge if requested. Each frame center is given a unique label
        (like "A3" or "B7") for identification. The resulting positions and labels are stored for later navigation.
        
        Parameters
        ----------
        horizontal_spacing_mm: float
            Horizontal spacing between neighboring grid spots in mm
        im_h_to_w_ratio: float,
            ratio height/width of frame dimensions
        center_pos: tuple(float, float)
            (x,y) coordinates of sample center, in stage coordinates
        randomize_frames : bool
            If True, shuffle the order of frame centers to avoid spatial bias.
        exclude_sample_margin : bool
            If True, apply a margin to avoid the rough border of the C tape (for EDX).
            If False, allow frames up to the tape edge (for size measurements).
        shape : str, optional
            'circle' to scan within a circular region (default), or 'rectangle' for a rectangular region.
    
        Sets
        ----
        self.frame_pos_mm : list of tuple
            List of (x, y) positions (in mm) corresponding to the center of each frame.
        self.frame_labels : list of str
            List of corresponding frame labels (e.g., 'A0', 'B3', ...).
        self.num_frames : int
            Number of analyzable frames.
            
        Notes
        -----
        Frame labelling follows a grid system:
        
            - Columns are labelled with letters (A, B, C, ...), rows with numbers (0, 1, 2, ...).
            - Each frame center is assigned a label such as 'A0', 'B3', etc.
            - For a circular sample, only frames within the circle are included.
        
        Example grid:
        
              y
              ^
              |        +-----------------------------+
              |        | A0 | B0 | C0 | D0 | E0 |    |
              |        +-----------------------------+
              |        | A1 | B1 | C1 | D1 | E1 |    |
              |        +-----------------------------+
              |        | A2 | B2 | C2 | D2 | E2 |    |
              |        +-----------------------------+
              |        | A3 | B3 | C3 | D3 | E3 |    |
              |        +-----------------------------+
              |        | A4 | B4 | C4 | D4 | E4 |    |
              |        +-----------------------------+
              +-------------------------------------> x

        '''
        cx, cy = center_pos
    
        # Determine the usable sample half width, optionally removing margin 
        if exclude_sample_margin:
            # margin defined as the size of a single frame
            margin = horizontal_spacing_mm * np.sqrt(1 + (im_h_to_w_ratio)**2)
            sample_hw_mm = self._sample_hw_mm - margin
        else:
            sample_hw_mm = self._sample_hw_mm
    
        frame_centers = []
        frame_labels = []
        alphabet_mapper = AlphabetMapper()
    
        # Define region checker function
        if self.sample_substrate_cfg.shape == cnst.CIRCLE_SUBSTRATE_SHAPE:
            def is_inside_region(x, y):
                return (x - cx) ** 2 + (y - cy) ** 2 < sample_hw_mm ** 2
        elif self.sample_substrate_cfg.shape == cnst.SQUARE_SUBSTRATE_SHAPE:
            rect_left   = cx - sample_hw_mm
            rect_right  = cx + sample_hw_mm
            rect_top    = cy - sample_hw_mm * im_h_to_w_ratio
            rect_bottom = cy + sample_hw_mm * im_h_to_w_ratio
            def is_inside_region(x, y):
                return (rect_left <= x <= rect_right) and (rect_top <= y <= rect_bottom)
        else:
            raise ValueError(f"Sample substrate shape must be one among {self.sample_substrate_cfg.ALLOWED_SHAPES}")
            
        half_n_frames_x = int(self._sample_hw_mm / horizontal_spacing_mm) + 1
        half_n_frames_y = int(self._sample_hw_mm / (horizontal_spacing_mm * im_h_to_w_ratio)) + 1
        
        frame_centers = []
        frame_labels = []
        alphabet_mapper = AlphabetMapper()
        
        for i in range(-half_n_frames_x, half_n_frames_x + 1):
            label_letter = alphabet_mapper.get_letter(i + half_n_frames_x)
            for j in range(-half_n_frames_y, half_n_frames_y + 1):
                label = label_letter + str(j + half_n_frames_y)
                x = cx + i * horizontal_spacing_mm
                y = cy + j * horizontal_spacing_mm * im_h_to_w_ratio
        
                if is_inside_region(x, y):
                    frame_centers.append((x, y))
                    frame_labels.append(label)
    
        if randomize_frames and frame_centers:
            frames = list(zip(frame_centers, frame_labels))
            np.random.shuffle(frames)
            frame_centers, frame_labels = zip(*frames)
    
        self.frame_pos_mm = list(frame_centers)
        self.frame_labels = list(frame_labels)
        self.num_frames = len(frame_centers)
 

    #%% X-ray spectra acquisition functions
    # =============================================================================
    def initialise_XS_analyzer(self, beam_voltage: float = None):
        """
        Initialize the EDS (Energy Dispersive X-ray Spectroscopy) analyzer and optionally set the electron beam voltage.
        
        If beam_voltage remains unspecified, uses current voltage initially set during EM initialization.
    
        Parameters
        ----------
        beam_voltage : float, optional
            Desired EM beam voltage (in kilovolts, kV) for EDS acquisition. If provided, sets the high tension
            to this value. If not provided, the current beam voltage is used.
    
        Notes
        -----
        - This method creates an EDS analyzer object via the EM driver.
        - The beam voltage should be specified in kV (e.g., 20.0 for 20 kV).
        """
        # Create EDS analyzer object
        self.analyzer = EM_driver.get_EDS_analyser_object()
        
        if beam_voltage is not None:
            # Set beam voltage (high tension) for EDS collection (expects volts)
            EM_driver.set_high_tension(beam_voltage)  # Convert kV to V     
    
            
    def get_XSp_coords(
        self, n_tot_sp_collected: int
    ) -> Tuple[bool, Optional[List[Tuple[float, float]]], Optional[int]]:
        """
        Determine X-ray spectrum acquisition coordinates for the next spectrum.
    
        Depending on the navigation mode, either prompts the user to select a spot
        manually, or automatically selects the next particle spot for 'powder' samples.
    
        Parameters
        ----------
        n_tot_sp_collected : int
            The current total number of spectra collected (used for labeling and selection).
    
        Returns
        -------
        success : bool
            True if a spot was selected/found, False if user stopped or no more particles.
        spots_xy_list : list of tuple[float, float] or None
            List of (x, y) coordinates for acquisition spots, or None if unsuccessful.
        particle_cntr : int or None
            The particle counter/index, or None if not applicable.
    
        Notes
        -----
        - For manual navigation, the user is prompted to select the spot.
        - For powder samples, the next particle is selected automatically.
        - Assumes required attributes (e.g., self.im_width, self.particle_finder) are initialized.
        """
    
        if self.measurement_cfg.is_manual_navigation:
            # Prompt user to select the X-ray acquisition spot
            prompt = Prompt_User(
                title="Select X-Ray Acquisition Spot",
                message="Center image on point where X-ray spectrum #{} is acquired.".format(n_tot_sp_collected)
            )
            prompt.run()
    
            if prompt.execution_stopped:
                print("Execution stopped by the user.")
                return False, None, None
    
            if prompt.ok_pressed:
                # Error handling if get_frame_width or im_width are not set or invalid
                try:
                    frame_width_mm = EM_driver.get_frame_width()
                    if not hasattr(self, 'im_width') or self.im_width == 0:
                        raise AttributeError("im_width attribute missing or zero.")
                    self.pixel_size_um = frame_width_mm / self.im_width * 1e3  # um
                except Exception as e:
                    print("Error determining pixel size: {}".format(e))
                    return False, None, None
    
                spots_xy_list = EM_driver.frame_pixel_to_rel_coords(
                    (int(self.im_width / 2), int(self.im_height / 2)),
                    self.im_width,
                    self.im_height
                ) # In manual mode, assume center of image
                particle_cntr = None  # Not applicable for manual mode
                return True, spots_xy_list, particle_cntr
            
        elif self.sample_cfg.type == cnst.S_BULK_SAMPLE_TYPE:
        
            # Center stage onto next acquisition spot
            movement_success = self.go_to_next_frame()
            if not movement_success:
                # Recalculate shifted grid of points
                recalc_success = self._calc_bulk_grid_acquisition_spots()
                if recalc_success:
                    movement_success = self.go_to_next_frame()
                    if not movement_success:
                        print("Error moving to next frame")
                        return False, None, None
                else:
                    return False, None, None 
                
            # Error handling if get_frame_width or im_width are not set or invalid
            try:
                frame_width_mm = EM_driver.get_frame_width()
                if not hasattr(self, 'im_width') or self.im_width == 0:
                    raise AttributeError("im_width attribute missing or zero.")
                self.pixel_size_um = frame_width_mm / self.im_width * 1e3  # um
            except Exception as e:
                print("Error determining pixel size: {}".format(e))
                return False, None, None

            spots_xy_list = EM_driver.frame_pixel_to_rel_coords(
                (int(self.im_width / 2), int(self.im_height / 2)),
                self.im_width,
                self.im_height
            )  # In bulk mode, only measures at center of image
            particle_cntr = None  # Not applicable for bulk mode
            return True, spots_xy_list, particle_cntr
        
        elif self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE:
            # Move to the next detected particle and get acquisition coordinates
            was_particle_found = self.particle_finder.go_to_next_particle()
            if not was_particle_found:
                if self.verbose:
                    print('No more particles could be found on the sample.')
                return False, None, None
    
            particle_cntr = self.particle_finder.tot_par_cntr
            try:
                spots_xy_list = self.particle_finder.get_XS_acquisition_spots_coord_list(
                    n_tot_sp_collected
                )
            except Exception as e:
                print("Error getting acquisition spot coordinates: {}".format(e))
                return False, None, None
    
            return True, spots_xy_list, particle_cntr
    
        print("Acquisition mode not implemented for sample type: {}".format(self.sample_cfg.type))
        return False, None, None
    
    
    def convert_XS_coords_to_pixels(self, xy_coords):
        """
        Convert XY coordinates from the XS coordinate system to pixel coordinates.
    
        This function uses the EM_driver to transform coordinates relative to the frame 
        into pixel positions based on the image width and height. The returned pixel 
        coordinates are integer values.
    
        Parameters
        ----------
        xy_coords : tuple(float, float)
            The XY coordinates in the XS coordinate system to be converted.
    
        Returns
        -------
        tuple(int, int)
            The corresponding (x, y) pixel coordinates as integers.
        """
        xy_coords_pixels = EM_driver.frame_rel_to_pixel_coords(
            xy_coords,
            self.im_width,
            self.im_height
        ).astype(int)[0]
    
        return xy_coords_pixels

            
    def acquire_XS_spot_spectrum(
        self,
        x: float,
        y: float,
        max_acquisition_time: float,
        target_acquisition_counts: int
    ) -> tuple:
        """
        Acquire an X-ray spectrum (EDS/WDS) at the specified (x, y) position.
    
        Parameters
        ----------
        x, y : float
            X, Y coordinates for spectrum acquisition (in relative image units, as required by EM_driver).
            
            Coordinate System
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
        max_acquisition_time : float
            Maximum allowed acquisition time in seconds.
        target_acquisition_counts : int
            Target total X-ray counts for the spectrum.
    
        Returns
        -------
        spectrum_data : Any
            The acquired spectrum data (format as returned by EM_driver).
        background_data : Any
            The measured background data (format as returned by EM_driver).
        real_time : float
            Real elapsed time for the acquisition (seconds).
        live_time : float
            Detector live time for the acquisition (seconds).
    
        Raises
        ------
        EMError
            If the spectrum acquisition fails.
        """
        try:
            spectrum_data, background_data, real_time, live_time = EM_driver.acquire_XS_spectral_data(
                self.analyzer, x, y, max_acquisition_time, target_acquisition_counts
            )
            return spectrum_data, background_data, real_time, live_time
        except Exception as e:
            raise EMError(f"Failed to acquire X-ray spectrum at ({x}, {y}): {e}") from e
    
    #%% Microscope Control Functions
    # =============================================================================
    def adjust_BCF(self) -> float:
        """
        Adjust brightness, contrast, and focus (BCF).
    
        If automatic brightness and contrast adjustment is enabled (self.microscope_cfg.is_auto_BC),
        calls the automatic adjustment method. Otherwise, sets brightness and contrast
        to fixed values and then performs autofocus.
    
        Returns
        -------
        None
            
        Updates
        -------
        self._last_EM_adjustment_time
        
        Raises
        ------
        EMError
            If an error occurs during brightness, contrast, or focus adjustment.
        """
        try:
            if self.microscope_cfg.is_auto_BC:
                # Automatically adjust brightness, contrast, and focus
                adj_time = self._autoadjust_BCF()
            else:
                # Set brightness and contrast to fixed values, then autofocus
                self._set_frame_BC()
                adj_time = self._auto_focus()
                time.wait(0.5)
            self._last_EM_adjustment_time = adj_time
        except Exception as e:
            # You may want to log e here with more context if needed.
            raise EMError(
                f"Failed to adjust brightness, contrast, and focus: {e}"
            ) from e
            
            
    def _set_frame_BC(self) -> None:
        """
        Set EM frame brightness and contrast to the values stored in 
        self.microscope_cfg.brightness and self.microscope_cfg.contrast.
    
        Raises
        ------
        EMError
            If the EM driver fails to set brightness or contrast.
    
        Notes
        -----
        - Assumes self.microscope_cfg.brightness and self.microscope_cfg.contrast are valid values. For maximum robustness,
            consider validating that self.microscope_cfg.brightness and self.microscope_cfg.contrast are within the allowed
            range for the hardware.
        """
        try:
            # Set brightness and contrast to fixed values to ensure C tape is below threshold
            EM_driver.set_brightness(self.microscope_cfg.brightness)
            EM_driver.set_contrastss(self.microscope_cfg.contrast)
        except Exception as e:
            # Wrap any driver exception in a custom error for clarity
            raise EMError(f"Failed to set brightness/contrast: {e}") from e
        

    def _autoadjust_BCF(self) -> float:
        """
        Automatically adjust brightness, contrast, and focus (BCF).
    
        This method:
            1. Calls the EM driver's automatic brightness and contrast adjustment routine.
            2. Calls autofocus function.
            3. Returns the current timestamp after adjustments.
    
        Returns
        -------
        float
            The timestamp (seconds since epoch) when the adjustment was completed.
    
        Raises
        ------
        EMError
            If an error occurs during brightness/contrast or focus adjustment.
        """
        try:
            # Automated adjustment of brightness and contrast
            EM_driver.auto_contrast_brightness()
            # Automated focus adjustment
            self._auto_focus()
            return time.time()
        except Exception as e:
            raise EMError(
                f"Failed to auto-adjust brightness/contrast/focus: {e}"
            ) from e
    
    
    def _auto_focus(self) -> float:
        """
        Automatically adjust focus, ensuring the working distance (WD) is within allowed limits.
    
        This method:
            1. Calls the EM driver's autofocus method, which returns the current working distance (WD) in millimeters.
            3. If the WD is out of allowed bounds (`self._min_wd`, `self._max_wd`), clips it and readjusts focus.
            4. Returns the current timestamp after adjustment.
    
        Returns
        -------
        float
            The timestamp (seconds since epoch) when the autofocus adjustment was completed.
    
        Raises
        ------
        EMError
            If an error occurs during autofocus or WD retrieval/adjustment.
        """
        try:
            # Perform autofocus and get current working distance
            wd = EM_driver.auto_focus()
    
            # If WD is out of allowed bounds, clip and readjust
            if not (self._min_wd < wd < self._max_wd):
                print(f"Working distance of {wd:.1f} mm obtained through autofocus was out of accepted limits.")
                wd = float(np.clip(wd, self._min_wd, self._max_wd))
                print(f"WD was set to {wd:.1f} mm")
                EM_driver.adjust_focus(wd)
    
            return time.time()
        except Exception as e:
            raise EMError(f"Failed to auto-focus EM: {e}") from e
    
    
    def set_frame_width(self, frame_width):
        """
        Set the frame width at the microscope and update the pixel size accordingly.
    
        Parameters
        ----------
        frame_width : float
            The desired frame width in millimeters.
    
        Updates
        -------
        self.pixel_size_um : float
            The image pixel size in micrometers, calculated as (frame_width / image width in pixels) * 1000.

        """
        try:
            # Set frame width at EM
            EM_driver.set_frame_width(frame_width)
            # Update pixel size (in um)
            self.pixel_size_um = frame_width / self.im_width * 1e3  # um
        except Exception as e:
            raise EMError(f"Failed to set frame width: {e}") from e
        
        
    def move_to_pos(self, pos):
        """
        Move the EM stage to the specified (x, y) position and update the current position.
    
        Parameters
        ----------
        pos : tuple of float
            Target (x, y) coordinates.
    
        Updates
        -------
        self._current_pos : tuple of float
            The current position is set to the specified (x, y) coordinates.
        """
        try:
            # Extract the x y coordinates
            x, y = pos
            # Move EM
            EM_driver.move_to(x, y)
            # Update current position
            self._current_pos = (x, y)
        except Exception as e:
            raise EMError(f"Failed to move to desired position: {e}") from e
        
        
    def convert_pixel_pos_to_mm(self, pos_pixels):
        """
        Convert a position from pixel coordinates in the current image to absolute stage coordinates in millimeters.
    
        Parameters
        ----------
        pos_pixels : array-like of float
            Position in pixel coordinates (x, y).
    
        Returns
        -------
        pos_abs_mm : ndarray of float
            Absolute position in millimeters (x, y), suitable for EM stage movement.
        """
        # Calculate the center of the image in pixels
        center_pixels = np.array([self.im_width, self.im_height]) / 2
    
        # Compute the shift from the image center in pixels
        shift_pixels = pos_pixels - center_pixels
    
        # Convert the shift to micrometers (um)
        shift_um = self.pixel_size_um * shift_pixels
    
        # Transform coordinates to match stage reference system
        shift_um_stage_coords = shift_um * EM_driver.image_to_stage_coords_transform
    
        # Calculate absolute position in mm for the EM
        pos_abs_mm = self._current_pos + shift_um_stage_coords * 1e-3
    
        return pos_abs_mm
    
    @staticmethod
    def standby():
        """
        Put microscope in standby mode
        """
        EM_driver.standby()
    
    #%% Frame and Particle Navigation Methods
    # =============================================================================
    def go_to_next_frame(self):
        '''
        Moves the microscope to the next frame position in the current sample.
    
        This function checks if there are remaining frames to analyze. If so, it moves the microscope 
        stage to the next frame center, adjusts the frame width if necessary, and (optionally) 
        updates EM brightness, contrast and focus.
        It also prints the current frame information if verbose mode is enabled.
    
        Returns
        -------
        bool
            True if moved to the next frame, False if no frames remain to be analysed.
        '''
        # Check if all frames have been analysed in the current sample.
        if self._frame_cntr >= self.num_frames:
            # Return False if there are no more frames available
            return False
    
        # Move to frame
        self.move_to_pos(self.frame_pos_mm[self._frame_cntr])
        self.current_frame_label = self.frame_labels[self._frame_cntr]
    
        # Set frame width
        if self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE:
            min_fw, max_fw = EM_driver.get_range_frame_width()
            # Check that microscope still accepts the current frame width
            self.grid_search_fw_mm = np.clip(self.grid_search_fw_mm, min_fw, max_fw)
            self.set_frame_width(self.grid_search_fw_mm)
    
        # Adjust EM settings (focus, contrast, brightness) if too long has passed since last adjustments
        if time.time() - self._last_EM_adjustment_time > self.refresh_time:
            self.adjust_BCF()
        
    
        if self.verbose:
            print_single_separator()
            print(f"Moved to frame {self.current_frame_label} (#{self._frame_cntr + 1}/{self.num_frames}).")
    
        # Update counter of particles in current frame
        self._frame_cntr += 1
    
        return True
    
    
    def save_frame_image(self, filename, im_annotations=None, frame_image = None, save_dir=None):
        """
        Save an annotated and raw electron microscopy (EM) frame as a multi-page TIFF.
        
        This function retrieves a raw grayscale EM image, generates an annotated
        RGB version with optional markers and a scale bar, and saves both images
        into a single multi-page TIFF file. The annotated image is stored as the
        first page, and the raw image as the second page.
        
        Parameters
        ----------
        filename : str
            Name used for saved .tif image file
            
        im_annotations : dict | list(dict) | None, optional
            A dictionary with the following keys:
                - 'text' : tuple(text, xy_coords)
                    text: str, text to print
                    xy_coords : (int, int), coordinates in pixels where to place the text
                - 'circle' : tuple(radius, xy_center, is_filled, border_thickness)
                    radius : int, radius of circle, in pixels
                    xy_center : (int, int), (x,y) coordinates of center, in pixels
                    thickness : int, thickness of circle border. Set to -1 for filling
            For multiple annotations, use a list of dictionaries.
            Do not include the relative key if that particular annotation is not desired.
            
        frame_image : np.array | None, opt
            Frame image to save. Captures the current frame image if not provided
        
        save_dir : str, optional
            Directory in which to save the TIFF file. Defaults to self.results_dir
            if available, otherwise it needs to be provided.
        
        Notes
        -----
        - Images are saved as RGB to maximize compatibility across platforms.

        """
        # Determine save directory
        if not save_dir:
            if self.results_dir:
                save_dir = self.results_dir
            else:
                warnings.warn(
                    "No directory specified for saving frame image.",
                    UserWarning
                )
                return

        if not isinstance(frame_image, np.ndarray):
            # Get raw grayscale image from EM (H, W), dtype: uint8 or uint16
            frame_image = EM_driver.get_image_data(self.im_width, self.im_height, 1)
    
        # Convert grayscale to RGB for annotation
        if len(frame_image.shape) == 2 or frame_image.shape[2] == 1:
            # Grayscale image
            color_image = cv2.cvtColor(frame_image, cv2.COLOR_GRAY2RGB)
        else:
            # Already color
            color_image = frame_image
    
        # Draw annotations if provided
        if im_annotations is not None:
            if isinstance(im_annotations, dict):
                im_annotations = list(im_annotations)
                
            for ann_dict in im_annotations:
                
                # Add dot/circles
                if self.an_circle_key in ann_dict.keys():
                    radius, xy_center, border_thickness = ann_dict[self.an_circle_key]
                
                    # Draw filled red circle
                    cv2.circle(color_image, tuple(xy_center), radius, (255, 0, 0), border_thickness)  # RGB red
                
                # Add label text
                if self.an_text_key in ann_dict.keys():
                    text, text_xy = ann_dict[self.an_text_key]
                    cv2.putText(
                        color_image,
                        text,
                        text_xy,
                        cv2.FONT_HERSHEY_SIMPLEX,
                        1,
                        (255, 0, 0),  # RGB red
                        2,
                        cv2.LINE_AA
                    )
        
        # Add scale bar
        color_image = draw_scalebar(color_image, self.pixel_size_um)
    
        # Prepare save path
        save_path = os.path.join(save_dir, f"{filename}.tif")
    
        # Ensure dtype consistency (convert to uint8 if needed)
        if frame_image.dtype != np.uint8:
            frame_image = (frame_image / frame_image.max() * 255).astype(np.uint8)
            color_image = (color_image / color_image.max() * 255).astype(np.uint8)
    
        # Convert grayscale to RGB for saving
        if frame_image.ndim == 2:
            frame_image = cv2.cvtColor(frame_image, cv2.COLOR_GRAY2RGB)
            
        image_description_d = {
            "sample_ID" : self.sample_cfg.ID,
            "microscope_ID" : self.microscope_cfg.ID,
            "microscope_type" : self.microscope_cfg.type,
            "detector" : self.microscope_cfg.detector_type,
            "resolution" : (self.im_width, self.im_height),
            "pixel_size_um": self.pixel_size_um # um
            } 

        # Convert dict to JSON string, ASCII-safe for Preview
        desc_str = json.dumps(image_description_d, ensure_ascii=True)
        
        # Convert numpy arrays to Pillow Image objects, force RGB mode
        im1 = Image.fromarray(color_image.astype('uint8'), mode='RGB')
        im2 = Image.fromarray(frame_image.astype('uint8'), mode='RGB')
        
        # Save as multi-page TIFF, with description on first page
        im1.save(
            save_path,
            format='TIFF',
            description=desc_str,
            save_all=True,
            append_images=[im2],
            compression=None
        )
        
    
#%% Electron Microscope Sample Finder class    
class EM_Sample_Finder:
    """
    Class for locating and managing samples in an electron microscope (EM).

    This class provides methods for detecting sample features (such as the center and size of a C tape)
    using the microscope's navigation camera.
    
    Attributes
    ----------
    microscope_ID : str
        Identifier for the target microscope.
    center_pos : np.ndarray
        Initial center position of sample as a NumPy array of shape (2,), representing [x, y] coordinates.
    results_dir : Optional[str]
        Directory path to save results, or None if not set.
    verbose : bool
        If True, enables detailed debug output.
    development_mode (bool): Whether class is being used for testing image processing functions, without real-time acquisition.
        Enables use of class outside microscope enviroment
    
    Example
    -------
    >>> import numpy as np
    >>> sample_center = np.array([0.0, 0.0])  # stage coordinates in mm
    >>> finder = EM_Sample_Finder(
    ...     microscope_ID='MySEM',
    ...     center_pos=sample_center,
    ...     sample_half_width_mm: 3,
    ...     substrate_width_mm: 12,
    ...     results_dir='./results',
    ...     verbose=True
    ... )
    >>> ctape_result = finder.detect_Ctape()
    Detecting position of C tape
    C tape detected.
    >>> if ctape_result is not None:
    ...     center_pos, sample_hw_mm = ctape_result
    ...     print("C tape center (mm):", center_pos)
    ...     print("C tape half-width (mm):", sample_hw_mm)
    ... else:
    ...     print("C tape detection failed.")
     
    Notes
    -----
    - The navigation camera image can be provided directly (for offline testing) or acquired live from the microscope.
    - For successful detection, the microscope calibration file must include the attributes:
      'navcam_im_w_mm', 'navcam_x_offset', 'navcam_y_offset'.
    - The detected center and half-width (radius) can be used to configure automated acquisition grids or for quality control.
    """

    def __init__(
        self,
        microscope_ID: str,
        center_pos: Union[np.ndarray, tuple, list],
        sample_half_width_mm: float,
        substrate_width_mm: float,
        development_mode: Optional[bool] = False, 
        results_dir: Optional[str] = None,
        verbose: bool = False
    ):
        """
        Initialize the EM_Sample_Finder object.

        Parameters:
            microscope_ID (str): Identifier for the target microscope. Used to load the appropriate instrument driver.
            center_pos (np.ndarray): Initial center position as a NumPy array of shape (2,), representing [x, y] coordinates.
            sample_half_width_mm (float): Half-width of sample
            results_dir (Optional[str], optional): Directory path to save results. If None, results are not saved to disk. Defaults to None.
            verbose (bool, optional): If True, enables detailed debug output. Defaults to False.
            development_mode (bool): Whether class is being used for testing image processing functions, without real-time acquisition.
                Enables use of class outside microscope enviroment
        """
        self.microscope_ID = microscope_ID
        # Load microscope driver for instrument microscope_ID
        EM_driver.load_microscope_driver(microscope_ID)
        if not development_mode:
            if not EM_driver.is_at_EM:
                raise EMError("Instrument driver could not be loaded")
        
        self._sample_half_width_mm = sample_half_width_mm
        self._substrate_width_mm = substrate_width_mm
        self._center_pos = np.array(center_pos, dtype=float)
        self.results_dir = results_dir
        self.development_mode = development_mode
        self.verbose = verbose


    def detect_Ctape(
        self, navcam_im: Optional[np.ndarray] = None,
    ) -> Optional[Tuple[np.ndarray, float]]:
        """
        Detects the effective center position and radius of the C tape using the navigation camera image.
        
        Uses navcam_im if provided (e.g., when testing function), or acquires it directly at microscope
        via EM_driver.get_navigation_camera_image()
        
        Returns
        -------
        (center_pos, sample_hw_mm) if successful, else None.
        """
        if self.verbose:
            print_single_separator()
            print('Detecting position of C tape...')
        
        # Collect NavCam image
        if navcam_im is None:
            navcam_im = EM_driver.get_navigation_camera_image()
            
        if navcam_im is None or not hasattr(navcam_im, "shape"):
            print("No valid navigation camera image provided or acquired. C-tape detection skipped")
            return None

        # Get size of image in pixels
        navcam_h, navcam_w, _ = navcam_im.shape
        
        # Load navigation camera calibrated parameters and check their presence
        required_attrs = ['navcam_im_w_mm', 'navcam_x_offset', 'navcam_y_offset']
        for attr in required_attrs:
            if not hasattr(EM_driver, attr):
                raise AttributeError(f"Microscope calibration file at {EM_driver.microscope_calib_dir} is missing required attribute '{attr}'")

        # Calculate pixel size
        navcam_pixel_size = EM_driver.navcam_im_w_mm / navcam_w

        # Calculate position of center of stub within navcam_im (in pixels)
        stub_c = ((self._center_pos / navcam_pixel_size + np.array([navcam_w, -navcam_h]) / 2) * np.array([1, -1])).astype(np.uint16)
        stub_c += np.array([EM_driver.navcam_x_offset, EM_driver.navcam_y_offset]).astype(np.uint16)

        # Calculate size of stub half-width in pixels
        stub_hw = int(self._substrate_width_mm / navcam_pixel_size / 2)

        # Crop image around stub
        y1 = max(stub_c[1] - stub_hw, 0)
        y2 = min(stub_c[1] + stub_hw + 1, navcam_h)
        x1 = max(stub_c[0] - stub_hw, 0)
        x2 = min(stub_c[0] + stub_hw + 1, navcam_w)
        stub_im = navcam_im[y1:y2, x1:x2]

        # Split RGB channels and detect edges
        # Normalize image to 8-bit and adjust brightness
        stub_im_normalized = self.normalise_img(stub_im)
        # if self.development_mode: cv2.imshow('Normalised', stub_im_normalized)
        channels = cv2.split(stub_im_normalized)
        edges_channels = [
            cv2.Canny(cv2.GaussianBlur(channel, (9, 9), 1.5), 100, 200)
            for channel in channels
        ]
        edges_combined = cv2.merge(edges_channels)
        if self.development_mode: cv2.imshow('Combined Edges (RGB)', edges_combined)
        gray = cv2.cvtColor(edges_combined, cv2.COLOR_BGR2GRAY)
        # if self.development_mode: cv2.imshow('Gray', gray)
        gray[gray < 100] = 0
        if self.development_mode:
            cv2.imshow('Thresholded Gray', gray)
            # cv2.imwrite(os.path.join(self.results_dir,'Thresholded Gray.png'), gray)
        
        
        # Detect circles in image
        min_radius = int(self._sample_half_width_mm / navcam_pixel_size / 1.5)
        max_radius = int(stub_hw * 0.9)
        circles = cv2.HoughCircles(
            gray, cv2.HOUGH_GRADIENT, dp=1.4, minDist=5,
            param1=50, param2=50, minRadius=min_radius, maxRadius=max_radius
        )
        
        if self.development_mode:
            output = stub_im.copy()
            if circles is not None:
                circles = np.uint16(np.around(circles))
                for (x, y, r) in circles[0, :]:
                    cv2.circle(output, (x, y), r, (0, 255, 0), 2)
                    cv2.circle(output, (x, y), 2, (0, 0, 255), 3)
            cv2.imshow('Detected Circles', output)
            # cv2.imwrite(os.path.join(self.results_dir,'Detected Circles.png'), output)

        
        # Average pixel intensity function
        def average_pixel_intensity(image: np.ndarray, center: tuple, radius: float) -> float:
            # Convert color to grayscale if needed
            if len(image.shape) == 3 and image.shape[2] == 3:
                gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
            else:
                gray_image = image
        
            x, y = int(center[0]), int(center[1])
            mask = np.zeros(gray_image.shape, dtype=np.uint8)
            cv2.circle(mask, (x, y), int(radius), 255, -1)  # white filled circle on black mask
        
            masked_pixels = gray_image[mask == 255]
            if masked_pixels.size == 0:
                return 0.0
            return float(np.mean(masked_pixels))

        stub_im_original = stub_im_normalized.copy()  # Don't modify for drawing

        # Filter circles by intensity
        filtered_circles = []
        if circles is not None:
            circles = np.round(circles[0, :]).astype("int")
            for (x, y, r) in circles:
                avg_intensity = average_pixel_intensity(stub_im_original, (x, y), r)
                if avg_intensity <= 100:
                    filtered_circles.append((x, y, r))

        if self.development_mode:
            filtered_debug = stub_im_original.copy()
            for (x, y, r) in filtered_circles:
                cv2.circle(filtered_debug, (x, y), r, (255, 0, 0), 2)  # Blue = accepted
            cv2.imshow('Filtered Circles', filtered_debug)

        # Find best circle
        x, y, r = None, None, None
        if len(filtered_circles) > 1:
            avg_x = np.mean([circle[0] for circle in filtered_circles])
            avg_y = np.mean([circle[1] for circle in filtered_circles])
            avg_center = (int(avg_x), int(avg_y))
            distances = [circle[2] - np.sqrt((circle[0] - avg_x) ** 2 + (circle[1] - avg_y) ** 2) for circle in filtered_circles]
            intersection_radius = int(min(distances))
            x, y = avg_center
            r = max(intersection_radius, min_radius)
        elif len(filtered_circles) > 0:
            x, y, r = filtered_circles[0]

        # Compute output
        sample_hw_mm_mm = None
        if len(filtered_circles) > 0 and x is not None and y is not None and r is not None:
            center_pos_eff = self._center_pos + (np.array([x, y]) - stub_hw) * navcam_pixel_size * np.array([1, -1])
            sample_hw_mm_mm = r * navcam_pixel_size * 0.9  # Add margin for safety
            Ctape_coords = (center_pos_eff, sample_hw_mm_mm)
            if self.verbose:
                print('C tape detected.')
        else:
            x, y = stub_hw, stub_hw
            r = int(self._sample_half_width_mm / navcam_pixel_size)
            Ctape_coords = None
            if self.verbose:
                print(f'The C tape could not be automatically detected. Using {tuple(self._center_pos)} instead.')

        # Draw detected circle on image
        cv2.circle(stub_im, (x, y), r, (0, 255, 0), 1)
        cv2.rectangle(stub_im, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)
        if self.development_mode: cv2.imshow('Detected region', stub_im)
        # Save result image if results_dir is set
        if self.results_dir:
            filename = os.path.join(self.results_dir, cnst.NAVCAM_IM_FILENAME + '.png')
            cv2.imwrite(filename, stub_im)

        return Ctape_coords
    
    
    def normalise_img(self, img: np.ndarray, target_brightness: float = 128.0) -> np.ndarray:
        """
        Normalize brightness of an RGB image to a target brightness level.
    
        Parameters:
            img (np.ndarray): Input RGB image (uint8).
            target_brightness (float): Desired average brightness (0–255).
    
        Returns:
            np.ndarray: Brightness-normalized RGB image (uint8).
        """
        if not isinstance(img, np.ndarray):
            raise TypeError(f"Expected image as np.ndarray, got {type(img)}")
    
        # Handle grayscale or 1-channel image
        if img.ndim == 2:
            img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
        elif img.ndim == 3:
            if img.shape[2] == 1:
                img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
            elif img.shape[2] == 4:
                img = img[:, :, :3]  # strip alpha channel
    
        # Final check
        if img.ndim != 3 or img.shape[2] != 3:
            raise ValueError("Input image must be an RGB image with 3 channels.")
    
        # Convert to grayscale to compute current brightness
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        current_brightness = np.mean(gray)
    
        # Avoid division by zero
        if current_brightness == 0:
            scale = 1.0
        else:
            scale = target_brightness / current_brightness
    
        # Scale image and clip to valid range
        img_float = img.astype(np.float32) * scale
        img_scaled = np.clip(img_float, 0, 255).astype(np.uint8)
    
        return img_scaled