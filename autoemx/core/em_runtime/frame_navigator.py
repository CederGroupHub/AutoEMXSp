#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Frame and grid navigation module for automated electron microscopy.

This module provides functionality for calculating scan frames, grid layouts,
and navigating between acquisition points in a sample.

Classes
-------
FrameNavigator
    Handles frame center calculation and frame-by-frame navigation.

Created on 2026
@author: Andrea
"""
from typing import List, Optional, Sequence, Tuple

import time
import numpy as np

from autoemx.core.em_runtime.xsp_spot_selection import XSpSpotSelectorCallback

import autoemx.utils.constants as cnst
from autoemx.utils.helper import AlphabetMapper, Prompt_User, print_single_separator
from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder

from autoemx._logging import get_logger
logger = get_logger(__name__)

# (frame_label, center_xy_mm, frame_width_mm)
FrameDescriptor = Tuple[str, Tuple[float, float], float]


class FrameNavigator:
    """
    Handles frame planning and navigation for electron microscopy samples.
    
    Calculates frame centers for systematic scanning across samples,
    supports multiple sample types (powder, bulk, manual), and manages
    frame-by-frame navigation during acquisition.
    
    Parameters
    ----------
    EM_controller : EM_Controller
        Reference to the parent EM_Controller for hardware access.
    sample_cfg : SampleConfig
        Sample configuration.
    sample_substrate_cfg : SampleSubstrateConfig
        Sample substrate configuration.
    measurement_cfg : MeasurementConfig
        Measurement configuration.
    powder_meas_cfg : PowderMeasurementConfig
        Powder measurement configuration.
    bulk_meas_cfg : BulkMeasurementConfig
        Bulk measurement configuration.
    center_pos : tuple
        Sample center position (x, y).
    sample_hw_mm : float
        Sample half-width in mm.
    im_width : int
        Image width in pixels.
    im_height : int
        Image height in pixels.
    EM_driver_obj : module
        EM driver for hardware access.
    results_dir : str or None
        Directory for saving results.
    verbose : bool
        If True, print progress information.
    development_mode : bool
        If True, enables offline testing mode.
    """
    
    def __init__(self, EM_controller, sample_cfg, sample_substrate_cfg, measurement_cfg,
                 powder_meas_cfg, bulk_meas_cfg, center_pos, sample_hw_mm,
                 im_width, im_height, EM_driver_obj, results_dir=None,
                 verbose=True, development_mode=False,
                 xsp_spot_selector: Optional[XSpSpotSelectorCallback] = None):
        self.EM_controller = EM_controller
        self.sample_cfg = sample_cfg
        self.sample_substrate_cfg = sample_substrate_cfg
        self.measurement_cfg = measurement_cfg
        self.powder_meas_cfg = powder_meas_cfg
        self.bulk_meas_cfg = bulk_meas_cfg
        self._center_pos = center_pos
        self._sample_hw_mm = sample_hw_mm
        self.im_width = im_width
        self.im_height = im_height
        self.EM_driver = EM_driver_obj
        self.results_dir = results_dir
        self.verbose = verbose
        self.development_mode = development_mode
        self.xsp_spot_selector = xsp_spot_selector
        
        # Frame tracking
        self._frame_cntr = 0
        self._bulk_offset_cntr = 0
        self.frame_pos_mm = []
        self.frame_labels = []
        self.num_frames = 0
        self.current_frame_label = ''
        self.particle_finder = None
        # Current frame center position in stage coordinates (mm).
        # Used by EM_Controller.convert_pixel_pos_to_mm.
        self._current_pos = center_pos
        self.grid_search_fw_mm = None
        self._band_frame_descriptors: List[FrameDescriptor] = []
        self._visited_band_frames: List[FrameDescriptor] = []
    
    
    def initialise_sample_navigator(self, grid_search_fw_mm, exclude_sample_margin=True):
        """
        Initialize the sample navigator for automated spectra collection.
        
        Supports 'powder' and 'bulk' samples with automated navigation,
        as well as manual navigation mode.
        
        Parameters
        ----------
        grid_search_fw_mm : float
            Frame width in mm for grid search.
        exclude_sample_margin : bool
            Whether to exclude sample margin when calculating frames.
            
        Raises
        ------
        NotImplementedError
            If the sample type is not supported.
        """
        self.grid_search_fw_mm = grid_search_fw_mm
        
        if self.sample_cfg.is_particle_acquisition:
            # Set frame width, and update current pixel size
            if self.EM_driver.is_microscope_connected():
                min_fw, max_fw = self.EM_driver.get_range_frame_width()
                self.grid_search_fw_mm = np.clip(
                    self.powder_meas_cfg.par_search_frame_width_um / 1000, 
                    min_fw, 
                    max_fw
                )
            self.EM_controller.set_frame_width(self.grid_search_fw_mm)
            
            # Calculate frame centers for particle search
            im_h_to_w_ratio = self.im_height / self.im_width
            self._calc_frame_centers(
                horizontal_spacing_mm=self.grid_search_fw_mm,
                im_h_to_w_ratio=im_h_to_w_ratio,
                center_pos=self._center_pos,
                randomize_frames=True,
                exclude_sample_margin=True
            )
            
            # Initialise particle finder for powder samples
            self.particle_finder = EM_Particle_Finder(
                self.EM_controller,
                powder_meas_cfg=self.powder_meas_cfg,
                xsp_spot_selector=self.xsp_spot_selector,
                results_dir=self.results_dir,
                verbose=self.verbose,
                development_mode=self.development_mode
            )
        
        elif self.sample_cfg.is_grid_acquisition:
            if self.EM_driver.is_microscope_connected():
                min_fw, max_fw = self.EM_driver.get_range_frame_width()
                self.grid_search_fw_mm = np.clip(
                    self.bulk_meas_cfg.image_frame_width_um / 1000, 
                    min_fw, 
                    max_fw
                )
            self.EM_controller.microscope_ctrl.set_frame_width(self.grid_search_fw_mm)
            # Construct grid of acquisition spots
            self._calc_bulk_grid_acquisition_spots()
        
        elif self.measurement_cfg.is_manual_navigation:
            self.frame_pos_mm = None
            self.frame_labels = None
            self.num_frames = np.inf
        
        else:
            raise NotImplementedError(
                f"Sample type '{self.sample_cfg.type}' is not supported for automated composition analysis. "
                "Use measurement_cfg.is_manual_navigation = True."
            )
    
    
    def _calc_bulk_grid_acquisition_spots(self) -> bool:
        """
        Calculate and apply offset to center position for bulk grid acquisition.
        
        Constructs a square grid of acquisition spots with optional offset
        in three directions: (x, 0), (x, x), and (0, x).
        
        Returns
        -------
        bool
            True if grid was constructed, False if offset exceeds grid spacing.
        """
        # Calculate offset distance in micrometers
        offset_dist_um = (
            self.bulk_meas_cfg.min_xsp_spots_distance_um * 
            np.ceil(self._bulk_offset_cntr / 3)
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
    
    
    def _calc_frame_centers(self, horizontal_spacing_mm, im_h_to_w_ratio, center_pos,
                            randomize_frames, exclude_sample_margin):
        """
        Generate evenly spaced scanning locations (frames) within a sample area.
        
        Calculates a grid of (x, y) positions covering a circular or rectangular 
        region, optionally avoiding sample edges. Assigns each frame a unique label.
        
        Parameters
        ----------
        horizontal_spacing_mm : float
            Horizontal spacing between grid spots in mm.
        im_h_to_w_ratio : float
            Aspect ratio (height/width) of frame dimensions.
        center_pos : tuple
            (x, y) coordinates of sample center in stage coordinates.
        randomize_frames : bool
            If True, shuffle frame order to reduce spatial bias.
        exclude_sample_margin : bool
            If True, apply margin to avoid rough border and contamination.
        """
        cx, cy = center_pos
        
        # Determine the usable sample half width
        if exclude_sample_margin:
            margin = 2 * horizontal_spacing_mm * np.sqrt(1 + (im_h_to_w_ratio) ** 2)
            sample_hw_mm = self._sample_hw_mm - margin
        else:
            sample_hw_mm = self._sample_hw_mm
        
        # Define region checker function based on substrate shape
        if self.sample_substrate_cfg.shape == cnst.CIRCLE_SUBSTRATE_SHAPE:
            def is_inside_region(x, y):
                return (x - cx) ** 2 + (y - cy) ** 2 < sample_hw_mm ** 2
        elif self.sample_substrate_cfg.shape == cnst.SQUARE_SUBSTRATE_SHAPE:
            rect_left = cx - sample_hw_mm
            rect_right = cx + sample_hw_mm
            rect_top = cy - sample_hw_mm * im_h_to_w_ratio
            rect_bottom = cy + sample_hw_mm * im_h_to_w_ratio
            
            def is_inside_region(x, y):
                return (rect_left <= x <= rect_right) and (rect_top <= y <= rect_bottom)
        else:
            raise ValueError(
                f"Sample substrate shape must be one among {self.sample_substrate_cfg.ALLOWED_SHAPES}"
            )
        
        half_n_frames_x = int(self._sample_hw_mm / horizontal_spacing_mm) + 1
        half_n_frames_y = int(self._sample_hw_mm / (horizontal_spacing_mm * im_h_to_w_ratio)) + 1
        
        frames = []
        alphabet_mapper = AlphabetMapper()
        for i in range(-half_n_frames_x, half_n_frames_x + 1):
            label_letter = alphabet_mapper.get_letter(i + half_n_frames_x)
            for j in range(-half_n_frames_y, half_n_frames_y + 1):
                label = label_letter + str(j + half_n_frames_y)
                x = cx + i * horizontal_spacing_mm
                y = cy + j * horizontal_spacing_mm * im_h_to_w_ratio
                if is_inside_region(x, y):
                    # Compute distance and angle for spiral ordering
                    dist = np.hypot(x - cx, y - cy)
                    angle = np.arctan2(y - cy, x - cx)
                    frames.append(((x, y), label, dist, angle))

        # Spiral order: sort by distance from center, then by angle
        if frames:
            frames.sort(key=lambda tup: (tup[2], tup[3]))
            frame_centers, frame_labels = zip(*[(f[0], f[1]) for f in frames])
        else:
            frame_centers, frame_labels = [], []

        # If randomize_frames is requested, override spiral with random
        if randomize_frames and frame_centers:
            import random
            zipped = list(zip(frame_centers, frame_labels))
            random.shuffle(zipped)
            frame_centers, frame_labels = zip(*zipped)

        self.frame_pos_mm = list(frame_centers)
        self.frame_labels = list(frame_labels)
        self.num_frames = len(frame_centers)
        self._sync_band_frame_descriptors()

    def reset_frame_cursor(self) -> None:
        """Reset the frame index so the next ``go_to_next_frame`` starts at frame 0."""
        self._frame_cntr = 0
        self.current_frame_label = ''
        self._visited_band_frames = []

    def _apply_search_frame_width(self, frame_width_mm: float) -> float:
        """Clip and apply the particle-search FOV; return the effective width in mm."""
        self.grid_search_fw_mm = float(frame_width_mm)
        if self.EM_driver.is_microscope_connected():
            min_fw, max_fw = self.EM_driver.get_range_frame_width()
            self.grid_search_fw_mm = float(np.clip(self.grid_search_fw_mm, min_fw, max_fw))
            self.EM_controller.set_frame_width(self.grid_search_fw_mm)
        return float(self.grid_search_fw_mm)

    def _sync_band_frame_descriptors(self) -> None:
        """Refresh planned-frame metadata for the active size band."""
        fw = float(self.grid_search_fw_mm) if self.grid_search_fw_mm is not None else 0.0
        self._band_frame_descriptors = [
            (str(label), (float(pos[0]), float(pos[1])), fw)
            for pos, label in zip(self.frame_pos_mm, self.frame_labels)
        ]

    def get_visited_band_frames(self) -> List[FrameDescriptor]:
        """Return frames visited during the current size band (in visit order)."""
        return list(self._visited_band_frames)

    def configure_top_level_particle_frames(
        self,
        frame_width_mm: float,
        randomize_frames: bool = True,
        exclude_sample_margin: bool = True,
    ) -> None:
        """Rebuild the sample-wide particle-search grid for the coarsest size band."""
        effective_fw = self._apply_search_frame_width(frame_width_mm)
        im_h_to_w_ratio = self.im_height / self.im_width
        self._calc_frame_centers(
            horizontal_spacing_mm=effective_fw,
            im_h_to_w_ratio=im_h_to_w_ratio,
            center_pos=self._center_pos,
            randomize_frames=randomize_frames,
            exclude_sample_margin=exclude_sample_margin,
        )
        self.reset_frame_cursor()

    def configure_particle_subframes(
        self,
        parent_frames: Sequence[FrameDescriptor],
        child_fw_mm: float,
        randomize_frames: bool = True,
    ) -> None:
        """Build hierarchical child frames inside each parent FOV.

        Child labels use a compact lowercase AlphabetMapper scheme (``a0``, ``b1``, …)
        and are joined to the parent id as ``{parent}_{local}``.
        """
        if not parent_frames:
            self.frame_pos_mm = []
            self.frame_labels = []
            self.num_frames = 0
            self._band_frame_descriptors = []
            self.reset_frame_cursor()
            return

        effective_fw = self._apply_search_frame_width(child_fw_mm)
        im_h_to_w_ratio = self.im_height / self.im_width
        alphabet_mapper = AlphabetMapper(alphabet='abcdefghijklmnopqrstuvwxyz')
        frames = []
        sample_cx, sample_cy = self._center_pos

        for parent_label, (cx, cy), parent_fw_mm in parent_frames:
            parent_fw_mm = float(parent_fw_mm)
            n_x = max(1, int(round(parent_fw_mm / effective_fw)))
            n_y = max(1, int(round((parent_fw_mm * im_h_to_w_ratio) / (effective_fw * im_h_to_w_ratio))))
            half_x = (n_x - 1) / 2.0
            half_y = (n_y - 1) / 2.0
            for i in range(n_x):
                letter = alphabet_mapper.get_letter(i)
                for j in range(n_y):
                    local_label = f"{letter}{j}"
                    x = float(cx) + (i - half_x) * effective_fw
                    y = float(cy) + (j - half_y) * effective_fw * im_h_to_w_ratio
                    label = f"{parent_label}_{local_label}"
                    dist = np.hypot(x - sample_cx, y - sample_cy)
                    angle = np.arctan2(y - sample_cy, x - sample_cx)
                    frames.append(((x, y), label, dist, angle))

        if frames:
            frames.sort(key=lambda tup: (tup[2], tup[3]))
            frame_centers, frame_labels = zip(*[(f[0], f[1]) for f in frames])
            if randomize_frames:
                import random
                zipped = list(zip(frame_centers, frame_labels))
                random.shuffle(zipped)
                frame_centers, frame_labels = zip(*zipped)
            self.frame_pos_mm = list(frame_centers)
            self.frame_labels = list(frame_labels)
        else:
            self.frame_pos_mm = []
            self.frame_labels = []
        self.num_frames = len(self.frame_pos_mm)
        self._sync_band_frame_descriptors()
        self.reset_frame_cursor()

    def _record_visited_frame(self) -> None:
        """Append the current frame to the visited list for hierarchical sub-banding."""
        if self._current_pos is None or self.grid_search_fw_mm is None:
            return
        descriptor = (
            str(self.current_frame_label),
            (float(self._current_pos[0]), float(self._current_pos[1])),
            float(self.grid_search_fw_mm),
        )
        self._visited_band_frames.append(descriptor)
    
    
    def go_to_next_frame(self, microscope_ctrl) -> bool:
        """
        Move the microscope to the next frame position.
        
        Checks for remaining frames, moves to next position, and adjusts
        frame width and EM settings as needed.
        
        Parameters
        ----------
        microscope_ctrl : MicroscopeController
            Reference to microscope controller.
            
        Returns
        -------
        bool
            True if moved to next frame, False if no frames remain.
        """
        is_particle_stats_measurement = (
            self.measurement_cfg.type == self.measurement_cfg.PARTICLE_STATS_MEAS_TYPE_KEY
        )
        
        if is_particle_stats_measurement and self.measurement_cfg.is_manual_navigation:
            prompt = Prompt_User(
                title="Select next Frame",
                message=f"Go to next frame to analyze (#{self._frame_cntr})."
            )
            prompt.run()
            
            if prompt.execution_stopped:
                logger.warning("⚠️ Execution stopped by the user.")
                return False
            
            if prompt.ok_pressed:
                self.current_frame_label = self._frame_cntr
                frame_width = self.EM_driver.get_frame_width()
                self.grid_search_fw_mm = frame_width
                # In manual mode, keep the last known stage position as frame reference.
                # If no explicit frame move occurred yet, use configured sample center.
                if self._current_pos is None:
                    self._current_pos = self._center_pos
        
        else:
            # Check if all frames have been analysed
            if self._frame_cntr >= self.num_frames:
                return False
            
            # Move to frame
            next_frame_pos = self.frame_pos_mm[self._frame_cntr]
            microscope_ctrl.move_to_pos(next_frame_pos)
            self._current_pos = next_frame_pos
            self.current_frame_label = self.frame_labels[self._frame_cntr]
            
            # Set frame width for particle acquisition
            if self.sample_cfg.is_particle_acquisition:
                min_fw, max_fw = self.EM_driver.get_range_frame_width()
                self.grid_search_fw_mm = np.clip(self.grid_search_fw_mm, min_fw, max_fw)
                microscope_ctrl.set_frame_width(self.grid_search_fw_mm)
            
            # Adjust EM settings if too long has passed since last adjustments
            if time.time() - microscope_ctrl._last_EM_adjustment_time > microscope_ctrl.refresh_time:
                microscope_ctrl.adjust_BCF()
        
        if self.verbose:
            print_single_separator()
            logger.info(f"▶️ Moved to frame {self.current_frame_label} (#{self._frame_cntr + 1}/{self.num_frames}).")

        self._record_visited_frame()
        
        # Update frame counter
        self._frame_cntr += 1
        
        return True
