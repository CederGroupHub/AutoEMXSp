#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Electron Microscope Sample Finder module.

Provides utilities for locating and managing samples within the EM, particularly
for detecting carbon tape position and dimensions.

Classes
-------
EM_Sample_Finder
    Utilities for locating samples and detecting carbon tape.

Created on 2026
@author: Andrea
"""
import os
from typing import Optional, Union, Tuple

import numpy as np
import cv2

from autoemx.utils.helper import EMError, print_single_separator
from autoemx.core.em_runtime.image_utilities import normalise_img
import autoemx.utils.constants as cnst
from autoemx import microscope_drivers as EM_driver

from autoemx._logging import get_logger
logger = get_logger(__name__)


class EM_Sample_Finder:
    """
    Class for locating and managing samples in an electron microscope (EM).
    
    Provides methods for detecting sample features (such as the center and size of 
    a C tape) using the microscope's navigation camera.
    
    Attributes
    ----------
    microscope_ID : str
        Identifier for the target microscope.
    center_pos : np.ndarray
        Initial center position of sample as [x, y] coordinates in mm.
    results_dir : Optional[str]
        Directory path to save results, or None if not set.
    verbose : bool
        If True, enables detailed debug output.
    development_mode : bool
        If True, enables offline testing without real-time acquisition.
    
    Example
    -------
    >>> import numpy as np
    >>> sample_center = np.array([0.0, 0.0])
    >>> finder = EM_Sample_Finder(
    ...     microscope_ID='MySEM',
    ...     center_pos=sample_center,
    ...     sample_half_width_mm=3,
    ...     substrate_width_mm=12,
    ...     results_dir='./results',
    ...     verbose=True
    ... )
    >>> ctape_result = finder.detect_Ctape()
    >>> if ctape_result is not None:
    ...     center_pos, sample_hw_mm = ctape_result
    ...     print("C tape center (mm):", center_pos)
    ...     print("C tape half-width (mm):", sample_hw_mm)
    
    Notes
    -----
    - The navigation camera image can be provided directly (for offline testing) or 
      acquired live from the microscope.
    - For successful detection, the microscope calibration file must include:
      'navcam_im_w_mm', 'navcam_x_offset', 'navcam_y_offset'.
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
        
        Parameters
        ----------
        microscope_ID : str
            Identifier for the target microscope.
        center_pos : array-like
            Initial center position as [x, y] coordinates in mm.
        sample_half_width_mm : float
            Half-width of sample in mm.
        substrate_width_mm : float
            Width of substrate (e.g., carbon tape) in mm.
        development_mode : bool, optional
            If True, enables offline testing mode (default: False).
        results_dir : Optional[str], optional
            Directory path to save results (default: None).
        verbose : bool, optional
            If True, enables detailed debug output (default: False).
        """
        self.microscope_ID = microscope_ID
        # Load microscope driver
        EM_driver.load_microscope_driver(microscope_ID)
        self.EM_driver = EM_driver
        if not development_mode:
            try:
                connect_result = False
                if hasattr(self.EM_driver, "connect_to_microscope"):
                    connect_result = bool(self.EM_driver.connect_to_microscope(warn_if_unavailable=True))

                is_connected = connect_result
                if hasattr(self.EM_driver, "is_microscope_connected"):
                    is_connected = bool(self.EM_driver.is_microscope_connected()) or is_connected
            except Exception as e:
                raise EMError("Instrument driver could not be loaded") from e

            if not is_connected:
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
        Detect the effective center position and radius of the C tape.
        
        Uses the provided navigation camera image or acquires it from the microscope.
        Applies edge detection and circle fitting to locate the tape.
        
        Parameters
        ----------
        navcam_im : Optional[np.ndarray]
            Navigation camera image. If None, acquires from microscope.
        
        Returns
        -------
        Optional[Tuple[np.ndarray, float]]
            Tuple of (center_pos, sample_hw_mm) if detection successful, else None.
        """
        if self.verbose:
            print_single_separator()
            logger.info('🔬 Detecting position of C tape...')
        
        # Collect NavCam image
        if navcam_im is None:
            navcam_im = self.EM_driver.get_navigation_camera_image()
        
        if navcam_im is None or not hasattr(navcam_im, "shape"):
            logger.warning("⚠️ No valid navigation camera image provided or acquired. C-tape detection skipped")
            return None
        
        # Get size of image in pixels
        navcam_h, navcam_w, _ = navcam_im.shape
        
        # Load navigation camera calibrated parameters
        required_attrs = ['navcam_im_w_mm', 'navcam_x_offset', 'navcam_y_offset']
        for attr in required_attrs:
            if not hasattr(self.EM_driver, attr):
                raise AttributeError(
                    f"Microscope calibration file at {self.EM_driver.microscope_calib_dir} "
                    f"is missing required attribute '{attr}'"
                )
        
        # Calculate pixel size
        navcam_pixel_size = self.EM_driver.navcam_im_w_mm / navcam_w
        
        # Calculate position of center of stub within navcam_im (in pixels)
        stub_c = (
            (self._center_pos / navcam_pixel_size + 
             np.array([navcam_w, -navcam_h]) / 2) * 
            np.array([1, -1])
        ).astype(np.uint16)
        stub_c += np.array([
            self.EM_driver.navcam_x_offset, 
            self.EM_driver.navcam_y_offset
        ]).astype(np.uint16)
        
        # Calculate size of stub half-width in pixels
        stub_hw = int(self._substrate_width_mm / navcam_pixel_size / 2)
        
        # Crop image around stub
        y1 = max(stub_c[1] - stub_hw, 0)
        y2 = min(stub_c[1] + stub_hw + 1, navcam_h)
        x1 = max(stub_c[0] - stub_hw, 0)
        x2 = min(stub_c[0] + stub_hw + 1, navcam_w)
        stub_im = navcam_im[y1:y2, x1:x2]
        
        # Detect edges
        stub_im_normalized = normalise_img(stub_im)
        channels = cv2.split(stub_im_normalized)
        edges_channels = [
            cv2.Canny(cv2.GaussianBlur(channel, (9, 9), 1.5), 100, 200)
            for channel in channels
        ]
        edges_combined = cv2.merge(edges_channels)
        if self.development_mode:
            cv2.imshow('Combined Edges (RGB)', edges_combined)
        
        gray = cv2.cvtColor(edges_combined, cv2.COLOR_BGR2GRAY)
        gray[gray < 100] = 0
        if self.development_mode:
            cv2.imshow('Thresholded Gray', gray)
        
        # Detect circles in image
        min_radius = int(self._sample_half_width_mm / navcam_pixel_size / 1.5)
        max_radius = int(stub_hw * 0.9)
        circles = cv2.HoughCircles(
            gray, cv2.HOUGH_GRADIENT, dp=1.4, minDist=5,
            param1=50, param2=50, minRadius=min_radius, maxRadius=max_radius
        )
        
        if self.development_mode and circles is not None:
            output = stub_im.copy()
            circles_uint16 = np.uint16(np.around(circles))
            for (x, y, r) in circles_uint16[0, :]:
                cv2.circle(output, (x, y), r, (0, 255, 0), 2)
                cv2.circle(output, (x, y), 2, (0, 0, 255), 3)
            cv2.imshow('Detected Circles', output)
        
        # Filter circles by intensity
        stub_im_original = normalise_img(stub_im)
        filtered_circles = self._filter_circles_by_intensity(
            circles, stub_im_original, min_radius
        )
        
        if self.development_mode:
            filtered_debug = stub_im_original.copy()
            for (x, y, r) in filtered_circles:
                cv2.circle(filtered_debug, (x, y), r, (255, 0, 0), 2)
            cv2.imshow('Filtered Circles', filtered_debug)
        
        # Find best circle
        x, y, r = self._find_best_circle(filtered_circles, min_radius)
        
        # Compute output
        if len(filtered_circles) > 0 and x is not None and y is not None and r is not None:
            center_pos_eff = (
                self._center_pos + 
                (np.array([x, y]) - stub_hw) * navcam_pixel_size * np.array([1, -1])
            )
            sample_hw_mm = r * navcam_pixel_size * 0.9
            Ctape_coords = (center_pos_eff, sample_hw_mm)
            if self.verbose:
                logger.info('✅ C tape detected.')
        else:
            x, y = stub_hw, stub_hw
            r = int(self._sample_half_width_mm / navcam_pixel_size)
            Ctape_coords = None
            if self.verbose:
                logger.warning(
                    f'⚠️ The C tape could not be automatically detected. '
                    f'Using {tuple(float(v) for v in self._center_pos)} instead.'
                )
        
        # Draw detected circle on image
        cv2.circle(stub_im, (x, y), r, (0, 255, 0), 1)
        cv2.rectangle(stub_im, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)
        if self.development_mode:
            cv2.imshow('Detected region', stub_im)
        
        # Save result image if results_dir is set
        if self.results_dir:
            filename = os.path.join(self.results_dir, cnst.NAVCAM_IM_FILENAME + '.png')
            cv2.imwrite(filename, stub_im)
        
        return Ctape_coords
    
    
    def _filter_circles_by_intensity(self, circles, stub_im, min_radius):
        """Filter detected circles by average pixel intensity."""
        filtered_circles = []
        if circles is not None:
            circles_int = np.round(circles[0, :]).astype("int")
            for (x, y, r) in circles_int:
                avg_intensity = self._average_pixel_intensity(stub_im, (x, y), r)
                if avg_intensity <= 100:
                    filtered_circles.append((x, y, r))
        return filtered_circles
    
    
    def _average_pixel_intensity(self, image: np.ndarray, center: tuple, 
                                 radius: float) -> float:
        """Calculate average pixel intensity within a circular region."""
        # Convert color to grayscale if needed
        if len(image.shape) == 3 and image.shape[2] == 3:
            gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        else:
            gray_image = image
        
        x, y = int(center[0]), int(center[1])
        mask = np.zeros(gray_image.shape, dtype=np.uint8)
        cv2.circle(mask, (x, y), int(radius), 255, -1)
        
        masked_pixels = gray_image[mask == 255]
        if masked_pixels.size == 0:
            return 0.0
        return float(np.mean(masked_pixels))
    
    
    def _find_best_circle(self, filtered_circles, min_radius):
        """Find the best circle from filtered candidates."""
        x, y, r = None, None, None
        if len(filtered_circles) > 1:
            avg_x = np.mean([circle[0] for circle in filtered_circles])
            avg_y = np.mean([circle[1] for circle in filtered_circles])
            avg_center = (int(avg_x), int(avg_y))
            distances = [
                circle[2] - np.sqrt((circle[0] - avg_x) ** 2 + (circle[1] - avg_y) ** 2)
                for circle in filtered_circles
            ]
            intersection_radius = int(min(distances))
            x, y = avg_center
            r = max(intersection_radius, min_radius)
        elif len(filtered_circles) > 0:
            x, y, r = filtered_circles[0]
        
        return x, y, r
