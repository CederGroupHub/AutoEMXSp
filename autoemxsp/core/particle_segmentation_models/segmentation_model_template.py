#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 09:34:39 2025

@author: Andrea
"""

import cv2
import numpy as np


def segment_particles(frame_image : np.array,
                      powder_meas_config : 'PowderMeasurementConfig' = None,
                      save_image : bool = False,
                      EM : 'EM_controller' = None):
    """
    Segments particles in the given frame image based on defined criteria.
    
    This function applies image processing techniques, such as thresholding 
    and contour filling, to generate a binary mask that represents the 
    detected particles in the input frame image. 

    Parameters
    ----------
    frame_image : ndarray
        A grayscale input image of the current frame containing particles to be detected.
        
    powder_meas_config : PowderMeasurementConfig object, optional
        A configuration object that contains relevant parameters for segmentation such as:
        - `par_brightness_thresh` (int): The brightness threshold used for binary segmentation.
        Additional parameters can be added as needed to fine-tune segmentation behavior.
    
    save_image : bool, optional
        Optionally save masked image, through EM_controller function. Default: False
    
    EM_controller : EM_controller object, optional
        Used to optionally save segmented image
        
    Returns
    -------
    par_mask : ndarray
        Choose either a binary mask or a labeled image:

        - **Binary mask (uint8 or bool)**
          Background pixels are set to ``0`` (black), and particle pixels
          are set to ``255`` (white).
          OK TO USE ONLY IF PARTICLE EDGES ARE NOT TOUCHING, because it will use
          cv2.ConnectedComponentsWithStats(img, connectivity=8) to discern particles

        - **Labeled image (uint8 or higher precision)**
            USE THIS IF DIFFERENT PARTICLES HAVE CONTIGOUS EDGES
          Background pixels are set to ``0`` (black). Each detected particle
          is assigned a unique positive integer label, starting from ``1`` and
          increasing consecutively (same format as the ``labels`` output from
          ``cv2.connectedComponents``).

    Note
    ----
    Particles may be segmented as contours or filled regions.
    A filling step is always performed later at the level of EM_Particle_Finder.
    """
    
    # Example thresholding to create a binary mask
    _, par_mask = cv2.threshold(frame_image, powder_meas_config.par_brightness_thresh, 255, cv2.THRESH_BINARY)
    
    # Use this code to masked image
    if save_image and EM:
        filename = f"{EM.sample_cfg.ID}_fr{EM.current_frame_label}_ml_mask"
        EM.save_frame_image(filename, frame_image = par_mask)
    
    return par_mask