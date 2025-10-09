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
        A binary mask of detected particles where:
        - Pixels corresponding to particles are set to 255 (white).
        - Background pixels are set to 0 (black).

    Note
    ----
    Particles may be segmented as contours or filled regions.
    A filling step is always performed later at the level of EM_Particle_Finder.
    """
    
    # Example thresholding to create a binary mask
    _, par_mask = cv2.threshold(frame_image, powder_meas_config.par_brightness_thresh, 255, cv2.THRESH_BINARY)
    
    if save_image and EM:
        filename = f"fr{EM.frame_labels[EM._frame_cntr-1]}_ml_mask"
        EM.save_frame_image(filename, frame_image = par_mask)
    
    return par_mask