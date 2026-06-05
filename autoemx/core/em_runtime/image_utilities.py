#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Image utilities and coordinate conversion module for electron microscopy.

This module provides image processing, saving, and coordinate transformation utilities
for scanning electron microscopy applications.

Functions
---------
convert_pixel_pos_to_mm
    Convert pixel coordinates to absolute stage coordinates.
convert_XS_coords_to_pixels
    Convert XS coordinates to pixel coordinates.
save_frame_image
    Save annotated and raw EM frame as multi-page TIFF.
normalise_img
    Normalize brightness of an RGB image.

Created on 2026
@author: Andrea
"""
import os
import json
import warnings
import numpy as np
import cv2
from PIL import Image

from autoemx.utils.helper import draw_scalebar
import autoemx.utils.constants as cnst


def convert_pixel_pos_to_mm(pos_pixels, im_width, im_height, pixel_size_um, 
                           current_pos, image_to_stage_coords_transform):
    """
    Convert a position from pixel coordinates to absolute stage coordinates in mm.
    
    Parameters
    ----------
    pos_pixels : array-like of float
        Position in pixel coordinates (x, y).
    im_width : int
        Image width in pixels.
    im_height : int
        Image height in pixels.
    pixel_size_um : float
        Pixel size in micrometers.
    current_pos : tuple
        Current position in mm (x, y).
    image_to_stage_coords_transform : array
        Transformation matrix for coordinate system conversion.
    
    Returns
    -------
    pos_abs_mm : ndarray of float
        Absolute position in millimeters (x, y).
    """
    # Calculate the center of the image in pixels
    center_pixels = np.array([im_width, im_height]) / 2
    
    # Compute the shift from the image center in pixels
    shift_pixels = pos_pixels - center_pixels
    
    # Convert the shift to micrometers
    shift_um = pixel_size_um * shift_pixels
    
    # Transform coordinates to match stage reference system
    shift_um_stage_coords = shift_um * image_to_stage_coords_transform
    
    # Calculate absolute position in mm for the EM
    pos_abs_mm = np.array(current_pos) + shift_um_stage_coords * 1e-3
    
    return pos_abs_mm


def convert_XS_coords_to_pixels(xy_coords, im_width, im_height, EM_driver):
    """
    Convert XY coordinates from the XS coordinate system to pixel coordinates.
    
    Parameters
    ----------
    xy_coords : tuple(float, float)
        The XY coordinates in the XS coordinate system.
    im_width : int
        Image width in pixels.
    im_height : int
        Image height in pixels.
    EM_driver : module
        EM driver module.
    
    Returns
    -------
    tuple(int, int)
        The corresponding (x, y) pixel coordinates as integers.
    """
    xy_coords_pixels = EM_driver.frame_rel_to_pixel_coords(
        xy_coords,
        im_width,
        im_height
    ).astype(int)[0]
    
    return xy_coords_pixels


def save_frame_image(frame_image, pixel_size_um, im_width, im_height,
                     sample_id, microscope_cfg, filename, results_dir,
                     im_annotations=None, scalebar=True,
                     image_extension: str = "tif",
                     save_raw_image: bool = True,
                     EM_driver=None,
                     auto_adjust_bc=True):
    """
    Save an annotated and optional raw electron microscopy (EM) frame.
    
    Generates a raw grayscale EM image and an annotated RGB version with optional 
    markers and scale bar. Both are saved into a single multi-page TIFF file.
    The annotated image is stored as the first page, and the raw image as the second page.
    
    Parameters
    ----------
    frame_image : np.ndarray
        Frame image to save. If None, acquires current frame.
    pixel_size_um : float
        Pixel size in micrometers.
    im_width : int
        Image width in pixels.
    im_height : int
        Image height in pixels.
    sample_id : str
        Canonical sample identifier used for metadata.
    microscope_cfg : MicroscopeConfig
        Microscope configuration.
    filename : str
        Name used for saved .tif image file (without extension).
    results_dir : str
        Directory in which to save the TIFF file.
    im_annotations : dict | list(dict) | None, optional
        Dictionary or list of dictionaries with annotations:
            - 'circle': (radius, xy_center, border_thickness)
            - 'text': (text, xy_coords)
    scalebar : bool, optional
        Whether to annotate the image with a scalebar (default: True).
    EM_driver : module, optional
        EM driver for image acquisition if frame_image is None.
    auto_adjust_bc : bool, optional
        Whether to auto-adjust brightness/contrast before saving (default: True).
    
    Notes
    -----
    - Images are saved as RGB to maximize compatibility.
    """
    # Determine save directory
    if not results_dir:
        warnings.warn(
            "No directory specified for saving frame image.",
            UserWarning
        )
        return
    
    if not isinstance(frame_image, np.ndarray):
        if EM_driver is None:
            raise ValueError("EM_driver required if frame_image is not provided")
        if auto_adjust_bc:
            EM_driver.auto_contrast_brightness()
        frame_image = EM_driver.get_image_data(im_width, im_height, 1)
    
    # Convert grayscale to RGB for annotation
    if len(frame_image.shape) == 2 or frame_image.shape[2] == 1:
        color_image = cv2.cvtColor(frame_image, cv2.COLOR_GRAY2RGB)
    else:
        color_image = frame_image.copy()
    
    # Draw annotations if provided
    an_circle_key = cnst.ANNOTATION_CIRCLE_KEY
    an_text_key = cnst.ANNOTATION_TEXT_KEY
    
    if im_annotations is not None:
        if isinstance(im_annotations, dict):
            im_annotations = [im_annotations]
        
        for ann_dict in im_annotations:
            # Add circles
            if an_circle_key in ann_dict.keys():
                radius, xy_center, border_thickness = ann_dict[an_circle_key]
                cv2.circle(color_image, tuple(xy_center), radius, (255, 0, 0), border_thickness)
            
            # Add label text
            if an_text_key in ann_dict.keys():
                text, text_xy = ann_dict[an_text_key]
                cv2.putText(
                    color_image,
                    text,
                    text_xy,
                    cv2.FONT_HERSHEY_SIMPLEX,
                    1,
                    (255, 0, 0),
                    2,
                    cv2.LINE_AA
                )
    
    # Add scale bar
    if scalebar:
        color_image = draw_scalebar(color_image, pixel_size_um)
    
    # Normalize extension and prepare paths.
    ext = str(image_extension).strip().lower().lstrip('.')
    if ext == "jpg":
        ext = "jpeg"

    save_path = os.path.join(results_dir, f"{filename}.{ext}")
    
    # Ensure dtype consistency (convert to uint8 if needed)
    if frame_image.dtype != np.uint8:
        frame_image_uint8 = (frame_image / frame_image.max() * 255).astype(np.uint8)
        if scalebar or im_annotations:
            color_image = (color_image / color_image.max() * 255).astype(np.uint8)
    else:
        frame_image_uint8 = frame_image
    
    # Convert grayscale to RGB for saving
    if frame_image_uint8.ndim == 2:
        frame_image_uint8 = cv2.cvtColor(frame_image_uint8, cv2.COLOR_GRAY2RGB)
    
    # Create image metadata
    image_description_d = {
        "sample_ID": sample_id,
        "microscope_ID": microscope_cfg.ID,
        "microscope_type": microscope_cfg.type,
        "detector": microscope_cfg.detector_type,
        "resolution": (im_width, im_height),
        "pixel_size_um": pixel_size_um
    }
    
    desc_str = json.dumps(image_description_d, ensure_ascii=True)
    
    # Convert numpy arrays to Pillow Image objects.
    if scalebar or im_annotations:
        im1 = Image.fromarray(color_image.astype('uint8'), mode='RGB')
    else:
        im1 = None
    im2 = Image.fromarray(frame_image_uint8.astype('uint8'), mode='RGB')

    annotated_image = im1 if im1 is not None else im2

    # TIFF supports multi-page saves. For non-TIFF formats, save raw in a sidecar file.
    if ext in {"tif", "tiff"}:
        if save_raw_image and im1 is not None:
            im1.save(
                save_path,
                format='TIFF',
                description=desc_str,
                save_all=True,
                append_images=[im2],
                compression=None,
            )
        else:
            annotated_image.save(save_path, format='TIFF', description=desc_str)
        return

    annotated_image.save(save_path, format=ext.upper())
    if save_raw_image:
        raw_path = os.path.join(results_dir, f"{filename}_raw.{ext}")
        im2.save(raw_path, format=ext.upper())


def normalise_img(img: np.ndarray, target_brightness: float = 128.0) -> np.ndarray:
    """
    Normalize brightness of an RGB image to a target brightness level.
    
    Parameters
    ----------
    img : np.ndarray
        Input RGB image (uint8).
    target_brightness : float, optional
        Desired average brightness (0–255), default 128.0.
    
    Returns
    -------
    np.ndarray
        Brightness-normalized RGB image (uint8).
    
    Raises
    ------
    TypeError
        If input is not a numpy array.
    ValueError
        If image does not have 3 channels after conversion.
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
