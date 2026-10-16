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
draw_annotations
    Draw circles and labels (e.g. X-ray spectrum spots) on an RGB image.
save_frame_image
    Save an EM frame, with optional annotations and scale bar.
normalise_img
    Normalize brightness of an RGB image.

Created on 2026
@author: Andrea
"""
import os
from pathlib import Path
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


def draw_annotations(color_image, im_annotations):
    """
    Draw circles and text labels on an RGB image, in place.

    Parameters
    ----------
    color_image : np.ndarray
        RGB image to annotate.
    im_annotations : dict | list(dict)
        Dictionary or list of dictionaries with annotations:
            - 'circle': (radius, xy_center, border_thickness), with border_thickness -1 for a filled dot
            - 'text': (text, xy_coords)

    Returns
    -------
    np.ndarray
        The annotated image.
    """
    an_circle_key = cnst.ANNOTATION_CIRCLE_KEY
    an_text_key = cnst.ANNOTATION_TEXT_KEY

    if isinstance(im_annotations, dict):
        im_annotations = [im_annotations]

    for ann_dict in im_annotations:
        # Add circles
        if an_circle_key in ann_dict.keys():
            radius, xy_center, border_thickness = ann_dict[an_circle_key]
            cv2.circle(color_image, tuple(int(c) for c in xy_center), radius, (255, 0, 0), border_thickness)

        # Add label text
        if an_text_key in ann_dict.keys():
            text, text_xy = ann_dict[an_text_key]
            cv2.putText(
                color_image,
                text,
                tuple(int(c) for c in text_xy),
                cv2.FONT_HERSHEY_SIMPLEX,
                1,
                (255, 0, 0),
                2,
                cv2.LINE_AA
            )
    return color_image


#: Suffix of the raw copy of an annotated image saved in a format without pages (e.g. ``<name>_raw.png``).
RAW_IMAGE_SUFFIX = "_raw"


def load_raw_saved_image(path):
    """
    RGB array of the un-annotated version of an image saved by :func:`save_frame_image`: its ``_raw`` copy if
    present, the second (raw) page of a two-page TIFF, or the image itself.
    """
    path = Path(path)
    raw_sidecar = path.with_name(f"{path.stem}{RAW_IMAGE_SUFFIX}{path.suffix}")
    source = raw_sidecar if raw_sidecar.exists() else path
    with Image.open(source) as im:
        # Two-page TIFFs saved during acquisition hold the annotated image first and the raw image second
        if source == path and getattr(im, "n_frames", 1) > 1:
            im.seek(1)
        return np.array(im.convert("RGB"))


def particle_mark_annotations(particle_marks, im_width, text_margin=30, text_pos_scale=0.9):
    """
    Annotations of the particles marked in a frame: a circle around each and, if it has one, its particle id at
    the top-right of the circle (at its centre near the top or right edge of the image).

    Parameters
    ----------
    particle_marks : sequence of (x, y, radius, particle_id or None)
        Particles, in pixels (see ``ParticleStatsFrame.particle_marks_px``).
    im_width : int
        Image width, in pixels.

    Returns
    -------
    list of dict
        Annotations for :func:`draw_annotations`.
    """
    annotations = []
    for x, y, radius, particle_id in particle_marks:
        annotation = {cnst.ANNOTATION_CIRCLE_KEY: (int(radius), (int(x), int(y)), 2)}
        if particle_id is not None:
            x_text = int(x + radius * text_pos_scale)
            y_text = int(y - radius * text_pos_scale)
            if x_text > im_width - text_margin or y_text < text_margin:
                x_text, y_text = int(x), int(y)
            annotation[cnst.ANNOTATION_TEXT_KEY] = (str(particle_id), (x_text, y_text))
        annotations.append(annotation)
    return annotations


def xsp_spot_annotation(spectrum_id, xy_center):
    """Annotation of an X-ray spectrum spot: red dot with the spectrum ID at its top-left."""
    x, y = int(xy_center[0]), int(xy_center[1])
    return {
        cnst.ANNOTATION_TEXT_KEY: (str(spectrum_id), (x - 30, y - 15)),
        cnst.ANNOTATION_CIRCLE_KEY: (10, (x, y), -1),
    }


def save_frame_image(frame_image, pixel_size_um, im_width, im_height,
                     sample_id, microscope_cfg, filename, results_dir,
                     im_annotations=None, scalebar=True,
                     image_extension: str = "tif",
                     save_raw_image: bool = False,
                     EM_driver=None,
                     auto_adjust_bc=True):
    """
    Save an electron microscopy (EM) frame, with optional markers and scale bar.

    If the image is annotated and ``save_raw_image`` is True, the raw image is saved as well:
    as second page of the TIFF file, or as ``<filename>_raw.<ext>`` for other formats.
    The image metadata (sample, microscope, resolution, pixel size) is stored in the
    TIFF/PNG image description, e.g. to draw scale bars on raw images afterwards.
    
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
    if im_annotations is not None:
        draw_annotations(color_image, im_annotations)
    
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
    
    # Convert numpy arrays to Pillow Image objects (RGB, from their 3 channels).
    is_annotated = bool(scalebar or im_annotations)
    raw_image = Image.fromarray(frame_image_uint8.astype('uint8'))
    image = Image.fromarray(color_image.astype('uint8')) if is_annotated else raw_image
    # A raw copy is only useful when the saved image is annotated
    save_raw_copy = save_raw_image and is_annotated

    if ext in {"tif", "tiff"}:
        if save_raw_copy:
            # Annotated image as first page, raw image as second page
            image.save(save_path, format='TIFF', description=desc_str, save_all=True,
                       append_images=[raw_image], compression=None)
        else:
            image.save(save_path, format='TIFF', description=desc_str)
        return

    if ext == "png":
        from PIL.PngImagePlugin import PngInfo

        png_info = PngInfo()
        png_info.add_text("Description", desc_str)
        image.save(save_path, format="PNG", pnginfo=png_info)
    else:
        image.save(save_path, format=ext.upper())
    if save_raw_copy:
        raw_path = os.path.join(results_dir, f"{filename}{RAW_IMAGE_SUFFIX}.{ext}")
        raw_image.save(raw_path, format=ext.upper())


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
