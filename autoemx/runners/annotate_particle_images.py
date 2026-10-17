#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate the particle and frame images of a sample with the X-ray spectrum spots and the particles found.

During acquisition, particle and frame images are saved without annotations (unless
``annotate_particle_images=True``), while the pixel position of every spectrum spot, and
of every particle found in a frame, is stored in the sample ledger. This runner draws each
spot as a dot with its spectrum ID, and each particle found in a frame (particle size
statistics, or frames searched for EDS particles) as a circle with its particle number, as
done during acquisition, plus a scale bar when the pixel size is known.

Annotated copies are written to ``SEM images/annotated``; the original images are not
modified.

Import this module in your own code and call `annotate_particle_images()`, or use the
script ``Annotate_Particle_Images.py``.

Workflow
--------
- Loads the sample ledger and the images ``SEM images/<sample>_par<N>_fr<F>_xyspots.<ext>``
  (or ``<sample>_fr<F>_xyspots.<ext>`` for grid/linescan acquisitions).
- Finds the spectra acquired on each image: same particle (or same frame, for images without particle).
- Draws a dot + spectrum ID at the pixel position of each spectrum, and a scale bar.
- Annotates the frame images recorded in the ledger (``particle_stats``, ``particle_stats_progress``,
  ``eds_size_scans``) with their particles (``particle_marks_px``).

Notes
-----
- The un-annotated image is used: the ``<name>_raw.<ext>`` copy if present, or the raw page of a
  two-page TIFF. Images saved with annotations and without raw copy are annotated again.
- The pixel size used for the scale bar is read from the ledger (``image_pixel_size_um`` of the
  particle) or, for older samples and images without particle, from the image metadata.
  Without it, no scale bar is drawn.

Created on Mon Oct  5 2026

@author: Andrea
"""

import json
import logging
from types import SimpleNamespace
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from PIL import Image

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.core.em_runtime.image_utilities import (
    RAW_IMAGE_SUFFIX,
    draw_annotations,
    load_raw_saved_image,
    particle_mark_annotations,
    xsp_spot_annotation,
)
from autoemx.utils.helper import draw_scalebar, get_sample_dir, parse_xsp_spots_image_name

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["annotate_particle_images"]

_IMAGE_EXTENSIONS = (".png", ".tif", ".tiff", ".jpg", ".jpeg", ".webp", ".bmp")


def _resolve_sample_dir(
    sample_path: Optional[str],
    sample_ID: Optional[str],
    results_path: Optional[str],
) -> Path:
    """Sample directory from its full path, or from its ID searched under ``results_path``."""
    if sample_path:
        sample_dir = Path(sample_path).expanduser().resolve()
    elif sample_ID:
        if results_path is None:
            results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)
        sample_dir = Path(get_sample_dir(results_path, sample_ID)).resolve()
    else:
        raise ValueError("Provide either sample_path, or sample_ID (and optionally results_path).")

    ledger_path = sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
    if not ledger_path.exists():
        raise FileNotFoundError(f"No {ledger_path.name} found in {sample_dir}. Is this an AutoEMX sample folder?")
    return sample_dir


def _pixel_size_from_metadata(path: Path) -> Optional[float]:
    """Pixel size (um) stored in the image description (TIFF tag or PNG text) by AutoEMX."""
    try:
        with Image.open(path) as im:
            desc = im.info.get("Description")
            if desc is None and hasattr(im, "tag_v2"):
                desc = im.tag_v2.get(270)  # TIFF ImageDescription
        value = json.loads(desc).get("pixel_size_um") if desc else None
        return float(value) if value else None
    except Exception:
        return None


def _spots_by_image(ledger) -> Tuple[Dict[int, list], Dict[str, list]]:
    """Spectrum spots grouped by particle ID and, for spots without particle, by frame ID."""
    by_particle: Dict[int, list] = {}
    by_frame: Dict[str, list] = {}
    for idx, sp in enumerate(ledger.spectra):
        details = sp.acquisition_details
        coords = details.spot_coordinates if details is not None else None
        if coords is None or coords.pixel_coordinates is None:
            continue
        spectrum_id = sp.spectrum_id if sp.spectrum_id is not None else str(idx)
        spot = (spectrum_id, coords.pixel_coordinates, details.frame_id)
        if details.particle_id is not None:
            by_particle.setdefault(int(details.particle_id), []).append(spot)
        elif details.frame_id is not None:
            by_frame.setdefault(str(details.frame_id), []).append(spot)
    return by_particle, by_frame


def _recorded_frames(ledger) -> list:
    """
    Frames of the ledger with an image, with the particles to mark (x, y, radius, index in the frame): those of the
    segmentation of the active statistics result (else the last segmentation of each frame).
    """
    result = ledger.size_result
    frames = []
    for frame in ledger.frames:
        if not frame.image_file or not frame.segmentations:
            continue
        preferred = None
        if result is not None:
            preferred = result.eds_segmentation_id if frame.purpose == "eds_search" else result.segmentation_id
        segmentation = frame.segmentation(preferred) or frame.segmentations[-1]
        marks = [(x, y, r, i) for i, (x, y, r) in enumerate(segmentation.particle_marks_px)]
        frames.append(SimpleNamespace(image_file=frame.image_file, pixel_size_um=frame.pixel_size_um,
                                      particle_marks_px=marks))
    return frames


def annotate_particle_images(
    sample_path: Optional[str] = None,
    sample_ID: Optional[str] = None,
    results_path: Optional[str] = None,
    output_dir: Optional[str] = None,
    scalebar: bool = True,
    overwrite: bool = True,
    verbose: bool = True,
) -> List[str]:
    """
    Write copies of the particle images of a sample annotated with the X-ray spectrum spots, and of its frame
    images annotated with the particles found in them.

    Parameters
    ----------
    sample_path : str, optional
        Full path to the sample folder (containing ledger.json).
    sample_ID : str, optional
        Sample identifier, used with ``results_path`` when ``sample_path`` is not given.
    results_path : str, optional
        Directory under which the ``sample_ID`` folder is searched. Default: ./Results
    output_dir : str, optional
        Folder for the annotated images. Default: ``SEM images/annotated`` in the sample folder.
    scalebar : bool, optional
        Whether to draw a scale bar (when the pixel size is known). Default: True.
    overwrite : bool, optional
        Whether to overwrite existing annotated images. Default: True.
    verbose : bool, optional
        Whether to log progress. Default: True.

    Returns
    -------
    written : list of str
        Paths of the annotated images.
    """
    sample_dir = _resolve_sample_dir(sample_path, sample_ID, results_path)
    images_dir = sample_dir / cnst.IMAGES_DIR
    out_dir = Path(output_dir).expanduser() if output_dir else images_dir / cnst.ANNOTATED_IMAGES_SUBDIR
    if not images_dir.is_dir():
        logging.warning(f"No '{cnst.IMAGES_DIR}' folder in {sample_dir}. Nothing to annotate.")
        return []

    ledger = load_sample_ledger(str(sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"))
    by_particle, by_frame = _spots_by_image(ledger)
    pixel_size_by_particle = {
        p.id: p.image_pixel_size_um for p in ledger.particles if p.image_pixel_size_um is not None
    }

    written: List[str] = []
    n_without_spots = 0
    n_without_scale = 0
    for path in sorted(images_dir.iterdir()):
        if path.suffix.lower() not in _IMAGE_EXTENSIONS or path.stem.endswith(RAW_IMAGE_SUFFIX):
            continue
        parsed = parse_xsp_spots_image_name(path.stem)
        if parsed is None:
            continue  # Not an image with spectrum spots (e.g. frame overview)
        particle_id, frame_id = parsed
        if particle_id is not None:
            spots = by_particle.get(particle_id, [])
            # Particle IDs are unique; check the frame only when it is recorded for both
            spots = [s for s in spots if s[2] is None or str(s[2]) == str(frame_id)] or spots
        else:
            spots = by_frame.get(str(frame_id), [])
        if not spots:
            n_without_spots += 1
            if verbose:
                logging.warning(f"No spectrum spots in the ledger for {path.name}. Skipped.")
            continue

        out_path = out_dir / path.name
        if out_path.exists() and not overwrite:
            continue

        image = load_raw_saved_image(path)
        draw_annotations(image, [xsp_spot_annotation(sid, xy) for sid, xy, _ in spots])
        if scalebar:
            pixel_size = pixel_size_by_particle.get(particle_id) or _pixel_size_from_metadata(path)
            if pixel_size:
                image = draw_scalebar(image, pixel_size)
            else:
                n_without_scale += 1

        out_dir.mkdir(parents=True, exist_ok=True)
        Image.fromarray(image.astype(np.uint8)).save(out_path)
        written.append(str(out_path))
        if verbose:
            logging.info(f"Annotated {path.name} with {len(spots)} spectrum spot(s).")

    # Frames: particles found in them
    for frame in _recorded_frames(ledger):
        path = images_dir / frame.image_file
        out_path = out_dir / frame.image_file
        if not path.exists() or (out_path.exists() and not overwrite):
            continue
        image = load_raw_saved_image(path)
        draw_annotations(image, particle_mark_annotations(frame.particle_marks_px, image.shape[1]))
        if scalebar:
            image = draw_scalebar(image, frame.pixel_size_um)
        out_dir.mkdir(parents=True, exist_ok=True)
        Image.fromarray(image.astype(np.uint8)).save(out_path)
        written.append(str(out_path))
        if verbose:
            logging.info(f"Annotated {path.name} with {len(frame.particle_marks_px)} particle(s).")

    if verbose:
        logging.info(f"{len(written)} annotated image(s) written to {out_dir}")
        if n_without_spots:
            logging.info(f"{n_without_spots} image(s) skipped: no spectrum spots recorded in the ledger.")
        if n_without_scale:
            logging.info(f"{n_without_scale} image(s) without scale bar: pixel size not found in ledger or image metadata.")
    return written
