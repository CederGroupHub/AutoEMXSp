#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Outlines of the particles found in a frame by a segmentation, to show the segmentation on the frame image (GUI:
Analysis > Particle sizes, frame browser).

One small JSON file per frame and segmentation, outside the ledger (loaded only when the frame is shown):
``<sample>/particle_outlines/segmentation_<id>/<frame_id>.json``, with the outline of each particle counted in the
statistics (in the order of the frame's particles: particle ``i`` of the frame) and of the other particles found
(touching the frame edge, or outside the diameters searched in the frame). Each outline is a polygon
``[x0, y0, x1, y1, ...]`` in pixels of the frame image.

Created on Sat Oct 17 2026
@author: Andrea
"""

import json
import os
import re
from typing import Dict, Iterable, List, Optional

import cv2
import numpy as np

__all__ = ["OUTLINES_DIR", "frame_outlines", "load_outlines", "outlines_path", "save_outlines"]

OUTLINES_DIR = "particle_outlines"


def outlines_path(sample_dir: str, segmentation_id: int, frame_id: str) -> str:
    """File of the particle outlines of a frame under a segmentation."""
    name = re.sub(r'[\\/:*?"<>|]', "_", str(frame_id))
    return os.path.join(sample_dir, OUTLINES_DIR, f"segmentation_{int(segmentation_id)}", f"{name}.json")


def _outline(mask: np.ndarray, x0: int, y0: int) -> List[int]:
    """Outline of a mask (its largest external contour, simplified to within a pixel) as [x0, y0, x1, y1, ...]."""
    padded = np.pad(mask.astype(np.uint8), 1)
    contours, _ = cv2.findContours(padded, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return []
    contour = cv2.approxPolyDP(max(contours, key=cv2.contourArea), 0.8, True).reshape(-1, 2)
    return (contour + (x0 - 1, y0 - 1)).astype(int).ravel().tolist()


def frame_outlines(labels: np.ndarray, counted: Iterable[int], others: Iterable[int] = ()) -> Dict[str, list]:
    """
    Outlines of the particles of a label image: ``counted`` (labels of the particles counted, in the order of the
    frame's particles) and ``others`` (labels of the other particles found).
    """
    from scipy import ndimage

    labels = np.asarray(labels).astype(np.int32, copy=False)
    slices = ndimage.find_objects(labels)

    def outline(label: int) -> List[int]:
        slc = slices[label - 1] if 0 < label <= len(slices) else None
        if slc is None:
            return []
        return _outline(labels[slc] == label, slc[1].start, slc[0].start)

    return {"counted": [outline(int(i)) for i in counted], "others": [outline(int(i)) for i in others]}


def save_outlines(sample_dir: str, segmentation_id: int, frame_id: str, outlines: Dict[str, list]) -> str:
    """Save the particle outlines of a frame under a segmentation; returns the file."""
    path = outlines_path(sample_dir, segmentation_id, frame_id)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"frame_id": frame_id, "segmentation_id": int(segmentation_id), **outlines}, fh,
                  separators=(",", ":"))
    return path


def load_outlines(sample_dir: str, segmentation_id: Optional[int], frame_id: str) -> Optional[Dict[str, list]]:
    """Particle outlines of a frame under a segmentation, or None if not saved (frames segmented before outlines
    were saved: segmenting them again with the same settings adds them)."""
    if segmentation_id is None:
        return None
    try:
        with open(outlines_path(sample_dir, segmentation_id, frame_id), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None
