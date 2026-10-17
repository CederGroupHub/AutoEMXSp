#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing of segmentation models that give one mask per object (SAM2, micro-sam): from the masks of a frame,
a label image of the particles whose size can be measured.

Promptable models such as SAM segment everything: particles, but also the substrate, agglomerates of particles
(as well as each particle in them), surface features of particles (pits, debris) and particles hidden behind
others. The masks are filtered as follows:

1. **Substrate**: masks not brighter on average than an Otsu threshold of the frame (particles are brighter than the
   substrate in SE/BSE images), or covering more than ``max_frame_fraction`` of the frame, are dropped.
2. **Masks inside a larger mask** (at least ``containment`` of their area in it, holes of the larger mask filled),
   from the smallest up: if the smaller masks cover at least ``agglomerate_coverage`` of the larger one, it is an
   agglomerate of particles and is dropped (its particles kept); otherwise the smaller masks are features of its
   surface and are dropped (the larger mask's holes filled).
3. **Occluded particles**: the masks are painted from the largest to the smallest (smaller particles in front of
   larger ones); a mask left with less than ``min_visible_fraction`` of its pixels is a particle partly hidden by
   others, of unknown size, and is dropped.
4. Optionally, masks with a solidity (area / convex hull area) below ``min_solidity`` (off by default: it would
   bias the statistics against irregular particles).

Particles touching the frame edge and outside the size range are left out afterwards by the particle finder.

Created on Sat Oct 17 2026
@author: Andrea
"""

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import cv2
import numpy as np

#: Options of the filters (set in PowderMeasurementConfig.par_segmentation_params) and their defaults
FILTER_DEFAULTS: Dict[str, Any] = {
    "min_visible_fraction": 0.85,
    "agglomerate_coverage": 0.7,
    "containment": 0.9,
    "drop_dark_masks": True,
    "max_frame_fraction": 0.5,
    "min_solidity": None,
}


@dataclass
class _Mask:
    """A mask cropped to its bounding box (rows y0:y1, columns x0:x1 of the frame)."""

    y0: int
    x0: int
    y1: int
    x1: int
    pixels: np.ndarray  # bool, shape (y1 - y0, x1 - x0)

    @property
    def area(self) -> int:
        return int(self.pixels.sum())

    def filled(self) -> "_Mask":
        """The mask with its holes filled."""
        from scipy import ndimage

        return _Mask(self.y0, self.x0, self.y1, self.x1, ndimage.binary_fill_holes(self.pixels))

    def overlap(self, other: "_Mask") -> int:
        """Number of pixels in both masks."""
        y0, x0 = max(self.y0, other.y0), max(self.x0, other.x0)
        y1, x1 = min(self.y1, other.y1), min(self.x1, other.x1)
        if y0 >= y1 or x0 >= x1:
            return 0
        a = self.pixels[y0 - self.y0:y1 - self.y0, x0 - self.x0:x1 - self.x0]
        b = other.pixels[y0 - other.y0:y1 - other.y0, x0 - other.x0:x1 - other.x0]
        return int(np.count_nonzero(a & b))


def _crop(mask: np.ndarray) -> Optional[_Mask]:
    ys = np.flatnonzero(mask.any(axis=1))
    if ys.size == 0:
        return None
    xs = np.flatnonzero(mask.any(axis=0))
    y0, y1, x0, x1 = int(ys[0]), int(ys[-1]) + 1, int(xs[0]), int(xs[-1]) + 1
    return _Mask(y0, x0, y1, x1, mask[y0:y1, x0:x1].astype(bool))


def _masks_of_labels(labels: np.ndarray) -> List[_Mask]:
    """The masks of the labels of a label image (0 = background), cropped."""
    from scipy import ndimage

    labels = labels.astype(np.int32, copy=False)
    masks = []
    for label, slc in enumerate(ndimage.find_objects(labels), start=1):
        if slc is not None:
            masks.append(_Mask(slc[0].start, slc[1].start, slc[0].stop, slc[1].stop, labels[slc] == label))
    return masks


def filter_options(params: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Filter options from model options (unknown keys ignored), with the defaults."""
    params = params or {}
    return {key: params.get(key, default) for key, default in FILTER_DEFAULTS.items()}


def _solidity(mask: _Mask) -> float:
    contours, _ = cv2.findContours(mask.pixels.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return 0.0
    hull_area = cv2.contourArea(cv2.convexHull(np.vstack(contours)))
    return mask.area / hull_area if hull_area > 0 else 1.0


def _substrate_threshold(image: np.ndarray) -> float:
    gray = image if image.dtype == np.uint8 else cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
    threshold, _ = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return float(threshold) if image.dtype == np.uint8 else float(
        image.min() + threshold / 255.0 * (float(image.max()) - float(image.min())))


def particle_labels(
    masks: Union[Iterable[np.ndarray], np.ndarray],
    image: np.ndarray,
    min_visible_fraction: float = FILTER_DEFAULTS["min_visible_fraction"],
    agglomerate_coverage: float = FILTER_DEFAULTS["agglomerate_coverage"],
    containment: float = FILTER_DEFAULTS["containment"],
    drop_dark_masks: bool = FILTER_DEFAULTS["drop_dark_masks"],
    max_frame_fraction: float = FILTER_DEFAULTS["max_frame_fraction"],
    min_solidity: Optional[float] = FILTER_DEFAULTS["min_solidity"],
) -> np.ndarray:
    """
    Label image (int32, background 0, particles 1, 2, ...) of the particles among the masks of a frame (see the
    module documentation).

    Parameters
    ----------
    masks : iterable of ndarray, or ndarray
        One bool mask per object, of the frame's shape (they may overlap), or a label image of the objects (integer,
        0 = background).
    image : ndarray
        The greyscale frame image.
    """
    if image.ndim == 3:
        image = image[..., 0]
    frame_area = image.shape[0] * image.shape[1]
    threshold = _substrate_threshold(image) if drop_dark_masks else None

    # 1. Substrate
    if isinstance(masks, np.ndarray) and masks.ndim == 2 and masks.dtype != np.bool_:
        cropped = _masks_of_labels(masks)
    else:
        cropped = (_crop(np.asarray(full, dtype=bool)) for full in masks)
    kept: List[_Mask] = []
    for mask in cropped:
        if mask is None or mask.area > max_frame_fraction * frame_area:
            continue
        if threshold is not None and image[mask.y0:mask.y1, mask.x0:mask.x1][mask.pixels].mean() <= threshold:
            continue
        kept.append(mask)

    # 2. Masks inside a larger mask (with its holes filled: in a label image, the masks inside another are holes in
    # it), from the smallest up (nested ones resolved first)
    filled = [m.filled() for m in kept]
    order = sorted(range(len(kept)), key=lambda k: filled[k].area)
    kept, filled = [kept[k] for k in order], [filled[k] for k in order]
    dropped = set()
    for i, outer in enumerate(filled):
        inner = [j for j in range(i) if j not in dropped and kept[j].area < outer.area
                 and kept[j].overlap(outer) >= containment * kept[j].area]
        if not inner:
            continue
        union = np.zeros_like(outer.pixels)
        for j in inner:
            m = kept[j]
            y0, x0 = max(m.y0, outer.y0), max(m.x0, outer.x0)
            y1, x1 = min(m.y1, outer.y1), min(m.x1, outer.x1)
            union[y0 - outer.y0:y1 - outer.y0, x0 - outer.x0:x1 - outer.x0] |= \
                m.pixels[y0 - m.y0:y1 - m.y0, x0 - m.x0:x1 - m.x0]
        if np.count_nonzero(union & outer.pixels) >= agglomerate_coverage * outer.area:
            dropped.add(i)  # agglomerate of the particles inside
        else:
            dropped.update(inner)  # features of the particle's surface, part of the particle
            kept[i] = outer
    kept = [m for i, m in enumerate(kept) if i not in dropped]

    # 4. Solidity
    if min_solidity:
        kept = [m for m in kept if _solidity(m) >= min_solidity]

    # 3. Occluded particles: painted from the largest to the smallest, twice (after dropping the hidden ones, those
    # they hid are visible again)
    for _ in range(2):
        kept.sort(key=lambda m: m.area, reverse=True)
        labels = _paint(kept, image.shape)
        visible = np.bincount(labels.ravel(), minlength=len(kept) + 1)[1:]
        keep = [visible[i] >= min_visible_fraction * m.area for i, m in enumerate(kept)]
        if all(keep):
            break
        kept = [m for m, k in zip(kept, keep) if k]
    labels = _paint(kept, image.shape)
    # Consecutive labels (a mask may be entirely hidden after the last painting)
    present = np.unique(labels)
    present = present[present > 0]
    lookup = np.zeros(len(kept) + 1, dtype=np.int32)
    lookup[present] = np.arange(1, present.size + 1, dtype=np.int32)
    return lookup[labels]


def _paint(masks: Sequence[_Mask], shape) -> np.ndarray:
    labels = np.zeros(shape[:2], dtype=np.int32)
    for i, m in enumerate(masks, start=1):
        labels[m.y0:m.y1, m.x0:m.x1][m.pixels] = i
    return labels
