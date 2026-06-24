#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Programmatic X-ray spot selection on automatically navigated particles.

When ``PowderMeasurementConfig.par_spot_selection_mode`` is ``'callback'``, the
framework navigates to each particle automatically and invokes ``xsp_spot_selector``
repeatedly on the same particle until the callback returns an empty list.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, List, Sequence, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder

XSpSpotSelectorCallback = Callable[["XSpSpotSelectionContext"], Sequence[Tuple[int, int]]]


@dataclass(frozen=True)
class XSpSpotSelectionContext:
    """
    Context passed to ``xsp_spot_selector`` after automatic particle navigation.

    The callback is invoked once per spot batch on the current particle. Return
    a non-empty list of pixel coordinates to measure, then wait for the next
    invocation to supply more spots. Return ``[]`` when finished with this
    particle (including to skip it without measuring).

    ``ref_image`` is the frame captured at spot-selection time. It is the same
    image used as the drift-correction reference during subsequent acquisitions
    on this particle (see ``EM_Controller.ref_image``).
    """

    particle_id: int
    n_tot_sp_collected: int
    ref_image: Any
    pixel_size_um: float
    frame_label: str
    im_width: int
    im_height: int
    particle_finder: "EM_Particle_Finder"


def validate_xsp_spot_pixels(
    spots: Sequence[Tuple[int, int]],
    *,
    im_width: int,
    im_height: int,
) -> List[Tuple[int, int]]:
    """
    Validate and normalize programmatic spot pixel coordinates.

    Raises
    ------
    ValueError
        If a spot is outside the image bounds.
    """
    if im_width <= 0 or im_height <= 0:
        raise ValueError("Image dimensions must be positive.")

    validated: List[Tuple[int, int]] = []
    seen = set()
    for spot in spots:
        if len(spot) != 2:
            raise ValueError(f"Each spot must be an (x, y) pair; got {spot!r}.")
        x, y = int(spot[0]), int(spot[1])
        if not (0 <= x < im_width and 0 <= y < im_height):
            raise ValueError(
                f"Spot ({x}, {y}) is outside image bounds "
                f"(width={im_width}, height={im_height})."
            )
        key = (x, y)
        if key not in seen:
            seen.add(key)
            validated.append(key)
    return validated
