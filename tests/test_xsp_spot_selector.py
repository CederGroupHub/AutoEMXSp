#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for programmatic X-ray spot selection callbacks."""
import pytest

from autoemx.core.em_runtime.xsp_spot_selection import validate_xsp_spot_pixels


def test_validate_xsp_spot_pixels_accepts_in_bounds_points():
    spots = validate_xsp_spot_pixels(
        [(4, 5), (8, 6)],
        im_width=12,
        im_height=10,
    )
    assert spots == [(4, 5), (8, 6)]


def test_validate_xsp_spot_pixels_rejects_out_of_bounds():
    with pytest.raises(ValueError, match="outside image bounds"):
        validate_xsp_spot_pixels([(12, 5)], im_width=12, im_height=10)


def test_validate_xsp_spot_pixels_deduplicates():
    spots = validate_xsp_spot_pixels(
        [(3, 4), (3, 4), (5, 6)],
        im_width=12,
        im_height=10,
    )
    assert spots == [(3, 4), (5, 6)]
