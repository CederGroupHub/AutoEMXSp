#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for multi-scale particle-statistics helpers and ledger-backed PSD export."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from autoemx.config.ledger_schemas import Coordinate2D, ParticleInfo
from autoemx.config.runtime_configs import PowderMeasurementConfig
from autoemx.core.em_runtime.frame_navigator import FrameNavigator
from autoemx.core.em_runtime.particle_finder import (
    EM_Particle_Finder,
    diameter_um_to_area_um2,
    frame_width_um_from_max_area,
    min_valid_diameter_um,
    plan_size_bands_um,
)
from autoemx.utils.helper import AlphabetMapper
import autoemx.utils.constants as cnst


def test_size_bands_from_minimum_particle_pixels():
    # A particle of 80 pixels has a diameter of ~10 pixels
    assert min_valid_diameter_um(0.5, 80) == pytest.approx(2 * 0.5 * np.sqrt(80 / np.pi))

    bands = plan_size_bands_um(0.1, 50.0, im_width_px=1920, min_par_pixels=80)
    assert len(bands) == 3
    assert bands[0][1] == 50.0 and bands[-1][0] == 0.1
    for i, (d_lo, d_hi, fw_um) in enumerate(bands):
        assert fw_um == pytest.approx(frame_width_um_from_max_area(diameter_um_to_area_um2(d_hi)))
        assert d_lo == pytest.approx(max(0.1, min_valid_diameter_um(fw_um / 1920, 80)))  # smallest valid particles
        if i:
            assert d_hi == bands[i - 1][0]  # contiguous, largest first
    # Frames 10 x the largest particle: each band spans 1920 / (20 * sqrt(80 / pi)) ~ 19 x in diameter
    assert bands[1][1] / bands[1][0] == pytest.approx(1920 / (20 * np.sqrt(80 / np.pi)))

    # More pixels per particle (or fewer per frame) -> narrower bands
    assert len(plan_size_bands_um(0.1, 50.0, im_width_px=1920, min_par_pixels=800)) > len(bands)
    assert plan_size_bands_um(2.0, 2.0, im_width_px=1920, min_par_pixels=80) == [(2.0, 2.0, 20.0)]
    with pytest.raises(ValueError, match="reduce min_par_pixels"):
        plan_size_bands_um(0.1, 50.0, im_width_px=100, min_par_pixels=80)


def test_diameter_area_and_frame_width_helpers():
    area = diameter_um_to_area_um2(2.0)
    assert area == pytest.approx(np.pi)
    # radius 1 um -> FOV 20 um
    assert frame_width_um_from_max_area(area) == pytest.approx(20.0)
    huge = diameter_um_to_area_um2(100.0)
    assert frame_width_um_from_max_area(huge) == pytest.approx(500.0)


def test_lowercase_alphabet_mapper():
    mapper = AlphabetMapper(alphabet="abcdefghijklmnopqrstuvwxyz")
    assert mapper.get_letter(0) == "a"
    assert mapper.get_letter(25) == "z"
    assert mapper.get_letter(26) == "aa"


def test_configure_particle_subframes_hierarchical_ids():
    nav = object.__new__(FrameNavigator)
    nav.EM_controller = SimpleNamespace(set_frame_width=lambda _fw: None)
    nav.EM_driver = SimpleNamespace(is_microscope_connected=lambda: False)
    nav.im_width = 100
    nav.im_height = 100
    nav._center_pos = (0.0, 0.0)
    nav._sample_hw_mm = 1.0
    nav.grid_search_fw_mm = None
    nav._frame_cntr = 0
    nav.current_frame_label = ""
    nav._visited_band_frames = []
    nav._band_frame_descriptors = []
    nav.frame_pos_mm = []
    nav.frame_labels = []
    nav.num_frames = 0

    parents = [
        ("A0", (0.0, 0.0), 0.10),  # 100 um parent FOV
    ]
    nav.configure_particle_subframes(parents, child_fw_mm=0.05, randomize_frames=False)

    assert nav.num_frames == 4  # 2x2 local grid
    assert all("_" in label for label in nav.frame_labels)
    assert "A0_a0" in nav.frame_labels
    assert all(label.startswith("A0_") for label in nav.frame_labels)
    # Local labels use lowercase letters
    locals_only = [label.split("_", 1)[1] for label in nav.frame_labels]
    assert all(local[0].islower() for local in locals_only)


def test_save_particle_statistics_from_particle_info(outputs_dir: Path):
    out = outputs_dir / "psd_from_particle_info"
    out.mkdir(parents=True, exist_ok=True)

    finder = object.__new__(EM_Particle_Finder)
    finder.results_dir = str(out)
    finder._sample_id = "demo"
    finder.verbose = False
    finder.powder_meas_cfg = PowderMeasurementConfig()
    finder.analyzed_pars = [
        ParticleInfo(
            id=0,
            area_um=float(np.pi),
            eq_diameter_um=2.0,
            coordinates=Coordinate2D(x=1.0, y=2.0),
            frame_id="A0",
        ),
        ParticleInfo(
            id=1,
            area_um=float(np.pi * 4.0),
            eq_diameter_um=4.0,
            coordinates=Coordinate2D(x=1.1, y=2.1),
            frame_id="A0_a1",
        ),
    ]

    stats_df = EM_Particle_Finder.save_particle_statistics(finder)
    assert int(stats_df.loc[0, "n_par_analysed"]) == 2
    assert float(stats_df.loc[0, "min"]) == pytest.approx(2.0)
    assert float(stats_df.loc[0, "max"]) == pytest.approx(4.0)

    sizes_path = out / f"demo_{cnst.PARTICLE_SIZES_FILENAME}.csv"
    assert sizes_path.is_file()
    content = sizes_path.read_text(encoding="utf-8")
    assert "A0_a1" in content
    assert cnst.PAR_EQ_D_KEY in content
