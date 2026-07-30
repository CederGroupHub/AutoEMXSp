#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for par_selection_mode='list' particle targeting."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from autoemx.config.ledger_schemas import (
    Coordinate2D,
    LedgerConfigs,
    ParticleInfo,
    SampleLedger,
)
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder


def _minimal_ledger(sample_dir: Path, particles: list[ParticleInfo]) -> SampleLedger:
    return SampleLedger(
        sample_id="demo",
        sample_path=str(sample_dir.resolve()),
        configs=LedgerConfigs(
            microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"),
            sample_cfg=SampleConfig(elements=["Pb", "Mo", "O"], type="powder"),
            measurement_cfg=MeasurementConfig(),
            sample_substrate_cfg=SampleSubstrateConfig(),
            plot_cfg=PlotConfig(),
        ),
        spectra=[],
        particles=particles,
    )


def _finder_with_list_mode(im_width: int = 1024) -> EM_Particle_Finder:
    powder_cfg = PowderMeasurementConfig(
        par_selection_mode="list",
        max_area_par=300.0,
        min_area_par=10.0,
    )
    moves = []
    frame_widths = []

    em = SimpleNamespace(
        is_initialized=True,
        sample_id="demo",
        im_width=im_width,
        im_height=im_width,
        pixel_size_um=1.0,
        move_to_pos=lambda pos: moves.append(tuple(pos)),
        set_frame_width=lambda fw: frame_widths.append(float(fw)),
        adjust_BCF=lambda: None,
        frame_navigator=SimpleNamespace(_current_pos=(0.0, 0.0)),
    )
    finder = EM_Particle_Finder(em, powder_meas_cfg=powder_cfg, verbose=False)
    finder._moves = moves
    finder._frame_widths = frame_widths
    return finder


def test_powder_config_accepts_list_mode():
    cfg = PowderMeasurementConfig(par_selection_mode="list")
    assert cfg.par_selection_mode == "list"


def test_set_particle_targets_rejects_both_and_neither(outputs_dir: Path):
    finder = _finder_with_list_mode()
    sample_dir = outputs_dir / "list_mode_validate"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ledger = _minimal_ledger(
        sample_dir,
        [ParticleInfo(id=1, coordinates=Coordinate2D(x=1.0, y=2.0))],
    )

    with pytest.raises(ValueError, match="exactly one"):
        finder.set_particle_targets(particle_ids=[1], particle_coords_mm=[(0.0, 0.0)], ledger=ledger)

    with pytest.raises(ValueError, match="exactly one"):
        finder.set_particle_targets(ledger=ledger)


def test_set_particle_targets_missing_id_raises(outputs_dir: Path):
    finder = _finder_with_list_mode()
    sample_dir = outputs_dir / "list_mode_missing_id"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ledger = _minimal_ledger(
        sample_dir,
        [ParticleInfo(id=1, coordinates=Coordinate2D(x=1.0, y=2.0))],
    )
    with pytest.raises(ValueError, match="not found in ledger"):
        finder.set_particle_targets(particle_ids=[1, 99], ledger=ledger)


def test_set_particle_targets_missing_coordinates_raises(outputs_dir: Path):
    finder = _finder_with_list_mode()
    sample_dir = outputs_dir / "list_mode_missing_coords"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ledger = _minimal_ledger(
        sample_dir,
        [ParticleInfo(id=1, area_um=1.0, eq_diameter_um=1.0)],
    )
    with pytest.raises(ValueError, match="missing absolute stage coordinates"):
        finder.set_particle_targets(particle_ids=[1], ledger=ledger)


def test_go_to_next_particle_by_id_uses_existing_particle(outputs_dir: Path):
    finder = _finder_with_list_mode()
    sample_dir = outputs_dir / "list_mode_by_id"
    sample_dir.mkdir(parents=True, exist_ok=True)
    particle = ParticleInfo(
        id=7,
        area_um=12.0,
        eq_diameter_um=4.0,
        coordinates=Coordinate2D(x=3.25, y=-1.5),
        frame_id="A0",
    )
    ledger = _minimal_ledger(sample_dir, [particle])
    finder.set_particle_targets(particle_ids=[7], ledger=ledger)

    assert finder.go_to_next_particle() is True
    assert finder.tot_par_cntr == 7
    assert finder.uses_absolute_particle_ids
    assert finder.current_par_coordinates_mm == pytest.approx((3.25, -1.5))
    assert finder.current_par_area_um2 == pytest.approx(12.0)
    assert finder._moves == [(3.25, -1.5)]
    assert finder.go_to_next_particle() is False


def test_go_to_next_particle_by_coordinates_increments_id():
    finder = _finder_with_list_mode()
    finder.set_particle_targets(particle_coords_mm=[(1.0, 2.0), (3.0, 4.0)])

    assert finder.go_to_next_particle() is True
    assert finder.tot_par_cntr == 1
    assert not finder.uses_absolute_particle_ids
    assert finder.current_par_coordinates_mm == pytest.approx((1.0, 2.0))
    assert finder.current_par_area_um2 is None

    assert finder.go_to_next_particle() is True
    assert finder.tot_par_cntr == 2
    assert finder.current_par_coordinates_mm == pytest.approx((3.0, 4.0))
    assert finder.go_to_next_particle() is False


def test_go_to_next_particle_explicit_coordinate_argument():
    finder = _finder_with_list_mode()
    # Explicit args work even without set_particle_targets when mode is list,
    # and also when coordinates are passed directly.
    assert finder.go_to_next_particle(coordinates=(5.5, 6.5)) is True
    assert finder.tot_par_cntr == 1
    assert finder.current_par_coordinates_mm == pytest.approx((5.5, 6.5))
    assert finder._moves == [(5.5, 6.5)]


def test_go_to_next_particle_rejects_both_args():
    finder = _finder_with_list_mode()
    with pytest.raises(ValueError, match="not both"):
        finder.go_to_next_particle(particle_id=1, coordinates=(0.0, 0.0))
