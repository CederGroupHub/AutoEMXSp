#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recalculating particle size statistics without some particles or frames (Process_Particle_Stats_Files.py)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import create_sample_ledger
from autoemx.config.ledger_schemas import LedgerConfigs, ParticleInfo, ParticleSizeBand, SampleLedger
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.particle_statistics import build_particle_stats_result, select_size_stats_particles
from autoemx.runners.process_particle_statistics import process_particle_statistics

SAMPLE_ID = "PSD_sample"
AREAS = {0: ("fr1", 3.0), 1: ("fr1", 12.0), 2: ("fr2", 7.0), 3: ("fr2", 50.0), 4: ("fr3", 20.0)}


@pytest.fixture
def results_path(tmp_path: Path) -> Path:
    sample_dir = tmp_path / SAMPLE_ID
    sample_dir.mkdir()
    pd.DataFrame({
        cnst.PAR_ID_DF_KEY: list(AREAS),
        cnst.FRAME_ID_DF_KEY: [frame for frame, _ in AREAS.values()],
        cnst.PAR_AREA_UM_KEY: [area for _, area in AREAS.values()],
        cnst.PAR_EQ_D_KEY: [2 * np.sqrt(area / np.pi) for _, area in AREAS.values()],
    }).to_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv", index=False)
    return tmp_path


def _diameters(*ids) -> np.ndarray:
    return np.array([2 * np.sqrt(AREAS[i][1] / np.pi) for i in ids])


def test_statistics_recalculated_without_filtered_particles_and_frames(results_path: Path):
    stats = process_particle_statistics(SAMPLE_ID, str(results_path), particle_ids_to_filter=[1],
                                        frame_ids_to_filter=["fr3"], verbose=False)
    kept = _diameters(0, 2, 3)
    assert stats["n_par_analysed"].iloc[0] == 3
    assert stats["mean"].iloc[0] == pytest.approx(kept.mean())
    assert stats["median"].iloc[0] == pytest.approx(np.median(kept))

    sample_dir = results_path / SAMPLE_ID
    processed = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}_processed.csv")
    assert list(processed[cnst.PAR_ID_DF_KEY]) == [0, 2, 3]
    assert (sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}_processed.csv").exists()
    assert list(sample_dir.glob("*hist*_processed.png"))


def test_ids_not_found_are_reported_and_nothing_left_returns_none(results_path: Path, caplog):
    stats = process_particle_statistics(SAMPLE_ID, str(results_path), particle_ids_to_filter=[99], verbose=False)
    assert stats["n_par_analysed"].iloc[0] == len(AREAS)
    assert "particle IDs [99]" in caplog.text
    assert process_particle_statistics(SAMPLE_ID, str(results_path), frame_ids_to_filter=["fr1", "fr2", "fr3"],
                                       verbose=False) is None


def test_missing_sample(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        process_particle_statistics("Nope", str(tmp_path))


LEDGER_PARTICLES = [
    ParticleInfo(id=0, area_um=10.0, frame_id="A0", selection="size_stats"),
    ParticleInfo(id=1, area_um=20.0, frame_id="A0", selection="size_stats"),
    ParticleInfo(id=2, area_um=30.0, frame_id="A1", selection="size_stats"),
    ParticleInfo(id=3, area_um=40.0, frame_id="A2", selection="auto"),
    ParticleInfo(id=4, area_um=50.0, frame_id="A2", selection="manual"),
]


@pytest.fixture
def ledger_sample(tmp_path: Path) -> Path:
    """Sample whose ledger holds a particle size statistics result (and no particle sizes file)."""
    sample_dir = tmp_path / SAMPLE_ID
    sample_dir.mkdir()
    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=[], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    ledger = create_sample_ledger(str(sample_dir), SAMPLE_ID, configs)
    ledger.particles = list(LEDGER_PARTICLES)
    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, frame_width_um=100.0, n_frames=3, scanned_area_um2=2e4,
                            particle_ids=[0, 1, 2])
    ledger.particle_stats = build_particle_stats_result(
        select_size_stats_particles(LEDGER_PARTICLES), PowderMeasurementConfig(min_area_par=5.0), n_par_target=3,
        bands=[band],
    )
    ledger.to_json_file(str(sample_dir / "ledger.json"))
    return tmp_path


def _ledger_result(results_path: Path):
    return SampleLedger.from_json_file(str(results_path / SAMPLE_ID / "ledger.json")).particle_stats


def test_ledger_statistics_recalculated_and_exclusions_recorded(ledger_sample: Path):
    original = _ledger_result(ledger_sample)
    stats = process_particle_statistics(SAMPLE_ID, str(ledger_sample), particle_ids_to_filter=[1, 99],
                                        frame_ids_to_filter=["A2"], verbose=False)
    assert stats["n_par_analysed"].iloc[0] == 2

    result = _ledger_result(ledger_sample)
    assert result.particle_ids == [0, 2]
    assert result.excluded_particle_ids == [1] and result.excluded_frame_ids == ["A2"]  # 99 was not found
    measurement_fields = {"d_min_um", "d_max_um", "frame_width_um", "n_frames", "scanned_area_um2", "particle_ids"}
    assert result.bands[0].model_dump(include=measurement_fields) == original.bands[0].model_dump(
        include=measurement_fields)  # the bands describe the measurement
    assert result.bands[0].number_density_per_mm2 == pytest.approx(2 / 2e4 * 1e6)  # 2 particles left in 20000 µm²
    assert result.n_par_target == 3 and result.powder_meas_cfg.min_area_par == 5.0
    sizes = pd.read_csv(ledger_sample / SAMPLE_ID / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}_processed.csv")
    assert list(sizes[cnst.PAR_ID_DF_KEY]) == [0, 2]

    # All particles stay in the ledger: recalculating without exclusions restores the statistics of the measurement
    process_particle_statistics(SAMPLE_ID, str(ledger_sample), verbose=False)
    restored = _ledger_result(ledger_sample)
    random_fields = {"datetime", "confidence_intervals_95"}  # bootstrap intervals vary slightly from run to run
    assert restored.model_dump(exclude=random_fields) == original.model_dump(exclude=random_fields)


def test_eds_auto_particles_follow_the_measurement_setting_unless_given(ledger_sample: Path):
    process_particle_statistics(SAMPLE_ID, str(ledger_sample), verbose=False)
    assert _ledger_result(ledger_sample).particle_ids == [0, 1, 2, 3]  # manual particle 4 never included

    process_particle_statistics(SAMPLE_ID, str(ledger_sample), include_eds_auto_particles=False, verbose=False)
    result = _ledger_result(ledger_sample)
    assert result.particle_ids == [0, 1, 2]
    assert result.powder_meas_cfg.include_eds_auto_particles is False
