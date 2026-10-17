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
from autoemx.config.ledger_schemas import (
    FrameRecord,
    FrameSegmentation,
    LedgerConfigs,
    ParticleInfo,
    SampleLedger,
    SegmentationConfig,
    SizeBand,
    SizeMeasurement,
)
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.particle_statistics import recalculate_particle_stats
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


# Particles measured by EDS (older ledgers, without frames searched by EDS): 3 counts (automatic selection), 4 not
LEDGER_PARTICLES = [
    ParticleInfo(id=3, area_um=40.0, frame_id="A2", selection="auto"),
    ParticleInfo(id=4, area_um=50.0, frame_id="A2", selection="manual"),
]
#: Particles of the frames of the measurement (µm²)
FRAME_AREAS = {"A0": [10.0, 20.0], "A1": [30.0], "A2": []}


@pytest.fixture
def ledger_sample(tmp_path: Path) -> Path:
    """Sample whose ledger holds a particle size statistics measurement and its result (no particle sizes file)."""
    sample_dir = tmp_path / SAMPLE_ID
    sample_dir.mkdir()
    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=[], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    ledger = create_sample_ledger(str(sample_dir), SAMPLE_ID, configs)
    ledger.particles = list(LEDGER_PARTICLES)
    cfg = PowderMeasurementConfig(min_area_par=5.0)
    segmentation = ledger.append_segmentation(SegmentationConfig.from_powder_meas_cfg(cfg, 0))
    band = SizeBand(d_min_um=1.0, d_max_um=10.0, frame_width_um=100.0, frame_ids=list(FRAME_AREAS))
    ledger.size_measurements.append(SizeMeasurement(measurement_id=0, status="completed", n_par_target=3,
                                                    powder_meas_cfg=cfg, segmentation_id=segmentation, bands=[band]))
    for frame_id, areas in FRAME_AREAS.items():
        ledger.add_frame(FrameRecord(
            frame_id=frame_id, label=frame_id, purpose="size_stats", measurement_id=0, frame_width_um=100.0,
            pixel_size_um=0.1, area_um2=2e4 / 3, search_window_um=(1.0, 10.0),
            segmentations=[FrameSegmentation(segmentation_id=segmentation, d_range_um=(1.0, 10.0),
                                             particle_areas_um2=areas, particle_edge_corrections=[1.0] * len(areas))],
        ))
    ledger.append_size_result(recalculate_particle_stats(ledger))
    ledger.to_json_file(str(sample_dir / "ledger.json"))
    return tmp_path


def _ledger(results_path: Path) -> SampleLedger:
    return SampleLedger.from_json_file(str(results_path / SAMPLE_ID / "ledger.json"))


def test_ledger_statistics_recalculated_and_exclusions_recorded(ledger_sample: Path):
    original = _ledger(ledger_sample).size_result
    assert original.n_par_by_selection == {"size_stats": 3, "auto": 1}
    stats = process_particle_statistics(SAMPLE_ID, str(ledger_sample), particles_to_filter=[("A0", 1), ("A9", 0)],
                                        frame_ids_to_filter=["A2"], verbose=False)
    assert stats["n_par_analysed"].iloc[0] == 2  # A0 particle 0 and A1 particle 0 (EDS particle 3 is in A2)

    ledger = _ledger(ledger_sample)
    result = ledger.size_result
    assert result.result_id == 1 and len(ledger.size_results) == 2  # the statistics of the measurement are kept
    assert result.excluded_particles == [("A0", 1)] and result.excluded_frame_ids == ["A2"]  # A9 was not found
    assert result.bands[0].frame_ids == original.bands[0].frame_ids == ["A0", "A1", "A2"]
    assert result.bands[0].number_density_per_mm2 == pytest.approx(2 / (2e4 * 2 / 3) * 1e6)  # A2 area left out
    assert result.n_par_target == 3 and result.powder_meas_cfg.min_area_par == 5.0
    sizes = pd.read_csv(ledger_sample / SAMPLE_ID / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}_processed.csv",
                        dtype={cnst.FRAME_ID_DF_KEY: str})
    assert list(zip(sizes[cnst.FRAME_ID_DF_KEY], sizes[cnst.PAR_ID_DF_KEY])) == [("A0", 0), ("A1", 0)]

    # Recalculating without exclusions gives back the statistics of the measurement (same result, now active)
    process_particle_statistics(SAMPLE_ID, str(ledger_sample), verbose=False)
    ledger = _ledger(ledger_sample)
    random_fields = {"datetime", "confidence_intervals_95"}  # bootstrap intervals vary slightly from run to run
    assert ledger.size_result.model_dump(exclude=random_fields) == original.model_dump(exclude=random_fields)
    assert len(ledger.size_results) == 2


def test_eds_auto_particles_follow_the_measurement_setting_unless_given(ledger_sample: Path):
    process_particle_statistics(SAMPLE_ID, str(ledger_sample), verbose=False)
    assert _ledger(ledger_sample).size_result.n_par_by_selection == {"size_stats": 3, "auto": 1}  # 4: manual

    process_particle_statistics(SAMPLE_ID, str(ledger_sample), include_eds_auto_particles=False, verbose=False)
    ledger = _ledger(ledger_sample)
    assert ledger.size_result.n_par_by_selection == {"size_stats": 3}
    assert ledger.size_result.powder_meas_cfg.include_eds_auto_particles is False
    assert len(ledger.size_results) == 2  # with and without what EDS measured
