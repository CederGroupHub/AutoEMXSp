#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ledgers saved before frames and segmentations were recorded (``particle_stats``, ``particle_stats_progress``,
``eds_size_scans``, made with AutoEMX 1.6.3: tests/inputs/legacy_size_stats) converted when loaded.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import SampleLedger
from autoemx.core.composition_analysis.ledger import LedgerModule
from autoemx.core.particle_statistics import recalculate_particle_stats
from ci_paths import INPUTS_DIR

LEGACY_DIR = INPUTS_DIR / "legacy_size_stats"
STATS_FIELDS = ("n_par_analysed", "mean", "stdev", "median", "min", "max", "D10", "D25", "D75", "D90", "D43", "Dv10",
                "Dv50", "Dv90")


def _legacy(name: str, tmp_path: Path) -> Path:
    """Copy of a legacy sample folder (loading a ledger rewrites it)."""
    sample_dir = tmp_path / name
    shutil.copytree(LEGACY_DIR / name, sample_dir)
    return sample_dir


def _migrated(name: str, tmp_path: Path) -> SampleLedger:
    return SampleLedger.model_validate(json.loads((_legacy(name, tmp_path) / "ledger.json").read_text()))


def test_measurement_with_eds_frames_migrated(tmp_path: Path):
    raw = json.loads((LEGACY_DIR / "measured_with_eds" / "ledger.json").read_text())
    ledger = _migrated("measured_with_eds", tmp_path)

    # One segmentation: same settings for the measurement and the EDS acquisition
    [segmentation] = ledger.segmentations
    assert (segmentation.segmentation_id, segmentation.par_brightness_thresh, segmentation.min_par_pixels) == (0, 100, 80)
    # Frames: those of the 2 size bands (6 each), then those searched by EDS, with their particles
    band_frames, eds_frames = ledger.frames_of("size_stats", 0), ledger.frames_of("eds_search")
    assert len(band_frames) == 12 and len(eds_frames) == 6
    assert [f.frame_id for f in eds_frames] == [f"EDS-B{i}" for i in range(6)] and eds_frames[0].label == "B0"
    legacy_frame = raw["particle_stats"]["bands"][0]["frames"][0]
    frame = ledger.frame(legacy_frame["frame_id"])
    assert frame.image_file == legacy_frame["image_file"] and frame.center_mm.x == legacy_frame["center_mm"]["x"]
    assert frame.segmentation(0).n_particles == legacy_frame["n_particles"]
    assert len(frame.segmentation(0).particle_marks_px) == legacy_frame["n_particles"]
    # Particles of the size bands moved to their frames; the particles measured by EDS stay
    assert [p.selection for p in ledger.particles] == ["auto"] * 3

    [measurement] = ledger.size_measurements
    assert (measurement.status, measurement.n_par_target, measurement.segmentation_id) == ("completed", 12, 0)
    assert [len(b.frame_ids) for b in measurement.bands] == [6, 6]
    result = ledger.size_result
    assert (result.result_id, result.measurement_id, result.segmentation_id, result.eds_segmentation_id) == (0, 0, 0, 0)
    assert result.excluded_frame_ids == ["A3"] and result.excluded_particles == [("A0", 0)]
    assert result.size_stats.model_dump() == raw["particle_stats"]["size_stats"]  # statistics kept

    # Recomputed from the migrated frames: the same statistics
    recomputed = recalculate_particle_stats(ledger)
    for key in STATS_FIELDS:
        assert getattr(recomputed.size_stats, key) == pytest.approx(getattr(result.size_stats, key)), key
    assert recomputed.n_par_by_selection == result.n_par_by_selection == raw["particle_stats"]["n_par_by_selection"]
    assert recomputed.inputs_key() == result.inputs_key()


def test_eds_only_sample_migrated(tmp_path: Path):
    ledger = _migrated("eds_only", tmp_path)
    assert len(ledger.frames_of("eds_search")) == 6 and ledger.size_measurements == []
    result = ledger.size_result
    assert (result.measurement_id, result.segmentation_id, result.eds_segmentation_id) == (None, None, 0)
    recomputed = recalculate_particle_stats(ledger)
    for key in STATS_FIELDS:
        assert getattr(recomputed.size_stats, key) == pytest.approx(getattr(result.size_stats, key)), key


def test_measurement_in_progress_migrated(tmp_path: Path):
    ledger = _migrated("in_progress", tmp_path)
    [measurement] = ledger.size_measurements
    assert (measurement.status, measurement.n_par_target) == ("in_progress", 50)
    assert ledger.segmentation(measurement.segmentation_id).par_brightness_thresh == 120  # its own settings
    assert len(measurement.bands) == 1 and len(ledger.frames) == 6 and ledger.size_results == []
    assert ledger.particles == []


def test_migrated_ledger_saved_in_the_new_format(tmp_path: Path):
    sample_dir = _legacy("measured_with_eds", tmp_path)
    ledger = load_sample_ledger(str(sample_dir / "ledger.json"))  # rewrites the file (sample path changed)
    ledger.to_json_file(str(sample_dir / "ledger.json"))
    saved = json.loads((sample_dir / "ledger.json").read_text())
    assert not {"particle_stats", "particle_stats_progress", "eds_size_scans"} & set(saved)
    assert len(saved["frames"]) == 18 and saved["active_size_result"] == 0
    assert SampleLedger.from_json_file(str(sample_dir / "ledger.json")) == ledger


def test_unreadable_ledger_never_replaced(tmp_path: Path):
    (tmp_path / "ledger.json").write_text('{"sample_id": "S", "broken": ')
    analyzer = SimpleNamespace(sample_result_dir=str(tmp_path))
    with pytest.raises(ValueError):
        LedgerModule._load_existing_ledger(analyzer)
    assert (tmp_path / "ledger.json").read_text() == '{"sample_id": "S", "broken": '
    (tmp_path / "ledger.json").write_text("")  # empty: as if there were none
    assert LedgerModule._load_existing_ledger(analyzer) is None
