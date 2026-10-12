#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Settings of experimental standards given to batch_fit_experimental_stds replace those in the ledger for that run."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from autoemx.runners.batch_fit_experimental_stds import batch_fit_experimental_stds
import autoemx.utils.constants as cnst

# Al standard of an SDD calibration (5 spectra, already fitted)
AL_STANDARD = (
    Path(__file__).resolve().parents[1] / "autoemx" / "calibrations" / "PhenomXL" / "Detector_channel_params_calibs"
    / "SDD calibrations" / "20260521_09h22m" / "Al"
)


@pytest.fixture
def stds_path(tmp_path: Path) -> Path:
    shutil.copytree(AL_STANDARD, tmp_path / "Al")
    return tmp_path


SETTINGS = ("els_to_use_for_mean_PB_calc", "min_acceptable_PB_ratio", "quant_flags_accepted", "generate_separate_std_dict")


def _saved_settings(stds_path: Path) -> dict:
    ledger = json.loads((stds_path / "Al" / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)).read_text())
    return {key: ledger["configs"]["measurement_cfg"]["exp_stds_cfg"][key] for key in SETTINGS}


def test_settings_given_replace_saved_ones_for_this_run(stds_path: Path):
    saved_before = _saved_settings(stds_path)
    [analyzer] = batch_fit_experimental_stds(
        std_ids=["Al"], stds_path=str(stds_path), update_std_library=False,
        els_to_use_for_mean_PB_calc=["none"], min_acceptable_PB_ratio=50, quant_flags_accepted=[0, -1],
        generate_separate_std_dict=True,
    )
    cfg = analyzer.exp_stds_cfg
    assert cfg.els_to_use_for_mean_PB_calc is None  # ["none"] normalized as when loaded from a ledger
    assert cfg.min_acceptable_PB_ratio == 50
    assert cfg.quant_flags_accepted == [0, -1]
    assert cfg.generate_separate_std_dict is True
    assert _saved_settings(stds_path) == saved_before  # the ledger keeps its own settings


def test_saved_settings_used_by_default(stds_path: Path):
    saved = _saved_settings(stds_path)
    [analyzer] = batch_fit_experimental_stds(std_ids=["Al"], stds_path=str(stds_path), update_std_library=False)
    assert analyzer.exp_stds_cfg.min_acceptable_PB_ratio == saved["min_acceptable_PB_ratio"]
    assert analyzer.exp_stds_cfg.quant_flags_accepted == saved["quant_flags_accepted"]


def test_invalid_setting_rejected(stds_path: Path):
    # Invalid values are rejected as when loaded from a ledger
    with pytest.raises(ValueError, match="not a recognized element"):
        batch_fit_experimental_stds(std_ids=["Al"], stds_path=str(stds_path), update_std_library=False,
                                    els_to_use_for_mean_PB_calc=["Xx"])


@pytest.mark.parametrize("min_pb_ratio", [10, 7000, 1e5])
def test_min_acceptable_pb_ratio_excludes_weak_peaks(stds_path: Path, min_pb_ratio):
    from autoemx.core.composition_analysis.standards import StandardsModule

    [analyzer] = batch_fit_experimental_stds(std_ids=["Al"], stds_path=str(stds_path), update_std_library=False,
                                             min_acceptable_PB_ratio=min_pb_ratio)
    all_pb_ratios = [
        peak.pb_ratio for record in analyzer.spectra_quant_records if record is not None and record.quant_flag == 0
        for key, peak in record.fit_result.fitted_peaks.items() if key == "Al_Ka1"
    ]
    lines = StandardsModule._assemble_std_PB_data_from_records(analyzer)
    expected = sorted(pb for pb in all_pb_ratios if pb >= min_pb_ratio)
    if expected:
        assert sorted(lines["Al_Ka1"].pb_ratios) == pytest.approx(expected)
    else:
        assert "Al_Ka1" not in lines  # no peak left: the line is not used
    assert (len(expected) == len(all_pb_ratios)) == (min_pb_ratio == 10)  # 7000 excludes some spectra
