#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Editing the standards file: deletions, exclusions from the means, backups, and the shared reference-mean calculation."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from autoemx.config.schema_models.standards import StandardLine
from autoemx.runners.modify_standards_file import modify_standards_file
import autoemx.utils.constants as cnst

STANDARDS_FILE = Path(__file__).resolve().parents[1] / "autoemx" / "calibrations" / "PhenomXL" / "EDS_Stds_15keV.json"


@pytest.fixture
def stds_file(tmp_path: Path) -> Path:
    path = tmp_path / STANDARDS_FILE.name
    shutil.copy(STANDARDS_FILE, path)
    return path


def _point_lines(path: Path) -> dict:
    return json.loads(path.read_text())["standards_by_mode"]["point"]


def _backups(path: Path) -> list:
    return sorted(path.parent.glob(f"{path.stem}_backup_*{path.suffix}"))


def test_no_change_requested_leaves_file_untouched(stds_file: Path):
    original = stds_file.read_bytes()
    summary = modify_standards_file(stds_file)
    assert summary["changed"] is False
    assert stds_file.read_bytes() == original
    assert _backups(stds_file) == []


def test_exclusion_recomputes_only_the_modified_line(stds_file: Path):
    # B_Ka1 has two entries used for the mean: LaB6 (5682.6) and BN (2123.1)
    original = _point_lines(stds_file)
    summary = modify_standards_file(stds_file, exclude_std_lines_from_mean_calc=[{"LaB6": ["B_Ka1"]}])

    lines = _point_lines(stds_file)
    assert summary["excluded_entries"] == [("LaB6", "B_Ka1")]
    assert summary["means_recomputed"] == ["B_Ka1"]
    lab6 = next(e for e in lines["B_Ka1"]["entries"] if e[cnst.STD_ID_KEY] == "LaB6")
    assert lab6[cnst.STD_USE_FOR_MEAN_KEY] is False
    assert lines["B_Ka1"]["reference_mean"][cnst.COR_PB_DF_KEY] == pytest.approx(2123.1)
    assert {k: v for k, v in lines.items() if k != "B_Ka1"} == {k: v for k, v in original.items() if k != "B_Ka1"}
    [backup] = _backups(stds_file)
    assert json.loads(backup.read_text())["standards_by_mode"]["point"] == original


def test_deleting_entries_and_empty_lines(stds_file: Path):
    # In_Mz1 has two entries (InP, In): deleting both removes the line
    summary = modify_standards_file(stds_file, delete_std_lines=[{"InP": ["In_Mz1"]}, {"In": ["In_Mz1"]}])
    assert sorted(summary["deleted_entries"]) == [("In", "In_Mz1"), ("InP", "In_Mz1")]
    assert summary["deleted_el_lines"] == ["In_Mz1"]
    assert "In_Mz1" not in _point_lines(stds_file)


def test_delete_el_lines_from_all_standards(stds_file: Path):
    assert "B_Ka1" in _point_lines(stds_file)
    summary = modify_standards_file(stds_file, delete_el_lines=["B_Ka1", "Xx_Ka1"])
    assert summary["deleted_el_lines"] == ["B_Ka1"]
    assert summary["not_found_el_lines"] == ["Xx_Ka1"]
    assert "B_Ka1" not in _point_lines(stds_file)


def test_without_backup(stds_file: Path):
    modify_standards_file(stds_file, delete_el_lines=["B_Ka1"], backup=False)
    assert _backups(stds_file) == []


def test_unknown_mode(stds_file: Path):
    with pytest.raises(ValueError, match="Available modes"):
        modify_standards_file(stds_file, delete_el_lines=["B_Ka1"], mode="line_scan")


def test_reference_mean_uses_only_entries_used_for_the_mean():
    def entry(standard_id, corrected_pb, use_for_mean_calc=None):
        payload = {cnst.STD_ID_KEY: standard_id, cnst.DATETIME_KEY: "2026-01-01", cnst.COR_PB_DF_KEY: corrected_pb,
                   cnst.STDEV_PB_DF_KEY: 1.0, cnst.REL_ER_PERCENT_PB_DF_KEY: 1.0}
        if use_for_mean_calc is not None:
            payload[cnst.STD_USE_FOR_MEAN_KEY] = use_for_mean_calc
        return payload

    line = StandardLine.model_validate({"entries": [
        entry("A", 100.0, True), entry("B", 200.0, True), entry("C", 900.0, False), entry("D", 900.0),  # no flag: not used
    ]})
    mean = line.compute_reference_mean("2026-10-11")
    assert mean.corrected_pb == pytest.approx(150.0)
    assert mean.stdev_pb == pytest.approx(np.std([100.0, 200.0]))
    assert mean.datetime == "2026-10-11"
    assert StandardLine.model_validate({"entries": [entry("C", 900.0, False)]}).compute_reference_mean("x") is None
