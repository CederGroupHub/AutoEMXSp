#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibrations are read from and written to the local folder ~/.autoemx/calibrations/<ID>/, never the package.

The local AutoEMX folder is a temporary folder in tests (autoemx_home fixture in conftest.py).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import autoemx._paths as paths
import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst

MICROSCOPE_ID = "PhenomXL"
STDS_FILENAME = "EDS_Stds_15keV.json"
PARAMS_DIR = cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR


def _tree_hashes(folder: Path) -> dict:
    return {
        path.relative_to(folder).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(folder.rglob("*")) if path.is_file() and "__pycache__" not in path.parts
    }


@pytest.fixture
def package_unchanged():
    """Fails the test if any file of the package calibration folder is added, removed or changed."""
    package_dir = paths.package_calibration_dir(MICROSCOPE_ID)
    before = _tree_hashes(package_dir)
    yield
    assert _tree_hashes(package_dir) == before


def _write_params(user_dir: Path, name: str, offset: float) -> Path:
    path = user_dir / PARAMS_DIR / f"{name}_{cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_FILENAME}.json"
    path.write_text(json.dumps({"point": {cnst.BEAM_CURRENT_KEY: 1.0, cnst.SCALE_KEY: 0.01, cnst.OFFSET_KEY: offset}}))
    return path


def test_local_folder_created_with_copy_of_package_files(autoemx_home: Path, package_unchanged):
    package_dir = paths.package_calibration_dir(MICROSCOPE_ID)
    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)

    assert user_dir == autoemx_home / "calibrations" / MICROSCOPE_ID
    assert (user_dir / STDS_FILENAME).read_bytes() == (package_dir / STDS_FILENAME).read_bytes()
    package_params = sorted(p.name for p in (package_dir / PARAMS_DIR).glob("*.json"))
    assert package_params
    assert sorted(p.name for p in (user_dir / PARAMS_DIR).glob("*.json")) == package_params
    # The detector model of the microscope model stays in the package
    for name in (cnst.DETECTOR_MODEL_FILENAME, cnst.DETECTOR_EFFICIENCY_FILENAME):
        assert not (user_dir / name).exists()
    assert not list(autoemx_home.glob("calibrations/.*_incomplete"))


def test_local_folder_is_filled_only_once(package_unchanged):
    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)
    (user_dir / STDS_FILENAME).write_text("local values")
    (user_dir / PARAMS_DIR / sorted((user_dir / PARAMS_DIR).glob("*.json"))[0].name).unlink()
    n_params = len(list((user_dir / PARAMS_DIR).glob("*.json")))

    assert paths.user_calibration_dir(MICROSCOPE_ID) == user_dir
    assert (user_dir / STDS_FILENAME).read_text() == "local values"
    assert len(list((user_dir / PARAMS_DIR).glob("*.json"))) == n_params  # deleted file not copied again


def test_interrupted_first_fill_is_redone(autoemx_home: Path):
    incomplete = autoemx_home / "calibrations" / f".{MICROSCOPE_ID}_incomplete"
    incomplete.mkdir(parents=True)
    (incomplete / "partial.json").write_text("{}")

    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)
    assert (user_dir / STDS_FILENAME).exists()
    assert not (user_dir / "partial.json").exists()
    assert not incomplete.exists()


def test_standards_file_of_new_beam_energy_copied_from_package(package_unchanged):
    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)
    (user_dir / STDS_FILENAME).unlink()
    path = paths.user_standards_file(MICROSCOPE_ID, "EDS", 15)
    assert path == user_dir / STDS_FILENAME
    assert path.read_bytes() == (paths.package_calibration_dir(MICROSCOPE_ID) / STDS_FILENAME).read_bytes()

    # An existing local file is never replaced
    path.write_text("local values")
    paths.user_standards_file(MICROSCOPE_ID, "EDS", 15)
    assert path.read_text() == "local values"

    # Beam energy without standards, in the package or locally
    assert not paths.user_standards_file(MICROSCOPE_ID, "EDS", 99).exists()


def test_detector_channel_params_read_from_local_folder(package_unchanged):
    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)
    _write_params(user_dir, "29991231_23h59m", offset=-0.123)

    calibs.load_microscope_calibrations(MICROSCOPE_ID, "point", load_detector_channel_params=True)
    assert calibs.detector_channel_params["point"][cnst.OFFSET_KEY] == -0.123
    assert Path(calibs.calibration_files_dir) == user_dir / PARAMS_DIR  # also where SDD sessions are saved


def test_new_detector_calibration_written_to_local_folder(package_unchanged):
    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)
    for path in (user_dir / PARAMS_DIR).glob("*.json"):
        path.unlink()
    _write_params(user_dir, "20000101_00h00m", offset=-0.1)
    calibs.load_microscope_calibrations(MICROSCOPE_ID, "point", load_detector_channel_params=True)
    loaded_params = calibs.detector_channel_params

    calibs.update_detector_channel_params("point", new_offset=-0.2, new_scale=0.02, verbose=False)

    latest = paths.latest_detector_channel_params_file(MICROSCOPE_ID)
    assert latest.parent == user_dir / PARAMS_DIR
    saved = json.loads(latest.read_text())["point"]
    assert (saved[cnst.OFFSET_KEY], saved[cnst.SCALE_KEY]) == (-0.2, 0.02)
    assert loaded_params["point"][cnst.OFFSET_KEY] == -0.1  # parameters loaded before are not changed


def test_standards_read_from_local_folder(package_unchanged):
    calibs.load_microscope_calibrations(MICROSCOPE_ID, "point")
    user_file = paths.user_standards_file(MICROSCOPE_ID, "EDS", 15)
    payload = json.loads(user_file.read_text())
    payload["standards_by_mode"]["point"] = {}
    user_file.write_text(json.dumps(payload))

    standards, standards_path = calibs.load_standards_model("EDS", 15)
    assert Path(standards_path) == user_file
    assert standards.standards_by_mode["point"] == {}


def test_load_standards_requires_loaded_calibrations(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(calibs, "loaded_microscope_ID", None)
    with pytest.raises(RuntimeError, match="load_microscope_calibrations"):
        calibs.load_standards_model("EDS", 15)


def test_gui_lists_local_and_package_standards(package_unchanged):
    from autoemx.gui import backend

    user_dir = paths.user_calibration_dir(MICROSCOPE_ID)
    payload = json.loads((user_dir / STDS_FILENAME).read_text())
    payload["beam_energy_keV"] = 20
    payload["standards_by_mode"] = {"point": {"Fe_Ka1": payload["standards_by_mode"]["point"]["Fe_Ka1"]}}
    (user_dir / "EDS_Stds_20keV.json").write_text(json.dumps(payload))
    (user_dir / STDS_FILENAME).unlink()  # still in the package: copied when first used

    assert backend.standards_beam_energies(MICROSCOPE_ID) == [15.0, 20.0]
    assert backend.quantifiable_elements(MICROSCOPE_ID, 20) == ["Fe"]
    assert backend.quantifiable_elements(MICROSCOPE_ID, 15) is not None
    assert (user_dir / STDS_FILENAME).exists()
