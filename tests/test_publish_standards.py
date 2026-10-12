#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Publishing the local standards file of an instrument to the package (on a copy of the package folder)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import autoemx._paths as paths
from autoemx.runners.publish_standards import publish_standards

MICROSCOPE_ID = "PhenomXL"
STDS_FILENAME = "EDS_Stds_15keV.json"


@pytest.fixture
def package_copy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Copy of the package calibrations in a git checkout, used instead of the package."""
    calibrations = tmp_path / "repo" / "autoemx" / "calibrations"
    (tmp_path / "repo" / ".git").mkdir(parents=True)
    (calibrations / MICROSCOPE_ID).mkdir(parents=True)
    shutil.copy2(paths.package_calibration_dir(MICROSCOPE_ID) / STDS_FILENAME, calibrations / MICROSCOPE_ID / STDS_FILENAME)
    paths.user_calibration_dir(MICROSCOPE_ID)  # local folder filled from the real package
    monkeypatch.setattr(paths, "PACKAGE_CALIBRATIONS_DIR", calibrations)
    return calibrations / MICROSCOPE_ID


def test_local_standards_published_with_report(package_copy: Path):
    local_file = paths.user_calibration_dir(MICROSCOPE_ID) / STDS_FILENAME
    local_file.write_text(local_file.read_text() + "\n")  # updated locally

    published = publish_standards(MICROSCOPE_ID, 15)

    assert published == package_copy / STDS_FILENAME
    assert published.read_bytes() == local_file.read_bytes()
    assert (package_copy / "experimental_standards_details_15keV.txt").exists()


def test_identical_standards_not_published_again(package_copy: Path):
    assert publish_standards(MICROSCOPE_ID, 15) is None
    assert not (package_copy / "experimental_standards_details_15keV.txt").exists()


def test_not_published_outside_git_checkout(package_copy: Path):
    shutil.rmtree(package_copy.parents[2] / ".git")
    with pytest.raises(RuntimeError, match="git checkout"):
        publish_standards(MICROSCOPE_ID, 15)


def test_missing_local_standards(package_copy: Path):
    with pytest.raises(FileNotFoundError, match="No local standards file"):
        publish_standards(MICROSCOPE_ID, 20)
