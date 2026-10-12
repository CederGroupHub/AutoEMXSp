#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Folders where AutoEMX keeps its local files.

AutoEMX writes only to ``~/.autoemx`` (AUTOEMX_HOME), never into the installed package.

The package folder of each microscope model, ``autoemx/calibrations/<microscope_ID>/``, holds its detector model
(detector_model.json) and detector efficiency, read from there, and the initial values of the calibrations of each
instrument.

The calibrations of each instrument (detector channel parameters, SDD calibration sessions, experimental standards
files) are read from and written to ``~/.autoemx/calibrations/<microscope_ID>/``. The first time the local folder is
needed, it is created and filled with a copy of the initial values of the package. After that, a package file is only
copied when the local folder does not have it (e.g. the standards file of a new beam energy). Detector channel
parameters are never copied again, since the most recent file is used: a newer one from the package would replace the
calibration of the local instrument.

Quantities computed from the calibrations (convolution matrices of the detector response) are cached in
``~/.autoemx/cache/<microscope_ID>/``.
"""

import shutil
from pathlib import Path
from typing import List, Optional

import autoemx.utils.constants as cnst

from autoemx._logging import get_logger
logger = get_logger(__name__)

AUTOEMX_HOME = Path.home() / ".autoemx"

PACKAGE_CALIBRATIONS_DIR = Path(__file__).resolve().parent / cnst.CALIBS_DIR


def package_calibration_dir(microscope_ID: str) -> Path:
    """Calibration folder of the microscope in the package (read only)."""
    return PACKAGE_CALIBRATIONS_DIR / microscope_ID


def standards_filename(meas_type: str, beam_energy: float) -> str:
    """Name of the experimental standards file, e.g. EDS_Stds_15keV.json."""
    return f"{meas_type}_{cnst.STD_FILENAME}_{int(beam_energy):d}keV.json"


def _initial_calibration_files(package_dir: Path) -> List[Path]:
    """Package files copied into a new local calibration folder (paths relative to package_dir)."""
    files = [path for path in package_dir.glob(f"*_{cnst.STD_FILENAME}_*keV.json") if path.is_file()]
    detector_params_dir = package_dir / cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR
    if detector_params_dir.is_dir():
        files.extend(
            path for path in detector_params_dir.rglob("*")
            if path.is_file() and path.name != ".DS_Store"
        )
    return sorted(path.relative_to(package_dir) for path in files)


def _describe_copied_files(files: List[Path]) -> str:
    """One line per standards file, and per kind of detector calibration file."""
    sessions_dir = Path(cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR, cnst.SDD_CALIBS_MEAS_DIR)
    standards = [path.name for path in files if len(path.parts) == 1]
    params = [path.name for path in files if path.parent == Path(cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR)]
    sessions = sorted({path.relative_to(sessions_dir).parts[0] for path in files if sessions_dir in path.parents})
    lines = [f"  - Standards file: {name}" for name in standards]
    if params:
        lines.append(f"  - Detector channel parameters: {len(params)} files (latest: {params[-1]})")
    if sessions:
        lines.append(f"  - SDD calibration sessions: {', '.join(sessions)}")
    return "\n".join(lines) or "  (no calibration files in the package)"


def user_calibration_dir(microscope_ID: str) -> Path:
    """
    Local calibration folder of the microscope, ``~/.autoemx/calibrations/<microscope_ID>/``.

    Created on first use and filled with a copy of the calibration files shipped in the package: the experimental
    standards files and the detector channel parameters (with any SDD calibration sessions). The copy is made in a
    temporary folder, renamed when complete, so an interrupted copy is never taken for a complete one.
    """
    user_dir = AUTOEMX_HOME / cnst.CALIBS_DIR / microscope_ID
    if user_dir.is_dir():
        return user_dir

    package_dir = package_calibration_dir(microscope_ID)
    files = _initial_calibration_files(package_dir) if package_dir.is_dir() else []
    tmp_dir = user_dir.with_name(f".{microscope_ID}_incomplete")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    (tmp_dir / cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR).mkdir(parents=True)
    for rel_path in files:
        (tmp_dir / rel_path).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(package_dir / rel_path, tmp_dir / rel_path)
    tmp_dir.rename(user_dir)

    logger.info(
        f"ℹ️ Created the local calibration folder of {microscope_ID} at '{user_dir}', "
        f"with a copy of the calibration files of the package:\n{_describe_copied_files(files)}\n"
        "AutoEMX reads and updates the calibrations of this microscope in this folder from now on."
    )
    return user_dir


def user_calibration_file(microscope_ID: str, filename: str) -> Path:
    """
    Path of a calibration file in the local calibration folder.

    If the local folder does not have the file but the package does, the package file is copied first (an existing
    local file is never replaced). The returned path may not exist, if neither has the file.
    """
    path = user_calibration_dir(microscope_ID) / filename
    package_path = package_calibration_dir(microscope_ID) / filename
    if not path.exists() and package_path.is_file():
        shutil.copy2(package_path, path)
        logger.info(f"ℹ️ Copied '{filename}' of the package to the local calibration folder '{path.parent}'.")
    return path


def user_standards_file(microscope_ID: str, meas_type: str, beam_energy: float) -> Path:
    """Path of the experimental standards file in the local calibration folder (see user_calibration_file)."""
    return user_calibration_file(microscope_ID, standards_filename(meas_type, beam_energy))


def user_cache_dir(microscope_ID: str) -> Path:
    """Cache of quantities computed from the calibrations of the microscope, ``~/.autoemx/cache/<microscope_ID>/``."""
    path = AUTOEMX_HOME / "cache" / microscope_ID
    path.mkdir(parents=True, exist_ok=True)
    return path


def latest_detector_channel_params_file(microscope_ID: str) -> Optional[Path]:
    """Most recent detector channel parameters file in the local calibration folder (None if there is none)."""
    params_dir = user_calibration_dir(microscope_ID) / cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_DIR
    files = sorted(params_dir.glob(f"*_{cnst.DETECTOR_CHANNEL_PARAMS_CALIBR_FILENAME}.json"))
    return files[-1] if files else None
