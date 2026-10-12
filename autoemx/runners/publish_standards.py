#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: Publish Experimental Standards to the Package

AutoEMX reads and updates the experimental standards of an instrument in its local calibration folder,
~/.autoemx/calibrations/<microscope_ID>/. The standards file in the package (autoemx/calibrations/<microscope_ID>/)
only holds the initial values copied to new local folders, and AutoEMX never changes it.

This runner copies the local standards file of a beam energy into the package, and regenerates its details report
next to it, so the new values can be committed and distributed as the initial values of the microscope.
It only runs in a git checkout of AutoEMX, where the package files can be committed.

Created on Mon Oct 12 2026
@author: Andrea
"""

import filecmp
import shutil
from pathlib import Path
from typing import Optional

import autoemx._paths as paths


def _git_checkout_root(path: Path) -> Optional[Path]:
    for parent in (path, *path.parents):
        if (parent / ".git").exists():
            return parent
    return None


def publish_standards(microscope_ID: str, beam_energy: float, meas_type: str = "EDS") -> Optional[Path]:
    """
    Copy the local standards file of a beam energy into the package, and regenerate its details report there.

    Parameters
    ----------
    microscope_ID : str
        Microscope whose standards are published (folder of autoemx/calibrations).
    beam_energy : float
        Beam energy (keV) of the standards file, e.g. 15 for EDS_Stds_15keV.json.
    meas_type : str, optional
        Measurement type of the standards file. Default is "EDS".

    Returns
    -------
    Path or None
        Path of the standards file in the package, or None if it was already identical to the local one.

    Raises
    ------
    RuntimeError
        If AutoEMX is not a git checkout (an installed package), where published files could not be committed.
    FileNotFoundError
        If the microscope has no folder in the package, or no local standards file for this beam energy.
    """
    package_dir = paths.package_calibration_dir(microscope_ID)
    if _git_checkout_root(paths.PACKAGE_CALIBRATIONS_DIR) is None:
        raise RuntimeError(
            f"AutoEMX at '{paths.PACKAGE_CALIBRATIONS_DIR.parent}' is not a git checkout: standards can only be "
            "published from a git checkout of AutoEMX, where they can be committed."
        )
    if not package_dir.is_dir():
        raise FileNotFoundError(f"No folder for microscope '{microscope_ID}' in the package: '{package_dir}'.")
    filename = paths.standards_filename(meas_type, beam_energy)
    local_path = paths.AUTOEMX_HOME / "calibrations" / microscope_ID / filename
    if not local_path.is_file():
        raise FileNotFoundError(f"No local standards file to publish: '{local_path}'.")

    package_path = package_dir / filename
    if package_path.is_file() and filecmp.cmp(local_path, package_path, shallow=False):
        print(f"The standards file in the package is already identical to the local one: {package_path}")
        return None
    shutil.copy2(local_path, package_path)
    print(f"Published '{local_path}'\n       to '{package_path}'.")

    from autoemx.runners.extract_experimental_standards_details import extract_experimental_standards_details
    extract_experimental_standards_details(
        microscope_ID=microscope_ID,
        voltage=beam_energy,
        standards_json_path=str(package_path),
        print_report=False,
    )
    print(f"Regenerated the standards details report in '{package_dir}'.\nCommit both files to distribute them.")
    return package_path
