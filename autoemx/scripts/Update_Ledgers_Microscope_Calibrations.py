#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update Ledgers with Latest Microscope Calibrations

Use this script to update the microscope_cfg.energy_zero and
microscope_cfg.bin_width fields in sample ledgers with the latest detector
channel calibration values loaded from the calibration files.

For each sample, the microscope ID and measurement mode are read directly from
the ledger, so calibration files are resolved automatically.

Created on Sat May 23 2026

@author: Andrea
"""

from pathlib import Path

# =============================================================================
# Input
# =============================================================================

# Directory containing the sample subfolders.
directory = '/Users/Andrea_1/Desktop/Work/Projects/SEM EDX automation/EDX standards/Auto measurements'

# Optional list of sample IDs (subfolder names) to process.
# If None or empty, all subfolders containing a ledger file are processed.
sample_ids = []
# Example:
# sample_ids = ['Sample_01', 'Sample_02']

# =============================================================================
# Run
# =============================================================================

import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import load_sample_ledger, save_sample_ledger
from autoemx.utils import print_single_separator

directory = Path(directory)
if not directory.is_dir():
    raise NotADirectoryError(f"Provided directory does not exist: {directory}")

# Resolve sample subfolders
if sample_ids is not None and len(sample_ids) > 0:
    sample_dirs = [directory / sid for sid in sample_ids] # type: ignore
else:
    sample_dirs = [p for p in sorted(directory.iterdir()) if p.is_dir()]

ledger_filename = cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT

updated = []
skipped = []

for sample_dir in sample_dirs:
    ledger_path = sample_dir / ledger_filename
    if not ledger_path.exists():
        print(f"⚠️  No ledger found for '{sample_dir.name}', skipping.")
        skipped.append(sample_dir.name)
        continue

    try:
        ledger = load_sample_ledger(ledger_path)
    except Exception as e:
        print(f"❌ Failed to load ledger for '{sample_dir.name}': {e}")
        skipped.append(sample_dir.name)
        continue

    microscope_ID = ledger.configs.microscope_cfg.ID
    meas_mode = ledger.configs.measurement_cfg.mode

    try:
        calibs.load_microscope_calibrations(
            microscope_ID,
            meas_mode,
            load_detector_channel_params=True,
        )
    except Exception as e:
        print(f"❌ Failed to load calibrations for '{sample_dir.name}' "
              f"(microscope='{microscope_ID}', mode='{meas_mode}'): {e}")
        skipped.append(sample_dir.name)
        continue

    energy_zero = calibs.detector_channel_params[meas_mode][cnst.OFFSET_KEY]
    bin_width = calibs.detector_channel_params[meas_mode][cnst.SCALE_KEY]

    ledger.configs.microscope_cfg.energy_zero = energy_zero
    ledger.configs.microscope_cfg.bin_width = bin_width

    try:
        save_sample_ledger(ledger, ledger_path)
    except Exception as e:
        print(f"❌ Failed to save ledger for '{sample_dir.name}': {e}")
        skipped.append(sample_dir.name)
        continue

    print(f"✅ '{sample_dir.name}': energy_zero={energy_zero}, bin_width={bin_width}")
    updated.append(sample_dir.name)

print_single_separator()
print(f"Done. Updated {len(updated)} ledger(s), skipped {len(skipped)}.")
if skipped:
    print(f"Skipped: {skipped}")
