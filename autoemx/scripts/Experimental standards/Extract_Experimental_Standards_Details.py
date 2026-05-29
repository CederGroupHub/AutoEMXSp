#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract experimental standards details.

Use this script to generate a text summary of standards available for a
microscope (or for a custom project standards file), including:
- voltages/modes/currents covered
- number of elements and peaks
- per-peak measured standard IDs and corrected PB values
- mean corrected PB per peak

The summary is printed in terminal and saved to a txt file.

Created on Tue May 12 2026

@author: Andrea
"""

from pathlib import Path

# =============================================================================
# Input
# =============================================================================

microscope_ID = 'PhenomXL'

# Optional voltage filter in keV. If None, include all voltages found.
voltage = None

# Optional path to a project-specific standards JSON file.
# If None, standards are loaded from autoemx/calibrations/<microscope_ID>/
# Example:
# standards_json_path = '/absolute/path/to/Results/<sample_ID>/EDS_Stds_15keV.json'
standards_json_path = None

# Directory where the txt report is saved. Default: this script's directory.
report_output_dir = str(Path(__file__).resolve().parent)

# =============================================================================
# Run
# =============================================================================

from autoemx.runners.extract_experimental_standards_details import extract_experimental_standards_details  # type: ignore

report_path = extract_experimental_standards_details(
    microscope_ID=microscope_ID,
    voltage=voltage,
    standards_json_path=standards_json_path,
    report_output_dir=report_output_dir,
)
