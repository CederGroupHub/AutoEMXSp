#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update EDS standards JSON file in place using the main runner from autoemx.runners.modify_standards_file.

Set the variables below as needed and run this script. The reference means of the modified lines are recomputed,
and a timestamped backup of the file is saved next to it before it is changed.
"""
import os

# Path to the EDS standards JSON file. The standards used by AutoEMX are in the local calibration folder of the
# microscope, ~/.autoemx/calibrations/<microscope_ID>/
standards_json_path = os.path.expanduser('~/.autoemx/calibrations/PhenomXL/EDS_Stds_15keV.json')  # <-- Set this

# List of dicts: {ID: el_line(s) or "all"} to delete
# Example: [{"BaSO4": "all"}, {"BaSO4": ["Ba_La1", "O_Ka1"]}]
delete_std_lines = None  # or []

# List of dicts: {ID: el_line(s) or "all"} to exclude from mean calc
# Example: [{"BaSO4": "all"}, {"BaSO4": ["Ba_La1", "O_Ka1"]}]
exclude_std_lines_from_mean_calc = None  # or []

# El_lines to delete from all standards (entries and reference mean)
# Example: ["Mo_La1"]
delete_el_lines = None  # or []

mode = "point"  # Measurement mode of the standards to modify
backup = True  # If True, saves a copy of the file (<name>_backup_<date>_<time>.json) before changing it

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.modify_standards_file import modify_standards_file

modify_standards_file(
    standards_json_path,
    delete_std_lines=delete_std_lines,
    exclude_std_lines_from_mean_calc=exclude_std_lines_from_mean_calc,
    delete_el_lines=delete_el_lines,
    mode=mode,
    backup=backup,
)
