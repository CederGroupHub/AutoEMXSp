#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update EDS standards JSON file in place using the main runner from autoemx.runners.modify_standards_file.

Set the variables below as needed and run this script.
"""


# Path to the EDS standards JSON file
standards_json_path = '/absolute/path/to/EDS_Stds_15keV.json'  # <-- Set this

# List of dicts: {ID: el_line(s) or "all"} to delete
# Example: [{"BaSO4": "all"}, {"BaSO4": ["Ba_La1", "O_Ka1"]}]
delete_std_lines = None  # or []

# List of dicts: {ID: el_line(s) or "all"} to exclude from mean calc
# Example: [{"BaSO4": "all"}, {"BaSO4": ["Ba_La1", "O_Ka1"]}]
exclude_std_lines_from_mean_calc = None  # or []

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.modify_standards_file import modify_standards_file

modify_standards_file(
    standards_json_path,
    delete_std_lines=delete_std_lines,
    exclude_std_lines_from_mean_calc=exclude_std_lines_from_mean_calc
)
