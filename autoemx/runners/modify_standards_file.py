#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: Modify EDS Standards JSON File In-Place

This runner updates an EDS standards JSON file in place. It supports:
    - Deleting standard entries by ID and el_line (or "all")
    - Deleting el_lines from all standards
    - Excluding entries from mean calculation by ID and el_line (or "all")
    - Recomputing the reference means of the modified el_lines from the entries used for the mean
      (Use_for_mean_calc = True)
    - Setting reference_mean to null if no entry used for the mean remains

Before the file is changed, a timestamped backup copy is saved next to it.

How to use:
    1. Set the variables in scripts/Experimental standards/Modify_EDS_Standards_File.py.
    2. Run the script to apply the changes.
    3. A summary of all actions taken (deletions, exclusions, means recomputed, etc.) is printed.

Created on Tue May 19 2026
@author: Andrea
"""

import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from autoemx.config.schema_models.standards import EDSStandardsFile


def _build_lookup(list_of_dicts):
    lookup = {}
    if not list_of_dicts:
        return lookup
    for d in list_of_dicts:
        for k, v in d.items():
            if v == "all":
                lookup[k] = "all"
            elif isinstance(v, list):
                if k not in lookup or lookup[k] == "all":
                    lookup[k] = set()
                lookup[k].update(v)
            else:
                if k not in lookup or lookup[k] == "all":
                    lookup[k] = set()
                lookup[k].add(v)
    return lookup


def _matches(lookup: Dict[str, Any], entry_id: str, el_line: str) -> bool:
    return entry_id in lookup and (lookup[entry_id] == "all" or el_line in lookup[entry_id])


def _not_found(lookup: Dict[str, Any], all_ids: set, all_el_lines: set) -> List[str]:
    not_found = []
    for id_, el_lines in lookup.items():
        if id_ not in all_ids:
            not_found.append(id_)
        elif el_lines != "all":
            not_found.extend(f"{id_}:{el}" for el in el_lines if el not in all_el_lines)
    return not_found


def _backup_path(file_path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return file_path.with_name(f"{file_path.stem}_backup_{stamp}{file_path.suffix}")


def modify_standards_file(
    file_path,
    delete_std_lines=None,
    exclude_std_lines_from_mean_calc=None,
    delete_el_lines: Optional[List[str]] = None,
    mode: str = "point",
    backup: bool = True,
) -> Dict[str, Any]:
    """
    Modify an EDS standards JSON file in place.

    Parameters
    ----------
    file_path : str or Path
        Path to the standards file (e.g. EDS_Stds_15keV.json).
    delete_std_lines : list of dict, optional
        Entries to delete, as {standard ID: el_line(s) or "all"}, e.g. [{"BaSO4": "all"}, {"BaSO4": ["Ba_La1"]}].
    exclude_std_lines_from_mean_calc : list of dict, optional
        Entries to exclude from the reference means, in the same format as delete_std_lines.
    delete_el_lines : list of str, optional
        El_lines deleted from all standards (entries and reference mean), e.g. ["Mo_La1"].
    mode : str, optional
        Measurement mode of the standards to modify (key of standards_by_mode). Default is "point".
    backup : bool, optional
        If True (default), a copy of the file is saved as <name>_backup_<date>_<time>.json before it is changed.

    Returns
    -------
    dict
        Summary of the changes: deleted_entries, excluded_entries, deleted_el_lines, means_recomputed,
        means_removed, not_found_deletions, not_found_exclusions, not_found_el_lines, backup_path, changed.
    """
    file_path = Path(file_path)
    with file_path.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    standards_file = EDSStandardsFile.from_payload(
        payload, meas_type=payload["measurement_type"], beam_energy_keV=payload["beam_energy_keV"]
    )
    if mode not in standards_file.standards_by_mode:
        raise ValueError(
            f"Mode '{mode}' not found in {file_path.name}. Available modes: {sorted(standards_file.standards_by_mode)}"
        )
    standards = standards_file.standards_by_mode[mode]

    delete_lookup = _build_lookup(delete_std_lines)
    exclude_lookup = _build_lookup(exclude_std_lines_from_mean_calc)
    all_el_lines = set(standards)
    all_ids = {entry.standard_id for line in standards.values() for entry in line.entries}

    summary: Dict[str, Any] = {
        "deleted_entries": [], "excluded_entries": [], "deleted_el_lines": [],
        "means_recomputed": [], "means_removed": [],
        "not_found_deletions": _not_found(delete_lookup, all_ids, all_el_lines),
        "not_found_exclusions": _not_found(exclude_lookup, all_ids, all_el_lines),
        "not_found_el_lines": [el_line for el_line in (delete_el_lines or []) if el_line not in all_el_lines],
        "backup_path": None,
    }
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    for el_line in list(standards):
        line = standards[el_line]
        if el_line in (delete_el_lines or []):
            del standards[el_line]
            summary["deleted_el_lines"].append(el_line)
            continue

        kept_entries = []
        is_line_modified = False
        for entry in line.entries:
            if _matches(delete_lookup, entry.standard_id, el_line):
                summary["deleted_entries"].append((entry.standard_id, el_line))
                is_line_modified = True
                continue
            if _matches(exclude_lookup, entry.standard_id, el_line) and entry.use_for_mean_calc:
                entry.use_for_mean_calc = False
                summary["excluded_entries"].append((entry.standard_id, el_line))
                is_line_modified = True
            kept_entries.append(entry)
        line.entries = kept_entries
        if not is_line_modified:
            # Untouched lines keep their reference mean (recomputing it from the entries, rounded to 0.1, could
            # differ slightly from the stored one, computed from unrounded values)
            continue

        if not line.entries:
            # No entry left for this line
            del standards[el_line]
            summary["deleted_el_lines"].append(el_line)
            continue

        new_mean = line.compute_reference_mean(now)
        old_mean = line.reference_mean
        if new_mean is None:
            if old_mean is not None:
                line.reference_mean = None
                summary["means_removed"].append(el_line)
        else:
            line.reference_mean = new_mean
            summary["means_recomputed"].append(el_line)

    summary["changed"] = any(
        summary[key] for key in ("deleted_entries", "excluded_entries", "deleted_el_lines", "means_recomputed",
                                 "means_removed")
    )
    if summary["changed"]:
        if backup:
            summary["backup_path"] = str(_backup_path(file_path))
            shutil.copy2(file_path, summary["backup_path"])
        standards_file.to_json_file(file_path, indent=2)

    _print_summary(summary, file_path)

    # Report on the modified standards, saved in the same folder
    if summary["changed"]:
        try:
            from autoemx.runners.extract_experimental_standards_details import extract_experimental_standards_details
            extract_experimental_standards_details(
                standards_json_path=str(file_path),
                report_output_dir=str(file_path.parent),
                print_report=False,
            )
            print(f"Experimental standards details report saved in: {file_path.parent}")
        except Exception as e:
            print(f"[Warning] Could not generate standards details report: {e}")

    return summary


def _print_summary(summary: Dict[str, Any], file_path: Path) -> None:
    print("\n--- EDS Standards Update Summary ---")
    sections = [
        ("deleted_entries", "Deleted entries:", "No entries deleted."),
        ("excluded_entries", "Excluded from mean calculation:", "No entries excluded from mean calculation."),
    ]
    for key, title, empty_message in sections:
        if summary[key]:
            print(title)
            for entry_id, el_line in summary[key]:
                print(f"  - ID: {entry_id}, el_line: {el_line}")
        else:
            print(empty_message)
    for key, title in (
        ("not_found_deletions", "Entries requested for deletion not found:"),
        ("not_found_exclusions", "Entries requested for exclusion not found:"),
        ("not_found_el_lines", "El_lines requested for deletion not found:"),
        ("deleted_el_lines", "Removed el_lines (deleted, or with no remaining entries):"),
        ("means_recomputed", "Recomputed means for el_lines:"),
        ("means_removed", "Removed reference_mean for el_lines (no entry used for the mean):"),
    ):
        if summary[key]:
            print(title)
            for item in summary[key]:
                print(f"  - {item}")
    print("--- End of Summary ---\n")
    if summary["changed"]:
        if summary["backup_path"]:
            print(f"Backup of the previous file: {summary['backup_path']}")
        print(f"Updated standards file: {file_path}")
    else:
        print(f"No changes: {file_path} was not modified.")
