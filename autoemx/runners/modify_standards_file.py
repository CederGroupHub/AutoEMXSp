#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: Modify EDS Standards JSON File In-Place

This runner script updates an EDS standards JSON file in place. It supports:
    - Deleting standard entries by ID and el_line (or "all")
    - Excluding entries from mean calculation by ID and el_line (or "all")
    - Recomputing the mean values using only entries with Use_for_mean_calc = True
    - Setting reference_mean to null if no valid entries remain

How to use:
    1. Set the variables at the top of this script (file path, delete/exclude lists).
    2. Run the script directly (as a runner) to apply the changes.
    3. The script prints a summary of all actions taken (deletions, exclusions, means recomputed, etc).

Intended for use as a command-line runner for standards file maintenance and batch updates.

Created on Tue May 19 2026
@author: Andrea
"""

import json
import numpy as np

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


def modify_standards_file(
    file_path,
    delete_std_lines=None,
    exclude_std_lines_from_mean_calc=None
):

    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    delete_lookup = _build_lookup(delete_std_lines)
    exclude_lookup = _build_lookup(exclude_std_lines_from_mean_calc)

    standards = data["standards_by_mode"]["point"]
    el_lines_to_remove = []
    deleted_entries = []
    excluded_entries = []
    not_found_deletions = set()
    not_found_exclusions = set()
    means_recomputed = []
    means_removed = []

    # Track all el_lines and IDs for reporting
    all_el_lines = set(standards.keys())
    all_ids = set()
    for std in standards.values():
        for entry in std["entries"]:
            if "ID" in entry:
                all_ids.add(entry["ID"])

    for el_line, std in list(standards.items()):
        new_entries = []
        for entry in std["entries"]:
            entry_id = entry.get("ID")
            deleted = False
            # Deletion
            if entry_id in delete_lookup:
                if delete_lookup[entry_id] == "all" or el_line in delete_lookup[entry_id]:
                    deleted_entries.append((entry_id, el_line))
                    deleted = True
            if deleted:
                continue  # skip (delete)
            # Exclusion from mean calc
            if entry_id in exclude_lookup:
                if exclude_lookup[entry_id] == "all" or el_line in exclude_lookup[entry_id]:
                    if entry.get("Use_for_mean_calc", True):
                        excluded_entries.append((entry_id, el_line))
                    entry["Use_for_mean_calc"] = False
            new_entries.append(entry)
        std["entries"] = new_entries

        # If no entries remain, mark for removal
        if len(std["entries"]) == 0:
            el_lines_to_remove.append(el_line)
            continue

        # Recompute means, but only mark as recomputed if changed
        mean_fields = [
            "Corrected_PB", "Stdev_PB", "Rel_stdev_PB (%)", "Measured_PB"
        ]
        mean_entry = {k: [] for k in mean_fields}
        datetimes = []
        for entry in std["entries"]:
            if entry.get("Use_for_mean_calc", True):
                for k in mean_fields:
                    if k in entry:
                        mean_entry[k].append(entry[k])
                if "datetime" in entry:
                    datetimes.append(entry["datetime"])
        reference_mean = {}
        if mean_entry["Corrected_PB"]:
            new_mean = float(np.mean(mean_entry["Corrected_PB"]))
            new_stdev = float(np.std(mean_entry["Corrected_PB"]))
            if new_mean != 0:
                new_rel_stdev = new_stdev / new_mean * 100
            else:
                new_rel_stdev = 0.0
            new_datetime = datetimes[-1] if datetimes else None
            reference_mean = {
                "Corrected_PB": new_mean,
                "Stdev_PB": new_stdev,
                "Rel_stdev_PB (%)": new_rel_stdev,
                "datetime": new_datetime
            }
            prev_mean = std.get("reference_mean", {})
            # Only append to means_recomputed if any value changed
            mean_changed = (
                not prev_mean or
                any(abs(reference_mean.get(k, None) - prev_mean.get(k, None)) > 1e-12 if isinstance(reference_mean.get(k, None), float) else reference_mean.get(k, None) != prev_mean.get(k, None) for k in reference_mean)
            )
            std["reference_mean"] = reference_mean
            if mean_changed:
                means_recomputed.append(el_line)
        else:
            # No valid entries for mean: set reference_mean to null (None in Python)
            if std.get("reference_mean", None) is not None:
                means_removed.append(el_line)
            std["reference_mean"] = None

    # Remove el_lines with no entries
    for el_line in el_lines_to_remove:
        del standards[el_line]

    # Check for not found deletions/exclusions
    if delete_lookup:
        for id_ in delete_lookup:
            if id_ not in all_ids:
                not_found_deletions.add(id_)
            elif delete_lookup[id_] != "all":
                for el in delete_lookup[id_]:
                    if el not in all_el_lines:
                        not_found_deletions.add(f"{id_}:{el}")
    if exclude_lookup:
        for id_ in exclude_lookup:
            if id_ not in all_ids:
                not_found_exclusions.add(id_)
            elif exclude_lookup[id_] != "all":
                for el in exclude_lookup[id_]:
                    if el not in all_el_lines:
                        not_found_exclusions.add(f"{id_}:{el}")

    # Validate and serialize using EDSStandardsFile to enforce schema
    from autoemx.config.schema_models.standards import EDSStandardsFile
    # Extract required fields for model
    schema_version = data.get("schema_version", 1)
    measurement_type = data["measurement_type"]
    beam_energy_keV = data["beam_energy_keV"]
    standards_by_mode = data["standards_by_mode"]
    # Build and validate model
    model = EDSStandardsFile.from_payload(
        payload={
            "schema_version": schema_version,
            "measurement_type": measurement_type,
            "beam_energy_keV": beam_energy_keV,
            "standards_by_mode": standards_by_mode,
        },
        meas_type=measurement_type,
        beam_energy_keV=beam_energy_keV,
    )
    # Write back using model's to_json_file (handles sorting and aliases)
    model.to_json_file(file_path, indent=2)

    # Print summary
    has_stds_file_changed = False
    print("\n--- EDS Standards Update Summary ---")
    if deleted_entries:
        has_stds_file_changed = True
        print(f"Deleted entries:")
        for entry_id, el_line in deleted_entries:
            print(f"  - ID: {entry_id}, el_line: {el_line}")
    else:
        print("No entries deleted.")
    if not_found_deletions:
        print("Entries requested for deletion not found:")
        for item in not_found_deletions:
            print(f"  - {item}")
    if excluded_entries:
        has_stds_file_changed = True
        print(f"Excluded from mean calculation:")
        for entry_id, el_line in excluded_entries:
            print(f"  - ID: {entry_id}, el_line: {el_line}")
    else:
        print("No entries excluded from mean calculation.")
    if not_found_exclusions:
        print("Entries requested for exclusion not found:")
        for item in not_found_exclusions:
            print(f"  - {item}")
    if el_lines_to_remove:
        print(f"Removed el_lines with no remaining entries:")
        for el_line in el_lines_to_remove:
            print(f"  - {el_line}")
    if means_recomputed:
        print(f"Recomputed means for el_lines:")
        for el_line in means_recomputed:
            print(f"  - {el_line}")
    if means_removed:
        print(f"Removed reference_mean for el_lines (no valid entries):")
        for el_line in means_removed:
            print(f"  - {el_line}")
    print("--- End of Summary ---\n")

    # After modifying standards, call extract_experimental_standards_details and save report in same folder
    try:
        from autoemx.runners.extract_experimental_standards_details import extract_experimental_standards_details
        from pathlib import Path
        standards_dir = str(Path(file_path).parent)
        extract_experimental_standards_details(
            standards_json_path=file_path,
            report_output_dir=standards_dir,
            print_report=False
        )
        print(f"Experimental standards details report saved in: {standards_dir}")
    except Exception as e:
        print(f"[Warning] Could not generate standards details report: {e}")

    if has_stds_file_changed:
        print(f"Updated standards file: {file_path}")