#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract and report details from experimental standards files.

This runner inspects standards JSON files and builds a human-readable report
with:
- coverage by voltage/mode/current
- number of elements represented per mode
- per-peak measured standards, corrected PB values, and mean corrected PB

The report is printed to terminal and written to a text file.
"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from statistics import mean
from typing import Any, Dict, List, Optional, Tuple

import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
import autoemx._paths as paths
from pymatgen.core import Element


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["extract_experimental_standards_details"]


_STD_FILENAME = cnst.STD_FILENAME
_STD_FILENAME_PATTERN = re.compile(
    rf"^(?P<meas_type>.+)_{re.escape(_STD_FILENAME)}_(?P<beam_kv>\d+)keV\.json$",
    re.IGNORECASE,
)
_BEAM_CURRENT_KEY = cnst.BEAM_CURRENT_KEY
_STD_ID_KEY = cnst.STD_ID_KEY
_STD_MEAN_ID = cnst.STD_MEAN_ID_KEY
_CORRECTED_PB_KEY = cnst.COR_PB_DF_KEY


def _extract_file_energy_and_type(std_path: Path) -> Tuple[int, str]:
    match = _STD_FILENAME_PATTERN.match(std_path.name)
    if match is None:
        raise ValueError(f"Invalid standards filename format: {std_path.name}")
    return int(match.group("beam_kv")), str(match.group("meas_type")).upper()


def _find_latest_detector_params(microscope_ID: str) -> Optional[Dict[str, Dict[str, float]]]:
    latest = paths.latest_detector_channel_params_file(microscope_ID)
    if latest is None:
        return None

    with latest.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if not isinstance(payload, dict):
        return None
    return payload


def _format_current(mode: str, detector_params: Optional[Dict[str, Dict[str, float]]]) -> str:
    if detector_params is None:
        return "unknown"
    mode_payload = detector_params.get(mode)
    if not isinstance(mode_payload, dict):
        return "unknown"
    current = mode_payload.get(_BEAM_CURRENT_KEY)
    if current is None:
        return "unknown"
    try:
        return f"{float(current):g}"
    except (TypeError, ValueError):
        return "unknown"


def _normalize_voltage(voltage: Optional[float]) -> Optional[int]:
    if voltage is None:
        return None
    if isinstance(voltage, bool):
        raise ValueError("voltage must be numeric")
    normalized = int(float(voltage))
    if normalized <= 0:
        raise ValueError("voltage must be positive")
    return normalized


def _normalize_line_payload(line_payload: Any) -> Tuple[List[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """Return (entries, reference_mean) from both legacy and schema formats."""
    if isinstance(line_payload, dict):
        entries = line_payload.get("entries", [])
        reference_mean = line_payload.get("reference_mean")
        if not isinstance(entries, list):
            entries = []
        if reference_mean is not None and not isinstance(reference_mean, dict):
            reference_mean = None
        return entries, reference_mean

    if isinstance(line_payload, list):
        entries: List[Dict[str, Any]] = []
        reference_mean: Optional[Dict[str, Any]] = None
        for row in line_payload:
            if not isinstance(row, dict):
                continue
            if str(row.get(_STD_ID_KEY, "")).strip() == _STD_MEAN_ID:
                reference_mean = row
            else:
                entries.append(row)
        return entries, reference_mean

    return [], None


def _to_float_or_none(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _sort_elements_by_atomic_number(elements: List[str]) -> List[str]:
    """Sort element symbols by atomic number using pymatgen."""
    def _sort_key(symbol: str) -> Tuple[int, str]:
        try:
            return Element(symbol).Z, symbol
        except Exception:
            return 10**9, symbol

    return sorted(elements, key=_sort_key)


def extract_experimental_standards_details(
    microscope_ID: str = dflt.microscope_ID,
    voltage: Optional[float] = None,
    standards_json_path: Optional[str] = None,
    report_output_dir: Optional[str] = None,
    print_report: bool = True,
    print_per_peak: bool = True,
) -> str:
    """Summarize standards coverage and per-peak corrected PB entries.

    Loads standards from either:
    - ``standards_json_path`` when provided, or
    - ``~/.autoemx/calibrations/<microscope_ID>/*_Stds_*keV.json`` (local calibration folder) by default.

    Writes a text report and prints it to stdout.
    """
    if not microscope_ID or not str(microscope_ID).strip():
        raise ValueError("microscope_ID cannot be empty")

    requested_voltage = _normalize_voltage(voltage)
    calib_dir = paths.user_calibration_dir(microscope_ID)

    source_label = "local_calibrations"
    output_dir = calib_dir

    if report_output_dir is not None:
        output_dir = Path(report_output_dir).expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

    if standards_json_path is not None:
        std_path = Path(standards_json_path).expanduser().resolve()
        if not std_path.exists() or not std_path.is_file():
            raise FileNotFoundError(f"Provided standards JSON file does not exist: '{std_path}'")
        all_std_files = [std_path]
        source_label = "custom_project_file"
        if report_output_dir is None:
            output_dir = std_path.parent
    else:
        all_std_files = sorted(calib_dir.glob("*_Stds_*keV.json"))
        if requested_voltage is not None:
            all_std_files = [
                path for path in all_std_files if f"_{requested_voltage}keV.json" in path.name
            ]

    if not all_std_files:
        if requested_voltage is None:
            raise FileNotFoundError(
                f"No standards files matching '*_Stds_*keV.json' were found in '{calib_dir}'."
            )
        raise FileNotFoundError(
            f"No standards file for {requested_voltage} keV was found in '{calib_dir}'."
        )

    detector_params = _find_latest_detector_params(microscope_ID)

    report_lines: List[str] = []
    report_lines.append("=" * 80)
    report_lines.append("Experimental standards summary")
    report_lines.append("-" * 80)
    report_lines.append(f"Microscope ID: {microscope_ID}")
    report_lines.append(f"Standards source: {source_label}")
    if standards_json_path is not None:
        report_lines.append(f"Custom standards file: {Path(standards_json_path).name}")
    report_lines.append(
        f"Voltage filter: {requested_voltage} keV" if requested_voltage is not None else "Voltage filter: all"
    )
    report_lines.append(f"Standards files found: {len(all_std_files)}")
    report_lines.append("")
    report_lines.append("Coverage by voltage and current")
    report_lines.append("-" * 80)

    per_file_sections: List[str] = []
    for std_path in all_std_files:
        with std_path.open("r", encoding="utf-8") as handle:
            raw_payload = json.load(handle)

        beam_kv_raw = raw_payload.get("beam_energy_keV")
        beam_kv: Optional[int] = None
        if beam_kv_raw is not None:
            try:
                beam_kv = int(beam_kv_raw)
            except (TypeError, ValueError):
                beam_kv = None

        try:
            beam_kv_from_name, file_meas_type = _extract_file_energy_and_type(std_path)
            if beam_kv is None:
                beam_kv = beam_kv_from_name
        except ValueError:
            file_meas_type = str(raw_payload.get("measurement_type", "EDS")).upper()

        if beam_kv is None:
            raise ValueError(
                f"Could not determine beam energy for standards file '{std_path}'. "
                "Include 'beam_energy_keV' in the JSON payload or use '*_Stds_<voltage>keV.json' naming."
            )

        if requested_voltage is not None and beam_kv != requested_voltage:
            continue

        measurement_type = str(raw_payload.get("measurement_type", file_meas_type))
        standards_by_mode_raw = raw_payload.get("standards_by_mode", {})
        if not isinstance(standards_by_mode_raw, dict):
            standards_by_mode_raw = {}

        report_lines.append(f"Voltage {beam_kv} keV ({std_path.name})")
        if not standards_by_mode_raw:
            report_lines.append("  - No modes found in standards file.")
            continue

        for mode, lines_by_peak_raw in standards_by_mode_raw.items():
            if not isinstance(lines_by_peak_raw, dict):
                continue
            current_txt = _format_current(mode, detector_params)


            # Extract unique element symbols in order of first appearance
            seen = set()
            unique_elements = []
            for peak_name in lines_by_peak_raw.keys():
                el = str(peak_name).split("_", maxsplit=1)[0]
                if el not in seen:
                    seen.add(el)
                    unique_elements.append(el)
            elements_txt = ", ".join(unique_elements) if unique_elements else "none"

            report_lines.append(
                "  - "
                f"Mode={mode}; current={current_txt}; "
                f"n_elements={len(unique_elements)}; n_peaks={len(lines_by_peak_raw)}; "
                f"elements=[{elements_txt}]"
            )

        section_lines: List[str] = []
        section_lines.append("")
        section_lines.append("-" * 80)
        section_lines.append(f"Per-peak standards at {beam_kv} keV ({std_path.name})")
        section_lines.append("-" * 80)

        section_lines.append(f"Measurement type: {measurement_type}")

        for mode, lines_by_peak_raw in standards_by_mode_raw.items():
            if not isinstance(lines_by_peak_raw, dict):
                continue
            current_txt = _format_current(mode, detector_params)
            section_lines.append(f"Mode: {mode} | current: {current_txt}")

            if not lines_by_peak_raw:
                section_lines.append("  (no peak data)")
                section_lines.append("")
                continue

            for peak_name, line_payload in lines_by_peak_raw.items():
                entries, reference_mean = _normalize_line_payload(line_payload)

                entry_ids_and_pb: List[Tuple[str, float]] = []
                entry_use_for_mean: Dict[str, str] = {}
                for row in entries:
                    std_id = str(row.get(_STD_ID_KEY, "")).strip()
                    corrected_pb = _to_float_or_none(row.get(_CORRECTED_PB_KEY))
                    if not std_id or corrected_pb is None:
                        continue
                    entry_ids_and_pb.append((std_id, corrected_pb))

                    use_for_mean_raw = row.get(cnst.STD_USE_FOR_MEAN_KEY)
                    if isinstance(use_for_mean_raw, bool):
                        entry_use_for_mean[std_id] = str(use_for_mean_raw)
                    else:
                        entry_use_for_mean[std_id] = "unknown"

                pb_values = [pb for _, pb in entry_ids_and_pb]

                if reference_mean is not None:
                    overall_mean = _to_float_or_none(reference_mean.get(_CORRECTED_PB_KEY))
                else:
                    overall_mean = None

                if overall_mean is None and pb_values:
                    overall_mean = float(mean(pb_values))

                section_lines.append(f"  Peak {peak_name}")
                if not entry_ids_and_pb:
                    section_lines.append("    - No measured standards entries")
                else:
                    for std_id, corrected_pb in entry_ids_and_pb:
                        use_for_mean_txt = entry_use_for_mean.get(std_id, "unknown")
                        section_lines.append(
                            f"    - ID={std_id}; Corrected_PB={corrected_pb:.6g}; "
                            f"Use_for_mean_calc={use_for_mean_txt}"
                        )

                if overall_mean is not None:
                    section_lines.append(f"    - Mean Corrected_PB={overall_mean:.6g}")
                else:
                    section_lines.append("    - Mean Corrected_PB=nan")

            section_lines.append("")

        per_file_sections.extend(section_lines)

    if not per_file_sections:
        raise FileNotFoundError(
            "No standards data matched the requested voltage filter in the selected standards source."
        )

    summary_report_lines = list(report_lines)
    report_lines.append("")
    report_lines.extend(per_file_sections)
    report_text = "\n".join(report_lines).rstrip() + "\n"

    suffix = f"_{requested_voltage}keV" if requested_voltage is not None else ""
    out_path = output_dir / f"experimental_standards_details{suffix}.txt"
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(report_text)

    if print_report:
        if print_per_peak:
            print(report_text)
        else:
            summary_text = "\n".join(summary_report_lines).rstrip() + "\n"
            print(summary_text)
    logging.info("Saved report to: %s", out_path)
    return str(out_path)
