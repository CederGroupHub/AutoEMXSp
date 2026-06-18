#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate relative X-ray line weights from repeated fits of an experimental standard.

Fits multiple spectra from one standard and records fitted peak areas for user-selected
lines (`free_area_el_lines`). Computes per-line mean/std and pairwise area correlations,
including a heuristic check for anti-correlated pairs with approximately constant area sum.
"""

from __future__ import annotations

import logging
import os
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from autoemx.config.ledger_io import load_sample_ledger
from autoemx.data.Xray_lines import get_el_xray_lines
from autoemx.runners.fit_and_quantify_spectrum_from_ledger import (
    fit_and_quantify_spectrum_from_ledger,
)
from autoemx.core.fitter import Background_Model
from autoemx.utils.helper import get_sample_dir, print_double_separator
import autoemx.utils.constants as cnst


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = [
    "calibrate_peak_weights_experimental_standard",
    "recompute_calibration",
]


def _resolve_results_root(path: Optional[str]) -> Path:
    if path is None:
        return Path(os.getcwd()) / cnst.RESULTS_DIR
    return Path(path).expanduser().resolve()


def _resolve_spectrum_ids(ledger, spectrum_ids_to_fit: Optional[Sequence[Any]]) -> List[Any]:
    if spectrum_ids_to_fit is not None:
        return list(spectrum_ids_to_fit)

    resolved = []
    for idx, spectrum_entry in enumerate(ledger.spectra):
        entry_id = spectrum_entry.spectrum_id
        if entry_id is None:
            resolved.append(idx)
            continue
        try:
            resolved.append(int(float(entry_id)))
        except (TypeError, ValueError):
            resolved.append(entry_id)
    return resolved


def _normalize_spectrum_id(value: Any) -> Any:
    """Normalize spectrum identifiers for robust matching (e.g., 3, 3.0, '3')."""
    if value is None:
        return None
    try:
        as_num = float(value)
    except (TypeError, ValueError):
        return str(value)

    if np.isfinite(as_num) and as_num.is_integer():
        return int(as_num)
    return as_num

def _load_abs_factors_from_table(
    ratios_abs_corrected_path: Path,
    free_area_el_lines: Sequence[str],
) -> Dict[str, float]:
    """
    Load per-line absorption-correction factors from a previously written
    'ratios_abs_corrected' table.

    The table has one row per line (keyed by the 'Line' column, e.g. 'Ma1') and a
    'Correction factor' column. We map each requested element-line back to its row by
    matching the line key, consuming rows positionally to disambiguate repeated keys.
    """
    if not ratios_abs_corrected_path.exists():
        raise FileNotFoundError(
            "Absorption-corrected ratios table not found for recompute. Run a normal "
            f"calibration first: {ratios_abs_corrected_path}"
        )

    table = pd.read_csv(ratios_abs_corrected_path)
    if "Line" not in table.columns or "Correction factor" not in table.columns:
        raise ValueError(
            "Absorption-corrected ratios table is missing required 'Line' and/or "
            f"'Correction factor' columns. File: {ratios_abs_corrected_path}"
        )

    line_key_series = table["Line"].astype(str)
    used_row_idx: set = set()
    factors: Dict[str, float] = {}
    for el_line in free_area_el_lines:
        line_key = _line_key(el_line)
        candidate_idx = [
            idx for idx in line_key_series[line_key_series == line_key].index
            if idx not in used_row_idx
        ]
        if len(candidate_idx) == 0:
            logging.warning(
                "Line '%s' (from '%s') not found in absorption-corrected ratios table; "
                "defaulting correction factor to 1.0. File: %s",
                line_key, el_line, ratios_abs_corrected_path,
            )
            factors[el_line] = 1.0
            continue
        row_idx = candidate_idx[0]
        used_row_idx.add(row_idx)

        factor_val = pd.to_numeric(table.loc[row_idx, "Correction factor"], errors="coerce")
        factors[el_line] = float(factor_val) if pd.notna(factor_val) else 1.0

    return factors

def _absorption_factor_per_line(
    free_area_el_lines: Sequence[str],
    quantifier,
    ref_peak: str,
    direction: str = "to_pure",
) -> Dict[str, float]:
    """
    Compute per-line absorption-correction factors for the standard, expressed as a
    *double ratio against ref_peak*, to match the convention used by
    _update_peak_weights in the quantifier.

    The quantifier expresses each dependent line's area as
        A(L) = A(R) * weight_NIST(L) * (a_pure(R)/a_pure(L)) * (a_sample(L)/a_sample(R))
    so the absorption-corrected ratio A(L)/A(R) -> pure-element ratio requires
        factor(L) = (a_pure(L)/a_sample(L)) / (a_pure(R)/a_sample(R))

    With this definition, factor(ref_peak) == 1 by construction, and the corrected
    ratio  (A_meas(L)/A_meas(R)) * factor(L)  equals the pure-element area ratio,
    identical to what the quantifier's dependent-line expression encodes.

    direction="to_pure":  ratio_corr = ratio_meas * factor(L)   [double ratio above]
    direction="strip":    factor(L) = a_sample(R) / a_sample(L)  [remove only self-abs]
    """
    fitted_params_for_calcs = quantifier.fit_result.params.copy()
    fitted_params_for_calcs['apply_det_response'].value = 0
    Background_Model(quantifier.is_particle,
                    beam_energy=quantifier.beam_energy,
                    emergence_angle=quantifier.emergence_angle)

    def _a_sample(el: str, line: str) -> float:
        en = np.array([get_el_xray_lines(el)[line]["energy (keV)"]])
        Background_Model._clear_cached_abs_att_variables()
        return float(Background_Model._abs_attenuation_phirho(en, **fitted_params_for_calcs)[0])

    def _a_pure(el: str, line: str) -> float:
        en = np.array([get_el_xray_lines(el)[line]["energy (keV)"]])
        Background_Model._clear_cached_abs_att_variables()
        pure_kwargs = {"det_angle": quantifier.emergence_angle, "adr_abs": 0, f"f_{el}": 1.0}
        return float(Background_Model._abs_attenuation_phirho(en, **pure_kwargs)[0])

    # Reference-peak absorption terms (denominator of the double ratio).
    if "_" not in ref_peak:
        raise ValueError(
            f"Reference peak '{ref_peak}' is not an element-line label; cannot build "
            "absorption double ratio against it."
        )
    ref_el, ref_line = ref_peak.split("_", 1)
    a_sample_ref = _a_sample(ref_el, ref_line)
    a_pure_ref = _a_pure(ref_el, ref_line)

    if a_sample_ref == 0 or a_pure_ref == 0:
        raise ValueError(
            f"Reference peak '{ref_peak}' has a zero absorption term "
            f"(a_pure={a_pure_ref}, a_sample={a_sample_ref}); cannot normalize."
        )

    factors: Dict[str, float] = {}
    for el_line in free_area_el_lines:
        if "_" not in el_line:
            factors[el_line] = 1.0
            continue
        el, line = el_line.split("_", 1)

        a_sample = _a_sample(el, line)

        if direction == "strip":
            # Remove self-absorption relative to the reference peak's self-absorption,
            # so the reference peak itself is left unchanged (factor == 1).
            factors[el_line] = a_sample_ref / a_sample if a_sample != 0 else 1.0
            continue

        a_pure = _a_pure(el, line)

        # Double ratio: (a_pure(L)/a_sample(L)) / (a_pure(R)/a_sample(R))
        factors[el_line] = (a_pure / a_sample) / (a_pure_ref / a_sample_ref)

    return factors


def _apply_absorption_factors(
    raw_df: pd.DataFrame,
    free_area_el_lines: Sequence[str],
    factors: Dict[str, float],
) -> pd.DataFrame:
    """Multiply each fitted area column by its constant absorption factor."""
    out = raw_df.copy()
    for el_line in free_area_el_lines:
        if el_line in out.columns:
            out[el_line] = pd.to_numeric(out[el_line], errors="coerce") * factors.get(el_line, 1.0)
    return out

def _collect_fit_row(
    std_id: str,
    spectrum_id: Any,
    quantifier,
    free_lines: Sequence[str],
) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        "std_ID": std_id,
        "spectrum_ID": spectrum_id,
        "redchi": getattr(getattr(quantifier, "fit_result", None), "redchi", np.nan),
    }

    fit_params = getattr(getattr(quantifier, "fit_result", None), "params", None)

    for el_line in free_lines:
        area_key = f"{el_line}_area"
        if fit_params is None or area_key not in fit_params:
            row[el_line] = np.nan
            continue
        row[el_line] = float(fit_params[area_key].value)

    return row

def _build_summary_table(
    std_id: str,
    areas_df: pd.DataFrame,
    free_lines: Sequence[str],
    line_energies: Dict[str, float],
    sum_cv_threshold: float,
    anti_corr_threshold: float = -0.7,
    split_cv_ratio: float = 0.5,
) -> pd.DataFrame:
    summary_rows: List[Dict[str, Any]] = []

    for el_line in free_lines:
        vals = pd.to_numeric(areas_df[el_line], errors="coerce").dropna()
        n = int(vals.shape[0])
        mean_v = float(vals.mean()) if n else np.nan
        std_v = float(vals.std(ddof=1)) if n > 1 else np.nan
        cv_v = float(std_v / mean_v) if n > 1 and np.isfinite(mean_v) and mean_v != 0 else np.nan

        summary_rows.append(
            {
                "record_type": "line_summary",
                "std_ID": std_id,
                "line_A": el_line,
                "line_B": "",
                "n": n,
                "mean_area": mean_v,
                "std_area": std_v,
                "cv_area": cv_v,
                "pearson_r": np.nan,
                "cv_A": np.nan,
                "cv_B": np.nan,
                "sum_mean": np.nan,
                "sum_std": np.nan,
                "sum_cv": np.nan,
                "is_anti_corr": False,
                "is_sum_locked_pair": False,
                "is_split_pair": False,
                "signal": "",
            }
        )

    # Pair diagnostics are computed only for adjacent-energy neighbors.
    for line_a, line_b in combinations(free_lines, 2):
        if not _are_adjacent_by_energy(line_a, line_b, line_energies):
            continue

        pair_df = areas_df[[line_a, line_b]].dropna()
        n = int(pair_df.shape[0])

        if n < 3:
            corr_v = sum_mean = sum_std = sum_cv = np.nan
            cv_a = cv_b = np.nan
            is_anti_corr = is_sum_locked = is_split_pair = False
        else:
            a = pair_df[line_a]
            b = pair_df[line_b]
            corr_v = float(a.corr(b))

            pair_sum = a + b
            sum_mean = float(pair_sum.mean())
            sum_std = float(pair_sum.std(ddof=1))
            sum_cv = float(sum_std / sum_mean) if sum_mean != 0 else np.nan

            # Per-peak CVs: the splitting signature is that each peak varies
            # much more than their sum.
            cv_a = float(a.std(ddof=1) / a.mean()) if a.mean() != 0 else np.nan
            cv_b = float(b.std(ddof=1) / b.mean()) if b.mean() != 0 else np.nan
            max_peak_cv = np.nanmax([cv_a, cv_b])

            is_anti_corr = bool(np.isfinite(corr_v) and corr_v <= anti_corr_threshold)
            is_sum_locked = bool(np.isfinite(sum_cv) and sum_cv <= sum_cv_threshold)
            # The sum must be substantially more stable than the noisier peak.
            is_variance_in_split = bool(
                np.isfinite(sum_cv)
                and np.isfinite(max_peak_cv)
                and max_peak_cv > 0
                and (sum_cv / max_peak_cv) <= split_cv_ratio
            )
            # A pair is a true "split" only if all three hold.
            is_split_pair = is_anti_corr and is_sum_locked and is_variance_in_split

        signal_tokens = []
        if is_split_pair:
            signal_tokens.append("SPLIT_PAIR")
        elif is_anti_corr:
            signal_tokens.append("ANTI_CORR_ONLY")
        elif is_sum_locked:
            signal_tokens.append("SUM_LOCKED_ONLY")

        summary_rows.append(
            {
                "record_type": "pair_correlation",
                "std_ID": std_id,
                "line_A": line_a,
                "line_B": line_b,
                "n": n,
                "mean_area": np.nan,
                "std_area": np.nan,
                "cv_area": np.nan,
                "pearson_r": corr_v,
                "cv_A": cv_a,
                "cv_B": cv_b,
                "sum_mean": sum_mean,
                "sum_std": sum_std,
                "sum_cv": sum_cv,
                "is_anti_corr": is_anti_corr,
                "is_sum_locked_pair": is_sum_locked,
                "is_split_pair": is_split_pair,
                "signal": "|".join(signal_tokens),
            }
        )

    return pd.DataFrame(summary_rows)


def _select_reference_peak(free_lines: Sequence[str]) -> str:
    """Choose a reference peak to compute relative ratios (e.g., La1, Ka1, Ma1)."""
    reference_suffixes = ("Ka1", "La1", "Ma1", "Mz1", "Ll")

    def suffix_of(line: str) -> str:
        return line.rsplit("_", 1)[-1]

    def element_of(line: str) -> str:
        return line.rsplit("_", 1)[0]

    # Filter down to actual reference candidates.
    candidates = [line for line in free_lines if suffix_of(line) in reference_suffixes]

    if len(candidates) == 0:
        raise ValueError(f"No reference line found among: {list(free_lines)}")

    if len(candidates) == 1:
        return candidates[0]

    # Beyond this point we have multiple candidates; only same-element pairs
    # are resolvable.
    elements = {element_of(line) for line in candidates}
    if len(elements) > 1:
        raise ValueError(
            f"Multiple elements present among reference lines: {list(candidates)}"
        )

    suffixes = {suffix_of(line): line for line in candidates}

    resolvable_pairs = {
        frozenset(("Mz1", "Ma1")): "Ma1",
        frozenset(("La1", "Ll")): "La1",
    }

    if len(candidates) == 2 and len(suffixes) == 2:
        key = frozenset(suffixes)
        if key in resolvable_pairs:
            return suffixes[resolvable_pairs[key]]

    raise ValueError(
        f"Multiple reference lines present and not a resolvable pair: {list(candidates)}"
    )


def _line_key(el_line: str) -> str:
    return el_line.split("_", 1)[1] if "_" in el_line else el_line


def _build_spectrum_column_labels(spectrum_ids: Sequence[Any]) -> List[str]:
    """Build per-spectrum column labels linked to ledger spectrum IDs."""
    return [f"Spectrum #{spectrum_id}" for spectrum_id in spectrum_ids]


def _parse_spectrum_id_from_column(col_name: str) -> Any:
    """Parse spectrum identifier from a column label like 'Spectrum #12'."""
    prefix = "Spectrum #"
    if not str(col_name).startswith(prefix):
        return None
    raw = str(col_name)[len(prefix):]
    return _normalize_spectrum_id(raw)

def _append_correction_factor_column(
    table_df: pd.DataFrame,
    line_order: Sequence[str],
    abs_factors: Dict[str, float],
) -> pd.DataFrame:
    """Insert a 'Correction factor' column matching each line row in the table."""
    out = table_df.copy()
    factor_col = [float(abs_factors.get(el_line, 1.0)) for el_line in line_order]
    # Place it just after the 'Original w' column when present, else at the end.
    if "Original w" in out.columns:
        loc = out.columns.get_loc("Original w")
        if not isinstance(loc, int):
            # Non-unique/duplicated column label; fall back to appending at the end.
            insert_at = len(out.columns)
        else:
            insert_at = loc + 1
    else:
        insert_at = len(out.columns)
    out.insert(insert_at, "Correction factor", factor_col)
    return out

def _resolve_output_dir(sample_dir: Path, std_id: str, ref_peak: str, free_area_el_lines: Sequence[str]) -> Path:
    """Resolve and create per-element calibration output directory."""

    output_dir = sample_dir / f"{ref_peak}_group_area_calib"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _build_output_paths(output_dir: Path, std_id: str, ref_peak: str, output_file_suffix: str) -> Dict[str, Path]:
    """Build deterministic output paths for all calibration CSV files."""
    return {
        "raw": output_dir / f"{std_id}_{ref_peak}_group_weights_raw{output_file_suffix}.csv",
        "areas_table": output_dir / f"{std_id}_{ref_peak}_group_weights_areas{output_file_suffix}.csv",
        "ratios_table": output_dir / f"{std_id}_{ref_peak}_group_weights_ratios{output_file_suffix}.csv",
        "ratios_abs_corrected": output_dir
        / f"{std_id}_{ref_peak}_group_weights_ratios_abs_corrected{output_file_suffix}.csv",
        "correlated_adjustments": output_dir
        / f"{std_id}_{ref_peak}_group_weights_correlated_adjustments{output_file_suffix}.csv",
        "correlated_diagnostics": output_dir
        / f"{std_id}_{ref_peak}_group_weights_correlated_diagnostics{output_file_suffix}.csv",
        "summary": output_dir / f"{std_id}_{ref_peak}_group_weights_summary{output_file_suffix}.csv",
    }

def _finalize_calibration_outputs(
    std_id: str,
    raw_df: pd.DataFrame,
    free_area_el_lines: Sequence[str],
    areas_by_line: pd.DataFrame,
    spectrum_ids_used: Sequence[Any],
    output_dir: Path,
    output_file_suffix: str,
    sum_cv_threshold: float,
    n_sig_figs: int,
    ref_peak: str,
    abs_factors: Optional[Dict[str, float]] = None,
    manual_groups: Optional[Sequence[Sequence[str]]] = None,
    log_verb: str = "Saved",
) -> Dict[str, Any]:
    """
    Shared post-fit pipeline: derive ratios, build tables, round, write CSVs, and
    assemble the return dict. Used by both calibrate and recompute entry points.
    """
    ref_vals = pd.to_numeric(areas_by_line[ref_peak], errors="coerce")

    # Raw (uncorrected) ratios.
    ratios_by_line = areas_by_line.div(ref_vals.replace(0, np.nan), axis=0)

    # Absorption-corrected ratios (separate frame; do not mutate the raw one).
    if abs_factors:
        ratios_corrected_by_line = _apply_absorption_factors(
            ratios_by_line, free_area_el_lines, abs_factors
        )
    else:
        ratios_corrected_by_line = ratios_by_line.copy()

    original_weights = _get_original_weights(free_area_el_lines)
    line_energies = _get_line_energies(free_area_el_lines)

    areas_table_df = _build_workbook_like_table(
        values_by_line=areas_by_line,
        line_order=free_area_el_lines,
        spectrum_ids=spectrum_ids_used,
        original_weights=original_weights,
        line_energies=line_energies,
    )
    ratios_table_df = _build_workbook_like_table(
        values_by_line=ratios_by_line,
        line_order=free_area_el_lines,
        spectrum_ids=spectrum_ids_used,
        original_weights=original_weights,
        line_energies=line_energies,
    )
    ratios_table_df.insert(0, "Reference peak", ref_peak)

    ratios_corrected_table_df = _build_workbook_like_table(
        values_by_line=ratios_corrected_by_line,
        line_order=free_area_el_lines,
        spectrum_ids=spectrum_ids_used,
        original_weights=original_weights,
        line_energies=line_energies,
    )
    ratios_corrected_table_df = _append_correction_factor_column(
        ratios_corrected_table_df, free_area_el_lines, abs_factors or {}
    )
    ratios_corrected_table_df.insert(0, "Reference peak", ref_peak)

    summary_df = _build_summary_table(
        std_id=std_id,
        areas_df=raw_df,
        free_lines=free_area_el_lines,
        line_energies=line_energies,
        sum_cv_threshold=sum_cv_threshold,
    )

    recommendations_df, diagnostics_df = _build_correlated_adjustments_table(
        summary_df=summary_df,
        ratios_by_line=ratios_corrected_by_line,
        areas_by_line=areas_by_line,
        line_energies=line_energies,
        original_weights=original_weights,
        ref_peak=ref_peak,
        sum_cv_threshold=sum_cv_threshold,
        spectrum_ids=spectrum_ids_used,
        manual_groups=manual_groups,
    )

    # Reduce exported numeric precision for easier manual use in calibration workbooks.
    raw_df_out = _round_dataframe_sig_figs(raw_df, n_sig_figs=n_sig_figs)
    areas_table_df_out = _round_dataframe_sig_figs(areas_table_df, n_sig_figs=n_sig_figs)
    ratios_table_df_out = _round_dataframe_sig_figs(ratios_table_df, n_sig_figs=n_sig_figs)
    ratios_corrected_table_df_out = _round_dataframe_sig_figs(
        ratios_corrected_table_df, n_sig_figs=n_sig_figs
    )
    recommendations_df_out = _round_dataframe_sig_figs(recommendations_df, n_sig_figs=n_sig_figs)
    diagnostics_df_out = _round_dataframe_sig_figs(diagnostics_df, n_sig_figs=n_sig_figs)
    summary_df_out = _round_dataframe_sig_figs(summary_df, n_sig_figs=n_sig_figs)

    output_paths = _build_output_paths(output_dir, std_id, ref_peak, output_file_suffix)

    raw_df_out.to_csv(output_paths["raw"], index=False)
    areas_table_df_out.to_csv(output_paths["areas_table"], index=False)
    ratios_table_df_out.to_csv(output_paths["ratios_table"], index=False)
    ratios_corrected_table_df_out.to_csv(output_paths["ratios_abs_corrected"], index=False)
    recommendations_df_out.to_csv(output_paths["correlated_adjustments"], index=False)
    diagnostics_df_out.to_csv(output_paths["correlated_diagnostics"], index=False)
    summary_df_out.to_csv(output_paths["summary"], index=False)

    print_double_separator()
    logging.info("%s raw peak areas to: %s", log_verb, output_paths["raw"])
    logging.info("%s workbook-like areas table to: %s", log_verb, output_paths["areas_table"])
    logging.info("%s workbook-like ratio table to: %s", log_verb, output_paths["ratios_table"])
    logging.info(
        "%s absorption-corrected ratio table to: %s",
        log_verb,
        output_paths["ratios_abs_corrected"],
    )
    logging.info(
        "%s correlated adjustments to: %s",
        log_verb,
        output_paths["correlated_adjustments"],
    )
    logging.info(
        "%s correlated diagnostics to: %s",
        log_verb,
        output_paths["correlated_diagnostics"],
    )
    logging.info("%s calibration summary to: %s", log_verb, output_paths["summary"])
    print_double_separator()

    return {
        "raw": raw_df_out,
        "areas_table": areas_table_df_out,
        "ratios_table": ratios_table_df_out,
        "ratios_abs_corrected": ratios_corrected_table_df_out,
        "correlated_adjustments": recommendations_df_out,
        "correlated_diagnostics": diagnostics_df_out,
        "summary": summary_df_out,
    }

def _build_workbook_like_table(
    values_by_line: pd.DataFrame,
    line_order: Sequence[str],
    spectrum_ids: Sequence[Any],
    original_weights: Optional[Dict[str, float]] = None,
    line_energies: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Build a line-by-row table with Spectrum # columns and summary stats."""
    if original_weights is None:
        original_weights = {}
    if line_energies is None:
        line_energies = {}

    spectrum_cols = _build_spectrum_column_labels(spectrum_ids)

    rows: List[Dict[str, Any]] = []
    for el_line in line_order:
        vals = pd.to_numeric(values_by_line[el_line], errors="coerce")
        row: Dict[str, Any] = {
            "Line": _line_key(el_line),
            "Energy (keV)": line_energies.get(el_line, np.nan),
            "Original w": original_weights.get(el_line, np.nan),
        }

        for idx, col_name in enumerate(spectrum_cols):
            row[col_name] = float(vals.iloc[idx]) if idx < len(vals) and pd.notna(vals.iloc[idx]) else np.nan

        mean_v = float(vals.mean()) if vals.notna().any() else np.nan
        std_v = float(vals.std(ddof=1)) if vals.notna().sum() > 1 else np.nan
        er_v = float(std_v / mean_v * 100) if np.isfinite(mean_v) and mean_v != 0 and np.isfinite(std_v) else np.nan

        row["Mean"] = mean_v
        row["Stdev"] = std_v
        row["Er (%)"] = er_v
        rows.append(row)

    return pd.DataFrame(rows)


def _get_original_weights(free_lines: Sequence[str]) -> Dict[str, float]:
    """Map element-line labels (e.g., Pr_Ma1) to tabulated X-ray line weights when available."""
    weights: Dict[str, float] = {}
    by_el: Dict[str, List[str]] = {}
    for el_line in free_lines:
        if "_" not in el_line:
            continue
        el, _ = el_line.split("_", 1)
        by_el.setdefault(el, []).append(el_line)

    for el, lines in by_el.items():
        try:
            line_data = get_el_xray_lines(el)
        except Exception:
            continue

        for el_line in lines:
            _, line = el_line.split("_", 1)
            w = line_data.get(line, {}).get("weight")
            if w is not None:
                weights[el_line] = float(w)

    return weights


def _get_line_energies(free_lines: Sequence[str]) -> Dict[str, float]:
    """Map element-line labels (e.g., Pr_Ma1) to tabulated X-ray line energies (keV)."""
    energies: Dict[str, float] = {}
    by_el: Dict[str, List[str]] = {}
    for el_line in free_lines:
        if "_" not in el_line:
            continue
        el, _ = el_line.split("_", 1)
        by_el.setdefault(el, []).append(el_line)

    for el, lines in by_el.items():
        try:
            line_data = get_el_xray_lines(el)
        except Exception:
            continue

        for el_line in lines:
            _, line = el_line.split("_", 1)
            en = line_data.get(line, {}).get("energy (keV)")
            if en is not None:
                energies[el_line] = float(en)

    return energies


def _are_adjacent_by_energy(line_a: str, line_b: str, line_energies: Dict[str, float]) -> bool:
    """Return True when two lines are adjacent in energy ordering among same-element free lines."""
    if "_" not in line_a or "_" not in line_b:
        return False
    el_a = line_a.split("_", 1)[0]
    el_b = line_b.split("_", 1)[0]
    if el_a != el_b:
        return False

    same_el_lines = sorted(
        [line for line in line_energies if line.startswith(f"{el_a}_") and np.isfinite(line_energies[line])],
        key=lambda line: line_energies[line],
    )
    if line_a not in same_el_lines or line_b not in same_el_lines:
        return False

    idx_a = same_el_lines.index(line_a)
    idx_b = same_el_lines.index(line_b)
    return abs(idx_a - idx_b) == 1

def _build_correlated_adjustments_table(
    summary_df: pd.DataFrame,
    ratios_by_line: pd.DataFrame,
    areas_by_line: pd.DataFrame,
    line_energies: Dict[str, float],
    original_weights: Dict[str, float],
    ref_peak: str,
    sum_cv_threshold: float,
    spectrum_ids: Sequence[Any],
    manual_groups: Optional[Sequence[Sequence[str]]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build correlated-peak adjustment tables for sum-locked groups.

    Returns
    -------
    (recommendations_df, group_diagnostics_df)
        recommendations_df : one row per peak, the actionable redistribution.
        group_diagnostics_df : one row per group, the supporting statistics.
    """
    manual_groups = [list(g) for g in (manual_groups or []) if len(g) >= 2]
    manually_locked_lines = {line for g in manual_groups for line in g}

    rec_columns = [
        "Group", "Lines in group", "Reference peak",
        "Line", "Energy (keV)",
        "Original w", "Original frac in group",
        "Recommended ratio", "n",
    ]
    diag_columns = [
        "Group", "Lines in group", "n_lines", "n",
        "Area sum mean", "Area sum stdev", "Area sum er (%)",
        "Meas total ratio mean", "Meas total ratio stdev", "Meas total ratio er (%)",
        "Source",
    ]

    pair_rows = summary_df[
        (summary_df["record_type"] == "pair_correlation")
        & (summary_df["is_split_pair"] == True)  # noqa: E712
    ]

    if not pair_rows.empty and manually_locked_lines:
        keep_mask = pair_rows.apply(
            lambda r: (str(r["line_A"]) not in manually_locked_lines)
            and (str(r["line_B"]) not in manually_locked_lines),
            axis=1,
        )
        pair_rows = pair_rows[keep_mask]

    if not pair_rows.empty:
        adjacency_mask = pair_rows.apply(
            lambda r: _are_adjacent_by_energy(
                str(r["line_A"]), str(r["line_B"]), line_energies
            ),
            axis=1,
        )
        pair_rows = pair_rows[adjacency_mask]

    if pair_rows.empty and not manual_groups:
        return pd.DataFrame(columns=rec_columns), pd.DataFrame(columns=diag_columns)

    # --- build connected groups (unchanged logic) ---
    adjacency: Dict[str, set] = {}
    for _, row in pair_rows.iterrows():
        a, b = str(row["line_A"]), str(row["line_B"])
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    connected_groups: List[List[str]] = []
    visited = set()
    for node in adjacency:
        if node in visited:
            continue
        stack, comp = [node], []
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            stack.extend(n for n in adjacency.get(cur, set()) if n not in visited)
        if len(comp) >= 2:
            connected_groups.append(sorted(comp))

    # Track provenance so the diagnostics table can show it.
    group_source: Dict[frozenset, str] = {}

    groups: List[List[str]] = []
    for comp in connected_groups:
        if len(comp) == 2:
            groups.append(comp)
            group_source[frozenset(comp)] = "auto-pair"
            continue
        if len(comp) > 3:
            sub_pairs = pair_rows[
                pair_rows["line_A"].isin(comp) & pair_rows["line_B"].isin(comp)
            ][["line_A", "line_B"]].drop_duplicates()
            for _, p in sub_pairs.iterrows():
                g = sorted([str(p["line_A"]), str(p["line_B"])])
                groups.append(g)
                group_source[frozenset(g)] = "auto-pair (split from >3)"
            continue

        comp_total = areas_by_line[comp].apply(pd.to_numeric, errors="coerce").sum(axis=1)
        comp_sum_cv = np.nan
        if comp_total.notna().sum() > 1:
            cmean = float(comp_total.mean())
            if np.isfinite(cmean) and cmean != 0:
                comp_sum_cv = float(comp_total.std(ddof=1)) / cmean

        if np.isfinite(comp_sum_cv) and comp_sum_cv <= sum_cv_threshold:
            groups.append(comp)
            group_source[frozenset(comp)] = "auto-triple"
        else:
            sub_pairs = pair_rows[
                pair_rows["line_A"].isin(comp) & pair_rows["line_B"].isin(comp)
            ][["line_A", "line_B"]].drop_duplicates()
            for _, p in sub_pairs.iterrows():
                g = sorted([str(p["line_A"]), str(p["line_B"])])
                groups.append(g)
                group_source[frozenset(g)] = "auto-pair (triple rejected)"

    for g in manual_groups:
        group_source[frozenset(g)] = "manual"
    groups = [list(g) for g in manual_groups] + groups

    unique_groups, seen = [], set()
    for g in groups:
        key = frozenset(g)
        if key in seen:
            continue
        seen.add(key)
        unique_groups.append(g)
    groups = unique_groups

    rec_rows: List[Dict[str, Any]] = []
    diag_rows: List[Dict[str, Any]] = []

    for g_idx, group_lines in enumerate(groups, start=1):
        missing = [
            line for line in group_lines
            if line not in ratios_by_line.columns or line not in areas_by_line.columns
        ]
        if missing:
            logging.warning("Skipping group %s: lines not in fitted data: %s",
                            group_lines, missing)
            continue

        area_total = areas_by_line[group_lines].apply(pd.to_numeric, errors="coerce").sum(axis=1)
        ratio_total = ratios_by_line[group_lines].apply(pd.to_numeric, errors="coerce").sum(axis=1)
        n = int(ratio_total.notna().sum())

        area_mean = float(area_total.mean()) if area_total.notna().any() else np.nan
        area_std = float(area_total.std(ddof=1)) if area_total.notna().sum() > 1 else np.nan
        area_er = (float(area_std / area_mean * 100)
                   if np.isfinite(area_mean) and area_mean != 0 and np.isfinite(area_std)
                   else np.nan)

        mean_total = float(ratio_total.mean()) if ratio_total.notna().any() else np.nan
        std_total = float(ratio_total.std(ddof=1)) if ratio_total.notna().sum() > 1 else np.nan
        er_total = (float(std_total / mean_total * 100)
                    if np.isfinite(mean_total) and mean_total != 0 and np.isfinite(std_total)
                    else np.nan)

        group_weights = {line: original_weights.get(line, np.nan) for line in group_lines}
        valid_w = {line: w for line, w in group_weights.items() if np.isfinite(w)}
        if len(valid_w) == len(group_lines) and sum(valid_w.values()) > 0:
            w_sum = float(sum(valid_w.values()))
            fracs = {line: float(valid_w[line] / w_sum) for line in group_lines}
        else:
            means = pd.to_numeric(ratios_by_line[group_lines].mean(axis=0), errors="coerce")
            mean_sum = float(means.sum()) if means.notna().any() else np.nan
            fracs = ({line: float(means[line] / mean_sum) for line in group_lines}
                     if np.isfinite(mean_sum) and mean_sum > 0
                     else {line: np.nan for line in group_lines})

        group_label = f"Group {g_idx}"
        group_lines_str = " + ".join(_line_key(line) for line in group_lines)

        diag_rows.append({
            "Group": group_label,
            "Lines in group": group_lines_str,
            "n_lines": len(group_lines),
            "n": n,
            "Area sum mean": area_mean,
            "Area sum stdev": area_std,
            "Area sum er (%)": area_er,
            "Meas total ratio mean": mean_total,
            "Meas total ratio stdev": std_total,
            "Meas total ratio er (%)": er_total,
            "Source": group_source.get(frozenset(group_lines), "auto"),
        })

        for line in group_lines:
            frac = fracs.get(line, np.nan)
            rec = (float(mean_total * frac)
                   if np.isfinite(mean_total) and np.isfinite(frac) else np.nan)
            rec_rows.append({
                "Group": group_label,
                "Lines in group": group_lines_str,
                "Reference peak": ref_peak,
                "Line": _line_key(line),
                "Energy (keV)": line_energies.get(line, np.nan),
                "Original w": group_weights.get(line, np.nan),
                "Original frac in group": frac,
                "Recommended ratio": rec,
                "n": n,
            })

    rec_df = pd.DataFrame(rec_rows, columns=rec_columns)
    diag_df = pd.DataFrame(diag_rows, columns=diag_columns)
    return rec_df, diag_df

def _round_to_sig_figs(val: Any, n_sig_figs: int) -> Any:
    """Round numeric scalars to a fixed number of significant figures."""
    if not isinstance(val, (int, float, np.integer, np.floating)):
        return val

    if not np.isfinite(val) or val == 0:
        return val

    return float(f"{float(val):.{n_sig_figs}g}")


def _round_dataframe_sig_figs(df: pd.DataFrame, n_sig_figs: int) -> pd.DataFrame:
    """Apply significant-figure rounding only to numeric columns."""
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(lambda x: _round_to_sig_figs(x, n_sig_figs))
    return out

def _load_ledger_for_standard(
    std_ID: str,
    results_path: Optional[str],
) -> Tuple[Path, Path, Path, Any]:
    """Resolve results root, sample dir, ledger path, and load the ledger."""
    results_root = _resolve_results_root(results_path)
    sample_dir = Path(get_sample_dir(str(results_root), std_ID))
    ledger_path = sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
    ledger = load_sample_ledger(str(ledger_path))
    return results_root, sample_dir, ledger_path, ledger

def calibrate_peak_weights_experimental_standard(
    std_ID: str,
    free_area_el_lines: Sequence[str],
    spectrum_IDs_to_fit: Optional[Sequence[Any]] = None,
    results_path: Optional[str] = None,
    spectrum_lims: Optional[Tuple[int, int]] = None,
    els_substrate: Optional[List[str]] = None,
    is_particle: bool = True,
    use_instrument_background: bool = False,
    fit_tol: float = 1e-4,
    sum_cv_threshold: float = 0.12,
    n_sig_figs: int = 5,
    output_file_suffix: str = "",
) -> Dict[str, Any]:
    """
    Fit spectra of one standard and export free-peak area calibration tables.

    Returns
    -------
    dict
        Contains output paths and generated DataFrames with keys:
        raw_path, summary_path, raw_df, summary_df
    """
    if len(free_area_el_lines) == 0:
        raise ValueError(
            "free_area_el_lines is empty. Provide at least one line to calibrate peak weights."
        )
    
    #Remove potential duplicate lines
    free_area_el_lines = list(dict.fromkeys(free_area_el_lines))

    results_root, sample_dir, _, ledger = _load_ledger_for_standard(std_ID, results_path)
    spectrum_ids = _resolve_spectrum_ids(ledger, spectrum_IDs_to_fit)

    print_double_separator(
        message=f"Calibrating peak weights for '{std_ID}' using {len(spectrum_ids)} spectra"
    )

    raw_rows: List[Dict[str, Any]] = []
    latest_valid_quantifier = None
    for spectrum_id in spectrum_ids:
        quantifier = fit_and_quantify_spectrum_from_ledger(
            sample_ID=std_ID,
            spectrum_ID=spectrum_id,
            ledger=ledger,
            is_standard=True,
            spectrum_lims=spectrum_lims,
            results_path=str(results_root),
            use_instrument_background=use_instrument_background,
            quantify_plot=False,
            plot_signal=False,
            zoom_plot=False,
            line_to_plot="",
            fit_tol=fit_tol,
            is_particle=is_particle,
            max_undetectable_w_fr=0,
            force_single_iteration=False,
            interrupt_fits_bad_spectra=False,
            els_substrate=els_substrate,
            print_results=False,
            quant_verbose=False,
            fitting_verbose=False,
            free_area_el_lines=list(free_area_el_lines),
        )

        if quantifier is None:
            logging.warning("Spectrum %s could not be fitted and will be skipped.", spectrum_id)
            continue
        
        latest_valid_quantifier = quantifier
        raw_rows.append(_collect_fit_row(std_ID, spectrum_id, quantifier, free_area_el_lines))

    if len(raw_rows) == 0:
        raise RuntimeError("No spectra were successfully fitted. No calibration output generated.")

    raw_df = pd.DataFrame(raw_rows)
    spectrum_ids_used = raw_df["spectrum_ID"].tolist()
    areas_by_line = raw_df[list(free_area_el_lines)].copy()

    ref_peak = _select_reference_peak(free_area_el_lines)

    output_dir = _resolve_output_dir(sample_dir, std_ID, ref_peak, free_area_el_lines)

    abs_factors = _absorption_factor_per_line(
        free_area_el_lines=free_area_el_lines,
        quantifier=latest_valid_quantifier,
        ref_peak = ref_peak,
        direction="to_pure",
    )

    return _finalize_calibration_outputs(
        std_id=std_ID,
        raw_df=raw_df,
        free_area_el_lines=free_area_el_lines,
        ref_peak=ref_peak,
        areas_by_line=areas_by_line,
        spectrum_ids_used=spectrum_ids_used,
        output_dir=output_dir,
        output_file_suffix=output_file_suffix,
        sum_cv_threshold=sum_cv_threshold,
        n_sig_figs=n_sig_figs,
        abs_factors=abs_factors,
        manual_groups=None,
        log_verb="Saved",
    )


def recompute_calibration(
    std_ID: str,
    free_area_el_lines: Sequence[str],
    spectrum_to_ignore: Optional[Sequence[Any]] = None,
    spectrum_IDs_to_fit: Optional[Sequence[Any]] = None,
    manual_groups: Optional[Sequence[Sequence[str]]] = None,
    results_path: Optional[str] = None,
    sum_cv_threshold: float = 0.12,
    n_sig_figs: int = 5,
) -> Dict[str, Any]:
    """
    Recompute calibration outputs excluding specific spectra, without refitting.

    Absorption-correction factors are not recomputed (no refit occurs); they are
    loaded from the 'ratios_abs_corrected' table produced by the original calibration.

    Parameters
    ----------
    spectrum_to_ignore
        Spectrum IDs to exclude before fitting. Matching is normalized so values like
        `3`, `3.0`, and "3" are treated as the same ID.
    manual_groups
        Optional user-defined groups of element-lines, each given as a sequence of
        element-line labels written like the entries of `free_area_el_lines`
        (e.g. ``[["Pr_Mz1", "Pr_Ma1"], ["La_Lb1", "La_Lb2"]]``). Each provided group
        is taken verbatim as a Group in the correlated-peak adjustments table, and any
        line appearing in a manual group is excluded from the automated detection of
        sum-locked groups, so those peaks cannot be reassigned automatically.
    """
    if len(free_area_el_lines) == 0:
        raise ValueError(
            "free_area_el_lines is empty. Provide at least one line to recalibrate peak weights."
        )
    # Remove potential duplicate lines.
    free_area_el_lines = list(dict.fromkeys(free_area_el_lines))
    ref_peak = _select_reference_peak(free_area_el_lines)

    _, sample_dir, _, ledger = _load_ledger_for_standard(std_ID, results_path)
    output_dir = _resolve_output_dir(sample_dir, std_ID, ref_peak, free_area_el_lines)

    # Normalize and validate manual groups against the available free lines.
    free_lines_set = set(free_area_el_lines)
    normalized_manual_groups: List[List[str]] = []
    for group in (manual_groups or []):
        group_lines = list(group)
        if len(group_lines) < 2:
            raise ValueError(
                f"Each manual group must contain at least two lines. Got: {group_lines}"
            )
        unknown = [line for line in group_lines if line not in free_lines_set]
        if unknown:
            raise ValueError(
                "Manual group references lines not present in free_area_el_lines: "
                f"{unknown}. Available lines: {sorted(free_lines_set)}"
            )
        normalized_manual_groups.append(group_lines)

    source_paths = _build_output_paths(output_dir, std_ID, ref_peak, output_file_suffix="")
    source_raw_path = source_paths["raw"]
    source_areas_table_path = source_paths["areas_table"]
    if not source_areas_table_path.exists():
        raise FileNotFoundError(
            "Areas table not found for recompute. Run a normal calibration first: "
            f"{source_areas_table_path}"
        )

    # Load absorption-correction factors from the original calibration's corrected
    # ratios table (recompute does not refit, so factors cannot be regenerated here).
    abs_factors = _load_abs_factors_from_table(
        ratios_abs_corrected_path=source_paths["ratios_abs_corrected"],
        free_area_el_lines=free_area_el_lines,
    )

    if len(free_area_el_lines) == 0:
        raise ValueError(
            "free_area_el_lines is empty. Provide at least one line to recalibrate peak weights."
        )
    #Remove potential duplicate lines
    free_area_el_lines = list(dict.fromkeys(free_area_el_lines))

    # Normalize and validate manual groups against the available free lines.
    free_lines_set = set(free_area_el_lines)
    normalized_manual_groups: List[List[str]] = []
    for group in (manual_groups or []):
        group_lines = list(group)
        if len(group_lines) < 2:
            raise ValueError(
                f"Each manual group must contain at least two lines. Got: {group_lines}"
            )
        unknown = [line for line in group_lines if line not in free_lines_set]
        if unknown:
            raise ValueError(
                "Manual group references lines not present in free_area_el_lines: "
                f"{unknown}. Available lines: {sorted(free_lines_set)}"
            )
        normalized_manual_groups.append(group_lines)

    source_paths = _build_output_paths(output_dir, std_ID, ref_peak, output_file_suffix="")
    source_raw_path = source_paths["raw"]
    source_areas_table_path = source_paths["areas_table"]
    if not source_areas_table_path.exists():
        raise FileNotFoundError(
            "Areas table not found for recompute. Run a normal calibration first: "
            f"{source_areas_table_path}"
        )

    areas_table_in = pd.read_csv(source_areas_table_path)
    spectrum_cols_all = [col for col in areas_table_in.columns if str(col).startswith("Spectrum #")]
    if len(spectrum_cols_all) == 0:
        raise ValueError(
            "Could not find any 'Spectrum #' columns in areas table. "
            f"File: {source_areas_table_path}"
        )

    spectrum_ids_all = [_parse_spectrum_id_from_column(col) for col in spectrum_cols_all]

    selected_norm: Optional[set] = None
    if spectrum_IDs_to_fit is not None:
        selected_norm = {_normalize_spectrum_id(spectrum_id) for spectrum_id in spectrum_IDs_to_fit}

    ignore_norm = {_normalize_spectrum_id(spectrum_id) for spectrum_id in (spectrum_to_ignore or [])}

    spectrum_cols_used: List[str] = []
    spectrum_ids_used: List[Any] = []
    for col_name, spectrum_id in zip(spectrum_cols_all, spectrum_ids_all):
        if selected_norm is not None and spectrum_id not in selected_norm:
            continue
        if spectrum_id in ignore_norm:
            continue
        spectrum_cols_used.append(col_name)
        spectrum_ids_used.append(spectrum_id)

    if len(spectrum_cols_used) == 0:
        raise ValueError(
            "After applying spectrum_IDs_to_fit and spectrum_to_ignore, no spectra remain in areas table."
        )

    if "Line" not in areas_table_in.columns:
        raise ValueError(
            "Areas table is missing required 'Line' column. "
            f"File: {source_areas_table_path}"
        )

    line_key_series = areas_table_in["Line"].astype(str)
    used_row_idx: set = set()
    areas_data: Dict[str, List[float]] = {}
    for el_line in free_area_el_lines:
        line_key = _line_key(el_line)
        candidate_idx = [idx for idx in line_key_series[line_key_series == line_key].index if idx not in used_row_idx]
        if len(candidate_idx) == 0:
            raise ValueError(
                f"Line '{line_key}' (from '{el_line}') not found in source areas table: {source_areas_table_path}"
            )
        row_idx = candidate_idx[0]
        used_row_idx.add(row_idx)

        vals = pd.to_numeric(areas_table_in.loc[row_idx, spectrum_cols_used], errors="coerce")
        areas_data[el_line] = [float(v) if pd.notna(v) else np.nan for v in vals.to_list()]

    areas_by_line = pd.DataFrame(areas_data)

    redchi_by_spectrum: Dict[Any, float] = {}
    if source_raw_path.exists():
        try:
            raw_in = pd.read_csv(source_raw_path)
            if "spectrum_ID" in raw_in.columns and "redchi" in raw_in.columns:
                for _, row in raw_in.iterrows():
                    sid_norm = _normalize_spectrum_id(row["spectrum_ID"])
                    redchi_val = pd.to_numeric(row["redchi"], errors="coerce")
                    if sid_norm not in redchi_by_spectrum:
                        redchi_by_spectrum[sid_norm] = float(redchi_val) if pd.notna(redchi_val) else np.nan
        except Exception:
            redchi_by_spectrum = {}

    raw_data: Dict[str, Any] = {
        "std_ID": [std_ID] * len(spectrum_ids_used),
        "spectrum_ID": spectrum_ids_used,
        "redchi": [redchi_by_spectrum.get(_normalize_spectrum_id(sid), np.nan) for sid in spectrum_ids_used],
    }
    for el_line in free_area_el_lines:
        raw_data[el_line] = areas_by_line[el_line].tolist()
    raw_df = pd.DataFrame(raw_data)

    return _finalize_calibration_outputs(
        std_id=std_ID,
        raw_df=raw_df,
        free_area_el_lines=free_area_el_lines,
        ref_peak=ref_peak,
        areas_by_line=areas_by_line,
        spectrum_ids_used=spectrum_ids_used,
        output_dir=output_dir,
        output_file_suffix="_recompute",
        sum_cv_threshold=sum_cv_threshold,
        n_sig_figs=n_sig_figs,
        abs_factors=abs_factors,
        manual_groups=normalized_manual_groups,
        log_verb="Recomputed",
    )