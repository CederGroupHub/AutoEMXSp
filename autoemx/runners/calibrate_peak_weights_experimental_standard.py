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
from datetime import datetime
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
from autoemx.utils.helper import get_sample_dir, print_double_separator
import autoemx.utils.constants as cnst


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["calibrate_peak_weights_experimental_standard"]


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
    abs_corr_threshold: float,
    sum_locked_neg_corr_threshold: float,
    sum_cv_threshold: float,
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
                "sum_mean": np.nan,
                "sum_std": np.nan,
                "sum_cv": np.nan,
                "is_high_corr": False,
                "is_sum_locked_pair": False,
                "signal": "",
            }
        )

    for line_a, line_b in combinations(free_lines, 2):
        pair_df = areas_df[[line_a, line_b]].dropna()
        n = int(pair_df.shape[0])

        if n < 3:
            corr_v = np.nan
            sum_mean = np.nan
            sum_std = np.nan
            sum_cv = np.nan
            is_high_corr = False
            is_sum_locked = False
        else:
            corr_v = float(pair_df[line_a].corr(pair_df[line_b]))
            pair_sum = pair_df[line_a] + pair_df[line_b]
            sum_mean = float(pair_sum.mean())
            sum_std = float(pair_sum.std(ddof=1)) if n > 1 else np.nan
            sum_cv = float(sum_std / sum_mean) if np.isfinite(sum_mean) and sum_mean != 0 else np.nan

            is_high_corr = bool(np.isfinite(corr_v) and abs(corr_v) >= abs_corr_threshold)
            is_sum_locked = bool(
                np.isfinite(corr_v)
                and np.isfinite(sum_cv)
                and corr_v <= sum_locked_neg_corr_threshold
                and sum_cv <= sum_cv_threshold
            )

        signal_tokens = []
        if is_high_corr:
            signal_tokens.append("HIGH_CORRELATION")
        if is_sum_locked:
            signal_tokens.append("SUM_LOCKED_PAIR")

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
                "sum_mean": sum_mean,
                "sum_std": sum_std,
                "sum_cv": sum_cv,
                "is_high_corr": is_high_corr,
                "is_sum_locked_pair": is_sum_locked,
                "signal": "|".join(signal_tokens),
            }
        )

    return pd.DataFrame(summary_rows)


def _select_reference_peak(free_lines: Sequence[str], areas_df: pd.DataFrame) -> str:
    """Choose a reference peak to compute relative ratios (e.g., La1, Ka1, Ma1)."""
    if len(free_lines) == 1:
        return free_lines[0]

    preferred_suffixes = ("Ka1", "La1", "Ma1", "Mz1", "Ll")
    for suffix in preferred_suffixes:
        candidates = [line for line in free_lines if line.endswith(f"_{suffix}")]
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            means = pd.to_numeric(areas_df[candidates].mean(axis=0), errors="coerce")
            return str(means.idxmax())

    means = pd.to_numeric(areas_df[list(free_lines)].mean(axis=0), errors="coerce")
    return str(means.idxmax())


def _line_key(el_line: str) -> str:
    return el_line.split("_", 1)[1] if "_" in el_line else el_line


def _build_workbook_like_table(
    values_by_line: pd.DataFrame,
    line_order: Sequence[str],
    original_weights: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """Build a line-by-row table with Meas columns and summary stats like calibration workbooks."""
    if original_weights is None:
        original_weights = {}

    n_meas = values_by_line.shape[0]
    meas_cols = [f"Meas {i}" for i in range(1, n_meas + 1)]

    rows: List[Dict[str, Any]] = []
    for el_line in line_order:
        vals = pd.to_numeric(values_by_line[el_line], errors="coerce")
        row: Dict[str, Any] = {
            "Line": _line_key(el_line),
            "Original w": original_weights.get(el_line, np.nan),
        }

        for idx, col_name in enumerate(meas_cols):
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


def _build_correlated_adjustments_table(
    summary_df: pd.DataFrame,
    ratios_by_line: pd.DataFrame,
    areas_by_line: pd.DataFrame,
    original_weights: Dict[str, float],
    ref_peak: str,
) -> pd.DataFrame:
    """
    For sum-locked correlated pairs, compute total measured ratio and recommended
    per-peak ratios preserving original relative proportions.
    """
    pair_rows = summary_df[
        (summary_df["record_type"] == "pair_correlation")
        & (summary_df["is_sum_locked_pair"] == True)
    ]

    if pair_rows.empty:
        return pd.DataFrame(
            columns=[
                "Group",
                "Reference peak",
                "Lines in group",
                "Line",
                "Original w",
                "Original frac in group",
                "Area sum mean",
                "Area sum stdev",
                "Area sum er (%)",
                "Meas total ratio mean",
                "Meas total ratio stdev",
                "Meas total ratio er (%)",
                "Recommended ratio",
                "n",
            ]
        )

    # Build connected groups from pairwise sum-locked relations.
    adjacency: Dict[str, set] = {}
    for _, row in pair_rows.iterrows():
        a = str(row["line_A"])
        b = str(row["line_B"])
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    groups: List[List[str]] = []
    visited = set()
    for node in adjacency:
        if node in visited:
            continue
        stack = [node]
        comp: List[str] = []
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp.append(cur)
            stack.extend([n for n in adjacency.get(cur, set()) if n not in visited])
        if len(comp) >= 2:
            groups.append(sorted(comp))

    out_rows: List[Dict[str, Any]] = []

    meas_cols = [f"Meas {i}" for i in range(1, ratios_by_line.shape[0] + 1)]

    for g_idx, group_lines in enumerate(groups, start=1):
        missing_cols = [line for line in group_lines if line not in ratios_by_line.columns or line not in areas_by_line.columns]
        if missing_cols:
            continue

        area_total = areas_by_line[group_lines].apply(pd.to_numeric, errors="coerce").sum(axis=1)
        ratio_total = ratios_by_line[group_lines].apply(pd.to_numeric, errors="coerce").sum(axis=1)

        n = int(ratio_total.notna().sum())

        area_mean = float(area_total.mean()) if area_total.notna().any() else np.nan
        area_std = float(area_total.std(ddof=1)) if area_total.notna().sum() > 1 else np.nan
        area_er = (
            float(area_std / area_mean * 100)
            if np.isfinite(area_mean) and area_mean != 0 and np.isfinite(area_std)
            else np.nan
        )

        mean_total = float(ratio_total.mean()) if ratio_total.notna().any() else np.nan
        std_total = float(ratio_total.std(ddof=1)) if ratio_total.notna().sum() > 1 else np.nan
        er_total = (
            float(std_total / mean_total * 100)
            if np.isfinite(mean_total) and mean_total != 0 and np.isfinite(std_total)
            else np.nan
        )

        # Fractions from original weights when available, fallback to measured mean ratios.
        group_weights = {line: original_weights.get(line, np.nan) for line in group_lines}
        valid_w = {line: w for line, w in group_weights.items() if np.isfinite(w)}

        if len(valid_w) == len(group_lines):
            w_sum = float(sum(valid_w.values()))
            if w_sum > 0:
                fracs = {line: float(valid_w[line] / w_sum) for line in group_lines}
            else:
                fracs = {line: np.nan for line in group_lines}
        else:
            means = pd.to_numeric(ratios_by_line[group_lines].mean(axis=0), errors="coerce")
            mean_sum = float(means.sum()) if means.notna().any() else np.nan
            if np.isfinite(mean_sum) and mean_sum > 0:
                fracs = {line: float(means[line] / mean_sum) for line in group_lines}
            else:
                fracs = {line: np.nan for line in group_lines}

        group_label = f"Group {g_idx}"
        group_lines_str = " + ".join(_line_key(line) for line in group_lines)

        for line in group_lines:
            frac = fracs.get(line, np.nan)
            rec = float(mean_total * frac) if np.isfinite(mean_total) and np.isfinite(frac) else np.nan

            out_row: Dict[str, Any] = {
                "Group": group_label,
                "Reference peak": ref_peak,
                "Lines in group": group_lines_str,
                "Line": _line_key(line),
                "Original w": group_weights.get(line, np.nan),
                "Original frac in group": frac,
                "Area sum mean": area_mean,
                "Area sum stdev": area_std,
                "Area sum er (%)": area_er,
                "Meas total ratio mean": mean_total,
                "Meas total ratio stdev": std_total,
                "Meas total ratio er (%)": er_total,
                "Recommended ratio": rec,
                "n": n,
            }

            for idx, col_name in enumerate(meas_cols):
                out_row[f"Area sum {col_name}"] = (
                    float(area_total.iloc[idx])
                    if idx < len(area_total) and pd.notna(area_total.iloc[idx])
                    else np.nan
                )
                out_row[f"Ratio sum {col_name}"] = (
                    float(ratio_total.iloc[idx])
                    if idx < len(ratio_total) and pd.notna(ratio_total.iloc[idx])
                    else np.nan
                )

            out_rows.append(out_row)

    return pd.DataFrame(out_rows)


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
    abs_corr_threshold: float = 0.90,
    sum_locked_neg_corr_threshold: float = -0.90,
    sum_cv_threshold: float = 0.12,
    n_sig_figs: int = 5,
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

    results_root = _resolve_results_root(results_path)
    sample_dir = Path(get_sample_dir(str(results_root), std_ID))
    ledger_path = sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"

    ledger = load_sample_ledger(str(ledger_path))
    spectrum_ids = _resolve_spectrum_ids(ledger, spectrum_IDs_to_fit)

    print_double_separator(
        message=f"Calibrating peak weights for '{std_ID}' using {len(spectrum_ids)} spectra"
    )

    raw_rows: List[Dict[str, Any]] = []
    for spectrum_id in spectrum_ids:
        logging.info("Fitting standard '%s' spectrum %s", std_ID, spectrum_id)
        quantifier = fit_and_quantify_spectrum_from_ledger(
            sample_ID=std_ID,
            spectrum_ID=spectrum_id,
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

        raw_rows.append(_collect_fit_row(std_ID, spectrum_id, quantifier, free_area_el_lines))

    if len(raw_rows) == 0:
        raise RuntimeError("No spectra were successfully fitted. No calibration output generated.")

    raw_df = pd.DataFrame(raw_rows)
    areas_by_line = raw_df[list(free_area_el_lines)].copy()

    ref_peak = _select_reference_peak(free_area_el_lines, areas_by_line)
    ref_vals = pd.to_numeric(areas_by_line[ref_peak], errors="coerce")
    ratios_by_line = areas_by_line.div(ref_vals.replace(0, np.nan), axis=0)

    original_weights = _get_original_weights(free_area_el_lines)

    areas_table_df = _build_workbook_like_table(
        values_by_line=areas_by_line,
        line_order=free_area_el_lines,
        original_weights=original_weights,
    )
    ratios_table_df = _build_workbook_like_table(
        values_by_line=ratios_by_line,
        line_order=free_area_el_lines,
        original_weights=original_weights,
    )

    ratios_table_df.insert(0, "Reference peak", ref_peak)

    summary_df = _build_summary_table(
        std_id=std_ID,
        areas_df=raw_df,
        free_lines=free_area_el_lines,
        abs_corr_threshold=abs_corr_threshold,
        sum_locked_neg_corr_threshold=sum_locked_neg_corr_threshold,
        sum_cv_threshold=sum_cv_threshold,
    )

    correlated_adjustments_df = _build_correlated_adjustments_table(
        summary_df=summary_df,
        ratios_by_line=ratios_by_line,
        areas_by_line=areas_by_line,
        original_weights=original_weights,
        ref_peak=ref_peak,
    )

    # Reduce exported numeric precision for easier manual use in calibration workbooks.
    raw_df_out = _round_dataframe_sig_figs(raw_df, n_sig_figs=n_sig_figs)
    areas_table_df_out = _round_dataframe_sig_figs(areas_table_df, n_sig_figs=n_sig_figs)
    ratios_table_df_out = _round_dataframe_sig_figs(ratios_table_df, n_sig_figs=n_sig_figs)
    correlated_adjustments_df_out = _round_dataframe_sig_figs(
        correlated_adjustments_df, n_sig_figs=n_sig_figs
    )
    summary_df_out = _round_dataframe_sig_figs(summary_df, n_sig_figs=n_sig_figs)

    # Save outputs in a dedicated per-element calibration folder.
    line_elements = sorted({el_line.split("_", 1)[0] for el_line in free_area_el_lines if "_" in el_line})
    if len(line_elements) == 1:
        el_tag = line_elements[0]
    elif len(line_elements) > 1:
        el_tag = "-".join(line_elements)
    else:
        el_tag = std_ID

    output_dir = sample_dir / f"{el_tag}_peaks_area_calib"
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H-%M")
    output_raw_path = output_dir / f"{std_ID}_peak_weights_raw_{timestamp}.csv"
    output_areas_table_path = output_dir / f"{std_ID}_peak_weights_areas_table_{timestamp}.csv"
    output_ratios_table_path = output_dir / f"{std_ID}_peak_weights_ratios_table_{timestamp}.csv"
    output_corr_adjust_path = output_dir / f"{std_ID}_peak_weights_correlated_adjustments_{timestamp}.csv"
    output_summary_path = output_dir / f"{std_ID}_peak_weights_summary_{timestamp}.csv"

    raw_df_out.to_csv(output_raw_path, index=False)
    areas_table_df_out.to_csv(output_areas_table_path, index=False)
    ratios_table_df_out.to_csv(output_ratios_table_path, index=False)
    correlated_adjustments_df_out.to_csv(output_corr_adjust_path, index=False)
    summary_df_out.to_csv(output_summary_path, index=False)

    print_double_separator()
    logging.info("Saved raw peak areas to: %s", output_raw_path)
    logging.info("Saved workbook-like areas table to: %s", output_areas_table_path)
    logging.info("Saved workbook-like ratio table to: %s", output_ratios_table_path)
    logging.info("Saved correlated-peak adjustments to: %s", output_corr_adjust_path)
    logging.info("Saved peak-weight summary to: %s", output_summary_path)

    return {
        "raw_path": str(output_raw_path),
        "areas_table_path": str(output_areas_table_path),
        "ratios_table_path": str(output_ratios_table_path),
        "correlated_adjustments_path": str(output_corr_adjust_path),
        "summary_path": str(output_summary_path),
        "reference_peak": ref_peak,
        "raw_df": raw_df_out,
        "areas_table_df": areas_table_df_out,
        "ratios_table_df": ratios_table_df_out,
        "correlated_adjustments_df": correlated_adjustments_df_out,
        "summary_df": summary_df_out,
    }
