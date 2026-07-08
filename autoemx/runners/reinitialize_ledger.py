#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reset quantification and clustering state for a sample ledger.

This runner removes performed quantifications (including their saved
quantification/clustering configs) and associated analysis outputs from the
ledger, so processing can start from scratch while preserving acquisition
data.

WARNING: This operation is destructive and cannot be undone.
"""

from __future__ import annotations

import shutil
from glob import glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import autoemx.utils.constants as cnst
from autoemx.config.ledger_schemas import SampleLedger


def _resolve_sample_dir(
    sample_path: Optional[str],
    sample_ID: Optional[str],
    project_path: Optional[str],
) -> Path:
    """Resolve sample directory from either full sample path or ID + project path."""
    if sample_path:
        resolved = Path(sample_path).expanduser().resolve()
    else:
        if not sample_ID or not project_path:
            raise ValueError(
                "Provide either sample_path, or both sample_ID and project_path."
            )
        resolved = (Path(project_path).expanduser().resolve() / sample_ID).resolve()

    if not resolved.exists() or not resolved.is_dir():
        msg = [f"Sample directory not found: {resolved}"]
        if sample_path is None:
            msg.append(
                "Check project_path. For repository examples, use project_path=<repo>/examples/Results."
            )
        raise FileNotFoundError(" ".join(msg))

    ledger_path = resolved / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
    if not ledger_path.exists():
        raise FileNotFoundError(
            f"Ledger file not found in sample directory: {ledger_path}"
        )

    return resolved


def _collect_analysis_targets(sample_dir: Path) -> List[Path]:
    """Collect analysis directories produced by quantification/clustering runs."""
    targets: List[Path] = []

    for path in glob(str(sample_dir / "analysis_quant*_clust*")):
        targets.append(Path(path))

    for dirname in (cnst.ANALYSIS_DIR, cnst.ANALYSIS_SUBDIR):
        candidate = sample_dir / dirname
        if candidate.exists():
            targets.append(candidate)

    return targets


def _collect_legacy_root_csv_targets(sample_dir: Path) -> List[Path]:
    """Collect legacy CSV result files from the sample root."""
    targets: List[Path] = []
    for stem in (cnst.COMPOSITIONS_FILENAME, cnst.CLUSTERS_FILENAME):
        for path in glob(str(sample_dir / f"{stem}*.csv")):
            targets.append(Path(path))
    return targets


def _reset_ledger_quant_history(ledger_path: Path, dry_run: bool) -> Tuple[bool, int, int]:
    """Clear quantification configs and per-spectrum quantification results in ledger."""
    ledger = SampleLedger.from_json_file(ledger_path)

    total_results_removed = 0
    for spectrum in ledger.spectra:
        total_results_removed += len(spectrum.quantification_results)
        spectrum.quantification_results = []

    quant_configs_removed = len(ledger.quantifications)
    had_quant_configs = quant_configs_removed > 0 or ledger.active_quant is not None
    ledger.quantifications = []
    ledger.active_quant = None

    if not dry_run:
        ledger.to_json_file(ledger_path)

    ledger_reset = had_quant_configs or total_results_removed > 0
    return ledger_reset, quant_configs_removed, total_results_removed


def _delete_path(path: Path, dry_run: bool) -> None:
    """Delete a file or directory, unless dry-run is active."""
    if dry_run:
        return

    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def reinitialize_ledger(
    sample_path: Optional[str] = None,
    sample_ID: Optional[str] = None,
    project_path: Optional[str] = None,
    dry_run: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> Dict[str, int | bool | str]:
    """
    Remove quantification/clustering state and outputs for one sample.

    Parameters
    ----------
    sample_path : str, optional
        Full path to the sample folder containing ledger.json.
    sample_ID : str, optional
        Sample identifier. Used only when sample_path is not provided.
    project_path : str, optional
        Parent directory containing sample folders. Used with sample_ID.
    dry_run : bool, optional
        If True, only print planned actions.
    force : bool, optional
        If True, skip interactive confirmation.
    verbose : bool, optional
        If True, print detailed progress.

    Returns
    -------
    summary : dict
        Dictionary with cleanup results and counts, including
        ``quant_configs_removed`` (number of saved quantification configs
        cleared from the ledger) and ``quant_results_removed`` (number of
        per-spectrum quantification result records cleared).
    """
    sample_dir = _resolve_sample_dir(sample_path, sample_ID, project_path)
    ledger_path = sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"

    analysis_targets = _collect_analysis_targets(sample_dir)
    root_csv_targets = _collect_legacy_root_csv_targets(sample_dir)

    if verbose:
        print("=" * 72)
        print("WARNING: This action removes quantification/clustering history and results.")
        print("WARNING: It cannot be undone.")
        print(f"Sample directory: {sample_dir}")
        print("=" * 72)
        print("Planned actions:")
        print(f"- Remove saved quantification configs and reset quantification history inside ledger: {ledger_path}")
        for path in analysis_targets:
            print(f"- Remove analysis directory: {path}")
        for path in root_csv_targets:
            print(f"- Remove legacy result file: {path}")

    if dry_run and verbose:
        print("\nDRY RUN: no files will be modified.")

    if not force and not dry_run:
        answer = input("Type YES to continue: ").strip()
        if answer != "YES":
            if verbose:
                print("Aborted by user.")
            return {
                "sample_dir": str(sample_dir),
                "aborted": True,
                "ledger_reset": False,
                "quant_configs_removed": 0,
                "quant_results_removed": 0,
                "analysis_dirs_removed": 0,
                "legacy_files_removed": 0,
            }

    ledger_reset, quant_configs_removed, result_records_removed = _reset_ledger_quant_history(
        ledger_path=ledger_path,
        dry_run=dry_run,
    )

    deleted_dirs = 0
    deleted_files = 0
    for target in analysis_targets + root_csv_targets:
        if target.exists():
            _delete_path(target, dry_run)
            if target.is_dir():
                deleted_dirs += 1
            else:
                deleted_files += 1

    summary: Dict[str, int | bool | str] = {
        "sample_dir": str(sample_dir),
        "aborted": False,
        "ledger_reset": ledger_reset,
        "quant_configs_removed": quant_configs_removed,
        "quant_results_removed": result_records_removed,
        "analysis_dirs_removed": deleted_dirs,
        "legacy_files_removed": deleted_files,
    }

    if verbose:
        print("\nSummary:")
        print(f"- Ledger reset applied: {ledger_reset}")
        print(f"- Quantification configs removed from ledger: {quant_configs_removed}")
        print(f"- Quantification result records removed from ledger: {result_records_removed}")
        print(f"- Analysis directories removed: {deleted_dirs}")
        print(f"- Legacy result CSV files removed: {deleted_files}")
        if dry_run:
            print("\nDry run complete. No changes were written.")
        else:
            print("\nReinitialization complete. The next quantification will start from scratch.")

    return summary
