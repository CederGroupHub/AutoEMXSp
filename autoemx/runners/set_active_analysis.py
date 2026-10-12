#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Choose the active quantification run and clustering analysis of a sample.

Every quantification run and clustering analysis is kept in the sample ledger. The active ones provide the
settings reused by later quantifications and analyses, and are shown by default in the GUI. This runner lists
them and makes a previous one active again, e.g. after declining a new analysis or to go back to earlier settings.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import get_active_quantification_config, load_sample_ledger
from autoemx.config.ledger_schemas import SampleLedger
from autoemx.runners.reinitialize_ledger import _resolve_sample_dir

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def _ledger_path(sample_dir: Path) -> Path:
    return sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"


def format_analyses(ledger: SampleLedger) -> str:
    """Table of the quantification runs and clustering analyses of a ledger, marking the active ones."""
    active_quant_config = get_active_quantification_config(ledger)
    lines = [f"Analyses of sample '{ledger.sample_id}' (* = active):"]
    for quant_config in ledger.quantifications:
        is_active_quant = active_quant_config is not None and quant_config is active_quant_config
        lines.append(
            f"{'*' if is_active_quant else ' '} Quantification {quant_config.quantification_id}"
            f" (AutoEMX {quant_config.autoemx_version})"
        )
        for index, analysis in enumerate(quant_config.clustering_analyses):
            cfg = analysis.config
            is_active = is_active_quant and quant_config.active_clustering_analysis_index == index
            k = len(analysis.result.centroids) if analysis.result is not None else None
            details = [cfg.method, cfg.geometry, f"k={k}" if k is not None else "not run",
                       f"max an. error {cfg.max_analytical_error_percent} w%",
                       f"min bckgrnd cnts {cfg.min_bckgrnd_cnts}"]
            lines.append(f"  {'*' if is_active else ' '}   Clustering {cfg.clustering_id}: {', '.join(details)}")
    if not ledger.quantifications:
        lines.append("  (no quantification yet)")
    return "\n".join(lines)


def list_analyses(
    sample_path: Optional[str] = None,
    sample_ID: Optional[str] = None,
    project_path: Optional[str] = None,
) -> str:
    """Print and return the table of the quantification runs and clustering analyses of a sample."""
    sample_dir = _resolve_sample_dir(sample_path, sample_ID, project_path)
    table = format_analyses(load_sample_ledger(_ledger_path(sample_dir)))
    logging.info(table)
    return table


def set_active_analysis(
    quantification_id: int,
    clustering_id: Optional[int] = None,
    sample_path: Optional[str] = None,
    sample_ID: Optional[str] = None,
    project_path: Optional[str] = None,
) -> SampleLedger:
    """
    Make a quantification run, and optionally one of its clustering analyses, the active ones.

    Parameters
    ----------
    quantification_id : int
        Id of the quantification run to make active.
    clustering_id : int, optional
        Id of the clustering analysis of that run to make active. If None, the run keeps its own active analysis.
    sample_path : str, optional
        Full path to the sample folder containing ledger.json.
    sample_ID : str, optional
        Sample identifier. Used only when sample_path is not provided.
    project_path : str, optional
        Parent directory containing sample folders. Used with sample_ID.

    Returns
    -------
    SampleLedger
        The updated ledger, saved to ledger.json.

    Raises
    ------
    ValueError
        If the quantification run or clustering analysis does not exist (the available ones are listed).
    """
    sample_dir = _resolve_sample_dir(sample_path, sample_ID, project_path)
    ledger_path = _ledger_path(sample_dir)
    ledger = load_sample_ledger(ledger_path)

    quant_config = next(
        (q for q in ledger.quantifications if q.quantification_id == int(quantification_id)), None
    )
    if quant_config is None:
        raise ValueError(f"Quantification {quantification_id} not found.\n{format_analyses(ledger)}")

    if clustering_id is not None:
        index = next(
            (i for i, a in enumerate(quant_config.clustering_analyses) if a.config.clustering_id == int(clustering_id)),
            None,
        )
        if index is None:
            raise ValueError(
                f"Clustering {clustering_id} not found in quantification {quantification_id}.\n{format_analyses(ledger)}"
            )
        quant_config.active_clustering_analysis_index = index

    ledger.active_quant = quant_config.quantification_id
    ledger.to_json_file(ledger_path)
    logging.info(f"✅ Active analysis updated.\n{format_analyses(ledger)}")
    return ledger
