#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Choosing the active clustering analysis: analyze_sample(activate_new_analysis=...) and set_active_analysis."""

from __future__ import annotations

import builtins
import json
import shutil
from pathlib import Path

import pytest

from autoemx.runners.analyze_sample import analyze_sample
from autoemx.runners.set_active_analysis import format_analyses, set_active_analysis
from autoemx.config.ledger_io import load_sample_ledger
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    # K-412 test ledger: one quantification run (0) with one clustering analysis (0), active
    sample_dir = tmp_path / K412_CLUSTER_MINI_ID
    shutil.copytree(INPUTS_DIR / K412_CLUSTER_MINI_ID, sample_dir)
    return sample_dir


def _active(sample_dir: Path) -> tuple:
    """(active quantification id, active clustering id, number of clustering analyses)."""
    ledger = json.loads((sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)).read_text())
    [quant_config] = [q for q in ledger["quantifications"] if q["quantification_id"] == ledger["active_quant"]]
    analyses = quant_config["clustering_analyses"]
    return ledger["active_quant"], analyses[quant_config["active_clustering_analysis_index"]]["config"]["clustering_id"], len(analyses)


def _analyse_with_new_settings(sample_dir: Path, **kwargs) -> None:
    # A different threshold creates a new clustering analysis (clustering 1)
    analyze_sample(K412_CLUSTER_MINI_ID, results_path=str(sample_dir.parent), show_plots=False,
                   min_bckgrnd_cnts=4.5, **kwargs)


def test_new_analysis_active_by_default(sample_dir: Path):
    _analyse_with_new_settings(sample_dir)
    assert _active(sample_dir) == (0, 1, 2)


def test_new_analysis_saved_but_not_active(sample_dir: Path):
    _analyse_with_new_settings(sample_dir, activate_new_analysis=False)
    assert _active(sample_dir) == (0, 0, 2)
    assert (sample_dir / "analysis_quant0_clust1").exists()  # results of the new analysis are kept


@pytest.mark.parametrize("answer, expected_active", [("n", 0), ("no", 0), ("y", 1), ("", 1)])
def test_ask_before_activating(sample_dir: Path, monkeypatch, answer, expected_active):
    questions = []
    monkeypatch.setattr(builtins, "input", lambda question: questions.append(question) or answer)
    _analyse_with_new_settings(sample_dir, activate_new_analysis="ask")
    assert len(questions) == 1 and "clustering analysis 1" in questions[0]
    assert _active(sample_dir) == (0, expected_active, 2)


def test_ask_without_console_activates(sample_dir: Path, monkeypatch):
    def no_console(_question):
        raise EOFError
    monkeypatch.setattr(builtins, "input", no_console)
    _analyse_with_new_settings(sample_dir, activate_new_analysis="ask")
    assert _active(sample_dir) == (0, 1, 2)


def test_no_question_when_settings_unchanged(sample_dir: Path, monkeypatch):
    monkeypatch.setattr(builtins, "input", lambda question: pytest.fail("no new analysis: nothing to ask"))
    analyze_sample(K412_CLUSTER_MINI_ID, results_path=str(sample_dir.parent), show_plots=False,
                   activate_new_analysis="ask")
    assert _active(sample_dir) == (0, 0, 1)


def test_set_active_analysis(sample_dir: Path):
    _analyse_with_new_settings(sample_dir)
    set_active_analysis(0, 0, sample_path=str(sample_dir))
    assert _active(sample_dir) == (0, 0, 2)
    table = format_analyses(load_sample_ledger(sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)))
    assert "*   Clustering 0" in table and "    Clustering 1" in table


@pytest.mark.parametrize("quantification_id, clustering_id", [(5, None), (0, 7)])
def test_set_active_analysis_unknown_ids(sample_dir: Path, quantification_id, clustering_id):
    with pytest.raises(ValueError, match="not found"):
        set_active_analysis(quantification_id, clustering_id, sample_path=str(sample_dir))


def test_gui_set_active_analysis(sample_dir: Path):
    from autoemx.gui import backend as be
    from autoemx.gui.app import create_app

    _analyse_with_new_settings(sample_dir)
    info = be.load_sample_info(str(sample_dir))
    first = next(a for a in info.analyses if a.clust_id == 0)
    assert be.set_active_analysis(info, first.key) == "quant0 / clust0 is now the active analysis."
    assert _active(sample_dir) == (0, 0, 2)
    assert "already the active analysis" in be.set_active_analysis(be.load_sample_info(str(sample_dir)), first.key)
    # The app (with the "Set as active" callback) still builds
    create_app(str(sample_dir.parent))
