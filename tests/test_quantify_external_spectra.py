#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Import and quantification of spectra acquired outside AutoEMX (2 spectra of the Wulfenite test sample)."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from ci_paths import INPUTS_DIR, WULFENITE_MINI_ID

from autoemx.runners.quantify_external_spectra import quantify_external_spectra


def _external_spectra(tmp_path: Path) -> Path:
    return Path(shutil.copytree(INPUTS_DIR / WULFENITE_MINI_ID / "spectra", tmp_path / "external"))


def _ledger(results: Path) -> dict:
    return json.loads((results / "Ext" / "ledger.json").read_text())


def test_ledger_built_without_quantification(tmp_path: Path):
    results = tmp_path / "results"
    analyzers = quantify_external_spectra(
        [{"ID": "Ext", "els": ["Pb", "Mo", "O"], "spectra_dir": str(_external_spectra(tmp_path))}],
        results_path=str(results), quantify=False, verbose=False,
    )
    assert analyzers == []
    ledger = _ledger(results)
    assert len(ledger["spectra"]) == 2
    assert ledger["quantifications"] == []
    assert ledger["configs"]["sample_cfg"]["elements"] == ["Pb", "Mo", "O"]
    assert ledger["configs"]["measurement_cfg"]["emergence_angle"] == 28.5  # from the detector model


def test_analysis_uses_the_settings_given(tmp_path: Path):
    results = tmp_path / "results"
    analyzers = quantify_external_spectra(
        [{"ID": "Ext", "els": ["Pb", "Mo", "O"], "spectra_dir": str(_external_spectra(tmp_path)), "cnd": ["PbMoO4"]}],
        results_path=str(results), max_n_clusters=3, quant_flags_accepted=[0, 1], num_CPU_cores=1, verbose=False,
    )
    assert len(analyzers) == 1
    [quantification] = _ledger(results)["quantifications"]
    active = quantification["clustering_analyses"][quantification["active_clustering_analysis_index"]]
    config = active.get("config", active)
    assert config["ref_formulae"] == ["PbMoO4"]
    assert config["max_k"] == 3
    assert config["quant_flags_accepted"] == [0, 1]
