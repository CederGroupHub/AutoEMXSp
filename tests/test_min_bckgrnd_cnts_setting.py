#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""min_bckgrnd_cnts (quant_flag 8) is a clustering setting, applied from the background counts stored with each
quantification result: changing it needs no re-quantification."""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from autoemx.config import ClusteringConfig
from autoemx.config.schema_models.quantification import QuantificationConfig, QuantificationResult
from autoemx.core.composition_analysis.clustering import effective_quant_flag
from autoemx.runners.analyze_sample import analyze_sample
from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze
import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID

OPTIONS = {"method": "PB", "spectrum_lims": [14, 1100], "fit_tolerance": 1e-4, "use_instrument_background": False}


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    # K-412: background counts under the reference peaks 4.0 (spectrum 14, quant_flag -1), 5.1, 5.1, 8.1, 5.1, 6.3;
    # quantified (and analysed) with min_bckgrnd_cnts = 5
    sample_dir = tmp_path / K412_CLUSTER_MINI_ID
    shutil.copytree(INPUTS_DIR / K412_CLUSTER_MINI_ID, sample_dir)
    return sample_dir


def _ledger(sample_dir: Path) -> dict:
    return json.loads((sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)).read_text())


def _stored_results(sample_dir: Path) -> list:
    return [
        [QuantificationResult.model_validate(r).model_dump() for r in s["quantification_results"]]
        for s in _ledger(sample_dir)["spectra"]
    ]


def _compositions(sample_dir: Path, analysis_id: int) -> pd.DataFrame:
    return pd.read_csv(sample_dir / f"analysis_quant0_clust{analysis_id}" / "Compositions.csv", index_col="Spectrum #")


def test_value_saved_by_1_4_1_moved_to_clustering_configs():
    config = QuantificationConfig.model_validate({
        "quantification_id": 0,
        "options": {**OPTIONS, "min_bckgrnd_cnts": 3.0},
        "clustering_analyses": [{"config": {"clustering_id": 0}, "result": None}],
    })
    assert "min_bckgrnd_cnts" not in config.options
    assert config.clustering_analyses[0].config.min_bckgrnd_cnts == 3.0


def test_value_is_part_of_clustering_fingerprint():
    assert ClusteringConfig(min_bckgrnd_cnts=5.0).fingerprint() != ClusteringConfig(min_bckgrnd_cnts=3.0).fingerprint()


def _analyzer(threshold: float, quantification_threshold: float = 5.0):
    quantification_config = SimpleNamespace(min_bckgrnd_cnts=quantification_threshold)
    return SimpleNamespace(
        clustering_cfg=SimpleNamespace(min_bckgrnd_cnts=threshold),
        current_quant_config=SimpleNamespace(clustering_analyses=[SimpleNamespace(config=quantification_config)]),
    )


def _record(flag, background, converged=True):
    return SimpleNamespace(
        quant_flag=flag, comment="stored", composition_atomic_fractions={"Si": 1.0},
        diagnostics=SimpleNamespace(min_background_ref_lines=background, converged=converged),
    )


@pytest.mark.parametrize("record, threshold, expected_flag", [
    (_record(-1, 4.0, converged=False), 5.0, -1),  # threshold used at quantification: stored flag
    (_record(0, 5.1), 6.0, 8),                     # higher threshold
    (_record(-1, 4.0, converged=False), 6.0, 8),
    (_record(8, 4.0), 3.0, 0),                     # lower threshold
    (_record(8, 4.0, converged=False), 3.0, -1),
    (_record(8, 4.0), 4.5, 8),
    (_record(4, 1.0), 6.0, 4),                     # other flags take precedence
    (_record(0, None), 6.0, 0),                    # no stored background counts
])
def test_effective_flag(record, threshold, expected_flag):
    assert effective_quant_flag(record, _analyzer(threshold))[0] == expected_flag


def test_changed_value_applied_at_analysis_without_requantification(sample_dir: Path):
    stored_results = _stored_results(sample_dir)

    analyze_sample(K412_CLUSTER_MINI_ID, results_path=str(sample_dir.parent), show_plots=False)
    analyze_sample(K412_CLUSTER_MINI_ID, results_path=str(sample_dir.parent), show_plots=False, min_bckgrnd_cnts=4.5)

    # Same threshold as at quantification: stored flags; threshold 4.5: spectrum 14 (4.0 counts) gets quant_flag 8
    assert _compositions(sample_dir, 0)["Quant_flag"].to_dict() == {14: -1, 25: 0, 26: 0, 64: 0, 78: 0, 87: 0}
    assert _compositions(sample_dir, 1)["Quant_flag"].to_dict() == {14: 8, 25: 0, 26: 0, 64: 0, 78: 0, 87: 0}
    [quant_config] = _ledger(sample_dir)["quantifications"]
    assert [a["config"]["min_bckgrnd_cnts"] for a in quant_config["clustering_analyses"]] == [5.0, 4.5]
    assert _stored_results(sample_dir) == stored_results  # quantification results not modified


def test_warning_for_results_without_stored_background(sample_dir: Path, caplog):
    ledger_path = sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
    payload = json.loads(ledger_path.read_text())
    payload["spectra"][1]["quantification_results"][0]["diagnostics"]["min_background_ref_lines"] = None
    ledger_path.write_text(json.dumps(payload))

    with caplog.at_level(logging.WARNING):
        analyze_sample(K412_CLUSTER_MINI_ID, results_path=str(sample_dir.parent), show_plots=False,
                       min_bckgrnd_cnts=4.5)
    assert "You should requantify first to apply this parameter" in caplog.text
    assert "cannot be applied to 1 spectra" in caplog.text


def test_legacy_ledger_with_default_value_is_not_requantified(sample_dir: Path):
    run_ids_before = [q["quantification_id"] for q in _ledger(sample_dir)["quantifications"]]
    batch_quantify_and_analyze(
        sample_IDs=[K412_CLUSTER_MINI_ID], results_path=str(sample_dir.parent), run_analysis=False,
        max_spectra_to_quantify=0,
    )
    saved = _ledger(sample_dir)
    assert [q["quantification_id"] for q in saved["quantifications"]] == run_ids_before
    assert saved["quantifications"][0]["clustering_analyses"][0]["config"]["min_bckgrnd_cnts"] == dflt.min_bckgrnd_cnts
