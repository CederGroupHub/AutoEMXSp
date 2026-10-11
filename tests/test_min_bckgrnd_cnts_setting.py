#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""min_bckgrnd_cnts is a quantification setting; ledgers that stored it in the clustering config are migrated."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from autoemx.config.schema_models.quantification import QuantificationConfig
from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze
import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID

OPTIONS = {"method": "PB", "spectrum_lims": [14, 1100], "fit_tolerance": 1e-4, "use_instrument_background": False}


def _legacy_quant_config(clustering_config: dict) -> dict:
    return {
        "quantification_id": 0,
        "options": dict(OPTIONS),
        "clustering_analyses": [{"config": clustering_config, "result": None}],
        "active_clustering_analysis_index": 0,
    }


def test_legacy_value_moved_from_first_clustering_analysis():
    config = QuantificationConfig.model_validate(_legacy_quant_config({"clustering_id": 0, "min_bckgrnd_cnts": 3.0}))
    assert config.options["min_bckgrnd_cnts"] == 3.0
    assert "min_bckgrnd_cnts" not in config.clustering_analyses[0].config.model_dump()


def test_legacy_run_without_value_gets_default():
    config = QuantificationConfig.model_validate(_legacy_quant_config({"clustering_id": 0}))
    assert config.options["min_bckgrnd_cnts"] == dflt.min_bckgrnd_cnts


def test_value_is_part_of_fingerprint():
    base = QuantificationConfig(quantification_id=0, options={**OPTIONS, "min_bckgrnd_cnts": 5.0})
    changed = QuantificationConfig(quantification_id=0, options={**OPTIONS, "min_bckgrnd_cnts": 3.0})
    assert base.fingerprint() != changed.fingerprint()


@pytest.mark.parametrize("drop_stored_value", [False, True])
def test_legacy_ledger_with_default_value_is_not_requantified(tmp_path: Path, drop_stored_value: bool):
    # The test ledger stores min_bckgrnd_cnts = 5.0 (the default) in its clustering config; older ledgers lack it
    sample_dir = tmp_path / K412_CLUSTER_MINI_ID
    shutil.copytree(INPUTS_DIR / K412_CLUSTER_MINI_ID, sample_dir)
    ledger_path = sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
    payload = json.loads(ledger_path.read_text())
    for quant_config in payload["quantifications"]:
        for analysis in quant_config["clustering_analyses"]:
            if drop_stored_value:
                analysis["config"].pop("min_bckgrnd_cnts", None)
            else:
                assert analysis["config"]["min_bckgrnd_cnts"] == dflt.min_bckgrnd_cnts
    ledger_path.write_text(json.dumps(payload))
    run_ids_before = [q["quantification_id"] for q in payload["quantifications"]]

    batch_quantify_and_analyze(
        sample_IDs=[K412_CLUSTER_MINI_ID], results_path=str(tmp_path), run_analysis=False, max_spectra_to_quantify=0
    )

    saved = json.loads(ledger_path.read_text())
    assert [q["quantification_id"] for q in saved["quantifications"]] == run_ids_before
    assert all(q["options"]["min_bckgrnd_cnts"] == dflt.min_bckgrnd_cnts for q in saved["quantifications"])
