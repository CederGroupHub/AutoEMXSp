#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The ledger and its quantification runs record the AutoEMX version; newer ledgers name the version to install."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import autoemx
from autoemx._version import LEGACY_VERSION
from autoemx.config.ledger_schemas import LedgerVersionError, SampleLedger
from autoemx.core.composition_analysis.analyser import EMXSp_Composition_Analyzer
from autoemx.core.composition_analysis.ledger import LedgerModule
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID


@pytest.fixture
def ledger_path(tmp_path: Path) -> Path:
    """Copy of a ledger saved before the version field existed."""
    sample_dir = tmp_path / K412_CLUSTER_MINI_ID
    shutil.copytree(INPUTS_DIR / K412_CLUSTER_MINI_ID, sample_dir)
    return sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)


def _write_as_newer_version(ledger_path: Path, with_unknown_field: bool) -> None:
    payload = json.loads(ledger_path.read_text())
    payload["autoemx_version"] = "99.0.0"
    if with_unknown_field:
        payload["field_added_in_99"] = True
    ledger_path.write_text(json.dumps(payload))


def test_ledger_without_version_gets_current_version(ledger_path: Path):
    assert "autoemx_version" not in json.loads(ledger_path.read_text())
    assert SampleLedger.from_json_file(ledger_path).autoemx_version == autoemx.__version__


def test_old_quantification_runs_get_legacy_version(ledger_path: Path):
    ledger = SampleLedger.from_json_file(ledger_path)
    assert ledger.quantifications
    assert all(q.autoemx_version == LEGACY_VERSION == "<=1.4.0" for q in ledger.quantifications)


def test_version_not_in_quantification_fingerprint(ledger_path: Path):
    # A new AutoEMX version alone must not trigger a new quantification run
    quant_config = SampleLedger.from_json_file(ledger_path).quantifications[0]
    same_science = quant_config.model_copy(update={"autoemx_version": "99.0.0"})
    assert same_science.fingerprint() == quant_config.fingerprint()


def test_saving_records_installed_version(ledger_path: Path):
    SampleLedger.from_json_file(ledger_path).to_json_file(ledger_path)
    assert json.loads(ledger_path.read_text())["autoemx_version"] == autoemx.__version__


def test_newer_unreadable_ledger_names_version_to_install(ledger_path: Path):
    _write_as_newer_version(ledger_path, with_unknown_field=True)
    with pytest.raises(LedgerVersionError, match=r"AutoEMX 99\.0\.0.*pip install autoemx==99\.0\.0"):
        SampleLedger.from_json_file(ledger_path)


def test_newer_readable_ledger_loads(ledger_path: Path):
    _write_as_newer_version(ledger_path, with_unknown_field=False)
    assert SampleLedger.from_json_file(ledger_path).autoemx_version == "99.0.0"


def test_analyser_does_not_replace_newer_ledger(ledger_path: Path):
    # A ledger that cannot be read because it is newer must not be treated as missing (it would be overwritten)
    _write_as_newer_version(ledger_path, with_unknown_field=True)
    analyzer = EMXSp_Composition_Analyzer.__new__(EMXSp_Composition_Analyzer)
    analyzer.sample_result_dir = str(ledger_path.parent)
    with pytest.raises(LedgerVersionError):
        LedgerModule._load_existing_ledger(analyzer)
