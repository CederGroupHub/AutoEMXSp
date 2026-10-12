#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""is_known_phase_mixture is a quantification setting; violin plots of the intermixing are a plot setting."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from autoemx.config.ledger_schemas import SampleLedger
from autoemx.config.schema_models.quantification import QuantificationConfig
from autoemx.core.composition_analysis.ledger import LedgerModule
from autoemx.core.composition_analysis.reference_matching import (
    _draw_intermixing_violin_plots,
    _is_known_phase_mixture,
)
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID

OPTIONS = {"method": "PB", "spectrum_lims": [14, 1100], "fit_tolerance": 1e-4, "use_instrument_background": False}


@pytest.mark.parametrize("legacy_value", [True, False])
def test_legacy_powder_setting_moved_to_quantification_runs(tmp_path: Path, legacy_value: bool):
    ledger_path = tmp_path / "ledger.json"
    shutil.copy(INPUTS_DIR / K412_CLUSTER_MINI_ID / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT), ledger_path)
    payload = json.loads(ledger_path.read_text())
    payload["configs"]["measurement_cfg"]["powder_meas_cfg"]["is_known_powder_mixture_meas"] = legacy_value
    ledger_path.write_text(json.dumps(payload))

    ledger = SampleLedger.from_json_file(ledger_path)
    assert [q.options["is_known_phase_mixture"] for q in ledger.quantifications] == [legacy_value]
    assert "is_known_powder_mixture_meas" not in ledger.configs.measurement_cfg.powder_meas_cfg.model_dump()


def test_known_phase_mixture_not_in_fingerprint():
    # The reference P/B values it selects are compared instead (reference_values_by_el_line)
    known = QuantificationConfig(quantification_id=0, options={**OPTIONS, "is_known_phase_mixture": True})
    regular = QuantificationConfig(quantification_id=0, options={**OPTIONS, "is_known_phase_mixture": False})
    assert known.fingerprint() == regular.fingerprint()


def test_reused_run_stores_current_value():
    reused = QuantificationConfig(quantification_id=0, options={**OPTIONS, "is_known_phase_mixture": False})
    candidate = QuantificationConfig(quantification_id=0, options={**OPTIONS, "is_known_phase_mixture": True})
    LedgerModule._update_non_fingerprint_options(reused, candidate)
    assert reused.options["is_known_phase_mixture"] is True


@pytest.mark.parametrize(
    "run_value, violin_setting, expected",
    [(True, None, True), (False, None, False), (False, True, True), (True, False, False)],
)
def test_violin_plots_follow_active_run_unless_set(run_value, violin_setting, expected):
    analyzer = SimpleNamespace(
        current_quant_config=SimpleNamespace(options={"is_known_phase_mixture": run_value}),
        quant_cfg=SimpleNamespace(is_known_phase_mixture=not run_value),  # the active run takes precedence
        plot_cfg=SimpleNamespace(intermixing_violin_plots=violin_setting),
    )
    assert _is_known_phase_mixture(analyzer) is run_value
    assert _draw_intermixing_violin_plots(analyzer) is expected
