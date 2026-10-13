#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EMXSp_Composition_Analyzer.from_ledger: the analyzer of an acquired sample, from its ledger."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from ci_paths import INPUTS_DIR, WULFENITE_MINI_ID

from autoemx.config.ledger_io import ledger_runtime_configs, load_sample_ledger
from autoemx.config.runtime_configs import ExpStandardsConfig
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer

CONFIG_ATTRS = ("microscope_cfg", "sample_cfg", "measurement_cfg", "sample_substrate_cfg", "quant_cfg",
                "clustering_cfg", "plot_cfg", "powder_meas_cfg", "bulk_meas_cfg", "exp_stds_cfg")


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    return Path(shutil.copytree(INPUTS_DIR / WULFENITE_MINI_ID, tmp_path / WULFENITE_MINI_ID))


def _state(analyzer: EMXSp_Composition_Analyzer) -> dict:
    state = {attr: getattr(analyzer, attr).model_dump() for attr in CONFIG_ATTRS}
    state.update(sample_id=analyzer.sample_id, sample_result_dir=str(Path(analyzer.sample_result_dir).resolve()),
                 is_acquisition=analyzer.is_acquisition, elements=list(analyzer.all_els_sample))
    return state


def test_same_analyzer_as_explicit_constructor(sample_dir: Path):
    ledger = load_sample_ledger(sample_dir / "ledger.json")
    configs = ledger_runtime_configs(ledger)
    explicit = EMXSp_Composition_Analyzer(
        microscope_cfg=configs.microscope_cfg, sample_cfg=configs.sample_cfg, measurement_cfg=configs.measurement_cfg,
        sample_substrate_cfg=configs.sample_substrate_cfg, quant_cfg=configs.quant_cfg,
        initial_clustering_cfg=configs.clustering_cfg, powder_meas_cfg=configs.powder_meas_cfg,
        bulk_meas_cfg=configs.bulk_meas_cfg, plot_cfg=configs.plot_cfg, is_acquisition=False,
        sample_id=ledger.sample_id, results_dir=str(sample_dir), verbose=False,
    )
    from_ledger = EMXSp_Composition_Analyzer.from_ledger(load_sample_ledger(sample_dir / "ledger.json"), verbose=False)
    assert _state(from_ledger) == _state(explicit)


def test_ledger_given_as_folder_file_or_model(sample_dir: Path):
    states = [
        _state(EMXSp_Composition_Analyzer.from_ledger(source, verbose=False))
        for source in (sample_dir, sample_dir / "ledger.json", str(sample_dir), load_sample_ledger(sample_dir / "ledger.json"))
    ]
    assert all(state == states[0] for state in states)


def test_configs_changed_by_the_caller_are_used(sample_dir: Path):
    ledger = load_sample_ledger(sample_dir / "ledger.json")
    configs = ledger_runtime_configs(ledger)
    configs.quant_cfg.fit_tolerance = 1e-3
    analyzer = EMXSp_Composition_Analyzer.from_ledger(ledger, configs, sample_id="Renamed", verbose=False)
    assert analyzer.quant_cfg.fit_tolerance == 1e-3
    assert analyzer.sample_id == "Renamed"
    assert analyzer.is_acquisition is False


def test_quantified_as_ordinary_sample_unless_standard_settings_given(sample_dir: Path):
    analyzer = EMXSp_Composition_Analyzer.from_ledger(sample_dir, verbose=False)
    assert analyzer.exp_stds_cfg == ExpStandardsConfig()

    standard = ExpStandardsConfig(is_exp_std_measurement=True, formula="PbMoO4")
    analyzer = EMXSp_Composition_Analyzer.from_ledger(sample_dir, exp_stds_cfg=standard, verbose=False)
    assert analyzer.exp_stds_cfg.is_exp_std_measurement and analyzer.exp_stds_cfg.formula == "PbMoO4"
