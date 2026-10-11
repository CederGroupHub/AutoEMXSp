#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility coverage for pymatgen Composition weight-fraction APIs."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from autoemx.config.runtime_configs import ExpStandardsConfig
from autoemx.config.schema_models.ledger import LedgerConfigs
from autoemx.utils.helper import composition_as_weight_dict
from pymatgen.core import Composition


def test_composition_as_weight_dict_uses_installed_pymatgen():
    weights = composition_as_weight_dict(Composition("PbMoO4"))
    assert set(weights) == {"Pb", "Mo", "O"}
    assert weights["Pb"] == pytest.approx(0.5644, abs=1e-3)
    assert sum(weights.values()) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize(
    "comp",
    [
        SimpleNamespace(as_weight_dict=lambda: {"Pb": 0.5, "O": 0.5}),
        SimpleNamespace(as_weight_dict={"Pb": 0.5, "O": 0.5}),
        SimpleNamespace(to_weight_dict=lambda: {"Pb": 0.5, "O": 0.5}),
        SimpleNamespace(to_weight_dict={"Pb": 0.5, "O": 0.5}),
    ],
)
def test_composition_as_weight_dict_accepts_legacy_and_current_shapes(comp):
    assert composition_as_weight_dict(comp) == {"Pb": 0.5, "O": 0.5}


def test_composition_as_weight_dict_prefers_as_weight_dict():
    comp = SimpleNamespace(
        as_weight_dict=lambda: {"Pb": 1.0},
        to_weight_dict={"O": 1.0},
    )
    assert composition_as_weight_dict(comp) == {"Pb": 1.0}


def test_exp_standards_config_computes_weight_fractions():
    cfg = ExpStandardsConfig(is_exp_std_measurement=True, formula="PbMoO4")
    w_frs = cfg.formula_w_frs()
    assert w_frs is not None
    assert w_frs["Pb"] == pytest.approx(0.5644, abs=1e-3)


def _std_ledger_configs_payload(exp_stds_cfg: dict) -> dict:
    return {
        "microscope_cfg": {"ID": "PhenomXL", "type": "SEM"},
        "sample_cfg": {"elements": ["Pb", "Mo", "O"], "type": "powder"},
        "measurement_cfg": {"exp_stds_cfg": exp_stds_cfg},
        "sample_substrate_cfg": {"elements": ["C"]},
        "plot_cfg": {},
    }


def test_legacy_exp_std_w_frs_moved_to_sample_cfg():
    # Old ledgers stored the standard's mass fractions in exp_stds_cfg (possibly overridden by the user)
    payload = _std_ledger_configs_payload(
        {"is_exp_std_measurement": True, "formula": "PbMoO4", "w_frs": {"Pb": 0.5, "Mo": 0.3, "O": 0.2}}
    )
    configs = LedgerConfigs.model_validate(payload)
    assert configs.sample_cfg.w_frs == {"Pb": 0.5, "Mo": 0.3, "O": 0.2}
    assert "w_frs" not in configs.measurement_cfg.exp_stds_cfg.model_dump()


def test_exp_std_w_frs_set_from_formula_when_missing():
    configs = LedgerConfigs.model_validate(
        _std_ledger_configs_payload({"is_exp_std_measurement": True, "formula": "PbMoO4"})
    )
    assert configs.sample_cfg.w_frs["Pb"] == pytest.approx(0.5644, abs=1e-3)
