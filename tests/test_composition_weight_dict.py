#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compatibility coverage for pymatgen Composition weight-fraction APIs."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from autoemx.config.runtime_configs import ExpStandardsConfig
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
    assert cfg.w_frs is not None
    assert cfg.w_frs["Pb"] == pytest.approx(0.5644, abs=1e-3)
