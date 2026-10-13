#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the configurable pre-fit total-counts threshold."""

from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

from autoemx.config.runtime_configs import QuantificationOptionsConfig
from autoemx.core.quantifier.quant_flags import detect_prefit_spectrum_issues
from autoemx.core.composition_analysis.ledger import LedgerModule
from autoemx.config.schema_models.quantification import QuantificationConfig


def _flat_spectrum(total_counts: float, n_channels: int = 100) -> np.ndarray:
    return np.full(n_channels, total_counts / n_channels, dtype=float)


def _prefit_kwargs(**overrides):
    kwargs = {
        "spectrum": _flat_spectrum(40000),
        "energy_vals": np.linspace(0.1, 10.0, 100),
        "sp_start": 0,
        "sp_end": 100,
        "target_acquisition_counts": 50000,
        "min_bckgrnd_cnts": None,
        "min_total_counts_fraction": 0.9,
    }
    kwargs.update(overrides)
    return kwargs


def test_default_fraction_flags_counts_below_90_percent():
    flag, comment = detect_prefit_spectrum_issues(**_prefit_kwargs())
    assert flag == 2
    assert comment == "Total counts too low"


def test_lower_fraction_allows_the_same_spectrum():
    flag, comment = detect_prefit_spectrum_issues(
        **_prefit_kwargs(min_total_counts_fraction=0.7)
    )
    assert flag != 2
    assert comment != "Total counts too low"


def test_zero_fraction_disables_the_total_counts_check():
    flag, _comment = detect_prefit_spectrum_issues(
        **_prefit_kwargs(
            spectrum=_flat_spectrum(1),
            min_total_counts_fraction=0.0,
        )
    )
    assert flag != 2


def test_empty_spectrum_still_flagged_as_missing():
    flag, comment = detect_prefit_spectrum_issues(
        **_prefit_kwargs(spectrum=np.array([]), min_total_counts_fraction=0.0)
    )
    assert flag == 1
    assert comment == "No spectral data present"


def test_quantification_options_default_is_90_percent():
    cfg = QuantificationOptionsConfig()
    assert cfg.min_total_counts_fraction == pytest.approx(0.9)


@pytest.mark.parametrize("bad_value", [-0.1, 1.1, float("nan"), float("inf")])
def test_quantification_options_rejects_invalid_fraction(bad_value):
    with pytest.raises(ValidationError):
        QuantificationOptionsConfig(min_total_counts_fraction=bad_value)


def _options_analyzer(min_total_counts_fraction: float) -> SimpleNamespace:
    return SimpleNamespace(
        quant_cfg=SimpleNamespace(
            method="PB",
            spectrum_lims=(14, 1100),
            fit_tolerance=1e-4,
            use_instrument_background=False,
            min_total_counts_fraction=min_total_counts_fraction,
        ),
        _apply_geom_factors=False,
        measurement_cfg=SimpleNamespace(beam_energy_keV=15.0, emergence_angle=35.0),
        det_ch_offset=0.0,
        det_ch_width=0.01,
        _min_total_counts_fraction=lambda: min_total_counts_fraction,
    )


def test_build_quantification_options_always_includes_fraction():
    options = LedgerModule._build_quantification_options(
        _options_analyzer(0.9)
    )
    assert options["min_total_counts_fraction"] == pytest.approx(0.9)


def test_build_quantification_options_includes_non_default_fraction():
    options = LedgerModule._build_quantification_options(
        _options_analyzer(0.5)
    )
    assert options["min_total_counts_fraction"] == pytest.approx(0.5)


def test_missing_fraction_in_old_ledger_options_defaults_to_90_percent():
    required_options = {
        "method": "PB",
        "spectrum_lims": [0, 2048],
        "fit_tolerance": 1e-3,
        "use_instrument_background": False,
    }
    legacy = QuantificationConfig(
        quantification_id=0,
        sample_elements=["Mo"],
        options=required_options,
    )
    explicit = QuantificationConfig(
        quantification_id=0,
        sample_elements=["Mo"],
        options={**required_options, "min_total_counts_fraction": 0.9},
    )
    assert legacy.options["min_total_counts_fraction"] == pytest.approx(0.9)
    assert legacy.fingerprint() == explicit.fingerprint()
    assert "min_total_counts_fraction" in legacy.fingerprint_payload()["options"]


def test_non_default_fraction_changes_quantification_fingerprint():
    required_options = {
        "method": "PB",
        "spectrum_lims": [0, 2048],
        "fit_tolerance": 1e-3,
        "use_instrument_background": False,
    }
    base = QuantificationConfig(
        quantification_id=0,
        sample_elements=["Mo"],
        options=required_options,
    )
    changed = QuantificationConfig(
        quantification_id=0,
        sample_elements=["Mo"],
        options={**required_options, "min_total_counts_fraction": 0.5},
    )
    assert base.fingerprint() != changed.fingerprint()
    assert "options.min_total_counts_fraction" in base.fingerprint_differences(changed)
