#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for quantification scientific-fingerprint / reference-value reuse."""

from types import SimpleNamespace
from typing import Any, Dict

import pytest

from autoemx.config.schema_models.quantification import QuantificationConfig
from autoemx.core.composition_analysis.analyser import EMXSp_Composition_Analyzer


def _minimal_options(**overrides: Any) -> Dict[str, Any]:
    options = {
        "method": "PB",
        "spectrum_lims": [0, 2048],
        "fit_tolerance": 1e-3,
        "use_instrument_background": False,
    }
    options.update(overrides)
    return options


def _quant_config(
    *,
    quantification_id: int = 0,
    sample_elements=None,
    substrate_elements=None,
    reference_values_by_el_line=None,
    reference_lines_by_element=None,
    options=None,
) -> QuantificationConfig:
    refs = dict(reference_values_by_el_line or {})
    lines = dict(reference_lines_by_element or {})
    if refs and not lines:
        lines = QuantificationConfig.derive_reference_lines_by_element(refs)
    return QuantificationConfig(
        quantification_id=quantification_id,
        sample_elements=list(sample_elements or ["Mo", "Pb"]),
        substrate_elements=list(substrate_elements or []),
        options=options or _minimal_options(),
        reference_values_by_el_line=refs,
        reference_lines_by_element=lines,
    )


def test_fingerprint_includes_sample_element_reference_values():
    payload = _quant_config(
        reference_values_by_el_line={"Mo_La1": 1.23, "Pb_Ma1": 4.56},
    ).fingerprint_payload()

    assert payload["reference_values_by_el_line"]["Mo_La1"] == pytest.approx(1.23)
    assert payload["reference_values_by_el_line"]["Pb_Ma1"] == pytest.approx(4.56)


def test_fingerprint_scopes_reference_values_to_sample_and_substrate_elements():
    payload = _quant_config(
        sample_elements=["Mo"],
        substrate_elements=["C"],
        reference_values_by_el_line={
            "Mo_La1": 1.0,
            "C_Ka1": 2.0,
            "Zn_Ka1": 9.0,  # unrelated element must not affect the signature
        },
    ).fingerprint_payload()

    assert set(payload["reference_values_by_el_line"]) == {"C_Ka1", "Mo_La1"}
    assert "Zn_Ka1" not in payload["reference_values_by_el_line"]


def test_fingerprint_changes_when_sample_element_reference_value_changes():
    base = _quant_config(reference_values_by_el_line={"Mo_La1": 1.0, "Pb_Ma1": 2.0})
    changed = _quant_config(reference_values_by_el_line={"Mo_La1": 1.5, "Pb_Ma1": 2.0})

    assert base.fingerprint() != changed.fingerprint()
    diffs = base.fingerprint_differences(changed)
    assert "reference_values_by_el_line.Mo_La1" in diffs


def test_fingerprint_ignores_unrelated_element_reference_value_changes():
    base = _quant_config(
        sample_elements=["Mo"],
        reference_values_by_el_line={"Mo_La1": 1.0, "Zn_Ka1": 3.0},
    )
    changed = _quant_config(
        sample_elements=["Mo"],
        reference_values_by_el_line={"Mo_La1": 1.0, "Zn_Ka1": 9.9},
    )

    assert base.fingerprint() == changed.fingerprint()


def test_reference_values_changed_respects_relevant_elements():
    current = {"Mo_La1": 1.001, "Zn_Ka1": 9.0}
    cached = {"Mo_La1": 1.004, "Zn_Ka1": 1.0}

    # Within 2-decimal tolerance for Mo; Zn ignored when not relevant.
    assert not EMXSp_Composition_Analyzer._reference_values_changed(
        current,
        cached,
        decimals=2,
        relevant_elements={"Mo"},
    )
    assert EMXSp_Composition_Analyzer._reference_values_changed(
        current,
        cached,
        decimals=2,
        relevant_elements={"Mo", "Zn"},
    )


def test_get_reference_values_loads_standards_when_standards_dict_is_none():
    """Default quant path leaves standards_dict=None; signature must still compare live refs."""
    live_refs = {"Mo_La1": 2.5, "Pb_Ma1": 3.5}
    active = _quant_config(reference_values_by_el_line={"Mo_La1": 1.0, "Pb_Ma1": 1.0})

    analyzer = SimpleNamespace(
        all_els_sample=["Mo", "Pb"],
        all_els_substrate=[],
        detectable_els_sample=["Mo", "Pb"],
        detectable_els_substrate=[],
        quant_cfg=SimpleNamespace(use_project_specific_std_dict=False, method="PB"),
        standards_dict=None,
        standards=None,
        measurement_cfg=SimpleNamespace(mode="TEM"),
        _extract_reference_values_from_standards=lambda load_if_missing: (
            live_refs if load_if_missing else {}
        ),
        _reference_values_changed=EMXSp_Composition_Analyzer._reference_values_changed,
    )

    with pytest.warns(UserWarning, match="Reference values in standards differ"):
        result = EMXSp_Composition_Analyzer._get_reference_values_by_el_line(
            analyzer,
            active_quant_config=active,
        )

    assert result == live_refs


def test_get_reference_values_reuses_cache_when_live_refs_match_within_tolerance():
    live_refs = {"Mo_La1": 1.004, "Pb_Ma1": 2.003}
    cached_refs = {"Mo_La1": 1.001, "Pb_Ma1": 2.001}
    active = _quant_config(reference_values_by_el_line=cached_refs)

    analyzer = SimpleNamespace(
        all_els_sample=["Mo", "Pb"],
        all_els_substrate=[],
        detectable_els_sample=["Mo", "Pb"],
        detectable_els_substrate=[],
        quant_cfg=SimpleNamespace(use_project_specific_std_dict=False, method="PB"),
        standards_dict=None,
        standards=None,
        measurement_cfg=SimpleNamespace(mode="TEM"),
        _extract_reference_values_from_standards=lambda load_if_missing: (
            live_refs if load_if_missing else {}
        ),
        _reference_values_changed=EMXSp_Composition_Analyzer._reference_values_changed,
    )

    result = EMXSp_Composition_Analyzer._get_reference_values_by_el_line(
        analyzer,
        active_quant_config=active,
    )
    assert result == cached_refs


def test_extract_reference_values_filters_to_quant_ref_lines():
    import autoemx.utils.constants as cnst

    standards_by_line = {
        "Mo_La1": [{cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY, cnst.COR_PB_DF_KEY: 1.2}],
        "Mo_Lb1": [{cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY, cnst.COR_PB_DF_KEY: 9.9}],  # not a quant ref line
        "Zn_Ka1": [{cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY, cnst.COR_PB_DF_KEY: 3.3}],  # not a sample element
    }

    analyzer = SimpleNamespace(
        quant_cfg=SimpleNamespace(use_project_specific_std_dict=False, method="PB"),
        standards_dict=standards_by_line,
        standards=None,
        measurement_cfg=SimpleNamespace(mode="TEM"),
        detectable_els_sample=["Mo"],
        detectable_els_substrate=[],
    )

    result = EMXSp_Composition_Analyzer._extract_reference_values_from_standards(
        analyzer,
        load_if_missing=False,
    )

    assert result == {"Mo_La1": 1.2}
