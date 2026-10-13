#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quantification flags: skipping fits for pre-fit issues, merging pre-fit and post-fit flags, post-fit flags."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from autoemx.core.quantifier.quant_flags import (
    check_fit_quant_validity,
    merge_prefit_quant_flag,
    should_interrupt_prefit_issue,
)


@pytest.mark.parametrize("flag, interrupt, expected", [
    (None, True, False), (1, False, True), (2, True, True), (2, False, False), (3, True, True), (3, False, False),
])
def test_should_interrupt_prefit_issue(flag, interrupt, expected):
    assert should_interrupt_prefit_issue(flag, interrupt) is expected


def test_prefit_flag_applied_only_to_otherwise_good_results():
    assert merge_prefit_quant_flag(0, "", 2, "Total counts too low") == (2, "Total counts too low")
    assert merge_prefit_quant_flag(4, "poor fit", 2, "Total counts too low") == (4, "poor fit")
    assert merge_prefit_quant_flag(0, "", None, None) == (0, "")


def _quantifier(substrate_counts: float = 0.0, total_counts: float = 1e5):
    return SimpleNamespace(tot_sp_counts=total_counts, fitted_peaks_info={"C_Ka1": {"peak_intensity": substrate_counts}})


@pytest.mark.parametrize("bad_quant_flag, is_valid, expected", [
    (1, True, 4), (2, True, 5), (3, False, 6), (None, False, 9), (None, True, 0), (-1, True, -1),
])
def test_post_fit_flags(bad_quant_flag, is_valid, expected):
    flag, _comment = check_fit_quant_validity(
        is_quant_fit_valid=is_valid, bad_quant_flag=bad_quant_flag, quantifier=_quantifier(),
        min_bckgrnd_ref_lines=100.0, detectable_els_substrate=["C"], min_bckgrnd_cnts=5.0,
    )
    assert flag == expected


def test_substrate_and_background_flags():
    flag, comment = check_fit_quant_validity(
        is_quant_fit_valid=True, bad_quant_flag=None, quantifier=_quantifier(substrate_counts=2e4),
        min_bckgrnd_ref_lines=100.0, detectable_els_substrate=["C"], min_bckgrnd_cnts=5.0,
    )
    assert flag == 7 and comment.startswith("C ")
    flag, _comment = check_fit_quant_validity(
        is_quant_fit_valid=True, bad_quant_flag=None, quantifier=_quantifier(),
        min_bckgrnd_ref_lines=2.0, detectable_els_substrate=["C"], min_bckgrnd_cnts=5.0,
    )
    assert flag == 8
