#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI tests for single-spectrum fitting and quantification."""

from autoemx.runners.fit_and_quantify_spectrum_from_ledger import (
    fit_and_quantify_spectrum_from_ledger,
)

from ci_paths import WULFENITE_MINI_ID


def test_single_spectrum_quantification(results_path):
    quantifier = fit_and_quantify_spectrum_from_ledger(
        sample_ID=WULFENITE_MINI_ID,
        spectrum_ID=0,
        results_path=results_path,
        is_standard=False,
        quantify_plot=False,
        plot_signal=False,
        zoom_plot=False,
        fit_tol=1e-4,
        is_particle=True,
        max_undetectable_w_fr=0,
        force_single_iteration=False,
        interrupt_fits_bad_spectra=False,
        print_results=False,
        quant_verbose=False,
        fitting_verbose=False,
    )

    assert quantifier is not None
    assert getattr(quantifier, "fit_result", None) is not None
