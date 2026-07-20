#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI tests for batch quantification on the mini Wulfenite ledger."""

from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze

from ci_paths import WULFENITE_MINI_ID


def test_batch_quantification_two_spectra(results_path):
    """Quantify both spectra in the 2-spectrum Wulfenite mini ledger."""
    analyzers = batch_quantify_and_analyze(
        sample_IDs=[WULFENITE_MINI_ID],
        quantification_method="PB",
        min_bckgrnd_cnts=5,
        results_path=results_path,
        max_analytical_error=5,
        run_analysis=False,
        interrupt_fits_bad_spectra=False,
        force_requantification=True,
    )

    assert analyzers
    analyzer = analyzers[0]
    assert analyzer is not None
    n_spectra = len(getattr(analyzer, "spectra_quant_records", []) or [])
    assert n_spectra >= 2
