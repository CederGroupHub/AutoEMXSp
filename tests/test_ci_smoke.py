#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI-safe smoke tests for core package workflows.

These tests intentionally avoid hardware-dependent paths and interactive plotting.
"""

import os
from pathlib import Path

from autoemx.runners.fit_and_quantify_spectrum_from_ledger import (
    fit_and_quantify_spectrum_from_ledger,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_PATH = REPO_ROOT / "examples" / "Results"


# Force a non-interactive backend for CI runners.
os.environ.setdefault("MPLBACKEND", "Agg")


def test_single_spectrum_quantification_smoke():
    quantifier = fit_and_quantify_spectrum_from_ledger(
        sample_ID="Wulfenite_example",
        spectrum_ID=4,
        results_path=str(RESULTS_PATH),
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
    assert hasattr(quantifier, "fit_result")
    assert quantifier.fit_result is not None
