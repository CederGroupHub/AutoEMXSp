#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paths and sample IDs for CI ledger fixtures under tests/inputs/."""

from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
INPUTS_DIR = TESTS_DIR / "inputs"
OUTPUTS_DIR = TESTS_DIR / "outputs"

WULFENITE_MINI_ID = "Wulfenite_mini"
K412_CLUSTER_MINI_ID = "K412_cluster_mini"
RESULTS_PATH = str(INPUTS_DIR)
