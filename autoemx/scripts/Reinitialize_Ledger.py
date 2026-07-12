#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reinitialize sample ledger state for quantification and clustering.

This script is a thin wrapper around the runner:
autoemx.runners.reinitialize_ledger.reinitialize_ledger

Use one of the following input modes:
1) sample_path: full path to the sample directory containing ledger.json
2) sample_ID + project_path: sample folder name and its parent directory

WARNING: This operation cannot be undone.
"""

# =============================================================================
# Sample Definition
# =============================================================================
import os
from pathlib import Path
import autoemx.utils.constants as cnst

# Option A: Full path to sample folder (contains ledger.json)
sample_path = None

# Option B: Sample identifier + project path (parent directory of sample_ID folder)
sample_ID = "Wulfenite_example"
# Defaults to examples/Results in the repository, so execution cwd does not matter.
project_path = str(Path(__file__).resolve().parents[2] / "examples" / cnst.RESULTS_DIR)

# =============================================================================
# Reset options
# =============================================================================
# Targeted removal: set to a quantification config id (int) or a list of ids to
# remove only those quantifications (their results, saved configs, and matching
# analysis_quant<id>_clust* folders). The active config is reverted to the latest
# remaining one. Leave as None to clear ALL quantification/clustering history.
quant_config_ids = None  # e.g. 2  or  [0, 3]

dry_run = False  # If True, prints the changes that would be made without modifying any files.
# Always recommended to run with dry_run=True first to verify the intended changes before applying them.

force = False  # If False, script asks for explicit YES confirmation.

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.reinitialize_ledger import reinitialize_ledger

summary = reinitialize_ledger(
    sample_path=sample_path,
    sample_ID=sample_ID,
    project_path=project_path,
    quant_config_ids=quant_config_ids,
    dry_run=dry_run,
    force=force,
    verbose=True,
)
