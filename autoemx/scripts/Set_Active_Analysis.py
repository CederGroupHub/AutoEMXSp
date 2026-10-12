#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Choose the active quantification run and clustering analysis of a sample.

All quantification runs and clustering analyses are kept in the sample ledger. The active ones provide the
settings reused by later quantifications and analyses, and are shown by default in the GUI.

This script is a thin wrapper around the runner:
autoemx.runners.set_active_analysis.set_active_analysis

Leave quantification_id as None to only list the analyses of the sample (the active ones are marked with *).

Use one of the following input modes:
1) sample_path: full path to the sample directory containing ledger.json
2) sample_ID + project_path: sample folder name and its parent directory
"""

# =============================================================================
# Sample Definition
# =============================================================================
from pathlib import Path
import autoemx.utils.constants as cnst

# Option A: Full path to sample folder (contains ledger.json)
sample_path = None

# Option B: Sample identifier + project path (parent directory of sample_ID folder)
sample_ID = "Wulfenite_example"
# Defaults to examples/Results in the repository, so execution cwd does not matter.
project_path = str(Path(__file__).resolve().parents[2] / "examples" / cnst.RESULTS_DIR)

# =============================================================================
# Analysis to make active
# =============================================================================
quantification_id = None  # Id of the quantification run to make active. None: only list the analyses
clustering_id = None  # Id of its clustering analysis to make active. None: the run keeps its active analysis

# =============================================================================
# Run
# =============================================================================
from autoemx.runners.set_active_analysis import list_analyses, set_active_analysis

if quantification_id is None:
    list_analyses(sample_path=sample_path, sample_ID=sample_ID, project_path=project_path)
else:
    set_active_analysis(
        quantification_id,
        clustering_id,
        sample_path=sample_path,
        sample_ID=sample_ID,
        project_path=project_path,
    )
