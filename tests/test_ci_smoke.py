#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI smoke suite entrypoint.

Imports the offline workflow and unit tests so GitHub Actions can run the full
suite via::

    pytest -q tests/test_ci_smoke.py

Hardware-only manual scripts remain capitalized as ``Test_*.py`` and are not
collected by pytest.
"""

from __future__ import annotations

import os

# Force a non-interactive backend for CI runners.
os.environ.setdefault("MPLBACKEND", "Agg")

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID, WULFENITE_MINI_ID

# Workflow tests (converted from former manual Test_*.py scripts)
from test_compositional_analysis import test_compositional_analysis_clustering
from test_ctape_detection_navcam import test_ctape_detection_navcam
from test_drift_correction import test_drift_correction_estimates_nonzero_shift
from test_quantification_and_analysis import test_batch_quantification_two_spectra
from test_quantification_single_spectrum import test_single_spectrum_quantification
from test_segmentation_particles_in_frame import test_segmentation_particles_in_frame
from test_segmentation_single_particle import test_segmentation_single_particle
from test_xsp_measurement_spot_selection import test_xsp_measurement_spot_selection

# Existing unit tests
from test_dbscan_clustering import (  # noqa: F401
    test_clustering_config_method_is_normalized,
    test_clustering_config_rejects_unknown_method,
    test_dbscan_fingerprint_includes_and_reacts_to_params,
    test_dbscan_params_defaults,
    test_dbscan_params_rejects_invalid_eps,
    test_dbscan_params_rejects_invalid_min_samples,
    test_kmeans_fingerprint_excludes_dbscan_block,
    test_normalize_dbscan_labels_makes_contiguous_and_keeps_noise,
    test_run_dbscan_clustering_all_noise_returns_zero_clusters,
    test_run_dbscan_clustering_single_cluster_has_nan_silhouette,
    test_run_dbscan_clustering_two_clusters_with_noise,
)
from test_quantification_fingerprint import (  # noqa: F401
    test_extract_reference_values_filters_to_quant_ref_lines,
    test_fingerprint_changes_when_sample_element_reference_value_changes,
    test_fingerprint_ignores_unrelated_element_reference_value_changes,
    test_fingerprint_includes_sample_element_reference_values,
    test_fingerprint_scopes_reference_values_to_sample_and_substrate_elements,
    test_get_reference_values_loads_standards_when_standards_dict_is_none,
    test_get_reference_values_reuses_cache_when_live_refs_match_within_tolerance,
    test_reference_values_changed_respects_relevant_elements,
)
from test_msa_writer_regressions import (  # noqa: F401
    test_write_spectrum_pointer_file_preserves_template_tail,
)
from test_xsp_spot_selector import (  # noqa: F401
    test_validate_xsp_spot_pixels_accepts_in_bounds_points,
    test_validate_xsp_spot_pixels_deduplicates,
    test_validate_xsp_spot_pixels_rejects_out_of_bounds,
)


def test_ci_fixtures_present():
    wulf_dir = INPUTS_DIR / WULFENITE_MINI_ID
    k412_dir = INPUTS_DIR / K412_CLUSTER_MINI_ID

    assert (wulf_dir / "ledger.json").is_file()
    assert (wulf_dir / "spectra" / "spectrum_0.msa").is_file()
    assert (wulf_dir / "spectra" / "spectrum_1.msa").is_file()
    assert (k412_dir / "ledger.json").is_file()
