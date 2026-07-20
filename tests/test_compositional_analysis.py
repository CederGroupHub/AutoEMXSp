#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CI tests for clustering / compositional analysis on the K-412 mini ledger."""

from autoemx.runners.analyze_sample import analyze_sample

from ci_paths import K412_CLUSTER_MINI_ID


def test_compositional_analysis_clustering(results_path):
    analyzer = analyze_sample(
        sample_ID=K412_CLUSTER_MINI_ID,
        results_path=results_path,
        k_forced=2,
        clustering_method="kmeans",
        max_analytical_error_percent=5,
        quant_flags_accepted=[0, -1],
        show_plots=False,
        plot_custom_plots=False,
        show_unused_compositions_cluster_plot=False,
        output_filename_suffix="_ci",
    )

    assert analyzer is not None
    clustering_cfg = getattr(analyzer, "clustering_cfg", None)
    assert clustering_cfg is not None
    resolved_k = getattr(clustering_cfg, "k_resolved", None) or getattr(
        clustering_cfg, "k_forced", None
    )
    assert resolved_k == 2
