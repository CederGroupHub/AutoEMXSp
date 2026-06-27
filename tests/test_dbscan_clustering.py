#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for DBSCAN clustering support in composition analysis."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from autoemx.config.ledger_schemas import ClusteringConfig, DBSCANParams
from autoemx.core.composition_analysis.clustering import ClusteringModule


# --- Config model -----------------------------------------------------------
def test_dbscan_params_defaults():
    params = DBSCANParams()
    assert params.eps > 0
    assert params.min_samples >= 1
    assert params.metric == "euclidean"


@pytest.mark.parametrize("bad_eps", [0, -0.1, float("nan"), float("inf")])
def test_dbscan_params_rejects_invalid_eps(bad_eps):
    with pytest.raises(ValueError):
        DBSCANParams(eps=bad_eps)


def test_dbscan_params_rejects_invalid_min_samples():
    with pytest.raises(ValueError):
        DBSCANParams(min_samples=0)


def test_clustering_config_rejects_unknown_method():
    with pytest.raises(ValueError):
        ClusteringConfig(method="not_a_method")


def test_clustering_config_method_is_normalized():
    assert ClusteringConfig(method="DBSCAN").method == "dbscan"


# --- Fingerprint behavior ---------------------------------------------------
def test_kmeans_fingerprint_excludes_dbscan_block():
    """Existing k-means configs must keep their fingerprint stable (no dbscan key)."""
    payload = ClusteringConfig(method="kmeans").fingerprint_payload()
    assert "dbscan" not in payload


def test_dbscan_fingerprint_includes_and_reacts_to_params():
    base = ClusteringConfig(method="dbscan", dbscan=DBSCANParams(eps=0.05))
    changed = ClusteringConfig(method="dbscan", dbscan=DBSCANParams(eps=0.04))
    assert "dbscan" in base.fingerprint_payload()
    assert base.fingerprint() != changed.fingerprint()


# --- Label normalization ----------------------------------------------------
def test_normalize_dbscan_labels_makes_contiguous_and_keeps_noise():
    raw = np.array([-1, 3, 3, 7, -1, 7])
    out = ClusteringModule._normalize_dbscan_labels(raw)
    assert out.tolist() == [-1, 0, 0, 1, -1, 1]


# --- End-to-end clustering helper -------------------------------------------
def _fake_analyzer(eps=0.05, min_samples=3, metric="euclidean"):
    return SimpleNamespace(
        clustering_cfg=SimpleNamespace(
            dbscan=DBSCANParams(eps=eps, min_samples=min_samples, metric=metric)
        )
    )


def test_run_dbscan_clustering_two_clusters_with_noise():
    rng = np.random.default_rng(0)
    cluster_a = rng.normal([0.1, 0.9], 0.01, size=(20, 2))
    cluster_b = rng.normal([0.9, 0.1], 0.01, size=(20, 2))
    outlier = np.array([[0.5, 0.5]])
    X = pd.DataFrame(np.vstack([cluster_a, cluster_b, outlier]), columns=["El1", "El2"])

    labels, centroids, k, sil_score, wcss = ClusteringModule._run_dbscan_clustering(
        _fake_analyzer(eps=0.05, min_samples=3), X
    )

    assert k == 2
    assert centroids.shape == (2, 2)
    assert int((labels == -1).sum()) == 1  # the outlier is noise
    assert set(labels.tolist()) == {-1, 0, 1}
    assert 0.0 < sil_score <= 1.0
    assert wcss >= 0.0


def test_run_dbscan_clustering_single_cluster_has_nan_silhouette():
    rng = np.random.default_rng(1)
    X = pd.DataFrame(rng.normal([0.3, 0.7], 0.01, size=(15, 2)), columns=["El1", "El2"])

    labels, centroids, k, sil_score, wcss = ClusteringModule._run_dbscan_clustering(
        _fake_analyzer(eps=0.1, min_samples=3), X
    )

    assert k == 1
    assert centroids.shape == (1, 2)
    assert np.isnan(sil_score)  # silhouette undefined for a single cluster


def test_run_dbscan_clustering_all_noise_returns_zero_clusters():
    # Far-apart points with min_samples high => everything is noise.
    X = pd.DataFrame([[0.0, 1.0], [1.0, 0.0], [0.5, 0.5]], columns=["El1", "El2"])

    labels, centroids, k, sil_score, wcss = ClusteringModule._run_dbscan_clustering(
        _fake_analyzer(eps=0.01, min_samples=3), X
    )

    assert k == 0
    assert centroids.shape == (0, 2)
    assert bool(np.all(labels == -1))
    assert np.isnan(sil_score)
