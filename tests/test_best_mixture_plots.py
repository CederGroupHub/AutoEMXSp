#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the plots of the best mixture of each cluster."""
import os
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import pytest

import autoemx.utils.constants as cnst
from autoemx.config.ledger_schemas import ClusteringConfig
from autoemx.config.runtime_configs import PlotConfig
from autoemx.core.composition_analysis.plotting import PlottingModule
from autoemx.core.composition_analysis.reference_matching import ReferenceMatchingModule


def _analyzer(tmp_path, els, refs, weights, els_to_plot=()):
    names = list(refs)
    return SimpleNamespace(
        sample_id='test',
        analysis_dir=str(tmp_path),
        ref_formulae=names,
        ref_phases_df=pd.DataFrame([refs[n] for n in names], columns=els),
        ref_weights_in_mixture=[weights[n] for n in names],
        detectable_els_sample=list(els),
        clustering_cfg=ClusteringConfig(),
        plot_cfg=PlotConfig(els_to_plot=list(els_to_plot)),
        quant_cfg=SimpleNamespace(is_known_phase_mixture=False),
    )


def _cluster(analyzer, phases, n=40, noise=0.005, seed=0):
    rng = np.random.default_rng(seed)
    H = analyzer.ref_phases_df.loc[[analyzer.ref_formulae.index(p) for p in phases]].to_numpy()
    X = np.clip(rng.dirichlet(np.ones(len(phases)), size=n) @ H + rng.normal(0, noise, (n, H.shape[1])), 0, None)
    return pd.DataFrame(X / X.sum(axis=1, keepdims=True), columns=analyzer.detectable_els_sample)


def _run(analyzer, compositions_df):
    labels = np.zeros(len(compositions_df), dtype=int)
    mixtures = ReferenceMatchingModule._assign_mixtures(
        analyzer, 1, labels, compositions_df, rms_dist_cluster=[1.0], max_raw_confs=[0.0],
        n_points_per_cluster=[len(compositions_df)],
    )
    PlottingModule._save_best_mixture_plots(analyzer, compositions_df, labels, mixtures)
    return sorted(os.listdir(analyzer.analysis_dir)), mixtures


CUALO = (['Cu', 'Al', 'O'], {'CuO': [0.5, 0, 0.5], 'Al2O3': [0, 0.4, 0.6], 'Cu': [1, 0, 0]},
         {'CuO': 2, 'Al2O3': 5, 'Cu': 1})


@pytest.mark.parametrize('phases', [['CuO', 'Al2O3'], ['CuO', 'Al2O3', 'Cu']])
def test_three_element_samples_get_ternary_plots(tmp_path, phases):
    analyzer = _analyzer(tmp_path, *CUALO)
    files, _ = _run(analyzer, _cluster(analyzer, phases))
    assert files == ['Mixture_plot_cl0.png', 'Mixture_plot_cl0_zoomed.png']


FOUR_ELS = (['Na', 'Al', 'Si', 'O'],
            {'NaAlSiO4': [1 / 7, 1 / 7, 1 / 7, 4 / 7], 'SiO2': [0, 0, 1 / 3, 2 / 3], 'Al2O3': [0, 0.4, 0, 0.6]},
            {'NaAlSiO4': 7, 'SiO2': 3, 'Al2O3': 5})


def test_binary_mixture_with_more_elements_gets_2d_plot(tmp_path):
    analyzer = _analyzer(tmp_path, *FOUR_ELS)
    files, _ = _run(analyzer, _cluster(analyzer, ['NaAlSiO4', 'SiO2']))
    assert files == ['Mixture_plot_cl0.png', 'Mixture_plot_cl0_zoomed.png']


def test_2d_plot_elements_differ_most_between_phases_or_follow_els_to_plot(tmp_path):
    els, refs, weights = FOUR_ELS
    H = np.array([refs['SiO2'], refs['Al2O3']])
    assert PlottingModule._mixture_plot_elements(_analyzer(tmp_path, *FOUR_ELS), H, els) == ['Al', 'Si']
    forced = _analyzer(tmp_path, *FOUR_ELS, els_to_plot=['Na', 'O'])
    assert PlottingModule._mixture_plot_elements(forced, H, els) == ['Na', 'O']


def test_ternary_mixture_with_more_elements_gets_phase_ternary(tmp_path):
    analyzer = _analyzer(tmp_path, *FOUR_ELS)
    files, mixtures = _run(analyzer, _cluster(analyzer, ['NaAlSiO4', 'SiO2', 'Al2O3']))
    assert len(ReferenceMatchingModule._sort_mixtures(mixtures[0])[0][cnst.REF_NAME_KEY]) == 3
    assert files == ['Mixture_plot_cl0.png']


def test_plot_best_mixture_defaults_to_true():
    assert PlotConfig().plot_best_mixture is True
    # Plot configs saved before the field existed
    assert PlotConfig.model_validate({'save_plots': True}).plot_best_mixture is True


# --- 2D clustering plot ----------------------------------------------------------
def _clustering_stub(tmp_path, els, features='at_fr', els_to_plot=()):
    return SimpleNamespace(
        all_els_sample=list(els), detectable_els_sample=list(els), sample_id='test', analysis_dir=str(tmp_path),
        clustering_cfg=ClusteringConfig(features=features), plot_cfg=PlotConfig(els_to_plot=list(els_to_plot)),
        ref_formulae=None, ref_phases_df=None,
    )


def _two_groups(els, seed=0):
    rng = np.random.default_rng(seed)
    d = len(els)
    a, b = np.full(d, 1 / d), np.r_[np.full(d - 1, 0.5 / (d - 1)), 0.5]
    X = np.r_[a + rng.normal(0, 0.005, (20, d)), b + rng.normal(0, 0.005, (20, d))]
    X = np.clip(X, 0, None); X /= X.sum(axis=1, keepdims=True)
    labels = np.r_[np.zeros(20, int), np.ones(20, int)]
    C = np.array([X[labels == j].mean(0) for j in range(2)])
    S = [X[labels == j].std(0) for j in range(2)]
    return pd.DataFrame(X, columns=els), C, labels, S


@pytest.mark.parametrize('els, features, expected, ternary', [
    (['Cu', 'Al', 'O'], 'at_fr', ['Cu', 'Al', 'O'], True),
    (['Cu', 'Al', 'O'], 'w_fr', None, False),
    (['Na', 'Al', 'Si', 'O'], 'at_fr', None, False),
])
def test_clustering_plot_2d(tmp_path, els, features, expected, ternary):
    stub = _clustering_stub(tmp_path, els, features)
    df, C, labels, S = _two_groups(els)
    chosen, is_ternary = PlottingModule._clustering_plot_2d_elements(stub, df, C)
    assert is_ternary == ternary
    assert chosen == expected if expected else len(chosen) == 2
    PlottingModule._save_clustering_plot_2d(stub, df, C, labels, S, [])
    assert sorted(os.listdir(tmp_path)) == ['Clustering_plot_2D.png', 'Clustering_plot_2D_zoomed.png']


def test_clustering_plot_2d_elements_follow_els_to_plot(tmp_path):
    stub = _clustering_stub(tmp_path, ['Na', 'Al', 'Si', 'O'], els_to_plot=['Na', 'O'])
    df, C, _, _ = _two_groups(['Na', 'Al', 'Si', 'O'])
    assert PlottingModule._clustering_plot_2d_elements(stub, df, C) == (['Na', 'O'], False)


def test_no_clustering_plot_2d_for_two_elements(tmp_path):
    stub = _clustering_stub(tmp_path, ['Ti', 'O'])
    df, C, labels, S = _two_groups(['Ti', 'O'])
    PlottingModule._save_clustering_plot_2d(stub, df, C, labels, S, [])
    assert os.listdir(tmp_path) == []


def test_plots_without_adjusttext(tmp_path, monkeypatch):
    """Labels keep their initial position when adjustText is not installed."""
    import builtins
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == 'adjustText':
            raise ImportError
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', fake_import)
    stub = _clustering_stub(tmp_path, ['Cu', 'Al', 'O'])
    df, C, labels, S = _two_groups(['Cu', 'Al', 'O'])
    PlottingModule._save_clustering_plot_2d(stub, df, C, labels, S, [])
    assert sorted(os.listdir(tmp_path)) == ['Clustering_plot_2D.png', 'Clustering_plot_2D_zoomed.png']
