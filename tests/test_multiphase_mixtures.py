#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the decomposition of clusters into mixtures of two or more candidate phases."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import autoemx.utils.constants as cnst
from autoemx.config.ledger_schemas import ClusteringConfig, MixtureParams
from autoemx.core.composition_analysis.reference_matching import ReferenceMatchingModule

ELS = ['Cu', 'Al', 'O']
# Atomic fractions of candidate phases and atoms per formula unit (weights for at_fr features)
REFS = {
    'CuO': ([0.5, 0.0, 0.5], 2),
    'Al2O3': ([0.0, 0.4, 0.6], 5),
    'Cu': ([1.0, 0.0, 0.0], 1),
    # Lies on the CuO-Al2O3 line (CuO + Al2O3)
    'CuAl2O4': ([1 / 7, 2 / 7, 4 / 7], 7),
}


def _make_analyzer(mixture=None, names=('CuO', 'Al2O3', 'Cu')):
    names = list(names)
    return SimpleNamespace(
        ref_formulae=names,
        ref_phases_df=pd.DataFrame([REFS[n][0] for n in names], columns=ELS),
        ref_weights_in_mixture=[REFS[n][1] for n in names],
        detectable_els_sample=ELS,
        clustering_cfg=ClusteringConfig(mixture=mixture or MixtureParams()),
        quant_cfg=SimpleNamespace(is_known_phase_mixture=False),
    )


def _mixture_cluster(phases, n=60, noise=0.005, seed=0):
    rng = np.random.default_rng(seed)
    H = np.array([REFS[p][0] for p in phases])
    W = rng.dirichlet(np.ones(len(phases)), size=n)
    X = np.clip(W @ H + rng.normal(0, noise, size=(n, len(ELS))), 0, None)
    return X / X.sum(axis=1, keepdims=True)


def _top_mixture(mixtures):
    return max(mixtures, key=lambda m: m[cnst.CONF_SCORE_KEY])


# --- Config model -----------------------------------------------------------
def test_mixture_params_defaults():
    params = ClusteringConfig().mixture
    assert params.max_n_phases == 4
    assert params.max_recon_error == pytest.approx(0.4)


@pytest.mark.parametrize("kwargs", [
    {"max_n_phases": 1}, {"max_recon_error": 0}, {"max_recon_error": float("inf")},
    {"conf_sigma": -1}, {"single_phase_min_ref_conf": 1.5}, {"nmf_min_mixture_conf": -0.1},
    {"equivalent_recon_error_tol": -0.01},
    {"max_reported_mixtures": 0}, {"report_within_conf_ratio": 0}, {"min_reported_conf_ratio": 1.2},
    {"equivalent_span_tol": 0},
])
def test_mixture_params_validation(kwargs):
    with pytest.raises(ValueError):
        MixtureParams(**kwargs)


def test_default_mixture_params_keep_fingerprint():
    """Existing configs must keep their fingerprint (no mixture block at default values)."""
    assert "mixture" not in ClusteringConfig().fingerprint_payload()
    binary_only = ClusteringConfig(mixture=MixtureParams(max_n_phases=2))
    assert binary_only.fingerprint_payload()["mixture"]["max_n_phases"] == 2
    assert binary_only.fingerprint() != ClusteringConfig().fingerprint()


# --- Mixture identification -------------------------------------------------
def test_binary_cluster_is_decomposed_into_two_phases():
    analyzer = _make_analyzer()
    X = _mixture_cluster(['CuO', 'Al2O3'])
    _, mixtures = ReferenceMatchingModule._identify_mixture_from_refs(analyzer, X)
    assert all(len(m[cnst.REF_NAME_KEY]) == 2 for m in mixtures)
    assert set(_top_mixture(mixtures)[cnst.REF_NAME_KEY]) == {'CuO', 'Al2O3'}


def test_ternary_cluster_is_decomposed_into_three_phases():
    analyzer = _make_analyzer()
    X = _mixture_cluster(['CuO', 'Al2O3', 'Cu'])
    max_conf, mixtures = ReferenceMatchingModule._identify_mixture_from_refs(analyzer, X)
    top = _top_mixture(mixtures)
    assert set(top[cnst.REF_NAME_KEY]) == {'CuO', 'Al2O3', 'Cu'}
    assert top[cnst.CONF_SCORE_KEY] == pytest.approx(max_conf)
    assert len(top[cnst.MOLAR_FRS_MEAN_KEY]) == 3
    assert sum(top[cnst.MOLAR_FRS_MEAN_KEY]) == pytest.approx(1)


def test_max_n_phases_two_restricts_to_binary_mixtures():
    analyzer = _make_analyzer(MixtureParams(max_n_phases=2))
    X = _mixture_cluster(['CuO', 'Al2O3', 'Cu'])
    _, mixtures = ReferenceMatchingModule._identify_mixture_from_refs(analyzer, X)
    assert all(len(m[cnst.REF_NAME_KEY]) == 2 for m in mixtures)


def test_molar_fractions_match_binary_formula():
    """The N-phase conversion reduces to the previous closed-form binary expression."""
    analyzer = _make_analyzer()
    rng = np.random.default_rng(1)
    c2 = rng.uniform(0, 1, 20)
    W = np.c_[1 - c2, c2]
    ref_weights = [2.0, 5.0]
    mix, _ = ReferenceMatchingModule._get_mixture_dict_with_conf(analyzer, W, ref_weights, 0.1, ['CuO', 'Al2O3'])
    r = ref_weights[0] / ref_weights[1]
    x1_legacy = 1 - c2 * r / (1 - c2 * (1 - r))
    assert mix[cnst.MOLAR_FR_MEAN_KEY] == pytest.approx(x1_legacy.mean())
    assert mix[cnst.MOLAR_FR_STDEV_KEY] == pytest.approx(x1_legacy.std())


def test_mixtures_df_reports_all_fractions_for_multiphase_mixtures():
    mixtures = [[
        {cnst.REF_NAME_KEY: ['CuO', 'Al2O3', 'Cu'], cnst.CONF_SCORE_KEY: 0.99, cnst.MOLAR_FR_MEAN_KEY: 0.5,
         cnst.MOLAR_FR_STDEV_KEY: 0.1, cnst.MOLAR_FRS_MEAN_KEY: [0.5, 0.3, 0.2], cnst.MOLAR_FRS_STDEV_KEY: [0.1, 0.1, 0.1]},
        {cnst.REF_NAME_KEY: ['CuO', 'Al2O3'], cnst.CONF_SCORE_KEY: 0.6, cnst.MOLAR_FR_MEAN_KEY: 0.6,
         cnst.MOLAR_FR_STDEV_KEY: 0.1},
    ]]
    df = ReferenceMatchingModule._build_mixtures_df(_make_analyzer(), mixtures)
    row = df.iloc[0]
    assert row[f'{cnst.MIX_DF_KEY}1'] == 'CuO, Al2O3, Cu'
    assert row[f'{cnst.MIX_ALL_COMPS_MEAN_DF_KEY}1'] == '0.50/0.30/0.20'
    assert np.isnan(row[f'{cnst.MIX_MOLAR_RATIO_DF_KEY}1'])
    assert row[f'{cnst.MIX_MOLAR_RATIO_DF_KEY}2'] == pytest.approx(1.5)


# --- Ranking of equivalent mixtures -----------------------------------------
def test_collinear_phases_prefer_closest_pair():
    """CuO+Al2O3 and CuAl2O4+Al2O3 reconstruct a CuAl2O4-Al2O3 mixture equally well; the closer pair is preferred."""
    analyzer = _make_analyzer(names=('CuO', 'Al2O3', 'CuAl2O4'))
    labels = np.zeros(60, dtype=int)
    X = _mixture_cluster(['CuAl2O4', 'Al2O3'])
    compositions_df = pd.DataFrame(X, columns=ELS)
    mixtures = ReferenceMatchingModule._assign_mixtures(
        analyzer, 1, labels, compositions_df, rms_dist_cluster=[0.1], max_raw_confs=[0.0], n_points_per_cluster=[60]
    )[0]
    top = ReferenceMatchingModule._sort_mixtures(mixtures)[0]
    assert set(top[cnst.REF_NAME_KEY]) == {'CuAl2O4', 'Al2O3'}
    assert top[cnst.MIX_RANK_KEY] == 1
    df = ReferenceMatchingModule._build_mixtures_df(analyzer, [mixtures])
    assert df.iloc[0][f'{cnst.MIX_DF_KEY}1'] in ('CuAl2O4, Al2O3', 'Al2O3, CuAl2O4')


def _mix(refs, recon_error, spread):
    return {cnst.REF_NAME_KEY: refs, cnst.CONF_SCORE_KEY: float(np.exp(-recon_error ** 2 / 0.5)),
            cnst.MIX_RECON_ERROR_KEY: recon_error, cnst.MIX_PHASES_SPREAD_KEY: spread}


def test_rank_mixtures_uses_tolerance():
    far = _mix(['Al2O3', 'Cu', 'CO'], 0.000, 3.2)
    close = _mix(['Al2O3', 'CuO', 'CO'], 0.003, 1.9)
    worse = _mix(['Al2O3', 'CuO', 'Cu'], 0.058, 2.6)
    ranked = ReferenceMatchingModule._rank_mixtures([far, worse, close], equivalent_recon_error_tol=0.01)
    assert [m[cnst.REF_NAME_KEY][1] for m in ranked] == ['CuO', 'Cu', 'CuO']
    assert [m[cnst.MIX_RANK_KEY] for m in ranked] == [1, 2, 3]
    # Without tolerance, the lowest reconstruction error wins
    ranked = ReferenceMatchingModule._rank_mixtures([far, worse, close], equivalent_recon_error_tol=0)
    assert ranked[0] is far


def test_rank_mixtures_prefers_fewer_phases_when_equivalent():
    ternary = _mix(['A', 'B', 'C'], 0.100, 0.5)
    binary = _mix(['A', 'B'], 0.105, 0.9)
    assert ReferenceMatchingModule._rank_mixtures([ternary, binary], 0.01)[0] is binary


def test_sort_mixtures_without_ranks_uses_confidence():
    """Results saved before mixtures were ranked are sorted by confidence score."""
    legacy = [{cnst.REF_NAME_KEY: ['A', 'B'], cnst.CONF_SCORE_KEY: 0.2},
              {cnst.REF_NAME_KEY: ['A', 'C'], cnst.CONF_SCORE_KEY: 0.9}]
    assert ReferenceMatchingModule._sort_mixtures(legacy)[0][cnst.CONF_SCORE_KEY] == 0.9


# --- Mixtures reported in Clusters.csv -----------------------------------------
def _conf_mixtures(confs):
    return [{cnst.REF_NAME_KEY: ['A', f'B{i}'], cnst.CONF_SCORE_KEY: c, cnst.MIX_RANK_KEY: i + 1,
             cnst.MOLAR_FR_MEAN_KEY: 0.5, cnst.MOLAR_FR_STDEV_KEY: 0.1} for i, c in enumerate(confs)]


def test_reported_mixtures_limited_to_five():
    mixtures = _conf_mixtures([1.0, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6])
    reported = ReferenceMatchingModule._select_reported_mixtures(mixtures, MixtureParams())
    assert [m[cnst.MIX_RANK_KEY] for m in reported] == [1, 2, 3, 4, 5]


def test_reported_mixtures_beyond_five_if_within_ten_percent():
    mixtures = _conf_mixtures([1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.91, 0.85])
    reported = ReferenceMatchingModule._select_reported_mixtures(mixtures, MixtureParams())
    assert len(reported) == 7


def test_reported_mixtures_never_below_half_of_best():
    mixtures = _conf_mixtures([0.9, 0.6, 0.44, 0.1])
    reported = ReferenceMatchingModule._select_reported_mixtures(mixtures, MixtureParams())
    assert [m[cnst.CONF_SCORE_KEY] for m in reported] == [0.9, 0.6]


def test_mixtures_df_notes_mixtures_only_in_ledger():
    analyzer = _make_analyzer()
    df = ReferenceMatchingModule._build_mixtures_df(analyzer, [_conf_mixtures([0.9, 0.6, 0.44, 0.1]), _conf_mixtures([0.9])])
    assert f'{cnst.MIX_DF_KEY}3' not in df.columns
    assert df.iloc[0][cnst.MIX_MORE_DF_KEY].startswith('2 more mixture(s)')
    assert '2 with confidence below 50% of the best (CS_mix < 0.45)' in df.iloc[0][cnst.MIX_MORE_DF_KEY]
    assert pd.isna(df.iloc[1][cnst.MIX_MORE_DF_KEY])


def test_mixtures_df_notes_mixtures_beyond_the_first_five():
    df = ReferenceMatchingModule._build_mixtures_df(_make_analyzer(), [_conf_mixtures([1.0, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6])])
    note = df.iloc[0][cnst.MIX_MORE_DF_KEY]
    assert note.startswith('2 more mixture(s)')
    assert '2 ranked after the first 5 with confidence more than 10% below the best (CS_mix < 0.90)' in note


# --- Equivalent decompositions ---------------------------------------------------
SR_TA_O = {  # at. fractions of (Sr, Ta, O): all on the SrO-Ta2O5 line
    'SrO': [0.5, 0, 0.5], 'Ta2O5': [0, 2 / 7, 5 / 7], 'Sr2Ta2O7': [2 / 11, 2 / 11, 7 / 11],
    'SrTa2O6': [1 / 9, 2 / 9, 6 / 9], 'Sr4Ta2O9': [4 / 15, 2 / 15, 9 / 15],
}


def test_collinear_pairs_collapse_to_one():
    comps = {f: np.array(v) for f, v in SR_TA_O.items()}
    mixtures = [{cnst.REF_NAME_KEY: list(p), cnst.CONF_SCORE_KEY: 0.95 - 0.01 * i, cnst.MIX_RANK_KEY: i + 1}
                for i, p in enumerate([('SrTa2O6', 'Sr4Ta2O9'), ('Ta2O5', 'SrO'), ('Sr2Ta2O7', 'SrTa2O6')])]
    kept, n_eq = ReferenceMatchingModule._collapse_equivalent_mixtures(mixtures, comps, 0.005)
    assert [m[cnst.MIX_RANK_KEY] for m in kept] == [1] and n_eq == 2


def test_ternaries_spanning_whole_space_are_kept():
    comps = {f: np.array(REFS[f][0]) for f in ('CuO', 'Al2O3', 'Cu')}
    comps['CO'] = np.array([0.0, 0.0, 1.0])
    mixtures = [{cnst.REF_NAME_KEY: ['Al2O3', 'Cu', 'CO'], cnst.CONF_SCORE_KEY: 1.0},
                {cnst.REF_NAME_KEY: ['CuO', 'Al2O3', 'CO'], cnst.CONF_SCORE_KEY: 0.99},
                {cnst.REF_NAME_KEY: ['CuO', 'Al2O3'], cnst.CONF_SCORE_KEY: 0.5}]
    kept, n_eq = ReferenceMatchingModule._collapse_equivalent_mixtures(mixtures, comps, 0.005)
    assert len(kept) == 3 and n_eq == 0


def test_mixtures_df_notes_equivalent_mixtures():
    els = ['Sr', 'Ta', 'O']
    names = list(SR_TA_O)
    analyzer = SimpleNamespace(ref_formulae=names, ref_phases_df=pd.DataFrame([SR_TA_O[n] for n in names], columns=els),
                               detectable_els_sample=els, clustering_cfg=ClusteringConfig())
    mixtures = [{cnst.REF_NAME_KEY: list(p), cnst.CONF_SCORE_KEY: 0.95, cnst.MIX_RANK_KEY: i + 1,
                 cnst.MOLAR_FR_MEAN_KEY: 0.5, cnst.MOLAR_FR_STDEV_KEY: 0.1}
                for i, p in enumerate([('SrTa2O6', 'Sr4Ta2O9'), ('Ta2O5', 'SrO'), ('Sr2Ta2O7', 'SrTa2O6')])]
    df = ReferenceMatchingModule._build_mixtures_df(analyzer, [mixtures])
    assert f'{cnst.MIX_DF_KEY}2' not in df.columns
    assert '2 equivalent to a listed mixture' in df.iloc[0][cnst.MIX_MORE_DF_KEY]



def test_violin_plots_only_for_reported_mixtures(monkeypatch):
    """Known powder mixtures: violin plots for the reported binary mixtures only; no per-point data left in the ledger."""
    from autoemx.core.composition_analysis.plotting import PlottingModule
    calls = []
    monkeypatch.setattr(PlottingModule, '_save_violin_plot_powder_mixture',
                        staticmethod(lambda self, W, names, cid: calls.append(tuple(names))))
    els = ['Sr', 'Ta', 'O']
    names = list(SR_TA_O)
    analyzer = SimpleNamespace(
        ref_formulae=names, ref_phases_df=pd.DataFrame([SR_TA_O[n] for n in names], columns=els),
        ref_weights_in_mixture=[2, 7, 11, 9, 15], detectable_els_sample=els, clustering_cfg=ClusteringConfig(),
        quant_cfg=SimpleNamespace(is_known_phase_mixture=True),
    )
    rng = np.random.default_rng(0)
    H = np.array([SR_TA_O['SrTa2O6'], SR_TA_O['Sr4Ta2O9']])
    X = rng.dirichlet(np.ones(2), size=40) @ H + rng.normal(0, 0.004, (40, 3))
    X = np.clip(X, 0, None); X /= X.sum(1, keepdims=True)
    mixtures = ReferenceMatchingModule._assign_mixtures(
        analyzer, 1, np.zeros(40, dtype=int), pd.DataFrame(X, columns=els),
        rms_dist_cluster=[1.0], max_raw_confs=[0.0], n_points_per_cluster=[40],
    )[0]
    assert len(mixtures) == 10  # all pairs are kept for known phase mixtures
    assert len(calls) == 1      # all pairs are collinear: only the top-ranked one is reported
    assert not any('_point_molar_frs' in m for m in mixtures)


def test_violin_plots_for_reported_three_phase_mixtures(monkeypatch):
    """Known phase mixtures: a reported mixture of three phases gets a violin plot per phase."""
    from autoemx.core.composition_analysis.plotting import PlottingModule
    calls = []
    monkeypatch.setattr(PlottingModule, '_save_violin_plot_powder_mixture',
                        staticmethod(lambda self, W, names, cid: calls.append(('binary', tuple(names)))))
    monkeypatch.setattr(PlottingModule, '_save_violin_plot_powder_mixture_multi',
                        staticmethod(lambda self, W, names, cid: calls.append(('multi', tuple(names), W.shape))))
    analyzer = _make_analyzer()
    analyzer.quant_cfg = SimpleNamespace(is_known_phase_mixture=True)
    X = _mixture_cluster(['CuO', 'Al2O3', 'Cu'])
    ReferenceMatchingModule._assign_mixtures(
        analyzer, 1, np.zeros(len(X), dtype=int), pd.DataFrame(X, columns=ELS),
        rms_dist_cluster=[1.0], max_raw_confs=[0.0], n_points_per_cluster=[len(X)],
    )
    assert ('multi', ('CuO', 'Al2O3', 'Cu'), (len(X), 3)) in calls
    assert not any(c[0] == 'binary' for c in calls)  # binaries fail and are not reported
