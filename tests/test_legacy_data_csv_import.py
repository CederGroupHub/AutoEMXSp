#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Results saved in the legacy format (Data.csv + Comp_analysis_configs.json) are converted to a ledger.

The fixture holds the first 7 spectra of a Wulfenite measurement (AutoEMX 2025): 6 accepted (quant_flag 0) and
1 flagged for a poor fit (quant_flag 4).
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from autoemx._version import LEGACY_VERSION
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.runners.analyze_sample import analyze_sample
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, WULFENITE_LEGACY_MINI_ID


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    sample_dir = tmp_path / WULFENITE_LEGACY_MINI_ID
    shutil.copytree(INPUTS_DIR / WULFENITE_LEGACY_MINI_ID, sample_dir)
    return sample_dir


def _ledger_path(sample_dir: Path) -> Path:
    return sample_dir / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)


def test_legacy_results_converted_to_ledger(sample_dir: Path):
    data = pd.read_csv(sample_dir / "Data.csv")
    ledger = load_sample_ledger(_ledger_path(sample_dir))

    assert _ledger_path(sample_dir).exists()
    # Configurations
    assert ledger.configs.sample_cfg.elements == ["Pb", "Mo", "O"]
    assert ledger.configs.sample_substrate_cfg.elements == ["C", "O", "Al"]
    assert ledger.configs.measurement_cfg.powder_meas_cfg.par_selection_mode == "auto"  # legacy key converted
    # Quantification run: computed by an AutoEMX version that did not record it, with the legacy settings
    [quant_config] = ledger.quantifications
    assert ledger.active_quant == quant_config.quantification_id == 0
    assert quant_config.autoemx_version == LEGACY_VERSION
    assert quant_config.options["fit_tolerance"] == 0.001
    assert quant_config.options["min_bckgrnd_cnts"] == 5.0  # stored in the legacy quant_cfg
    assert quant_config.options["is_known_phase_mixture"] is False  # stored in the legacy powder settings
    clustering_config = quant_config.clustering_analyses[0].config
    assert clustering_config.ref_formulae == ["PbMoO4"]
    assert clustering_config.geometry == "euclidean" and clustering_config.auto_merge_clusters is False

    # Spectra: one file and entry per row, with the results of Data.csv
    assert [s.spectrum_id for s in ledger.spectra] == [str(i) for i in data["Spectrum #"]]
    assert [s.acquisition_details.particle_id for s in ledger.spectra] == list(data["Particle #"])
    for spectrum, (_, row) in zip(ledger.spectra, data.iterrows()):
        assert (sample_dir / spectrum.spectrum_relpath).exists()
        [result] = spectrum.quantification_results
        assert result.quant_flag == row["Quant_flag"]
        assert result.analytical_error == pytest.approx(row["An er w%"] / 100)
        for el in ("Pb", "Mo", "O"):
            assert result.composition_atomic_fractions[el] == pytest.approx(row[f"{el}_at%"] / 100)
    assert ledger.spectra[6].quantification_results[0].comment == "poor fit"
    # Accepted spectra: the comment "<value> min. ref. bckgrnd counts" becomes a diagnostic
    assert ledger.spectra[0].quantification_results[0].diagnostics.min_background_ref_lines == pytest.approx(34.1)


def test_converted_ledger_is_loaded_unchanged(sample_dir: Path):
    load_sample_ledger(_ledger_path(sample_dir))
    first = json.loads(_ledger_path(sample_dir).read_text())
    load_sample_ledger(_ledger_path(sample_dir))
    assert json.loads(_ledger_path(sample_dir).read_text()) == first


def test_converted_ledger_can_be_analysed(sample_dir: Path):
    analyze_sample(WULFENITE_LEGACY_MINI_ID, results_path=str(sample_dir.parent), show_plots=False)
    [analysis_dir] = sample_dir.glob("analysis*")
    clusters = pd.read_csv(analysis_dir / "Clusters.csv", index_col="Metric")
    assert list(clusters.columns) == ["Cluster_0"]
    assert int(clusters.loc["n_points", "Cluster_0"]) == 6  # the spectrum flagged for a poor fit is excluded
