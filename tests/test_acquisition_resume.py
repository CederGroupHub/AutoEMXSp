#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Restarting an acquisition (quantify=True) on the same sample, with a fake microscope and a stubbed quantification.

The resumed acquisition must keep the existing spectra and results, continue spectrum and particle ids, reuse the
quantification run and quantify only the new spectra. Deleted spectrum files are dropped from the ledger only if not
yet quantified.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import autoemx.core.composition_analysis.analyser as analyser_mod
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.config.schema_models.quantification import QuantificationResult
import autoemx.utils.constants as cnst

from test_acquisition_with_quantification import _FakeEMController

SAMPLE_ID = "resume_test"
SPECTRA_PER_PARTICLE = 2


class _FakeMicroscope(_FakeEMController):
    """Fake microscope finding particles one after the other; like the particle finder, ids start at 1."""

    def get_XSp_coords(self, next_spectrum_id):
        self._particle_idx += 1
        self._low_next = False
        return True, [(100 + 10 * j, 100) for j in range(SPECTRA_PER_PARTICLE)], self._particle_idx


@pytest.fixture
def quantified_indices(monkeypatch):
    """Patch the microscope and the quantification; record the spectra quantified in each run."""
    runs = []

    def _stub_worker(payload):
        runs[-1].append(int(payload["spectrum_index"]))
        comp = {"Si": 0.4, "O": 0.6}
        record = QuantificationResult(
            quantification_id=payload["quantification_id"], quant_flag=0, composition_atomic_fractions=comp
        )
        return int(payload["spectrum_index"]), {cnst.COMP_AT_FR_KEY: comp, cnst.AN_ER_KEY: 0.0}, record, None

    monkeypatch.setattr(analyser_mod, "EM_Controller", _FakeMicroscope)
    monkeypatch.setattr(analyser_mod, "_quantify_spectrum_worker", _stub_worker)
    monkeypatch.setattr(analyser_mod, "_resolve_parallel_mode", lambda requested_cores: (1, None, False))
    monkeypatch.setattr(
        analyser_mod.EMXSp_Composition_Analyzer,
        "_evaluate_clustering_convergence",
        lambda self, tot_n_spectra, n_spectra_to_collect: (tot_n_spectra, False, False),
    )
    return runs


def _acquire(results_dir: Path, n_spectra: int, runs: list) -> None:
    runs.append([])
    analyzer = analyser_mod.EMXSp_Composition_Analyzer(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"),
        sample_cfg=SampleConfig(elements=["Si", "O"], type="powder"),
        measurement_cfg=MeasurementConfig(min_n_spectra=n_spectra, max_n_spectra=n_spectra),
        sample_substrate_cfg=SampleSubstrateConfig(elements=["C"], auto_detection=False),
        powder_meas_cfg=PowderMeasurementConfig(max_spectra_per_par=SPECTRA_PER_PARTICLE),
        is_acquisition=True,
        verbose=False,
        results_dir=str(results_dir),
        sample_id=SAMPLE_ID,
    )
    analyzer.run_collection_and_quantification(quantify=True)


def _ledger_path(results_dir: Path) -> Path:
    return results_dir / SAMPLE_ID / (cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)


def test_resumed_acquisition_continues_where_it_stopped(tmp_path: Path, quantified_indices):
    _acquire(tmp_path, n_spectra=4, runs=quantified_indices)
    _acquire(tmp_path, n_spectra=6, runs=quantified_indices)

    ledger = json.loads(_ledger_path(tmp_path).read_text())
    spectra = ledger["spectra"]
    assert [s["spectrum_id"] for s in spectra] == ["0", "1", "2", "3", "4", "5"]
    assert [s["acquisition_details"]["particle_id"] for s in spectra] == [1, 1, 2, 2, 3, 3]
    # Same quantification run, one result per spectrum, and only the new spectra quantified when resuming
    assert [q["quantification_id"] for q in ledger["quantifications"]] == [0]
    assert all([r["quantification_id"] for r in s["quantification_results"]] == [0] for s in spectra)
    assert quantified_indices == [[0, 1, 2, 3], [4, 5]]


def test_resume_drops_deleted_spectra_only_if_not_quantified(tmp_path: Path, quantified_indices):
    _acquire(tmp_path, n_spectra=4, runs=quantified_indices)

    # Spectrum 1: deleted before being quantified (e.g. a wrong acquisition). Spectrum 2: deleted after quantification.
    ledger_path = _ledger_path(tmp_path)
    payload = json.loads(ledger_path.read_text())
    payload["spectra"][1]["quantification_results"] = []
    ledger_path.write_text(json.dumps(payload))
    for spectrum_id in ("1", "2"):
        (tmp_path / SAMPLE_ID / cnst.SPECTRA_DIR / f"{cnst.SPECTRUM_FILENAME_PREFIX}{spectrum_id}.msa").unlink()

    _acquire(tmp_path, n_spectra=5, runs=quantified_indices)

    spectra = json.loads(ledger_path.read_text())["spectra"]
    assert [s["spectrum_id"] for s in spectra] == ["0", "2", "3", "4", "5"]
    spectrum_2 = spectra[1]
    assert [r["quant_flag"] for r in spectrum_2["quantification_results"]] == [0]  # result kept, not re-flagged
    assert quantified_indices[1] == [3, 4]  # the two new spectra (indices in the updated ledger)
