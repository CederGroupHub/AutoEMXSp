#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Acquisition loop with quantify=True, using a fake microscope and a stubbed quantification worker.

Quantification itself is covered elsewhere; this checks the acquisition bookkeeping:
ledger entries, per-particle quantification, low-count particle skip, and the
element-identification pool carried over from earlier particles.
"""

from __future__ import annotations

import json
from pathlib import Path

import autoemx.core.composition_analysis.analyser as analyser_mod
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PowderMeasurementConfig,
    QuantificationOptionsConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.config.schema_models.quantification import (
    QuantificationDiagnostics,
    QuantificationResult,
)
import autoemx.utils.constants as cnst

from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID

SAMPLE_ID = "acq_quant_test"
TARGET_COUNTS = 50000
# (particle id, number of spots, first spectrum has low counts)
PARTICLES = [(0, 2, False), (1, 2, True), (2, 2, False)]


def _read_msa(path: Path) -> tuple[list[str], list[float]]:
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("##SPECTRUM"))
    header = lines[: start + 1]
    counts = [float(line.split(",")[0]) for line in lines[start + 1:] if not line.startswith("#")]
    return header, counts


class _FakeEMController:
    """Serves one real K-412 spectrum, rescaled to the target (or a low) number of counts."""

    def __init__(self, *args, **kwargs):
        header, counts = _read_msa(sorted((INPUTS_DIR / K412_CLUSTER_MINI_ID / "spectra").glob("*.msa"))[0])
        self._header = header
        self._counts = [c * TARGET_COUNTS / sum(counts) for c in counts]
        self.current_frame_label = "A1"
        self.particle_finder = None
        self.pixel_size_um = 0.1
        self.measurement_cfg = args[0].measurement_cfg
        self._particle_idx = 0
        self._low_next = False

    def initialise_SEM(self): pass
    def initialise_sample_navigator(self, **kwargs): pass
    def initialise_XS_analyzer(self): pass
    def save_frame_image(self, *args, **kwargs): pass

    def get_XSp_coords(self, next_spectrum_id):
        if self._particle_idx >= len(PARTICLES):
            return False, [], None
        particle_id, n_spots, low_counts = PARTICLES[self._particle_idx]
        self._particle_idx += 1
        self._low_next = low_counts
        return True, [(100 + 10 * j, 100) for j in range(n_spots)], particle_id

    def acquire_XS_spot_spectrum(self, x, y, max_time, target_counts, elements=None, msa_file_path=None):
        scale = 0.3 if self._low_next else 1.0
        self._low_next = False
        counts = [c * scale for c in self._counts]
        Path(msa_file_path).write_text(
            "\n".join(self._header) + "\n" + "\n".join(map(str, counts)) + "\n#ENDOFDATA   :\n"
        )
        return counts, None


def test_acquisition_with_quantification(tmp_path: Path, monkeypatch):
    received_pools = {}

    def _stub_worker(payload):
        idx = int(payload["spectrum_index"])
        received_pools[idx] = list(payload["identification_pool"])
        comp = {"Si": 0.4, "O": 0.6}
        record = QuantificationResult(
            quantification_id=payload["quantification_id"],
            quant_flag=0,
            composition_atomic_fractions=comp,
            # The first spectrum "identifies" Ti, which later particles must test first
            diagnostics=QuantificationDiagnostics(added_elements=["Ti"] if idx == 0 else []),
        )
        return idx, {cnst.COMP_AT_FR_KEY: comp, cnst.AN_ER_KEY: 0.0}, record, None

    monkeypatch.setattr(analyser_mod, "EM_Controller", _FakeEMController)
    monkeypatch.setattr(analyser_mod, "_quantify_spectrum_worker", _stub_worker)
    monkeypatch.setattr(analyser_mod, "_resolve_parallel_mode", lambda requested_cores: (1, None, False))
    monkeypatch.setattr(
        analyser_mod.EMXSp_Composition_Analyzer,
        "_evaluate_clustering_convergence",
        lambda self, tot_n_spectra, n_spectra_to_collect: (tot_n_spectra, False, False),
    )

    analyzer = analyser_mod.EMXSp_Composition_Analyzer(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"),
        sample_cfg=SampleConfig(elements=["Si", "O"], type="powder"),
        measurement_cfg=MeasurementConfig(
            min_n_spectra=4,
            max_n_spectra=5,
            target_acquisition_counts=TARGET_COUNTS,
        ),
        sample_substrate_cfg=SampleSubstrateConfig(elements=["C"], auto_detection=False),
        quant_cfg=QuantificationOptionsConfig(identification_seed_spectra=3),
        powder_meas_cfg=PowderMeasurementConfig(max_spectra_per_par=2),
        is_acquisition=True,
        verbose=False,
        results_dir=str(tmp_path),
        sample_id=SAMPLE_ID,
    )
    analyzer.run_collection_and_quantification(quantify=True)

    ledger = json.loads((tmp_path / SAMPLE_ID / "ledger.json").read_text())
    spectra = ledger["spectra"]

    # Particle 1 is skipped after its first (low-count) spectrum
    assert [s["acquisition_details"]["particle_id"] for s in spectra] == [0, 0, 1, 2, 2]
    # Every spectrum is quantified exactly once, in the same quantification run
    assert all(len(s["quantification_results"]) == 1 for s in spectra)
    assert {s["quantification_results"][0]["quantification_id"] for s in spectra} == {0}
    # Pool is built from the first 3 spectra (across particles) and passed to the later ones
    assert received_pools == {0: [], 1: [], 2: [], 3: ["Ti"], 4: ["Ti"]}
