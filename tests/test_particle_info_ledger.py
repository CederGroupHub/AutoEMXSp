#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for ledger ParticleInfo geometry and clustering summary fields."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from autoemx.config.ledger_schemas import (
    AcquisitionDetails,
    LedgerConfigs,
    ParticleInfo,
    SampleLedger,
    SpectrumEntry,
)
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.config.schema_models.clustering import ClusteringResult
from autoemx.core.composition_analysis.analyser import EMXSp_Composition_Analyzer
import autoemx.utils.constants as cnst


def _minimal_ledger_configs() -> LedgerConfigs:
    return LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"),
        sample_cfg=SampleConfig(elements=["Pb", "Mo", "O"], type="powder"),
        measurement_cfg=MeasurementConfig(),
        sample_substrate_cfg=SampleSubstrateConfig(),
        plot_cfg=PlotConfig(),
    )


def _spectrum(particle_id: int, spectrum_id: str = "0") -> SpectrumEntry:
    return SpectrumEntry(
        spectrum_id=spectrum_id,
        total_counts=1000,
        live_acquisition_time=1.0,
        acquisition_details=AcquisitionDetails(particle_id=particle_id, frame_id="A1"),
        spectrum_relpath=f"spectra/spectrum_{spectrum_id}.msa",
    )


def test_particle_info_round_trip(outputs_dir: Path):
    sample_dir = outputs_dir / "particle_info_round_trip"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ledger = SampleLedger(
        sample_id="demo",
        sample_path=str(sample_dir.resolve()),
        configs=_minimal_ledger_configs(),
        spectra=[_spectrum(1)],
        particles=[
            ParticleInfo(
                id=1,
                area_um=3.14,
                eq_diameter_um=2.0,
                clusters=[0, 2],
                composition="Candidate phases: PbMoO4 (CScnd: 0.90)",
            )
        ],
    )
    path = sample_dir / "ledger.json"
    ledger.to_json_file(path)

    loaded = SampleLedger.from_json_file(path)
    assert len(loaded.particles) == 1
    particle = loaded.particles[0]
    assert particle.id == 1
    assert particle.area_um == pytest.approx(3.14)
    assert particle.eq_diameter_um == pytest.approx(2.0)
    assert particle.clusters == [0, 2]
    assert "PbMoO4" in (particle.composition or "")


def test_legacy_ledger_loads_with_empty_particles():
    legacy_path = Path(__file__).parent / "inputs" / "Wulfenite_mini" / "ledger.json"
    ledger = SampleLedger.from_json_file(legacy_path)
    assert ledger.particles == []


def test_duplicate_particle_ids_rejected(outputs_dir: Path):
    sample_dir = outputs_dir / "particle_info_dup"
    sample_dir.mkdir(parents=True, exist_ok=True)
    with pytest.raises(ValueError, match="ParticleInfo.id values must be unique"):
        SampleLedger(
            sample_id="demo",
            sample_path=str(sample_dir.resolve()),
            configs=_minimal_ledger_configs(),
            spectra=[_spectrum(1)],
            particles=[
                ParticleInfo(id=1, area_um=1.0, eq_diameter_um=1.0),
                ParticleInfo(id=1, area_um=2.0, eq_diameter_um=1.5),
            ],
        )


def test_eq_diameter_from_area():
    area = math.pi  # diameter 2
    assert EMXSp_Composition_Analyzer._eq_diameter_um_from_area(area) == pytest.approx(2.0)


def test_composition_string_template():
    result = ClusteringResult(
        centroids=[[0.5], [0.5]],
        wcss=1.0,
        sil_score=0.5,
        tot_n_points=4,
        refs_assigned_rows=[
            {"Cnd1": "PbMoO4", "CS_cnd1": 0.9},
            {"Cnd1": "PbMoO3", "CS_cnd1": 0.5},
        ],
        clusters_assigned_mixtures=[
            [{"refs": ["PbMoO4", "SiO2"], "conf_score": 0.8, "mean": 0.6, "stddev": 0.01}],
            [],
        ],
    )
    text = EMXSp_Composition_Analyzer._format_particle_composition_string([0, 1], result)
    assert text == (
        "Candidate phases: PbMoO4 (CScnd: 0.90), PbMoO3 (CScnd: 0.50); "
        "candidate mixtures: PbMoO4+SiO2 (CSmix: 0.80)"
    )


def test_extract_spectrum_info_inserts_size_columns():
    analyser = object.__new__(EMXSp_Composition_Analyzer)
    spectrum = _spectrum(7, spectrum_id="3")
    particle = ParticleInfo(id=7, area_um=4.0129, eq_diameter_um=2.256758334)
    row = analyser._extract_spectrum_info(
        spectrum, 3, particles_by_id={7: particle}
    )
    keys = list(row.keys())
    assert keys[:4] == [
        cnst.SP_ID_DF_KEY,
        cnst.PAR_ID_DF_KEY,
        cnst.PAR_AREA_UM_KEY,
        cnst.PAR_EQ_D_KEY,
    ]
    assert keys[4] == cnst.FRAME_ID_DF_KEY
    assert row[cnst.PAR_AREA_UM_KEY] == 4.013
    assert row[cnst.PAR_EQ_D_KEY] == 2.257


def test_extract_spectrum_info_omits_size_columns_without_particles():
    analyser = object.__new__(EMXSp_Composition_Analyzer)
    row = analyser._extract_spectrum_info(_spectrum(1), 0, particles_by_id=None)
    assert cnst.PAR_AREA_UM_KEY not in row
    assert cnst.PAR_EQ_D_KEY not in row


def test_update_particles_from_clustering(outputs_dir: Path):
    analyser = object.__new__(EMXSp_Composition_Analyzer)
    sample_dir = outputs_dir / "particle_info_clustering"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ledger = SampleLedger(
        sample_id="demo",
        sample_path=str(sample_dir.resolve()),
        configs=_minimal_ledger_configs(),
        spectra=[
            _spectrum(1, "0"),
            _spectrum(1, "1"),
            _spectrum(2, "2"),
        ],
        particles=[
            ParticleInfo(id=1, area_um=1.0, eq_diameter_um=1.128),
            ParticleInfo(id=2, area_um=2.0, eq_diameter_um=1.595),
        ],
    )
    result = ClusteringResult(
        centroids=[[0.5], [0.5]],
        wcss=1.0,
        sil_score=0.4,
        tot_n_points=3,
        refs_assigned_rows=[
            {"Cnd1": "PbMoO4", "CS_cnd1": 0.9},
            {"Cnd1": "SiO2", "CS_cnd1": 0.7},
        ],
        clusters_assigned_mixtures=[[], []],
    )
    # spectrum 0 -> cluster 0, spectrum 1 -> cluster 1, spectrum 2 -> cluster 0
    changed = analyser._update_particles_from_clustering(
        ledger,
        labels=[0, 1, 0],
        df_indices=[0, 1, 2],
        clustering_result=result,
    )
    assert changed
    assert ledger.particles[0].clusters == [0, 1]
    assert "PbMoO4" in (ledger.particles[0].composition or "")
    assert "SiO2" in (ledger.particles[0].composition or "")
    assert ledger.particles[1].clusters == [0]
    assert ledger.particles[1].composition == "Candidate phases: PbMoO4 (CScnd: 0.90)"


def test_update_particles_noop_when_empty(outputs_dir: Path):
    analyser = object.__new__(EMXSp_Composition_Analyzer)
    sample_dir = outputs_dir / "particle_info_noop"
    sample_dir.mkdir(parents=True, exist_ok=True)
    ledger = SampleLedger(
        sample_id="demo",
        sample_path=str(sample_dir.resolve()),
        configs=_minimal_ledger_configs(),
        spectra=[_spectrum(1)],
        particles=[],
    )
    result = ClusteringResult(centroids=[], wcss=0.0, sil_score=0.0, tot_n_points=0)
    assert not analyser._update_particles_from_clustering(
        ledger, labels=[0], df_indices=[0], clustering_result=result
    )
