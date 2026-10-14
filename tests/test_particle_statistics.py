#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Particle size statistics and their export (autoemx.core.particle_statistics)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import create_sample_ledger
from autoemx.config.ledger_schemas import (
    LedgerConfigs,
    ParticleInfo,
    ParticleSizeBand,
    ParticleStatsFrame,
    ParticleStatsProgress,
    SampleLedger,
)
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.particle_statistics import (
    _weighted_percentile,
    bootstrap_confidence_intervals,
    save_particle_size_histogram,
    band_particle_weights,
    build_particle_stats_result,
    edge_correction_factor,
    compute_particle_size_stats,
    save_particle_statistics,
    select_size_stats_particles,
)


def test_summary_statistics():
    stats = compute_particle_size_stats([1.0, 2.0, 3.0, 4.0, 10.0])
    assert stats["n_par_analysed"] == 5
    assert stats["mean"] == pytest.approx(4.0)
    assert stats["stdev"] == pytest.approx(np.std([1.0, 2.0, 3.0, 4.0, 10.0]))  # population
    assert stats["median"] == 3.0 and stats["min"] == 1.0 and stats["max"] == 10.0
    assert stats["D10"] == pytest.approx(1.4) and stats["D90"] == pytest.approx(7.6)
    with pytest.raises(ValueError):
        compute_particle_size_stats([])


def test_saved_files(tmp_path: Path):
    particles = [
        ParticleInfo(id=0, area_um=float(np.pi), eq_diameter_um=2.0, frame_id="A0"),
        ParticleInfo(id=1, area_um=float(4 * np.pi), frame_id=None),  # diameter from the area
    ]
    stats = save_particle_statistics(particles, str(tmp_path), "S", output_file_suffix="_x")
    assert list(stats.columns)[:2] == ["Measure", "n_par_analysed"]
    assert stats["mean"].iloc[0] == pytest.approx(3.0)

    sizes = pd.read_csv(tmp_path / f"S_{cnst.PARTICLE_SIZES_FILENAME}_x.csv", keep_default_na=False)
    assert list(sizes[cnst.FRAME_ID_DF_KEY]) == ["A0", ""]
    assert list(sizes[cnst.PAR_EQ_D_KEY]) == pytest.approx([2.0, 4.0])
    assert (tmp_path / f"S_{cnst.PARTICLE_STATS_FILENAME}_x.csv").exists()
    assert (tmp_path / f"S_{cnst.PARTICLE_STAT_HIST_FILENAME}_x.png").exists()
    with pytest.raises(ValueError):
        save_particle_statistics([], str(tmp_path), "S")


PARTICLES = [
    ParticleInfo(id=0, area_um=10.0, selection="size_stats"),
    ParticleInfo(id=1, area_um=20.0, selection="auto"),
    ParticleInfo(id=2, area_um=30.0, selection="manual"),
    ParticleInfo(id=3, area_um=40.0, selection="list"),
    ParticleInfo(id=4, area_um=50.0),  # selection not recorded (older ledgers)
    ParticleInfo(id=5, selection="auto"),  # no area
]


def test_particles_selected_for_size_statistics():
    assert [p.id for p in select_size_stats_particles(PARTICLES)] == [0, 1]
    assert [p.id for p in select_size_stats_particles(PARTICLES, include_eds_auto_particles=False)] == [0]


def test_result_stored_in_the_ledger(tmp_path: Path):
    cfg = PowderMeasurementConfig(min_area_par=5.0)
    result = build_particle_stats_result(select_size_stats_particles(PARTICLES), cfg, n_par_target=100)
    assert result.particle_ids == [0, 1]
    assert result.n_par_by_selection == {"size_stats": 1, "auto": 1}
    assert result.n_par_target == 100 and result.powder_meas_cfg.min_area_par == 5.0
    assert result.size_stats.median == pytest.approx(np.mean([2 * np.sqrt(a / np.pi) for a in (10.0, 20.0)]))
    cfg.min_area_par = 1.0
    assert result.powder_meas_cfg.min_area_par == 5.0  # settings copied

    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=[], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    ledger = create_sample_ledger(str(tmp_path), "S", configs)
    ledger.particles = list(PARTICLES)
    ledger.particle_stats = result
    ledger.to_json_file(str(tmp_path / "ledger.json"))
    reloaded = SampleLedger.from_json_file(str(tmp_path / "ledger.json"))
    assert reloaded.particle_stats == result
    assert [p.selection for p in reloaded.particles] == ["size_stats", "auto", "manual", "list", None, "auto"]

    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, frame_width_um=50.0, n_frames=2, scanned_area_um2=3e3,
                            particle_ids=[0])
    [weighted_band] = build_particle_stats_result(PARTICLES[:2], cfg, bands=[band]).bands
    assert weighted_band.model_copy(update={"number_density_per_mm2": None, "n_eds_particles": 0}) == band
    assert weighted_band.n_eds_particles == 1  # particle 1 (EDS, 5 µm) is in the band
    assert weighted_band.number_density_per_mm2 == pytest.approx(1 / 3e3 * 1e6)  # 1 particle found in 3000 µm²
    with pytest.raises(ValueError, match="one size band"):
        build_particle_stats_result(PARTICLES[:2], cfg, bands=[band, band])

    with pytest.raises(ValueError, match="not in the ledger"):
        SampleLedger.model_validate({**ledger.model_dump(), "particles": PARTICLES[1:]})
    ledger.particle_stats = build_particle_stats_result(PARTICLES[:2], cfg, bands=[band.model_copy(update={"particle_ids": [0, 9]})])
    with pytest.raises(ValueError, match=r"not in the ledger: \[9\]"):
        SampleLedger.model_validate(ledger.model_dump())



# --- Weighting by size band ----------------------------------------------------------------------------------------

def _psd_particles(first_id: int, n: int, diameter_um: float) -> list:
    return [ParticleInfo(id=i, area_um=float(np.pi * diameter_um**2 / 4), eq_diameter_um=diameter_um,
                         selection="size_stats") for i in range(first_id, first_id + n)]


def _two_band_measurement(n_large: int, n_small: int):
    """
    Sample with 1000 particles of 5 µm and 10 particles of 50 µm per mm²: the large band scans the area holding
    n_large large particles, the small band the (100 times smaller) area holding n_small small ones.
    """
    large, small = _psd_particles(0, n_large, 50.0), _psd_particles(n_large, n_small, 5.0)
    bands = [
        ParticleSizeBand(d_min_um=10.0, d_max_um=100.0, n_frames=1, scanned_area_um2=n_large / 10 * 1e6,
                         particle_ids=[p.id for p in large]),
        ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=1, scanned_area_um2=n_small / 1000 * 1e6,
                         particle_ids=[p.id for p in small]),
    ]
    return large + small, bands


def test_band_weights_give_the_sample_distribution_whatever_the_number_collected():
    particles, bands = _two_band_measurement(n_large=10, n_small=100)
    result = build_particle_stats_result(particles, PowderMeasurementConfig(), bands=bands)
    # 1 % of the particles are large: D10-D90 are small particles, the mean is weighted accordingly
    assert result.size_stats.median == result.size_stats.D90 == 5.0
    assert result.size_stats.mean == pytest.approx((1000 * 5.0 + 10 * 50.0) / 1010)
    assert result.unweighted_size_stats.mean == pytest.approx((100 * 5.0 + 10 * 50.0) / 110)  # 9 % large, as before
    assert [b.number_density_per_mm2 for b in result.bands] == pytest.approx([10.0, 1000.0])

    # Collecting 3 times more large particles (over 3 times the area) does not change the statistics
    particles_3x, bands_3x = _two_band_measurement(n_large=30, n_small=100)
    other = build_particle_stats_result(particles_3x, PowderMeasurementConfig(), bands=bands_3x)
    for key in ("mean", "stdev", "median", "D10", "D25", "D75", "D90"):
        assert getattr(other.size_stats, key) == pytest.approx(getattr(result.size_stats, key))


def test_single_band_statistics_are_unweighted():
    particles = [ParticleInfo(id=i, area_um=float(a), selection="size_stats") for i, a in enumerate([3, 8, 20, 50])]
    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=4, scanned_area_um2=1e4, particle_ids=[0, 1, 2, 3])
    result = build_particle_stats_result(particles, PowderMeasurementConfig(), bands=[band])
    assert result.size_stats == result.unweighted_size_stats
    d = np.sort(np.random.default_rng(0).lognormal(1.0, 0.5, 50))
    for q in (10, 25, 50, 75, 90):  # equal weights: as numpy, method 'hazen'
        assert _weighted_percentile(d, np.full(d.size, 0.3), q) == pytest.approx(np.percentile(d, q, method="hazen"))
    # A heavy last value counts in proportion to its weight
    assert _weighted_percentile(np.array([1.0, 10.0]), np.array([1.0, 9.0]), 60) == 10.0  # its middle is at 55 %
    assert _weighted_percentile(np.array([1.0, 10.0]), np.array([1.0, 9.0]), 50) == pytest.approx(1 + 9 * 0.45 / 0.5)


def test_eds_particles_share_the_weight_of_their_band():
    particles, bands = _two_band_measurement(n_large=10, n_small=100)
    eds = [
        ParticleInfo(id=1000, area_um=float(np.pi * 30.0**2 / 4), selection="auto"),  # 30 µm: large band
        ParticleInfo(id=1001, area_um=float(np.pi * 200.0**2 / 4), selection="auto"),  # 200 µm: no band
    ]
    weights, weighted_bands = band_particle_weights(particles + eds, bands)
    assert [b.n_eds_particles for b in weighted_bands] == [1, 0]
    large_band = [k for k, p in enumerate(particles + eds) if p.id < 10 or p.id == 1000]
    assert weights[large_band].sum() == pytest.approx(10 / 1e6)  # band weight unchanged: its density per µm²
    assert weights[-1] == 0

    result = build_particle_stats_result(particles + eds, PowderMeasurementConfig(), bands=bands)
    assert result.n_par_without_weight == 1 and result.size_stats.n_par_analysed == 111
    assert result.size_stats.max == 50.0  # the 200 µm particle is left out



def test_frame_records_and_progress_in_the_ledger(tmp_path: Path):
    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=0, scanned_area_um2=0.0)
    for i, ids in enumerate([[0], []]):  # frames without particles count too
        band.add_frame(ParticleStatsFrame(frame_id=f"A{i}", frame_width_um=100.0, pixel_size_um=0.1,
                                          area_um2=6000.0, n_particles=len(ids)), ids)
    assert (band.n_frames, band.scanned_area_um2, band.particle_ids) == (2, 12000.0, [0])
    assert band.frame_area_um2(["A1", "B9"]) == 6000.0
    with pytest.raises(ValueError, match="n_frames"):
        ParticleSizeBand.model_validate({**band.model_dump(), "n_frames": 3})

    # Leaving out the empty frame A1 halves the area: the density doubles
    weights, [weighted] = band_particle_weights(PARTICLES[:1], [band], excluded_frame_ids=["A1"])
    assert weighted.number_density_per_mm2 == pytest.approx(1 / 6000 * 1e6)

    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=[], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    ledger = create_sample_ledger(str(tmp_path), "S", configs)
    ledger.particle_stats_progress = ParticleStatsProgress(n_par_target=10, powder_meas_cfg=PowderMeasurementConfig(),
                                                           bands=[band])
    with pytest.raises(ValueError, match=r"particle_stats_progress references particle ids not in the ledger: \[0\]"):
        SampleLedger.model_validate(ledger.model_dump())
    ledger.particles = list(PARTICLES[:1])
    ledger.to_json_file(str(tmp_path / "ledger.json"))
    assert SampleLedger.from_json_file(str(tmp_path / "ledger.json")).particle_stats_progress == ledger.particle_stats_progress



def test_edge_correction():
    # 100 x 50 px frame, 3 px margin: a 1 px particle has 92 x 42 allowed positions out of 100 x 50
    assert edge_correction_factor(1, 1, 100, 50, 3) == pytest.approx(100 * 50 / (92 * 42))
    assert edge_correction_factor(40, 20, 100, 50, 3) > edge_correction_factor(10, 5, 100, 50, 3) > 1
    with pytest.raises(ValueError):
        edge_correction_factor(93, 10, 100, 50, 3)

    # In a band, particles count (and share its weight) in proportion to their edge correction
    small = ParticleInfo(id=0, area_um=float(np.pi), eq_diameter_um=2.0, selection="size_stats", edge_correction=1.0)
    large = ParticleInfo(id=1, area_um=float(np.pi * 16), eq_diameter_um=8.0, selection="size_stats",
                         edge_correction=3.0)
    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=1, scanned_area_um2=1e3, particle_ids=[0, 1])
    weights, [weighted] = band_particle_weights([small, large], [band])
    assert weights[1] / weights[0] == pytest.approx(3.0)
    assert weighted.number_density_per_mm2 == pytest.approx(4 / 1e3 * 1e6)
    result = build_particle_stats_result([small, large], PowderMeasurementConfig(), bands=[band])
    assert result.size_stats.mean == pytest.approx((2.0 + 3 * 8.0) / 4)
    assert result.unweighted_size_stats.mean == pytest.approx(5.0)



def test_volume_weighted_statistics():
    # 999 particles of 1 µm and 1 of 10 µm: the large one holds half of the volume
    d = np.array([1.0] * 999 + [10.0])
    stats = compute_particle_size_stats(d)
    assert stats["D43"] == pytest.approx((999 * 1 + 10**4) / (999 * 1 + 10**3))
    assert stats["Dv10"] == 1.0 and stats["Dv90"] == 10.0
    assert stats["median"] == 1.0  # number-weighted
    # With number weights, the volume weights are multiplied by them
    weighted = compute_particle_size_stats([1.0, 10.0], weights=[999.0, 1.0])
    assert weighted["D43"] == pytest.approx(stats["D43"])
    assert compute_particle_size_stats([3.0, 3.0])["D43"] == pytest.approx(3.0)


def test_size_distribution_plot(tmp_path: Path):
    wide = np.geomspace(0.5, 80.0, 200)  # more than a factor 10: logarithmic bins
    stats = compute_particle_size_stats(wide, weights=np.linspace(1, 2, 200))
    path = save_particle_size_histogram(wide, str(tmp_path), "S", weights=np.linspace(1, 2, 200), size_stats=stats,
                                        image_format="svg")
    assert path.endswith("S_Par_size_distribution_hist.svg") and Path(path).stat().st_size > 0
    # A single size, unweighted, without statistics
    assert Path(save_particle_size_histogram([2.0, 2.0], str(tmp_path), "T")).exists()



def _band_with_frames(n_frames: int, rng: np.random.Generator):
    """Band of n_frames frames of 1000 µm², each with 5-15 particles of lognormal diameters."""
    band = ParticleSizeBand(d_min_um=0.1, d_max_um=100.0, n_frames=0, scanned_area_um2=0.0)
    particles = []
    for i in range(n_frames):
        diameters = rng.lognormal(1.0, 0.4, int(rng.integers(5, 16)))
        frame_particles = [ParticleInfo(id=len(particles) + k, area_um=float(np.pi * d**2 / 4), eq_diameter_um=float(d),
                                        frame_id=f"A{i}", selection="size_stats") for k, d in enumerate(diameters)]
        particles += frame_particles
        band.add_frame(ParticleStatsFrame(frame_id=f"A{i}", frame_width_um=40.0, pixel_size_um=0.04, area_um2=1000.0,
                                          n_particles=len(frame_particles)), [p.id for p in frame_particles])
    return particles, band


def test_bootstrap_confidence_intervals():
    particles, band = _band_with_frames(10, np.random.default_rng(1))
    cis = bootstrap_confidence_intervals(particles, [band], rng=np.random.default_rng(2))
    assert set(cis) == {"D10", "D50", "D90", "D43", "Dv10", "Dv50", "Dv90"}
    stats = compute_particle_size_stats([p.eq_diameter_um for p in particles])
    assert cis["D50"][0] <= stats["median"] <= cis["D50"][1]
    assert cis == bootstrap_confidence_intervals(particles, [band], rng=np.random.default_rng(2))  # reproducible

    # More frames: narrower intervals
    more_particles, more_frames = _band_with_frames(100, np.random.default_rng(1))
    wide, narrow = cis["D50"], bootstrap_confidence_intervals(more_particles, [more_frames])["D50"]
    assert narrow[1] - narrow[0] < (wide[1] - wide[0]) / 2

    # Without bands, particles are resampled
    assert bootstrap_confidence_intervals(particles)["D50"][0] <= stats["median"]
    assert bootstrap_confidence_intervals(particles[:4]) == {}  # too few to resample
    few_particles, few_frames = _band_with_frames(4, np.random.default_rng(1))
    assert bootstrap_confidence_intervals(few_particles, [few_frames]) == {}



def test_tiny_weights_are_not_taken_as_equal():
    # Band weights are particles per µm² (~1e-8): their relative differences must count
    d = np.array([1.0, 2.0, 4.0])
    assert compute_particle_size_stats(d, weights=[1e-9, 1e-9, 3e-9])["mean"] == pytest.approx((1 + 2 + 12) / 5)
