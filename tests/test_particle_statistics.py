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
    FrameRecord,
    FrameSegmentation,
    LedgerConfigs,
    ParticleInfo,
    ParticleSizeBand,
    ParticleStatsFrame,
    SampleLedger,
    SegmentationConfig,
    SizeBand,
    SizeMeasurement,
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
    build_size_result,
    compute_particle_size_stats,
    edge_correction_factor,
    recalculate_particle_stats,
    save_particle_size_histogram,
    save_particle_statistics,
    select_size_stats_particles,
    weighted_size_observations,
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
    ParticleInfo(id=0, area_um=10.0, selection="size_stats"),  # older ledgers: now particles of frames
    ParticleInfo(id=1, area_um=20.0, selection="auto"),
    ParticleInfo(id=2, area_um=30.0, selection="manual"),
    ParticleInfo(id=3, area_um=40.0, selection="list"),
    ParticleInfo(id=4, area_um=50.0),  # selection not recorded (older ledgers)
    ParticleInfo(id=5, selection="auto"),  # no area
]


def test_particles_selected_for_size_statistics():
    # Besides the particles of frames: the particles measured by EDS with automatic selection
    assert [p.id for p in select_size_stats_particles(PARTICLES)] == [1]
    assert select_size_stats_particles(PARTICLES, include_eds_auto_particles=False) == []


def _area(d_um: float) -> float:
    return float(np.pi * d_um**2 / 4)


def _frame(frame_id: str, diameters_um, area_um2: float = 1000.0, corrections=None) -> ParticleStatsFrame:
    """Working frame with particles of the given diameters."""
    return ParticleStatsFrame(frame_id=frame_id, frame_width_um=40.0, pixel_size_um=0.04, area_um2=area_um2,
                              n_particles=len(diameters_um), particle_areas_um2=[_area(d) for d in diameters_um],
                              particle_edge_corrections=list(corrections or [1.0] * len(diameters_um)))


def _band(d_min_um: float, d_max_um: float, frames) -> ParticleSizeBand:
    band = ParticleSizeBand(d_min_um=d_min_um, d_max_um=d_max_um, n_frames=0, scanned_area_um2=0.0)
    for frame in frames:
        band.add_frame(frame)
    return band


def _ledger(tmp_path: Path) -> SampleLedger:
    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=[], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    return create_sample_ledger(str(tmp_path), "S", configs)


def test_result_stored_in_the_ledger(tmp_path: Path):
    cfg = PowderMeasurementConfig(min_area_par=5.0)
    ledger = _ledger(tmp_path)
    segmentation = ledger.append_segmentation(SegmentationConfig.from_powder_meas_cfg(cfg, 0))
    assert ledger.append_segmentation(SegmentationConfig.from_powder_meas_cfg(cfg, 7)) == segmentation  # reused
    other = ledger.append_segmentation(SegmentationConfig.from_powder_meas_cfg(
        cfg.model_copy(update={"par_brightness_thresh": 150}), 0))
    assert (segmentation, other) == (0, 1) and ledger.segmentations[1].label == "threshold_bright · 150 · 80 px"
    measurement = SizeMeasurement(measurement_id=0, status="completed", n_par_target=100, powder_meas_cfg=cfg,
                                  segmentation_id=segmentation, bands=[SizeBand(d_min_um=1.0, d_max_um=10.0)])
    ledger.size_measurements.append(measurement)
    for i, diameters in enumerate([[2.0, 3.0], [4.0]]):
        frame = ledger.add_frame(FrameRecord(
            frame_id="A0", label="A0", purpose="size_stats", measurement_id=0, frame_width_um=40.0, pixel_size_um=0.04,
            area_um2=1000.0, search_window_um=(1.0, 10.0),
            segmentations=[FrameSegmentation(segmentation_id=segmentation, d_range_um=(1.0, 10.0),
                                             particle_areas_um2=[_area(d) for d in diameters],
                                             particle_edge_corrections=[1.0] * len(diameters))],
        ))
        measurement.bands[0].frame_ids.append(frame.frame_id)
    assert measurement.bands[0].frame_ids == ["A0", "A0#2"]  # frame ids unique in the ledger

    result = recalculate_particle_stats(ledger)
    assert (result.measurement_id, result.segmentation_id, result.n_par_target) == (0, 0, 100)
    assert result.n_par_by_selection == {"size_stats": 3} and result.bands[0].frame_ids == ["A0", "A0#2"]
    assert result.size_stats.median == 3.0 and result.powder_meas_cfg.min_area_par == 5.0
    stored = ledger.append_size_result(result)
    assert stored.result_id == 0 and ledger.size_result == stored
    ledger.to_json_file(str(tmp_path / "ledger.json"))
    reloaded = SampleLedger.from_json_file(str(tmp_path / "ledger.json"))
    assert reloaded.size_result == stored and reloaded.frames == ledger.frames

    # Same inputs: the result is replaced; other exclusions: a new result (the active one), the first one kept
    assert ledger.append_size_result(recalculate_particle_stats(ledger)).result_id == 0
    excluded = ledger.append_size_result(recalculate_particle_stats(ledger, [("A0", 1)], ["A0#2"]))
    assert excluded.result_id == 1 and ledger.active_size_result == 1 and len(ledger.size_results) == 2
    assert excluded.size_stats.n_par_analysed == 1 and excluded.excluded_particles == [("A0", 1)]

    # Results of a segmentation are added to frames, never replaced
    ledger.add_frame_segmentation("A0", FrameSegmentation(segmentation_id=other, d_range_um=(1.0, 10.0)))
    with pytest.raises(ValueError, match="already segmented"):
        ledger.add_frame_segmentation("A0", FrameSegmentation(segmentation_id=other, d_range_um=(1.0, 10.0)))
    # Frames of A0 without result of segmentation 'other' are left out of its statistics
    by_other = recalculate_particle_stats(ledger, [], [], segmentation_id=other)
    assert by_other is None  # A0 has no particle under it, A0#2 is not segmented with it

    data = ledger.model_dump()
    for bad, match in (
        ({"size_measurements": [measurement.model_copy(update={"segmentation_id": 9}).model_dump()]},
         "unknown segmentation"),
        ({"active_size_result": 5}, "active_size_result"),
        ({"frames": data["frames"] + [data["frames"][0]]}, "frame_id values must be unique"),
    ):
        with pytest.raises(ValueError, match=match):
            SampleLedger.model_validate({**data, **bad})



# --- Weighting by size band ----------------------------------------------------------------------------------------

def _two_band_measurement(n_large: int, n_small: int):
    """
    Sample with 1000 particles of 5 µm and 10 particles of 50 µm per mm²: the large band scans the area holding
    n_large large particles, the small band the (100 times smaller) area holding n_small small ones.
    """
    return [
        _band(10.0, 100.0, [_frame("L", [50.0] * n_large, area_um2=n_large / 10 * 1e6)]),
        _band(1.0, 10.0, [_frame("S", [5.0] * n_small, area_um2=n_small / 1000 * 1e6)]),
    ]


def test_band_weights_give_the_sample_distribution_whatever_the_number_collected():
    result = build_size_result(PowderMeasurementConfig(), bands=_two_band_measurement(n_large=10, n_small=100))
    # 1 % of the particles are large: D10-D90 are small particles, the mean is weighted accordingly
    assert result.size_stats.median == pytest.approx(5.0) and result.size_stats.D90 == pytest.approx(5.0)
    assert result.size_stats.mean == pytest.approx((1000 * 5.0 + 10 * 50.0) / 1010)
    assert result.unweighted_size_stats.mean == pytest.approx((100 * 5.0 + 10 * 50.0) / 110)  # 9 % large, as before
    assert [b.number_density_per_mm2 for b in result.bands] == pytest.approx([10.0, 1000.0])
    assert [b.n_particles for b in result.bands] == [10, 100] and result.n_par_by_selection == {"size_stats": 110}

    # Collecting 3 times more large particles (over 3 times the area) does not change the statistics
    other = build_size_result(PowderMeasurementConfig(), bands=_two_band_measurement(n_large=30, n_small=100))
    for key in ("mean", "stdev", "median", "D10", "D25", "D75", "D90"):
        assert getattr(other.size_stats, key) == pytest.approx(getattr(result.size_stats, key))


def test_single_band_statistics_are_unweighted():
    band = _band(1.0, 10.0, [_frame(f"A{i}", [2 * np.sqrt(a / np.pi)], area_um2=2500.0)
                             for i, a in enumerate([3, 8, 20, 50])])
    result = build_size_result(PowderMeasurementConfig(), bands=[band])
    assert result.size_stats == result.unweighted_size_stats
    d = np.sort(np.random.default_rng(0).lognormal(1.0, 0.5, 50))
    for q in (10, 25, 50, 75, 90):  # equal weights: as numpy, method 'hazen'
        assert _weighted_percentile(d, np.full(d.size, 0.3), q) == pytest.approx(np.percentile(d, q, method="hazen"))
    # A heavy last value counts in proportion to its weight
    assert _weighted_percentile(np.array([1.0, 10.0]), np.array([1.0, 9.0]), 60) == 10.0  # its middle is at 55 %
    assert _weighted_percentile(np.array([1.0, 10.0]), np.array([1.0, 9.0]), 50) == pytest.approx(1 + 9 * 0.45 / 0.5)


def test_eds_particles_share_the_weight_of_their_band():
    # Samples without frames searched by EDS: the particles measured by EDS share the weight of their band
    bands = _two_band_measurement(n_large=10, n_small=100)
    eds = [
        ParticleInfo(id=1000, area_um=_area(30.0), selection="auto"),  # 30 µm: large band
        ParticleInfo(id=1001, area_um=_area(200.0), selection="auto"),  # 200 µm: no band
    ]
    obs, weights, weighted_bands, _ = weighted_size_observations(eds, bands)
    assert [b.n_eds_particles for b in weighted_bands] == [1, 0]
    large_band = (obs.eq_diameters_um > 20) & (obs.eq_diameters_um < 100)
    assert weights[large_band].sum() == pytest.approx(10 / 1e6)  # band weight unchanged: its density per µm²
    assert weights[obs.eq_diameters_um > 100].sum() == 0

    result = build_size_result(PowderMeasurementConfig(), bands=bands, particles=eds)
    assert result.n_par_without_weight == 1 and result.size_stats.n_par_analysed == 111
    assert result.size_stats.max == pytest.approx(50.0)  # the 200 µm particle is left out
    assert result.n_par_by_selection == {"size_stats": 110, "auto": 2}


def test_frames_and_particles_left_out():
    band = _band(1.0, 10.0, [_frame("A0", [2.0, 3.0], area_um2=6000.0), _frame("A1", [], area_um2=6000.0)])
    assert (band.n_frames, band.scanned_area_um2, band.n_particles) == (2, 12000.0, 2)
    assert band.frame_area_um2(["A1", "B9"]) == 6000.0
    with pytest.raises(ValueError, match="n_frames"):
        ParticleSizeBand.model_validate({**band.model_dump(), "n_frames": 3})

    # Leaving out the empty frame A1 halves the area: the density doubles
    _, _, [weighted], _ = weighted_size_observations([], [band], excluded_frame_ids=["A1"])
    assert weighted.number_density_per_mm2 == pytest.approx(2 / 6000 * 1e6)
    # A particle left out (frame, index)
    obs, _, [weighted], _ = weighted_size_observations([], [band], excluded_particles=[("A0", 1)])
    assert list(obs.eq_diameters_um) == pytest.approx([2.0]) and obs.particle_keys == [("A0", 0)]


def test_edge_correction():
    # 100 x 50 px frame, 3 px margin: a 1 px particle has 92 x 42 allowed positions out of 100 x 50
    assert edge_correction_factor(1, 1, 100, 50, 3) == pytest.approx(100 * 50 / (92 * 42))
    assert edge_correction_factor(40, 20, 100, 50, 3) > edge_correction_factor(10, 5, 100, 50, 3) > 1
    with pytest.raises(ValueError):
        edge_correction_factor(93, 10, 100, 50, 3)

    # In a band, particles count (and share its weight) in proportion to their edge correction
    band = _band(1.0, 10.0, [_frame("A0", [2.0, 8.0], area_um2=1e3, corrections=[1.0, 3.0])])
    _, weights, [weighted], _ = weighted_size_observations([], [band])
    assert weights[1] / weights[0] == pytest.approx(3.0)
    assert weighted.number_density_per_mm2 == pytest.approx(4 / 1e3 * 1e6)
    result = build_size_result(PowderMeasurementConfig(), bands=[band])
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
    """Band of n_frames frames of 1000 µm², each with 5-15 particles of lognormal diameters, and the diameters."""
    frames = [_frame(f"A{i}", rng.lognormal(1.0, 0.4, int(rng.integers(5, 16)))) for i in range(n_frames)]
    band = _band(0.1, 100.0, frames)
    diameters = [2 * np.sqrt(a / np.pi) for f in frames for a in f.particle_areas_um2]
    return diameters, band


def test_bootstrap_confidence_intervals():
    diameters, band = _band_with_frames(10, np.random.default_rng(1))
    cis = bootstrap_confidence_intervals([], [band], rng=np.random.default_rng(2))
    assert set(cis) == {"D10", "D50", "D90", "D43", "Dv10", "Dv50", "Dv90"}
    stats = compute_particle_size_stats(diameters)
    assert cis["D50"][0] <= stats["median"] <= cis["D50"][1]
    assert cis == bootstrap_confidence_intervals([], [band], rng=np.random.default_rng(2))  # reproducible

    # More frames: narrower intervals
    _, more_frames = _band_with_frames(100, np.random.default_rng(1))
    wide, narrow = cis["D50"], bootstrap_confidence_intervals([], [more_frames])["D50"]
    assert narrow[1] - narrow[0] < (wide[1] - wide[0]) / 2

    # Without bands, particles (measured by EDS) are resampled
    particles = [ParticleInfo(id=i, area_um=_area(d), selection="auto") for i, d in enumerate(diameters)]
    assert bootstrap_confidence_intervals(particles)["D50"][0] <= stats["median"]
    assert bootstrap_confidence_intervals(particles[:4]) == {}  # too few to resample
    _, few_frames = _band_with_frames(4, np.random.default_rng(1))
    assert bootstrap_confidence_intervals([], [few_frames]) == {}



def test_tiny_weights_are_not_taken_as_equal():
    # Band weights are particles per µm² (~1e-8): their relative differences must count
    d = np.array([1.0, 2.0, 4.0])
    assert compute_particle_size_stats(d, weights=[1e-9, 1e-9, 3e-9])["mean"] == pytest.approx((1 + 2 + 12) / 5)


def test_size_distribution_reproduces_the_statistics(tmp_path: Path):
    from _psd_helpers import make_eds_frames_sample, make_measurement_sample

    from autoemx.config.ledger_schemas import SampleLedger
    from autoemx.core.particle_statistics import compute_particle_size_stats, recalculate_particle_stats, size_distribution

    measured = SampleLedger.from_json_file(str(make_measurement_sample(tmp_path / "M") / "ledger.json"))
    eds_only = SampleLedger.from_json_file(str(make_eds_frames_sample(tmp_path / "E") / "ledger.json"))
    cases = [
        (measured, measured.size_result),  # stored result
        (eds_only, recalculate_particle_stats(eds_only)),  # computed from the frames searched by EDS
        (measured, recalculate_particle_stats(measured, frame_ids_to_filter=["B0F0"])),  # frame left out
    ]
    for ledger, result in cases:
        d, w = size_distribution(ledger, result)
        assert np.all(w > 0)
        stats = compute_particle_size_stats(d, w)
        assert stats["n_par_analysed"] == result.size_stats.n_par_analysed
        assert stats["median"] == pytest.approx(result.size_stats.median)
        assert stats["D90"] == pytest.approx(result.size_stats.D90)
    # Default: the statistics saved in the ledger
    assert np.array_equal(size_distribution(measured)[0], size_distribution(measured, measured.size_result)[0])
