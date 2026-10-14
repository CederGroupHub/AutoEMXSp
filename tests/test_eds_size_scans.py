#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Frames searched by EDS counted in particle size statistics (SampleLedger.eds_size_scans)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import cv2
import numpy as np
import pytest

from _image_helpers import load_grayscale_image, make_dev_controller
from ci_paths import INPUTS_DIR

from autoemx.config.ledger_io import create_sample_ledger
from autoemx.config.ledger_schemas import (
    Coordinate2D,
    LedgerConfigs,
    ParticleInfo,
    ParticleSizeBand,
    ParticleStatsFrame,
)
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.composition_analysis.ledger import LedgerModule
from autoemx.core.em_runtime.particle_finder import EM_Particle_Finder, min_valid_diameter_um
from autoemx.core.particle_statistics import (
    overlapping_frame_ids,
    recalculate_particle_stats,
    weighted_size_observations,
)


def _area(d_um: float) -> float:
    return float(np.pi * d_um**2 / 4)


def _eds_frame(frame_id: str, diameters_um, area_um2: float = 2e3) -> ParticleStatsFrame:
    return ParticleStatsFrame(frame_id=frame_id, frame_width_um=50.0, pixel_size_um=0.05, area_um2=area_um2,
                              n_particles=len(diameters_um), particle_areas_um2=[_area(d) for d in diameters_um],
                              particle_edge_corrections=[1.0] * len(diameters_um))


def _eds_scan(frames) -> ParticleSizeBand:
    scan = ParticleSizeBand(d_min_um=3.0, d_max_um=30.0, frame_width_um=50.0, n_frames=0, scanned_area_um2=0.0)
    for frame in frames:
        scan.add_frame(frame, [])
    return scan


def test_frames_searched_by_eds_are_recorded(tmp_path: Path):
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    em, cfg = make_dev_controller(str(tmp_path))
    cfg.min_area_par, cfg.max_area_par, cfg.min_par_pixels = 1.0, 100.0, 10
    finder = EM_Particle_Finder(em, cfg, results_dir=str(tmp_path), development_mode=False)
    em.pixel_size_um, em.frame_navigator.current_frame_label = 1.0, "A3"
    # As in the frame search of EDS: segmentation of the frame, then record
    par_mask, _ = finder._get_particles_on_substrate_mask(cv2.GaussianBlur(image, (5, 5), 0))
    num_labels, _, stats, _ = finder._get_connected_components_with_stats(par_mask)
    finder._record_eds_search_frame(image.shape, num_labels, stats)

    [((d_min, d_max, frame_width), frame)] = finder.pending_eds_frames
    assert d_min == pytest.approx(max(2 * np.sqrt(1.0 / np.pi), min_valid_diameter_um(1.0, 10)))
    assert d_max == pytest.approx(2 * np.sqrt(100.0 / np.pi))
    assert frame_width == image.shape[1] and frame.area_um2 == image.shape[0] * image.shape[1]
    assert frame.frame_id == "EDS-A3"
    # All the valid particles of the frame are counted, not only those EDS then visits
    assert frame.n_particles > cfg.max_n_par_per_frame
    assert all(10 <= a <= 100 for a in frame.particle_areas_um2)
    assert all(c >= 1 for c in frame.particle_edge_corrections)


def test_overlapping_size_ranges_weighted_by_the_area_able_to_observe_them():
    # Size band 1-10 µm over 1000 µm², EDS frames 3-30 µm over 2 x 2000 µm²
    psd = [ParticleInfo(id=i, area_um=_area(d), eq_diameter_um=d, selection="size_stats") for i, d in enumerate([2.0, 5.0])]
    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=1, scanned_area_um2=1e3, particle_ids=[0, 1])
    eds = _eds_scan([_eds_frame("EDS-A0", [5.0, 20.0]), _eds_frame("EDS-A1", [])])
    obs, weights, [weighted_band], [weighted_eds] = weighted_size_observations(psd, [band], (), [eds])
    expected = {2.0: 1 / 1e3, 5.0: 1 / (1e3 + 4e3), 20.0: 1 / 4e3}  # sizes seen by the band, by both, by EDS only
    assert weights == pytest.approx([expected[d] for d in obs.eq_diameters_um])
    assert weighted_band.number_density_per_mm2 == pytest.approx(2 / 1e3 * 1e6)
    assert weighted_eds.number_density_per_mm2 == pytest.approx(2 / 4e3 * 1e6)

    # Leaving out a frame searched by EDS removes its particles and area
    obs, weights, _, _ = weighted_size_observations(psd, [band], ["EDS-A0"], [eds])
    assert list(obs.eq_diameters_um) == [2.0, 5.0]
    assert weights == pytest.approx([1 / 1e3, 1 / (1e3 + 2e3)])


def _ledger(tmp_path: Path):
    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=["Si"], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    return create_sample_ledger(str(tmp_path), "S", configs)


def test_eds_frames_stored_in_the_ledger(tmp_path: Path):
    ledger = _ledger(tmp_path)
    finder = SimpleNamespace(pending_eds_frames=[
        ((3.0, 30.0, 50.0), _eds_frame("EDS-A0", [5.0])),
        ((3.0, 30.0, 50.0), _eds_frame("EDS-A1", [6.0, 7.0])),
        ((1.0, 10.0, 20.0), _eds_frame("EDS-B0", [2.0])),  # other settings: another scan
    ])
    analyzer = SimpleNamespace(EM_controller=SimpleNamespace(particle_finder=finder))
    assert LedgerModule._store_eds_search_frames(analyzer, ledger)
    assert finder.pending_eds_frames == []
    assert [[f.frame_id for f in scan.frames] for scan in ledger.eds_size_scans] == [["EDS-A0", "EDS-A1"], ["EDS-B0"]]
    assert ledger.eds_size_scans[0].scanned_area_um2 == 4e3
    assert not LedgerModule._store_eds_search_frames(analyzer, ledger)  # nothing new


def test_statistics_of_a_sample_measured_only_by_eds(tmp_path: Path):
    ledger = _ledger(tmp_path)
    rng = np.random.default_rng(0)
    ledger.eds_size_scans = [_eds_scan([_eds_frame(f"EDS-A{i}", rng.lognormal(2.0, 0.3, 10)) for i in range(6)])]
    # Particles measured by EDS are counted by the frames: left out
    ledger.particles = [ParticleInfo(id=1, area_um=_area(8.0), selection="auto", frame_id="A0")]

    result = recalculate_particle_stats(ledger)
    assert result.particle_ids == [] and result.n_par_by_selection == {"eds_frames": 60}
    assert result.size_stats.n_par_analysed == 60 and len(result.eds_scans) == 1
    assert "D50" in result.confidence_intervals_95  # 6 frames: enough for confidence intervals
    assert recalculate_particle_stats(ledger, frame_ids_to_filter=["EDS-A0"]).size_stats.n_par_analysed == 50
    # Without what EDS measured: nothing left
    assert recalculate_particle_stats(ledger, include_eds_auto_particles=False) is None



def test_eds_frame_searched_again_counted_once(tmp_path: Path):
    ledger = _ledger(tmp_path)

    def frame(frame_id, x_mm):
        return ((3.0, 30.0, 50.0), _eds_frame(frame_id, [5.0]).model_copy(update={"center_mm": Coordinate2D(x=x_mm, y=1.0)}))

    finder = SimpleNamespace(pending_eds_frames=[frame("EDS-A0", 0.0), frame("EDS-A1", 0.05)])
    analyzer = SimpleNamespace(EM_controller=SimpleNamespace(particle_finder=finder))
    LedgerModule._store_eds_search_frames(analyzer, ledger)
    # Resumed acquisition: A1 searched again (same position, 0.1 µm off), and a new frame
    finder.pending_eds_frames = [frame("EDS-A1", 0.0501), frame("EDS-A2", 0.10)]
    LedgerModule._store_eds_search_frames(analyzer, ledger)
    [scan] = ledger.eds_size_scans
    assert [f.frame_id for f in scan.frames] == ["EDS-A0", "EDS-A1", "EDS-A1#2", "EDS-A2"]  # all recorded, ids unique
    assert overlapping_frame_ids(ledger.eds_size_scans) == ["EDS-A1#2"]  # ... but counted once
    obs, _, _, [weighted] = weighted_size_observations([], [], (), ledger.eds_size_scans)
    assert len(obs.eq_diameters_um) == 3 and weighted.number_density_per_mm2 == pytest.approx(3 / 6e3 * 1e6)



def _frame(frame_id, x_mm, width_um, diameters=()):
    return _eds_frame(frame_id, list(diameters), area_um2=width_um * width_um * 0.625).model_copy(
        update={"center_mm": Coordinate2D(x=x_mm, y=0.0), "frame_width_um": width_um})


def test_overlapping_frames_counted_once_keeping_the_largest():
    # Size band 1-10 µm: frames of 50 µm; EDS frames of 100 µm, 3-30 µm, one overlapping band frame B1
    band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=0, scanned_area_um2=0.0)
    for frame_id, x_mm, ids in (("B0", 0.0, [0]), ("B1", 0.2, [1])):
        band.add_frame(_frame(frame_id, x_mm, 50.0).model_copy(update={"particle_areas_um2": None,
                                                                         "particle_edge_corrections": None,
                                                                         "n_particles": 1}), ids)
    eds = _eds_scan([_frame("EDS-A0", 0.23, 100.0, [5.0]), _frame("EDS-A1", 0.5, 100.0, [6.0])])
    assert overlapping_frame_ids([band, eds]) == ["B1"]  # the smaller frame of the overlapping pair

    psd = [ParticleInfo(id=i, area_um=_area(4.0), eq_diameter_um=4.0, frame_id=f, selection="size_stats")
           for i, f in enumerate(["B0", "B1"])]
    obs, _, [weighted_band], _ = weighted_size_observations(psd, [band], (), [eds])
    assert [obs.frame_ids[k] for k in range(len(obs.frame_ids))] == ["B0", "EDS-A0", "EDS-A1"]
    assert weighted_band.number_density_per_mm2 == pytest.approx(1 / (50 * 50 * 0.625) * 1e6)  # B1 area removed

    # Frames of different size bands (a subframe inside its parent) count different particles: both kept
    finer = ParticleSizeBand(d_min_um=0.1, d_max_um=1.0, n_frames=0, scanned_area_um2=0.0)
    finer.add_frame(_frame("B0_a0", 0.0, 5.0).model_copy(update={"particle_areas_um2": None,
                                                                   "particle_edge_corrections": None,
                                                                   "n_particles": 0}), [])
    assert overlapping_frame_ids([band, finer]) == []


def test_frame_position_only_when_known_and_unique_frame_ids(tmp_path: Path):
    em, cfg = make_dev_controller(str(tmp_path))
    finder = EM_Particle_Finder(em, cfg, results_dir=str(tmp_path), development_mode=False)
    navigator = em.frame_navigator
    navigator._current_pos, navigator.current_pos_is_measured = (1.0, 2.0), True
    assert finder._current_frame_position() == (1.0, 2.0)
    navigator.current_pos_is_measured = False  # stage position not read after a manual move
    assert finder._current_frame_position() is None

    # Manual frames of a resumed measurement restart their numbering: ids made unique in the band
    navigator.current_frame_label = "0"
    finder._stats_band = ParticleSizeBand(d_min_um=1.0, d_max_um=10.0, n_frames=0, scanned_area_um2=0.0)
    assert finder._stats_frame_id() == "0"
    finder._stats_band.add_frame(_frame("0", 0.0, 50.0), [])
    assert finder._stats_frame_id() == "0#2"
