#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Particle size statistics collected end to end, with a fake microscope serving the example frame image."""

from __future__ import annotations

import json
import sys
import time
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import autoemx.core.composition_analysis.analyser as analyser_mod
import autoemx.core.em_runtime.particle_finder as particle_finder_mod
import autoemx.utils.constants as cnst
import test_acquisition_with_quantification as eds_test
from autoemx.config.ledger_io import create_sample_ledger
from autoemx.config.ledger_schemas import LedgerConfigs, ParticleInfo, ParticleStatsFrame, SampleLedger
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.em_runtime.particle_finder import (
    EM_Particle_Finder,
    min_valid_diameter_um,
    plan_size_bands_um,
)
from autoemx.runners.collect_particle_statistics import collect_particle_statistics
from autoemx.runners.process_particle_statistics import process_particle_statistics

from _image_helpers import load_grayscale_image
from ci_paths import INPUTS_DIR

SAMPLE_ID = "PSD_test"
N_FRAMES = 3


class _FakeMicroscopeCtrl:
    _last_EM_adjustment_time = time.time()
    refresh_time = 1e9


class _FakeFrameNavigator:
    """Size bands: N_FRAMES top-level frames for the largest particles, then one subframe per visited parent."""

    def __init__(self):
        self.grid_search_fw_mm = None
        self.current_frame_label = ""
        self._current_pos = None
        self._labels = []
        self._cursor = 0

    def _set_frames(self, labels):
        self._labels, self._cursor = list(labels), 0

    def configure_top_level_particle_frames(self, frame_width_mm):
        self.grid_search_fw_mm = frame_width_mm
        self._set_frames(f"A{i}" for i in range(N_FRAMES))

    def configure_particle_subframes(self, parent_frames, child_fw_mm):
        self.grid_search_fw_mm = child_fw_mm
        self._set_frames(f"{label}_a0" for label, _, _ in parent_frames)

    def reset_frame_cursor(self):
        self._cursor = 0

    def skip_frames(self, frame_labels):
        self._set_frames(label for label in self._labels if label not in set(frame_labels))

    def get_visited_band_frames(self):
        return [(label, (0.0, 0.0), self.grid_search_fw_mm) for label in self._labels[:self._cursor]]

    def go_to_next_frame(self):
        if self._cursor == len(self._labels):
            return False
        self.current_frame_label = self._labels[self._cursor]
        self._current_pos = (float(self._cursor), 0.0)
        self._cursor += 1
        return True


class _FakeEMController:
    """Visits N_FRAMES frames (or the frames of _FakeFrameNavigator, if used), each showing the example frame image."""

    use_frame_navigator = False

    def __init__(self, configs, sample_id, results_dir, **kwargs):
        image, pixel_size_um = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
        self.image = image
        self.sample_id = sample_id
        self.im_height, self.im_width = image.shape
        self.pixel_size_um = pixel_size_um
        self.measurement_cfg = configs.measurement_cfg
        self.is_initialized = True
        self.microscope_ctrl = _FakeMicroscopeCtrl()
        self.frame_navigator = _FakeFrameNavigator() if self.use_frame_navigator else None
        self.current_frame_label = None
        self.n_frames_visited = 0
        self.saved_images = []
        _FakeEMController.last = self
        self.particle_finder = EM_Particle_Finder(self, configs.powder_meas_cfg, results_dir=results_dir, verbose=False,
                                                  development_mode=False)

    def initialise_SEM(self): pass
    def initialise_sample_navigator(self, **kwargs): pass
    def save_frame_image(self, filename, **kwargs):
        self.saved_images.append((filename, kwargs.get("im_annotations"), kwargs.get("scalebar")))

    def adjust_BCF(self): pass

    def go_to_next_frame(self):
        if self.frame_navigator is not None:
            moved = self.frame_navigator.go_to_next_frame()
            self.current_frame_label = self.frame_navigator.current_frame_label
            # The pixel size follows the frame width set by the navigation
            self.pixel_size_um = self.frame_navigator.grid_search_fw_mm * 1000 / self.im_width
            return moved
        if self.n_frames_visited == N_FRAMES:
            return False
        self.current_frame_label = f"A{self.n_frames_visited}"
        self.n_frames_visited += 1
        return True

    def convert_pixel_pos_to_mm(self, pixel_pos):
        return (self.n_frames_visited + pixel_pos[0] * 1e-4, pixel_pos[1] * 1e-4)


def _frame_area_um2() -> float:
    image, pixel_size_um = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    return image.shape[0] * image.shape[1] * pixel_size_um**2


def _d_um(area_um2: float) -> float:
    return float(2 * np.sqrt(area_um2 / np.pi))


@pytest.fixture
def fake_microscope(monkeypatch):
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    monkeypatch.setattr(analyser_mod, "EM_Controller", _FakeEMController)
    monkeypatch.setattr(particle_finder_mod.EM_driver, "is_microscope_connected", lambda: True, raising=False)
    monkeypatch.setattr(particle_finder_mod.EM_driver, "get_image_data", lambda *args: image, raising=False)


def _collect(results_dir: Path, n_par_target: int, min_area_par: float = 1.0, max_area_par: float = 100.0,
             min_par_pixels: int = 10, **powder_meas_cfg_kwargs) -> Path:
    collect_particle_statistics(
        [{"ID": SAMPLE_ID, "pos": (0.0, 0.0)}], n_par_target=n_par_target,
        powder_meas_cfg_kwargs={"min_area_par": min_area_par, "max_area_par": max_area_par,
                                "min_par_pixels": min_par_pixels, **powder_meas_cfg_kwargs},
        verbose=False, results_dir=str(results_dir),
    )
    return results_dir / SAMPLE_ID


def test_particle_statistics_saved_to_csv_and_ledger(tmp_path: Path, fake_microscope):
    sample_dir = _collect(tmp_path, n_par_target=10**6)  # more than available: all frames are visited

    sizes = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv", dtype={cnst.FRAME_ID_DF_KEY: str})
    stats = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}.csv")
    assert list(sample_dir.glob(f"{SAMPLE_ID}_{cnst.PARTICLE_STAT_HIST_FILENAME}*.png"))

    # Every frame shows the same image, so it contributes the same particles
    per_frame = sizes.groupby(cnst.FRAME_ID_DF_KEY).size()
    assert list(per_frame.index) == [f"A{i}" for i in range(N_FRAMES)]
    assert per_frame.nunique() == 1 and per_frame.iloc[0] > 5
    assert list(sizes[cnst.PAR_ID_DF_KEY]) == list(range(len(sizes)))
    assert sizes[cnst.PAR_AREA_UM_KEY].between(1.0, 100.0).all()

    diameters = sizes[cnst.PAR_EQ_D_KEY].to_numpy()
    assert stats["n_par_analysed"].iloc[0] == len(sizes)

    # The particles are also in the ledger
    ledger = json.loads((sample_dir / "ledger.json").read_text())

    # Particles touching the frame edge are left out, the others weighted by their edge correction
    corrections = np.array([p["edge_correction"] for p in ledger["particles"]])
    assert np.all(corrections >= 1) and corrections.max() > corrections.min()
    assert stats["mean"].iloc[0] == pytest.approx(np.sum(corrections * diameters) / np.sum(corrections))
    assert ledger["particle_stats"]["unweighted_size_stats"]["mean"] == pytest.approx(diameters.mean())
    assert ledger["particle_stats"]["unweighted_size_stats"]["median"] == pytest.approx(np.median(diameters))
    assert [p["id"] for p in ledger["particles"]] == list(sizes[cnst.PAR_ID_DF_KEY])
    assert [p["frame_id"] for p in ledger["particles"]] == list(sizes[cnst.FRAME_ID_DF_KEY])
    assert np.allclose([p["eq_diameter_um"] for p in ledger["particles"]], diameters)
    assert ledger["spectra"] == []
    assert {p["selection"] for p in ledger["particles"]} == {"size_stats"}

    # ... with the statistics and the settings of the measurement
    result = ledger["particle_stats"]
    assert result["particle_ids"] == list(sizes[cnst.PAR_ID_DF_KEY])
    assert result["n_par_by_selection"] == {"size_stats": len(sizes)}
    assert result["n_par_target"] == 10**6
    assert result["powder_meas_cfg"]["min_area_par"] == 1.0 and result["powder_meas_cfg"]["include_eds_auto_particles"]
    for key in ("n_par_analysed", "mean", "stdev", "median", "max", "min", "D10", "D25", "D75", "D90"):
        assert result["size_stats"][key] == pytest.approx(stats[key].iloc[0])

    # One size band (1-100 µm² is less than a decade in diameter), no frame width set without automatic navigation
    [band] = result["bands"]
    assert (band["d_min_um"], band["d_max_um"]) == pytest.approx((_d_um(1.0), _d_um(100.0)))
    assert band["frame_width_um"] is None
    assert band["n_frames"] == N_FRAMES
    assert band["scanned_area_um2"] == pytest.approx(N_FRAMES * _frame_area_um2())
    assert band["particle_ids"] == result["particle_ids"]
    # ... with a record of each frame
    assert [f["frame_id"] for f in band["frames"]] == [f"A{i}" for i in range(N_FRAMES)]
    assert sum(f["area_um2"] for f in band["frames"]) == pytest.approx(band["scanned_area_um2"])
    assert [f["n_particles"] for f in band["frames"]] == list(per_frame)
    assert ledger["particle_stats_progress"] is None  # measurement completed

    # Frame images saved without annotations; the particles to mark on them are in the ledger
    saved = _FakeEMController.last.saved_images
    assert [name for name, _, _ in saved] == [f"{SAMPLE_ID}_frA{i}_{kind}" for i in range(N_FRAMES)
                                              for kind in ("particles", "par_mask")]
    assert all(annotations is None and scalebar is False for _, annotations, scalebar in saved)
    for frame in band["frames"]:
        assert frame["image_file"] == f"{SAMPLE_ID}_fr{frame['frame_id']}_particles.png"
        frame_particles = [p for p in ledger["particles"] if p["frame_id"] == frame["frame_id"]]
        assert [mark[3] for mark in frame["particle_marks_px"]] == [p["id"] for p in frame_particles]


def test_collection_stops_at_the_frame_reaching_the_target(tmp_path: Path, fake_microscope):
    sample_dir = _collect(tmp_path, n_par_target=1)
    ledger = json.loads((sample_dir / "ledger.json").read_text())
    assert {p["frame_id"] for p in ledger["particles"]} == {"A0"}


def test_size_bands_recorded(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    im_height, im_width = image.shape
    min_area, max_area, min_pixels = 1e-3, 1e4, 20
    sample_dir = _collect(tmp_path, n_par_target=10**6, min_area_par=min_area, max_area_par=max_area,
                          min_par_pixels=min_pixels, size_stats_plot_format="pdf")

    ledger = json.loads((sample_dir / "ledger.json").read_text())
    particles = {p["id"]: p for p in ledger["particles"]}
    result = ledger["particle_stats"]
    bands = result["bands"]
    expected_bands = plan_size_bands_um(_d_um(min_area), _d_um(max_area), im_width, min_pixels)
    assert len(bands) == len(expected_bands) >= 2
    assert [(b["d_min_um"], b["d_max_um"], b["frame_width_um"]) for b in bands] == pytest.approx(expected_bands)

    for i, (band, (d_lo, d_hi, fw_um)) in enumerate(zip(bands, expected_bands)):
        pixel_size_um = fw_um / im_width
        assert band["n_frames"] == N_FRAMES  # top-level frames, then one subframe in each parent
        assert band["scanned_area_um2"] == pytest.approx(N_FRAMES * im_height * im_width * pixel_size_um**2)
        assert band["particle_ids"], f"no particles in band {i}"
        for particle_id in band["particle_ids"]:
            particle = particles[particle_id]
            assert d_lo - 1e-9 <= particle["eq_diameter_um"] <= d_hi + 1e-9
            assert particle["area_um"] / pixel_size_um**2 >= min_pixels - 1e-6  # enough pixels for a valid area
            assert ("_" in particle["frame_id"]) == (i > 0)  # finer bands: subframes
        if band["d_min_um"] > _d_um(min_area):
            assert band["d_min_um"] == pytest.approx(min_valid_diameter_um(pixel_size_um, min_pixels))
    assert [pid for band in bands for pid in band["particle_ids"]] == result["particle_ids"] == list(particles)

    # Statistics weighted by band: each band by its particles per unit area scanned
    for band in bands:  # particles corrected for those touching the frame edge, per area scanned
        n_corrected = sum(particles[pid]["edge_correction"] for pid in band["particle_ids"])
        assert band["number_density_per_mm2"] == pytest.approx(n_corrected / band["scanned_area_um2"] * 1e6)
    stats = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}.csv")
    sizes = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv")
    assert stats["n_bands"].iloc[0] == len(bands)
    for key in ("D43", "Dv10", "Dv50", "Dv90"):  # volume-weighted values
        assert stats[key].iloc[0] == pytest.approx(result["size_stats"][key])
    assert result["size_stats"]["Dv50"] > result["size_stats"]["median"]
    assert (sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STAT_HIST_FILENAME}.pdf").exists()
    assert stats["mean"].iloc[0] == pytest.approx(result["size_stats"]["mean"])
    assert result["size_stats"]["mean"] != pytest.approx(result["unweighted_size_stats"]["mean"])
    weights = sizes[cnst.PAR_WEIGHT_KEY].to_numpy()
    assert weights.mean() == pytest.approx(1.0)
    assert np.sum(weights * sizes[cnst.PAR_EQ_D_KEY]) / np.sum(weights) == pytest.approx(result["size_stats"]["mean"])


def test_measured_statistics_recalculated_from_the_ledger(tmp_path: Path, fake_microscope):
    sample_dir = _collect(tmp_path, n_par_target=10**6)
    (sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv").unlink()  # read from the ledger

    stats = process_particle_statistics(SAMPLE_ID, str(tmp_path), frame_ids_to_filter=["A1"], verbose=False)
    ledger = json.loads((sample_dir / "ledger.json").read_text())
    kept = [p for p in ledger["particles"] if p["frame_id"] != "A1"]
    result = ledger["particle_stats"]
    assert stats["n_par_analysed"].iloc[0] == len(kept) == len(ledger["particles"]) * 2 // 3
    assert result["particle_ids"] == [p["id"] for p in kept]
    assert result["bands"][0]["particle_ids"] == [p["id"] for p in ledger["particles"]]
    assert result["excluded_frame_ids"] == ["A1"]
    band = result["bands"][0]
    frame_area = next(f["area_um2"] for f in band["frames"] if f["frame_id"] == "A1")
    corrections = np.array([p["edge_correction"] for p in kept])
    assert band["number_density_per_mm2"] == pytest.approx(
        corrections.sum() / (band["scanned_area_um2"] - frame_area) * 1e6)
    diameters = np.array([p["eq_diameter_um"] for p in kept])
    assert result["size_stats"]["mean"] == pytest.approx(np.sum(corrections * diameters) / corrections.sum())


# --- Particle size statistics and EDS on the same sample ---------------------------------------------------------

def _eds_ledger(sample_dir: Path) -> None:
    """Ledger of a sample measured by EDS: particles 1-3 (auto, auto, manual)."""
    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"), sample_cfg=SampleConfig(elements=["Si", "O"], type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    ledger = create_sample_ledger(str(sample_dir), SAMPLE_ID, configs)
    ledger.particles = [
        ParticleInfo(id=1, area_um=30.0, frame_id="B0", selection="auto"),
        ParticleInfo(id=2, area_um=40.0, frame_id="B0", selection="auto"),
        ParticleInfo(id=3, area_um=50.0, frame_id="B1", selection="manual"),
    ]
    sample_dir.mkdir(parents=True)
    ledger.to_json_file(str(sample_dir / "ledger.json"))


def test_size_statistics_after_eds(tmp_path: Path, fake_microscope):
    _eds_ledger(tmp_path / SAMPLE_ID)
    sample_dir = _collect(tmp_path, n_par_target=1)  # the 3 EDS particles do not count towards the target

    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    eds, psd = ledger.particles[:3], ledger.particles[3:]
    assert [(p.id, p.selection, p.area_um) for p in eds] == [(1, "auto", 30.0), (2, "auto", 40.0), (3, "manual", 50.0)]
    assert [p.id for p in psd] == list(range(4, 4 + len(psd)))
    assert psd and {p.frame_id for p in psd} == {"A0"} and {p.selection for p in psd} == {"size_stats"}
    assert ledger.configs.sample_cfg.elements == ["Si", "O"]  # configs of the EDS measurement kept
    assert ledger.configs.measurement_cfg.type == "EDS"

    result = ledger.particle_stats
    assert result.particle_ids == [1, 2] + [p.id for p in psd]  # EDS particles chosen automatically included
    assert result.n_par_by_selection == {"auto": 2, "size_stats": len(psd)}
    assert result.bands[0].particle_ids == [p.id for p in psd]

    # A new measurement keeps the particles of the earlier one, but leaves them out of its statistics
    _collect(tmp_path, n_par_target=1)
    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    new_psd = ledger.particles[3 + len(psd):]
    assert [p.id for p in ledger.particles[:3 + len(psd)]] == [1, 2, 3] + [p.id for p in psd]
    assert [p.id for p in new_psd] == list(range(4 + len(psd), 4 + 2 * len(psd)))
    assert ledger.particle_stats.particle_ids == [1, 2] + [p.id for p in new_psd]
    # ... also when recalculated
    process_particle_statistics(SAMPLE_ID, str(tmp_path), verbose=False)
    assert SampleLedger.from_json_file(str(sample_dir / "ledger.json")).particle_stats.particle_ids == (
        [1, 2] + [p.id for p in new_psd]
    )


def _run_eds(tmp_path: Path, monkeypatch, eds_frames=()) -> None:
    """EDS acquisition (5 spectra, 3 particles of 20 µm²) on the sample, the frames searched by EDS being eds_frames."""

    class _FakeEDSController(eds_test._FakeEMController):
        """Particle counters start at 1, as those of the particle finder; particles of 20 µm²."""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.particle_finder = SimpleNamespace(current_par_area_um2=20.0, current_par_coordinates_mm=None,
                                                   pending_eds_frames=list(eds_frames))

        def get_XSp_coords(self, next_spectrum_id):
            success, spots, particle_cntr = super().get_XSp_coords(next_spectrum_id)
            return success, spots, None if particle_cntr is None else particle_cntr + 1

    monkeypatch.setattr(analyser_mod, "EM_Controller", _FakeEDSController)
    monkeypatch.setattr(
        analyser_mod.EMXSp_Composition_Analyzer, "_evaluate_clustering_convergence",
        lambda self, tot_n_spectra, n_spectra_to_collect: (tot_n_spectra, False, False),
    )
    analyzer = analyser_mod.EMXSp_Composition_Analyzer(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"),
        sample_cfg=SampleConfig(elements=["Si", "O"], type="powder"),
        measurement_cfg=MeasurementConfig(min_n_spectra=4, max_n_spectra=5,
                                          target_acquisition_counts=eds_test.TARGET_COUNTS),
        sample_substrate_cfg=SampleSubstrateConfig(elements=["C"], auto_detection=False),
        powder_meas_cfg=PowderMeasurementConfig(max_spectra_per_par=2),
        is_acquisition=True, verbose=False, results_dir=str(tmp_path), sample_id=SAMPLE_ID,
    )
    analyzer.run_collection_and_quantification(quantify=False)


def test_eds_after_size_statistics(tmp_path: Path, fake_microscope, monkeypatch):
    sample_dir = _collect(tmp_path, n_par_target=1)
    psd_ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    max_psd_id = max(p.id for p in psd_ledger.particles)
    assert psd_ledger.configs.measurement_cfg.type == MeasurementConfig.PARTICLE_STATS_MEAS_TYPE_KEY

    _run_eds(tmp_path, monkeypatch)

    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    # EDS particles get the ids after the particle size statistics ones, which are kept with their statistics
    eds_ids = [max_psd_id + 1, max_psd_id + 1, max_psd_id + 2, max_psd_id + 3, max_psd_id + 3]
    assert [s.acquisition_details.particle_id for s in ledger.spectra] == eds_ids
    assert ledger.particles[:len(psd_ledger.particles)] == psd_ledger.particles
    assert [(p.id, p.selection) for p in ledger.particles[len(psd_ledger.particles):]] == [
        (pid, "auto") for pid in sorted(set(eds_ids))
    ]
    # The statistics are recalculated to include the EDS particles, with the settings and bands of the measurement
    result, psd_result = ledger.particle_stats, psd_ledger.particle_stats
    assert result.particle_ids == psd_result.particle_ids + sorted(set(eds_ids))
    assert result.n_par_by_selection == {**psd_result.n_par_by_selection, "auto": 3}
    assert result.powder_meas_cfg == psd_result.powder_meas_cfg
    [band], [psd_band] = result.bands, psd_result.bands
    assert band.model_copy(update={"n_eds_particles": 0}) == psd_band  # same band and density
    assert band.n_eds_particles == 3  # the EDS particles (5 µm) share its weight
    # The configs now describe the spectra (read by quantification); the size statistics settings stay in the result
    assert ledger.configs.measurement_cfg.type == "EDS"
    assert ledger.configs.sample_cfg.elements == ["Si", "O"]



# --- Saving progress and resuming ------------------------------------------------------------------------------------

TWO_BANDS = dict(min_area_par=1e-3, max_area_par=1e4, min_par_pixels=20)


def _failing_after(n_images: int, monkeypatch):
    """Microscope images that fail after n_images (e.g. a crash or a lost connection)."""
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    n_calls = {"n": 0}

    def get_image_data(*args):
        n_calls["n"] += 1
        if n_calls["n"] > n_images:
            raise RuntimeError("Lost connection to the microscope")
        return image

    monkeypatch.setattr(particle_finder_mod.EM_driver, "get_image_data", get_image_data, raising=False)


def _bands(sample_dir: Path) -> list:
    return SampleLedger.from_json_file(str(sample_dir / "ledger.json")).particle_stats.bands


def test_interrupted_measurement_saved_and_resumed(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    monkeypatch.setattr(particle_finder_mod, "STATS_SAVE_EVERY_N_FRAMES", 1)
    complete = _bands(_collect(tmp_path / "complete", n_par_target=10**6, **TWO_BANDS))

    # Interrupted at the 2nd frame of the 2nd band: the frames analysed so far are saved
    _failing_after(N_FRAMES + 1, monkeypatch)
    sample_dir = _collect(tmp_path / "resumed", n_par_target=10**6, **TWO_BANDS)
    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    assert ledger.particle_stats is None
    progress = ledger.particle_stats_progress
    assert [band.n_frames for band in progress.bands] == [N_FRAMES, 1]
    assert {p.id for p in ledger.particles} == {pid for band in progress.bands for pid in band.particle_ids}

    # Resumed: the interrupted band continues in the frames not analysed yet
    _failing_after(10**6, monkeypatch)
    _collect(tmp_path / "resumed", n_par_target=10**6, **TWO_BANDS)
    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    assert ledger.particle_stats_progress is None
    resumed = ledger.particle_stats.bands
    for band, complete_band in zip(resumed, complete):
        assert sorted(f.frame_id for f in band.frames) == sorted(f.frame_id for f in complete_band.frames)
        assert len(band.particle_ids) == len(complete_band.particle_ids)
        assert band.scanned_area_um2 == pytest.approx(complete_band.scanned_area_um2)
    assert len(ledger.particles) == sum(len(band.particle_ids) for band in complete)  # no particle counted twice


def test_interrupted_measurement_discarded_with_other_settings(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    monkeypatch.setattr(particle_finder_mod, "STATS_SAVE_EVERY_N_FRAMES", 1)
    _failing_after(N_FRAMES + 1, monkeypatch)
    sample_dir = _collect(tmp_path, n_par_target=10**6, **TWO_BANDS)
    interrupted_ids = {p.id for p in SampleLedger.from_json_file(str(sample_dir / "ledger.json")).particles}

    _failing_after(10**6, monkeypatch)
    _collect(tmp_path, n_par_target=10**6, **{**TWO_BANDS, "min_par_pixels": 30})
    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    assert ledger.particle_stats_progress is None
    assert interrupted_ids <= {p.id for p in ledger.particles}  # kept in the ledger ...
    assert not interrupted_ids & set(ledger.particle_stats.particle_ids)  # ... but left out of the new statistics
    assert ledger.particle_stats.powder_meas_cfg.min_par_pixels == 30



# --- Stopping on precision -----------------------------------------------------------------------------------------

def _varying_images(monkeypatch, seed: int = 0):
    """Frames that differ: the example image with a random vertical stripe (40 % of its width) blanked."""
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    rng = np.random.default_rng(seed)

    def get_image_data(*args):
        frame = image.copy()
        start = int(rng.integers(0, int(image.shape[1] * 0.6)))
        frame[:, start:start + int(image.shape[1] * 0.4)] = 0
        print("VARYING", start, flush=True)
        return frame

    monkeypatch.setattr(particle_finder_mod.EM_driver, "get_image_data", get_image_data, raising=False)


def test_measurement_stops_when_precise(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(sys.modules[__name__], "N_FRAMES", 60)
    _varying_images(monkeypatch)
    # D43 (the fake frames repeat the same particles, so that percentiles barely vary between resamples)
    sample_dir = _collect(tmp_path, n_par_target=5, size_stats_ci_tolerance=0.02, size_stats_ci_quantities=["D43"],
                          size_stats_max_n_frames=60)
    result = SampleLedger.from_json_file(str(sample_dir / "ledger.json")).particle_stats
    [band] = result.bands
    assert result.precision_reached is True
    assert 5 < band.n_frames < 60  # more than the minimum for confidence intervals, fewer than available
    low, high = result.confidence_intervals_95["D43"]
    assert low <= result.size_stats.D43 <= high
    stats = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}.csv")
    assert (stats["D43_CI95_low"].iloc[0], stats["D43_CI95_high"].iloc[0]) == pytest.approx((low, high))


def test_measurement_stopping_on_precision_limited_by_frames(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    monkeypatch.setattr(sys.modules[__name__], "N_FRAMES", 20)
    _varying_images(monkeypatch)
    sample_dir = _collect(tmp_path, n_par_target=1, size_stats_ci_tolerance=0.001, size_stats_max_n_frames=8,
                          **TWO_BANDS)
    result = SampleLedger.from_json_file(str(sample_dir / "ledger.json")).particle_stats
    assert result.precision_reached is False
    assert sum(band.n_frames for band in result.bands) == 8
    # Both bands were extended; the finer one in subframes of the frames of the coarser one
    assert all(band.n_frames > 1 for band in result.bands)
    parents = {frame.frame_id for frame in result.bands[0].frames}
    assert all(frame.frame_id.rsplit("_", 1)[0] in parents for frame in result.bands[1].frames)



def test_frames_searched_by_eds_counted_in_the_statistics(tmp_path: Path, fake_microscope, monkeypatch):
    sample_dir = _collect(tmp_path, n_par_target=1)
    psd_result = SampleLedger.from_json_file(str(sample_dir / "ledger.json")).particle_stats
    eds_frames = [
        ((3.0, 20.0, 200.0), ParticleStatsFrame(
            frame_id=f"EDS-B{i}", frame_width_um=200.0, pixel_size_um=0.1, area_um2=2.4e4, n_particles=4,
            particle_areas_um2=[10.0, 20.0, 40.0, 80.0], particle_edge_corrections=[1.0, 1.0, 1.1, 1.2],
        ))
        for i in range(6)
    ]
    _run_eds(tmp_path, monkeypatch, eds_frames)

    ledger = SampleLedger.from_json_file(str(sample_dir / "ledger.json"))
    [scan] = ledger.eds_size_scans  # saved with the spectra
    assert [f.frame_id for f in scan.frames] == [f"EDS-B{i}" for i in range(6)]
    result = ledger.particle_stats
    # The EDS frames count the particles measured by EDS, which are left out
    assert result.n_par_by_selection == {"size_stats": len(psd_result.particle_ids), "eds_frames": 24}
    assert result.particle_ids == psd_result.particle_ids
    assert len(result.eds_scans) == 1 and result.eds_scans[0].number_density_per_mm2 > 0
    assert result.size_stats.n_par_analysed == len(psd_result.particle_ids) + 24
    assert any(p.selection == "auto" for p in ledger.particles)  # still in the ledger
