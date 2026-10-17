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
from autoemx.config.ledger_schemas import (
    FrameSegmentation,
    LedgerConfigs,
    ParticleInfo,
    ParticleSizeBand,
    ParticleStatsFrame,
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
from autoemx.core import particle_statistics
from autoemx.core.composition_analysis.ledger import LedgerModule
from autoemx.core.em_runtime import image_utilities
from autoemx.core.particle_statistics import frames_overlap
from autoemx.core.em_runtime.particle_finder import (
    EM_Particle_Finder,
    min_valid_diameter_um,
    plan_size_bands_um,
)
from autoemx.runners.collect_particle_statistics import collect_particle_statistics
from autoemx.runners.process_particle_statistics import process_particle_statistics
from autoemx.runners.reanalyze_particle_statistics import reanalyze_particle_statistics

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
    write_images = False  # save the frame images (as during acquisition), e.g. to reanalyse them

    def __init__(self, configs, sample_id, results_dir, **kwargs):
        image, pixel_size_um = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
        self.image = image
        self.sample_id = sample_id
        self.im_height, self.im_width = image.shape
        self.pixel_size_um = pixel_size_um
        self.measurement_cfg = configs.measurement_cfg
        self.microscope_cfg = configs.microscope_cfg
        self.results_dir = results_dir
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
        if self.write_images and kwargs.get("frame_image") is not None:
            image_utilities.save_frame_image(
                kwargs["frame_image"], self.pixel_size_um, self.im_width, self.im_height, self.sample_id,
                self.microscope_cfg, filename, self.results_dir, im_annotations=kwargs.get("im_annotations"),
                scalebar=bool(kwargs.get("scalebar")), image_extension=self.measurement_cfg.saved_images_extension,
            )

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


def _ledger(sample_dir: Path) -> SampleLedger:
    return SampleLedger.from_json_file(str(sample_dir / "ledger.json"))


def _frame_particles(ledger: SampleLedger, segmentation_id: int = 0):
    """Diameters, corrections and frame ids of the particles of the size frames under a segmentation."""
    rows = [(2 * np.sqrt(a / np.pi), c, f.frame_id)
            for f in ledger.frames_of("size_stats") if f.segmentation(segmentation_id)
            for a, c in zip(f.segmentation(segmentation_id).particle_areas_um2,
                            f.segmentation(segmentation_id).particle_edge_corrections)]
    return (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]), [r[2] for r in rows])


def test_particle_statistics_saved_to_csv_and_ledger(tmp_path: Path, fake_microscope):
    sample_dir = _collect(tmp_path, n_par_target=10**6)  # more than available: all frames are visited

    sizes = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv", dtype={cnst.FRAME_ID_DF_KEY: str})
    stats = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}.csv")
    assert list(sample_dir.glob(f"{SAMPLE_ID}_{cnst.PARTICLE_STAT_HIST_FILENAME}*.png"))

    # Every frame shows the same image, so it contributes the same particles (numbered in their frame)
    per_frame = sizes.groupby(cnst.FRAME_ID_DF_KEY).size()
    assert list(per_frame.index) == [f"A{i}" for i in range(N_FRAMES)]
    assert per_frame.nunique() == 1 and per_frame.iloc[0] > 5
    assert list(sizes[cnst.PAR_ID_DF_KEY]) == list(sizes.groupby(cnst.FRAME_ID_DF_KEY).cumcount())
    assert sizes[cnst.PAR_AREA_UM_KEY].between(1.0, 100.0).all()

    diameters = sizes[cnst.PAR_EQ_D_KEY].to_numpy()
    assert stats["n_par_analysed"].iloc[0] == len(sizes)

    # The frames, with their particles, are in the ledger (no ParticleInfo for them)
    ledger = _ledger(sample_dir)
    assert ledger.particles == [] and ledger.spectra == []
    [segmentation] = ledger.segmentations
    assert (segmentation.par_segmentation_model, segmentation.min_par_pixels) == ("threshold_bright", 10)
    assert [(f.frame_id, f.purpose, f.measurement_id) for f in ledger.frames] == [
        (f"A{i}", "size_stats", 0) for i in range(N_FRAMES)]
    frame_d, corrections, frame_ids = _frame_particles(ledger)
    assert frame_d == pytest.approx(diameters) and frame_ids == list(sizes[cnst.FRAME_ID_DF_KEY])

    # Particles touching the frame edge are left out, the others weighted by their edge correction
    assert np.all(corrections >= 1) and corrections.max() > corrections.min()
    assert stats["mean"].iloc[0] == pytest.approx(np.sum(corrections * diameters) / np.sum(corrections))

    # ... with the measurement and its statistics
    [measurement] = ledger.size_measurements
    assert (measurement.status, measurement.n_par_target, measurement.segmentation_id) == ("completed", 10**6, 0)
    assert measurement.powder_meas_cfg.min_area_par == 1.0
    result = ledger.size_result
    assert (result.measurement_id, result.segmentation_id, result.eds_segmentation_id) == (0, 0, None)
    assert result.n_par_by_selection == {"size_stats": len(sizes)}
    assert result.n_par_target == 10**6 and result.powder_meas_cfg.include_eds_auto_particles
    assert result.unweighted_size_stats.mean == pytest.approx(diameters.mean())
    assert result.unweighted_size_stats.median == pytest.approx(np.median(diameters))
    for key in ("n_par_analysed", "mean", "stdev", "median", "max", "min", "D10", "D25", "D75", "D90"):
        assert getattr(result.size_stats, key) == pytest.approx(stats[key].iloc[0])

    # One size band (1-100 µm² is less than a decade in diameter), no frame width set without automatic navigation
    [band] = measurement.bands
    assert (band.d_min_um, band.d_max_um) == pytest.approx((_d_um(1.0), _d_um(100.0)))
    assert band.frame_width_um is None and band.frame_ids == [f"A{i}" for i in range(N_FRAMES)]
    [summary] = result.bands
    assert summary.n_frames == N_FRAMES and summary.n_particles == len(sizes)
    assert summary.scanned_area_um2 == pytest.approx(N_FRAMES * _frame_area_um2())

    # Frame images saved without annotations; the particles to mark on them are in the ledger
    saved = _FakeEMController.last.saved_images
    assert [name for name, _, _ in saved] == [f"{SAMPLE_ID}_frA{i}_{kind}" for i in range(N_FRAMES)
                                              for kind in ("particles", "par_mask")]
    assert all(annotations is None and scalebar is False for _, annotations, scalebar in saved)
    for frame in ledger.frames:
        assert frame.image_file == f"{SAMPLE_ID}_fr{frame.frame_id}_particles.png"
        result_0 = frame.segmentation(0)
        assert len(result_0.particle_marks_px) == result_0.n_particles == per_frame[frame.frame_id]


def test_collection_stops_at_the_frame_reaching_the_target(tmp_path: Path, fake_microscope):
    sample_dir = _collect(tmp_path, n_par_target=1)
    assert [f.frame_id for f in _ledger(sample_dir).frames] == ["A0"]


def test_size_bands_recorded(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    im_height, im_width = image.shape
    min_area, max_area, min_pixels = 1e-3, 1e4, 20
    sample_dir = _collect(tmp_path, n_par_target=10**6, min_area_par=min_area, max_area_par=max_area,
                          min_par_pixels=min_pixels, size_stats_plot_format="pdf")

    ledger = _ledger(sample_dir)
    [measurement] = ledger.size_measurements
    result = ledger.size_result
    bands = measurement.bands
    expected_bands = plan_size_bands_um(_d_um(min_area), _d_um(max_area), im_width, min_pixels)
    assert len(bands) == len(expected_bands) >= 2
    assert [(b.d_min_um, b.d_max_um, b.frame_width_um) for b in bands] == pytest.approx(expected_bands)

    for i, (band, summary, (d_lo, d_hi, fw_um)) in enumerate(zip(bands, result.bands, expected_bands)):
        pixel_size_um = fw_um / im_width
        assert len(band.frame_ids) == summary.n_frames == N_FRAMES  # top-level frames, then a subframe in each parent
        assert summary.scanned_area_um2 == pytest.approx(N_FRAMES * im_height * im_width * pixel_size_um**2)
        assert summary.n_particles, f"no particles in band {i}"
        for frame_id in band.frame_ids:
            frame = ledger.frame(frame_id)
            assert ("_" in frame.label) == (i > 0)  # finer bands: subframes
            assert frame.search_window_um == pytest.approx((band.d_min_um, band.d_max_um))
            for area in frame.segmentation(0).particle_areas_um2:
                assert d_lo - 1e-9 <= 2 * np.sqrt(area / np.pi) <= d_hi + 1e-9
                assert area / pixel_size_um**2 >= min_pixels - 1e-6  # enough pixels for a valid area
        assert min_valid_diameter_um(pixel_size_um, min_pixels) <= band.d_min_um * (1 + 1e-9)

    # Statistics weighted by band: each band by its particles per unit area scanned
    for band, summary in zip(bands, result.bands):  # particles corrected for those touching the frame edge
        corrected = sum(sum(ledger.frame(f).segmentation(0).particle_edge_corrections) for f in band.frame_ids)
        assert summary.number_density_per_mm2 == pytest.approx(corrected / summary.scanned_area_um2 * 1e6)
    stats = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}.csv")
    sizes = pd.read_csv(sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv")
    assert stats["n_bands"].iloc[0] == len(bands)
    for key in ("D43", "Dv10", "Dv50", "Dv90"):  # volume-weighted values
        assert stats[key].iloc[0] == pytest.approx(getattr(result.size_stats, key))
    assert result.size_stats.Dv50 > result.size_stats.median
    assert (sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STAT_HIST_FILENAME}.pdf").exists()
    assert stats["mean"].iloc[0] == pytest.approx(result.size_stats.mean)
    assert result.size_stats.mean != pytest.approx(result.unweighted_size_stats.mean)
    weights = sizes[cnst.PAR_WEIGHT_KEY].to_numpy()
    assert weights.mean() == pytest.approx(1.0)
    assert np.sum(weights * sizes[cnst.PAR_EQ_D_KEY]) / np.sum(weights) == pytest.approx(result.size_stats.mean)


def test_measured_statistics_recalculated_from_the_ledger(tmp_path: Path, fake_microscope):
    sample_dir = _collect(tmp_path, n_par_target=10**6)
    (sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv").unlink()  # read from the ledger

    stats = process_particle_statistics(SAMPLE_ID, str(tmp_path), frame_ids_to_filter=["A1"], verbose=False)
    ledger = _ledger(sample_dir)
    diameters, corrections, frame_ids = _frame_particles(ledger)
    kept = np.array(frame_ids) != "A1"
    result = ledger.size_result
    assert stats["n_par_analysed"].iloc[0] == kept.sum() == len(diameters) * 2 // 3
    assert result.excluded_frame_ids == ["A1"] and len(ledger.size_results) == 2  # the first result kept
    [band] = result.bands
    assert band.frame_ids == ["A0", "A1", "A2"]
    frame_area = ledger.frame("A1").area_um2
    assert band.number_density_per_mm2 == pytest.approx(
        corrections[kept].sum() / (band.scanned_area_um2 - frame_area) * 1e6)
    assert result.size_stats.mean == pytest.approx(
        np.sum(corrections[kept] * diameters[kept]) / corrections[kept].sum())


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

    ledger = _ledger(sample_dir)
    assert [(p.id, p.selection, p.area_um) for p in ledger.particles] == [
        (1, "auto", 30.0), (2, "auto", 40.0), (3, "manual", 50.0)]  # particles measured by EDS kept
    assert [f.frame_id for f in ledger.frames] == ["A0"]
    assert ledger.configs.sample_cfg.elements == ["Si", "O"]  # configs of the EDS measurement kept
    assert ledger.configs.measurement_cfg.type == "EDS"

    result = ledger.size_result
    n_frame = ledger.frames[0].segmentation(0).n_particles
    # Without frames searched by EDS (older ledgers), the EDS particles chosen automatically are included
    assert result.n_par_by_selection == {"auto": 2, "size_stats": n_frame}

    # A new measurement keeps the frames of the earlier one, but leaves them out of its statistics
    _collect(tmp_path, n_par_target=1)
    ledger = _ledger(sample_dir)
    assert [f.frame_id for f in ledger.frames] == ["A0", "A0#2"]
    assert [m.measurement_id for m in ledger.size_measurements] == [0, 1]
    assert ledger.size_result.measurement_id == 1 and ledger.size_result.bands[0].frame_ids == ["A0#2"]
    assert len(ledger.size_results) == 2  # the statistics of the first measurement kept
    # ... also when recalculated
    process_particle_statistics(SAMPLE_ID, str(tmp_path), verbose=False)
    assert _ledger(sample_dir).size_result.bands[0].frame_ids == ["A0#2"]


def _run_eds(tmp_path: Path, monkeypatch, eds_frames=()) -> "analyser_mod.EMXSp_Composition_Analyzer":
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
        powder_meas_cfg=PowderMeasurementConfig(max_spectra_per_par=2, min_par_pixels=10),
        is_acquisition=True, verbose=False, results_dir=str(tmp_path), sample_id=SAMPLE_ID,
    )
    analyzer.run_collection_and_quantification(quantify=False)
    return analyzer


def test_eds_after_size_statistics(tmp_path: Path, fake_microscope, monkeypatch):
    sample_dir = _collect(tmp_path, n_par_target=1)
    psd_ledger = _ledger(sample_dir)
    assert psd_ledger.configs.measurement_cfg.type == MeasurementConfig.PARTICLE_STATS_MEAS_TYPE_KEY

    _run_eds(tmp_path, monkeypatch)

    ledger = _ledger(sample_dir)
    eds_ids = [s.acquisition_details.particle_id for s in ledger.spectra]
    assert [(p.id, p.selection) for p in ledger.particles] == [(pid, "auto") for pid in sorted(set(eds_ids))]
    assert ledger.frames == psd_ledger.frames  # the frames of the measurement are kept
    # The statistics are recalculated to include the EDS particles (no frames searched by EDS recorded here), with the
    # settings and bands of the measurement: same inputs, the result is replaced
    result, psd_result = ledger.size_result, psd_ledger.size_result
    assert len(ledger.size_results) == 1 and result.result_id == psd_result.result_id
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
    """Bands of the active statistics result."""
    return _ledger(sample_dir).size_result.bands


def test_interrupted_measurement_saved_and_resumed(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    monkeypatch.setattr(particle_finder_mod, "STATS_SAVE_EVERY_N_FRAMES", 1)
    complete = _bands(_collect(tmp_path / "complete", n_par_target=10**6, **TWO_BANDS))

    # Interrupted at the 2nd frame of the 2nd band: the frames analysed so far are saved
    _failing_after(N_FRAMES + 1, monkeypatch)
    sample_dir = _collect(tmp_path / "resumed", n_par_target=10**6, **TWO_BANDS)
    ledger = _ledger(sample_dir)
    assert ledger.size_result is None
    [measurement] = ledger.size_measurements
    assert measurement.status == "in_progress"
    assert [len(band.frame_ids) for band in measurement.bands] == [N_FRAMES, 1]
    assert len(ledger.frames) == N_FRAMES + 1

    # Resumed: the interrupted band continues in the frames not analysed yet
    _failing_after(10**6, monkeypatch)
    _collect(tmp_path / "resumed", n_par_target=10**6, **TWO_BANDS)
    ledger = _ledger(sample_dir)
    [measurement] = ledger.size_measurements
    assert measurement.status == "completed"
    resumed = ledger.size_result.bands
    for band, complete_band in zip(resumed, complete):
        assert sorted(band.frame_ids) == sorted(complete_band.frame_ids)
        assert band.n_particles == complete_band.n_particles
        assert band.scanned_area_um2 == pytest.approx(complete_band.scanned_area_um2)
    assert len(ledger.frames) == sum(b.n_frames for b in complete)  # no frame counted twice


def test_interrupted_measurement_discarded_with_other_settings(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    monkeypatch.setattr(particle_finder_mod, "STATS_SAVE_EVERY_N_FRAMES", 1)
    _failing_after(N_FRAMES + 1, monkeypatch)
    sample_dir = _collect(tmp_path, n_par_target=10**6, **TWO_BANDS)
    interrupted_frames = {f.frame_id for f in _ledger(sample_dir).frames}

    _failing_after(10**6, monkeypatch)
    _collect(tmp_path, n_par_target=10**6, **{**TWO_BANDS, "min_par_pixels": 30})
    ledger = _ledger(sample_dir)
    assert [(m.measurement_id, m.status) for m in ledger.size_measurements] == [(0, "interrupted"), (1, "completed")]
    assert interrupted_frames <= {f.frame_id for f in ledger.frames}  # kept in the ledger ...
    result = ledger.size_result
    assert result.measurement_id == 1 and not interrupted_frames & {f for b in result.bands for f in b.frame_ids}
    assert result.powder_meas_cfg.min_par_pixels == 30
    assert [s.min_par_pixels for s in ledger.segmentations] == [20, 30]  # one segmentation per settings



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
    result = _ledger(sample_dir).size_result
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
    ledger = _ledger(sample_dir)
    result = ledger.size_result
    assert result.precision_reached is False and ledger.size_measurements[0].precision_reached is False
    assert sum(band.n_frames for band in result.bands) == 8
    # Both bands were extended; the finer one in subframes of the frames of the coarser one
    assert all(band.n_frames > 1 for band in result.bands)
    parents = {ledger.frame(f).label for f in result.bands[0].frame_ids}
    assert all(ledger.frame(f).label.rsplit("_", 1)[0] in parents for f in result.bands[1].frame_ids)



def test_frames_searched_by_eds_counted_in_the_statistics(tmp_path: Path, fake_microscope, monkeypatch):
    sample_dir = _collect(tmp_path, n_par_target=1)
    psd_result = _ledger(sample_dir).size_result
    n_psd = psd_result.size_stats.n_par_analysed
    eds_frames = [
        ((3.0, 20.0, 200.0), ParticleStatsFrame(
            frame_id=f"EDS-B{i}", frame_width_um=200.0, pixel_size_um=0.1, area_um2=2.4e4, n_particles=4,
            particle_areas_um2=[10.0, 20.0, 40.0, 80.0], particle_edge_corrections=[1.0, 1.0, 1.1, 1.2],
        ))
        for i in range(6)
    ]
    _run_eds(tmp_path, monkeypatch, eds_frames)

    ledger = _ledger(sample_dir)
    eds = ledger.frames_of("eds_search")  # saved with the spectra, segmented with the EDS settings
    assert [f.frame_id for f in eds] == [f"EDS-B{i}" for i in range(6)]
    eds_segmentation = eds[0].segmentations[0].segmentation_id
    assert eds_segmentation == psd_result.segmentation_id  # same segmentation settings as the measurement
    result = ledger.size_result
    # The EDS frames count the particles measured by EDS, which are left out; the statistics without them are kept
    assert (result.measurement_id, result.eds_segmentation_id) == (0, eds_segmentation)
    assert result.n_par_by_selection == {"size_stats": n_psd, "eds_frames": 24}
    assert len(result.eds_scans) == 1 and result.eds_scans[0].number_density_per_mm2 > 0
    assert result.size_stats.n_par_analysed == n_psd + 24
    assert len(ledger.size_results) == 2
    assert any(p.selection == "auto" for p in ledger.particles)  # still in the ledger



def test_interrupted_measurement_with_other_bands_restarted(tmp_path: Path, fake_microscope, monkeypatch):
    # E.g. interrupted with an earlier version of AutoEMX, which planned other size bands
    monkeypatch.setattr(_FakeEMController, "use_frame_navigator", True)
    monkeypatch.setattr(particle_finder_mod, "STATS_SAVE_EVERY_N_FRAMES", 1)
    _failing_after(N_FRAMES + 1, monkeypatch)
    sample_dir = _collect(tmp_path, n_par_target=10**6, **TWO_BANDS)
    ledger_path = sample_dir / "ledger.json"
    ledger = SampleLedger.from_json_file(str(ledger_path))
    ledger.size_measurements[0].bands[0].d_min_um *= 0.5  # same largest size, but wider
    ledger.to_json_file(str(ledger_path))

    _failing_after(10**6, monkeypatch)
    _collect(tmp_path, n_par_target=10**6, **TWO_BANDS)
    ledger = SampleLedger.from_json_file(str(ledger_path))
    assert [m.status for m in ledger.size_measurements] == ["interrupted", "completed"]
    result = ledger.size_result
    expected = plan_size_bands_um(_d_um(1e-3), _d_um(1e4), 1920, 20)
    assert [(b.d_min_um, b.d_max_um) for b in result.bands] == pytest.approx([(lo, hi) for lo, hi, _ in expected])
    assert [b.n_frames for b in result.bands] == [N_FRAMES, N_FRAMES]  # all frames analysed again



# --- Reanalysis from the saved frame images ------------------------------------------------------------------------

def _images_with_an_empty_frame(monkeypatch, empty_call: int = 2):
    """The example image, except an empty frame (no particle) at the empty_call-th image."""
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    n_calls = {"n": 0}

    def get_image_data(*args):
        n_calls["n"] += 1
        return np.zeros_like(image) if n_calls["n"] == empty_call else image

    monkeypatch.setattr(particle_finder_mod.EM_driver, "get_image_data", get_image_data, raising=False)


def test_every_frame_saved_and_reanalysed(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "write_images", True)
    _images_with_an_empty_frame(monkeypatch)
    sample_dir = _collect(tmp_path, n_par_target=10**6)
    images_dir = sample_dir / cnst.IMAGES_DIR
    ledger = _ledger(sample_dir)
    measured = ledger.size_result
    frames = ledger.frames
    # Every frame saved, the one without particles too
    assert [f.segmentation(0).n_particles == 0 for f in frames] == [False, True, False]
    assert all((images_dir / f.image_file).exists() for f in frames)
    # With the outlines of their particles (shown in the GUI)
    from autoemx.core.particle_outlines import load_outlines

    assert [len(load_outlines(str(sample_dir), 0, f.frame_id)["counted"]) for f in frames] == [
        f.segmentation(0).n_particles for f in frames]

    # Same settings: the segmentation of the measurement, same statistics
    reanalyze_particle_statistics(SAMPLE_ID, str(tmp_path), verbose=False)
    ledger = _ledger(sample_dir)
    assert len(ledger.segmentations) == 1 and len(ledger.size_results) == 1
    assert ledger.size_result.size_stats.model_dump() == pytest.approx(measured.size_stats.model_dump())
    assert (sample_dir / f"{SAMPLE_ID}_{cnst.PARTICLE_STATS_FILENAME}_reanalysed.csv").exists()

    # Brighter threshold: a new segmentation, added to the frames (the first one kept), and its statistics
    reanalyze_particle_statistics(SAMPLE_ID, str(tmp_path), par_brightness_thresh=200, verbose=False)
    ledger = _ledger(sample_dir)
    result = ledger.size_result
    segmentation = ledger.segmentation(result.segmentation_id)
    assert (segmentation.segmentation_id, segmentation.source, segmentation.par_brightness_thresh) == (
        1, "reanalysis", 200)
    assert all(f.segmentation(0) and f.segmentation(1) for f in ledger.frames)
    assert result.size_stats.n_par_analysed < measured.size_stats.n_par_analysed
    assert len(ledger.size_results) == 2 and ledger.size_results[0].size_stats == measured.size_stats

    # A frame whose image is missing cannot be segmented again: left out (listed in the result)
    (images_dir / frames[2].image_file).unlink()
    reanalyze_particle_statistics(SAMPLE_ID, str(tmp_path), par_brightness_thresh=150, verbose=False)
    result = _ledger(sample_dir).size_result
    assert result.bands[0].frame_ids == [frames[0].frame_id, frames[1].frame_id]
    assert result.missing_frame_ids == [frames[2].frame_id]


def test_eds_frames_reanalysed(tmp_path: Path, fake_microscope, monkeypatch):
    monkeypatch.setattr(_FakeEMController, "write_images", True)
    sample_dir = _collect(tmp_path, n_par_target=1)
    ledger_path = sample_dir / "ledger.json"
    ledger = SampleLedger.from_json_file(str(ledger_path))
    frame = ledger.frames[0]
    # An EDS search frame recorded without particles, whose saved image (the example image) has some
    ledger.add_frame(frame.model_copy(update={
        "frame_id": "EDS-B0", "label": "B0", "purpose": "eds_search", "measurement_id": None,
        "search_window_um": (_d_um(10.0), _d_um(100.0)),
        "segmentations": [FrameSegmentation(segmentation_id=0, d_range_um=(_d_um(10.0), _d_um(100.0)))],
    }))
    ledger.to_json_file(str(ledger_path))

    reanalyze_particle_statistics(SAMPLE_ID, str(tmp_path), par_brightness_thresh=110, verbose=False)
    ledger = SampleLedger.from_json_file(str(ledger_path))
    eds_frame = ledger.frame("EDS-B0")
    assert eds_frame.segmentation(0).n_particles == 0  # the recorded result is kept
    resegmented = eds_frame.segmentation(1)
    assert resegmented.n_particles > 0 and all(a >= 10.0 for a in resegmented.particle_areas_um2)
    assert ledger.size_result.eds_segmentation_id == 1
    assert ledger.size_result.n_par_by_selection["eds_frames"] == resegmented.n_particles


class _PositionedFrameNavigator(_FakeFrameNavigator):
    """As _FakeFrameNavigator, with frame positions: top-level frames side by side (A{i} at x = i frame widths),
    subframes at the centre of their parent. Frames overlapping areas already searched can be left out."""

    def __init__(self, im_h_to_w):
        super().__init__()
        self.im_h_to_w = im_h_to_w
        self._pos, self.parents = {}, []

    def configure_top_level_particle_frames(self, frame_width_mm):
        super().configure_top_level_particle_frames(frame_width_mm)
        self._pos = {f"A{i}": (i * frame_width_mm, 0.0) for i in range(N_FRAMES)}

    def configure_particle_subframes(self, parent_frames, child_fw_mm):
        super().configure_particle_subframes(parent_frames, child_fw_mm)
        self.parents = [label for label, _, _ in parent_frames]
        self._pos = {f"{label}_a0": tuple(center) for label, center, _ in parent_frames}

    def skip_frames_overlapping(self, areas_mm):
        areas = list(areas_mm)
        half = (self.grid_search_fw_mm / 2, self.grid_search_fw_mm * self.im_h_to_w / 2)
        skipped = [label for label in self._labels if areas and np.any(frames_overlap(
            self._pos[label], half, [center for center, _ in areas], [h for _, h in areas]))]
        self.skip_frames(skipped)
        return len(skipped)

    def go_to_next_frame(self):
        moved = super().go_to_next_frame()
        if moved:
            self._current_pos = self._pos[self.current_frame_label]
        return moved


def _size_stats_after_eds(tmp_path: Path, monkeypatch, eds_frames, n_par_target: int, size_cfg_kwargs=None,
                          **powder_meas_cfg_kwargs):
    """EDS acquisition, then particle size statistics counting the frames searched by EDS, on the same analyzer
    (whose microscope is replaced by one finding particles in positioned frames)."""
    analyzer = _run_eds(tmp_path, monkeypatch, eds_frames)
    image, _ = load_grayscale_image(INPUTS_DIR / "example_frame_image.tif")
    controller = _FakeEMController(LedgerModule._build_ledger_configs(analyzer), SAMPLE_ID,
                                   str(tmp_path / SAMPLE_ID / cnst.IMAGES_DIR))
    controller.frame_navigator = _PositionedFrameNavigator(image.shape[0] / image.shape[1])
    cfg = controller.particle_finder.powder_meas_cfg
    for name, value in {"min_area_par": 1.0, "max_area_par": 100.0, "min_par_pixels": 10,
                        **powder_meas_cfg_kwargs}.items():
        setattr(cfg, name, value)
    analyzer.EM_controller = controller
    saved_cfg = cfg.model_copy(deep=True)
    analyzer.collect_particle_stats_after_eds(n_par_target, size_cfg_kwargs)
    assert cfg == saved_cfg and controller.measurement_cfg.type == "EDS"  # settings of the EDS acquisition restored
    return SampleLedger.from_json_file(str(tmp_path / SAMPLE_ID / "ledger.json")), controller


def _eds_frame(frame_id, center_mm, width_um, diameters_um, d_range=(_d_um(1.0), _d_um(100.0))):
    areas = [float(np.pi * d**2 / 4) for d in diameters_um]
    return ((d_range[0], d_range[1], width_um), ParticleStatsFrame(
        frame_id=frame_id, center_mm={"x": center_mm[0], "y": center_mm[1]}, frame_width_um=width_um, pixel_size_um=width_um / 1000,
        area_um2=width_um**2 * 0.75, n_particles=len(areas), particle_areas_um2=areas,
        particle_edge_corrections=[1.0] * len(areas)))


def test_size_statistics_continuing_eds(tmp_path: Path, fake_microscope, monkeypatch):
    # Particles found per frame (the same image in every frame), without frames searched by EDS
    ref, _ = _size_stats_after_eds(tmp_path / "ref", monkeypatch, [], n_par_target=1)
    assert ref.size_result.bands[0].frame_ids == ["A0"]
    n_frame = ref.frame("A0").segmentation(ref.size_result.segmentation_id).n_particles
    # EDS searched a small frame at the position of frame A0, finding 3 particles of the band (and one larger)
    eds_frames = [_eds_frame("EDS-A0", (0.0, 0.0), 10.0, [5.0, 5.0, 6.0, 50.0])]
    ledger, _ = _size_stats_after_eds(tmp_path, monkeypatch, eds_frames, n_par_target=n_frame + 2)

    result = ledger.size_result
    [band] = result.bands
    # A0, searched by EDS, is not searched again; with the 3 EDS particles of the band, one frame reaches the target
    # (two would be needed without them)
    assert band.frame_ids == ["A1"]
    assert result.eds_scans[0].frame_ids == ["EDS-A0"] and result.eds_segmentation_id == result.segmentation_id
    assert result.overlapping_frame_ids == []
    assert result.n_par_by_selection == {"size_stats": n_frame, "eds_frames": 4}  # every particle of the EDS frame
    [measurement] = ledger.size_measurements
    assert measurement.count_eds_frames and measurement.status == "completed"
    assert ledger.configs.measurement_cfg.type == "EDS"


def test_size_statistics_after_eds_target_reached_by_eds(tmp_path: Path, fake_microscope, monkeypatch):
    # EDS frames already hold enough particles of the band: no frame is searched
    eds_frames = [_eds_frame(f"EDS-B{i}", (5.0 + i, 5.0), 10.0, [5.0] * 4) for i in range(3)]
    ledger, _ = _size_stats_after_eds(tmp_path, monkeypatch, eds_frames, n_par_target=10)
    result = ledger.size_result
    assert [b.n_frames for b in result.bands] == [0] and ledger.frames_of("size_stats") == []
    assert result.n_par_by_selection == {"eds_frames": 12} and result.size_stats.n_par_analysed == 12


def test_frames_searched_by_eds_are_parents_of_the_second_band(tmp_path: Path):
    finder = EM_Particle_Finder.__new__(EM_Particle_Finder)
    finder.verbose = False
    finder._eds_frames = [particle_statistics.EdsSearchFrame("EDS-A0", (2000.0, 0.0), (100.0, 75.0), (5.0, 50.0), 200.0)]
    nav = _PositionedFrameNavigator(0.75)
    band0 = ParticleSizeBand(d_min_um=5.0, d_max_um=50.0, frame_width_um=1000.0, n_frames=0, scanned_area_um2=0)
    band0.add_frame(ParticleStatsFrame(frame_id="A0", center_mm={"x": 0.0, "y": 0.0}, frame_width_um=1000.0, pixel_size_um=1.0,
                                       area_um2=7.5e5, n_particles=0))
    progress = particle_finder_mod.StatsProgress(n_par_target=5, powder_meas_cfg=PowderMeasurementConfig(), bands=[band0])
    # Second band, below the sizes searched by EDS: searched in the EDS frame too
    assert finder._plan_band_frames(nav, progress, 1, 0.1, d_range_um=(1.0, 5.0))
    assert nav.parents == ["A0", "EDS-A0"] and nav._labels == ["A0_a0", "EDS-A0_a0"]
    # Sizes also searched by EDS: subframes inside the EDS frame left out
    assert finder._plan_band_frames(nav, progress, 1, 0.1, d_range_um=(3.0, 6.0))
    assert nav._labels == ["A0_a0"]


class _FakeRunnerAnalyzer:
    """Composition analyzer without microscope recording what runners call; IDs with 'fail' fail."""

    calls: list = []

    def __init__(self, sample_id, **_):
        self.sample_id = sample_id
        finder = SimpleNamespace(get_particle_stats=lambda n: self._call("stats", n))
        self.EM_controller = SimpleNamespace(particle_finder=finder)

    def _call(self, name, *args):
        if "fail" in self.sample_id:
            raise RuntimeError("no particles found")
        _FakeRunnerAnalyzer.calls.append((self.sample_id, name, *args))

    def run_collection_and_quantification(self, **_):
        _FakeRunnerAnalyzer.calls.append((self.sample_id, "eds"))

    def collect_particle_stats_after_eds(self, n_par_target, powder_meas_cfg_kwargs=None):
        self._call("stats_after_eds", n_par_target, *([powder_meas_cfg_kwargs] if powder_meas_cfg_kwargs else []))


def test_runners_measure_size_statistics_after_eds_and_report_each_sample(tmp_path: Path, monkeypatch):
    import autoemx.runners.batch_acquire_and_analyze as batch_runner
    import autoemx.runners.collect_particle_statistics as psd_runner

    monkeypatch.setattr(batch_runner, "EMXSp_Composition_Analyzer", _FakeRunnerAnalyzer)
    monkeypatch.setattr(psd_runner, "EMXSp_Composition_Analyzer", _FakeRunnerAnalyzer)
    samples = [{"ID": sid, "els": ["Si"], "pos": (0, 0)} for sid in ("A", "B_fail")]
    done = []
    _FakeRunnerAnalyzer.calls = []
    batch_runner.batch_acquire_and_analyze(samples, development_mode=True, verbose=False, results_dir=str(tmp_path),
                                           size_stats_n_par_target=50, size_stats_powder_meas_cfg_kwargs={"x": 1},
                                           on_sample_done=lambda s, e: done.append((s, e)))
    assert _FakeRunnerAnalyzer.calls == [("A", "eds"), ("A", "stats_after_eds", 50, {"x": 1}), ("B_fail", "eds")]
    assert done[0] == ("A", None) and isinstance(done[1][1], RuntimeError)  # size statistics failed: reported

    _FakeRunnerAnalyzer.calls, done = [], []
    batch_runner.batch_acquire_and_analyze(samples[:1], development_mode=True, verbose=False,
                                           results_dir=str(tmp_path))
    assert _FakeRunnerAnalyzer.calls == [("A", "eds")]  # EDS only by default

    _FakeRunnerAnalyzer.calls = []
    psd_runner.collect_particle_statistics(samples, n_par_target=20, development_mode=True, verbose=False,
                                           results_dir=str(tmp_path), on_sample_done=lambda s, e: done.append((s, e)))
    assert _FakeRunnerAnalyzer.calls == [("A", "stats", 20)]
    assert done[0] == ("A", None) and isinstance(done[1][1], RuntimeError)


def test_size_statistics_after_eds_with_their_own_particle_settings(tmp_path: Path, fake_microscope, monkeypatch):
    # Size statistics of smaller particles than those of the EDS acquisition (settings restored after)
    ledger, _ = _size_stats_after_eds(tmp_path, monkeypatch, [], n_par_target=1, size_cfg_kwargs={"max_area_par": 20.0})
    [measurement] = ledger.size_measurements
    assert measurement.bands[0].d_max_um == pytest.approx(_d_um(20.0))
    assert measurement.powder_meas_cfg.max_area_par == ledger.size_result.powder_meas_cfg.max_area_par == 20.0


def test_eds_frames_segmented_again_with_the_settings_of_the_size_statistics(tmp_path: Path, fake_microscope,
                                                                            monkeypatch):
    import shutil

    from autoemx.core.saved_frames import SavedFrames, segment_frames

    # A saved frame image (the example image), as EDS would have saved it
    monkeypatch.setattr(_FakeEMController, "write_images", True)
    ref_frame = _ledger(_collect(tmp_path / "ref", n_par_target=1)).frames[0]
    images_dir = tmp_path / SAMPLE_ID / cnst.IMAGES_DIR
    images_dir.mkdir(parents=True)
    shutil.copy(tmp_path / "ref" / SAMPLE_ID / cnst.IMAGES_DIR / ref_frame.image_file, images_dir)
    monkeypatch.setattr(_FakeEMController, "write_images", False)
    # EDS frames: one with the saved image (recorded with other particles), one without image
    with_image = _eds_frame("EDS-C0", (50.0, 50.0), ref_frame.frame_width_um, [5.0] * 4)
    with_image = (with_image[0], with_image[1].model_copy(update={
        "image_file": ref_frame.image_file, "pixel_size_um": ref_frame.pixel_size_um}))
    no_image = _eds_frame("EDS-C1", (60.0, 50.0), ref_frame.frame_width_um, [5.0] * 4)

    ledger, _ = _size_stats_after_eds(tmp_path, monkeypatch, [with_image, no_image], n_par_target=1,
                                      size_cfg_kwargs={"par_brightness_thresh": 150})
    result = ledger.size_result
    size_segmentation = ledger.segmentation(result.segmentation_id)
    assert size_segmentation.par_brightness_thresh == 150 and result.eds_segmentation_id == result.segmentation_id
    # Segmented again with the threshold of the size statistics (the EDS result kept); the frame without image left
    # out of the statistics
    frame = ledger.frame("EDS-C0")
    assert frame.segmentation(0).particle_areas_um2 == with_image[1].particle_areas_um2
    resegmented = frame.segmentation(size_segmentation.segmentation_id)
    expected, _ = segment_frames([frame], SavedFrames(SAMPLE_ID, "PhenomXL"),
                                 PowderMeasurementConfig(min_par_pixels=10, par_brightness_thresh=150), str(images_dir),
                                 size_segmentation.segmentation_id)
    assert resegmented.particle_areas_um2 == expected["EDS-C0"].particle_areas_um2 and resegmented.n_particles > 0
    assert resegmented.particle_areas_um2 != with_image[1].particle_areas_um2
    assert ledger.frame("EDS-C1").segmentation(size_segmentation.segmentation_id) is None
    assert "EDS-C1" in result.missing_frame_ids
    assert result.n_par_by_selection["eds_frames"] == resegmented.n_particles


def test_eds_frames_not_segmented_again_with_the_same_settings(tmp_path: Path, fake_microscope, monkeypatch):
    eds = _eds_frame("EDS-C0", (50.0, 50.0), 10.0, [5.0] * 4)
    # Same segmentation settings as EDS: the frames keep the particles recorded by EDS
    ledger, _ = _size_stats_after_eds(tmp_path, monkeypatch, [eds], n_par_target=1,
                                      size_cfg_kwargs={"par_segmentation_model": "threshold_bright"})
    [frame] = ledger.frames_of("eds_search")
    assert len(frame.segmentations) == 1 and frame.segmentations[0].particle_areas_um2 == eds[1].particle_areas_um2
    assert ledger.size_result.eds_segmentation_id == frame.segmentations[0].segmentation_id
