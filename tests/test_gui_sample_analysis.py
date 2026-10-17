#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the sample-analysis GUI backend and figures (no quantification is run)."""

from __future__ import annotations

import inspect
import shutil
from pathlib import Path

import pytest

import autoemx.utils.constants as cnst

pytest.importorskip("plotly")

from autoemx.config.ledger_schemas import (
    AitchisonParams,
    ClusterMergeParams,
    DBSCANParams,
    MixtureParams,
)
from autoemx.gui import backend as be
from autoemx.gui import plots as pl
from autoemx.runners.analyze_sample import analyze_sample
from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze
from ci_paths import INPUTS_DIR, K412_CLUSTER_MINI_ID, TESTS_DIR, WULFENITE_MINI_ID

POWDER_EXAMPLE = TESTS_DIR.parent / "examples" / "Results" / "Powder_mixture_example"


@pytest.fixture()
def results_dir(tmp_path: Path) -> Path:
    for sample in (WULFENITE_MINI_ID, K412_CLUSTER_MINI_ID):
        shutil.copytree(INPUTS_DIR / sample, tmp_path / sample)
    return tmp_path


@pytest.fixture()
def powder_dir(tmp_path: Path) -> Path:
    if not (POWDER_EXAMPLE / "spectra").is_dir():
        pytest.skip("Powder_mixture_example spectra not available")
    dest = tmp_path / "res" / POWDER_EXAMPLE.name
    shutil.copytree(POWDER_EXAMPLE, dest)
    return dest


def test_find_samples(results_dir: Path):
    samples = be.find_samples(str(results_dir))
    assert [s.sample_id for s in samples] == [K412_CLUSTER_MINI_ID, WULFENITE_MINI_ID]
    with pytest.raises(FileNotFoundError):
        be.find_samples(str(results_dir / "missing"))


def test_param_specs_cover_every_clustering_submodel_field():
    keys = set(be.SPECS_BY_KEY)
    for section, model in (("dbscan", DBSCANParams), ("aitchison", AitchisonParams),
                           ("merge", ClusterMergeParams), ("mixture", MixtureParams)):
        for name in be._SUBMODEL_FIELDS.get(section) or model.model_fields:
            assert f"{section}.{name}" in keys
    # Mixture decomposition: only the main parameters are exposed
    assert {k for k in keys if k.startswith("mixture.")} == {
        "mixture.max_n_phases", "mixture.collapse_equivalent_mixtures", "mixture.equivalent_span_tol"}


def test_kwargs_match_runner_signatures(results_dir: Path):
    from autoemx.runners.fit_and_quantify_spectrum_from_ledger import fit_and_quantify_spectrum_from_ledger

    info = be.load_sample_info(str(results_dir / WULFENITE_MINI_ID))
    values = be.coerce_values(be.sample_param_values(info))
    assert set(be.analysis_kwargs(values)) <= set(inspect.signature(analyze_sample).parameters)

    qvalues = be.coerce_values({s.key: s.default for s in be.QUANT_SPECS}, be.QUANT_SPECS)
    qkwargs = be.quantification_kwargs(qvalues)
    assert set(qkwargs) <= set(inspect.signature(batch_quantify_and_analyze).parameters)
    # Empty / 'saved' options keep each sample's saved value
    assert qkwargs["spectrum_lims"] is None and qkwargs["use_project_specific_std_dict"] is None
    assert "use_instrument_background" not in qkwargs and qkwargs["run_analysis"] is True
    assert qkwargs["max_spectra_to_quantify"] is None

    svalues = be.coerce_values(be.single_param_values(info), be.SINGLE_SPECS)
    skwargs = be.single_fit_kwargs(svalues)
    skwargs.pop("quantify")
    assert set(skwargs) <= set(inspect.signature(fit_and_quantify_spectrum_from_ledger).parameters)
    assert skwargs["els_sample"] == ["Pb", "Mo", "O"]


def test_coerce_values_reports_all_errors(results_dir: Path):
    info = be.load_sample_info(str(results_dir / WULFENITE_MINI_ID))
    raw = be.sample_param_values(info)
    raw.update({"dbscan.eps": "-1", "clust.ref_formulae": "PbMoO4\nXx2"})
    with pytest.raises(ValueError) as exc:
        be.coerce_values(raw)
    msg = str(exc.value)
    assert "eps" in msg and "Xx2" in msg
    single = be.single_param_values(info)
    single.update({"single.els_sample": "Pb, Qq", "single.std_formula": "Zz3"})
    with pytest.raises(ValueError) as exc:
        be.coerce_values(single, be.SINGLE_SPECS)
    assert "Qq" in str(exc.value) and "Zz3" in str(exc.value)


def test_k_forced_mapping(results_dir: Path):
    info = be.load_sample_info(str(results_dir / WULFENITE_MINI_ID))
    values = be.coerce_values({**be.sample_param_values(info), "clust.k_forced": 3})
    kwargs = be.analysis_kwargs(values)
    assert kwargs["k_forced"] == 3 and kwargs["k_finding_method"] is None
    values = be.coerce_values({**be.sample_param_values(info), "clust.k_forced": None})
    kwargs = be.analysis_kwargs(values)
    assert kwargs["k_forced"] is False and kwargs["k_finding_method"] == values["clust.k_finding_method"]


def test_analysis_and_figures(powder_dir: Path):
    info = be.load_sample_info(str(powder_dir))
    values = be.coerce_values(be.sample_param_values(info))
    values["clust.ref_formulae"] = ["MgO", "Al2O3", "MgAl2O4"]
    values["clust.k_forced"] = 2
    analyzer = analyze_sample(
        sample_ID=powder_dir.name, results_path=str(powder_dir.parent), **be.analysis_kwargs(values)
    )
    assert analyzer is not None

    info = be.load_sample_info(str(powder_dir))
    data = be.load_analysis(info, info.active_analysis)
    assert data.summary["n_clusters"] == 2
    assert data.summary["n_clustered"] == sum(data.n_points)
    assert set(data.comps["status"]) <= {"clustered", "discarded", "not quantified"}
    assert list(data.ref_comps.index) == ["MgO", "Al2O3", "MgAl2O4"]
    assert data.ref_comps.loc["MgAl2O4", "Mg"] == pytest.approx(1 / 7)

    sid = data.comps.loc[data.comps["status"] == "clustered", "spectrum"].iloc[0]
    for mode, axes in (("3d", ["Mg", "Al", "O"]), ("ternary", ["Mg", "Al", "O"]), ("2d", ["Mg", "Al"])):
        fig = pl.clustering_figure(data, axes, mode, highlight=sid)
        names = {t.name for t in fig.data}
        assert {"Cluster 0", "Centroids", "Candidate phases", f"Spectrum {sid}"} <= {
            n.split(" (")[0] if n else n for n in names
        }
    assert pl.distribution_figure(data).data

    raw = be.load_raw_spectrum(info, sid)
    assert len(raw["energy"]) == len(raw["counts"]) > 0
    assert pl.spectrum_figure(raw).data


def test_clustering_zoom(powder_dir: Path):
    info = be.load_sample_info(str(powder_dir))
    data = be.load_analysis(info, info.active_analysis)
    assert [o["value"] for o in pl.zoom_options(None)] == [pl.ZOOM_FULL, pl.ZOOM_DATA]
    options = pl.zoom_options(data)
    assert [o["value"] for o in options][2:] == [f"{pl.ZOOM_CLUSTER}{i}" for i in range(len(data.centroids))]
    assert options[2]["label"] == f"Cluster 0 ({data.n_points[0]})"

    axes = ["Mg", "Al", "O"]
    members = data.comps[(data.comps["status"] == "clustered") & (data.comps["cluster"] == 0)]
    pts = members[axes].to_numpy(dtype=float) * 100
    spans = {}
    for zoom in (pl.ZOOM_FULL, pl.ZOOM_DATA, f"{pl.ZOOM_CLUSTER}0"):
        scene = pl.clustering_figure(data, axes, "3d", zoom=zoom).layout.scene
        ranges = [sorted(scene[f"{a}axis"].range) for a in "xyz"]
        assert all(lo <= pts[:, i].min() and pts[:, i].max() <= hi for i, (lo, hi) in enumerate(ranges))
        spans[zoom] = [hi - lo for lo, hi in ranges]
    assert spans[pl.ZOOM_FULL] == [106, 106, 106]
    assert all(c < d <= f for c, d, f in zip(spans[f"{pl.ZOOM_CLUSTER}0"], spans[pl.ZOOM_DATA],
                                             spans[pl.ZOOM_FULL]))

    # Ternary: the zoomed triangle contains the (normalised) cluster; 2D: the axes contain it
    tern = pl.clustering_figure(data, axes, "ternary", zoom=f"{pl.ZOOM_CLUSTER}0").layout.ternary
    norm = pts / pts.sum(axis=1, keepdims=True) * 100
    mins = [tern.baxis.min, tern.caxis.min, tern.aaxis.min]  # b: x, c: y, a: z (as in _point_trace)
    assert sum(mins) > 0 and all(m <= norm[:, i].min() for i, m in enumerate(mins))
    lay = pl.clustering_figure(data, axes[:2], "2d", zoom=f"{pl.ZOOM_CLUSTER}0").layout
    for i, ax in enumerate((lay.xaxis, lay.yaxis)):
        assert ax.range[0] <= pts[:, i].min() and pts[:, i].max() <= ax.range[1] and ax.range[1] - ax.range[0] < 100
    # Unknown cluster: zoomed on all spectra
    assert pl.clustering_figure(data, axes, "3d", zoom=f"{pl.ZOOM_CLUSTER}99").layout.scene.xaxis.range == \
        pl.clustering_figure(data, axes, "3d", zoom=pl.ZOOM_DATA).layout.scene.xaxis.range


def test_low_confidence_mixtures_hidden_by_default(powder_dir: Path):
    info = be.load_sample_info(str(powder_dir))
    data = be.load_analysis(info, info.active_analysis)
    data.ref_formulae, data.ref_comps = ["MgO", "Al2O3"], be.reference_compositions(
        ["MgO", "Al2O3"], data.elements, data.features, ["H", "He", "Li"])
    data.mixtures = [[{"refs": ["MgO", "Al2O3"], "conf_score": 0.5}]]

    def n_lines(**kw):
        fig = pl.clustering_figure(data, ["Mg", "Al", "O"], "3d", **kw)
        return sum(1 for t in fig.data if t.legendgroup == "mixtures")

    assert n_lines() == 0
    assert n_lines(min_mixture_conf=0.4) == 1


def test_create_launcher(tmp_path: Path):
    import os
    import sys

    from autoemx.gui.__main__ import create_launcher

    ext = ".bat" if os.name == "nt" else ".command" if sys.platform == "darwin" else ".sh"
    old_launcher = tmp_path / f"AutoEMX GUI{ext}"
    old_launcher.write_text('python -m autoemx.gui\n')
    other_file = tmp_path / f"AutoEMX GUI notes{ext}"
    other_file.write_text("not a launcher")

    path = create_launcher(str(tmp_path), results_folder=str(tmp_path))
    assert path.name == f"AutoEMX{ext}"
    assert not old_launcher.exists() and other_file.exists()  # only the old launcher is replaced
    text = path.read_text()
    assert sys.executable in text and "-m autoemx.gui" in text and str(tmp_path) in text
    if os.name != "nt":
        assert os.access(path, os.X_OK)


def test_sem_images_linked_to_spectra(tmp_path: Path):
    from PIL import Image

    images_dir = tmp_path / "S" / "SEM images"
    images_dir.mkdir(parents=True)
    for name in ("S_par10_frE1_xyspots.png", "S_par2_frE1_xyspots.png", "S_frE1_particles.png",
                 "S_fr3_xyspots.tif", "S_par2_frE1_xyspots_raw.png"):
        Image.new("RGB", (40, 30)).save(images_dir / name)
    images = be.find_sem_images(str(tmp_path / "S"))
    assert [im.name for im in images] == [
        "S_fr3_xyspots.tif", "S_par2_frE1_xyspots.png", "S_par10_frE1_xyspots.png", "S_frE1_particles.png"]
    assert images[3].kind == "frame" and images[3].label == "Frame E1"
    assert be.sem_image_for_spectrum(images, 10, "E1").name == "S_par10_frE1_xyspots.png"
    assert be.sem_image_for_spectrum(images, None, "3").name == "S_fr3_xyspots.tif"
    assert be.sem_image_for_spectrum(images, 7, "E1") is None
    png, size = be.read_image_png(images[0].path, 20)
    assert size == (40, 30) and png[:4] == b"\x89PNG"

    comps = be.pd.DataFrame({"spectrum": ["0", "1", "2"], "particle": [2, 10, 10], "frame": ["E1"] * 3,
                             "px": [5, 6, None], "py": [5, 6, 7]})
    assert be.spectra_on_image(comps, images[2])["spectrum"].tolist() == ["1"]


def test_load_analysis_without_clustering(results_dir: Path):
    info = be.load_sample_info(str(results_dir / K412_CLUSTER_MINI_ID))
    data = be.load_analysis(info, None)
    assert len(data.comps) == info.n_spectra
    assert data.summary["n_clusters"] == 0
    fig = pl.clustering_figure(data, data.elements[:3], "3d")
    assert fig.data


def test_active_params():
    from autoemx.gui.app import active_params

    base = {"clust.method": "kmeans", "clust.geometry": "euclidean", "clust.k_forced": None,
            "clust.auto_merge_clusters": True, "clust.do_matrix_decomposition": True}
    sections, params = active_params(base)
    assert not sections["dbscan"] and not sections["aitchison"]
    assert sections["merge"] and sections["mixture"] and params["clust.max_k"]
    for geometry in ("aitchison", "auto"):
        assert active_params({**base, "clust.geometry": geometry})[0]["aitchison"]
    sections, params = active_params({**base, "clust.k_forced": 3})
    assert not sections["merge"] and not params["clust.k_finding_method"]
    sections, params = active_params({**base, "clust.method": "dbscan"})
    assert sections["dbscan"] and not sections["merge"] and not params["clust.k_forced"]
    assert not active_params({**base, "clust.do_matrix_decomposition": False})[0]["mixture"]


def test_sample_summary_and_quant_progress(results_dir: Path):
    summary = be.sample_summary(str(results_dir / K412_CLUSTER_MINI_ID))
    assert summary["sample"] == K412_CLUSTER_MINI_ID and summary["n_spectra"] == 6
    assert summary["elements"] == "Fe, Mg, Ca, Al, Si, O" and len(summary["date"]) == 16
    runs = be.quantification_runs(be.load_sample_info(str(results_dir / K412_CLUSTER_MINI_ID)))
    assert runs and runs[-1]["n_spectra"] == 6

    log = (f"{be.QUANT_SAMPLE_MARKER} 1/2: A\nINFO: Starting quantification of 3 spectra on up to 6 cores.\n"
           f" Spectrum #1/5:\n  x\n Spectrum #0/5:\n\n{be.QUANT_SAMPLE_MARKER} 2/2: B\n"
           "INFO: Starting quantification of 4 spectra\n Spectrum #3/9:\n")
    progress = be.quant_progress(log, ["A", "B"])
    assert progress[0] == {"sample": "A", "state": "done", "done": 2, "total": 3}
    assert progress[1] == {"sample": "B", "state": "running", "done": 1, "total": 4}


def test_acquisition_settings():
    from autoemx.gui.tab_acquisition import active_acq_params
    from autoemx.runners.batch_acquire_and_analyze import batch_acquire_and_analyze

    defaults = {s.key: s.default for s in be.ACQ_SPECS}
    values = be.coerce_values(defaults, be.ACQ_SPECS)
    kwargs = be.acquisition_kwargs(values)
    assert set(kwargs) <= set(inspect.signature(batch_acquire_and_analyze).parameters)
    assert kwargs["max_XSp_acquisition_time"] == 25  # 50000 counts / 10000 * 5 s, as in Run_Acquisition.py
    assert kwargs["els_substrate"] == ["C"]
    # Without quantification, "number of spectra" sets the spectra collected; with it, min and max
    # (also as the minimum: the analyser raises max_n_spectra to min_n_spectra)
    no_quant = be.acquisition_kwargs({**values, "aacq.n_spectra": 5})
    assert (no_quant["min_n_spectra"], no_quant["max_n_spectra"]) == (5, 5)
    quant = be.acquisition_kwargs({**values, "aacq.quantify_spectra": True, "aacq.n_spectra": 5})
    assert (quant["min_n_spectra"], quant["max_n_spectra"]) == (50, 100)
    assert kwargs["powder_meas_cfg_kwargs"]["max_area_par"] == 10000.0
    assert "par_spot_selection_mode" not in kwargs["powder_meas_cfg_kwargs"]

    samples = be.acquisition_samples([
        {"ID": "Anorthite", "els": "Ca, Al, Si, O", "x": "-37.5", "y": -37.5, "cnd": "CaAl2Si2O8"},
        {"ID": "", "els": "", "x": None, "y": None, "cnd": ""},  # empty rows are ignored
    ])
    assert samples == [{"ID": "Anorthite", "els": ["Ca", "Al", "Si", "O"], "pos": (-37.5, -37.5),
                        "cnd": ["CaAl2Si2O8"]}]
    with pytest.raises(ValueError) as exc:
        be.acquisition_samples([{"ID": "a/b", "els": "Qq", "x": "s", "y": 1, "cnd": "Xx2"},
                                {"ID": "a/b", "els": "O", "x": 0, "y": 0}])
    assert all(word in str(exc.value) for word in ("cannot contain", "Qq", "numbers", "Xx2", "duplicated"))

    script = be.acquisition_script(samples, kwargs, "/data/results")
    compile(script, "Run_Acquisition_GUI.py", "exec")
    assert "batch_acquire_and_analyze(" in script and "'Anorthite'" in script

    active = active_acq_params(defaults)
    assert (active["apowder"], active["abulk"], active["aquant"]) == ("on", "off", "off")
    assert (active["aacq.n_spectra"], active["aacq.min_n_spectra"], active["aacq.contrast"]) == ("on", "hidden", "hidden")
    active = active_acq_params({**defaults, "asample.sample_type": "bulk", "aacq.quantify_spectra": True,
                                "aacq.auto_adjust_brightness_contrast": False})
    assert (active["apowder"], active["abulk"], active["aquant"]) == ("off", "on", "on")
    assert (active["aacq.n_spectra"], active["aacq.max_n_spectra"], active["aacq.brightness"]) == ("hidden", "on", "on")

    log = "Sample 'Anorthite'\n🔬 Acquiring spectrum #0...\n🔬 Acquiring spectrum #1...\nSample 'B'\n🔬 Acquiring spectrum #2..."
    progress = be.acquisition_progress(log, ["Anorthite", "B", "C"], 100)
    assert [(p["state"], p["done"]) for p in progress] == [("done", 2), ("running", 1), ("waiting", 0)]
    log = ("Sample 'A'\n🔬 Acquiring spectrum #0...\n12:00:00 ERROR: Sample 'A': acquisition/quantification failed: "
           "no particles found\nTraceback...\nSample 'B'\n🔬 Acquiring spectrum #0...")
    progress = be.acquisition_progress(log, ["A", "B"], 100)
    assert [(p["state"], p["error"]) for p in progress] == [("failed", "no particles found"), ("running", None)]


def test_import_spectra_folder(tmp_path: Path):
    source = INPUTS_DIR / WULFENITE_MINI_ID / cnst.SPECTRA_DIR
    scan = be.inspect_spectra_folder(str(source))
    n_files = len([p for p in source.iterdir() if p.suffix.lower() in cnst.EMSA_SPECTRUM_EXTENSIONS])
    assert scan["n_files"] == n_files > 0 and scan["beam_energies"] == [15.0] and scan["calibration"]

    results = tmp_path / "results"
    results.mkdir()
    (results / "taken").mkdir()
    args = dict(folder=str(source), results_dir=str(results), elements="Pb, Mo, O", substrate="C, O, Al",
                sample_type="powder", microscope_id=be.dflt.microscope_ID, beam_energy=15)
    for bad, msg in ((dict(sample_id="taken"), "already exists"), (dict(sample_id="a b"), "sample ID"),
                     (dict(sample_id="ok", elements=""), "elements"), (dict(sample_id="ok", beam_energy=None), "beam")):
        with pytest.raises(ValueError, match=msg):
            be.import_kwargs(**{**args, **bad})

    kwargs = be.import_kwargs(sample_id="Imported", **args)
    be._run_import("", {"kwargs": kwargs})
    info = be.load_sample_info(str(results / "Imported"))
    assert info.n_spectra == n_files and info.elements == ["Pb", "Mo", "O"]
    cfg = info.ledger.configs
    assert (cfg.microscope_cfg.energy_zero, cfg.microscope_cfg.bin_width) == pytest.approx(scan["calibration"])
    assert cfg.sample_cfg.type == "powder" and cfg.measurement_cfg.beam_energy_keV == 15
    assert len(list(source.iterdir())) == n_files  # source untouched
    # Failed import: no partial sample folder left
    with pytest.raises(Exception):
        be._run_import("", {"kwargs": {**kwargs, "samples": [{**kwargs["samples"][0], "ID": "Empty",
                                                                 "spectra_dir": str(tmp_path)}]}})
    assert not (results / "Empty").exists()


def test_import_requires_one_energy_calibration(tmp_path: Path):
    source = INPUTS_DIR / WULFENITE_MINI_ID / cnst.SPECTRA_DIR
    files = sorted(p for p in source.iterdir() if p.suffix.lower() in cnst.EMSA_SPECTRUM_EXTENSIONS)[:2]
    mixed = tmp_path / "mixed"
    mixed.mkdir()
    for i, f in enumerate(files):
        text = f.read_text(encoding="utf-8", errors="replace")
        if i:  # second file: another channel width
            text = "\n".join("#XPERCHAN    : 5.0" if line.upper().startswith("#XPERCHAN") else line
                              for line in text.splitlines())
        (mixed / f.name).write_text(text, encoding="utf-8")
    scan = be.inspect_spectra_folder(str(mixed))
    assert scan["calibration"] is None and "different energy calibrations" in scan["calibration_error"]

    results = tmp_path / "results"
    results.mkdir()
    kwargs = be.import_kwargs(str(mixed), str(results), "Mixed", "Pb, Mo, O", "C", "powder",
                              be.dflt.microscope_ID, 15)
    with pytest.raises(RuntimeError):
        be._run_import("", {"kwargs": kwargs})
    assert not (results / "Mixed").exists()


def test_windows_launcher_shortcut():
    from PIL import Image

    from autoemx.gui.__main__ import LAUNCHER_ICON_WINDOWS, _windows_shortcut_script

    # Windows icon with the sizes shown on the Desktop, up to 256 px
    assert (256, 256) in Image.open(LAUNCHER_ICON_WINDOWS).info["sizes"]
    script = _windows_shortcut_script(Path("C:/Users/O'Neil/Desktop/AutoEMX.lnk"), Path("C:/x/AutoEMX_1.bat"),
                                      LAUNCHER_ICON_WINDOWS, Path("C:/Users/O'Neil"))
    assert "'C:/Users/O''Neil/Desktop/AutoEMX.lnk'" in script.replace("\\", "/")  # quotes escaped
    assert f"{LAUNCHER_ICON_WINDOWS},0" in script and script.endswith("$s.Save()")


def test_remove_unfinished_sample(tmp_path: Path):
    unfinished = tmp_path / "unfinished"
    (unfinished / cnst.SPECTRA_DIR).mkdir(parents=True)
    finished = tmp_path / "finished"
    finished.mkdir()
    (finished / be.LEDGER_NAME).write_text("{}", encoding="utf-8")
    for path in (unfinished, finished, tmp_path / "missing", None):
        be.remove_unfinished_sample(str(path) if path else None)
    assert not unfinished.exists() and finished.exists()


# ----------------------------------------------------------------------------- particle sizes
@pytest.fixture()
def psd_dir(tmp_path: Path) -> Path:
    from _psd_helpers import make_eds_frames_sample, make_empty_sample, make_measurement_sample

    make_measurement_sample(tmp_path / "PSD_only")
    make_eds_frames_sample(tmp_path / "EDS_frames_only")
    make_empty_sample(tmp_path / "Empty")
    return tmp_path


def test_particle_sizes_of_samples_without_spectra(psd_dir: Path):
    from autoemx.gui.common import get_analysis, get_particle_sizes

    ledgers = {p.parent.name: p.stat().st_mtime for p in psd_dir.glob("*/ledger.json")}
    measured = be.load_sample_info(str(psd_dir / "PSD_only"))
    assert (measured.has_eds, measured.has_sizes) == (False, True)
    data = be.load_particle_sizes(measured)
    assert data.stored and data.source == "measurement" and data.eds_d_min_um is None
    assert len(data.diameters_um) == data.result.size_stats.n_par_analysed

    eds = get_particle_sizes(str(psd_dir / "EDS_frames_only"))  # no statistics in the ledger: computed
    assert not eds.stored and eds.source == "eds_frames" and eds.eds_d_min_um == 3.0
    assert eds.result.size_stats.n_par_analysed == 60

    empty = be.load_sample_info(str(psd_dir / "Empty"))
    assert not empty.has_sizes and be.load_particle_sizes(empty) is None
    # EDS views have no data for these samples (no error)
    for name in ("PSD_only", "EDS_frames_only", "Empty"):
        assert get_analysis(str(psd_dir / name), None)[1] is None
    # Read-only: ledgers not rewritten
    assert {p.parent.name: p.stat().st_mtime for p in psd_dir.glob("*/ledger.json")} == ledgers

    for d in (data, eds, None):
        assert pl.psd_figure(d).data is not None


def test_particle_size_callbacks(psd_dir: Path, results_dir: Path):
    from autoemx.gui.app import create_app

    app = create_app("")
    client = app.server.test_client()

    def call(prefix, outputs, sample_dir, state=(), extra_inputs=()):
        key = next(k for k in app.callback_map if k.startswith(prefix))
        body = {"output": key, "outputs": [{"id": i, "property": p} for i, p in outputs],
                "inputs": [{"id": "sample-dd", "property": "value", "value": str(sample_dir)},
                           {"id": "data-version", "property": "data", "value": 0}, *extra_inputs],
                "state": list(state), "changedPropIds": ["sample-dd.value"]}
        res = client.post("/_dash-update-component", json=body)
        assert res.status_code == 200, res.data
        return res.get_json()["response"]

    eds_tabs = ["plot-tab", "clusters-tab", "spectra-tab", "dist-tab", "figures-tab", "settings-tab"]
    grey_outputs = [(t, "disabled") for t in eds_tabs + ["psd-tab"]] + [
        ("tabs", "value"), ("analysis-sidebar", "className"), ("spectrum-panel", "className"),
        ("analysis-dd", "placeholder"), ("set-active-btn", "disabled")]

    def grey(sample_dir, tab):
        return call("..plot-tab.disabled", grey_outputs, sample_dir,
                    [{"id": "tabs", "property": "value", "value": tab}])

    res = grey(psd_dir / "PSD_only", "plot")
    assert all(res[t]["disabled"] for t in eds_tabs) and not res["psd-tab"]["disabled"]
    assert res["tabs"]["value"] == "psd" and "no-eds" in res["analysis-sidebar"]["className"]
    assert res["set-active-btn"]["disabled"]
    res = grey(psd_dir / "Empty", "clusters")
    assert res["psd-tab"]["disabled"] and res["tabs"]["value"] == "sem"
    res = grey(results_dir / WULFENITE_MINI_ID, "psd")  # spectra, no particle sizes
    assert not any(res[t]["disabled"] for t in eds_tabs) and res["psd-tab"]["disabled"]
    assert res["tabs"]["value"] == "plot" and res["analysis-sidebar"]["className"] == "sidebar"

    psd_outputs = [("psd-banner", "children"), ("psd-graph", "figure"), ("psd-tables", "children"),
                   ("psd-compare-controls", "hidden"), ("psd-compare-samples", "hidden")]
    mode = [{"id": i, "property": "value", "value": v}
            for i, v in (("psd-mode", "sample"), ("psd-compare", []), ("psd-kind", "cumulative"),
                         ("psd-basis", "number"), ("psd-result-dd", None))]
    options = [{"id": "sample-dd", "property": "options", "value": []}]
    res = call("..psd-banner.children", psd_outputs, psd_dir / "EDS_frames_only", options, mode)
    assert "frames searched by EDS only" in str(res["psd-banner"]["children"])
    assert "not saved" in str(res["psd-banner"]["children"])
    res = call("..psd-banner.children", psd_outputs, psd_dir / "Empty", options, mode)
    assert "No particle size data" in str(res["psd-banner"]["children"])


@pytest.fixture()
def series_dir(tmp_path: Path) -> Path:
    """Results folder: a temperature series and an EDS-only sample, a sub-folder of two samples, a sample without
    particle sizes."""
    from _psd_helpers import make_eds_frames_sample, make_empty_sample, make_measurement_sample

    for i, (name, median) in enumerate((("S_900C", 1.5), ("S_1000C", 2.5), ("S_1100C", 4.0))):
        make_measurement_sample(tmp_path / name, seed=i, median_um=median, sigma=0.4)
    make_eds_frames_sample(tmp_path / "S_EDS")
    make_measurement_sample(tmp_path / "batch2" / "H_A", seed=10, median_um=0.8)
    make_measurement_sample(tmp_path / "batch2" / "H_B", seed=11, median_um=1.2)
    make_empty_sample(tmp_path / "No_sizes")
    return tmp_path


def test_size_comparison_table_and_figure(series_dir: Path):
    from autoemx.gui.common import get_particle_sizes, has_sizes

    assert has_sizes(str(series_dir / "S_EDS")) and not has_sizes(str(series_dir / "No_sizes"))
    assert not be.sample_has_sizes(str(series_dir / "missing"))
    items = [(name, get_particle_sizes(str(series_dir / name))) for name in ("S_900C", "S_1100C", "S_EDS")]
    df = be.size_comparison_table(items)
    assert list(df["Sample"]) == ["S_900C", "S_1100C", "S_EDS"]
    assert list(df["Source"]) == ["size measurement", "size measurement", "EDS frames"]
    assert list(df["Saved"]) == [True, True, False]
    assert df.loc[0, "D50"] == pytest.approx(items[0][1].result.size_stats.median)
    assert df.loc[0, "D50 CI low"] <= df.loc[0, "D50"] <= df.loc[0, "D50 CI high"]
    for kind in ("cumulative", "histogram"):
        for basis in ("number", "volume"):
            fig = pl.psd_comparison_figure(items, kind, basis)
            assert {t.name for t in fig.data if t.showlegend is not False} == {"S_900C", "S_1100C", "S_EDS"}
    assert pl.psd_comparison_figure([]).layout.annotations  # empty-state message


def test_size_comparison_callbacks(series_dir: Path):
    from autoemx.gui.app import create_app

    app = create_app("")
    client = app.server.test_client()
    samples = [{"label": s.label, "value": s.sample_dir} for s in be.find_samples(str(series_dir))]
    by_label = {o["label"]: o["value"] for o in samples}

    def post(prefix, outputs, inputs, state=(), changed=()):
        key = next(k for k in app.callback_map if k.startswith(prefix))
        outs = [{"id": i, "property": p} for i, p in outputs]
        body = {"output": key, "outputs": outs if len(outs) > 1 else outs[0],  # one output: not a list
                "inputs": [{"id": i, "property": p, "value": v} for i, p, v in inputs],
                "state": [{"id": i, "property": p, "value": v} for i, p, v in state],
                "changedPropIds": list(changed)}
        res = client.post("/_dash-update-component", json=body)
        assert res.status_code == 200, res.data
        return res.get_json()["response"]

    def choices(sample_label):
        res = post("..psd-compare.options", [("psd-compare", "options"), ("psd-compare", "value")],
                   [("sample-dd", "options", samples), ("sample-dd", "value", by_label[sample_label]),
                    ("data-version", "data", 0)], changed=["sample-dd.value"])
        return [o["label"] for o in res["psd-compare"]["options"]], res["psd-compare"]["value"]

    labels, value = choices("S_1000C")
    # Samples with particle sizes, in natural order; by default those of the folder of the shown sample
    assert labels == ["batch2/H_A", "batch2/H_B", "S_900C", "S_1000C", "S_1100C", "S_EDS"]
    assert sorted(value) == sorted(by_label[s] for s in ("S_900C", "S_1000C", "S_1100C", "S_EDS"))
    assert sorted(choices("batch2/H_A")[1]) == sorted([by_label["batch2/H_A"], by_label["batch2/H_B"]])

    outputs = [("psd-banner", "children"), ("psd-graph", "figure"), ("psd-tables", "children"),
               ("psd-compare-controls", "hidden"), ("psd-compare-samples", "hidden")]

    def view(mode, compared):
        return post("..psd-banner.children", outputs,
                    [("sample-dd", "value", by_label["S_1000C"]), ("data-version", "data", 0),
                     ("psd-mode", "value", mode), ("psd-compare", "value", compared),
                     ("psd-kind", "value", "cumulative"), ("psd-basis", "value", "number"),
                     ("psd-result-dd", "value", None)],
                    [("sample-dd", "options", samples)], ["psd-mode.value"])

    res = view("compare", value)
    assert res["psd-compare-controls"]["hidden"] is False
    assert "4 samples compared" in str(res["psd-banner"]["children"])
    assert "S_EDS" in str(res["psd-banner"]["children"])  # EDS frames only: flagged
    rows = res["psd-tables"]["children"][0]["props"]["children"][1]["props"]["data"]
    assert [r["Sample"] for r in rows] == ["S_900C", "S_1000C", "S_1100C", "S_EDS"]
    assert view("sample", value)["psd-compare-controls"]["hidden"] is True

    res = post("psd-download.data", [("psd-download", "data")], [("psd-export-btn", "n_clicks", 1)],
               [("psd-compare", "value", value), ("sample-dd", "options", samples), ("psd-mode", "value", "compare"),
                ("sample-dd", "value", by_label["S_1000C"])], ["psd-export-btn.n_clicks"])
    assert res["psd-download"]["data"]["filename"] == "Particle_size_comparison.csv"
    assert "D50 CI low" in res["psd-download"]["data"]["content"]


def test_size_frames_and_saving_exclusions(tmp_path: Path):
    from _psd_helpers import make_eds_frames_sample, make_measurement_sample

    from autoemx.config.ledger_schemas import SampleLedger
    from autoemx.gui.common import get_info, get_particle_sizes

    sample = str(make_measurement_sample(tmp_path / "M", images=True))
    data = get_particle_sizes(sample)
    frames = be.size_frames(get_info(sample), data)
    assert [f.frame_id for f in frames][:2] == ["B0F0", "B0F1"] and len(frames) == 12
    assert all(f.image_path and Path(f.image_path).exists() and f.status == "included" for f in frames)
    assert len(frames[0].marks) == frames[0].n_particles and frames[0].scan_label.startswith("Band 1")
    assert "psd" in be._LEDGER_WRITERS  # never runs together with an acquisition or quantification

    # Saved through the job function: frames and particles (frame, index) left out recorded in a new result
    key = ["B0F0", frames[0].marks[0][3]]
    assert be._run_particle_stats(sample, {"frame_ids": ["B0F1"], "particles": [key]}) == {"ok": True}
    ledger = SampleLedger.from_json_file(str(tmp_path / "M" / "ledger.json"))
    result = ledger.size_result
    assert result.excluded_frame_ids == ["B0F1"] and [list(k) for k in result.excluded_particles] == [key]
    assert result.size_stats.n_par_analysed == data.result.size_stats.n_par_analysed - frames[1].n_particles - 1
    assert len(ledger.size_results) == 2  # the statistics without exclusions are kept
    assert be.size_frames(get_info(sample), get_particle_sizes(sample))[1].status == "left out"
    assert not be._run_particle_stats(sample, {"frame_ids": [f.frame_id for f in frames]})["ok"]

    # EDS frames only: statistics computed, then saved
    eds = str(make_eds_frames_sample(tmp_path / "E", images=True))
    assert be.can_save_particle_sizes(get_particle_sizes(eds)) and not get_particle_sizes(eds).stored
    assert be._run_particle_stats(eds, {})["ok"]
    assert get_particle_sizes(eds).stored


def test_size_frames_callback(tmp_path: Path):
    from _psd_helpers import make_eds_frames_sample, make_measurement_sample

    from autoemx.gui.app import create_app

    sample = str(make_measurement_sample(tmp_path / "M", images=True))
    eds = str(make_eds_frames_sample(tmp_path / "E", images=True))
    app = create_app("")
    client = app.server.test_client()
    outputs = [("psd-frame-dd", "options"), ("psd-frame-dd", "value"), ("psd-frame-graph", "figure"),
               ("psd-frame-caption", "children"), ("psd-frame-out", "value"), ("psd-excl", "data"),
               ("psd-pending", "children"), ("psd-recalc-btn", "disabled"), ("psd-recalc-btn", "children"),
               ("psd-reset-btn", "disabled")]

    def view(sample_dir, changed, frame=None, out=(), click=None, excl=None):
        key = next(k for k in app.callback_map if k.startswith("..psd-frame-dd.options"))
        values = {"sample-dd.value": sample_dir, "data-version.data": 0, "psd-frame-dd.value": frame,
                  "psd-frame-prev.n_clicks": None, "psd-frame-next.n_clicks": 1 if changed == "psd-frame-next.n_clicks"
                  else None, "psd-frame-out.value": list(out), "psd-frame-graph.clickData": click,
                  "psd-reset-btn.n_clicks": None, "psd-result-dd.value": None}
        body = {"output": key, "outputs": [{"id": i, "property": p} for i, p in outputs],
                "inputs": [{"id": k.split(".")[0], "property": k.split(".")[1], "value": v} for k, v in values.items()],
                "state": [{"id": "psd-excl", "property": "data", "value": excl}], "changedPropIds": [changed]}
        res = client.post("/_dash-update-component", json=body)
        assert res.status_code == 200, res.data
        return {f"{i}.{p}": v for i, d in res.get_json()["response"].items() for p, v in d.items()}

    res = view(sample, "sample-dd.value")
    assert res["psd-frame-dd.value"] == "B0F0" and res["psd-excl.data"]["frames"] == []
    assert res["psd-recalc-btn.disabled"] is True  # nothing changed, statistics saved
    excl = res["psd-excl.data"]
    res = view(sample, "psd-frame-next.n_clicks", frame="B0F0", excl=excl)
    assert res["psd-frame-dd.value"] == "B0F1"
    res = view(sample, "psd-frame-out.value", frame="B0F1", out=["out"], excl=excl)
    assert res["psd-excl.data"]["frames"] == ["B0F1"] and res["psd-recalc-btn.disabled"] is False
    assert "not saved" in str(res["psd-pending.children"]) and res["psd-reset-btn.disabled"] is False
    assert "left out" in res["psd-frame-dd.options"][1]["label"]
    # Clicking a particle leaves it out (frame, index), clicking it again counts it again
    pid = res["psd-frame-graph.figure"]["data"][0]["customdata"][0]
    click = {"points": [{"customdata": pid}]}
    res = view(sample, "psd-frame-graph.clickData", frame="B0F1", out=["out"], click=click, excl=res["psd-excl.data"])
    assert res["psd-excl.data"]["particles"] == [["B0F1", pid]]
    res = view(sample, "psd-frame-graph.clickData", frame="B0F1", out=["out"], click=click, excl=res["psd-excl.data"])
    assert res["psd-excl.data"]["particles"] == []

    # EDS frames: statistics not saved yet; their particles can be left out one by one too
    res = view(eds, "sample-dd.value")
    assert res["psd-recalc-btn.disabled"] is False and res["psd-recalc-btn.children"] == "Save statistics"
    assert {t.get("name") for t in res["psd-frame-graph.figure"]["data"]} == {"Counted (frame searched by EDS)"}
    res = view(eds, "psd-frame-graph.clickData", frame="EDS-A0", click={"points": [{"customdata": 2}]},
               excl=res["psd-excl.data"])
    assert res["psd-excl.data"]["particles"] == [["EDS-A0", 2]]
    res = view(eds, "psd-frame-graph.clickData", frame="EDS-A0", click={"points": [{"customdata": None}]},
               excl=res["psd-excl.data"])
    assert res["psd-excl.data"]["particles"] == [["EDS-A0", 2]]  # no particle clicked


def test_acquisition_settings_of_particle_size_statistics():
    from autoemx.config.runtime_configs import PowderMeasurementConfig
    from autoemx.gui.tab_acquisition import active_acq_params
    from autoemx.runners.batch_acquire_and_analyze import batch_acquire_and_analyze
    from autoemx.runners.collect_particle_statistics import collect_particle_statistics

    defaults = {s.key: s.default for s in be.ACQ_SPECS}
    values = be.coerce_values(defaults, be.ACQ_SPECS)
    # Settings of size statistics in their own section, passed with the powder settings (valid for the config)
    assert {s.name for s in be.ACQ_SPECS if s.section == "asize"} == {
        "n_par_target", "min_par_diameter_um", "max_par_diameter_um", *be.ACQ_SIZE_FIELDS, *be.ACQ_PSD_POWDER_FIELDS}
    # ... right below the acquisition settings
    sections = [key for key, _ in be.ACQ_SECTIONS]
    assert sections.index("asize") == sections.index("aacq") + 1
    assert be.SPECS_BY_KEY["asize.par_segmentation_model"].choices == be.SPECS_BY_KEY[
        "apowder.par_segmentation_model"].choices
    assert values["asize.size_stats_ci_quantities"] == ["D50"]
    values["asize.size_stats_ci_quantities"] = ["D90", "D10", "Dx"]  # checkboxes: choices only, in their order
    values = be.coerce_values(values, be.ACQ_SPECS)
    assert values["asize.size_stats_ci_quantities"] == ["D10", "D90"]
    eds = be.acquisition_kwargs(values)
    PowderMeasurementConfig(**eds["powder_meas_cfg_kwargs"])
    assert eds["powder_meas_cfg_kwargs"]["size_stats_ci_quantities"] == ["D10", "D90"]
    assert "size_stats_n_par_target" not in eds and "n_par_target" not in eds["powder_meas_cfg_kwargs"]
    assert eds["powder_meas_cfg_kwargs"]["max_area_par"] == values["apowder.max_area_par"]  # EDS: powder settings

    # EDS + particle sizes
    combined = be.acquisition_kwargs({**values, "asize.n_par_target": 300, "asize.size_stats_ci_tolerance": 5.0},
                                     "eds+sizes")
    assert set(combined) <= set(inspect.signature(batch_acquire_and_analyze).parameters)
    assert combined["size_stats_n_par_target"] == 300
    assert combined["powder_meas_cfg_kwargs"]["size_stats_ci_tolerance"] == 0.05  # entered in %
    # Particles measured by the size statistics: their own settings, diameters as particle areas
    size_cfg = combined["size_stats_powder_meas_cfg_kwargs"]
    assert size_cfg["min_area_par"] == pytest.approx(3.14159 * 0.15**2, rel=1e-4)
    assert size_cfg["max_area_par"] == pytest.approx(3.14159 * 10**2, rel=1e-4)
    assert size_cfg["par_segmentation_model"] == values["asize.par_segmentation_model"]
    PowderMeasurementConfig(**{**combined["powder_meas_cfg_kwargs"], **size_cfg})
    with pytest.raises(ValueError, match="diameters"):
        be.acquisition_kwargs({**values, "asize.min_par_diameter_um": 30.0}, "eds+sizes")
    with pytest.raises(ValueError, match="automatic navigation.*automatic particle selection"):
        be.acquisition_kwargs({**values, "aacq.is_manual_navigation": True,
                               "apowder.par_selection_mode": "manual"}, "eds+sizes")

    # Particle sizes only
    psd = be.psd_acquisition_kwargs({**values, "asize.n_par_target": 150, "apowder.max_area_par": 1.0})
    assert set(psd) <= set(inspect.signature(collect_particle_statistics).parameters)
    assert psd["n_par_target"] == 150 and "contrast" not in psd  # automatic brightness/contrast
    cfg = PowderMeasurementConfig(**psd["powder_meas_cfg_kwargs"])
    # Particle settings of the size statistics section (particles of 0.3–20 µm by default)
    assert cfg.par_selection_mode == "auto" and cfg.max_area_par == pytest.approx(3.14159 * 10**2, rel=1e-4)
    assert cfg.min_area_par == pytest.approx(3.14159 * 0.15**2, rel=1e-3)
    assert cfg.min_par_pixels == values["asize.min_par_pixels"] and "save_raw_images" not in psd
    assert "max_spectra_per_par" not in psd["powder_meas_cfg_kwargs"]  # EDS settings left out
    manual_bc = be.psd_acquisition_kwargs({**values, "aacq.auto_adjust_brightness_contrast": False,
                                           "aacq.contrast": 4.0, "aacq.brightness": 0.4})
    assert (manual_bc["contrast"], manual_bc["brightness"]) == (4.0, 0.4)

    samples = be.acquisition_samples([{"ID": "P1", "els": "", "x": 1, "y": 2}], require_elements=False)
    assert samples[0]["els"] == []
    with pytest.raises(ValueError, match="no elements"):
        be.acquisition_samples([{"ID": "P1", "els": "", "x": 1, "y": 2}])
    for mode, kwargs, runner in (("sizes", psd, "collect_particle_statistics("),
                                 ("eds+sizes", combined, "comp_analyzer = batch_acquire_and_analyze(")):
        script = be.acquisition_script([{"ID": "P1", "pos": (1.0, 2.0)}], kwargs, "/data/results", mode)
        compile(script, "script.py", "exec")
        assert runner in script and f"from autoemx.runners.{runner.split('= ')[-1][:-1]} import" in script

    # Fields shown for each measurement
    active = active_acq_params(defaults, "eds")
    assert active["asize.n_par_target"] == "hidden" and active.get("asize.min_par_pixels", "on") == "on"
    assert active["asize.max_par_diameter_um"] == "hidden"  # EDS: the powder settings are used
    assert active["aimg.save_raw_images"] == "off"  # raw copies of annotated images only
    assert active_acq_params({**defaults, "aimg.annotate_particle_images": ["on"]})["aimg.save_raw_images"] == "on"
    active = active_acq_params(defaults, "eds+sizes")
    assert active.get("asize.n_par_target", "on") == "on" and active["asize.size_stats_max_n_frames"] == "off"
    # Same size settings as for particle sizes only
    assert all(active.get(f"asize.{s.name}", "on") == "on" for s in be.ACQ_SPECS if s.section == "asize"
               and not s.name.startswith("size_stats_m") and s.name != "size_stats_ci_quantities")
    active = active_acq_params({**defaults, "asize.size_stats_ci_tolerance": 0.05}, "sizes")
    assert active.get("asize.size_stats_max_n_frames", "on") == "on"
    # Particle sizes only: particle settings in the size statistics section, no EDS settings
    assert (active["aacq.beam_energy"], active["aquant"], active["abulk"], active["apowder"]) == (
        "hidden", "hidden", "hidden", "hidden")
    assert active.get("asize.max_par_diameter_um", "on") == "on" and active["aimg.save_raw_images"] == "hidden"
    # Brightness threshold: only for the threshold segmentation model
    other_models = [m for m in be.SPECS_BY_KEY["apowder.par_segmentation_model"].choices if m != "threshold_bright"]
    assert active_acq_params(defaults)["apowder.par_brightness_thresh"] == "on"
    for model in other_models:
        for section in ("apowder", "asize"):
            assert active_acq_params({**defaults, f"{section}.par_segmentation_model": model}, "sizes")[
                f"{section}.par_brightness_thresh"] == "off"

    # Progress: spectra, then particles of the size band being measured
    log = ("Sample 'A'\n🔬 Acquiring spectrum #0...\nSample 'A': particle size statistics\nSize band 1: 2–40 µm\n"
           "A total of 12 particles have now been analyzed in this size band. 38 more to go.\n"
           "Size band 2: 0.4–2 µm\nA total of 1 particle has now been analyzed in this size band. 49 more to go.")
    [p] = be.acquisition_progress(log, ["A"], 100)
    assert p["done"] == 1 and p["sizes"] == {"band": 2, "done": 1, "total": 50}
    log = "Sample 'A'\nERROR: Sample 'A': particle stats collection failed: no frame\nSample 'B'\nSize band 1: 2–40 µm"
    progress = be.acquisition_progress(log, ["A", "B"], 0)
    assert [(p["state"], p["error"]) for p in progress] == [("failed", "no frame"), ("running", None)]
    assert progress[1]["sizes"] == {"band": 1, "done": 0, "total": None}


def test_size_results_segmentations_and_active(tmp_path: Path):
    from _psd_helpers import make_eds_frames_sample, make_measurement_sample

    from autoemx.config.ledger_schemas import SampleLedger
    from autoemx.gui.common import get_info, get_particle_sizes

    sample = str(make_measurement_sample(tmp_path / "M", images=True))
    info = get_info(sample)
    n_particles = info.ledger.size_result.size_stats.n_par_analysed
    assert be.size_result_options(info) == [
        {"label": f"Result 0 · measurement 0 · threshold_bright · 100 · 80 px · {n_particles} particles · active",
         "value": 0}]

    # Segmented again (another threshold): a new segmentation and result, made the active one; the first kept
    assert be._run_resegmentation(sample, {"par_brightness_thresh": 200, "par_segmentation_model": "threshold_bright",
                                           "min_par_pixels": 80}) == {"ok": True}
    ledger = SampleLedger.from_json_file(str(tmp_path / "M" / "ledger.json"))
    assert [s.par_brightness_thresh for s in ledger.segmentations] == [100, 200]
    assert ledger.active_size_result == 1 and ledger.size_result.segmentation_id == 1
    info = get_info(sample)
    labels = [o["label"] for o in be.size_result_options(info)]
    assert labels[0].startswith("Result 1 · measurement 0 · threshold_bright · 200") and labels[0].endswith("active")
    first = get_particle_sizes(sample, 0)
    assert first.result.result_id == 0 and first.result.segmentation_id == 0 and first.stored
    assert get_particle_sizes(sample).result.result_id == 1  # default: the active one
    assert not be.mixed_segmentations(first)

    # Recalculating a result that is not the active one (frame left out): its segmentation is kept
    assert be._run_particle_stats(sample, {"frame_ids": ["B0F0"], "result_id": 0}) == {"ok": True}
    ledger = SampleLedger.from_json_file(str(tmp_path / "M" / "ledger.json"))
    assert (ledger.size_result.segmentation_id, ledger.size_result.excluded_frame_ids) == (0, ["B0F0"])

    assert be.set_active_size_result(get_info(sample), 1) == "Result 1 is now the active one."
    assert SampleLedger.from_json_file(str(tmp_path / "M" / "ledger.json")).active_size_result == 1
    with pytest.raises(ValueError, match="No result 9"):
        be.set_active_size_result(get_info(sample), 9)

    # Sample measured only by EDS: its frames segmented again
    eds = str(make_eds_frames_sample(tmp_path / "E", images=True))
    assert be._run_resegmentation(eds, {"par_brightness_thresh": 150})["ok"]
    result = get_particle_sizes(eds).result
    assert (result.measurement_id, result.eds_segmentation_id) == (None, 1) and get_particle_sizes(eds).stored


def test_compare_segmentations_and_result_choices(tmp_path: Path):
    from _psd_helpers import make_measurement_sample

    from autoemx.gui.app import create_app

    sample = str(make_measurement_sample(tmp_path / "M", images=True))
    be._run_resegmentation(sample, {"par_brightness_thresh": 200})
    app = create_app("")
    client = app.server.test_client()

    def post(prefix, outputs, inputs, state=(), changed=()):
        key = next(k for k in app.callback_map if k.startswith(prefix))
        outs = [{"id": i, "property": p} for i, p in outputs]
        body = {"output": key, "outputs": outs if len(outs) > 1 else outs[0],
                "inputs": [{"id": i, "property": p, "value": v} for i, p, v in inputs],
                "state": [{"id": i, "property": p, "value": v} for i, p, v in state], "changedPropIds": list(changed)}
        res = client.post("/_dash-update-component", json=body)
        assert res.status_code == 200, res.data
        return res.get_json()["response"]

    res = post("..psd-result-dd.options", [("psd-result-dd", "options"), ("psd-result-dd", "value"),
                                           ("psd-seg-model", "value"), ("psd-seg-thresh", "value"),
                                           ("psd-seg-px", "value"), ("psd-seg-params", "value")],
               [("sample-dd", "value", sample), ("data-version", "data", 0)], changed=["sample-dd.value"])
    assert [o["value"] for o in res["psd-result-dd"]["options"]] == [1, 0] and res["psd-result-dd"]["value"] == 1
    assert (res["psd-seg-model"]["value"], res["psd-seg-thresh"]["value"]) == ("threshold_bright", 200)
    assert res["psd-seg-params"]["value"] == ""  # the threshold model has no options

    outputs = [("psd-banner", "children"), ("psd-graph", "figure"), ("psd-tables", "children"),
               ("psd-compare-controls", "hidden"), ("psd-compare-samples", "hidden")]
    res = post("..psd-banner.children", outputs,
               [("sample-dd", "value", sample), ("data-version", "data", 0), ("psd-mode", "value", "results"),
                ("psd-compare", "value", []), ("psd-kind", "value", "cumulative"), ("psd-basis", "value", "number"),
                ("psd-result-dd", "value", None)], [("sample-dd", "options", [])], ["psd-mode.value"])
    assert "2 results of this sample compared" in str(res["psd-banner"]["children"])
    rows = res["psd-tables"]["children"][0]["props"]["children"][1]["props"]["data"]
    assert [r["Result"].split(" · ")[0] for r in rows] == ["Result 0", "Result 1"]
    assert res["psd-compare-samples"]["hidden"] is True  # the samples list does not apply
    assert len([t for t in res["psd-graph"]["figure"]["data"] if t.get("showlegend") is not False]) == 2

    # Shown result: its segmentation named in the banner
    res = post("..psd-banner.children", outputs,
               [("sample-dd", "value", sample), ("data-version", "data", 0), ("psd-mode", "value", "sample"),
                ("psd-compare", "value", []), ("psd-kind", "value", "cumulative"), ("psd-basis", "value", "number"),
                ("psd-result-dd", "value", 0)], [("sample-dd", "options", [])], ["psd-result-dd.value"])
    assert "size frames threshold_bright · 100 · 80 px" in str(res["psd-banner"]["children"])


def test_segmentation_model_not_installed(tmp_path: Path, monkeypatch):
    """A model not installed is listed as such; selecting it shows how to install it, and nothing runs with it."""
    from _psd_helpers import make_measurement_sample

    from autoemx.core.em_runtime.particle_segmentation_models import SAM2
    from autoemx.gui.app import create_app
    from autoemx.gui.common import _pid, segmentation_note_id

    monkeypatch.setattr(SAM2, "is_available", lambda: False)
    assert "SAM2" in be.par_segmentation_models()
    assert be.segmentation_model_label("SAM2") == "SAM2 (not installed)"
    assert be.segmentation_model_label("threshold_bright") == "threshold_bright"
    assert 'pip install "autoemx[segmentation]"' in be.segmentation_install_note("SAM2")
    assert be.segmentation_install_note("threshold_bright") is None
    assert be._validate_models({"asize.par_segmentation_model": "SAM2"}) == [be.segmentation_install_note("SAM2")]

    sample = str(make_measurement_sample(tmp_path / "M", images=True))
    app = create_app("")
    client = app.server.test_client()

    def post(output, inputs, state=()):
        key = next(k for k in app.callback_map if k.startswith(output))
        outs = [{"id": i, "property": p} for i, p in [(o.rsplit(".", 1)[0], o.rsplit(".", 1)[1])
                                                     for o in key.removeprefix("..").removesuffix("..").split("...")]]
        body = {"output": key, "outputs": outs if len(outs) > 1 else outs[0],
                "inputs": [{"id": i, "property": p, "value": v} for i, p, v in inputs],
                "state": [{"id": i, "property": p, "value": v} for i, p, v in state],
                "changedPropIds": [f"{inputs[0][0]}.{inputs[0][1]}"]}
        res = client.post("/_dash-update-component", json=body)
        assert res.status_code == 200, res.data
        return res.get_json()["response"]

    # Segment again: the note only for a model not installed, and no job started with it
    res = post("..psd-seg-thresh-box.hidden", [("psd-seg-model", "value", "SAM2")])
    assert "not installed" in res["psd-seg-note"]["children"]
    # Fields of the selected model only: options for SAM2, the threshold for the threshold model
    assert (res["psd-seg-thresh-box"]["hidden"], res["psd-seg-params-box"]["hidden"]) == (True, False)
    res = post("..psd-seg-thresh-box.hidden", [("psd-seg-model", "value", "threshold_bright")])
    assert res["psd-seg-note"]["children"] is None
    assert (res["psd-seg-thresh-box"]["hidden"], res["psd-seg-params-box"]["hidden"]) == (False, True)
    res = post("..psd-job-store.data", [("psd-reseg-btn", "n_clicks", 1)],
               [("sample-dd", "value", sample), ("psd-seg-model", "value", "SAM2"), ("psd-seg-thresh", "value", 100),
                ("psd-seg-px", "value", 80), ("psd-seg-params", "value", "")])
    assert "not installed" in str(res["psd-recalc-msg"]["children"]) and "psd-job-store" not in res

    # Acquisition: the same note below the model setting
    note_id, model_id = segmentation_note_id("asize.par_segmentation_model"), _pid("asize.par_segmentation_model")
    key = next(k for k in app.callback_map if k.startswith("{") and '"seg-note"' in k and "asize" in k)
    body = {"output": key, "outputs": {"id": note_id, "property": "children"},
            "inputs": [{"id": model_id, "property": "value", "value": "SAM2"}], "changedPropIds": []}
    res = client.post("/_dash-update-component", json=body).get_json()["response"]
    assert "not installed" in str(res)
