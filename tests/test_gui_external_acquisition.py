#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for launching acquisitions from an external program: Acquisition tab of the GUI opened with samples
given in its URL, and report of the run (no microscope: the composition analyzer is replaced by a fake one).
"""

from __future__ import annotations

import json
import threading
import time
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlsplit

import pytest

pytest.importorskip("dash")

from autoemx.gui import backend as be
from autoemx.gui import run_report as rr
from autoemx.gui import tab_acquisition as ta
from autoemx.gui.__main__ import _read_samples_file, acquisition_url
from autoemx.gui.app import create_app

SAMPLES = [
    {"ID": "LaNbO4_A", "els": ["La", "Nb", "O"], "pos": (1.5, -2), "cnd": ["LaNbO4"]},
    {"ID": "Hematite", "els": "Fe, O", "pos": [3, 4]},
]
SAVED_ROW = {"uid": "saved1", "ID": "Saved", "els": "Si, O", "x": 0, "y": 0, "cnd": ""}


def _query(url: str) -> str:
    return "?" + urlsplit(url).query


def _rows(search):
    return [{k: v for k, v in r.items() if k != "uid"} for r in ta.prefill_from_query(search)["rows"]]


# ----------------------------------------------------------------------------- URL
def test_acquisition_url_round_trip(tmp_path: Path):
    url = acquisition_url(SAMPLES, str(tmp_path), port=8123)
    assert url.startswith("http://127.0.0.1:8123/?acq=")
    data = ta.prefill_from_query(_query(url))
    assert data["folder"] == str(tmp_path)
    assert _rows(_query(url)) == [
        {"ID": "LaNbO4_A", "els": "La, Nb, O", "x": 1.5, "y": -2, "cnd": "LaNbO4"},
        {"ID": "Hematite", "els": "Fe, O", "x": 3, "y": 4, "cnd": ""},
    ]
    assert len({r["uid"] for r in data["rows"]}) == 2


def test_acquisition_url_without_prefill():
    assert acquisition_url(port=8123) == "http://127.0.0.1:8123/"
    assert ta.prefill_from_query(None) is None
    assert ta.prefill_from_query("") is None
    assert ta.prefill_from_query("?other=1") is None


def test_prefill_folder_or_samples_only(tmp_path: Path):
    only_folder = ta.prefill_from_query(_query(acquisition_url(None, str(tmp_path))))
    assert only_folder == {"folder": str(tmp_path), "rows": [], "run_id": None}
    only_samples = ta.prefill_from_query(_query(acquisition_url(SAMPLES[:1])))
    assert only_samples["folder"] is None and len(only_samples["rows"]) == 1


@pytest.mark.parametrize("acq, error", [
    ("{bad", "not valid JSON"),
    ("[]", "JSON object"),
    ('{"folder": 3}', "'folder'"),
    ('{"samples": {"ID": "A"}}', "'samples' must be a list"),
    ('{"samples": ["A"]}', "each sample must be a dict"),
    ('{"samples": [{"ID": "A", "pos": [1]}]}', "'pos' must be [x, y]"),
])
def test_prefill_invalid(acq, error):
    from urllib.parse import urlencode

    with pytest.raises(ValueError, match=error.replace("[", r"\[").replace("]", r"\]")):
        ta.prefill_from_query("?" + urlencode({"acq": acq}))


# ----------------------------------------------------------------------------- samples file
def test_read_samples_file(tmp_path: Path):
    as_list = tmp_path / "list.json"
    as_list.write_text(json.dumps(SAMPLES), encoding="utf-8")
    assert _read_samples_file(str(as_list)) == (json.loads(json.dumps(SAMPLES)), None, None)

    as_dict = tmp_path / "dict.json"
    as_dict.write_text(json.dumps({"folder": "/res", "samples": SAMPLES}), encoding="utf-8")
    assert _read_samples_file(str(as_dict))[1] == "/res"

    bad = tmp_path / "bad.json"
    bad.write_text('{"oops": 1}', encoding="utf-8")
    with pytest.raises(ValueError):
        _read_samples_file(str(bad))


# ----------------------------------------------------------------------------- callbacks
@pytest.fixture(scope="module")
def app():
    return create_app("")


@pytest.fixture(scope="module")
def client(app):
    return app.server.test_client()


def _cell_values(children) -> list:
    """Values of the cells of the samples list rendered by render_rows."""
    out = []
    for row in children[1:]:  # first row: header
        out.append([c["props"].get("value") for c in row["props"]["children"]
                    if c["props"]["id"].get("type") == "a-cell"])
    return out


def _render_rows(client, search, cells=(), saved=None, changed=("url.search",)):
    body = {
        "output": "a-rows.children",
        "outputs": {"id": "a-rows", "property": "children"},
        "inputs": [
            {"id": "a-add-row", "property": "n_clicks", "value": None},
            [], [],
            {"id": "url", "property": "search", "value": search},
        ],
        "state": [list(cells), {"id": "a-saved-rows", "property": "data", "value": saved}],
        "changedPropIds": list(changed),
    }
    res = client.post("/_dash-update-component", json=body)
    if res.status_code == 204:
        return None
    return res.get_json()["response"].get("a-rows", {}).get("children")  # None: rows not updated


def _cells(row):
    return [{"id": {"type": "a-cell", "uid": row["uid"], "col": col}, "property": "value", "value": row.get(col)}
            for col in ("ID", "els", "x", "y", "cnd")]


def test_render_rows_uses_samples_of_the_url(client):
    search = _query(acquisition_url(SAMPLES))
    for changed in ((), ("url.search",)):  # page load, or URL read just after
        assert _cell_values(_render_rows(client, search, saved=[SAVED_ROW], changed=changed)) == [
            ["LaNbO4_A", "La, Nb, O", 1.5, -2, "LaNbO4"],
            ["Hematite", "Fe, O", 3, 4, ""],
        ]


def test_render_rows_without_url_samples(client):
    saved = [SAVED_ROW["ID"], SAVED_ROW["els"], SAVED_ROW["x"], SAVED_ROW["y"], SAVED_ROW["cnd"]]
    # Page load: rows saved in the browser, or the example
    assert _cell_values(_render_rows(client, "", saved=[SAVED_ROW], changed=())) == [saved]
    assert _cell_values(_render_rows(client, "", changed=("url.search",)))[0][0] == ta._EXAMPLE_ROW["ID"]
    # URL read once the rows are drawn (or invalid): rows kept
    assert _render_rows(client, "", cells=_cells(SAVED_ROW), saved=[SAVED_ROW]) is None
    assert _render_rows(client, "?acq=%7Bbad", cells=_cells(SAVED_ROW), saved=[SAVED_ROW]) is None


def _prefill(app, client, search, stored_link=None):
    # Outputs with allow_duplicate carry a hash in the callback key: take it from the callback map
    key = next(k for k in app.callback_map if k.startswith("..a-prefill-msg.children"))
    body = {
        "output": key,
        "outputs": [{"id": "a-prefill-msg", "property": "children"}, {"id": "folder", "property": "value"},
                    {"id": "main-tabs", "property": "value"}, {"id": "a-link", "property": "data"}],
        "inputs": [{"id": "url", "property": "search", "value": search}],
        "state": [{"id": "a-link", "property": "data", "value": stored_link}],
        "changedPropIds": ["url.search"],
    }
    res = client.post("/_dash-update-component", json=body)
    return None if res.status_code == 204 else res.get_json()["response"] or None  # {}: nothing updated


def test_prefill_callback(app, client, tmp_path: Path):
    res = _prefill(app, client, _query(acquisition_url(SAMPLES, str(tmp_path))))
    assert res["folder"]["value"] == str(tmp_path)
    assert res["main-tabs"]["value"] == "acq"
    assert "2 samples loaded" in json.dumps(res["a-prefill-msg"])
    assert "a-link" not in res or res["a-link"]["data"] is None  # no run ID: runs are not reported

    res = _prefill(app, client, _query(acquisition_url(SAMPLES, str(tmp_path), run_id="run-1")))
    assert res["a-link"]["data"] == {"run_id": "run-1", "requested": ["LaNbO4_A", "Hematite"], "folder": str(tmp_path)}
    assert "run-1" in json.dumps(res["a-prefill-msg"])

    # Page reloaded (no link in the URL any more): the stored link gives the folder back, still reported
    res = _prefill(app, client, "", stored_link={"run_id": "run-1", "requested": ["A"], "folder": str(tmp_path)})
    assert res["folder"]["value"] == str(tmp_path) and "a-link" not in res and "main-tabs" not in res
    assert "run-1" in json.dumps(res["a-prefill-msg"])

    res = _prefill(app, client, "?acq=%7Bbad")
    assert res["main-tabs"]["value"] == "acq" and "folder" not in res
    assert "could not be read" in json.dumps(res["a-prefill-msg"])

    assert _prefill(app, client, "") is None


# ----------------------------------------------------------------------------- run report
@pytest.fixture()
def reports_dir(tmp_path: Path, monkeypatch) -> Path:
    path = tmp_path / "reports"
    monkeypatch.setattr(rr, "REPORTS_DIR", path)
    return path


def _ledger(sample_dir: Path, n: int) -> None:
    sample_dir.mkdir(parents=True, exist_ok=True)
    (sample_dir / rr.LEDGER_NAME).write_text(json.dumps({"spectra": [{}] * n}), encoding="utf-8")


def _states(report):
    return [(s["ID"], s["status"]) for s in report["samples"]]


def test_run_id_checked(reports_dir: Path):
    for bad in ("../x", "a/b", "", "a b", 3, "x" * 101):
        with pytest.raises(ValueError):
            rr.clean_run_id(bad)
        with pytest.raises(ValueError):
            acquisition_url(SAMPLES, run_id=bad)
    assert ta.prefill_from_query(_query(acquisition_url(SAMPLES, run_id="Run_1-a")))["run_id"] == "Run_1-a"
    with pytest.raises(ValueError, match="run_id"):
        ta.prefill_from_query("?acq=" + json.dumps({"run_id": "../x"}))


def test_report_of_a_run(tmp_path: Path, reports_dir: Path):
    results = tmp_path / "res"
    _ledger(results / "A", 4)  # spectra from an earlier run
    report = rr.new_report("r1", ["A", "B", "C"], str(results), requested=["A", "B"])
    assert report["status"] == "running" and _states(report) == [("A", "running"), ("B", "waiting"), ("C", "waiting")]
    assert report["samples"][0]["path"] == str(results / "A") and report["samples"][0]["n_spectra_before"] == 4

    _ledger(results / "A", 10)
    rr.sample_done(report, "A", None)
    rr.sample_done(report, "B", RuntimeError("no particles found"))
    assert _states(report) == [("A", "done"), ("B", "failed"), ("C", "running")]
    assert report["samples"][0]["n_new_spectra"] == 6
    assert report["samples"][1]["error"] == "RuntimeError: no particles found"

    rr.write_report(report)
    assert rr.acquisition_report("r1") == report
    assert not list(reports_dir.glob("*.tmp"))
    rr.finish_report(report, "stopped", "Acquisition stopped")
    assert _states(report) == [("A", "done"), ("B", "failed"), ("C", "stopped")]
    assert report["samples"][2]["error"] == "Acquisition stopped" and report["finished"]


def test_report_of_a_failed_run():
    report = rr.new_report("r2", ["A", "B"], "/res")
    rr.finish_report(report, "failed", "EMError: Instrument driver could not be loaded")
    assert _states(report) == [("A", "failed"), ("B", "not_run")]
    assert {s["error"] for s in report["samples"]} == {"EMError: Instrument driver could not be loaded"}


def test_wait_for_acquisition(reports_dir: Path):
    assert rr.acquisition_report("r3") is None
    with pytest.raises(TimeoutError, match="not started"):
        rr.wait_for_acquisition("r3", timeout=0.2, poll_interval=0.05)
    report = rr.new_report("r3", ["A"], "/res")
    rr.write_report(report)
    with pytest.raises(TimeoutError, match="still running"):
        rr.wait_for_acquisition("r3", timeout=0.2, poll_interval=0.05)

    def finish():
        time.sleep(0.3)
        rr.sample_done(report, "A", None)
        rr.finish_report(report, "finished")
        rr.write_report(report)

    threading.Thread(target=finish).start()
    final = rr.wait_for_acquisition("r3", timeout=10, poll_interval=0.05)
    assert final["status"] == "finished" and _states(final) == [("A", "done")]
    assert not rr.report_path("r3").exists()  # deleted once read


class _FakeAnalyzer:
    """Composition analyzer without microscope: writes a ledger with 3 more spectra, fails for IDs with 'fail'."""

    def __init__(self, sample_id, results_dir=None, **_):
        self.sample_dir, self.EM_controller = Path(results_dir) / sample_id, None

    def run_collection_and_quantification(self, **_):
        if "fail" in self.sample_dir.name:
            raise RuntimeError("no particles found")
        _ledger(self.sample_dir, rr.count_spectra(str(self.sample_dir)) + 3)


def test_acquisition_job_reports_each_sample(tmp_path: Path, reports_dir: Path, monkeypatch):
    import autoemx.runners.batch_acquire_and_analyze as runner

    monkeypatch.setattr(runner, "EMXSp_Composition_Analyzer", _FakeAnalyzer)
    results = str(tmp_path / "res")
    samples = [{"ID": sid, "els": ["Si", "O"], "pos": (0, 0), "cnd": []} for sid in ("A", "B_fail", "C")]
    kwargs = be.acquisition_kwargs(be.coerce_values({s.key: s.default for s in be.ACQ_SPECS}, be.ACQ_SPECS))
    report = rr.new_report("r4", [s["ID"] for s in samples], results)
    seen = []
    monkeypatch.setattr(rr, "write_report", lambda rep: seen.append(_states(rep)))
    be._run_acquisition("", {"samples": samples, "kwargs": kwargs, "results_dir": results, "report": report})
    assert seen == [  # after each sample, then final
        [("A", "done"), ("B_fail", "running"), ("C", "waiting")],
        [("A", "done"), ("B_fail", "failed"), ("C", "running")],
        [("A", "done"), ("B_fail", "failed"), ("C", "done")],
        [("A", "done"), ("B_fail", "failed"), ("C", "done")],
    ]
    assert report["status"] == "finished" and [s["n_new_spectra"] for s in report["samples"]] == [3, 0, 3]


def _job(cancelled=False, result=None):
    return SimpleNamespace(process=SimpleNamespace(join=lambda: None), cancelled=cancelled, result=lambda: result)


def test_follow_job_makes_report_final(reports_dir: Path):
    # Stopped: the sample being acquired is 'stopped', the next ones 'not_run'
    report = rr.new_report("r5", ["A", "B"], "/res")
    rr.write_report(report)
    rr.follow_job(_job(cancelled=True, result={"ok": False, "error": "Cancelled."}), report)
    assert rr.acquisition_report("r5")["status"] == "stopped"
    assert _states(rr.acquisition_report("r5")) == [("A", "stopped"), ("B", "not_run")]

    # Process died
    report = rr.new_report("r6", ["A"], "/res")
    rr.write_report(report)
    rr.follow_job(_job(result={"ok": False, "error": "Process exited unexpectedly (code -9)."}), report)
    assert rr.acquisition_report("r6")["error"] == "Process exited unexpectedly (code -9)."

    # Made final by the job process: kept as is
    report = rr.new_report("r7", ["A"], "/res")
    final = json.loads(json.dumps(report))
    rr.sample_done(final, "A", None)
    rr.finish_report(final, "finished")
    rr.write_report(final)
    rr.follow_job(_job(result={"ok": True}), report)
    assert rr.acquisition_report("r7") == final

    # Report left by an earlier run with the same ID: replaced
    old = rr.new_report("r8", ["OLD"], "/res")
    old["started"] = "2000-01-01T00:00:00"
    rr.write_report(old)
    report = rr.new_report("r8", ["A"], "/res")
    rr.follow_job(_job(result={"ok": False, "error": "boom"}), report)
    assert _states(rr.acquisition_report("r8")) == [("A", "failed")]


def test_runner_on_sample_done(tmp_path: Path, monkeypatch):
    import autoemx.runners.batch_acquire_and_analyze as runner

    monkeypatch.setattr(runner, "EMXSp_Composition_Analyzer", _FakeAnalyzer)
    calls = []

    def on_done(sid, error):
        calls.append((sid, repr(error)))
        raise KeyError("errors of the callback do not stop the batch")

    samples = [{"ID": sid, "els": ["Si", "O"], "pos": (0, 0)} for sid in ("A", "B_fail", "C")]
    runner.batch_acquire_and_analyze(samples, development_mode=True, verbose=False, results_dir=str(tmp_path),
                                     on_sample_done=on_done)
    assert calls == [("A", "None"), ("B_fail", "RuntimeError('no particles found')"), ("C", "None")]


def test_particle_size_statistics_job_reports_each_sample(tmp_path: Path, reports_dir: Path, monkeypatch):
    import autoemx.runners.collect_particle_statistics as psd_runner

    measured = []

    class _FakePsdAnalyzer:
        def __init__(self, sample_id, **_):
            def stats(n):
                if "fail" in sample_id:
                    raise RuntimeError("no particles found")
                measured.append((sample_id, n))
            self.EM_controller = SimpleNamespace(particle_finder=SimpleNamespace(get_particle_stats=stats))

    monkeypatch.setattr(psd_runner, "EMXSp_Composition_Analyzer", _FakePsdAnalyzer)
    results = str(tmp_path / "res")
    values = be.coerce_values({s.key: s.default for s in be.ACQ_SPECS}, be.ACQ_SPECS)
    samples = [{"ID": sid, "pos": (0, 0)} for sid in ("A", "B_fail")]
    report = rr.new_report("r9", [s["ID"] for s in samples], results)
    seen = []
    monkeypatch.setattr(rr, "write_report", lambda rep: seen.append(_states(rep)))
    be._run_acquisition("", {"samples": samples, "kwargs": be.psd_acquisition_kwargs(values), "mode": "sizes",
                             "results_dir": results, "report": report})
    assert measured == [("A", values["asize.n_par_target"])]
    assert seen[-1] == [("A", "done"), ("B_fail", "failed")] and report["status"] == "finished"
