#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acquisition tab of the AutoEMX GUI: acquire (and optionally quantify) the spectra of a list of
samples with the electron microscope (``batch_acquire_and_analyze``), with every option of
``Run_Acquisition.py``, possibly continued by particle size statistics; or measure only the particle
size statistics of the samples (``collect_particle_statistics``). Settings and the sample list are kept
in the browser between sessions.

The results folder and the samples can be given in the page URL (``?acq=<JSON>``), e.g. by an external
script; see ``prefill_from_query``. With a run ID, the runs started from that page are reported to it
(see ``run_report``).
"""

from __future__ import annotations

import json
import os
import threading
import uuid
from typing import Any, Dict, List, Optional
from urllib.parse import parse_qs

from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from autoemx.gui import backend as be
from autoemx.gui import run_report
from autoemx.gui.common import JOBS, _chip, _pid, _row_class, _spec, form_values, param_sections

_NOTES = {
    "apowder": "used for powder samples",
    "asize": "used for particle sizes",
    "abulk": "used for bulk-type samples",
    "aquant": "used when quantifying during acquisition",
}
# Columns of the samples list: (key, header, placeholder)
_COLUMNS = [
    ("ID", "Sample ID", "e.g. Anorthite_mineral"),
    ("els", "Elements", "e.g. Ca, Al, Si, O"),
    ("x", "x (mm)", "x"),
    ("y", "y (mm)", "y"),
    ("cnd", "Candidate phases (optional)", "e.g. CaAl2Si2O8, CaO"),
]
_EXAMPLE_ROW = {"ID": "Anorthite_mineral", "els": "Ca, Al, Si, O", "x": -37.5, "y": -37.5, "cnd": "CaAl2Si2O8"}


def _sample_rows(rows: List[Dict[str, Any]]) -> List[Any]:
    """Editable rows of the samples list (plain text fields, so that their content can be edited in place)."""
    header = html.Div([html.Span(name, className=f"a-col a-col-{key}") for key, name, _ in _COLUMNS]
                      + [html.Span("", className="a-col a-col-del"), html.Span("", className="a-col a-col-del")],
                      className="a-row a-head")
    out = [header]
    for row in rows:
        uid = row["uid"]
        cells = [
            dcc.Input(id={"type": "a-cell", "uid": uid, "col": key}, value=row.get(key),
                      type="number" if key in ("x", "y") else "text", step="any", placeholder=placeholder,
                      className=f"param-input a-col a-col-{key}")
            for key, _, placeholder in _COLUMNS
        ]
        cells.append(html.Button("⧉", id={"type": "a-copy", "uid": uid}, n_clicks=0,
                                 title="Copy this sample (inserted below)", className="btn btn-small a-col a-col-del"))
        cells.append(html.Button("×", id={"type": "a-del", "uid": uid}, n_clicks=0, title="Remove this sample",
                                 className="btn btn-small a-col a-col-del"))
        out.append(html.Div(cells, className="a-row"))
    return out


def _rows_from_cells(cells: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Rows of the samples list from the ``ctx.states_list``/``inputs_list`` entry of its ALL-pattern cells."""
    rows: Dict[str, Dict[str, Any]] = {}
    for cell in cells or []:
        uid, col = cell["id"]["uid"], cell["id"]["col"]
        rows.setdefault(uid, {"uid": uid})[col] = cell.get("value")
    return list(rows.values())


def _row_from_sample(smp: Any) -> Dict[str, Any]:
    """Row of the samples list from a sample of ``Run_Acquisition.py`` (``{'ID', 'els', 'pos', 'cnd'}``)."""
    if not isinstance(smp, dict):
        raise ValueError("each sample must be a dict with keys 'ID', 'els', 'pos', 'cnd'")
    pos = smp.get("pos")
    if pos is None:
        x = y = None
    elif isinstance(pos, (list, tuple)) and len(pos) == 2:
        x, y = pos
    else:
        raise ValueError(f"sample {smp.get('ID')}: 'pos' must be [x, y] (mm)")

    def text(value):
        return ", ".join(str(v) for v in value) if isinstance(value, (list, tuple)) else str(value or "")

    return {"uid": uuid.uuid4().hex[:8], "ID": str(smp.get("ID") or ""), "els": text(smp.get("els")),
            "x": x, "y": y, "cnd": text(smp.get("cnd"))}


def prefill_from_query(search: Optional[str]) -> Optional[Dict[str, Any]]:
    """
    Results folder and sample rows given in the query string of the page URL, or None if there are none.

    The query is ``?acq=<URL-encoded JSON>``, with the JSON ``{"folder": ..., "samples": [...], "run_id": ...}``
    (all keys optional) and the samples as in ``Run_Acquisition.py``: ``{"ID", "els", "pos": [x, y], "cnd"}``,
    with ``els`` and ``cnd`` as lists or comma-separated text. Raises ``ValueError`` if the value is not valid.
    """
    values = parse_qs((search or "").lstrip("?")).get("acq")
    if not values:
        return None
    try:
        data = json.loads(values[0])
    except ValueError:
        raise ValueError("not valid JSON") from None
    if not isinstance(data, dict):
        raise ValueError("expected a JSON object with keys 'folder' and 'samples'")
    folder, samples = data.get("folder"), data.get("samples")
    if folder is not None and not isinstance(folder, str):
        raise ValueError("'folder' must be a path")
    if samples is not None and not isinstance(samples, list):
        raise ValueError("'samples' must be a list")
    run_id = data.get("run_id")
    return {"folder": folder or None, "rows": [_row_from_sample(s) for s in samples or []],
            "run_id": run_report.clean_run_id(run_id) if run_id is not None else None}


def _prefill_rows(search: Optional[str]) -> List[Dict[str, Any]]:
    try:
        prefill = prefill_from_query(search)
    except ValueError:
        return []  # reported by the prefill callback
    return prefill["rows"] if prefill else []


# =============================================================================
# Layout
# =============================================================================
def layout() -> List[Any]:
    sidebar = html.Div(
        [
            html.Div([
                html.Span("Measurement", className="side-title"),
                dcc.RadioItems(id="a-mode", value="eds", className="seg a-mode", persistence="v1",
                               persistence_type="local",
                               options=[{"label": label, "value": key} for key, label in be.ACQ_MODES.items()]),
                html.Div(id="a-mode-note", className="hint a-mode-note"),
            ], className="side-block"),
            html.Div(
                [
                    html.Div([html.Span("Settings", className="side-title"),
                              html.Button("Reset to defaults", id="a-reset", className="btn btn-small",
                                          title="Restore the default settings of Run_Acquisition.py")],
                             className="side-title-row"),
                    html.Div(param_sections(be.ACQ_SPECS, be.ACQ_SECTIONS, _NOTES,
                                            open_sections=("amicro", "asample", "aacq", "asize"), persist=True),
                             className="param-sections"),
                ],
                className="side-block side-params",
            ),
            html.Div(
                [
                    html.Div([
                        dcc.ConfirmDialogProvider(
                            html.Button("Start acquisition", id="a-run-btn", className="btn btn-primary"),
                            id="a-run-confirm",
                            message="Start the acquisition? The microscope stage will move to each sample "
                                    "and the beam will be used. Make sure the sample holder is loaded as defined.",
                        ),
                        dcc.ConfirmDialogProvider(
                            html.Button("Stop", id="a-cancel-btn", className="btn btn-danger", disabled=True),
                            id="a-cancel-confirm",
                            message="Stop the acquisition now? Spectra already acquired are kept. The microscope "
                                    "is left in its current state (stage position, beam).",
                        ),
                    ], className="run-buttons"),
                    html.Div(id="a-run-msg", className="run-msg"),
                    html.Div(id="a-progress", className="q-progress"),
                    html.Details([html.Summary("Log"), html.Pre(id="a-log", className="log")],
                                 id="a-log-details", className="log-box"),
                ],
                className="side-block run-block",
            ),
        ],
        className="sidebar",
    )
    main = html.Div(
        [
            html.Div(id="a-status", className="a-status"),
            html.Div(
                [
                    html.Div([
                        html.B("Samples to acquire"),
                        html.Span(" · one row per sample · position of the sample centre on the stage, in mm · "
                                  "elements quantified (substrate elements are set in the settings)",
                                  className="muted"),
                    ], className="q-table-title"),
                    html.Div([
                        html.Button("Add sample", id="a-add-row", className="btn btn-small"),
                        html.Button("Quantifiable elements", id="a-pt-btn", className="btn btn-small",
                                    title="Periodic table of the elements with P/B standards for the chosen "
                                          "microscope and beam energy"),
                        html.Button("Save as script…", id="a-script-btn", className="btn btn-small",
                                    title="Download a Python script running this acquisition, e.g. to run it "
                                          "on the microscope computer"),
                        html.Span(id="a-folder-msg", className="muted"),
                    ], className="q-table-actions"),
                    html.Div(id="a-prefill-msg", className="run-msg"),
                    html.Div(id="a-rows", className="a-rows"),
                    html.Div(id="a-table-msg", className="run-msg"),
                ],
                className="q-table-box",
            ),
            html.Div([
                html.B("Notes"),
                html.Ul([
                    html.Li("Each sample is saved in a sub-folder named after its ID, in the results folder "
                            "chosen at the top. Acquiring an existing sample adds spectra to it."),
                    html.Li("With manual navigation, or manual particle selection, a window asks you to centre "
                            "each spot or particle on the microscope."),
                    html.Li("After the acquisition, quantify the spectra in the Quantification tab "
                            "(unless they were quantified during the acquisition)."),
                ]),
            ], className="q-details a-notes"),
            dcc.Download(id="a-download"),
        ],
        className="main q-main",
    )
    return [
        dcc.Store(id="a-job-store"),
        # Run ID and samples given in the page URL by an external program, to report the runs to it. Kept for
        # the browser tab (the URL no longer has them after loading), so that runs started after a reload are
        # reported too
        dcc.Store(id="a-link", storage_type="session"),
        # Sample list, kept in the browser between sessions
        dcc.Store(id="a-saved-rows", storage_type="local"),
        dcc.Interval(id="a-poll", interval=2000, disabled=True),
        sidebar,
        main,
    ]


ON, OFF, HIDDEN = "on", "off", "hidden"  # parameter shown, greyed out, or not shown


# Settings of EDS spectra, not used by particle size statistics only
_EDS_ONLY_PARAMS = ("asample.sample_type", "asample.els_substrate", "asample.working_distance_tolerance",
                    "aacq.beam_energy", "aacq.target_Xsp_counts", "aacq.max_XSp_acquisition_time",
                    "aacq.quantify_spectra", "aacq.n_spectra", "aacq.min_n_spectra", "aacq.max_n_spectra",
                    "aimg.annotate_particle_images", "aimg.save_raw_images")
# Settings of a particle size statistics measurement (not of the EDS frames counted in the statistics)
_SIZE_MEASUREMENT_PARAMS = ("n_par_target", "min_par_diameter_um", "max_par_diameter_um", *be.ACQ_PSD_POWDER_FIELDS,
                            "size_stats_ci_tolerance", "size_stats_ci_quantities", "size_stats_max_n_frames",
                            "size_stats_max_time_min", "size_stats_plot_format")


def active_acq_params(values: Dict[str, Any], mode: str = "eds") -> Dict[str, str]:
    """State (on / off / hidden) of the acquisition sections ('section') and parameters ('section.name')
    with the current settings, for the measurement ``mode`` (see ``backend.ACQ_MODES``). Not listed: on."""
    sample_type = values.get("asample.sample_type")
    quantify = bool(values.get("aacq.quantify_spectra"))
    auto_bc = bool(values.get("aacq.auto_adjust_brightness_contrast"))
    active = {
        "apowder": ON if sample_type in be.ACQ_PARTICLE_TYPES else OFF,
        "abulk": ON if sample_type in be.ACQ_GRID_TYPES else OFF,
        "aquant": ON if quantify else OFF,
        "asample.is_auto_substrate_detection":
            ON if values.get("asample.sample_substrate_type") == be.cnst.CTAPE_SUBSTRATE_TYPE else OFF,
        "aacq.contrast": HIDDEN if auto_bc else ON,
        "aacq.brightness": HIDDEN if auto_bc else ON,
        "aacq.n_spectra": HIDDEN if quantify else ON,
        "aacq.min_n_spectra": ON if quantify else HIDDEN,
        "aacq.max_n_spectra": ON if quantify else HIDDEN,
        # Raw copies: of annotated images only
        "aimg.save_raw_images": ON if values.get("aimg.annotate_particle_images") else OFF,
    }
    # Brightness threshold: used by the threshold segmentation only
    for section in ("apowder", "asize"):
        model = values.get(f"{section}.par_segmentation_model")
        active[f"{section}.par_brightness_thresh"] = ON if model in (None, "threshold_bright") else OFF
    if mode == "eds":
        active.update({f"asize.{name}": HIDDEN for name in _SIZE_MEASUREMENT_PARAMS})
    elif not values.get("asize.size_stats_ci_tolerance"):
        # Limits of a measurement stopping at precision
        active.update({f"asize.{name}": OFF for name in ("size_stats_ci_quantities", "size_stats_max_n_frames",
                                                         "size_stats_max_time_min")})
    if mode == "sizes":
        # Particle settings in the size statistics section
        active.update({key: HIDDEN for key in _EDS_ONLY_PARAMS})
        active.update({"apowder": HIDDEN, "abulk": HIDDEN, "aquant": HIDDEN})
    return active


def _progress_view(progress: List[Dict[str, Any]]) -> List[Any]:
    items = []
    for p in progress:
        frac = min(1.0, p["done"] / p["total"]) if p["total"] else 0.0
        label = {"waiting": "waiting", "running": f"{p['done']} / {p['total']} spectra",
                 "done": f"done · {p['done']} spectra" if p["total"] else "done", "failed": "failed"}[p["state"]]
        sizes = p.get("sizes")
        if sizes and p["state"] == "running":  # particle size statistics (after the spectra, if any)
            total = sizes["total"]
            label = f"size band {sizes['band']}: {sizes['done']}" + (f" / {total} particles" if total else " particles")
            frac = min(1.0, sizes["done"] / total) if total else 0.0
        error = [html.Div(p["error"], className="err q-prog-error")] if p.get("error") else []
        items.append(html.Div([
            html.Div([html.Span(p["sample"], className="q-prog-name"), html.Span(label, className="muted")],
                     className="q-prog-head"),
            html.Div(html.Div(className=f"q-bar-fill q-{p['state']}",
                              style={"width": f"{100 * (1.0 if p['state'] == 'done' else frac):.0f}%"}),
                     className="q-bar"),
        ] + error, className="q-prog-item"))
    return items


# =============================================================================
# Callbacks
# =============================================================================
def register(app) -> None:
    @app.callback(
        Output("a-status", "children"),
        Input("main-tabs", "value"),
        Input(_pid("amicro.microscope_ID"), "value"),
        Input(_pid("aacq.beam_energy"), "value"),
    )
    def status(_, microscope, beam_kv):
        microscope = microscope or be.dflt.microscope_ID
        ok, text = be.microscope_status(microscope)
        chips = [_chip("microscope", text, "chip-ok" if ok else "chip-warn")]
        energies = be.standards_beam_energies(microscope)
        try:
            has_stds = be.quantifiable_elements(microscope, float(beam_kv)) is not None
        except (TypeError, ValueError):
            has_stds = False
        chips.append(_chip("P/B standards",
                           f"{float(beam_kv):g} kV available" if has_stds else
                           f"none at {beam_kv} kV (available: {', '.join(f'{e:g} kV' for e in energies) or 'none'})",
                           "chip-ok" if has_stds else "chip-warn"))
        return chips

    @app.callback(
        Output("a-folder-msg", "children"),
        Input("folder", "value"),
    )
    def folder_msg(folder):
        return f"Saved in: {folder}" if folder else "Choose the results folder at the top."

    @app.callback(
        Output("a-rows", "children"),
        Input("a-add-row", "n_clicks"),
        Input({"type": "a-del", "uid": ALL}, "n_clicks"),
        Input({"type": "a-copy", "uid": ALL}, "n_clicks"),
        Input("url", "search"),
        State({"type": "a-cell", "uid": ALL, "col": ALL}, "value"),
        State("a-saved-rows", "data"),
    )
    def render_rows(_add, deletes, _copies, search, _cells, saved):
        # Rows are re-drawn only when a sample is added or removed, not while typing in them
        trig = ctx.triggered_id
        if trig in (None, "url"):
            # Page load: samples given in the URL, else rows saved in the browser, or an example
            rows = _prefill_rows(search)
            if not rows:
                if ctx.states_list[0]:  # rows already drawn: keep them
                    return no_update
                rows = [r for r in (saved or []) if isinstance(r, dict)] or [dict(_EXAMPLE_ROW)]
        else:
            rows = _rows_from_cells(ctx.states_list[0])
            if trig == "a-add-row":
                rows.append({"ID": "", "els": "", "x": None, "y": None, "cnd": ""})
            elif isinstance(trig, dict) and trig.get("type") in ("a-del", "a-copy"):
                # New rows also trigger this callback (with n_clicks 0): only act on a real click
                if not ctx.triggered[0].get("value"):
                    return no_update
                i = next((k for k, r in enumerate(rows) if r["uid"] == trig["uid"]), None)
                if i is None:
                    return no_update
                if trig["type"] == "a-del":
                    rows.pop(i)
                else:
                    ids = {str(r.get("ID") or "") for r in rows}
                    copy_id, n = f"{rows[i].get('ID') or 'sample'}_copy", 2
                    while copy_id in ids:
                        copy_id, n = f"{rows[i].get('ID') or 'sample'}_copy{n}", n + 1
                    rows.insert(i + 1, {**rows[i], "uid": uuid.uuid4().hex[:8], "ID": copy_id})
        for r in rows:
            r.setdefault("uid", uuid.uuid4().hex[:8])
        return _sample_rows(rows)

    @app.callback(
        Output("a-prefill-msg", "children"),
        Output("folder", "value", allow_duplicate=True),
        Output("main-tabs", "value", allow_duplicate=True),
        Output("a-link", "data"),
        Input("url", "search"),
        State("a-link", "data"),
        prevent_initial_call="initial_duplicate",
    )
    def prefill(search, stored_link):
        # Results folder and samples given in the URL (the rows are drawn by render_rows)
        try:
            data = prefill_from_query(search)
        except ValueError as exc:
            return html.Span(f"The samples given in the link could not be read: {exc}", className="err"), \
                no_update, "acq", None
        if data is None:
            if not (stored_link and stored_link.get("run_id")):
                return no_update, no_update, no_update, no_update
            # Page reloaded (the URL no longer has the link): its results folder again, still reported
            return _reported_msg(stored_link["run_id"]), stored_link.get("folder") or no_update, no_update, no_update
        n = len(data["rows"])
        msg = [html.Span(f"{n} sample{'s' if n != 1 else ''} loaded from the link.", className="ok")] if n else []
        link = None
        if data["run_id"]:
            link = {"run_id": data["run_id"], "requested": [r["ID"] for r in data["rows"]], "folder": data["folder"]}
            msg.append(_reported_msg(data["run_id"]))
        return msg, data["folder"] or no_update, "acq", link

    def _reported_msg(run_id):
        return html.Span(f" Runs are reported to the program that opened this page (run {run_id}).",
                         className="muted")

    # Remove the samples from the address bar once loaded, so that reloading the page keeps the edits
    app.clientside_callback(
        """
        function(_) {
            if (window.location.search.indexOf('acq=') >= 0) {
                window.history.replaceState(window.history.state, '', window.location.pathname + window.location.hash);
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output("a-prefill-msg", "title"),
        Input("a-prefill-msg", "children"),
        prevent_initial_call=True,
    )

    @app.callback(
        Output({"type": "aparam", "key": ALL}, "value"),
        Input("a-reset", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset(_):
        from autoemx.gui.common import ui_value

        return [ui_value(o["id"], {s.key: s.default for s in be.ACQ_SPECS}) for o in ctx.outputs_list]

    @app.callback(
        Output({"type": "aparam-section", "key": ALL}, "className"),
        Output({"type": "aparam-row", "key": ALL}, "className"),
        Input({"type": "aparam", "key": ALL}, "value"),
        Input("a-mode", "value"),
    )
    def grey_out(_, mode):
        active = active_acq_params(form_values(ctx.inputs_list[0]), mode)
        extra = {ON: "", OFF: " disabled", HIDDEN: " hidden-param"}
        sections = ["param-section" + extra[active.get(o["id"]["key"], ON)] for o in ctx.outputs_list[0]]
        rows = [_row_class(_spec(o["id"])) + extra[active.get(_spec(o["id"]).key, ON)] for o in ctx.outputs_list[1]]
        return sections, rows

    @app.callback(Output("a-mode-note", "children"), Input("a-mode", "value"))
    def mode_note(mode):
        return {
            "eds": "",
            "eds+sizes": "After the spectra of each sample, frames are searched for particle size statistics: the "
                         "particles found in the frames searched by EDS count, and only more frames are searched to "
                         "reach the target. With another segmentation model or threshold than EDS, the frames "
                         "searched by EDS are segmented again with them, from their saved images. Powder samples, "
                         "automatic navigation and particle selection only.",
            "sizes": "Particle sizes only: no spectra are acquired, the elements of the samples are not needed.",
        }.get(mode, "")

    @app.callback(
        Output("a-table-msg", "children"),
        Output("a-saved-rows", "data"),
        Input({"type": "a-cell", "uid": ALL, "col": ALL}, "value"),
        Input(_pid("amicro.microscope_ID"), "value"),
        Input(_pid("aacq.beam_energy"), "value"),
        Input("a-mode", "value"),
        State("sample-dd", "options"),
    )
    def check_table(_, microscope, beam_kv, mode, options):
        rows = _rows_from_cells(ctx.inputs_list[0])
        if not ctx.inputs_list[0]:
            return html.Span("Add the samples to acquire.", className="muted"), no_update
        return _check_message(rows, options, microscope or be.dflt.microscope_ID, beam_kv, mode), rows

    def _check_message(rows, options, microscope, beam_kv, mode="eds"):
        try:
            samples = be.acquisition_samples(rows, require_elements=mode != "sizes")
        except ValueError as exc:
            return html.Span(str(exc), className="err")
        existing = {os.path.basename(o["value"]) for o in (options or [])}
        again = [s["ID"] for s in samples if s["ID"] in existing]
        out = [html.Span(f"{len(samples)} sample{'s' if len(samples) != 1 else ''} to acquire.", className="ok")]
        if again:
            out.append(html.Span(f" Already in the results folder (spectra will be added): {', '.join(again)}.",
                                 className="warn"))
        missing = []
        for smp in samples if mode != "sizes" else []:
            no_std = be.elements_without_standards(smp["els"], microscope, beam_kv)
            if no_std:
                missing.append(f"{smp['ID']}: {', '.join(no_std)}")
        if missing:
            out.append(html.Span(f" No P/B standard (these elements cannot be quantified): {'; '.join(missing)}.",
                                 className="err"))
        return out

    def _settings(states_form, rows, folder, mode):
        values = be.coerce_values(form_values(states_form), be.ACQ_SPECS)
        samples = be.acquisition_samples(rows, require_elements=mode != "sizes")
        if not folder:
            raise ValueError("Choose the results folder at the top")
        if mode == "sizes":
            return [{"ID": s["ID"], "pos": s["pos"]} for s in samples], be.psd_acquisition_kwargs(values)
        return samples, be.acquisition_kwargs(values, mode)

    @app.callback(
        Output("a-download", "data"),
        Output("a-run-msg", "children", allow_duplicate=True),
        Input("a-script-btn", "n_clicks"),
        State({"type": "aparam", "key": ALL}, "value"),
        State({"type": "a-cell", "uid": ALL, "col": ALL}, "value"),
        State("folder", "value"),
        State("a-mode", "value"),
        prevent_initial_call=True,
    )
    def script(_, __, ___, folder, mode):
        try:
            samples, kwargs = _settings(ctx.states_list[0], _rows_from_cells(ctx.states_list[1]), folder, mode)
        except ValueError as exc:
            return no_update, html.Div([html.B("Invalid settings: "), str(exc)], className="err")
        name = "Collect_Particle_Statistics_GUI.py" if mode == "sizes" else "Run_Acquisition_GUI.py"
        return dict(content=be.acquisition_script(samples, kwargs, folder, mode), filename=name), \
            html.Span("Script downloaded.", className="ok")

    @app.callback(
        Output("a-job-store", "data"),
        Output("a-poll", "disabled"),
        Output("a-run-msg", "children"),
        Input("a-run-confirm", "submit_n_clicks"),
        State({"type": "aparam", "key": ALL}, "value"),
        State({"type": "a-cell", "uid": ALL, "col": ALL}, "value"),
        State("folder", "value"),
        State("a-link", "data"),
        State("a-mode", "value"),
        prevent_initial_call=True,
    )
    def run(_, __, ___, folder, link, mode):
        try:
            samples, kwargs = _settings(ctx.states_list[0], _rows_from_cells(ctx.states_list[1]), folder, mode)
        except ValueError as exc:
            return no_update, no_update, html.Div([html.B("Invalid settings: "), str(exc)], className="err")
        payload = {"samples": samples, "kwargs": kwargs, "results_dir": os.path.abspath(os.path.expanduser(folder)),
                   "mode": mode}
        what = "Particle size statistics" if mode == "sizes" else "Acquisition"
        desc = f"{what} of {len(samples)} sample{'s' if len(samples) > 1 else ''}"
        report = None
        if link and link.get("run_id"):
            report = run_report.new_report(link["run_id"], [s["ID"] for s in samples], payload["results_dir"],
                                           link.get("requested"))
            payload["report"] = report
            # Written before the job process updates it, so that a report of a previous run is never read
            run_report.write_report(report)
        try:
            job = JOBS.start("acquisition", "", payload, desc)
        except RuntimeError as exc:
            if report:
                run_report.finish_report(report, "failed", str(exc))
                run_report.write_report(report)
            return no_update, no_update, html.Span(str(exc), className="err")
        if report:
            threading.Thread(target=run_report.follow_job, args=(job, report), name=f"report-{report['run_id']}",
                             daemon=False).start()
        return ({"id": job.job_id, "done": False, "samples": [s["ID"] for s in samples],
                 "max_n": kwargs.get("max_n_spectra", 0)}, False, html.Span(f"{desc} started…", className="running"))

    @app.callback(
        Output("a-run-msg", "children", allow_duplicate=True),
        Output("a-progress", "children"),
        Output("a-log", "children"),
        Output("a-run-btn", "disabled"),
        Output("a-cancel-btn", "disabled"),
        Output("a-job-store", "data", allow_duplicate=True),
        Output("a-poll", "disabled", allow_duplicate=True),
        Output("a-log-details", "open"),
        Output("scan-btn", "n_clicks", allow_duplicate=True),
        Input("a-poll", "n_intervals"),
        State("a-job-store", "data"),
        State("scan-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def poll(_, job_data, scans):
        job = JOBS.get(job_data["id"]) if job_data else None
        if job is None or job_data.get("done"):
            return (no_update,) * 6 + (True, no_update, no_update)
        log = job.log_tail(200000)
        progress = be.acquisition_progress(log, job_data["samples"], job_data["max_n"])
        res = job.result()
        if res is None:
            msg = html.Span(f"{job.description} running… {job.elapsed() / 60:.1f} min", className="running")
            return msg, _progress_view(progress), log[-20000:], True, False, no_update, False, True, no_update
        for p in progress:
            if p["state"] == "running":
                p["state"] = "done" if res.get("ok") else "failed"
        failed = [p["sample"] for p in progress if p["state"] == "failed"]
        if res.get("ok") and failed:
            msg = html.Span(f"{job.description} finished in {job.elapsed() / 60:.1f} min. "
                            f"{len(failed)} of {len(progress)} samples failed: {', '.join(failed)} "
                            "(see below, and the log).", className="warn")
        elif res.get("ok"):
            msg = html.Span(f"{job.description} finished in {job.elapsed() / 60:.1f} min.", className="ok")
        else:
            msg = html.Span(f"{job.description} stopped: {res.get('error')}", className="err")
        # Rescan the results folder so that the new samples appear in the other tabs
        return (msg, _progress_view(progress), log[-20000:], False, True, {**job_data, "done": True}, True,
                not res.get("ok") or bool(failed), (scans or 0) + 1)

    @app.callback(
        Output("a-run-msg", "children", allow_duplicate=True),
        Input("a-cancel-confirm", "submit_n_clicks"),
        State("a-job-store", "data"),
        prevent_initial_call=True,
    )
    def cancel(_, job_data):
        if job_data:
            JOBS.cancel(job_data["id"])
        return html.Span("Stopping…", className="warn")
