#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quantification tab of the AutoEMX GUI: quantify the spectra of the selected samples
(``batch_quantify_and_analyze``), and inspect the quantification runs of each sample.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import plotly.graph_objects as go
from dash import ALL, Input, Output, State, ctx, dash_table, dcc, html, no_update

from autoemx.gui import backend as be
from autoemx.gui import plots as pl
from autoemx.gui.common import (
    GRAPH_CONFIG,
    JOBS,
    SUMMARIES,
    _info_cache,
    _pick_folder,
    _pid,
    form_values,
    get_info,
    param_sections,
    peak_overlaps_view,
)

_NOTES = {
    "quant": "empty / 'saved' keeps each sample's value",
    "qafter": "optional",
}
_COLUMNS = [
    {"name": "Sample", "id": "sample", "editable": False},
    {"name": "Acquired", "id": "date", "editable": False},
    {"name": "Spectra", "id": "n_spectra", "type": "numeric", "editable": False},
    {"name": "Quantified", "id": "n_quantified", "type": "numeric", "editable": False},
    {"name": "Runs", "id": "n_runs", "type": "numeric", "editable": False},
    {"name": "Known sample elements", "id": "elements", "editable": True},
    {"name": "Known substrate elements", "id": "substrate", "editable": True},
]


# =============================================================================
# Layout
# =============================================================================
def layout() -> List[Any]:
    sidebar = html.Div(
        [
            html.Div(
                [
                    html.Div([html.Span("Settings", className="side-title")], className="side-title-row"),
                    html.Div(param_sections(be.QUANT_SPECS, be.QUANT_SECTIONS, _NOTES,
                                            open_sections=("quant", "qrun", "qafter")),
                             className="param-sections"),
                ],
                className="side-block side-params",
            ),
            html.Div(
                [
                    html.Div(id="q-selected-msg", className="hint"),
                    html.Div([
                        html.Button("Quantify selected samples", id="q-run-btn", className="btn btn-primary"),
                        html.Button("Cancel", id="q-cancel-btn", className="btn btn-danger", disabled=True),
                    ], className="run-buttons"),
                    html.Div(id="q-run-msg", className="run-msg"),
                    html.Div(id="q-progress", className="q-progress"),
                    html.Details([html.Summary("Log"), html.Pre(id="q-log", className="log")],
                                 id="q-log-details", className="log-box"),
                ],
                className="side-block run-block",
            ),
        ],
        className="sidebar",
    )
    main = html.Div(
        [
            html.Div(
                [
                    html.Div([
                        html.B("Samples"),
                        html.Span(" · tick the samples to quantify · click a column header to sort (e.g. by "
                                  "acquisition date) · edit the elements in the table to change them · "
                                  "substrate elements are not quantified, nor identified as sample elements: "
                                  "to quantify one (e.g. O in oxides), add it to the known sample elements",
                                  className="muted"),
                    ], className="q-table-title"),
                    html.Div([
                        dcc.Checklist(id="q-select-all", options=[{"label": "Select/Unselect all", "value": "all"}],
                                      value=[], className="param-check q-select-all",
                                      inputClassName="q-select-all-input"),
                        html.Button("Quantifiable elements", id="q-pt-btn", className="btn btn-small",
                                    title="Periodic table of the elements with P/B standards for the microscope "
                                          "and beam energy of the ticked samples"),
                        html.Button("Import spectra folder…", id="q-imp-open", className="btn btn-small",
                                    title="Create a sample from a folder of EMSA spectra (.msa/.emsa/.msg), "
                                          "e.g. acquired without AutoEMX, to quantify and analyse it here"),
                        html.Span(id="q-loading-msg", className="muted"),
                    ], className="q-table-actions"),
                    dash_table.DataTable(
                        id="q-table", columns=_COLUMNS, data=[], row_selectable="multi", selected_rows=[],
                        sort_action="native", sort_by=[], filter_action="native", page_action="none",
                        fixed_rows={"headers": True}, editable=True,
                        style_table={"maxHeight": "42vh", "overflowY": "auto"},
                        style_cell={"fontSize": 13, "padding": "4px 8px", "textAlign": "left",
                                    "minWidth": 60, "maxWidth": 380, "overflow": "hidden",
                                    "textOverflow": "ellipsis"},
                        style_cell_conditional=[
                            {"if": {"column_id": c}, "textAlign": "right", "width": "96px"}
                            for c in ("n_spectra", "n_quantified", "n_runs")
                        ] + [{"if": {"column_id": "date"}, "width": "130px"}],
                        style_header={"fontWeight": 600, "backgroundColor": "#f3f5f8"},
                        style_data_conditional=[
                            {"if": {"column_editable": True}, "backgroundColor": "#fbfcfe"},
                            {"if": {"state": "active"}, "backgroundColor": "#fff3e0", "border": "1px solid #ff9800"},
                        ],
                        tooltip_delay=300,
                    ),
                ],
                className="q-table-box",
            ),
            html.Div([
                html.Div([html.Div(id="q-details-title"),
                          html.Button("Analyse this sample →", id="q-to-analysis", className="btn btn-small")],
                         className="q-details-head"),
                html.Div(id="q-overlap-msg"),
                html.Div(id="q-details"),
            ], className="q-details"),
        ],
        className="main q-main",
    )
    return [
        dcc.Store(id="q-job-store"),
        dcc.Store(id="q-edits", data={}),
        dcc.Interval(id="q-poll", interval=1500, disabled=True),
        dcc.Interval(id="q-summary-poll", interval=1500, disabled=True),
        dcc.Store(id="q-imp-job"),
        dcc.Store(id="q-imp-scan"),
        dcc.Interval(id="q-imp-poll", interval=1000, disabled=True),
        _import_modal(),
        sidebar,
        main,
    ]


def _field(label: str, control: Any, hint: str = "") -> html.Div:
    return html.Div([html.Label(label, className="imp-label"), html.Div([control] + (
        [html.Div(hint, className="muted imp-hint")] if hint else []), className="imp-ctrl")], className="imp-row")


def _import_modal() -> html.Div:
    """Window to create a sample (spectra folder + ledger) from a folder of EMSA spectra."""
    micros = be.available_microscopes()
    return html.Div(
        html.Div([
            html.Div([html.B("Import a folder of EMSA spectra"),
                      html.Button("×", id="q-imp-close", className="btn btn-small", title="Close")],
                     className="pt-head"),
            html.Div("The spectra are copied into a new sample folder of the results folder, with a ledger, so "
                     "that the sample can be quantified and analysed like those acquired with AutoEMX. "
                     "Their energy calibration is read from the file headers. "
                     "The original files are not modified.", className="hint"),
            _field("Spectra folder", html.Div([
                dcc.Input(id="q-imp-folder", type="text", debounce=True, placeholder="/path/to/spectra",
                          className="imp-input"),
                html.Button("Browse…", id="q-imp-browse", className="btn btn-small"),
            ], className="imp-inline")),
            html.Div(id="q-imp-info", className="imp-info"),
            _field("Sample ID", dcc.Input(id="q-imp-id", type="text", className="imp-input"),
                   "name of the new sample folder"),
            _field("Known sample elements", dcc.Input(id="q-imp-els", type="text", placeholder="e.g. Pb, Mo, O",
                                                className="imp-input")),
            _field("Known substrate elements", dcc.Input(id="q-imp-sub", type="text", value="C",
                                                   className="imp-input"),
                   "elements of the substrate (e.g. carbon tape), empty if none; they are not quantified: "
                   "to quantify one (e.g. O in oxides), list it among the sample elements"),
            _field("Sample type", dcc.Dropdown(id="q-imp-type", value=be.cnst.S_POWDER_SAMPLE_TYPE, clearable=False,
                                               options=[{"label": t, "value": t} for t in be.IMPORT_SAMPLE_TYPES],
                                               className="imp-dd")),
            _field("Microscope", dcc.Dropdown(id="q-imp-micro", value=be.dflt.microscope_ID, clearable=False,
                                              options=[{"label": m, "value": m} for m in micros],
                                              className="imp-dd"),
                   "its P/B standards are used for the quantification"),
            _field("Beam energy", html.Div([dcc.Input(id="q-imp-kv", type="number", min=1, max=40, step="any",
                                                      className="imp-num"),
                                            html.Span("keV", className="param-unit")], className="imp-inline")),
            html.Div([html.Button("Import", id="q-imp-run", className="btn btn-primary"),
                      html.Button("Cancel", id="q-imp-cancel", className="btn btn-danger", disabled=True),
                      html.Div(id="q-imp-msg", className="run-msg")], className="imp-actions"),
        ], className="pt-window imp-window"),
        id="q-imp-modal", className="pt-backdrop", hidden=True,
    )


# =============================================================================
# Views
# =============================================================================
def _flag_figure(run: Dict[str, Any]) -> go.Figure:
    flags = sorted(run["flags"].items())
    n_missing = run["n_spectra"] - run["n_processed"]
    labels = [pl.flag_label(f) for f, _ in flags] + (["not processed"] if n_missing else [])
    values = [n for _, n in flags] + ([n_missing] if n_missing else [])
    colors = ["#2e7d32" if f in (0, -1) else "#c77700" if f == 8 else "#b0b6bf" for f, _ in flags]
    colors += ["#dfe3ea"] if n_missing else []
    fig = go.Figure(go.Bar(x=values, y=labels, orientation="h", marker_color=colors, text=values,
                           textposition="outside", cliponaxis=False))
    fig.update_layout(
        template="plotly_white", height=max(160, 34 * len(labels) + 60),
        margin=dict(l=10, r=40, t=30, b=10), title=dict(text=f"Quant flags · run {run['id']}", font=dict(size=13)),
        yaxis=dict(autorange="reversed"), xaxis=dict(title="spectra"),
    )
    return fig


def _details_view(info: be.SampleInfo) -> Tuple[List[Any], List[Any]]:
    runs = be.quantification_runs(info)
    head = [html.B(info.sample_id),
            html.Span(f" · {info.n_spectra} spectra · {len(runs)} quantification run(s)", className="muted")]
    if not runs:
        return head, [html.Div("Not quantified yet.", className="hint")]
    header = html.Tr([html.Th(h) for h in ("Run", "Elements", "Substrate", "Quantified", "Processed",
                                           "Fit tol.", "Channels", "Min counts", "Analyses", "Quant flags")])
    rows = []
    for r in reversed(runs):
        flags = " · ".join(f"{f}: {n}" for f, n in sorted(r["flags"].items()))
        lims = r["spectrum_lims"]
        rows.append(html.Tr([
            html.Td([str(r["id"]), html.Span(" active", className="badge badge-ok")] if r["active"] else str(r["id"])),
            html.Td(r["elements"]), html.Td(r["substrate"]),
            html.Td(f"{r['n_quantified']} / {r['n_spectra']}"),
            html.Td(f"{r['n_processed']} / {r['n_spectra']}"),
            html.Td(f"{r['fit_tolerance']:.0e}" if r["fit_tolerance"] else "—"),
            html.Td(f"{int(lims[0])}–{int(lims[1])}" if lims else "—"),
            html.Td(f"{r['min_total_counts_fraction']:.0%}" if r["min_total_counts_fraction"] is not None else "—"),
            html.Td(str(r["n_analyses"])), html.Td(flags or "—"),
        ], className="active-run" if r["active"] else ""))
    active = next((r for r in runs if r["active"]), runs[-1])
    return head, [
        html.Div([
            html.Div(html.Table([header] + rows, className="mix-table runs-table"), className="q-runs"),
            dcc.Graph(figure=_flag_figure(active), config=GRAPH_CONFIG, className="q-flags"),
        ], className="q-details-row"),
    ]


def _progress_view(progress: List[Dict[str, Any]], results: Optional[List[Dict[str, Any]]] = None) -> List[Any]:
    ok_by_sample = {r["sample"]: r["ok"] for r in (results or [])}
    items = []
    for p in progress:
        state = p["state"]
        if p["sample"] in ok_by_sample:
            state = "done" if ok_by_sample[p["sample"]] else "failed"
        total, done = p["total"], p["done"]
        frac = (done / total) if total else (1.0 if state in ("done", "failed") else 0.0)
        label = {"waiting": "waiting", "running": f"{done} / {total} spectra" if total else "starting…",
                 "done": f"done · {done} / {total} spectra" if total else "done (nothing to quantify)",
                 "failed": "failed"}[state]
        items.append(html.Div([
            html.Div([html.Span(p["sample"], className="q-prog-name"), html.Span(label, className="muted")],
                     className="q-prog-head"),
            html.Div(html.Div(className=f"q-bar-fill q-{state}", style={"width": f"{100 * frac:.0f}%"}),
                     className="q-bar"),
        ], className="q-prog-item"))
    return items


def _selected_ids(selected_rows: Optional[List[int]], rows: Optional[List[Dict[str, Any]]]) -> List[str]:
    """Sample folders ticked in the table (its selected_rows are indices in its data)."""
    rows = rows or []
    return [rows[i]["id"] for i in (selected_rows or []) if 0 <= i < len(rows)]


# =============================================================================
# Callbacks
# =============================================================================
def register(app) -> None:
    @app.callback(
        Output("q-table", "data"),
        Output("q-summary-poll", "disabled"),
        Output("q-loading-msg", "children"),
        Input("sample-dd", "options"),
        Input("q-summary-poll", "n_intervals"),
        Input("data-version", "data"),
        State("q-edits", "data"),
    )
    def table_data(options, _, __, edits):
        dirs = [o["value"] for o in (options or [])]
        summaries, n_loading = SUMMARIES.get(dirs)
        rows = []
        for d in dirs:
            s = summaries.get(d)
            row = {"id": d, "sample": Path(d).name, "date": "…", "n_spectra": None, "n_quantified": None,
                   "n_runs": None, "elements": "", "substrate": ""}
            if s:
                row.update({k: s[k] for k in ("date", "n_spectra", "n_quantified", "n_runs", "elements", "substrate")})
            row.update((edits or {}).get(d, {}))
            rows.append(row)
        msg = f"reading {n_loading} ledger(s)…" if n_loading else ""
        return rows, n_loading == 0, msg

    @app.callback(
        Output("q-edits", "data"),
        Input("q-table", "data_timestamp"),
        State("q-table", "data"),
        State("q-edits", "data"),
        prevent_initial_call=True,
    )
    def record_edits(_, rows, edits):
        # Elements edited in the table are kept until the next refresh of the summaries
        summaries, _n = SUMMARIES.get([r["id"] for r in rows or []])
        edits = dict(edits or {})
        for r in rows or []:
            s = summaries.get(r["id"])
            if not s:
                continue
            changed = {k: r[k] for k in ("elements", "substrate") if (r.get(k) or "") != (s.get(k) or "")}
            if changed:
                edits[r["id"]] = changed
            else:
                edits.pop(r["id"], None)
        return edits

    @app.callback(
        Output("q-table", "selected_rows"),
        Output("q-select-all", "value"),
        Input("q-select-all", "value"),
        Input("q-table", "selected_rows"),
        Input("q-table", "derived_virtual_row_ids"),
        State("q-table", "data"),
        prevent_initial_call=True,
    )
    def select_all(check, selected_rows, shown_ids, rows):
        # The checkbox ticks/unticks all listed samples (after filtering), and is ticked when they are all ticked.
        # The ticks drawn in the table are its selected_rows (indices in its data): the only source of truth.
        shown = set(shown_ids or [])
        if ctx.triggered_id == "q-select-all":
            if not check:
                return [], no_update
            return [i for i, r in enumerate(rows or []) if r["id"] in shown], no_update
        all_selected = bool(shown) and shown <= set(_selected_ids(selected_rows, rows))
        return no_update, (["all"] if all_selected else [])

    # Options shown as 'saved' display the value saved for the ticked samples when they all share it
    _number_keys = ["quant.min_bckgrnd_cnts", "quant.min_total_counts_fraction", "quant.fit_tolerance"]
    _choice_keys = ["quant.use_project_specific_std_dict", "quant.is_known_phase_mixture"]

    def _fmt_saved(key: str, value: Any) -> str:
        if key == "quant.fit_tolerance":
            return f"{value:.0e}"
        return f"{value:g}"

    @app.callback(
        [Output(_pid(k), "placeholder") for k in _number_keys]
        + [Output(_pid("quant.spectrum_lims", "__min"), "placeholder"),
           Output(_pid("quant.spectrum_lims", "__max"), "placeholder")]
        + [Output(_pid(k), "options") for k in _choice_keys],
        Input("q-table", "selected_rows"),
        Input("q-table", "data"),
    )
    def saved_placeholders(selected_rows, rows):
        selected = _selected_ids(selected_rows, rows)
        summaries, _n = SUMMARIES.get(selected)
        common = be.common_saved_values(list(summaries.values())) if selected and len(summaries) == len(selected) else {}
        out: List[Any] = [_fmt_saved(k, common[k]) if k in common else "saved" for k in _number_keys]
        lims = common.get("quant.spectrum_lims")
        out += [str(lims[0]) if lims else "min", str(lims[1]) if lims else "max"]
        for k in _choice_keys:
            saved_label = "saved" + (f" ({'yes' if common[k] else 'no'})" if k in common else "")
            out.append([{"label": saved_label if c == be.SAVED else c, "value": c}
                        for c in be.SPECS_BY_KEY[k].choices])
        return out

    # ---------------------------------------------------------------- import a folder of spectra
    @app.callback(
        Output("q-imp-modal", "hidden"),
        Input("q-imp-open", "n_clicks"),
        Input("q-imp-close", "n_clicks"),
        prevent_initial_call=True,
    )
    def import_window(*_):
        return ctx.triggered_id == "q-imp-close"

    @app.callback(
        Output("q-imp-folder", "value"),
        Input("q-imp-browse", "n_clicks"),
        prevent_initial_call=True,
    )
    def import_browse(_):
        return _pick_folder("Select the folder of EMSA spectra") or no_update

    @app.callback(
        Output("q-imp-info", "children"),
        Output("q-imp-scan", "data"),
        Output("q-imp-id", "value"),
        Output("q-imp-kv", "value"),
        Input("q-imp-folder", "value"),
        State("q-imp-id", "value"),
        prevent_initial_call=True,
    )
    def import_inspect(folder, sample_id):
        if not folder:
            return "", None, no_update, no_update
        try:
            scan = be.inspect_spectra_folder(folder)
        except Exception as exc:
            return html.Span(str(exc), className="err"), None, no_update, no_update
        name = Path(scan["folder"]).name
        if name.lower() == be.cnst.SPECTRA_DIR.lower():
            name = Path(scan["folder"]).parent.name  # <sample>/spectra
        if not scan["n_files"]:
            return html.Span("No EMSA spectra (.msa, .emsa, .msg) in this folder.", className="err"), None, \
                no_update, no_update
        parts: List[Any] = [f"{scan['n_files']} spectra"]
        beams = scan["beam_energies"]
        if beams:
            parts.append(f" · beam energy in the headers: {', '.join(f'{b:g}' for b in beams)} keV")
        if scan["calibration"]:
            off, width = scan["calibration"]
            parts.append(f" · energy calibration: offset {off * 1000:.1f} eV, {width * 1000:.3f} eV/channel")
        else:
            parts.append(html.Span(f" · cannot import: {scan['calibration_error']}.", className="err"))
        if len(beams) > 1:
            parts.append(html.Span(" · the files have different beam energies: check them.", className="err"))
        return parts, scan, sample_id or name, beams[0] if len(beams) == 1 else no_update

    @app.callback(
        Output("q-imp-job", "data"),
        Output("q-imp-poll", "disabled"),
        Output("q-imp-msg", "children"),
        Output("q-imp-cancel", "disabled"),
        Input("q-imp-run", "n_clicks"),
        State("q-imp-folder", "value"),
        State("folder", "value"),
        State("q-imp-id", "value"),
        State("q-imp-els", "value"),
        State("q-imp-sub", "value"),
        State("q-imp-type", "value"),
        State("q-imp-micro", "value"),
        State("q-imp-kv", "value"),
        State("q-imp-scan", "data"),
        prevent_initial_call=True,
    )
    def import_run(_, folder, results_dir, sample_id, els, sub, sample_type, micro, kv, scan):
        if scan and not scan.get("calibration"):
            return no_update, no_update, html.Span(f"Cannot import: {scan.get('calibration_error')}.",
                                                   className="err"), no_update
        try:
            kwargs = be.import_kwargs(folder, results_dir, sample_id, els, sub, sample_type, micro, kv)
            job = JOBS.start("import", "", {"kwargs": kwargs}, f"Import of {sample_id.strip()}")
        except (ValueError, RuntimeError) as exc:
            return no_update, no_update, html.Span(str(exc), className="err"), no_update
        sample_dir = str(Path(kwargs["results_path"]) / kwargs["samples"][0]["ID"])
        return {"id": job.job_id, "sample": sample_id.strip(), "dir": sample_dir}, False, \
            html.Span("Importing the spectra…", className="running"), False

    @app.callback(
        Output("q-imp-msg", "children", allow_duplicate=True),
        Output("q-imp-poll", "disabled", allow_duplicate=True),
        Output("scan-btn", "n_clicks", allow_duplicate=True),
        Output("q-run-msg", "children", allow_duplicate=True),
        Output("q-imp-cancel", "disabled", allow_duplicate=True),
        Input("q-imp-poll", "n_intervals"),
        State("q-imp-job", "data"),
        State("scan-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def import_poll(_, job_data, n_scan):
        job = JOBS.get(job_data["id"]) if job_data else None
        if job is None:
            return no_update, True, no_update, no_update, True
        res = job.result()
        if res is None:
            return no_update, False, no_update, no_update, False
        if job.cancelled:
            return html.Span("Import cancelled.", className="warn"), True, no_update, no_update, True
        if not res.get("ok"):
            log = job.log_tail(3000).strip().splitlines()
            detail = res.get("error") or (log[-1] if log else "")
            return html.Span(f"Import failed: {detail}", className="err"), True, no_update, no_update, True
        text = f"Imported as sample {job_data['sample']}: tick it in the table to quantify it."
        return html.Span(text, className="ok"), True, (n_scan or 0) + 1, html.Span(text, className="ok"), True

    @app.callback(
        Output("q-imp-msg", "children", allow_duplicate=True),
        Input("q-imp-cancel", "n_clicks"),
        State("q-imp-job", "data"),
        prevent_initial_call=True,
    )
    def import_cancel(_, job_data):
        if not job_data:
            return no_update
        JOBS.cancel(job_data["id"])
        # The import removes its partly built sample folder only when it ends by itself
        be.remove_unfinished_sample(job_data.get("dir"))
        return html.Span("Import cancelled.", className="warn")

    @app.callback(
        Output("q-run-btn", "children"),
        Input(_pid("qafter.run_analysis"), "value"),
    )
    def run_label(run_analysis):
        return "Quantify and analyse selected samples" if run_analysis else "Quantify selected samples"

    @app.callback(
        Output("q-selected-msg", "children"),
        Input("q-table", "selected_rows"),
        Input("q-table", "data"),
    )
    def selected_msg(selected_rows, rows):
        selected = _selected_ids(selected_rows, rows)
        n = len(selected)
        if not n:
            return "Tick samples in the table to quantify them."
        out: List[Any] = [f"{n} sample{'s' if n > 1 else ''} selected."]
        # Elements (as edited in the table) without P/B standards for the microscope and beam energy of each sample
        summaries, _n = SUMMARIES.get(selected)
        by_id = {r["id"]: r for r in rows or []}
        missing, no_stds = [], []
        for d in selected:
            summ = summaries.get(d)
            if not summ:
                continue
            try:
                els = be._elements_list(by_id.get(d, {}).get("elements") or summ.get("elements"))
            except ValueError:
                continue
            no_std = be.elements_without_standards(els, summ.get("microscope") or be.dflt.microscope_ID,
                                                   summ.get("beam_energy_keV"), summ.get("meas_mode") or "point")
            if no_std is None:
                no_stds.append(summ["sample"])
            elif no_std:
                missing.append(f"{summ['sample']}: {', '.join(no_std)}")
        if missing:
            out.append(html.Span(f" No P/B standard for {'; '.join(missing)}.", className="err"))
        if no_stds:
            out.append(html.Span(f" No standards at the beam energy of {', '.join(no_stds)}.", className="err"))
        return out

    @app.callback(
        Output("sample-dd", "value", allow_duplicate=True),
        Input("q-table", "active_cell"),
        prevent_initial_call=True,
    )
    def show_sample(cell):
        return cell["row_id"] if cell and cell.get("row_id") else no_update

    @app.callback(
        Output("q-details-title", "children"),
        Output("q-details", "children"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
    )
    def details(sample_dir, _):
        if not sample_dir:
            return "", html.Div("Click a sample in the table to see its quantification runs.", className="hint")
        try:
            return _details_view(get_info(sample_dir))
        except Exception as exc:
            return Path(sample_dir).name, html.Div(f"Cannot read the ledger: {exc}", className="err")

    @app.callback(
        Output("q-overlap-msg", "children"),
        Input("sample-dd", "value"),
        Input("q-table", "data"),
        Input("data-version", "data"),
    )
    def overlap_msg(sample_dir, rows, _):
        """Peak overlaps of the clicked sample, with its elements as edited in the table."""
        if not sample_dir:
            return None
        try:
            info = get_info(sample_dir)
        except Exception:
            return None
        if info is None:
            return None
        row = next((r for r in rows or [] if r.get("id") == sample_dir), {})
        try:
            els = be._elements_list(row["elements"]) if row.get("elements") else None
            substrate = be._elements_list(row["substrate"]) if "substrate" in row else None
        except ValueError:
            return None
        return peak_overlaps_view(be.sample_peak_overlaps(info, els, substrate))

    @app.callback(
        Output("main-tabs", "value", allow_duplicate=True),
        Input("q-to-analysis", "n_clicks"),
        prevent_initial_call=True,
    )
    def to_analysis(n):
        return "analysis" if n else no_update

    @app.callback(
        Output("q-job-store", "data"),
        Output("q-poll", "disabled"),
        Output("q-run-msg", "children"),
        Input("q-run-btn", "n_clicks"),
        State({"type": "qparam", "key": ALL}, "value"),
        State("q-table", "selected_rows"),
        State("q-table", "derived_virtual_row_ids"),
        State("q-table", "data"),
        prevent_initial_call=True,
    )
    def run(_, __, selected_rows, order, rows):
        selected = _selected_ids(selected_rows, rows)
        if not selected:
            return no_update, no_update, html.Span("Tick at least one sample in the table.", className="err")
        try:
            values = be.coerce_values(form_values(ctx.states_list[0]), be.QUANT_SPECS)
        except ValueError as exc:
            return no_update, no_update, html.Div([html.B("Invalid settings: "), str(exc)], className="err")
        summaries, _n = SUMMARIES.get(selected)
        by_id = {r["id"]: r for r in rows or []}
        ordered = [d for d in (order or []) if d in set(selected)] + [d for d in selected if d not in set(order or [])]
        samples, errors = [], []
        for d in ordered:
            row, s = by_id.get(d, {}), summaries.get(d) or {}
            entry: Dict[str, Any] = {"dir": d}
            for key, arg in (("elements", "els_sample"), ("substrate", "els_substrate")):
                text = (row.get(key) or "").strip()
                if s and text == (s.get(key) or ""):
                    continue  # unchanged: keep the ledger value
                try:
                    entry[arg] = be._elements_list(text)
                except ValueError as exc:
                    errors.append(f"{Path(d).name}: {exc}")
            if entry.get("els_sample") == []:
                errors.append(f"{Path(d).name}: no known sample elements")
            samples.append(entry)
        if errors:
            return no_update, no_update, html.Div([html.B("Invalid elements: "), "; ".join(errors)], className="err")
        payload = {"samples": samples, "kwargs": be.quantification_kwargs(values)}
        desc = f"Quantification of {len(samples)} sample{'s' if len(samples) > 1 else ''}"
        try:
            job = JOBS.start("quant", "", payload, desc)
        except RuntimeError as exc:
            return no_update, no_update, html.Span(str(exc), className="err")
        names = [Path(e["dir"]).name for e in samples]
        return {"id": job.job_id, "done": False, "samples": names, "dirs": [e["dir"] for e in samples]}, False, \
            html.Span(f"{desc} started…", className="running")

    @app.callback(
        Output("q-run-msg", "children", allow_duplicate=True),
        Output("q-progress", "children"),
        Output("q-log", "children"),
        Output("q-run-btn", "disabled"),
        Output("q-cancel-btn", "disabled"),
        Output("q-job-store", "data", allow_duplicate=True),
        Output("q-poll", "disabled", allow_duplicate=True),
        Output("q-log-details", "open"),
        Output("data-version", "data", allow_duplicate=True),
        Input("q-poll", "n_intervals"),
        State("q-job-store", "data"),
        State("data-version", "data"),
        prevent_initial_call=True,
    )
    def poll(_, job_data, version):
        job = JOBS.get(job_data["id"]) if job_data else None
        if job is None or job_data.get("done"):
            return (no_update,) * 6 + (True, no_update, no_update)
        log = job.log_tail(200000)
        progress = be.quant_progress(log, job_data["samples"])
        res = job.result()
        if res is None:
            n_done = sum(p["state"] == "done" for p in progress)
            msg = html.Span(f"{job.description} running… sample {min(n_done + 1, len(progress))} of "
                            f"{len(progress)} · {job.elapsed() / 60:.1f} min", className="running")
            # Ledgers of finished samples changed: refresh the tables every few polls
            return (msg, _progress_view(progress), log[-20000:], True, False, no_update, False, True, no_update)
        for d in job_data.get("dirs", []):
            _info_cache.pop(d, None)
        if res.get("ok"):
            text = f"{job.description} finished in {job.elapsed() / 60:.1f} min."
            msg = html.Span(text + (" " + res["warning"] if res.get("warning") else ""),
                            className="warn" if res.get("warning") else "ok")
        else:
            msg = html.Span(f"{job.description} failed: {res.get('error')}", className="err")
        return (msg, _progress_view(progress, res.get("samples")), log[-20000:], False, True,
                {**job_data, "done": True}, True, not res.get("ok"), (version or 0) + 1)

    @app.callback(
        Output("q-run-msg", "children", allow_duplicate=True),
        Input("q-cancel-btn", "n_clicks"),
        State("q-job-store", "data"),
        prevent_initial_call=True,
    )
    def cancel(_, job_data):
        if job_data:
            JOBS.cancel(job_data["id"])
        return html.Span("Cancelling… (spectra already quantified are kept in the ledgers)", className="warn")
