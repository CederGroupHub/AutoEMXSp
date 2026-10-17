#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-spectrum tab of the AutoEMX GUI: fit and quantify one spectrum of a sample
(``fit_and_quantify_spectrum_from_ledger``) or an external EMSA file (``fit_and_quantify_spectrum``),
with every option of ``Fit_Quant_Single_AutoEMX_Spectrum.py``. Results are not saved to the ledger.
"""

from __future__ import annotations

import base64
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from dash import ALL, Input, Output, State, ctx, dash_table, dcc, html, no_update

from autoemx.gui import backend as be
from autoemx.gui import plots as pl
from autoemx.gui.common import (
    GRAPH_CONFIG,
    JOBS,
    _chip,
    _fmt,
    _pick_file,
    form_values,
    get_analysis,
    get_info,
    param_sections,
    peak_overlaps_view,
    ui_value,
)

_UPLOAD_DIR = tempfile.mkdtemp(prefix="autoemx_gui_uploads_")
_SOURCES = [
    {"label": "Spectrum of the current sample", "value": "sample"},
    {"label": "External file (.msa, .emsa, .msg)", "value": "file"},
]


def _source_key(source: str, sample_dir: Optional[str], spectrum_id: Optional[str], path: Optional[str]) -> Optional[str]:
    if source == "file":
        return f"file:{path}" if path else None
    return f"sample:{sample_dir}:{spectrum_id}" if sample_dir and spectrum_id is not None else None


# =============================================================================
# Layout
# =============================================================================
def layout() -> List[Any]:
    sidebar = html.Div(
        [
            html.Div(
                [
                    html.Label("Spectrum", className="side-label"),
                    dcc.RadioItems(id="s-source", options=_SOURCES, value="sample", className="s-source"),
                    html.Div([
                        html.Button("◀", id="s-prev", className="btn btn-small", title="Previous spectrum"),
                        dcc.Dropdown(id="s-spectrum", options=[], clearable=False, placeholder="Spectrum",
                                     className="s-spectrum-dd"),
                        html.Button("▶", id="s-next", className="btn btn-small", title="Next spectrum"),
                    ], id="s-sample-block", className="s-row"),
                    html.Div([
                        html.Div([
                            dcc.Input(id="s-file-path", type="text", debounce=True, className="param-input",
                                      placeholder="Path of an EMSA spectrum file"),
                            html.Button("Browse…", id="s-browse", className="btn btn-small"),
                        ], className="s-row"),
                        dcc.Upload(id="s-upload", children=html.Div("or drop a file here"), className="s-upload",
                                   multiple=False),
                    ], id="s-file-block", hidden=True),
                ],
                className="side-block",
            ),
            html.Div(
                [
                    html.Div([
                        html.Span("Settings", className="side-title"),
                        html.Button("Use sample settings", id="s-reset", className="btn btn-small",
                                    title="Fill the form with the settings of the current sample's active quantification"),
                    ], className="side-title-row"),
                    html.Div(param_sections(be.SINGLE_SPECS, be.SINGLE_SECTIONS, open_sections=("single", "sfit"),
                                            extra={"single": [html.Div(id="s-overlap-msg")]}),
                             className="param-sections"),
                ],
                className="side-block side-params",
            ),
            html.Div(
                [
                    html.Div([html.Button("Fit and quantify", id="s-run-btn", className="btn btn-primary"),
                              html.Button("Cancel", id="s-cancel-btn", className="btn btn-danger", disabled=True)],
                             className="run-buttons"),
                    html.Div(id="s-run-msg", className="run-msg"),
                    html.Details([html.Summary("Log"), html.Pre(id="s-log", className="log")],
                                 id="s-log-details", className="log-box"),
                ],
                className="side-block run-block",
            ),
        ],
        className="sidebar",
    )
    main = html.Div(
        [
            html.Div([
                html.Div(id="s-title", className="s-title"),
                html.Div(id="s-metrics", className="summary"),
            ], className="analysis-bar"),
            html.Div([
                dcc.Checklist(id="s-log-y", options=[{"label": "Log scale", "value": "on"}], value=[],
                              className="param-check"),
                dcc.Checklist(id="s-bars", options=[{"label": "Background counts", "value": "on"}], value=["on"],
                              className="param-check"),
                html.Div([html.Span("Zoom to line", className="axis-tag"),
                          dcc.Dropdown(id="s-zoom", options=[], placeholder="full spectrum", className="color-dd")],
                         className="axis-box"),
            ], className="plot-controls"),
            dcc.Loading(dcc.Graph(id="s-graph", config=GRAPH_CONFIG, className="s-graph",
                                  figure=pl.single_spectrum_figure(None)),
                        type="circle", delay_show=400, parent_className="s-graph-box"),
            html.Div([
                html.Div([html.Div([html.B("Composition"),
                                    html.Button("Download CSV", id="s-dl-comp-btn", className="btn btn-small")],
                                   className="s-card-head"),
                          html.Div(id="s-comp")], className="card s-card"),
                html.Div([html.Div([html.B("Fitted peaks"),
                                    html.Button("Download CSV", id="s-dl-peaks-btn", className="btn btn-small")],
                                   className="s-card-head"),
                          dash_table.DataTable(
                              id="s-peaks", data=[], page_action="none", sort_action="native",
                              style_table={"maxHeight": "300px", "overflowY": "auto"},
                              style_cell={"fontSize": 12, "padding": "3px 8px", "textAlign": "right"},
                              style_header={"fontWeight": 600, "backgroundColor": "#f3f5f8"},
                              style_data_conditional=[{"if": {"filter_query": '{ref} = "yes"'}, "fontWeight": 600}],
                          )], className="card s-card"),
            ], className="s-cards"),
            dcc.Download(id="s-download"),
        ],
        className="main s-main",
    )
    return [
        dcc.Store(id="s-job-store"),
        dcc.Store(id="s-result"),
        dcc.Interval(id="s-poll", interval=1000, disabled=True),
        sidebar,
        main,
    ]


# =============================================================================
# Views
# =============================================================================
def _stored_record(sample_dir: Optional[str], spectrum_id: Optional[str]) -> Optional[pd.Series]:
    """Composition of a spectrum in the sample's active quantification."""
    if not sample_dir or spectrum_id is None:
        return None
    try:
        _, data = get_analysis(sample_dir, None)
    except Exception:
        return None
    rows = data.comps[data.comps["spectrum"].astype(str) == str(spectrum_id)] if data is not None else None
    return rows.iloc[0] if rows is not None and not rows.empty else None


def _composition_view(fit: Optional[Dict[str, Any]], stored: Optional[pd.Series]) -> List[Any]:
    if not fit or not fit.get("comp_at"):
        if fit is not None:
            return [html.Div("Fitted only (quantification off or interrupted).", className="hint")]
        return [html.Div("Run the fit to quantify the spectrum.", className="hint")]
    els = list(fit["comp_at"])
    head = [html.Th("")] + [html.Th(el) for el in els]
    rows = [
        html.Tr([html.Td("at%")] + [html.Td(_fmt(fit["comp_at"][el] * 100, 2)) for el in els]),
        html.Tr([html.Td("w%")] + [html.Td(_fmt((fit["comp_w"].get(el) or 0) * 100, 2)) for el in els]),
    ]
    if stored is not None:
        rows.append(html.Tr([html.Td("ledger at%", className="muted")] + [
            html.Td(_fmt(stored.get(f"{el} at%"), 2), className="muted") for el in els]))
        rows.append(html.Tr([html.Td("ledger w%", className="muted")] + [
            html.Td(_fmt(stored.get(f"{el} w%"), 2), className="muted") for el in els]))
    return [html.Table([html.Tr(head)] + rows, className="comp-table s-comp-table")]


def _metrics_view(fit: Optional[Dict[str, Any]]) -> List[Any]:
    if not fit:
        return []
    from autoemx.web.pipeline import QUANT_BEAM_KV, beam_energy_supports_quantification

    an_err = fit.get("analytical_error")
    chips = [
        _chip("R²", _fmt(fit.get("r_squared"), 5)),
        _chip("χ²ᵣ", _fmt(fit.get("redchi_sq"))),
        _chip("an. err", f"{_fmt(an_err * 100 if an_err is not None else None)} w%"),
        _chip("flag", pl.flag_label(fit.get("quant_flag"))),
    ]
    check = fit.get("element_check")
    if check is not None:
        found = check["added_quantified"] or check["added_not_quantified"] or check["possible"]
        chips.append(_chip("element check", be.element_check_message(check), "chip-warn" if found else ""))
    kv = fit.get("beam_energy_kV")
    if kv is not None:
        ok = beam_energy_supports_quantification(kv)
        chips.append(_chip("beam", f"{kv:g} kV" + ("" if ok else f" (standards are for {QUANT_BEAM_KV:g} kV)"),
                           "" if ok else "chip-warn"))
    return chips


def _peaks_payload(fit: Optional[Dict[str, Any]]):
    cols = [("line", "Line"), ("ref", "Ref."), ("center", "Centre (keV)"), ("th_energy", "Theory (keV)"),
            ("fwhm", "FWHM (keV)"), ("area", "Area"), ("height", "Height"), ("pb_ratio", "P/B")]
    columns = []
    for cid, name in cols:
        col = {"name": name, "id": cid}
        if cid not in ("line", "ref"):
            col.update(type="numeric", format={"specifier": ".4f" if cid in ("center", "th_energy", "fwhm") else ",.1f"})
        columns.append(col)
    rows = []
    for p in (fit or {}).get("peaks", []):
        rows.append({**{k: p.get(k) for k, _ in cols if k not in ("ref",)}, "ref": "yes" if p.get("reference") else ""})
    return rows, columns


# =============================================================================
# Callbacks
# =============================================================================
def register(app) -> None:
    @app.callback(
        Output("s-sample-block", "hidden"),
        Output("s-file-block", "hidden"),
        Input("s-source", "value"),
    )
    def toggle_source(source):
        return source != "sample", source != "file"

    @app.callback(
        Output("s-spectrum", "options"),
        Output("s-spectrum", "value"),
        Output("s-source", "value"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
        Input("single-request", "data"),
        Input("s-prev", "n_clicks"),
        Input("s-next", "n_clicks"),
        State("s-spectrum", "value"),
    )
    def spectra(sample_dir, _, request, _p, _n, current):
        source = no_update
        if not sample_dir:
            return [], None, source
        try:
            _, data = get_analysis(sample_dir, None)
        except Exception:
            return [], None, source
        if data is None:  # no spectra
            return [], None, source
        comps = data.comps
        options = []
        for _i, r in comps.iterrows():
            parts = [f"Spectrum {r['spectrum']}"]
            if r["particle"] not in (None, "", "None") and pd.notna(r["particle"]):
                parts.append(f"particle {r['particle']}")
            parts.append(str(r["status"]))
            options.append({"label": " · ".join(parts), "value": str(r["spectrum"])})
        ids = [o["value"] for o in options]
        trig = ctx.triggered_id
        value = current if current in ids else (ids[0] if ids else None)
        if trig == "single-request" and request and request.get("sample_dir") == sample_dir:
            value, source = request["spectrum_id"], "sample"
        elif trig in ("s-prev", "s-next") and ids:
            i = ids.index(value) if value in ids else 0
            value = ids[(i + (1 if trig == "s-next" else -1)) % len(ids)]
        return options, value, source

    @app.callback(
        Output({"type": "sparam", "key": ALL}, "value"),
        Input("sample-dd", "value"),
        Input("s-reset", "n_clicks"),
    )
    def defaults(sample_dir, _):
        info = None
        if sample_dir:
            try:
                info = get_info(sample_dir)
            except Exception:
                info = None
        values = be.single_param_values(info)
        return [ui_value(o["id"], values) for o in ctx.outputs_list]

    @app.callback(
        Output("s-file-path", "value"),
        Input("s-browse", "n_clicks"),
        Input("s-upload", "contents"),
        State("s-upload", "filename"),
        prevent_initial_call=True,
    )
    def choose_file(_, contents, filename):
        if ctx.triggered_id == "s-browse":
            path = _pick_file()
            return path if path else no_update
        if not contents:
            return no_update
        data = base64.b64decode(contents.split(",", 1)[1])
        path = os.path.join(_UPLOAD_DIR, Path(filename or "spectrum.msa").name)
        with open(path, "wb") as fh:
            fh.write(data)
        return path

    @app.callback(
        Output("s-graph", "figure"),
        Output("s-title", "children"),
        Output("s-metrics", "children"),
        Output("s-comp", "children"),
        Output("s-peaks", "data"),
        Output("s-peaks", "columns"),
        Output("s-zoom", "options"),
        Input("s-source", "value"),
        Input("s-spectrum", "value"),
        Input("s-file-path", "value"),
        Input("s-result", "data"),
        Input("s-log-y", "value"),
        Input("s-zoom", "value"),
        Input("s-bars", "value"),
        Input("sample-dd", "value"),
        Input({"type": "sparam", "key": "sfit__spectrum_lims__min"}, "value"),
        Input({"type": "sparam", "key": "sfit__spectrum_lims__max"}, "value"),
    )
    def view(source, sid, path, result, log_y, zoom, bars, sample_dir, lim_min, lim_max):
        key = _source_key(source, sample_dir, sid, path)
        fit = result if result and result.get("key") == key else None
        raw, title = None, ""
        try:
            if source == "file" and path:
                from autoemx.utils import load_msa

                energy, counts, _meta = load_msa(path)
                raw = {"energy": energy, "counts": counts}
                title = Path(path).name
            elif source == "sample" and sample_dir and sid is not None:
                info = get_info(sample_dir)
                raw = be.load_raw_spectrum(info, sid)
                title = f"{info.sample_id} · spectrum {sid}"
        except Exception as exc:
            return (pl.single_spectrum_figure(None), html.Span(f"Cannot read the spectrum: {exc}", className="err"),
                    [], [], [], [], [])
        if fit:
            title += " · fitted" + (" and quantified" if fit.get("quantified") else "")
        stored = _stored_record(sample_dir, sid) if source == "sample" else None
        rows, columns = _peaks_payload(fit)
        zoom_options = [{"label": p["line"].replace("_", " "), "value": p["line"]} for p in (fit or {}).get("peaks", [])]
        # Before fitting, show the range that will be fitted (spectrum limits of the form, or the defaults)
        try:
            lims = be.coerce_values({"sfit.spectrum_lims": [lim_min, lim_max]}, [be.SPECS_BY_KEY["sfit.spectrum_lims"]])
            lims = lims["sfit.spectrum_lims"] or list(be.dflt.spectrum_lims)
        except ValueError:
            lims = list(be.dflt.spectrum_lims)
        fig = pl.single_spectrum_figure(raw, fit, log_y=bool(log_y), zoom_line=zoom if fit else None,
                                        show_bckgrnd_cnts=bool(bars), title=title, channel_lims=lims)
        return (fig, html.B(title or "No spectrum"), _metrics_view(fit), _composition_view(fit, stored),
                rows, columns, zoom_options)

    @app.callback(
        Output("s-overlap-msg", "children"),
        Input({"type": "sparam", "key": "single__els_sample"}, "value"),
        Input({"type": "sparam", "key": "single__els_substrate"}, "value"),
        Input("s-source", "value"),
        Input("s-file-path", "value"),
        Input("sample-dd", "value"),
        Input({"type": "sparam", "key": "sfit__spectrum_lims__min"}, "value"),
        Input({"type": "sparam", "key": "sfit__spectrum_lims__max"}, "value"),
        Input({"type": "sparam", "key": "sfit__free_area_el_lines"}, "value"),
    )
    def overlap_msg(els_sample, els_substrate, source, path, sample_dir, lim_min, lim_max, free_area_el_lines):
        """Warn, as the elements are edited, about peak overlaps that may compromise the quantification."""
        try:
            els_sample = be._elements_list(els_sample)
            els_substrate = be._elements_list(els_substrate)
            fit_values = be.coerce_values(
                {"sfit.spectrum_lims": [lim_min, lim_max], "sfit.free_area_el_lines": free_area_el_lines},
                [be.SPECS_BY_KEY["sfit.spectrum_lims"], be.SPECS_BY_KEY["sfit.free_area_el_lines"]],
            )
            lims, free_area_el_lines = fit_values["sfit.spectrum_lims"], fit_values["sfit.free_area_el_lines"]
        except ValueError:
            return None
        from autoemx.web.pipeline import substrate_elements_note

        note = substrate_elements_note(els_sample, els_substrate)
        note = html.Div(note, className="hint") if note else None
        try:
            if source == "file":
                if not path:
                    return note
                from autoemx.web.pipeline import load_uploaded_spectrum

                _, _, geometry = load_uploaded_spectrum(path)
                overlaps = be.peak_overlaps(be.dflt.measurement_type, els_sample, els_substrate,
                                            geometry["beam_energy"], geometry["det_ch_offset"],
                                            geometry["det_ch_width"], lims, None, free_area_el_lines)
            else:
                info = get_info(sample_dir)
                if info is None:
                    return note
                overlaps = be.sample_peak_overlaps(info, els_sample, els_substrate, lims, free_area_el_lines)
        except Exception:
            return note
        return [note, peak_overlaps_view(overlaps)]

    @app.callback(
        Output("s-run-btn", "children"),
        Input({"type": "sparam", "key": "sfit__quantify"}, "value"),
    )
    def run_label(quantify):
        return "Fit and quantify" if quantify else "Fit spectrum"

    @app.callback(
        Output("s-job-store", "data"),
        Output("s-poll", "disabled"),
        Output("s-run-msg", "children"),
        Input("s-run-btn", "n_clicks"),
        State({"type": "sparam", "key": ALL}, "value"),
        State("s-source", "value"),
        State("s-spectrum", "value"),
        State("s-file-path", "value"),
        State("sample-dd", "value"),
        prevent_initial_call=True,
    )
    def run(_, __, source, sid, path, sample_dir):
        key = _source_key(source, sample_dir, sid, path)
        if key is None:
            return no_update, no_update, html.Span("Choose a spectrum first.", className="err")
        if source == "file" and not os.path.isfile(path):
            return no_update, no_update, html.Span(f"File not found: {path}", className="err")
        try:
            values = be.coerce_values(form_values(ctx.states_list[0]), be.SINGLE_SPECS)
        except ValueError as exc:
            return no_update, no_update, html.Div([html.B("Invalid settings: "), str(exc)], className="err")
        payload = {"source": source, "spectrum_id": sid, "path": path, "fit_kwargs": be.single_fit_kwargs(values)}
        desc = f"Fit of {Path(path).name}" if source == "file" else f"Fit of spectrum {sid}"
        try:
            job = JOBS.start("single", sample_dir or "", payload, desc)
        except RuntimeError as exc:
            return no_update, no_update, html.Span(str(exc), className="err")
        return {"id": job.job_id, "done": False, "key": key}, False, html.Span(f"{desc}…", className="running")

    @app.callback(
        Output("s-run-msg", "children", allow_duplicate=True),
        Output("s-log", "children"),
        Output("s-result", "data"),
        Output("s-job-store", "data", allow_duplicate=True),
        Output("s-poll", "disabled", allow_duplicate=True),
        Output("s-run-btn", "disabled"),
        Output("s-cancel-btn", "disabled"),
        Input("s-poll", "n_intervals"),
        State("s-job-store", "data"),
        prevent_initial_call=True,
    )
    def poll(_, job_data):
        job = JOBS.get(job_data["id"]) if job_data else None
        if job is None or job_data.get("done"):
            return no_update, no_update, no_update, no_update, True, False, True
        log = job.log_tail()
        res = job.result()
        if res is None:
            return (html.Span(f"{job.description} running… {job.elapsed():.0f} s", className="running"), log,
                    no_update, no_update, False, True, False)
        if job.cancelled:
            msg, result = html.Span(f"{job.description} cancelled.", className="warn"), no_update
        elif res.get("ok"):
            msg = html.Span(f"{job.description} done in {job.elapsed():.0f} s.", className="ok")
            result = {**res, "key": job_data["key"]}
        else:
            msg, result = html.Span(f"{job.description} failed: {res.get('error')}", className="err"), no_update
        return msg, log, result, {**job_data, "done": True}, True, False, True

    @app.callback(
        Output("s-run-msg", "children", allow_duplicate=True),
        Input("s-cancel-btn", "n_clicks"),
        State("s-job-store", "data"),
        prevent_initial_call=True,
    )
    def cancel(_, job_data):
        if job_data:
            JOBS.cancel(job_data["id"])
        return html.Span("Cancelling…", className="warn")

    @app.callback(
        Output("s-download", "data"),
        Input("s-dl-comp-btn", "n_clicks"),
        Input("s-dl-peaks-btn", "n_clicks"),
        State("s-result", "data"),
        prevent_initial_call=True,
    )
    def download(_c, _p, fit):
        if not fit:
            return no_update
        name = (fit.get("key") or "spectrum").split(":")[-1]
        name = Path(name).stem if "/" in name or "." in name else f"spectrum_{name}"
        if ctx.triggered_id == "s-dl-comp-btn":
            els = list(fit.get("comp_at") or {})
            df = pd.DataFrame({
                "element": els,
                "at%": [fit["comp_at"][el] * 100 for el in els],
                "w%": [(fit["comp_w"].get(el) or 0) * 100 for el in els],
            })
            an_err = fit.get("analytical_error")
            df["analytical error (w%)"] = an_err * 100 if an_err is not None else np.nan
            df["R2"] = fit.get("r_squared")
            df["reduced chi2"] = fit.get("redchi_sq")
            df["quant flag"] = fit.get("quant_flag")
            return dcc.send_data_frame(df.to_csv, f"{name}_composition.csv", index=False)
        return dcc.send_data_frame(pd.DataFrame(fit.get("peaks") or []).to_csv, f"{name}_peaks.csv", index=False)
