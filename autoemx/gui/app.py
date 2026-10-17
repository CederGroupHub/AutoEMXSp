#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoEMX GUI (Dash), with three tabs:

- Acquisition: acquire the spectra of a list of samples with the microscope (``batch_acquire_and_analyze``).
- Quantification: quantify the spectra of one or more samples (``batch_quantify_and_analyze``).
- Analysis: set every clustering parameter, run ``analyze_sample``, and explore the results in
  interactive 3D / ternary / 2D clustering plots linked to the spectra and SEM images.
- Single spectrum: fit and quantify one spectrum of a sample, or an external EMSA file.

Launch with ``python -m autoemx.gui``.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from dash import ALL, Dash, Input, Output, State, ctx, dash_table, dcc, get_asset_url, html, no_update
from urllib.parse import quote

from flask import Response, abort, send_file

from autoemx.config.runtime_configs import PowderMeasurementConfig
from autoemx.core.particle_outlines import load_outlines
from autoemx.gui import backend as be
from autoemx.gui import plots as pl
from autoemx.gui import tab_acquisition, tab_quantification, tab_single
from autoemx.gui.periodic_table import periodic_table
from autoemx.gui.common import (
    GRAPH_CONFIG as _GRAPH_CONFIG,
    JOBS,
    SUMMARIES,
    _chip,
    _fmt,
    _info_cache,
    _mtime,
    _open_path,
    _pick_folder,
    _pid,
    _point_spectrum,
    _row_class,
    _spec,
    form_values,
    get_analysis,
    get_info,
    get_particle_sizes,
    has_sizes,
    param_sections,
    peak_overlaps_view,
    register_clear_buttons,
    ui_value,
)

_CITATION = (
    "A. Giunto et al., Accurate SEM-EDS Quantification, Automation, and Machine Learning Enable "
    "High-Throughput Compositional Characterization of Powders, Nature Communications 17, 9735 (2026)."
)
_SECTION_NOTES = {
    "dbscan": "used when method = dbscan",
    "aitchison": "used when geometry = aitchison or auto",
    "merge": "k-means with automatic k and merging on",
    "mixture": "used when mixture decomposition is on",
    "plot": "saved PNG figures",
}
_OPEN_SECTIONS = {"filter", "clust"}
_PLOT_OPTIONS = [
    {"label": "Discarded", "value": pl.OPT_DISCARDED},
    {"label": "±1σ", "value": pl.OPT_STD},
    {"label": "Candidates", "value": pl.OPT_REFS},
    {"label": "Best mixtures", "value": pl.OPT_MIXTURES},
    {"label": "Centroids", "value": pl.OPT_CENTROIDS},
]
_COLOR_BY = [
    {"label": "Cluster", "value": "cluster"},
    {"label": "Analytical error", "value": "an_err"},
    {"label": "Particle", "value": "particle"},
    {"label": "Quant flag", "value": "quant_flag"},
    {"label": "R²", "value": "r_squared"},
]
def active_params(values: Dict[str, Any]) -> Tuple[Dict[str, bool], Dict[str, bool]]:
    """Which analysis parameter sections and single parameters are used with the current settings."""
    method = values.get("clust.method")
    kmeans = method == "kmeans"
    auto_k = values.get("clust.k_forced") in (None, "")
    sections = {
        "dbscan": method == "dbscan",
        "aitchison": values.get("clust.geometry") in ("aitchison", "auto"),
        "merge": kmeans and auto_k and bool(values.get("clust.auto_merge_clusters")),
        "mixture": bool(values.get("clust.do_matrix_decomposition")),
    }
    params = {
        "clust.k_forced": kmeans,
        "clust.k_finding_method": kmeans and auto_k,
        "clust.max_k": kmeans and auto_k,
        "clust.auto_merge_clusters": kmeans and auto_k,
    }
    return sections, params


# =============================================================================
# Layout
# =============================================================================
def _header(results_dir: str) -> html.Div:
    return html.Div(
        [
            html.Div([html.Img(src=get_asset_url("autoemx-logo-dark.svg"), alt="AutoEMX", className="brand-logo")],
                     className="brand-box"),
            html.Div(
                [
                    dcc.Input(id="folder", type="text", value=results_dir, debounce=True,
                              placeholder="Results folder (contains sample folders with ledger.json)",
                              className="folder-input"),
                    html.Button("Browse…", id="browse-btn", className="btn"),
                    html.Button("Rescan", id="scan-btn", className="btn"),
                ],
                className="folder-box",
            ),
            html.Div(id="scan-msg", className="scan-msg"),
            html.Div(
                [html.Span("Sample", className="header-label"),
                 dcc.Dropdown(id="sample-dd", options=[], placeholder="Select a sample", clearable=False,
                              className="header-sample")],
                className="header-sample-box",
            ),
        ],
        className="header",
    )


#: Particle segmentation models (the first one uses a brightness threshold), installed or not
PAR_SEGMENTATION_MODELS = be.par_segmentation_models()


PAGES = [("acq", "Acquisition"), ("quant", "Quantification"), ("analysis", "Analysis"),
         ("single", "Single spectrum")]


def _nav() -> dcc.Tabs:
    # Open on Acquisition on the microscope computer, else on Quantification
    first = "acq" if be.microscope_status()[0] else "quant"
    return dcc.Tabs(
        id="main-tabs", value=first, className="main-tabs",
        children=[dcc.Tab(label=label, value=value, className="main-tab", selected_className="main-tab--selected")
                  for value, label in PAGES],
    )


def _sidebar() -> html.Div:
    return html.Div(
        [
            html.Div(
                [
                    html.Label("Sample", className="side-label"),
                    html.Div(id="sample-info", className="sample-info"),
                ],
                className="side-block",
            ),
            html.Div(
                [
                    html.Div([
                        html.Span("Parameters", className="side-title"),
                        html.Button("Load settings of shown analysis", id="load-settings-btn", className="btn btn-small",
                                    title="Fill the form with the settings of the analysis selected above the plot"),
                    ], className="side-title-row"),
                    html.Div(param_sections(be.PARAM_SPECS, be.SECTIONS, _SECTION_NOTES, _OPEN_SECTIONS),
                             className="param-sections"),
                ],
                className="side-block side-params",
            ),
            html.Div(
                [
                    html.Div([
                        html.Button("Run analysis", id="run-btn", className="btn btn-primary"),
                        html.Button("Cancel", id="cancel-btn", className="btn btn-danger", disabled=True),
                    ], className="run-buttons"),
                    html.Div(id="run-msg", className="run-msg"),
                    html.Details([html.Summary("Log"), html.Pre(id="run-log", className="log")],
                                 id="log-details", className="log-box"),
                ],
                className="side-block run-block",
            ),
        ],
        id="analysis-sidebar",
        className="sidebar",
    )


def _plot_controls() -> html.Div:
    return html.Div(
        [
            dcc.RadioItems(id="plot-mode", value="3d", inline=True, className="seg",
                           options=[{"label": "3D", "value": "3d"}, {"label": "Ternary", "value": "ternary"},
                                    {"label": "2D", "value": "2d"}]),
            html.Div([html.Span("X", className="axis-tag"),
                      dcc.Dropdown(id="axis-x", clearable=False, className="axis-dd")], className="axis-box"),
            html.Div([html.Span("Y", className="axis-tag"),
                      dcc.Dropdown(id="axis-y", clearable=False, className="axis-dd")], className="axis-box"),
            html.Div([html.Span("Z", className="axis-tag"),
                      dcc.Dropdown(id="axis-z", clearable=False, className="axis-dd")], className="axis-box"),
            html.Div([html.Span("Colour", className="axis-tag"),
                      dcc.Dropdown(id="color-by", options=_COLOR_BY, value="cluster", clearable=False,
                                   className="color-dd")], className="axis-box"),
            html.Div([html.Span("Zoom on", className="axis-tag"),
                      dcc.Dropdown(id="zoom-on", options=pl.zoom_options(None), value=pl.ZOOM_FULL,
                                   clearable=False, className="color-dd")],
                     className="axis-box", title="Zoom the axes on all spectra or on one cluster. In 3D and "
                     "ternary, the mouse wheel zooms further and right-drag pans (3D); double-click resets."),
            dcc.Checklist(id="plot-opts", options=_PLOT_OPTIONS, value=pl.DEFAULT_OPTIONS, inline=True,
                          className="plot-opts"),
            html.Div([html.Span("mix. conf ≥", className="axis-tag"),
                      dcc.Input(id="mix-conf", type="number", min=0, max=1, step=0.05, debounce=True,
                                value=pl.DEFAULT_MIN_MIXTURE_CONF, className="param-input conf-input")],
                     className="axis-box", title="Best mixtures are drawn only above this confidence score"),
            html.Div([
                html.Button("Save plot as HTML", id="save-html-btn", className="btn btn-small",
                            title="Save the interactive (rotatable) plot in the analysis folder"),
                html.Button("Open analysis folder", id="open-folder-btn", className="btn btn-small"),
            ], className="summary-actions"),
        ],
        className="plot-controls",
    )


def _spectrum_panel() -> html.Div:
    return html.Div(
        [
            html.Div([
                html.Button("◀", id="prev-btn", className="btn btn-small", title="Previous spectrum"),
                html.Div(id="spectrum-title", className="spectrum-title"),
                html.Button("▶", id="next-btn", className="btn btn-small", title="Next spectrum"),
            ], className="spectrum-nav"),
            html.Div(id="spectrum-info", className="spectrum-info"),
            dcc.Graph(id="spectrum-graph", config=_GRAPH_CONFIG, className="spectrum-graph",
                      figure=pl.spectrum_figure(None)),
            html.Div([
                html.Button("Fit spectrum", id="fit-btn", className="btn",
                            title="Re-fit and quantify this spectrum to show the fitted model (≈10 s)"),
                html.Button("Cancel", id="fit-cancel-btn", className="btn btn-danger", hidden=True),
                html.Button("Show in SEM image", id="show-image-btn", className="btn",
                            title="Show the image of the particle where this spectrum was collected"),
                html.Button("Open in Single spectrum", id="to-single-btn", className="btn",
                            title="Fit and quantify this spectrum with custom settings in the Single spectrum tab"),
                dcc.Checklist(id="log-y", options=[{"label": "Log scale", "value": "on"}], value=[],
                              className="param-check"),
                html.Span(id="fit-msg", className="fit-msg"),
            ], className="spectrum-actions"),
        ],
        id="spectrum-panel",
        className="spectrum-panel",
    )


def _tabs() -> dcc.Tabs:
    return dcc.Tabs(
        id="tabs", value="plot", className="tabs",
        children=[
            dcc.Tab(id="plot-tab", label="Clustering", value="plot", children=[
                _plot_controls(),
                # Mouse-wheel zoom: Plotly's in 2D; in 3D and ternary, assets/cluster_zoom.js
                dcc.Loading(dcc.Graph(id="cluster-graph", config={**_GRAPH_CONFIG, "scrollZoom": "cartesian"},
                                      className="cluster-graph"),
                            type="circle", delay_show=400, parent_className="graph-fill"),
            ]),
            dcc.Tab(id="clusters-tab", label="Clusters", value="clusters", children=[html.Div(id="clusters-view", className="tab-body")]),
            dcc.Tab(id="spectra-tab", label="Spectra", value="spectra", children=[html.Div([
                html.Div("Click a row to show the spectrum. Filter with e.g. “> 5”, “< -3” or text.",
                         className="hint"),
                dash_table.DataTable(
                    id="spectra-table", page_size=25, sort_action="native", filter_action="native",
                    style_table={"overflowX": "auto"}, style_as_list_view=True,
                    style_cell={"fontSize": 13, "padding": "4px 8px",
                                "textAlign": "right", "minWidth": 56},
                    style_header={"fontWeight": 600, "backgroundColor": "#f3f5f8"},
                ),
            ], className="tab-body")]),
            dcc.Tab(id="dist-tab", label="Distributions", value="dist", children=[html.Div(
                dcc.Graph(id="dist-graph", config=_GRAPH_CONFIG, className="dist-graph"),
                className="tab-body graph-fill")]),
            dcc.Tab(id="psd-tab", label="Particle sizes", value="psd", children=[html.Div([
                html.Div([
                    dcc.RadioItems(id="psd-mode", value="sample", inline=True, className="seg",
                                   options=[{"label": "This sample", "value": "sample"},
                                            {"label": "Compare samples", "value": "compare"},
                                            {"label": "Compare segmentations", "value": "results"}]),
                    html.Div([
                        dcc.Dropdown(id="psd-result-dd", clearable=False, placeholder="Computed now (not saved)",
                                     className="psd-result-dd"),
                        html.Button("Set as active", id="psd-set-active-btn", className="btn btn-small",
                                    title="Make the shown statistics the active ones of the sample: shown and used "
                                          "by default (e.g. recalculated when EDS adds frames)"),
                        html.Details([
                            html.Summary("Segment again…"),
                            html.Div([
                                html.Span("model", className="axis-tag"),
                                dcc.Dropdown(id="psd-seg-model", clearable=False, className="psd-seg-model",
                                             options=[{"label": be.segmentation_model_label(m), "value": m}
                                                      for m in PAR_SEGMENTATION_MODELS],
                                             value=PAR_SEGMENTATION_MODELS[0]),
                                html.Span([  # threshold model only
                                    html.Span("threshold", className="axis-tag"),
                                    dcc.Input(id="psd-seg-thresh", type="number", min=0, max=255, step=1,
                                              className="param-input psd-seg-input"),
                                ], id="psd-seg-thresh-box", className="psd-seg-field"),
                                html.Span("min pixels", className="axis-tag"),
                                dcc.Input(id="psd-seg-px", type="number", min=1, step=1,
                                          className="param-input psd-seg-input"),
                                html.Span([  # other models only
                                    html.Span("options", className="axis-tag"),
                                    dcc.Input(id="psd-seg-params", type="text", debounce=True,
                                              placeholder='model defaults, or e.g. {"pred_iou_thresh": 0.7}',
                                              className="param-input psd-seg-params"),
                                ], id="psd-seg-params-box", className="psd-seg-field", hidden=True,
                                    title="Options of the segmentation model, as JSON (see the model's "
                                          "documentation); empty for its defaults"),
                                html.Button("Segment again", id="psd-reseg-btn", className="btn btn-small",
                                            title="Segment the frames again from their saved images with these "
                                                  "settings: a new segmentation and its statistics (the earlier ones "
                                                  "are kept)"),
                            ], className="psd-seg-row"),
                            html.Div(id="psd-seg-note", className="param-note"),
                        ], className="psd-seg-box", open=True),
                        html.Span(id="psd-result-msg"),
                    ], id="psd-result-controls", className="psd-result-controls"),
                    html.Div([
                        html.Div(dcc.Dropdown(id="psd-compare", multi=True, placeholder="Samples with particle sizes"),
                                 id="psd-compare-samples", className="psd-compare-dd"),
                        dcc.RadioItems(id="psd-kind", value="cumulative", inline=True, className="seg",
                                       options=[{"label": "Cumulative", "value": "cumulative"},
                                                {"label": "Histogram", "value": "histogram"}]),
                        dcc.RadioItems(id="psd-basis", value="number", inline=True, className="seg",
                                       options=[{"label": "Number", "value": "number"},
                                                {"label": "Volume", "value": "volume"}]),
                        html.Button("Export CSV", id="psd-export-btn", className="btn btn-small",
                                    title="Download the statistics of the compared samples"),
                        dcc.Download(id="psd-download"),
                    ], id="psd-compare-controls", className="psd-compare-controls", hidden=True),
                ], className="psd-controls"),
                html.Div(id="psd-banner", className="psd-banner"),
                dcc.Graph(id="psd-graph", config=_GRAPH_CONFIG, className="psd-graph", figure=pl.psd_figure(None)),
                html.Div(id="psd-tables", className="psd-tables"),
                html.Div([
                    html.Div("Frames", className="psd-section-title"),
                    html.Div([
                        html.Button("◀", id="psd-frame-prev", className="btn btn-small", title="Previous frame"),
                        dcc.Dropdown(id="psd-frame-dd", clearable=False, placeholder="No frames",
                                     className="psd-frame-dd"),
                        html.Button("▶", id="psd-frame-next", className="btn btn-small", title="Next frame"),
                        dcc.Checklist(id="psd-frame-out", options=[{"label": "Leave this frame out",
                                                                    "value": "out"}],
                                      value=[], className="param-check"),
                        html.Span(id="psd-frame-caption", className="muted"),
                    ], className="psd-frame-bar"),
                    html.Div("Click a particle to leave it out of the statistics (or count it again).",
                             className="hint"),
                    dcc.Graph(id="psd-frame-graph", className="psd-frame-graph",
                              config={**_GRAPH_CONFIG, "scrollZoom": True}),
                    html.Div([
                        html.Span(id="psd-pending", className="psd-pending"),
                        html.Button("Recalculate and save", id="psd-recalc-btn", className="btn btn-primary",
                                    title="Recalculate the statistics without the frames and particles left out, "
                                          "and save them in the sample ledger"),
                        html.Button("Reset", id="psd-reset-btn", className="btn",
                                    title="Back to the frames and particles left out of the saved statistics"),
                        html.Span(id="psd-recalc-msg"),
                    ], className="psd-recalc-bar"),
                ], id="psd-frames-section", className="psd-frames-section"),
            ], className="tab-body")]),
            dcc.Tab(id="sem-tab", label="SEM images", value="sem", children=[html.Div([
                html.Div([
                    html.Button("◀", id="sem-prev-btn", className="btn btn-small", title="Previous image"),
                    html.Span(id="sem-counter", className="muted sem-counter"),
                    html.Button("▶", id="sem-next-btn", className="btn btn-small", title="Next image"),
                    html.Div(id="sem-caption", className="sem-caption"),
                    dcc.Checklist(id="sem-spots", options=[{"label": "Show spectrum spots", "value": "on"}],
                                  value=["on"], className="param-check"),
                    html.Button("Open image", id="sem-open-btn", className="btn btn-small",
                                title="Open the full-resolution image in your image viewer"),
                ], className="sem-bar"),
                dcc.Graph(id="sem-graph", className="sem-graph",
                          config={**_GRAPH_CONFIG, "scrollZoom": True}),
                html.Div(id="sem-gallery", className="sem-gallery"),
            ], className="tab-body sem-body")]),
            dcc.Tab(id="figures-tab", label="Saved figures", value="figures", children=[html.Div(id="figures-view", className="tab-body")]),
            dcc.Tab(id="settings-tab", label="Settings used", value="settings", children=[html.Div(id="settings-view", className="tab-body")]),
        ],
    )


def _main() -> html.Div:
    return html.Div(
        [
            html.Div(
                [
                    html.Div([
                        html.Label("Analysis", className="side-label"),
                        dcc.Dropdown(id="analysis-dd", options=[], placeholder="No analysis yet", clearable=False,
                                     className="analysis-dd"),
                        html.Button("Set as active", id="set-active-btn", className="btn btn-small",
                                    title="Make the shown analysis the active one: its settings are reused by "
                                          "later analyses and it is shown by default"),
                    ], className="analysis-select"),
                    html.Div(id="summary", className="summary"),
                ],
                className="analysis-bar",
            ),
            html.Div(id="action-msg", className="action-msg"),
            html.Div([html.Div(_tabs(), className="tabs-col"), _spectrum_panel()], className="main-row"),
        ],
        className="main",
    )


def _periodic_table_modal() -> html.Div:
    """Window with the periodic table of the elements quantifiable with the standards of a microscope."""
    micros = be.available_microscopes()
    return html.Div(
        html.Div([
            html.Div([html.B("Quantifiable elements"),
                      html.Button("×", id="pt-close", className="btn btn-small", title="Close")],
                     className="pt-head"),
            html.Div([
                html.Span("Microscope", className="axis-tag"),
                dcc.Dropdown(id="pt-micro", options=[{"label": m, "value": m} for m in micros],
                             value=be.dflt.microscope_ID, clearable=False, className="pt-dd"),
                html.Span("Beam energy", className="axis-tag"),
                dcc.Dropdown(id="pt-kv", options=[], clearable=False, className="pt-dd"),
            ], className="pt-controls"),
            html.Div(id="pt-info", className="pt-info"),
            html.Div(id="pt-table"),
            html.Div([
                html.Span([html.Span(className="pt-swatch pt-quant"), " quantifiable"]),
                html.Span([html.Span(className="pt-swatch"), " no standard"]),
                html.Span([html.Span(className="pt-swatch pt-sample pt-quant"), " sample element"]),
                html.Span([html.Span(className="pt-swatch pt-missing"), " sample element without standard"]),
            ], className="pt-legend"),
        ], className="pt-window"),
        id="pt-modal", className="pt-backdrop", hidden=True,
    )


def build_layout(results_dir: str) -> html.Div:
    return html.Div(
        [
            dcc.Store(id="data-version", data=0),
            dcc.Store(id="selected-spectrum"),
            dcc.Store(id="job-store"),
            dcc.Store(id="fit-job-store"),
            dcc.Store(id="fit-result"),
            dcc.Store(id="folder-store"),
            dcc.Store(id="sem-image"),
            dcc.Store(id="mode-forced", data=False),
            dcc.Store(id="sem-paths"),
            dcc.Interval(id="poll", interval=1000, disabled=True),
            dcc.Store(id="psd-excl"),
            dcc.Store(id="psd-job-store"),
            dcc.Interval(id="psd-poll", interval=1000, disabled=True),
            dcc.Store(id="single-request"),
            dcc.Store(id="pt-elements"),
            # Query string of the page URL: samples to acquire given by an external script (?acq=...)
            dcc.Location(id="url", refresh=False),
            _header(results_dir),
            _nav(),
            html.Div(tab_acquisition.layout(), id="page-acq", className="body page", hidden=True),
            html.Div(tab_quantification.layout(), id="page-quant", className="body page", hidden=True),
            html.Div([_sidebar(), _main()], id="page-analysis", className="body page", hidden=True),
            html.Div(tab_single.layout(), id="page-single", className="body page", hidden=True),
            html.Div(["If you use AutoEMX, please cite: ", _CITATION], className="footer"),
            _periodic_table_modal(),
        ],
        className="app",
    )


# =============================================================================
# Views
# =============================================================================
def _summary_view(info: be.SampleInfo, data: be.AnalysisData) -> List[Any]:
    s = data.summary
    chips = []
    if s.get("has_labels"):
        chips.append(_chip("clusters", s.get("n_clusters", 0)))
        if s.get("silhouette") is not None and s.get("n_clusters", 0) > 1:
            chips.append(_chip("silhouette", f"{s['silhouette']:.2f}"))
        chips.append(_chip("clustered", s["n_clustered"], "chip-ok"))
        if s["n_noise"]:
            chips.append(_chip("noise", s["n_noise"]))
        chips.append(_chip("discarded", s["n_discarded"], "chip-warn"))
    else:
        chips.append(_chip("status", "not analysed", "chip-warn"))
        chips.append(_chip("quantified", s["n_spectra"] - s["n_not_quantified"]))
    chips.append(_chip("not quantified", s["n_not_quantified"], "chip-muted"))
    chips.append(_chip("method", f"{data.config.method} · {data.config.geometry} · {data.features}"))
    return chips


def _clusters_view(data: be.AnalysisData) -> List[Any]:
    if not len(data.centroids):
        return [html.Div("No clusters for this analysis. Run the analysis to cluster the spectra.", className="hint")]
    unit = data.unit
    cards = []
    for i, (c, sd) in enumerate(zip(data.centroids, data.stdevs)):
        n = data.n_points[i] if i < len(data.n_points) else "?"
        comp = html.Table(
            [html.Tr([html.Th(el) for el in data.elements])]
            + [html.Tr([html.Td(f"{v * 100:.1f} ± {e * 100:.1f}") for v, e in zip(c, sd)])],
            className="comp-table",
        )
        body = [html.Div([html.Span(className="swatch", style={"background": pl.cluster_color(i)}),
                          html.B(f"Cluster {i}"), html.Span(f"  {n} spectra · {unit}", className="muted")],
                         className="card-title"), comp]
        cand = data.candidates[i] if i < len(data.candidates) else None
        if cand:
            body.append(html.Div([html.B("Candidate phases: "), ", ".join(
                f"{name} (conf. {conf:.2f})" for name, conf in cand)], className="card-line"))
        mixes = data.mixtures[i] if i < len(data.mixtures) else []
        if mixes:
            header = html.Tr([html.Th("#"), html.Th("Phases"), html.Th("Confidence"), html.Th("Recon. error"),
                              html.Th("Molar fractions (mean ± σ)")])
            rows = []
            for rank, m in enumerate(mixes, start=1):
                refs = m.get("refs") or []
                means = m.get("means") or ([m.get("mean"), 1 - m["mean"]] if m.get("mean") is not None and len(refs) == 2 else [])
                stds = m.get("stddevs") or ([m.get("stddev")] * len(means) if m.get("stddev") is not None else [])
                fr = ", ".join(
                    f"{r}: {mu:.2f} ± {sg:.2f}" for r, mu, sg in zip(refs, means, stds + [np.nan] * len(means))
                ) if means else "—"
                rows.append(html.Tr([html.Td(m.get("rank", rank)), html.Td(" + ".join(refs)),
                                     html.Td(_fmt(m.get("conf_score"), 2)), html.Td(_fmt(m.get("recon_error"), 3)),
                                     html.Td(fr)]))
            body.append(html.Div([html.B("Mixtures of candidate phases"),
                                  html.Table([header] + rows, className="mix-table")], className="card-line"))
        cards.append(html.Div(body, className="card"))
    out: List[Any] = [html.Div(cards, className="cards")]
    if data.clusters_table is not None:
        df = data.clusters_table
        out.append(html.Details([
            html.Summary("Clusters.csv"),
            dash_table.DataTable(
                data=df.astype(str).to_dict("records"), columns=[{"name": c, "id": c} for c in df.columns],
                style_table={"overflowX": "auto"}, style_cell={"fontSize": 12},
                style_header={"fontWeight": 600},
            ),
        ], className="raw-csv"))
    return out


def _table_payload(data: be.AnalysisData):
    df = data.comps.copy()
    cols = ["spectrum", "particle", "cluster", "status"]
    cols += [f"{el} at%" for el in data.elements] + [f"{el} w%" for el in data.elements]
    cols += ["an_err", "quant_flag", "r_squared", "redchi_sq", "total_counts", "comment"]
    names = {"an_err": "An. err w%", "quant_flag": "Flag", "r_squared": "R²", "redchi_sq": "Red. χ²",
             "total_counts": "Counts"}
    df["cluster"] = df["cluster"].astype("Int64")
    df["id"] = df["spectrum"].astype(str)
    numeric = set(cols) - {"spectrum", "status", "comment", "particle", "cluster"}
    columns = []
    for c in cols:
        col = {"name": names.get(c, c), "id": c}
        if c in numeric:
            col["type"] = "numeric"
            col["format"] = {"specifier": ".5f" if c == "r_squared" else ".2f"}
        if c == "quant_flag":
            col["format"] = {"specifier": ".0f"}
        if c == "total_counts":
            col["format"] = {"specifier": ",.0f"}
        columns.append(col)
    for c in ("particle",):
        df[c] = df[c].astype(str).replace({"None": "", "nan": ""})
    records = df[cols + ["id"]].replace({np.nan: None}).astype(object).where(df[cols + ["id"]].notna(), None)
    records = records.to_dict("records")
    for r in records:
        if r.get("cluster") is not None:
            r["cluster"] = int(r["cluster"])
    styles = [
        {"if": {"filter_query": '{status} = "discarded"'}, "color": "#8a8f98"},
        {"if": {"filter_query": '{status} = "not quantified"'}, "color": "#b5b9c0", "fontStyle": "italic"},
        {"if": {"state": "active"}, "backgroundColor": "#fff3e0", "border": "1px solid #ff9800"},
    ]
    for k in range(len(data.centroids)):
        styles.append({"if": {"filter_query": f"{{cluster}} = {k}", "column_id": "cluster"},
                       "backgroundColor": pl.cluster_color(k), "color": "white", "fontWeight": 600})
    for c in ("status", "comment", "spectrum", "particle"):
        styles.append({"if": {"column_id": c}, "textAlign": "left"})
    return records, columns, styles


def _figures_view(data: be.AnalysisData) -> List[Any]:
    if not data.images:
        return [html.Div("No saved figures for this analysis.", className="hint")]
    items = []
    for path in data.images:
        name = Path(path).name
        items.append(html.Figure([
            html.A(html.Img(src=f"/analysis-file?path={path}", className="saved-img"),
                   href=f"/analysis-file?path={path}", target="_blank"),
            html.Figcaption(name),
        ], className="saved-fig"))
    return [html.Div(f"Folder: {data.folder}", className="hint"), html.Div(items, className="gallery")]


def _sample_info_view(info: be.SampleInfo) -> List[Any]:
    return [
        html.Div([html.Span("Elements", className="k"), html.Span(", ".join(info.elements), className="v")]),
        html.Div([html.Span("Substrate", className="k"), html.Span(", ".join(info.substrate) or "—", className="v")]),
        html.Div([html.Span("Spectra", className="k"),
                  html.Span(f"{info.n_spectra} ({info.n_quantified} quantified)", className="v")]),
        html.Div([html.Span("Analyses", className="k"), html.Span(str(len(info.analyses)), className="v")]),
        html.Div([html.Span("Sizes", className="k"),
                  html.Span("particle size data" if info.has_sizes else "—", className="v")]),
        None if info.has_eds else html.Div("No EDS spectra: the analysis settings below do not apply.",
                                           className="hint"),
        peak_overlaps_view(be.sample_peak_overlaps(info)),
    ]


# Rows of the particle size statistics table: (label, key in ParticleSizeStats, key of the confidence interval)
_PSD_STAT_ROWS = [
    ("Particles", "n_par_analysed", None), ("Mean", "mean", None), ("Std. dev.", "stdev", None),
    ("Min", "min", None), ("D10", "D10", "D10"), ("D25", "D25", None), ("D50", "median", "D50"),
    ("D75", "D75", None), ("D90", "D90", "D90"), ("Max", "max", None), ("D[4,3]", "D43", "D43"),
    ("Dv10", "Dv10", "Dv10"), ("Dv50", "Dv50", "Dv50"), ("Dv90", "Dv90", "Dv90"),
]


# Origin of the particles counted (keys of ParticleStatsResult.n_par_by_selection)
_PSD_SELECTION_LABELS = {"size_stats": "size measurement", "auto": "measured by EDS",
                         "eds_frames": "in EDS frames", "manual": "chosen manually", "list": "from list"}


def _psd_banner(data: be.ParticleSizeData, info: Optional[be.SampleInfo] = None) -> List[Any]:
    """Where the statistics come from (frames and segmentations), and what they can miss."""
    r = data.result
    n_eds_frames = sum(s.n_frames for s in r.eds_scans)
    if data.source.startswith("measurement"):
        text = f"From the particle size measurement ({len(r.bands)} size band{'s' if len(r.bands) != 1 else ''})"
        if n_eds_frames:
            text += f" and {n_eds_frames} frames searched by EDS"
        elif data.source == "measurement+eds":
            text += " and the particles measured by EDS"
        text += "."
    elif data.source == "eds_frames":
        text = (f"From the {n_eds_frames} frames searched by EDS only. They observe particles from "
                f"{data.eds_d_min_um:.3g} µm: smaller particles need a particle size measurement.")
    else:
        text = ("From the particles measured by EDS only (acquired before EDS frames were recorded): sizes are "
                "not weighted by the area scanned.")
    if data.stored:
        text += f" Statistics of {r.datetime.replace('T', ' ')}."
    else:
        text += " Computed now from the sample data, not saved."
    chips = [_chip("particles", r.size_stats.n_par_analysed)]
    if r.precision_reached is not None:
        chips.append(_chip("target precision", "reached" if r.precision_reached else "not reached",
                           "chip-ok" if r.precision_reached else "chip-warn"))
    chips += [_chip(_PSD_SELECTION_LABELS.get(sel, sel), n, "chip-muted") for sel, n in r.n_par_by_selection.items()]
    if r.excluded_frame_ids or r.excluded_particles:
        chips.append(_chip("left out", f"{len(r.excluded_frame_ids)} frames, {len(r.excluded_particles)} "
                                       "particles", "chip-warn"))
    if r.overlapping_frame_ids:
        chips.append(_chip("overlapping frames", len(r.overlapping_frame_ids), "chip-muted"))
    out = [html.Div(text, className="psd-source")]
    if info is not None:
        segmentation = info.ledger.segmentation(r.segmentation_id)
        eds = info.ledger.segmentation(r.eds_segmentation_id)
        names = ([f"size frames {segmentation.label}"] if segmentation is not None else []) + (
            [f"frames searched by EDS {eds.label}"] if eds is not None else [])
        if names:
            out.append(html.Div("Segmentation: " + "; ".join(names) + ".", className="muted"))
        if be.mixed_segmentations(data):
            out.append(html.Div("The size frames and the frames searched by EDS are segmented differently: segment "
                                "them again with the same settings for consistent statistics.", className="warn"))
    if r.missing_frame_ids:
        out.append(html.Div(f"{len(r.missing_frame_ids)} frame(s) without result of this segmentation left out "
                            f"(no saved image): {', '.join(r.missing_frame_ids)}.", className="warn"))
    return out + [html.Div(chips, className="summary")]


def _psd_tables(data: be.ParticleSizeData) -> List[Any]:
    r = data.result
    w, u, ci = r.size_stats, r.unweighted_size_stats, r.confidence_intervals_95

    def val(stats, key):
        v = getattr(stats, key, None) if stats is not None else None
        return "—" if v is None else (str(v) if key == "n_par_analysed" else f"{v:.3g}")

    rows = [html.Tr([html.Th(""), html.Th("Weighted"), html.Th("95 % CI"), html.Th("Unweighted")])]
    for label, key, ci_key in _PSD_STAT_ROWS:
        interval = ci.get(ci_key) if ci_key else None
        rows.append(html.Tr([html.Td(label), html.Td(val(w, key)),
                             html.Td(f"{interval[0]:.3g}–{interval[1]:.3g}" if interval else ""),
                             html.Td(val(u, key))]))
    stats_table = html.Div([html.Div("Equivalent diameter (µm)", className="psd-table-title"),
                            html.Table(rows, className="comp-table")])

    def scans_table(title: str, scans) -> Optional[html.Div]:
        if not scans:
            return None
        head = html.Tr([html.Th(h) for h in ("Diameters (µm)", "Frame width (µm)", "Frames", "Area (mm²)",
                                             "Particles", "Particles/mm²")])
        body = [html.Tr([
            html.Td(f"{s.d_min_um:.3g}–{s.d_max_um:.3g}"),
            html.Td(f"{s.frame_width_um:.3g}" if s.frame_width_um else "—"),
            html.Td(str(s.n_frames)),
            html.Td(f"{s.scanned_area_um2 / 1e6:.3g}"),
            html.Td(str(s.n_particles)),
            html.Td(f"{s.number_density_per_mm2:.3g}" if s.number_density_per_mm2 is not None else "—"),
        ]) for s in scans]
        return html.Div([html.Div(title, className="psd-table-title"),
                         html.Table([head] + body, className="comp-table")])

    return [stats_table, scans_table("Size bands", r.bands), scans_table("Frames searched by EDS", r.eds_scans)]


def _natural_key(text: str) -> List[Any]:
    """Sort key putting e.g. "S_900C" before "S_1000C"."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


def _compared_sizes(sample_dirs, samples) -> Tuple[List[Tuple[str, be.ParticleSizeData]], List[str]]:
    """Particle sizes of the samples to compare, labelled as in the sample list, and the samples that failed."""
    labels = {o["value"]: o["label"] for o in samples or []}
    items, errors = [], []
    for sample_dir in sorted(sample_dirs or [], key=lambda d: _natural_key(labels.get(d, d))):
        label = labels.get(sample_dir, os.path.basename(sample_dir))
        try:
            data = get_particle_sizes(sample_dir)
        except Exception as exc:
            errors.append(f"{label}: {exc}")
            continue
        if data is not None:
            items.append((label, data))
    return items, errors


def _compared_results(sample_dir) -> List[Tuple[str, be.ParticleSizeData]]:
    """Statistics results of a sample, labelled (oldest first)."""
    info = get_info(sample_dir) if sample_dir else None
    if info is None:
        return []
    items = []
    for result in info.ledger.size_results:
        data = get_particle_sizes(sample_dir, result.result_id)
        if data is not None:
            items.append((be.size_result_label(info.ledger, result), data))
    return items


def _comparison_banner(items, errors) -> List[Any]:
    out = [html.Div(f"{len(items)} sample{'s' if len(items) != 1 else ''} compared. Number distributions count "
                    "particles; volume distributions weight them by their volume (d³).", className="psd-source")]
    eds_only = [label for label, data in items if data.source == "eds_frames"]
    if eds_only:
        out.append(html.Div("From frames searched by EDS only, missing the particles too small for them: "
                            + ", ".join(eds_only) + ".", className="hint"))
    out += [html.Div(f"Cannot compute {e}", className="err") for e in errors]
    return out


def _comparison_table(items, first_column: str = "Sample") -> List[Any]:
    if not items:
        return []
    df = be.size_comparison_table(items).rename(columns={"Sample": first_column})
    for col, _key, ci_key in be.COMPARISON_STATS:
        if ci_key:
            df[f"{col} 95 % CI"] = [f"{lo:.3g}–{hi:.3g}" if pd.notna(lo) and pd.notna(hi) else ""
                                    for lo, hi in zip(df[f"{col} CI low"], df[f"{col} CI high"])]
    df["Saved"] = np.where(df["Saved"], "yes", "computed")
    shown = [first_column, "Source", "Particles", "D10", "D50", "D50 95 % CI", "D90", "Span (D90-D10)/D50",
             "D[4,3]", "Dv50", "Saved"]
    numeric = {"D10", "D50", "D90", "Span (D90-D10)/D50", "D[4,3]", "Dv50"}
    columns = [{"name": c, "id": c, **({"type": "numeric", "format": {"specifier": ".3~g"}} if c in numeric else {})}
               for c in shown]
    return [html.Div([
        html.Div("Equivalent diameter (µm)", className="psd-table-title"),
        dash_table.DataTable(
            data=df[shown].to_dict("records"), columns=columns, sort_action="native",
            style_as_list_view=True, style_table={"overflowX": "auto"},
            style_cell={"fontSize": 13, "padding": "4px 8px", "textAlign": "right"},
            style_cell_conditional=[{"if": {"column_id": c}, "textAlign": "left"} for c in (first_column, "Source")],
            style_header={"fontWeight": 600, "backgroundColor": "#f3f5f8"},
        ),
    ], className="psd-compare-table")]


def _spectrum_info(data: be.AnalysisData, sid: str, fit: Optional[Dict[str, Any]]) -> List[Any]:
    rows = data.comps[data.comps["spectrum"].astype(str) == str(sid)]
    if rows.empty:
        return [html.Div("Spectrum not found", className="hint")]
    r = rows.iloc[0]
    status = r["status"]
    cl = r["cluster"]
    badge_style = {"background": pl.cluster_color(cl)} if status == "clustered" and pd.notna(cl) else {}
    badge = f"Cluster {int(cl)}" if status == "clustered" and pd.notna(cl) else status
    head = html.Div([
        html.Span(badge, className="badge", style=badge_style),
        html.Span(f"Particle {r['particle']}" if r["particle"] not in (None, "", "None") and pd.notna(r["particle"]) else "",
                  className="muted"),
    ], className="spectrum-badges")
    comp_rows = [html.Tr([html.Th("")] + [html.Th(el) for el in data.elements])]
    comp_rows.append(html.Tr([html.Td("at%")] + [html.Td(_fmt(r.get(f"{el} at%"))) for el in data.elements]))
    comp_rows.append(html.Tr([html.Td("w%")] + [html.Td(_fmt(r.get(f"{el} w%"))) for el in data.elements]))
    if fit and fit.get("comp_at"):
        comp_rows.append(html.Tr([html.Td("refit at%", className="muted")] + [
            html.Td(_fmt((fit["comp_at"].get(el) or 0) * 100), className="muted") for el in data.elements]))
    metrics = html.Div([
        _chip("an. err", f"{_fmt(r['an_err'])} w%"),
        _chip("flag", pl.flag_label(r["quant_flag"])),
        _chip("R²", _fmt(r["r_squared"], 5)),
        _chip("χ²ᵣ", _fmt(r["redchi_sq"])),
        _chip("counts", f"{int(r['total_counts']):,}" if pd.notna(r["total_counts"]) else "—"),
    ], className="metrics")
    out = [head, html.Table(comp_rows, className="comp-table"), metrics]
    if isinstance(r["comment"], str) and r["comment"]:
        out.append(html.Div(r["comment"], className="comment"))
    return out


# =============================================================================
# Helpers
# =============================================================================
def _default_axes(data: be.AnalysisData, current: List[Optional[str]]) -> List[str]:
    """Keep the current axes if valid; else use the elements spreading most between clusters."""
    els = data.detectable or data.elements
    cur = [c for c in current if c in els]
    if len(cur) == len(current) and len(set(cur)) == len(cur) and len(cur) == min(3, len(els)):
        return cur
    if len(els) <= 3:
        return list(els)
    X = data.comps[data.comps["status"] != "not quantified"][els].to_numpy(dtype=float)
    spread = data.centroids[:, [data.elements.index(e) for e in els]].std(axis=0) if len(data.centroids) > 1 \
        else np.nanstd(X, axis=0) if len(X) else np.zeros(len(els))
    top = sorted(np.argsort(-np.nan_to_num(spread), kind="stable")[:3])
    return [els[i] for i in top]


# =============================================================================
# App
# =============================================================================
def _log_callback_error(err: Exception) -> None:
    """Print callback errors to the terminal (Dash only reports them in the browser)."""
    import traceback

    traceback.print_exception(type(err), err, err.__traceback__)
    raise err


def create_app(results_dir: Optional[str] = None) -> Dash:
    app = Dash(
        __name__,
        title="AutoEMX",
        assets_folder=str(Path(__file__).with_name("assets")),
        suppress_callback_exceptions=False,
        on_error=_log_callback_error,
    )
    app.layout = build_layout(results_dir or "")
    allowed_roots: Dict[str, bool] = {}

    @app.server.route("/analysis-file")
    def analysis_file():  # pragma: no cover - exercised through the browser
        from flask import request

        path = os.path.realpath(request.args.get("path", ""))
        if not path.lower().endswith(".png") or not _allowed(path):
            abort(404)
        if not os.path.exists(path):
            abort(404)
        return send_file(path, mimetype="image/png", max_age=0)

    def _allowed(path: str) -> bool:
        return any(path.startswith(os.path.realpath(root) + os.sep) for root in allowed_roots)

    _png_cache: Dict[Tuple[str, float, int], bytes] = {}

    @app.server.route("/sem-image")
    def sem_image():  # pragma: no cover - exercised through the browser
        """SEM image as PNG (TIFFs converted), downscaled to ``max`` pixels."""
        from flask import request

        path = os.path.realpath(request.args.get("path", ""))
        if not path.lower().endswith(be.SEM_IMAGE_EXTENSIONS) or not _allowed(path) or not os.path.exists(path):
            abort(404)
        max_size = int(request.args.get("max", 0) or 0)
        key = (path, _mtime(path), max_size)
        png = _png_cache.get(key)
        if png is None:
            png = be.read_image_png(path, max_size or None)[0]
            if len(_png_cache) > 400:
                _png_cache.clear()
            _png_cache[key] = png
        return Response(png, mimetype="image/png")

    # ---------------------------------------------------------------- pages
    @app.callback(
        Output("page-acq", "hidden"),
        Output("page-quant", "hidden"),
        Output("page-analysis", "hidden"),
        Output("page-single", "hidden"),
        Input("main-tabs", "value"),
    )
    def show_page(page):
        # Pages stay in the layout (hidden), so they keep their state when switching
        return tuple(page != value for value, _ in PAGES)

    register_clear_buttons(app)

    # ---------------------------------------------------------------- periodic table of quantifiable elements
    @app.callback(
        Output("pt-modal", "hidden"),
        Output("pt-micro", "value"),
        Output("pt-kv", "value"),
        Output("pt-elements", "data"),
        Input("a-pt-btn", "n_clicks"),
        Input("q-pt-btn", "n_clicks"),
        Input("pt-close", "n_clicks"),
        State(_pid("amicro.microscope_ID"), "value"),
        State(_pid("aacq.beam_energy"), "value"),
        State({"type": "a-cell", "uid": ALL, "col": ALL}, "value"),
        State("q-table", "selected_rows"),
        State("q-table", "data"),
        State("sample-dd", "value"),
        prevent_initial_call=True,
    )
    def open_periodic_table(_a, _q, _c, micro, kv, _cells, q_selected, q_rows, sample_dir):
        trig = ctx.triggered_id
        if trig == "pt-close" or not ctx.triggered[0].get("value"):
            return True, no_update, no_update, no_update
        elements: List[str] = []
        if trig == "a-pt-btn":
            # Microscope and beam energy of the acquisition settings; elements of the samples to acquire
            for cell in ctx.states_list[2]:
                if cell["id"]["col"] == "els":
                    try:
                        elements += be._elements_list(cell.get("value"))
                    except ValueError:
                        pass
        else:
            # Microscope and beam energy of the ticked samples (or of the current sample)
            dirs = [q_rows[i]["id"] for i in (q_selected or []) if q_rows and i < len(q_rows)] or \
                ([sample_dir] if sample_dir else [])
            summaries, _n = SUMMARIES.get(dirs)
            micro, kv = be.dflt.microscope_ID, 15.0
            if summaries:
                first = next(iter(summaries.values()))
                micro, kv = first.get("microscope") or micro, first.get("beam_energy_keV") or kv
                for summ in summaries.values():
                    elements += [e.strip() for e in (summ.get("elements") or "").split(",") if e.strip()]
        return False, micro or be.dflt.microscope_ID, float(kv) if kv else 15.0, sorted(set(elements))

    @app.callback(
        Output("pt-kv", "options"),
        Output("pt-table", "children"),
        Output("pt-info", "children"),
        Input("pt-micro", "value"),
        Input("pt-kv", "value"),
        Input("pt-elements", "data"),
    )
    def show_periodic_table(micro, kv, elements):
        energies = be.standards_beam_energies(micro) if micro else []
        values = sorted(set(energies) | ({float(kv)} if kv else set()))
        options = [{"label": f"{e:g} kV" + ("" if e in energies else " (no standards)"), "value": e} for e in values]
        if not micro or not kv:
            return options, periodic_table(None), ""
        available = be.quantifiable_elements(micro, kv)
        elements = elements or []
        if available is None:
            info = html.Span(f"No P/B standards for {micro} at {kv:g} kV: spectra can be fitted but not quantified. "
                             f"Standards are available at: {', '.join(f'{e:g} kV' for e in energies) or 'none'}.",
                             className="err")
        else:
            missing = [e for e in elements if e not in available and e not in be._DEFAULT_UNDETECTABLE_ELS]
            info = [html.Span(f"{len(available)} elements quantifiable with the standards of {micro} at {kv:g} kV "
                              f"({be.dflt.measurement_mode} mode).")]
            if missing:
                info.append(html.Span(f" No standard for the sample element(s): {', '.join(missing)}.", className="err"))
            elif elements:
                info.append(html.Span(" All sample elements can be quantified.", className="ok"))
        # Elements undetectable by EDS (e.g. Li) are not quantified anyway: not marked
        detectable = [e for e in elements if e not in be._DEFAULT_UNDETECTABLE_ELS]
        return options, periodic_table(available, detectable), info
    tab_acquisition.register(app)
    tab_quantification.register(app)
    tab_single.register(app)

    # ---------------------------------------------------------------- folder / samples
    @app.callback(
        Output("folder", "value"),
        Input("browse-btn", "n_clicks"),
        prevent_initial_call=True,
    )
    def browse(_):
        path = _pick_folder()
        return path if path else no_update

    @app.callback(
        Output("sample-dd", "options"),
        Output("sample-dd", "value"),
        Output("scan-msg", "children"),
        Input("folder", "value"),
        Input("scan-btn", "n_clicks"),
        State("sample-dd", "value"),
    )
    def scan(folder, _, current):
        if not folder:
            return [], None, "Choose the folder containing your AutoEMX sample folders."
        try:
            samples = be.find_samples(folder)
        except FileNotFoundError as exc:
            return [], None, html.Span(str(exc), className="err")
        allowed_roots[os.path.realpath(os.path.expanduser(folder))] = True
        if not samples:
            return [], None, html.Span("No sample folders (with ledger.json) found here.", className="err")
        options = [{"label": s.label, "value": s.sample_dir} for s in samples]
        dirs = [s.sample_dir for s in samples]
        value = current if current in dirs else dirs[0]
        return options, value, f"{len(samples)} sample{'s' if len(samples) != 1 else ''}"

    # ---------------------------------------------------------------- sample / analyses
    @app.callback(
        Output("sample-info", "children"),
        Output("analysis-dd", "options"),
        Output("analysis-dd", "value"),
        Output(_pid("plot.els_to_plot"), "options"),
        Output(_pid("plot.els_excluded_clust_plot"), "options"),
        Output("sem-tab", "label"),
        Output("show-image-btn", "children"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
    )
    def sample_changed(sample_dir, _):
        try:
            info = get_info(sample_dir)
        except Exception as exc:
            return html.Span(f"Cannot read ledger: {exc}", className="err"), [], None, [], [], no_update, no_update
        if info is None:
            return [], [], None, [], [], no_update, no_update
        # Microscope type saved in the ledger (e.g. "SEM") names the image tab.
        em_type = str(getattr(info.ledger.configs.microscope_cfg, "type", "") or "SEM")
        options = [{"label": a.label, "value": a.key} for a in reversed(info.analyses)]
        # After a run (or on sample change) show the active analysis, i.e. the latest one run.
        active = info.active_analysis
        value = active.key if active else None
        els = [{"label": e, "value": e} for e in info.elements]
        return (_sample_info_view(info), options, value, els, els, f"{em_type} images",
                f"Show in {em_type} image")

    # Sub-tabs showing EDS data: greyed out for samples without spectra (e.g. only particle sizes)
    eds_tabs = ["plot-tab", "clusters-tab", "spectra-tab", "dist-tab", "figures-tab", "settings-tab"]

    @app.callback(
        *[Output(t, "disabled") for t in eds_tabs],
        Output("psd-tab", "disabled"),
        Output("tabs", "value", allow_duplicate=True),
        Output("analysis-sidebar", "className"),
        Output("spectrum-panel", "className"),
        Output("analysis-dd", "placeholder"),
        Output("set-active-btn", "disabled"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
        State("tabs", "value"),
        prevent_initial_call=True,  # the sample is set by the folder scan
    )
    def grey_out_missing_data(sample_dir, _, tab):
        try:
            info = get_info(sample_dir)
        except Exception:
            info = None
        has_eds = info is None or info.has_eds
        has_sizes = info is not None and info.has_sizes
        disabled = {t.removesuffix("-tab"): not has_eds for t in eds_tabs}
        disabled["psd"] = not has_sizes
        if disabled.get(tab):
            tab = "plot" if has_eds else ("psd" if has_sizes else "sem")
        return (*[disabled[t.removesuffix("-tab")] for t in eds_tabs], disabled["psd"], tab,
                "sidebar" + ("" if has_eds else " no-eds"), "spectrum-panel" + ("" if has_eds else " no-eds"),
                "No analysis yet" if has_eds else "No EDS spectra", not has_eds)

    @app.callback(
        Output("psd-compare", "options"),
        Output("psd-compare", "value"),
        Input("sample-dd", "options"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
    )
    def compare_choices(samples, sample_dir, _):
        # Samples of the results folder with particle sizes; by default those in the folder of the shown sample
        options = sorted((o for o in samples or [] if has_sizes(o["value"])), key=lambda o: _natural_key(o["label"]))
        parent = os.path.dirname(sample_dir) if sample_dir else None
        value = [o["value"] for o in options if os.path.dirname(o["value"]) == parent]
        return options, value

    @app.callback(
        Output("psd-banner", "children"),
        Output("psd-graph", "figure"),
        Output("psd-tables", "children"),
        Output("psd-compare-controls", "hidden"),
        Output("psd-compare-samples", "hidden"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
        Input("psd-mode", "value"),
        Input("psd-compare", "value"),
        Input("psd-kind", "value"),
        Input("psd-basis", "value"),
        Input("psd-result-dd", "value"),
        State("sample-dd", "options"),
    )
    def particle_sizes(sample_dir, _, mode, compared, kind, basis, result_id, samples):
        if mode == "compare":
            items, errors = _compared_sizes(compared, samples)
            return (_comparison_banner(items, errors), pl.psd_comparison_figure(items, kind, basis),
                    _comparison_table(items), False, False)
        if mode == "results":
            items = _compared_results(sample_dir)
            banner = [html.Div(f"{len(items)} result{'s' if len(items) != 1 else ''} of this sample compared: "
                               "other segmentations, or frames and particles left out.", className="psd-source")]
            return (banner, pl.psd_comparison_figure(items, kind, basis), _comparison_table(items, "Result"), False,
                    True)  # the samples list does not apply: hidden
        try:
            data = get_particle_sizes(sample_dir, result_id)
        except Exception as exc:
            return html.Div(f"Cannot compute the particle size statistics: {exc}", className="err"), \
                pl.psd_figure(None), [], True, True
        if data is None:
            return (html.Div("No particle size data for this sample.", className="hint"), pl.psd_figure(None), [],
                    True, True)
        return _psd_banner(data, get_info(sample_dir)), pl.psd_figure(data), _psd_tables(data), True, True

    @app.callback(
        Output("psd-download", "data"),
        Input("psd-export-btn", "n_clicks"),
        State("psd-compare", "value"),
        State("sample-dd", "options"),
        State("psd-mode", "value"),
        State("sample-dd", "value"),
        prevent_initial_call=True,
    )
    def export_comparison(_, compared, samples, mode, sample_dir):
        if mode == "results":
            items, name = _compared_results(sample_dir), f"{Path(sample_dir).name}_size_results.csv"
        else:
            items, name = _compared_sizes(compared, samples)[0], "Particle_size_comparison.csv"
        if not items:
            return no_update
        return dcc.send_data_frame(be.size_comparison_table(items).to_csv, name, index=False)

    @app.callback(Output("psd-frames-section", "hidden"), Output("psd-result-controls", "hidden"),
                  Input("psd-mode", "value"))
    def frames_section_shown(mode):
        return mode != "sample", mode != "sample"

    @app.callback(
        Output("psd-result-dd", "options"),
        Output("psd-result-dd", "value"),
        Output("psd-seg-model", "value"),
        Output("psd-seg-thresh", "value"),
        Output("psd-seg-px", "value"),
        Output("psd-seg-params", "value"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
    )
    def result_choices(sample_dir, _):
        # Results of the sample (the active one shown); segmentation settings of the active one to start from
        info = get_info(sample_dir) if sample_dir else None
        if info is None:
            return [], None, no_update, no_update, no_update, no_update
        ledger = info.ledger
        result = ledger.size_result
        segmentation = ledger.segmentation(result.segmentation_id if result is not None else None) or ledger.segmentation(
            result.eds_segmentation_id if result is not None else None)
        if segmentation is None:
            cfg = ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig()
            model, threshold, px = cfg.par_segmentation_model, cfg.par_brightness_thresh, cfg.min_par_pixels
            params = cfg.par_segmentation_params
        else:
            model, threshold, px = (segmentation.par_segmentation_model, segmentation.par_brightness_thresh,
                                    segmentation.min_par_pixels)
            params = segmentation.model_params
        return (be.size_result_options(info), ledger.active_size_result,
                model if model in PAR_SEGMENTATION_MODELS else no_update, threshold, px,
                json.dumps(params) if params else "")

    @app.callback(Output("psd-seg-thresh-box", "hidden"), Output("psd-seg-params-box", "hidden"),
                  Output("psd-seg-note", "children"), Input("psd-seg-model", "value"))
    def threshold_used(model):
        # Fields of the selected model: a brightness threshold for the threshold model, options for the others;
        # how to install a model not installed
        return (model != PAR_SEGMENTATION_MODELS[0], model == PAR_SEGMENTATION_MODELS[0],
                be.segmentation_install_note(model))

    @app.callback(
        Output("psd-result-msg", "children"),
        Output("data-version", "data", allow_duplicate=True),
        Input("psd-set-active-btn", "n_clicks"),
        State("sample-dd", "value"),
        State("psd-result-dd", "value"),
        State("data-version", "data"),
        prevent_initial_call=True,
    )
    def set_active_result(_, sample_dir, result_id, version):
        info = get_info(sample_dir) if sample_dir else None
        if info is None:
            return html.Span("Select a sample first.", className="err"), no_update
        try:
            message = be.set_active_size_result(info, result_id)
        except Exception as exc:
            return html.Span(str(exc), className="err"), no_update
        _info_cache.pop(sample_dir, None)
        return html.Span(message, className="ok"), (version or 0) + 1

    @app.callback(
        Output("psd-job-store", "data", allow_duplicate=True),
        Output("psd-poll", "disabled", allow_duplicate=True),
        Output("psd-recalc-msg", "children", allow_duplicate=True),
        Input("psd-reseg-btn", "n_clicks"),
        State("sample-dd", "value"),
        State("psd-seg-model", "value"),
        State("psd-seg-thresh", "value"),
        State("psd-seg-px", "value"),
        State("psd-seg-params", "value"),
        prevent_initial_call=True,
    )
    def segment_again(_, sample_dir, model, threshold, min_pixels, params_text):
        if not sample_dir:
            return no_update, no_update, html.Span("Select a sample first.", className="err")
        if be.segmentation_install_hint(model):
            return no_update, no_update, html.Span(be.segmentation_install_note(model), className="err")
        params = {}
        if model != PAR_SEGMENTATION_MODELS[0] and (params_text or "").strip():
            try:
                params = json.loads(params_text)
            except ValueError:
                params = None
            if not isinstance(params, dict):
                return no_update, no_update, html.Span(
                    'Options must be a JSON object, e.g. {"pred_iou_thresh": 0.7}.', className="err")
        payload = {"par_segmentation_model": model, "par_brightness_thresh": threshold, "min_par_pixels": min_pixels,
                   "par_segmentation_params": params}
        try:
            job = JOBS.start("psd_reseg", sample_dir, payload, f"Segmentation of {Path(sample_dir).name}")
        except RuntimeError as exc:
            return no_update, no_update, html.Span(str(exc), className="err")
        return ({"id": job.job_id, "done": False}, False,
                html.Span("Segmenting the frames again from their saved images…", className="running"))

    @app.callback(
        Output("psd-frame-dd", "options"),
        Output("psd-frame-dd", "value"),
        Output("psd-frame-graph", "figure"),
        Output("psd-frame-caption", "children"),
        Output("psd-frame-out", "value"),
        Output("psd-excl", "data"),
        Output("psd-pending", "children"),
        Output("psd-recalc-btn", "disabled"),
        Output("psd-recalc-btn", "children"),
        Output("psd-reset-btn", "disabled"),
        Input("sample-dd", "value"),
        Input("data-version", "data"),
        Input("psd-frame-dd", "value"),
        Input("psd-frame-prev", "n_clicks"),
        Input("psd-frame-next", "n_clicks"),
        Input("psd-frame-out", "value"),
        Input("psd-frame-graph", "clickData"),
        Input("psd-reset-btn", "n_clicks"),
        Input("psd-result-dd", "value"),
        State("psd-excl", "data"),
    )
    def size_frames_view(sample_dir, _version, frame_id, _prev, _next, frame_out, click, _reset, result_id, excl):
        # Frames and particles to leave out ("psd-excl") are edited here, then saved by "Recalculate and save"
        trigger = ctx.triggered_id
        try:
            data = get_particle_sizes(sample_dir, result_id)
        except Exception:
            data = None
        if data is None:
            return ([], None, pl.size_frame_figure(None, None, []), "", [], None, "", True, "Recalculate and save",
                    True)
        frames = be.size_frames(get_info(sample_dir), data)
        saved = {"frames": sorted(data.result.excluded_frame_ids),
                 "particles": sorted([str(f), int(i)] for f, i in data.result.excluded_particles)}
        shown = data.result.result_id if data.stored else None
        if (trigger in (None, "sample-dd", "data-version", "psd-reset-btn", "psd-result-dd")
                or (excl or {}).get("sample") != sample_dir or (excl or {}).get("result") != shown):
            excl = {"sample": sample_dir, "result": shown, **saved}
        ids = [f.frame_id for f in frames]
        if frame_id not in ids:
            frame_id = ids[0] if ids else None
        if trigger in ("psd-frame-prev", "psd-frame-next") and ids:
            frame_id = ids[(ids.index(frame_id) + (1 if trigger == "psd-frame-next" else -1)) % len(ids)]
        frame = next((f for f in frames if f.frame_id == frame_id), None)
        # Particles are left out one by one, as (frame_id, index in the frame)
        out_frames, out_particles = set(excl["frames"]), {(str(f), int(i)) for f, i in excl["particles"]}
        if trigger == "psd-frame-out" and frame is not None:
            (out_frames.add if "out" in (frame_out or []) else out_frames.discard)(frame.frame_id)
        if trigger == "psd-frame-graph" and click and click.get("points"):
            pid = click["points"][0].get("customdata")
            if frame is not None and pid is not None and 0 <= int(pid) < frame.n_particles:
                key = (frame.frame_id, int(pid))
                (out_particles.discard if key in out_particles else out_particles.add)(key)
        excl = {"sample": sample_dir, "result": shown, "frames": sorted(out_frames),
                "particles": sorted(map(list, out_particles))}

        def option(f):
            state = "left out" if f.frame_id in out_frames else ("overlapping" if f.status == "overlapping" else "")
            label = f"{f.frame_id} · {f.scan_label.split(' (')[0]} · {f.n_particles} particles"
            return {"label": label + (f" · {state}" if state else ""), "value": f.frame_id}

        if frame is not None and frame.image_path:
            from PIL import Image

            with Image.open(frame.image_path) as im:
                size = im.size
            figure = pl.size_frame_figure(
                f"/sem-image?max=1600&path={quote(frame.image_path)}", size, frame.marks,
                {i for f, i in out_particles if f == frame.frame_id}, frame.is_eds,
                frame.frame_id in out_frames, uirevision=frame.image_path,
                outlines=load_outlines(sample_dir, frame.segmentation_id, frame.frame_id),
            )
        else:
            figure = pl.size_frame_figure(None, None, [])
        caption = ""
        if frame is not None:
            caption = (f"{frame.scan_label} · {frame.n_particles} particles · {frame.area_um2:.4g} µm²"
                       + (" · overlaps a larger frame: not counted" if frame.status == "overlapping" else ""))
        changed = excl["frames"] != saved["frames"] or excl["particles"] != saved["particles"]
        left_out = f"Left out: {len(out_frames)} frame(s), {len(out_particles)} particle(s)"
        pending = html.Span(left_out + (" — not saved: recalculate to apply" if changed else ""),
                            className="warn" if changed else "muted")
        can_save = be.can_save_particle_sizes(data) and (changed or not data.stored)
        label = "Save statistics" if not data.stored and not changed else "Recalculate and save"
        return ([option(f) for f in frames], frame_id, figure, caption,
                ["out"] if frame is not None and frame.frame_id in out_frames else [], excl, pending, not can_save,
                label, not changed)

    @app.callback(
        Output("psd-job-store", "data"),
        Output("psd-poll", "disabled"),
        Output("psd-recalc-msg", "children"),
        Input("psd-recalc-btn", "n_clicks"),
        State("sample-dd", "value"),
        State("psd-excl", "data"),
        prevent_initial_call=True,
    )
    def recalculate_sizes(_, sample_dir, excl):
        if not sample_dir or not excl or excl.get("sample") != sample_dir:
            return no_update, no_update, no_update
        payload = {"frame_ids": excl["frames"], "particles": excl["particles"], "result_id": excl.get("result")}
        try:
            job = JOBS.start("psd", sample_dir, payload, f"Particle size statistics of {Path(sample_dir).name}")
        except RuntimeError as exc:
            return no_update, no_update, html.Span(str(exc), className="err")
        return {"id": job.job_id, "done": False}, False, html.Span("Recalculating…", className="running")

    @app.callback(
        Output("psd-recalc-msg", "children", allow_duplicate=True),
        Output("psd-job-store", "data", allow_duplicate=True),
        Output("psd-poll", "disabled", allow_duplicate=True),
        Output("data-version", "data", allow_duplicate=True),
        Input("psd-poll", "n_intervals"),
        State("psd-job-store", "data"),
        State("data-version", "data"),
        prevent_initial_call=True,
    )
    def poll_sizes_job(_, job_data, version):
        job = JOBS.get(job_data["id"]) if job_data else None
        if job is None or job_data.get("done"):
            return no_update, no_update, True, no_update
        res = job.result()
        if res is None:
            doing = "Segmenting the frames again" if job.kind == "psd_reseg" else "Recalculating"
            return html.Span(f"{doing}… {job.elapsed():.0f} s", className="running"), no_update, False, no_update
        _info_cache.pop(job.sample_dir, None)
        msg = (html.Span(f"Saved in the sample ledger ({job.elapsed():.0f} s).", className="ok") if res.get("ok")
               else html.Span(f"Not saved: {res.get('error')}", className="err"))
        return msg, {**job_data, "done": True}, True, (version or 0) + 1

    @app.callback(
        Output({"type": "param", "key": ALL}, "value"),
        Input("sample-dd", "value"),
        Input("load-settings-btn", "n_clicks"),
        State("analysis-dd", "value"),
    )
    def load_params(sample_dir, _, analysis_key):
        outputs = ctx.outputs_list
        info = get_info(sample_dir) if sample_dir else None
        if info is None:
            return [no_update] * len(outputs)
        ref = be.find_analysis(info, analysis_key) if ctx.triggered_id == "load-settings-btn" else None
        values = be.sample_param_values(info, ref)
        return [ui_value(o["id"], values) for o in outputs]

    @app.callback(
        Output("summary", "children"),
        Output("clusters-view", "children"),
        Output("spectra-table", "data"),
        Output("spectra-table", "columns"),
        Output("spectra-table", "style_data_conditional"),
        Output("figures-view", "children"),
        Output("settings-view", "children"),
        Output("dist-graph", "figure"),
        Output("folder-store", "data"),
        Output("axis-x", "options"),
        Output("axis-y", "options"),
        Output("axis-z", "options"),
        Output("axis-x", "value"),
        Output("axis-y", "value"),
        Output("axis-z", "value"),
        Output("plot-mode", "value"),
        Output("mode-forced", "data"),
        Input("analysis-dd", "value"),
        Input("data-version", "data"),
        State("sample-dd", "value"),
        State("axis-x", "value"),
        State("axis-y", "value"),
        State("axis-z", "value"),
        State("plot-mode", "value"),
        State("mode-forced", "data"),
    )
    def analysis_changed(key, _, sample_dir, ax, ay, az, mode, mode_forced):
        empty = ([], [html.Div("Select a sample.", className="hint")], [], [], [], [], [], pl.distribution_figure(None),
                 None, [], [], [], None, None, None, mode, mode_forced)
        try:
            info, data = get_analysis(sample_dir, key)
        except Exception as exc:
            msg = html.Div(f"Cannot load this analysis: {exc}", className="err")
            return ([msg],) + empty[1:]
        if data is None:
            return empty
        records, columns, styles = _table_payload(data)
        settings = [html.Pre(data.config_summary or "No Analysis_config_summary.txt for this analysis.",
                             className="settings-pre")]
        els = data.detectable or data.elements
        opts = [{"label": e, "value": e} for e in els]
        axes = _default_axes(data, [ax, ay, az][:min(3, len(els))])
        axes = (axes + [None, None, None])[:3]
        # Samples with 2 elements can only be shown in 2D; restore 3D afterwards if 2D was forced.
        if len(els) < 3:
            mode, mode_forced = "2d", mode != "2d" or bool(mode_forced)
        elif mode_forced:
            mode, mode_forced = "3d", False
        return (
            _summary_view(info, data), _clusters_view(data), records, columns, styles,
            _figures_view(data), settings, pl.distribution_figure(data), data.folder or info.sample_dir,
            opts, opts, opts, axes[0], axes[1], axes[2], mode, mode_forced,
        )

    # ---------------------------------------------------------------- clustering plot
    @app.callback(
        Output("axis-z", "disabled"),
        Input("plot-mode", "value"),
    )
    def z_disabled(mode):
        return mode == "2d"

    @app.callback(
        Output("cluster-graph", "figure"),
        Input("plot-mode", "value"),
        Input("axis-x", "value"),
        Input("axis-y", "value"),
        Input("axis-z", "value"),
        Input("color-by", "value"),
        Input("plot-opts", "value"),
        Input("mix-conf", "value"),
        Input("selected-spectrum", "data"),
        Input("analysis-dd", "value"),
        Input("data-version", "data"),
        Input("zoom-on", "value"),
        State("sample-dd", "value"),
    )
    def cluster_plot(mode, ax, ay, az, color_by, opts, mix_conf, selected, key, version, zoom, sample_dir):
        try:
            _, data = get_analysis(sample_dir, key)
        except Exception:
            data = None
        axes = [ax, ay] if mode == "2d" else [ax, ay, az]
        # A new zoom choice resets the view; other changes (e.g. clicking a point) keep the user's zoom
        uirev = f"{sample_dir}|{key}|{mode}|{'-'.join(map(str, axes))}|{version}|{zoom}"
        min_conf = float(mix_conf) if mix_conf not in (None, "") else pl.DEFAULT_MIN_MIXTURE_CONF
        return pl.clustering_figure(data, axes, mode, opts, color_by, highlight=selected, uirevision=uirev,
                                    min_mixture_conf=min_conf, zoom=zoom or pl.ZOOM_FULL)

    @app.callback(
        Output("zoom-on", "options"),
        Output("zoom-on", "value"),
        Input("analysis-dd", "value"),
        Input("data-version", "data"),
        State("sample-dd", "value"),
        State("zoom-on", "value"),
    )
    def zoom_choices(key, _, sample_dir, current):
        try:
            _, data = get_analysis(sample_dir, key)
        except Exception:
            data = None
        options = pl.zoom_options(data)
        keep = current in {o["value"] for o in options}
        return options, current if keep else pl.ZOOM_FULL

    # ---------------------------------------------------------------- spectrum selection
    @app.callback(
        Output("fit-msg", "children", allow_duplicate=True),
        Input("selected-spectrum", "data"),
        prevent_initial_call=True,
    )
    def clear_fit_msg(_):
        return ""

    @app.callback(
        Output("selected-spectrum", "data"),
        Input("cluster-graph", "clickData"),
        Input("dist-graph", "clickData"),
        Input("sem-graph", "clickData"),
        Input("spectra-table", "active_cell"),
        Input("prev-btn", "n_clicks"),
        Input("next-btn", "n_clicks"),
        Input("sample-dd", "value"),
        State("selected-spectrum", "data"),
        State("analysis-dd", "value"),
        prevent_initial_call=True,
    )
    def select_spectrum(click, dist_click, sem_click, cell, _p, _n, sample_dir, current, key):
        trig = ctx.triggered_id
        if trig == "sample-dd":
            return None
        if trig == "cluster-graph":
            return _point_spectrum(click) or no_update
        if trig == "dist-graph":
            return _point_spectrum(dist_click) or no_update
        if trig == "sem-graph":
            return _point_spectrum(sem_click) or no_update
        if trig == "spectra-table":
            return str(cell["row_id"]) if cell and cell.get("row_id") is not None else no_update
        _, data = get_analysis(sample_dir, key)
        if data is None or data.comps.empty:
            return no_update
        ids = data.comps["spectrum"].astype(str).tolist()
        if current not in ids:
            return ids[0]
        i = ids.index(current) + (1 if trig == "next-btn" else -1)
        return ids[i % len(ids)]

    @app.callback(
        Output("spectrum-title", "children"),
        Output("spectrum-info", "children"),
        Output("spectrum-graph", "figure"),
        Input("selected-spectrum", "data"),
        Input("fit-result", "data"),
        Input("log-y", "value"),
        Input("analysis-dd", "value"),
        State("sample-dd", "value"),
    )
    def spectrum_panel(sid, fit, log_y, key, sample_dir):
        if not sid or not sample_dir:
            return "No spectrum selected", [], pl.spectrum_figure(None)
        info, data = get_analysis(sample_dir, key)
        if info is None or data is None:
            return "No spectrum selected", [], pl.spectrum_figure(None)
        fit = fit if fit and fit.get("sample_dir") == sample_dir and str(fit.get("spectrum_id")) == str(sid) else None
        try:
            raw = be.load_raw_spectrum(info, sid)
        except Exception as exc:
            return f"Spectrum {sid}", [html.Div(str(exc), className="err")], pl.spectrum_figure(None)
        title = f"Spectrum {sid}"
        fig = pl.spectrum_figure(raw, fit, log_y=bool(log_y), title=f"{info.sample_id} · {title}" + (" · fitted" if fit else ""))
        return title, _spectrum_info(data, sid, fit), fig

    # ---------------------------------------------------------------- greyed-out parameters
    @app.callback(
        Output({"type": "param-section", "key": ALL}, "className"),
        Output({"type": "param-row", "key": ALL}, "className"),
        Input(_pid("clust.method"), "value"),
        Input(_pid("clust.geometry"), "value"),
        Input(_pid("clust.k_forced"), "value"),
        Input(_pid("clust.auto_merge_clusters"), "value"),
        Input(_pid("clust.do_matrix_decomposition"), "value"),
    )
    def grey_out_unused(method, geometry, k_forced, auto_merge, do_mixture):
        values = {
            "clust.method": method, "clust.geometry": geometry, "clust.k_forced": k_forced,
            "clust.auto_merge_clusters": bool(auto_merge), "clust.do_matrix_decomposition": bool(do_mixture),
        }
        sections, params = active_params(values)
        section_classes = [
            "param-section" + ("" if sections.get(o["id"]["key"], True) else " disabled")
            for o in ctx.outputs_list[0]
        ]
        row_classes = []
        for o in ctx.outputs_list[1]:
            spec = _spec(o["id"])
            row_classes.append(_row_class(spec) + ("" if params.get(spec.key, True) else " disabled"))
        return section_classes, row_classes

    # ---------------------------------------------------------------- SEM images
    @app.callback(
        Output("sem-gallery", "children"),
        Output("sem-paths", "data"),
        Input("analysis-dd", "value"),
        Input("data-version", "data"),
        State("sample-dd", "value"),
    )
    def sem_gallery(key, _, sample_dir):
        if not sample_dir:
            return [], []
        images = be.find_sem_images(sample_dir)
        if not images:
            return [html.Div(f"No SEM images in {os.path.join(sample_dir, be.cnst.IMAGES_DIR)}.", className="hint")], []
        _, data = get_analysis(sample_dir, key)
        items = []
        for i, im in enumerate(images):
            caption = [html.B(im.label)]
            if data is not None and im.kind == "spots":
                ids = be.spectra_on_image(data.comps, im)["spectrum"].astype(str).tolist()
                if ids:
                    caption.append(html.Span(" · spectra " + ", ".join(ids), className="muted"))
            if im.kind == "other":
                caption = [html.Span(im.name, className="muted")]
            items.append(html.Div(
                [html.Img(src=f"/sem-image?max=420&path={quote(im.path)}", className="saved-img"),
                 html.Div(caption, className="sem-thumb-caption")],
                id={"type": "sem-thumb", "index": i}, n_clicks=0, className="saved-fig sem-thumb",
                title="Show in the viewer",
            ))
        return items, [im.path for im in images]

    @app.callback(
        Output("sem-image", "data"),
        Input("selected-spectrum", "data"),
        Input({"type": "sem-thumb", "index": ALL}, "n_clicks"),
        Input("sem-paths", "data"),
        Input("sem-prev-btn", "n_clicks"),
        Input("sem-next-btn", "n_clicks"),
        State("sample-dd", "value"),
        State("analysis-dd", "value"),
        State("sem-image", "data"),
    )
    def choose_sem_image(sid, clicks, paths, _prev, _next, sample_dir, key, current):
        trig = ctx.triggered_id
        if trig in ("sem-prev-btn", "sem-next-btn"):
            if not paths:
                return no_update
            i = paths.index(current) if current in paths else -1
            return paths[(i + (1 if trig == "sem-next-btn" else -1)) % len(paths)]
        if isinstance(trig, dict) and trig.get("type") == "sem-thumb":
            i = trig["index"]
            if not clicks or i >= len(clicks) or not clicks[i] or not paths or i >= len(paths):
                return no_update
            return paths[i]
        if not sample_dir:
            return None
        images = be.find_sem_images(sample_dir)
        if sid:
            _, data = get_analysis(sample_dir, key)
            row = data.comps[data.comps["spectrum"].astype(str) == str(sid)] if data is not None else None
            if row is not None and not row.empty:
                im = be.sem_image_for_spectrum(images, row.iloc[0]["particle"], row.iloc[0]["frame"])
                if im is not None:
                    return im.path
        if current and paths and current in paths:
            return no_update
        spots = [im for im in images if im.kind == "spots"]
        return (spots or images)[0].path if images else None

    @app.callback(
        Output("sem-graph", "figure"),
        Output("sem-caption", "children"),
        Output("sem-counter", "children"),
        Output({"type": "sem-thumb", "index": ALL}, "className"),
        Input("sem-image", "data"),
        Input("selected-spectrum", "data"),
        Input("sem-spots", "value"),
        Input("analysis-dd", "value"),
        Input("data-version", "data"),
        State("sample-dd", "value"),
        State("sem-paths", "data"),
    )
    def sem_view(path, sid, show_spots, key, _, sample_dir, paths):
        classes = ["saved-fig sem-thumb" + (" active" if p == path else "") for p in (paths or [])]
        n_thumbs = len(ctx.outputs_list[3])
        classes = classes[:n_thumbs] + ["saved-fig sem-thumb"] * (n_thumbs - len(classes))
        counter = f"{paths.index(path) + 1} / {len(paths)}" if paths and path in paths else ""
        if not path or not sample_dir or not os.path.exists(path):
            return pl.sem_image_figure(None, None), "No images for this sample.", counter, classes
        images = be.find_sem_images(sample_dir)
        image = next((im for im in images if im.path == path), None)
        _, data = get_analysis(sample_dir, key)
        spots = be.spectra_on_image(data.comps, image) if (data is not None and image is not None) else None
        from PIL import Image

        with Image.open(path) as im:
            size = im.size
        fig = pl.sem_image_figure(
            f"/sem-image?max=1600&path={quote(path)}", size, spots, data, selected=sid,
            show_spots=bool(show_spots), uirevision=path,
        )
        caption = [html.B(image.label if image else Path(path).name)]
        if spots is not None and len(spots):
            caption.append(html.Span(f" · {len(spots)} spectra · click a spot to select it", className="muted"))
        if sid and (spots is None or str(sid) not in spots["spectrum"].astype(str).tolist()):
            caption.append(html.Span(f" · spectrum {sid} is not in this image", className="warn"))
        return fig, caption, counter, classes

    @app.callback(
        Output("single-request", "data"),
        Output("main-tabs", "value"),
        Input("to-single-btn", "n_clicks"),
        State("selected-spectrum", "data"),
        State("sample-dd", "value"),
        prevent_initial_call=True,
    )
    def open_in_single(n, sid, sample_dir):
        if not sid or not sample_dir:
            return no_update, no_update
        return {"sample_dir": sample_dir, "spectrum_id": str(sid), "n": n}, "single"

    @app.callback(
        Output("tabs", "value"),
        Input("show-image-btn", "n_clicks"),
        State("selected-spectrum", "data"),
        prevent_initial_call=True,
    )
    def show_spectrum_image(_, sid):
        # The image of the selected spectrum is already chosen by choose_sem_image.
        return "sem" if sid else no_update

    app.clientside_callback(
        """
        function(classes, tab) {
            // Keep the image shown in the viewer visible in the thumbnail strip (also right after
            // switching to the tab, when the strip has just been drawn).
            [80, 400].forEach(function(delay) {
                setTimeout(function() {
                    const el = document.querySelector('.sem-thumb.active');
                    if (el) { el.scrollIntoView({behavior: 'smooth', block: 'nearest', inline: 'center'}); }
                }, delay);
            });
            return window.dash_clientside.no_update;
        }
        """,
        Output("sem-gallery", "title"),
        Input({"type": "sem-thumb", "index": ALL}, "className"),
        Input("tabs", "value"),
    )

    @app.callback(
        Output("action-msg", "children", allow_duplicate=True),
        Output("data-version", "data", allow_duplicate=True),
        Input("set-active-btn", "n_clicks"),
        State("sample-dd", "value"),
        State("analysis-dd", "value"),
        State("data-version", "data"),
        prevent_initial_call=True,
    )
    def set_active(_, sample_dir, analysis_key, version):
        info = get_info(sample_dir) if sample_dir else None
        if info is None:
            return html.Span("Select a sample first.", className="err"), no_update
        try:
            message = be.set_active_analysis(info, analysis_key)
        except Exception as exc:
            return html.Span(f"Could not change the active analysis: {exc}", className="err"), no_update
        _info_cache.pop(sample_dir, None)
        return html.Span(message, className="ok"), (version or 0) + 1

    @app.callback(
        Output("action-msg", "children", allow_duplicate=True),
        Input("sem-open-btn", "n_clicks"),
        State("sem-image", "data"),
        prevent_initial_call=True,
    )
    def open_sem_image(_, path):
        if path and os.path.exists(path):
            _open_path(path)
        return no_update

    # ---------------------------------------------------------------- jobs
    @app.callback(
        Output("job-store", "data"),
        Output("poll", "disabled"),
        Output("run-msg", "children"),
        Input("run-btn", "n_clicks"),
        State({"type": "param", "key": ALL}, "value"),
        State("sample-dd", "value"),
        prevent_initial_call=True,
    )
    def run(_, __, sample_dir):
        if not sample_dir:
            return no_update, no_update, html.Span("Select a sample first.", className="err")
        try:
            values = be.coerce_values(form_values(ctx.states_list[0]))
        except ValueError as exc:
            return no_update, no_update, html.Div([html.B("Invalid parameters: "), str(exc)], className="err")
        payload = {"quantify": False, "analyse": True, "analysis_kwargs": be.analysis_kwargs(values)}
        desc = f"Analysis of {Path(sample_dir).name}"
        try:
            job = JOBS.start("analysis", sample_dir, payload, desc)
        except RuntimeError as exc:
            return no_update, no_update, html.Span(str(exc), className="err")
        return {"id": job.job_id, "done": False}, False, html.Span(f"{desc} started…", className="running")

    @app.callback(
        Output("fit-job-store", "data"),
        Output("poll", "disabled", allow_duplicate=True),
        Output("fit-msg", "children"),
        Input("fit-btn", "n_clicks"),
        State("selected-spectrum", "data"),
        State("sample-dd", "value"),
        prevent_initial_call=True,
    )
    def fit(_, sid, sample_dir):
        # Fitted with the settings of the sample's active quantification
        if not sid or not sample_dir:
            return no_update, no_update, "Select a spectrum first."
        try:
            job = JOBS.start("fit", sample_dir, {"spectrum_id": sid}, f"Fit of spectrum {sid}")
        except RuntimeError as exc:
            return no_update, no_update, str(exc)
        return {"id": job.job_id, "done": False}, False, f"Fitting spectrum {sid}…"

    @app.callback(
        Output("run-msg", "children", allow_duplicate=True),
        Output("run-log", "children"),
        Output("run-btn", "disabled"),
        Output("cancel-btn", "disabled"),
        Output("data-version", "data"),
        Output("job-store", "data", allow_duplicate=True),
        Output("fit-job-store", "data", allow_duplicate=True),
        Output("fit-result", "data"),
        Output("fit-msg", "children", allow_duplicate=True),
        Output("poll", "disabled", allow_duplicate=True),
        Output("log-details", "open"),
        Output("fit-cancel-btn", "hidden"),
        Input("poll", "n_intervals"),
        State("job-store", "data"),
        State("fit-job-store", "data"),
        State("data-version", "data"),
        prevent_initial_call=True,
    )
    def poll(_, job_data, fit_data, version):
        run_msg = log = data_version = job_out = fit_job_out = fit_result = fit_msg = no_update
        run_disabled = cancel_disabled = log_open = fit_cancel_hidden = no_update
        any_running = False

        job = JOBS.get(job_data["id"]) if job_data else None
        if job is not None and not job_data.get("done"):
            log = job.log_tail()
            res = job.result()
            if res is None:
                any_running = True
                run_msg = html.Span(f"{job.description} running… {job.elapsed():.0f} s", className="running")
                run_disabled, cancel_disabled, log_open = True, False, True
            else:
                job_out = {**job_data, "done": True}
                run_disabled, cancel_disabled = False, True
                if res.get("ok"):
                    text = f"{job.description} finished in {job.elapsed():.0f} s."
                    if res.get("n_clusters") is not None:
                        text += f" {res['n_clusters']} cluster(s)."
                    run_msg = html.Span(text, className="ok")
                    if res.get("warning"):
                        run_msg = html.Span(text + " " + res["warning"], className="warn")
                    log_open = False
                else:
                    run_msg = html.Span(f"{job.description} failed: {res.get('error')}", className="err")
                _info_cache.pop(job.sample_dir, None)
                data_version = (version or 0) + 1

        fjob = JOBS.get(fit_data["id"]) if fit_data else None
        if fjob is not None and not fit_data.get("done"):
            res = fjob.result()
            fit_cancel_hidden = res is not None
            if res is None:
                any_running = True
                fit_msg = f"Fitting… {fjob.elapsed():.0f} s"
            else:
                fit_job_out = {**fit_data, "done": True}
                if fjob.cancelled:
                    fit_msg = "Fit cancelled."
                elif res.get("ok"):
                    fit_result = {**res, "sample_dir": fjob.sample_dir}
                    r2 = res.get("r_squared")
                    fit_msg = f"Fitted in {fjob.elapsed():.0f} s" + (f" · R² {r2:.5f}" if r2 is not None else "")
                else:
                    fit_msg = f"Fit failed: {res.get('error')}"
        return (run_msg, log, run_disabled, cancel_disabled, data_version, job_out, fit_job_out, fit_result,
                fit_msg, not any_running, log_open, fit_cancel_hidden)

    @app.callback(
        Output("run-msg", "children", allow_duplicate=True),
        Input("cancel-btn", "n_clicks"),
        State("job-store", "data"),
        prevent_initial_call=True,
    )
    def cancel(_, job_data):
        if job_data:
            JOBS.cancel(job_data["id"])
        return html.Span("Cancelling…", className="warn")

    @app.callback(
        Output("fit-msg", "children", allow_duplicate=True),
        Input("fit-cancel-btn", "n_clicks"),
        State("fit-job-store", "data"),
        prevent_initial_call=True,
    )
    def fit_cancel(_, fit_data):
        if fit_data:
            JOBS.cancel(fit_data["id"])
        return "Cancelling…"

    # ---------------------------------------------------------------- actions
    @app.callback(
        Output("action-msg", "children"),
        Input("open-folder-btn", "n_clicks"),
        Input("save-html-btn", "n_clicks"),
        State("folder-store", "data"),
        State("cluster-graph", "figure"),
        State("plot-mode", "value"),
        prevent_initial_call=True,
    )
    def actions(_o, _s, folder, figure, mode):
        if not folder or not os.path.isdir(folder):
            return html.Span("No analysis folder.", className="err")
        if ctx.triggered_id == "open-folder-btn":
            _open_path(folder)
            return no_update
        import plotly.graph_objects as go

        path = os.path.join(folder, f"Clustering_plot_interactive_{mode}.html")
        go.Figure(figure).write_html(path, include_plotlyjs=True, config=_GRAPH_CONFIG)
        return html.Span(f"Saved {path}", className="ok")

    return app
