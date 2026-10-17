#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Interactive Plotly figures of the AutoEMX GUI."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from autoemx.core.composition_analysis.clustering_plot_axes import (
    CLUSTERING_3D_VIEW_AZIM,
    CLUSTERING_3D_VIEW_ELEV,
)
from autoemx.gui.backend import QUANT_FLAG_MEANINGS, AnalysisData, ParticleSizeData

CLUSTER_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
]
STATUS_STYLE = {
    "noise": dict(color="#b0b0b0", name="Noise (unclustered)"),
    "discarded": dict(color="#9a9a9a", name="Discarded spectra"),
    "not analysed": dict(color="#4a6fa5", name="Quantified spectra"),
}
REF_COLOR = "#1a237e"
CENTROID_COLOR = "#111111"
HIGHLIGHT_COLOR = "#ff1744"

# Plot options (values of the "show" checklist)
OPT_DISCARDED = "discarded"
OPT_STD = "std"
OPT_REFS = "refs"
OPT_MIXTURES = "mixtures"
OPT_CENTROIDS = "centroids"
DEFAULT_OPTIONS = [OPT_DISCARDED, OPT_STD, OPT_REFS, OPT_MIXTURES, OPT_CENTROIDS]
# Zoom of the clustering plot (values of the "Zoom on" menu): whole 0-100 % range, all spectra, or one
# cluster (ZOOM_CLUSTER + its number)
ZOOM_FULL = "full"
ZOOM_DATA = "data"
ZOOM_CLUSTER = "cl"
# Best mixtures with a lower confidence score are not drawn (unless the threshold is lowered in the GUI)
DEFAULT_MIN_MIXTURE_CONF = 0.8


def cluster_color(i: int) -> str:
    return CLUSTER_COLORS[int(i) % len(CLUSTER_COLORS)]


def _empty_figure(message: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font=dict(size=15, color="#666"),
                       xref="paper", yref="paper", x=0.5, y=0.5)
    fig.update_layout(xaxis=dict(visible=False), yaxis=dict(visible=False),
                      template="plotly_white", margin=dict(l=10, r=10, t=10, b=10))
    return fig


def _hover_template(data: AnalysisData) -> str:
    """Hover text: spectrum, particle, composition in both units, quality metrics."""
    lines = ["<b>Spectrum %{customdata[0]}</b>", "Particle: %{customdata[1]}", "Status: %{customdata[2]}"]
    n = 3
    comp_parts = []
    for el in data.elements:
        comp_parts.append(f"{el} %{{customdata[{n}]:.1f}}")
        n += 1
    lines.append(f"{data.unit}: " + ", ".join(comp_parts))
    lines.append(f"An. error: %{{customdata[{n}]:.1f}} w%")
    lines.append(f"Quant flag: %{{customdata[{n + 1}]}}")
    lines.append(f"R²: %{{customdata[{n + 2}]:.5f}}")
    return "<br>".join(lines) + "<extra></extra>"


def _customdata(data: AnalysisData, df: pd.DataFrame) -> np.ndarray:
    cols = [df["spectrum"].astype(str), df["particle"].astype(str), df["status"]]
    cols += [df[el] * 100 for el in data.elements]
    cols += [df["an_err"], df["quant_flag"].astype(str), df["r_squared"]]
    return np.column_stack([np.asarray(c, dtype=object) for c in cols])


def _project(values: np.ndarray, mode: str) -> List[np.ndarray]:
    """Fractions (n, len(axes)) to plot coordinates: percent, or normalised ternary components."""
    values = np.atleast_2d(np.asarray(values, dtype=float))
    if mode == "ternary":
        total = values.sum(axis=1, keepdims=True)
        total[total == 0] = np.nan
        values = values / total
    return [values[:, i] * 100 for i in range(values.shape[1])]


def _point_trace(plot_mode: str, coords: List[np.ndarray], **kw) -> Any:
    """Trace for the plot type; ``kw`` may hold the trace's own ``mode`` (markers, lines, ...)."""
    if plot_mode == "3d":
        return go.Scatter3d(x=coords[0], y=coords[1], z=coords[2], **kw)
    if plot_mode == "ternary":
        return go.Scatterternary(a=coords[2], b=coords[0], c=coords[1], **kw)
    return go.Scatter(x=coords[0], y=coords[1], **kw)


def _marker(mode: str, color: Any, symbol: str = "circle", size: Optional[int] = None, **kw) -> Dict[str, Any]:
    if size is None:
        size = 4 if mode == "3d" else 8
    if mode == "3d" and symbol not in ("circle", "circle-open", "cross", "diamond", "diamond-open",
                                       "square", "square-open", "x"):
        symbol = {"triangle-up": "diamond-open", "star": "diamond"}.get(symbol, "circle")
    marker = dict(color=color, symbol=symbol, size=size, **kw)
    return marker


def _camera_eye(elev_deg: float, azim_deg: float, distance: float = 2.2) -> Dict[str, float]:
    """Plotly camera position matching a matplotlib ``view_init(elev, azim)``."""
    elev, azim = np.radians(elev_deg), np.radians(azim_deg)
    return dict(
        x=float(distance * np.cos(elev) * np.cos(azim)),
        y=float(distance * np.cos(elev) * np.sin(azim)),
        z=float(distance * np.sin(elev)),
    )


def _ellipsoid(center: np.ndarray, radii: np.ndarray, n: int = 24):
    u = np.linspace(0, 2 * np.pi, n)
    v = np.linspace(0, np.pi, n)
    x = center[0] + radii[0] * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radii[1] * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radii[2] * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z


def clustering_figure(
    data: Optional[AnalysisData],
    axes: Sequence[str],
    mode: str = "3d",
    options: Iterable[str] = DEFAULT_OPTIONS,
    color_by: str = "cluster",
    highlight: Optional[str] = None,
    uirevision: Optional[str] = None,
    min_mixture_conf: float = DEFAULT_MIN_MIXTURE_CONF,
    zoom: str = ZOOM_FULL,
) -> go.Figure:
    """Clustering plot: 3D (3 elements), ternary (3 elements, normalised) or 2D (2 elements).

    *zoom*: axis ranges on the whole 0-100 % range (``ZOOM_FULL``), on all spectra (``ZOOM_DATA``), or on
    one cluster (``f"{ZOOM_CLUSTER}{i}"``).
    """
    if data is None:
        return _empty_figure("Select a sample")
    options = set(options or [])
    axes = [a for a in axes if a in data.elements]
    needed = 2 if mode == "2d" else 3
    if len(axes) != needed or len(set(axes)) != needed:
        return _empty_figure(f"Choose {needed} different elements for the {mode.upper()} plot")
    comps = data.comps
    quantified = comps[comps["status"] != "not quantified"]
    if quantified.empty:
        return _empty_figure("No quantified spectra. Run the quantification first.")
    idx = [data.elements.index(a) for a in axes]
    unit = data.unit
    hover = _hover_template(data)
    fig = go.Figure()

    def add_points(df: pd.DataFrame, name: str, marker: Dict[str, Any], legendgroup: str, **kw) -> None:
        if df.empty:
            return
        fig.add_trace(_point_trace(
            mode, _project(df[axes].to_numpy(dtype=float), mode),
            mode="markers", name=name, legendgroup=legendgroup, marker=marker,
            customdata=_customdata(data, df), hovertemplate=hover, **kw,
        ))

    show_discarded = OPT_DISCARDED in options
    if color_by == "cluster":
        if show_discarded:
            add_points(quantified[quantified["status"] == "discarded"], STATUS_STYLE["discarded"]["name"],
                       _marker(mode, STATUS_STYLE["discarded"]["color"], "triangle-up", opacity=0.55),
                       "discarded")
        add_points(quantified[quantified["status"] == "not analysed"], STATUS_STYLE["not analysed"]["name"],
                   _marker(mode, STATUS_STYLE["not analysed"]["color"]), "not analysed")
        add_points(quantified[quantified["status"] == "noise"], STATUS_STYLE["noise"]["name"],
                   _marker(mode, STATUS_STYLE["noise"]["color"]), "noise")
        clustered = quantified[quantified["status"] == "clustered"]
        for cl in sorted(clustered["cluster"].dropna().unique()):
            df = clustered[clustered["cluster"] == cl]
            add_points(df, f"Cluster {int(cl)} ({len(df)})",
                       _marker(mode, cluster_color(cl), line=dict(width=0.5, color="white") if mode != "3d" else None),
                       f"cl{int(cl)}")
    else:
        df = quantified if show_discarded else quantified[quantified["status"] != "discarded"]
        if color_by == "an_err":
            values, title, scale = df["an_err"].astype(float), "An. error (w%)", "RdBu"
            cmid = 0.0
        elif color_by == "particle":
            values, title, scale, cmid = pd.to_numeric(df["particle"], errors="coerce"), "Particle #", "Turbo", None
        elif color_by == "quant_flag":
            values, title, scale, cmid = pd.to_numeric(df["quant_flag"], errors="coerce"), "Quant flag", "Viridis", None
        else:
            values, title, scale, cmid = pd.to_numeric(df["r_squared"], errors="coerce"), "R²", "Viridis", None
        marker = _marker(mode, values, colorscale=scale, showscale=True,
                         colorbar=dict(title=title, thickness=14, len=0.6))
        if cmid is not None:
            marker["cmid"] = cmid
        add_points(df, "Spectra", marker, "spectra")

    # Standard deviation around each centroid
    centroids = data.centroids[:, idx] if len(data.centroids) else np.zeros((0, len(axes)))
    stdevs = data.stdevs[:, idx] if len(data.stdevs) else np.zeros((0, len(axes)))
    if OPT_STD in options and mode != "ternary":
        for i, (c, s) in enumerate(zip(centroids, stdevs)):
            if np.any(np.isnan(s)):
                continue
            color = cluster_color(i)
            if mode == "3d":
                x, y, z = _ellipsoid(c * 100, s * 100)
                fig.add_trace(go.Surface(
                    x=x, y=y, z=z, opacity=0.18, showscale=False, hoverinfo="skip",
                    colorscale=[[0, color], [1, color]], name=f"Cluster {i} ±1σ",
                    legendgroup=f"cl{i}", showlegend=False,
                ))
            else:
                t = np.linspace(0, 2 * np.pi, 80)
                fig.add_trace(go.Scatter(
                    x=c[0] * 100 + s[0] * 100 * np.cos(t), y=c[1] * 100 + s[1] * 100 * np.sin(t),
                    mode="lines", fill="toself", line=dict(color=color, dash="dash", width=1),
                    opacity=0.25, hoverinfo="skip", legendgroup=f"cl{i}", showlegend=False,
                ))

    # Mixture lines: best-ranked mixture of each cluster
    ref_comps = data.ref_comps
    if OPT_MIXTURES in options and len(ref_comps) and data.mixtures:
        first = True
        for i, mixes in enumerate(data.mixtures):
            if not mixes:
                continue
            best = mixes[0]
            refs = [r for r in best.get("refs", []) if r in ref_comps.index]
            conf = best.get("conf_score")
            if len(refs) < 2 or conf is None or conf < min_mixture_conf:
                continue
            pts = ref_comps.loc[refs, axes].to_numpy(dtype=float)
            if len(refs) > 2:
                pts = np.vstack([pts, pts[:1]])
            label = f"Cluster {i}: {' + '.join(refs)} (conf. {conf:.2f})"
            fig.add_trace(_point_trace(
                mode, _project(pts, mode), mode="lines",
                line=dict(color=cluster_color(i), dash="dash", width=4 if mode == "3d" else 2),
                name="Best mixtures" if first else label, legendgroup="mixtures", showlegend=first,
                hovertemplate=label + "<extra></extra>",
            ))
            first = False

    if OPT_CENTROIDS in options and len(centroids):
        text = [f"C{i}" for i in range(len(centroids))]
        hover_c = []
        for i, c in enumerate(data.centroids):
            n = data.n_points[i] if i < len(data.n_points) else "?"
            comp = ", ".join(f"{el} {v * 100:.1f}" for el, v in zip(data.elements, c))
            sd = ", ".join(f"{el} {v * 100:.1f}" for el, v in zip(data.elements, data.stdevs[i])) if i < len(data.stdevs) else ""
            hover_c.append(f"<b>Cluster {i}</b> ({n} spectra)<br>{unit}: {comp}<br>σ: {sd}")
        fig.add_trace(_point_trace(
            mode, _project(centroids, mode), mode="markers+text", text=text,
            textposition="top center", textfont=dict(size=13, color=CENTROID_COLOR),
            marker=_marker(mode, CENTROID_COLOR, "x", size=6 if mode == "3d" else 13),
            name="Centroids", legendgroup="centroids", hovertext=hover_c, hoverinfo="text",
        ))

    if OPT_REFS in options and len(ref_comps):
        pts = ref_comps[axes].to_numpy(dtype=float)
        hover_r = [
            f"<b>{f}</b><br>{unit}: " + ", ".join(f"{el} {v * 100:.1f}" for el, v in zip(data.elements, row))
            for f, row in zip(ref_comps.index, ref_comps[data.elements].to_numpy(dtype=float))
        ]
        fig.add_trace(_point_trace(
            mode, _project(pts, mode), mode="markers+text", text=list(ref_comps.index),
            textposition="top center" if mode == "3d" else "top right", textfont=dict(size=13, color=REF_COLOR),
            marker=_marker(mode, REF_COLOR, "star" if mode != "3d" else "diamond", size=6 if mode == "3d" else 15),
            name="Candidate phases", legendgroup="refs", hovertext=hover_r, hoverinfo="text",
            **({} if mode == "3d" else {"cliponaxis": False}),
        ))

    if highlight is not None:
        sel = comps[comps["spectrum"].astype(str) == str(highlight)]
        sel = sel[sel["status"] != "not quantified"]
        if not sel.empty:
            fig.add_trace(_point_trace(
                mode, _project(sel[axes].to_numpy(dtype=float), mode), mode="markers",
                marker=_marker(mode, HIGHLIGHT_COLOR, "circle-open", size=12 if mode == "3d" else 20,
                               line=dict(color=HIGHLIGHT_COLOR, width=3)),
                name=f"Spectrum {highlight}", hoverinfo="skip", showlegend=True,
            ))

    _layout_clustering(fig, data, axes, mode, unit, zoom, uirevision, OPT_REFS in options)
    return fig


def zoom_options(data: Optional[AnalysisData]) -> List[Dict[str, str]]:
    """Choices of the "Zoom on" menu for this analysis."""
    options = [{"label": "Whole range", "value": ZOOM_FULL}, {"label": "All spectra", "value": ZOOM_DATA}]
    if data is not None:
        for i in range(len(data.centroids)):
            n = data.n_points[i] if i < len(data.n_points) else None
            options.append({"label": f"Cluster {i}" + (f" ({n})" if n is not None else ""),
                            "value": f"{ZOOM_CLUSTER}{i}"})
    return options


def _zoom_points(data: AnalysisData, axes: Sequence[str], zoom: str, show_refs: bool, mode: str) -> np.ndarray:
    """Points (percent) the axes are zoomed onto: all spectra and centroids, or the spectra and centroid of
    one cluster."""
    comps = data.comps[data.comps["status"] != "not quantified"]
    idx = [data.elements.index(a) for a in axes]
    if zoom.startswith(ZOOM_CLUSTER) and zoom[len(ZOOM_CLUSTER):].isdigit():
        i = int(zoom[len(ZOOM_CLUSTER):])
        members = comps[(comps["status"] == "clustered") & (pd.to_numeric(comps["cluster"], errors="coerce") == i)]
        pts = members[list(axes)].to_numpy(dtype=float) * 100
        if i < len(data.centroids):
            pts = np.vstack([pts, data.centroids[i:i + 1, idx] * 100])
        if len(pts):
            return pts
    pts = comps[list(axes)].to_numpy(dtype=float) * 100
    if len(data.centroids):
        pts = np.vstack([pts, data.centroids[:, idx] * 100])
    if mode == "3d" and show_refs and len(data.ref_comps):
        pts = np.vstack([pts, data.ref_comps[list(axes)].to_numpy(dtype=float) * 100])
    return pts


def _zoom_range(values: np.ndarray, pad: float = 0.08) -> List[float]:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return [0, 100]
    lo, hi = float(vals.min()), float(vals.max())
    span = max(hi - lo, 0.5)
    return [max(0.0, lo - pad * span), min(100.0, hi + pad * span)]


# Gap (percent) between the 0 / 100 % planes and the walls of the 3D box, so that markers and labels
# at 0 % (e.g. candidate phases without one of the elements) are not hidden by the walls.
_WALL_GAP_3D = 3.0


def _range_3d(values: np.ndarray, zoom: bool) -> List[float]:
    """Axis range of the 3D plot, leaving a gap between the data and the walls."""
    if not zoom:
        return [-_WALL_GAP_3D, 100 + _WALL_GAP_3D]
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return [-_WALL_GAP_3D, 100 + _WALL_GAP_3D]
    lo, hi = float(vals.min()), float(vals.max())
    gap = 0.15 * max(hi - lo, 0.5)  # room for the labels above the markers
    return [max(-_WALL_GAP_3D, lo - gap), min(100 + _WALL_GAP_3D, hi + gap)]


def _ticks(lo: float, hi: float, n: int = 6) -> List[float]:
    """Round tick values within [lo, hi], restricted to 0-100 %."""
    lo, hi = max(lo, 0.0), min(hi, 100.0)
    step = next((s for s in (0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 25) if (hi - lo) / s <= n), 25)
    start = np.ceil(lo / step) * step
    return [round(float(v), 2) for v in np.arange(start, hi + 1e-9, step)]


def _ternary_zoom_mins(pts: np.ndarray, pad: float = 0.15) -> List[float]:
    """Axis minima (percent) of a ternary diagram zoomed onto the normalised points, with a margin of
    *pad* times their spread."""
    pts = np.asarray(pts, dtype=float)
    total = pts.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        norm = pts / total * 100
    norm = norm[np.all(np.isfinite(norm), axis=1)]
    if norm.size == 0:
        return [0.0, 0.0, 0.0]
    spread = max(float((norm.max(axis=0) - norm.min(axis=0)).max()), 0.5)
    mins = np.clip(norm.min(axis=0) - pad * spread, 0, None)
    if mins.sum() > 99.5:  # keep a visible triangle
        mins *= 99.5 / mins.sum()
    return [float(np.floor(m * 10) / 10) for m in mins]


def _layout_clustering(fig: go.Figure, data: AnalysisData, axes: Sequence[str], mode: str, unit: str,
                       zoom: str, uirevision: Optional[str], show_refs: bool = False) -> None:
    title = f"{data.config.method} clustering · {data.sample_id}"
    pts = _zoom_points(data, axes, zoom or ZOOM_FULL, show_refs, mode)
    zoom = (zoom or ZOOM_FULL) != ZOOM_FULL
    common = dict(
        template="plotly_white",
        title=dict(text=title, x=0.01, font=dict(size=15)),
        legend=dict(itemsizing="constant", bgcolor="rgba(255,255,255,0.75)", font=dict(size=12)),
        margin=dict(l=10, r=10, t=45, b=10),
        uirevision=uirevision,
        hoverlabel=dict(font_size=12),
    )
    if mode == "3d":
        def ax(i: int) -> Dict[str, Any]:
            rng = _range_3d(pts[:, i], zoom)
            ticks = _ticks(*rng)
            if i < 2:
                rng = rng[::-1]  # x and y reversed, as in the saved plot: (0,0,0) at the back
            # No opaque wall panes: they hide markers and labels lying on or crossing the walls.
            return dict(title=f"{axes[i]} ({unit})", range=rng, tickvals=ticks, showbackground=False,
                        gridcolor="#d3d8df", zeroline=False, showspikes=False)
        # Same starting view as the saved matplotlib clustering plot.
        eye = _camera_eye(CLUSTERING_3D_VIEW_ELEV, CLUSTERING_3D_VIEW_AZIM)
        fig.update_layout(
            scene=dict(xaxis=ax(0), yaxis=ax(1), zaxis=ax(2), aspectmode="cube",
                       camera=dict(eye=eye, up=dict(x=0, y=0, z=1)), uirevision=uirevision),
            **common,
        )
    elif mode == "ternary":
        mins = _ternary_zoom_mins(pts) if zoom else [0.0, 0.0, 0.0]
        fig.update_layout(
            ternary=dict(
                sum=100,
                aaxis=dict(title=f"{axes[2]}", ticksuffix="%", min=mins[2]),
                baxis=dict(title=f"{axes[0]}", ticksuffix="%", min=mins[0]),
                caxis=dict(title=f"{axes[1]}", ticksuffix="%", min=mins[1]),
                bgcolor="#fbfbfc",
            ),
            **common,
        )
        fig.update_layout(title_text=f"{title} · normalised {'-'.join(axes)} ({unit})",
                          margin=dict(l=60, r=60, t=70, b=50))
    else:
        rx = _zoom_range(pts[:, 0]) if zoom else [0, 100]
        ry = _zoom_range(pts[:, 1]) if zoom else [0, 100]
        fig.update_layout(
            xaxis=dict(title=f"{axes[0]} ({unit})", range=rx, zeroline=False, constrain="domain"),
            yaxis=dict(title=f"{axes[1]} ({unit})", range=ry, zeroline=False, scaleanchor="x", scaleratio=1,
                       constrain="domain"),
            **common,
        )


def distribution_figure(data: Optional[AnalysisData], elements: Optional[Sequence[str]] = None) -> go.Figure:
    """Strip/box plots of each element's fraction, by cluster."""
    if data is None:
        return _empty_figure("Select a sample")
    df = data.comps[data.comps["status"].isin(["clustered", "noise", "not analysed", "discarded"])]
    if df.empty:
        return _empty_figure("No quantified spectra")
    elements = [e for e in (elements or data.elements) if e in data.elements]
    groups = []
    for status in ("clustered",):
        for cl in sorted(df.loc[df["status"] == status, "cluster"].dropna().unique()):
            groups.append((f"Cluster {int(cl)}", df[(df["status"] == status) & (df["cluster"] == cl)], cluster_color(cl)))
    for status, style in STATUS_STYLE.items():
        part = df[df["status"] == status]
        if not part.empty:
            groups.append((style["name"], part, style["color"]))
    fig = go.Figure()
    for name, part, color in groups:
        long = part.melt(id_vars=["spectrum"], value_vars=elements, var_name="element", value_name="v")
        fig.add_trace(go.Box(
            x=long["element"], y=long["v"] * 100, name=name, marker_color=color,
            boxpoints="all", jitter=0.4, pointpos=0, marker=dict(size=4, opacity=0.7),
            line=dict(width=1), customdata=long["spectrum"],
            hovertemplate="Spectrum %{customdata}<br>%{x}: %{y:.1f}<extra>" + name + "</extra>",
        ))
    fig.update_layout(
        boxmode="group", template="plotly_white", yaxis_title=f"Fraction ({data.unit})",
        margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.08),
    )
    return fig


def spectrum_figure(
    spectrum: Optional[Dict[str, Any]],
    fit: Optional[Dict[str, Any]] = None,
    log_y: bool = False,
    title: str = "",
) -> go.Figure:
    """Raw spectrum, or data/fit/background with residuals and peak labels when ``fit`` is given."""
    if spectrum is None and fit is None:
        return _empty_figure("Click a point or a table row")
    fig = go.Figure()
    if fit is not None:
        e = np.asarray(fit["energy"], dtype=float)
        counts = np.asarray(fit["counts"], dtype=float)
        model = np.asarray(fit["fit"], dtype=float)
        bkg = np.asarray(fit["background"], dtype=float)
        fig.add_trace(go.Scatter(x=e, y=counts, mode="markers", name="Data",
                                 marker=dict(size=3, color="#1f77b4")))
        fig.add_trace(go.Scatter(x=e, y=model, mode="lines", name="Fit", line=dict(color="#ff7f0e", width=1.6)))
        fig.add_trace(go.Scatter(x=e, y=bkg, mode="lines", name="Background",
                                 line=dict(color="#2ca02c", width=1.2, dash="dash")))
        fig.add_trace(go.Scatter(x=e, y=counts - model, mode="lines", name="Residuals",
                                 line=dict(color="#7f7f7f", width=1), yaxis="y2", visible="legendonly"))
        for energy, label in fit.get("peak_labels", []) or []:
            i = int(np.argmin(np.abs(e - energy))) if e.size else 0
            y = float(max(model[i], counts[i])) if e.size else 0
            fig.add_annotation(x=energy, y=np.log10(max(y, 1)) if log_y else y, text=label.replace("_", " "),
                               showarrow=True, arrowhead=0, ay=-25, font=dict(size=11))
        # Background counts under the reference lines, as in plot_quantified_spectrum. They are taken
        # from the background without detector response, so they need not reach the plotted background.
        bars = fit.get("bckgrnd_cnts") or []
        if bars:
            base = 1.0 if log_y else 0.0
            xs, ys, text = [], [], []
            for line, energy, value in bars:
                xs += [energy, energy, None]
                ys += [base, max(value, base), None]
                text += [None, f"{line.replace('_', ' ')}: {value:.1f} background counts", None]
            fig.add_trace(go.Scatter(
                x=xs, y=ys, mode="lines", name="Background counts", line=dict(color="#d62728", width=2.5),
                hovertext=text, hoverinfo="text",
            ))
        fig.update_layout(yaxis2=dict(overlaying="y", side="right", showgrid=False, showticklabels=False,
                                      zeroline=False))
        xmax = float(e.max()) if e.size else 10
    else:
        e = np.asarray(spectrum["energy"], dtype=float)
        counts = np.asarray(spectrum["counts"], dtype=float)
        fig.add_trace(go.Scatter(x=e, y=counts, mode="lines", name="Counts", line=dict(color="#1f77b4", width=1.2)))
        nz = np.nonzero(counts > 0)[0]
        xmax = float(e[nz[-1]]) if nz.size else float(e.max())
    fig.update_layout(
        template="plotly_white",
        xaxis=dict(title="Energy (keV)", range=[0, min(xmax, 20)]),
        yaxis=dict(title="Counts", type="log" if log_y else "linear"),
        margin=dict(l=10, r=10, t=30, b=10),
        legend=dict(orientation="h", y=1.0, yanchor="bottom", x=0, font=dict(size=11)),
        uirevision=title,
    )
    return fig


def flag_label(flag: Any) -> str:
    try:
        f = int(flag)
    except (TypeError, ValueError):
        return "—"
    return f"{f} · {QUANT_FLAG_MEANINGS.get(f, 'unknown')}"


def sem_image_figure(
    image_src: Optional[str],
    size: Optional[Sequence[int]],
    spots: Optional[pd.DataFrame] = None,
    data: Optional[AnalysisData] = None,
    selected: Optional[str] = None,
    show_spots: bool = True,
    uirevision: Optional[str] = None,
) -> go.Figure:
    """SEM image (pixel coordinates) with the spectrum spots as clickable markers coloured by cluster."""
    if not image_src or not size:
        return _empty_figure("No SEM image for this spectrum")
    width, height = int(size[0]), int(size[1])
    fig = go.Figure()
    fig.add_layout_image(dict(
        source=image_src, xref="x", yref="y", x=0, y=0, sizex=width, sizey=height,
        xanchor="left", yanchor="top", sizing="stretch", layer="below",
    ))
    if show_spots and spots is not None and not spots.empty and data is not None:
        hover = _hover_template(data)
        groups = []
        clustered = spots[spots["status"] == "clustered"]
        for cl in sorted(clustered["cluster"].dropna().unique()):
            groups.append((f"Cluster {int(cl)}", clustered[clustered["cluster"] == cl], cluster_color(cl)))
        for status, style in (("noise", STATUS_STYLE["noise"]), ("discarded", STATUS_STYLE["discarded"]),
                              ("not analysed", STATUS_STYLE["not analysed"])):
            part = spots[spots["status"] == status]
            if not part.empty:
                groups.append((style["name"], part, style["color"]))
        part = spots[spots["status"] == "not quantified"]
        if not part.empty:
            groups.append(("Not quantified", part, "#ffffff"))
        for name, part, color in groups:
            outline = "#555555" if color == "#ffffff" else "white"
            fig.add_trace(go.Scatter(
                x=part["px"], y=part["py"], mode="markers", name=name,
                marker=dict(color=color, size=13, line=dict(color=outline, width=1.5), opacity=0.9),
                customdata=_customdata(data, part), hovertemplate=hover,
            ))
        if selected is not None:
            sel = spots[spots["spectrum"].astype(str) == str(selected)]
            if not sel.empty:
                fig.add_trace(go.Scatter(
                    x=sel["px"], y=sel["py"], mode="markers", name=f"Spectrum {selected}",
                    marker=dict(color="rgba(0,0,0,0)", size=26, line=dict(color=HIGHLIGHT_COLOR, width=3)),
                    hoverinfo="skip",
                ))
    fig.update_layout(
        template="plotly_white",
        xaxis=dict(range=[0, width], visible=False, constrain="domain"),
        yaxis=dict(range=[height, 0], visible=False, scaleanchor="x", scaleratio=1, constrain="domain"),
        margin=dict(l=0, r=0, t=0, b=0),
        legend=dict(orientation="h", y=-0.01, yanchor="top", x=0, font=dict(size=11)),
        uirevision=uirevision,
        dragmode="pan",
    )
    return fig


def single_spectrum_figure(
    spectrum: Optional[Dict[str, Any]],
    fit: Optional[Dict[str, Any]] = None,
    log_y: bool = False,
    zoom_line: Optional[str] = None,
    show_bckgrnd_cnts: bool = True,
    title: str = "",
    channel_lims: Optional[Sequence[int]] = None,
) -> go.Figure:
    """
    Large spectrum plot: raw spectrum, or data/fit/background with a residuals panel underneath.

    As in ``XSp_Quantifier.plot_quantified_spectrum``, only the fitted energy range is shown: the
    channels ``channel_lims`` of the raw spectrum, or the energies of the fit.
    """
    from plotly.subplots import make_subplots

    if spectrum is None and fit is None:
        return _empty_figure("Choose a spectrum")
    if fit is None:
        e = np.asarray(spectrum["energy"], dtype=float)
        counts = np.asarray(spectrum["counts"], dtype=float)
        if channel_lims is not None and len(channel_lims) == 2:
            lo, hi = max(0, int(channel_lims[0])), min(len(e), int(channel_lims[1]))
            if hi > lo:
                e, counts = e[lo:hi], counts[lo:hi]
        fig = go.Figure(go.Scatter(x=e, y=counts, mode="lines", name="Counts", line=dict(color="#1f77b4", width=1.2)))
        fig.update_layout(
            template="plotly_white",
            title=dict(text=title, x=0.01, font=dict(size=14)),
            xaxis=dict(title="Energy (keV)", range=[float(e[0]), float(e[-1])] if e.size else None),
            yaxis=dict(title="Counts", type="log" if log_y else "linear"),
            margin=dict(l=10, r=10, t=40, b=10),
            uirevision=f"{title}|raw|{log_y}",
        )
        return fig
    e = np.asarray(fit["energy"], dtype=float)
    counts = np.asarray(fit["counts"], dtype=float)
    model = np.asarray(fit["fit"], dtype=float)
    bkg = np.asarray(fit["background"], dtype=float)
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.78, 0.22], vertical_spacing=0.03)
    fig.add_trace(go.Scatter(x=e, y=counts, mode="markers", name="Data", marker=dict(size=3.5, color="#1f77b4")),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=e, y=model, mode="lines", name="Fit", line=dict(color="#ff7f0e", width=1.6)),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=e, y=bkg, mode="lines", name="Background",
                             line=dict(color="#2ca02c", width=1.3, dash="dash")), row=1, col=1)
    bars = fit.get("bckgrnd_cnts") or []
    if show_bckgrnd_cnts and bars:
        base = 1.0 if log_y else 0.0
        xs, ys, text = [], [], []
        for line, energy, value in bars:
            xs += [energy, energy, None]
            ys += [base, max(value, base), None]
            text += [None, f"{line.replace('_', ' ')}: {value:.1f} background counts", None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name="Background counts",
                                 line=dict(color="#d62728", width=2.5), hovertext=text, hoverinfo="text"),
                      row=1, col=1)
    fig.add_trace(go.Scatter(x=e, y=counts - model, mode="lines", name="Residuals", showlegend=False,
                             line=dict(color="#7f7f7f", width=1)), row=2, col=1)
    fig.add_hline(y=0, line=dict(color="#b0b6bf", width=1), row=2, col=1)
    for energy, label in fit.get("peak_labels", []) or []:
        i = int(np.argmin(np.abs(e - energy))) if e.size else 0
        y = float(max(model[i], counts[i])) if e.size else 0
        fig.add_annotation(x=energy, y=np.log10(max(y, 1)) if log_y else y, text=label.replace("_", " "),
                           showarrow=True, arrowhead=0, ay=-28, font=dict(size=12), row=1, col=1)
    x_range = [float(e.min()), float(e.max())] if e.size else None  # fitted range
    y_range = None
    peak = next((p for p in fit.get("peaks", []) if p["line"] == zoom_line), None) if zoom_line else None
    if peak and peak.get("center"):
        half = max(3 * float(peak.get("fwhm") or 0.1), 0.15)
        x_range = [peak["center"] - half, peak["center"] + half]
        sel = (e >= x_range[0]) & (e <= x_range[1])
        if sel.any() and not log_y:
            y_range = [0, float(max(counts[sel].max(), model[sel].max())) * 1.15]
    fig.update_layout(
        template="plotly_white",
        title=dict(text=title, x=0.01, font=dict(size=14)),
        margin=dict(l=10, r=10, t=40, b=10),
        legend=dict(orientation="h", y=1.02, yanchor="bottom", x=1, xanchor="right"),
        uirevision=f"{title}|{zoom_line}|{log_y}",
        hovermode="x unified",
    )
    fig.update_yaxes(title_text="Counts", type="log" if log_y else "linear", row=1, col=1,
                     **({"range": y_range} if y_range else {}))
    fig.update_yaxes(title_text="Residuals", row=2, col=1)
    fig.update_xaxes(range=x_range, row=1, col=1)
    fig.update_xaxes(title_text="Energy (keV)", range=x_range, row=2, col=1)
    return fig


PSD_COLOR = "#1f77b4"
PSD_MARK_COLOR = "#d62728"


def psd_figure(data: Optional[ParticleSizeData], bins: int = 20) -> go.Figure:
    """Particle size distribution: number fraction by equivalent diameter (weighted as the statistics), with
    logarithmic bins when the diameters span more than a factor 10, cumulative number and volume distributions,
    and D10, D50 and D90 (as ``particle_statistics.save_particle_size_histogram``)."""
    if data is None:
        return _empty_figure("No particle size data for this sample")
    d, w = np.asarray(data.diameters_um, dtype=float), np.asarray(data.weights, dtype=float)
    if d.size == 0:
        return _empty_figure("No particle with a weight in the statistics")
    d_min, d_max = float(d.min()), float(d.max())
    is_log = d_max / d_min > 10
    if np.isclose(d_min, d_max):
        edges = np.array([d_min * 0.95, d_max * 1.05])
    elif is_log:
        edges = np.geomspace(d_min, d_max, bins + 1)
    else:
        edges = np.linspace(d_min, d_max, bins + 1)
    fractions, _ = np.histogram(d, bins=edges, weights=w / w.sum() * 100)

    fig = go.Figure()
    # Histogram as a filled step line: bars keep their width on a logarithmic axis
    fig.add_trace(go.Scatter(
        x=np.repeat(edges, 2)[1:-1], y=np.repeat(fractions, 2), mode="lines", fill="tozeroy", name="Number fraction",
        line=dict(color=PSD_COLOR, width=1), fillcolor="rgba(31,119,180,0.45)", hoverinfo="skip",
    ))
    centres = np.sqrt(edges[:-1] * edges[1:]) if is_log else (edges[:-1] + edges[1:]) / 2
    fig.add_trace(go.Scatter(
        x=centres, y=fractions, mode="markers", marker=dict(size=14, opacity=0), showlegend=False,
        customdata=np.column_stack([edges[:-1], edges[1:]]),
        hovertemplate="%{customdata[0]:.3g}–%{customdata[1]:.3g} µm: %{y:.1f} %<extra></extra>",
    ))
    order = np.argsort(d, kind="stable")
    for dist_w, dash, name in ((w, "solid", "Cumulative, number"), (w * d**3, "dash", "Cumulative, volume")):
        cumulative = np.cumsum(dist_w[order]) / np.sum(dist_w) * 100
        fig.add_trace(go.Scatter(
            x=d[order], y=cumulative, mode="lines", line=dict(color="#222", width=1.4, dash=dash, shape="hv"),
            name=name, yaxis="y2", hovertemplate="%{x:.3g} µm: %{y:.1f} %<extra>" + name + "</extra>",
        ))
    stats = data.result.size_stats
    y_top = float(fractions.max()) * 1.12
    for value, label in ((stats.D10, "D10"), (stats.median, "D50"), (stats.D90, "D90")):
        fig.add_trace(go.Scatter(
            x=[value, value], y=[0, y_top], mode="lines+text", text=["", f" {label}"], textposition="bottom right",
            textfont=dict(color=PSD_MARK_COLOR, size=11), line=dict(color=PSD_MARK_COLOR, width=1.2, dash="dot"),
            showlegend=False, hovertemplate=f"{label}: {value:.3g} µm<extra></extra>",
        ))
    fig.update_layout(
        template="plotly_white", margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.08),
        xaxis=dict(title="Equivalent diameter (µm)", type="log" if is_log else "linear",
                   dtick="D2" if is_log else None),
        yaxis=dict(title="Number fraction (%)", range=[0, y_top * 1.02], rangemode="tozero"),
        yaxis2=dict(title="Cumulative (%)", overlaying="y", side="right", range=[0, 100], showgrid=False),
        hovermode="closest",
    )
    return fig


def psd_comparison_figure(items: Sequence[tuple], kind: str = "cumulative", basis: str = "number",
                          bins: int = 25) -> go.Figure:
    """Size distributions of several samples (``items``: (label, ParticleSizeData)), by number or volume:
    cumulative curves with their median marked, or histograms on common bins (logarithmic when the diameters
    span more than a factor 10)."""
    items = [(label, data) for label, data in items if data is not None and len(data.diameters_um)]
    if not items:
        return _empty_figure("Select samples with particle size data to compare")
    all_d = np.concatenate([np.asarray(data.diameters_um, dtype=float) for _, data in items])
    d_min, d_max = float(all_d.min()), float(all_d.max())
    is_log = d_max / d_min > 10
    if np.isclose(d_min, d_max):
        edges = np.array([d_min * 0.95, d_max * 1.05])
    else:
        edges = np.geomspace(d_min, d_max, bins + 1) if is_log else np.linspace(d_min, d_max, bins + 1)
    fig = go.Figure()
    for i, (label, data) in enumerate(items):
        color = cluster_color(i)
        d = np.asarray(data.diameters_um, dtype=float)
        w = np.asarray(data.weights, dtype=float) * (d**3 if basis == "volume" else 1)
        if kind == "histogram":
            fractions, _ = np.histogram(d, bins=edges, weights=w / w.sum() * 100)
            fig.add_trace(go.Scatter(
                x=np.repeat(edges, 2)[1:-1], y=np.repeat(fractions, 2), mode="lines", name=label,
                line=dict(color=color, width=1.8), hovertemplate="%{x:.3g} µm: %{y:.1f} %<extra>" + label + "</extra>",
            ))
            continue
        order = np.argsort(d, kind="stable")
        cumulative = np.cumsum(w[order]) / w.sum() * 100
        fig.add_trace(go.Scatter(
            x=d[order], y=cumulative, mode="lines", name=label, legendgroup=label,
            line=dict(color=color, width=1.8, shape="hv"),
            hovertemplate="%{x:.3g} µm: %{y:.1f} %<extra>" + label + "</extra>",
        ))
        stats = data.result.size_stats
        median = stats.Dv50 if basis == "volume" else stats.median
        if median is not None:
            fig.add_trace(go.Scatter(
                x=[median], y=[50], mode="markers", legendgroup=label, showlegend=False,
                marker=dict(color=color, size=9, line=dict(color="white", width=1)),
                hovertemplate=f"{label}<br>{'Dv50' if basis == 'volume' else 'D50'}: {median:.3g} µm<extra></extra>",
            ))
    y_title = f"{'Volume' if basis == 'volume' else 'Number'} fraction (%)"
    fig.update_layout(
        template="plotly_white", margin=dict(l=10, r=10, t=30, b=10), legend=dict(orientation="h", y=1.08),
        xaxis=dict(title="Equivalent diameter (µm)", type="log" if is_log else "linear", dtick="D2" if is_log else None),
        yaxis=dict(title=("Cumulative " + y_title[0].lower() + y_title[1:]) if kind == "cumulative" else y_title,
                   range=[0, 102] if kind == "cumulative" else None, rangemode="tozero"),
        hovermode="closest",
    )
    return fig


FRAME_MARK_COLORS = {"in": "#2ca02c", "out": "#d62728", "eds": "#17becf", "other": "#ffb000"}


def size_frame_figure(image_src: Optional[str], size: Optional[Sequence[int]], marks: Sequence[Sequence[Any]],
                      excluded_particles: Iterable[int] = (), is_eds: bool = False,
                      frame_excluded: bool = False, uirevision: Optional[str] = None,
                      outlines: Optional[Dict[str, list]] = None) -> go.Figure:
    """Frame image of a particle size measurement with its particles outlined (``outlines``, see
    autoemx.core.particle_outlines; circled if not saved): counted (green, cyan in frames searched by EDS), left out
    (red), or not counted (amber, dashed: touching the frame edge, or outside the diameters searched); ``excluded_particles``
    are their indices in the frame. Clicking the centre of a particle gives its index as customdata."""
    if not image_src or not size:
        return _empty_figure("No image saved for this frame")
    width, height = int(size[0]), int(size[1])
    excluded = set(excluded_particles)
    fig = go.Figure()
    fig.add_layout_image(dict(source=image_src, xref="x", yref="y", x=0, y=0, sizex=width, sizey=height,
                              xanchor="left", yanchor="top", sizing="stretch", layer="below",
                              opacity=0.45 if frame_excluded else 1))
    groups: Dict[str, List[Sequence[Any]]] = {"in": [], "out": [], "eds": []}
    for x, y, r, index in marks:
        key = "out" if index in excluded or frame_excluded else ("eds" if is_eds else "in")
        groups[key].append((x, y, r, index))
    names = {"in": "Counted", "out": "Left out", "eds": "Counted (frame searched by EDS)"}
    angles = np.linspace(0, 2 * np.pi, 33)
    counted_outlines = (outlines or {}).get("counted") or []

    def shape(x, y, r, index):
        # Outline of the particle if saved, else a circle of its area
        poly = counted_outlines[index] if index is not None and index < len(counted_outlines) else None
        if poly:
            xy = np.asarray(poly, dtype=float).reshape(-1, 2)
            xy = np.vstack([xy, xy[:1], [np.nan, np.nan]])
            return xy[:, 0], xy[:, 1]
        return np.append(x + r * np.cos(angles), np.nan), np.append(y + r * np.sin(angles), np.nan)

    for key, pts in groups.items():
        if not pts:
            continue
        fig.add_trace(go.Scatter(
            x=[p[0] for p in pts], y=[p[1] for p in pts], mode="markers", name=names[key], legendgroup=key,
            customdata=[p[3] for p in pts],
            marker=dict(color=FRAME_MARK_COLORS[key], size=8, opacity=0.25),  # click target, circle stays visible
            hovertemplate="Particle %{customdata}<extra></extra>",
        ))
        # Outlines of a group as one line trace (polygons separated by gaps): fast with hundreds of particles
        shapes = [shape(*p) for p in pts]
        xs, ys = np.concatenate([sx for sx, _ in shapes]), np.concatenate([sy for _, sy in shapes])
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines", name=names[key], legendgroup=key, hoverinfo="skip", showlegend=False,
            line=dict(color=FRAME_MARK_COLORS[key], width=2 if outlines is None else 1.5,
                      dash="dot" if key == "out" else "solid"),
        ))
    others = [np.asarray(p, dtype=float).reshape(-1, 2) for p in (outlines or {}).get("others") or [] if p]
    if others:
        xy = np.vstack([np.vstack([o, o[:1], [np.nan, np.nan]]) for o in others])
        fig.add_trace(go.Scatter(
            x=xy[:, 0], y=xy[:, 1], mode="lines", hoverinfo="skip",
            name="Not counted (cut by frame edge or outside specified size range)",
            line=dict(color=FRAME_MARK_COLORS["other"], width=1, dash="dash"),
        ))
    fig.update_layout(
        template="plotly_white", margin=dict(l=0, r=0, t=0, b=0), dragmode="pan", uirevision=uirevision,
        xaxis=dict(range=[0, width], visible=False, constrain="domain"),
        yaxis=dict(range=[height, 0], visible=False, scaleanchor="x", scaleratio=1, constrain="domain"),
        legend=dict(orientation="h", y=-0.01, yanchor="top", x=0, font=dict(size=11)),
    )
    return fig
