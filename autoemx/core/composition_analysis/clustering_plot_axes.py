#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared axis formatting helpers for composition clustering plots."""

from __future__ import annotations

from typing import Any, List, Optional, Sequence

import numpy as np
from matplotlib.ticker import FuncFormatter, Locator, MaxNLocator

CLUSTERING_3D_VIEW_ELEV = 15.0
CLUSTERING_3D_VIEW_AZIM = -135.0


def _clamp_fraction_limits(low: float, high: float) -> tuple[float, float]:
    """Clamp composition axis limits to the valid [0, 1] fraction range."""
    low = max(0.0, min(1.0, float(low)))
    high = max(0.0, min(1.0, float(high)))
    if high < low:
        low, high = high, low
    return low, high


def _clamp_axis_limits(ax: Any, axis_name: str) -> None:
    """Ensure a single axis never displays values outside [0, 1]."""
    get_lim = getattr(ax, f"get_{axis_name}lim")
    set_lim = getattr(ax, f"set_{axis_name}lim")
    lo, hi = get_lim()
    data_lo, data_hi = _clamp_fraction_limits(lo, hi)
    if lo > hi:
        if lo != data_hi or hi != data_lo:
            set_lim(data_hi, data_lo, emit=False)
    elif lo != data_lo or hi != data_hi:
        set_lim(data_lo, data_hi, emit=False)


def clamp_composition_axis_limits(ax: Any, *, is_3d: bool = False) -> None:
    """Clamp all composition axes to [0, 1], preserving reversed display order."""
    _clamp_axis_limits(ax, "x")
    _clamp_axis_limits(ax, "y")
    if is_3d:
        _clamp_axis_limits(ax, "z")


def composition_axis_ticks(lo: float, hi: float, *, nbins: int = 6) -> np.ndarray:
    """Return non-negative tick positions for a composition axis segment."""
    data_lo, data_hi = _clamp_fraction_limits(min(lo, hi), max(lo, hi))
    if data_hi <= data_lo:
        return np.array([0.0])

    ticks = MaxNLocator(nbins=nbins).tick_values(data_lo, data_hi)
    ticks = ticks[(ticks >= 0.0) & (ticks <= 1.0)]

    if data_lo <= 0.0 and not np.any(np.isclose(ticks, 0.0)):
        ticks = np.sort(np.unique(np.concatenate(([0.0], ticks))))

    return ticks


class _CompositionTickLocator(Locator):
    """Locator that only proposes tick positions within the [0, 1] range.

    Unlike calling ``Axes.set_xticks``/``set_yticks``/``set_zticks``
    reactively (whose ``_set_tick_locations`` implicitly *expands* the axis
    view limits to keep every tick visible, via ``Axis.set_view_interval``),
    a ``Locator`` subclass is queried fresh on every draw without ever
    mutating the view limits. This lets the view zoom/pan completely freely
    while the ticks/labels still never show anything outside 0-100%.
    """

    def __init__(self, nbins: int = 6) -> None:
        super().__init__()
        self._nbins = nbins

    def tick_values(self, vmin: float, vmax: float) -> np.ndarray:
        return composition_axis_ticks(vmin, vmax, nbins=self._nbins)

    def __call__(self) -> np.ndarray:
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)


def compute_zoom_limits(
    values: np.ndarray,
    margin_ratio: float = 0.08,
    min_span: float = 0.10,
    central_fraction: float = 0.90,
) -> tuple[float, float]:
    """Return axis limits that frame the central fraction of sample values."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0

    lower_q = max(0.0, (1.0 - central_fraction) / 2.0)
    upper_q = min(1.0, 1.0 - lower_q)
    v_min = float(np.quantile(values, lower_q))
    v_max = float(np.quantile(values, upper_q))
    span = max(v_max - v_min, min_span)
    margin = span * margin_ratio
    low = max(0.0, v_min - margin)
    high = min(1.0, v_max + margin)

    if high <= low:
        center = 0.5 * (v_min + v_max)
        half = max(min_span * 0.5, 0.02)
        low = max(0.0, center - half)
        high = min(1.0, center + half)

    return low, high


def expand_limits(
    low: float,
    high: float,
    margin_ratio: float = 0.20,
) -> tuple[float, float]:
    """Expand axis limits by a fraction of their span, clamped to [0, 1]."""
    span = high - low
    margin = span * margin_ratio
    return max(0.0, low - margin), min(1.0, high + margin)


def int_percent_formatter() -> FuncFormatter:
    """Format fractional composition values as integer percentages."""
    return FuncFormatter(lambda value, _: f"{int(round(value * 100.0))}")


def apply_dynamic_axis_formatting(ax: Any, *, is_3d: bool = False, nbins: int = 6) -> None:
    """Install self-updating tick locators/formatters for the composition axes.

    The locator recomputes tick positions from the *current* view limits on
    every draw (e.g. after an interactive pan/zoom/rotate) and always clamps
    them to [0, 1], so labels never appear outside the 0-100% range, while
    the view itself remains free to zoom/pan without limit.
    """
    formatter = int_percent_formatter()
    ax.xaxis.set_major_locator(_CompositionTickLocator(nbins=nbins))
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_locator(_CompositionTickLocator(nbins=nbins))
    ax.yaxis.set_major_formatter(formatter)
    if is_3d:
        ax.zaxis.set_major_locator(_CompositionTickLocator(nbins=nbins))
        ax.zaxis.set_major_formatter(formatter)


def apply_fixed_full_range_ticks(ax: Any, *, is_3d: bool = False) -> None:
    """Use fixed 10%-step ticks for the full 0-100% composition range."""
    ticks = np.arange(0, 1, 0.1)
    tick_labels = [f"{x * 100:.0f}" for x in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels)
    if is_3d:
        ax.set_zticks(ticks)
        ax.set_zticklabels(tick_labels)


def gather_clustering_zoom_points(
    els_comps_list: np.ndarray,
    centroids: np.ndarray,
    unused_compositions_list: Sequence,
    elements: List[str],
    ref_phases_df: Optional[Any] = None,
) -> np.ndarray:
    """Collect points used to compute data-driven clustering plot axis limits."""
    n_dims = len(elements)
    zoom_points = [np.asarray(els_comps_list, dtype=float).T]
    if unused_compositions_list:
        zoom_points.append(np.asarray(unused_compositions_list, dtype=float))
    zoom_points.append(np.asarray(centroids, dtype=float))

    if ref_phases_df is not None and unused_compositions_list:
        ref_phases_df_zoom = ref_phases_df[elements]
        unused_arr = np.asarray(unused_compositions_list, dtype=float)
        threshold = 0.20
        for _, row in ref_phases_df_zoom.iterrows():
            ref_point = np.array(row.values, dtype=float)
            dists = np.linalg.norm(unused_arr - ref_point, axis=1)
            if np.any(dists < threshold):
                zoom_points.append(ref_point.reshape(1, -1))

    if not zoom_points:
        return np.empty((0, n_dims))
    return np.vstack([pts for pts in zoom_points if pts.size > 0])


def point_within_composition_limits(
    point: Sequence[float],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    zlim: Optional[tuple[float, float]] = None,
) -> bool:
    """Return True when a composition-space point lies inside the given axis limits."""
    x_lo, x_hi = _clamp_fraction_limits(min(xlim), max(xlim))
    y_lo, y_hi = _clamp_fraction_limits(min(ylim), max(ylim))
    coords = np.asarray(point, dtype=float)
    if not (x_lo <= coords[0] <= x_hi and y_lo <= coords[1] <= y_hi):
        return False
    if zlim is None:
        return True
    z_lo, z_hi = _clamp_fraction_limits(min(zlim), max(zlim))
    return z_lo <= coords[2] <= z_hi


def compute_data_driven_axis_limits(
    all_points: np.ndarray,
    *,
    is_3d: bool = False,
    margin_ratio: float = 0.20,
) -> tuple[tuple[float, float], tuple[float, float], Optional[tuple[float, float]]]:
    """Return data-driven composition axis limits for clustering zoom plots."""
    default = (0.0, 1.0)
    if all_points.size == 0:
        return default, default, default if is_3d else None

    x_low, x_high = compute_zoom_limits(all_points[:, 0])
    y_low, y_high = compute_zoom_limits(all_points[:, 1])
    x_low, x_high = expand_limits(x_low, x_high, margin_ratio=margin_ratio)
    y_low, y_high = expand_limits(y_low, y_high, margin_ratio=margin_ratio)
    x_low, x_high = _clamp_fraction_limits(x_low, x_high)
    y_low, y_high = _clamp_fraction_limits(y_low, y_high)

    zlim = None
    if is_3d:
        z_low, z_high = compute_zoom_limits(all_points[:, 2])
        z_low, z_high = expand_limits(z_low, z_high, margin_ratio=margin_ratio)
        zlim = _clamp_fraction_limits(z_low, z_high)

    return (x_low, x_high), (y_low, y_high), zlim


def apply_data_driven_axis_limits(
    ax: Any,
    all_points: np.ndarray,
    *,
    is_3d: bool = False,
    reversed_xy: bool = False,
    margin_ratio: float = 0.20,
) -> None:
    """Set axis limits from clustered data and apply adaptive tick formatting."""
    xlim, ylim, zlim = compute_data_driven_axis_limits(
        all_points,
        is_3d=is_3d,
        margin_ratio=margin_ratio,
    )
    x_low, x_high = xlim
    y_low, y_high = ylim

    if is_3d and reversed_xy:
        ax.set_xlim(x_high, x_low)
        ax.set_ylim(y_high, y_low)
    else:
        ax.set_xlim(x_low, x_high)
        ax.set_ylim(y_low, y_high)

    if is_3d and zlim is not None:
        ax.set_zlim(*zlim)

    clamp_composition_axis_limits(ax, is_3d=is_3d)
    apply_dynamic_axis_formatting(ax, is_3d=is_3d)


def configure_interactive_clustering_axes(
    ax: Any,
    *,
    is_3d: bool = False,
    all_points: Optional[np.ndarray] = None,
    reversed_xy: bool = False,
) -> None:
    """Configure axes for interactive clustering plots with zoom-aware ticks."""
    if all_points is not None and all_points.size > 0:
        apply_data_driven_axis_limits(
            ax,
            all_points,
            is_3d=is_3d,
            reversed_xy=reversed_xy,
        )
    else:
        clamp_composition_axis_limits(ax, is_3d=is_3d)
        apply_dynamic_axis_formatting(ax, is_3d=is_3d)
