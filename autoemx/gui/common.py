#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared state and helpers of the AutoEMX GUI tabs: data caches, job manager, parameter forms."""

from __future__ import annotations

import os
import subprocess
import sys
import threading
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from dash import dcc, html

from autoemx._paths import AUTOEMX_HOME
from autoemx.gui import backend as be

GRAPH_CONFIG = {
    "displaylogo": False,
    "toImageButtonOptions": {"format": "png", "scale": 3},
    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
}

JOBS = be.JobManager()

# =============================================================================
# Cached loading (keyed on file modification times)
# =============================================================================
_info_cache: Dict[str, Tuple[float, be.SampleInfo]] = {}
_analysis_cache: Dict[Tuple[Any, ...], be.AnalysisData] = {}
# Callbacks run in parallel threads; loading ledgers and calibration modules is not thread-safe.
_load_lock = threading.RLock()


def _mtime(path: str) -> float:
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


def get_info(sample_dir: Optional[str]) -> Optional[be.SampleInfo]:
    if not sample_dir:
        return None
    with _load_lock:
        return _get_info(sample_dir)


def _get_info(sample_dir: str) -> be.SampleInfo:
    m = _mtime(os.path.join(sample_dir, be.LEDGER_NAME))
    cached = _info_cache.get(sample_dir)
    if cached and cached[0] == m:
        return cached[1]
    info = be.load_sample_info(sample_dir)
    _info_cache[sample_dir] = (m, info)
    return info


_sizes_cache: Dict[Tuple[str, Optional[int]], Tuple[float, Optional[be.ParticleSizeData]]] = {}


def get_particle_sizes(sample_dir: Optional[str], result_id: Optional[int] = None) -> Optional[be.ParticleSizeData]:
    """Particle size statistics of a sample: result ``result_id`` (default: the active one), cached (computing them
    may run a bootstrap)."""
    if not sample_dir:
        return None
    with _load_lock:
        m = _mtime(os.path.join(sample_dir, be.LEDGER_NAME))
        cached = _sizes_cache.get((sample_dir, result_id))
        if cached and cached[0] == m:
            return cached[1]
        data = be.load_particle_sizes(_get_info(sample_dir), result_id)
        if len(_sizes_cache) > 200:
            _sizes_cache.clear()
        _sizes_cache[(sample_dir, result_id)] = (m, data)
        return data


_has_sizes_cache: Dict[str, Tuple[float, bool]] = {}


def has_sizes(sample_dir: str) -> bool:
    """Whether a sample has particle size data (cached, from its ledger file)."""
    m = _mtime(os.path.join(sample_dir, be.LEDGER_NAME))
    cached = _has_sizes_cache.get(sample_dir)
    if cached and cached[0] == m:
        return cached[1]
    value = be.sample_has_sizes(sample_dir)
    _has_sizes_cache[sample_dir] = (m, value)
    return value


def get_analysis(sample_dir: Optional[str], key: Optional[str]) -> Tuple[Optional[be.SampleInfo], Optional[be.AnalysisData]]:
    if not sample_dir:
        return None, None
    with _load_lock:
        return _get_analysis(sample_dir, key)


def _get_analysis(sample_dir: str, key: Optional[str]) -> Tuple[be.SampleInfo, be.AnalysisData]:
    info = _get_info(sample_dir)
    if not info.has_eds:  # e.g. samples with only particle size data
        return info, None
    ref = be.find_analysis(info, key)
    folder = ref.folder if ref else ""
    cache_key = (sample_dir, ref.key if ref else None, _mtime(os.path.join(sample_dir, be.LEDGER_NAME)),
                 _mtime(os.path.join(folder, "Compositions.csv")) if folder else 0)
    data = _analysis_cache.get(cache_key)
    if data is None:
        data = be.load_analysis(info, ref)
        if len(_analysis_cache) > 20:
            _analysis_cache.clear()
        _analysis_cache[cache_key] = data
    return info, data


# =============================================================================
# Parameter form
# =============================================================================
# Component-id type of the parameters of each form (Dash ids cannot contain dots in values)
FORM_ID_TYPES = {
    **{key: "param" for key, _ in be.SECTIONS},
    **{key: "qparam" for key, _ in be.QUANT_SECTIONS},
    **{key: "sparam" for key, _ in be.SINGLE_SECTIONS},
    **{key: "aparam" for key, _ in be.ACQ_SECTIONS},
}


_PAIR_PARTS = ("__min", "__max")  # suffixes of the two fields of a (min, max) parameter


def _to_ui(spec: be.ParamSpec, value: Any) -> Any:
    if spec.kind == "bool":
        return ["on"] if value else []
    if spec.kind == "formulae":
        return "\n".join(value or [])
    if spec.kind in ("elements", "flags", "multi"):
        return list(value or [])
    return value


def _from_ui(spec: be.ParamSpec, value: Any) -> Any:
    if spec.kind == "bool":
        return bool(value) and "on" in value
    return value


def _pid(key: str, part: str = "") -> Dict[str, str]:
    """Component id of a parameter (``part``: '__min' / '__max' for the fields of (min, max) parameters)."""
    section = key.split(".", 1)[0]
    return {"type": FORM_ID_TYPES[section], "key": key.replace(".", "__") + part}


def _pair_part(component_id: Dict[str, str]) -> Optional[int]:
    """Index (0: min, 1: max) of the field of a (min, max) parameter, or None."""
    for i, suffix in enumerate(_PAIR_PARTS):
        if component_id["key"].endswith(suffix):
            return i
    return None


def _spec(component_id: Dict[str, str]) -> be.ParamSpec:
    key = component_id["key"]
    if _pair_part(component_id) is not None:
        key = key[: -len("__min")]
    return be.SPECS_BY_KEY[key.replace("__", ".")]


def form_values(states: List[Dict[str, Any]]) -> Dict[str, Any]:
    """GUI values of a form, from the ``ctx.states_list`` entry of its ALL-pattern State."""
    values: Dict[str, Any] = {}
    for s in states:
        spec = _spec(s["id"])
        part = _pair_part(s["id"])
        if part is None:
            values[spec.key] = _from_ui(spec, s.get("value"))
        else:
            pair = values.setdefault(spec.key, [None, None])
            pair[part] = s.get("value")
    return values


def ui_value(component_id: Dict[str, str], values: Dict[str, Any]) -> Any:
    """Value to show in a form component, from parameter values (as from ``backend.*_param_values``)."""
    spec = _spec(component_id)
    value = values.get(spec.key)
    part = _pair_part(component_id)
    if part is not None:
        return value[part] if value else None
    return _to_ui(spec, value)


def _row_class(spec: be.ParamSpec) -> str:
    if spec.kind == "bool":
        return "param-row param-bool"
    return "param-row" + (" param-stacked" if spec.kind in ("flags", "multi", "formulae", "elements") else "")


def _row_id(key: str) -> Dict[str, str]:
    return {"type": _pid(key)["type"] + "-row", "key": key.replace(".", "__")}


def _clear_button(spec: be.ParamSpec) -> html.Button:
    """↺ button emptying an optional parameter, i.e. back to its saved/automatic value."""
    pid = _pid(spec.key)
    return html.Button("↺", id={"type": pid["type"] + "-clear", "key": pid["key"]}, n_clicks=0,
                       className="clear-btn", title=f"Clear: use the {spec.placeholder or 'default'} value")


def _param_control(spec: be.ParamSpec, persist: bool = False) -> html.Div:
    """Control of one parameter. With ``persist``, values set by the user are kept in the browser between sessions."""
    pid = _pid(spec.key)
    value = _to_ui(spec, spec.default)
    keep = {"persistence": "v1", "persistence_type": "local"} if persist else {}
    # Clicking "?" shows the help text below the parameter (assets/help.js); hovering shows it as a tooltip
    help_icon = (html.Span("?", className="help", title=spec.help, role="button", tabIndex="0",
                           **{"data-help": spec.help}) if spec.help else None)
    if spec.kind == "bool":
        ctrl = dcc.Checklist(id=pid, options=[{"label": spec.label, "value": "on"}], value=value,
                             className="param-check", **keep)
        return html.Div([ctrl, help_icon], id=_row_id(spec.key), className=_row_class(spec), title=spec.help)
    placeholder = spec.placeholder or ("none" if spec.kind.endswith("_opt") else "")
    if spec.kind in ("int", "int_opt", "float", "float_opt"):
        ctrl = dcc.Input(id=pid, type="number", value=value, step="any", debounce=True,
                         placeholder=placeholder, className="param-input", **keep)
    elif spec.kind in ("pair", "pair_opt"):
        lo, hi = (value or [None, None])[:2]
        ctrl = html.Div([
            dcc.Input(id=_pid(spec.key, "__min"), type="number", value=lo, step=1, debounce=True,
                      placeholder="min", className="param-input", **keep),
            dcc.Input(id=_pid(spec.key, "__max"), type="number", value=hi, step=1, debounce=True,
                      placeholder="max", className="param-input", **keep),
        ], className="pair-inputs")
    elif spec.kind == "choice":
        ctrl = dcc.Dropdown(id=pid, options=[{"label": str(c), "value": c} for c in spec.choices],
                            value=value, clearable=False, className="param-dd", **keep)
    elif spec.kind == "elements":
        ctrl = dcc.Dropdown(id=pid, options=[], value=value, multi=True, placeholder="default",
                            className="param-dd")
    elif spec.kind == "flags":
        ctrl = dcc.Checklist(
            id=pid, value=value, className="flag-list",
            options=[{"label": f"{f} · {m}", "value": f} for f, m in be.QUANT_FLAG_MEANINGS.items()], **keep,
        )
    elif spec.kind == "multi":
        ctrl = dcc.Checklist(id=pid, value=value, className="flag-list", inline=True,
                             options=[{"label": str(c), "value": c} for c in spec.choices], **keep)
    elif spec.kind == "formulae":
        ctrl = dcc.Textarea(id=pid, value=value, placeholder="MgO\nAl2O3\nMgAl2O4", className="param-text",
                            rows=3)
    else:
        ctrl = dcc.Input(id=pid, type="text", value=value, debounce=True, placeholder=placeholder,
                         className="param-input", **keep)
    if spec.unit or spec.kind in ("int_opt", "float_opt", "pair_opt"):
        parts = [ctrl]
        if spec.unit:
            parts.append(html.Span(spec.unit, className="param-unit"))
        if spec.kind in ("int_opt", "float_opt", "pair_opt"):
            parts.append(_clear_button(spec))
        ctrl = html.Div(parts, className="param-ctrl")
    return html.Div(
        [html.Label([spec.label, help_icon], className="param-label", title=spec.help), ctrl],
        id=_row_id(spec.key), className=_row_class(spec),
    )


def register_clear_buttons(app) -> None:
    """Callbacks of the ↺ buttons of optional parameters, one per form."""
    from dash import ALL, Input, Output, ctx, no_update

    for id_type in sorted(set(FORM_ID_TYPES.values())):
        @app.callback(
            Output({"type": id_type, "key": ALL}, "value", allow_duplicate=True),
            Input({"type": id_type + "-clear", "key": ALL}, "n_clicks"),
            prevent_initial_call=True,
        )
        def clear(_clicks):
            trig = ctx.triggered_id
            if not isinstance(trig, dict) or not ctx.triggered[0].get("value"):
                return [no_update] * len(ctx.outputs_list)
            return [None if o["id"]["key"] in (trig["key"], trig["key"] + "__min", trig["key"] + "__max")
                    else no_update for o in ctx.outputs_list]


def param_sections(
    specs: Sequence[be.ParamSpec],
    sections: Sequence[Tuple[str, str]],
    notes: Optional[Dict[str, str]] = None,
    open_sections: Sequence[str] = (),
    persist: bool = False,
    extra: Optional[Dict[str, List[Any]]] = None,
) -> List[html.Details]:
    """Collapsible sections of a parameter form (``extra``: components appended to a section, by section key)."""
    out = []
    for key, title in sections:
        section_specs = [s for s in specs if s.section == key]
        note = (notes or {}).get(key)
        summary = [html.Span(title)]
        if note:
            summary.append(html.Span(note, className="section-note"))
        out.append(html.Details(
            [html.Summary(summary)] + [_param_control(s, persist) for s in section_specs]
            + list((extra or {}).get(key, [])),
            id={"type": FORM_ID_TYPES[key] + "-section", "key": key}, open=key in open_sections,
            className="param-section",
        ))
    return out


def _chip(label: str, value: Any, cls: str = "") -> html.Span:
    return html.Span([html.Span(label, className="chip-label"), html.Span(str(value), className="chip-value")],
                     className="chip " + cls)


def peak_overlaps_view(overlaps: Sequence[Dict[str, Any]], title: str = "") -> Optional[html.Div]:
    """Warning listing the peak overlaps that may compromise the quantification (None if there are none)."""
    if not overlaps:
        return None
    from autoemx.core.quantifier.peak_overlaps import affected_elements, format_overlap

    items = [html.Li(format_overlap(o)) for o in overlaps]
    header = f"{title}: p" if title else "P"
    return html.Div(
        [html.Div(f"⚠️ {header}eak overlaps may cause inaccurate quantification of "
                  f"{', '.join(affected_elements(overlaps))}:"),
         html.Ul(items, className="overlap-list")],
        className="warn overlap-warn",
    )


def _fmt(v: Any, nd: int = 1) -> str:
    try:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "—"
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return str(v)


def _pick_file(prompt: str = "Select a spectrum file") -> Optional[str]:
    """Native file dialog (in a subprocess, so it never blocks the server thread)."""
    try:
        if sys.platform == "darwin":
            res = subprocess.run(["osascript", "-e", f'POSIX path of (choose file with prompt "{prompt}")'],
                                 capture_output=True, text=True, timeout=600)
        else:
            code = ("import tkinter as tk; from tkinter import filedialog; r = tk.Tk(); r.withdraw(); "
                    "r.attributes('-topmost', True); print(filedialog.askopenfilename(title='" + prompt + "', "
                    "filetypes=[('EMSA spectra', '*.msa *.emsa *.msg'), ('All files', '*')]))")
            res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600)
    except Exception:
        return None
    return res.stdout.strip() or None


def _pick_folder(prompt: str = "Select the AutoEMX results folder") -> Optional[str]:
    """Native folder dialog (in a subprocess, so it never blocks the server thread)."""
    try:
        if sys.platform == "darwin":
            res = subprocess.run(
                ["osascript", "-e", f'POSIX path of (choose folder with prompt "{prompt}")'],
                capture_output=True, text=True, timeout=600,
            )
        else:
            code = ("import tkinter as tk; from tkinter import filedialog; r = tk.Tk(); r.withdraw(); "
                    "r.attributes('-topmost', True); print(filedialog.askdirectory(title='" + prompt + "'))")
            res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=600)
    except Exception:
        return None
    path = res.stdout.strip()
    return path.rstrip("/") if path else None


def _open_path(path: str) -> None:
    if sys.platform == "darwin":
        subprocess.Popen(["open", path])
    elif os.name == "nt":
        os.startfile(path)  # type: ignore[attr-defined]
    else:
        subprocess.Popen(["xdg-open", path])


def _point_spectrum(click: Optional[Dict[str, Any]]) -> Optional[str]:
    if not click or not click.get("points"):
        return None
    cd = click["points"][0].get("customdata")
    if isinstance(cd, list) and cd:
        return str(cd[0])
    if isinstance(cd, (str, int, float)):
        return str(cd)
    return None




# =============================================================================
# Sample summaries (samples table), loaded in the background and cached on disk
# =============================================================================
class SummaryLoader:
    """
    Reads the summary of every sample of a folder in background threads.

    Ledgers on cloud drives can take seconds each to download, so summaries are cached on disk
    (keyed on ledger path, size and modification time) and only changed ledgers are re-read.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._summaries: Dict[str, Dict[str, Any]] = {}
        self._pending: set = set()
        self._cache_path = str(AUTOEMX_HOME / "gui_sample_summaries.json")
        self._disk: Dict[str, Any] = self._read_disk()

    def _read_disk(self) -> Dict[str, Any]:
        import json

        try:
            with open(self._cache_path, encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            return {}

    def _write_disk(self) -> None:
        import json

        try:
            os.makedirs(os.path.dirname(self._cache_path), exist_ok=True)
            tmp = self._cache_path + ".tmp"
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(self._disk, fh)
            os.replace(tmp, self._cache_path)
        except Exception:
            pass

    @staticmethod
    def _stamp(sample_dir: str) -> Optional[List[float]]:
        try:
            st = os.stat(os.path.join(sample_dir, be.LEDGER_NAME))
            return [st.st_mtime, st.st_size]
        except OSError:
            return None

    def get(self, sample_dirs: Sequence[str]) -> Tuple[Dict[str, Dict[str, Any]], int]:
        """Summaries available now (by sample dir), and the number still loading. Starts loading the others."""
        out: Dict[str, Dict[str, Any]] = {}
        to_load = []
        with self._lock:
            for d in sample_dirs:
                stamp = self._stamp(d)
                cached = self._disk.get(d)
                if cached and stamp is not None and cached.get("stamp") == stamp and "meas_type" in cached["summary"]:
                    out[d] = cached["summary"]
                elif d in self._pending:
                    if d in self._summaries:
                        out[d] = self._summaries[d]  # previous version while reloading
                else:
                    self._pending.add(d)
                    to_load.append(d)
            n_pending = len([d for d in sample_dirs if d in self._pending])
        if to_load:
            threading.Thread(target=self._load, args=(to_load,), daemon=True).start()
        return out, n_pending

    def _load(self, sample_dirs: List[str]) -> None:
        from concurrent.futures import ThreadPoolExecutor

        def one(d: str) -> None:
            try:
                stamp = self._stamp(d)
                summary = be.sample_summary(d)
                with self._lock:
                    self._summaries[d] = summary
                    self._disk[d] = {"stamp": stamp, "summary": summary}
            except Exception:
                pass
            finally:
                with self._lock:
                    self._pending.discard(d)

        with ThreadPoolExecutor(max_workers=6) as pool:
            list(pool.map(one, sample_dirs))
        with self._lock:
            self._write_disk()


SUMMARIES = SummaryLoader()
