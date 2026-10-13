#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Backend of the AutoEMX GUI (no Dash code).

- Finds AutoEMX sample folders (folders containing ``ledger.json``).
- Describes every quantification/analysis parameter, with its value saved in the ledger.
- Runs ``batch_quantify_and_analyze`` / ``analyze_sample`` / single-spectrum fits in a
  separate process, capturing their log.
- Loads the results of any clustering analysis saved in a sample, for interactive plotting.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import re
import signal
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import types
import typing
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import autoemx.utils.constants as cnst
import autoemx._paths as paths
from autoemx.utils.helper import parse_xsp_spots_image_name
from autoemx.config.ledger_io import get_active_quantification_config, load_sample_ledger
from autoemx.config.ledger_schemas import (
    AitchisonParams,
    ClusterMergeParams,
    ClusteringConfig,
    DBSCANParams,
    MixtureParams,
)
import autoemx.config.defaults as dflt
from autoemx.config.runtime_configs import PlotConfig, QuantificationOptionsConfig

LEDGER_NAME = f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
_DEFAULT_UNDETECTABLE_ELS = ("H", "He", "Li")
_SKIPPED_DIRS = {"spectra", "__pycache__", ".git"}

QUANT_FLAG_MEANINGS: Dict[int, str] = {
    0: "OK",
    -1: "OK, not converged",
    1: "Acquisition error / no data",
    2: "Total counts too low",
    3: "Low-energy background too low",
    4: "Poor fit",
    5: "Analytical error > 50 w%",
    6: "Excessive X-ray absorption",
    7: "Substrate contamination",
    8: "Background under reference peak too low",
    9: "Fit interrupted (unknown)",
}


# =============================================================================
# Sample discovery and ledger summaries
# =============================================================================
@dataclass
class SampleRef:
    sample_id: str
    sample_dir: str
    label: str


def find_samples(root: str, max_depth: int = 5) -> List[SampleRef]:
    """Return every folder under ``root`` (included) that contains a ``ledger.json``."""
    root_path = Path(root).expanduser()
    if not root_path.is_dir():
        raise FileNotFoundError(f"Folder not found: {root}")
    samples: List[SampleRef] = []
    root_depth = len(root_path.parts)
    for dirpath, dirnames, filenames in os.walk(root_path):
        depth = len(Path(dirpath).parts) - root_depth
        dirnames[:] = sorted(
            d for d in dirnames
            if d not in _SKIPPED_DIRS and not d.startswith(("analysis_quant", "."))
        )
        if depth >= max_depth:
            dirnames[:] = []
        if LEDGER_NAME in filenames:
            path = Path(dirpath)
            rel = path.relative_to(root_path)
            label = path.name if str(rel) == "." else str(rel)
            samples.append(SampleRef(sample_id=path.name, sample_dir=str(path), label=label))
    return sorted(samples, key=lambda s: s.label.lower())


_summary_cache: Dict[str, Tuple[float, Dict[str, Any]]] = {}


def _acquisition_date(sample_dir: Path, ledger: Dict[str, Any]) -> str:
    """Date/time of the first spectrum (EMSA header), else of the ledger file, as 'YYYY-MM-DD HH:MM'."""
    from datetime import datetime

    for sp in ledger.get("spectra", [])[:1]:
        rel = sp.get("spectrum_relpath")
        if not rel:
            continue
        try:
            header: Dict[str, str] = {}
            with open(sample_dir / rel, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.upper().startswith("#SPECTRUM"):
                        break
                    if ":" in line:
                        k, v = line[1:].split(":", 1)
                        header[k.strip().upper()] = v.strip()
            date = datetime.strptime(header["DATE"].title(), "%d-%b-%Y")
            hh, mm = (header.get("TIME", "00:00").split(":") + ["00"])[:2]
            return date.replace(hour=int(hh), minute=int(mm[:2])).strftime("%Y-%m-%d %H:%M")
        except Exception:
            pass
    stat = (sample_dir / LEDGER_NAME).stat()
    ts = getattr(stat, "st_birthtime", stat.st_mtime)
    return datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M")


def _saved_quant_options(cfgs: Dict[str, Any], active: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Values of the quantification-form options saved for a sample (keys of QUANT_SPECS)."""
    defaults = QuantificationOptionsConfig()
    opts = dict((active or {}).get("options") or {})
    analyses = (active or {}).get("clustering_analyses") or []
    idx = (active or {}).get("active_clustering_analysis_index")
    cluster_cfg = (analyses[idx] if idx is not None and idx < len(analyses) else (analyses[-1] if analyses else {}))
    min_cnts = (cluster_cfg.get("config") or {}).get("min_bckgrnd_cnts", ClusteringConfig().min_bckgrnd_cnts)
    lims = opts.get("spectrum_lims", defaults.spectrum_lims)
    return {
        "quant.min_bckgrnd_cnts": min_cnts,
        "quant.min_total_counts_fraction": opts.get("min_total_counts_fraction", defaults.min_total_counts_fraction),
        "quant.spectrum_lims": [int(round(v)) for v in lims] if lims else None,
        "quant.fit_tolerance": opts.get("fit_tolerance", defaults.fit_tolerance),
        "quant.use_project_specific_std_dict": bool(opts.get("use_project_specific_std_dict",
                                                            defaults.use_project_specific_std_dict)),
        "quant.is_known_phase_mixture": bool(opts.get("is_known_phase_mixture", False)),
    }


CALIBRATIONS_DIR = Path(__file__).resolve().parents[1] / "calibrations"


def available_microscopes() -> List[str]:
    """Microscopes with a detector model (folders of autoemx/calibrations with a detector_model.json)."""
    ids = [d.name for d in sorted(CALIBRATIONS_DIR.iterdir()) if (d / cnst.DETECTOR_MODEL_FILENAME).exists()]
    return ids or [dflt.microscope_ID]


def standards_beam_energies(microscope_id: str, meas_type: str = dflt.measurement_type) -> List[float]:
    """Beam energies (keV) for which the microscope has a P/B standards file.

    Those in the local calibration folder, and those of the package (copied to the local folder when first used).
    """
    energies = set()
    for folder in (paths.user_calibration_dir(microscope_id), paths.package_calibration_dir(microscope_id)):
        for path in folder.glob(f"{meas_type}_{cnst.STD_FILENAME}_*keV.json"):
            match = re.search(r"_(\d+(?:\.\d+)?)keV$", path.stem)
            if match:
                energies.add(float(match.group(1)))
    return sorted(energies)


def quantifiable_elements(microscope_id: str, beam_energy_keV: Any, meas_mode: str = dflt.measurement_mode,
                          meas_type: str = dflt.measurement_type) -> Optional[List[str]]:
    """Elements with at least one P/B standard for this microscope, beam energy and mode (None: no standards file)."""
    try:
        kv = float(beam_energy_keV)
    except (TypeError, ValueError):
        return None
    if abs(kv - round(kv)) > 0.25:
        return None
    path = paths.user_standards_file(microscope_id, meas_type, round(kv))
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            payload = json.load(fh)
        modes = payload.get("standards_by_mode", payload)
        lines = modes.get(meas_mode) or next(iter(modes.values()), {})
        return sorted({key.split("_")[0] for key in lines})
    except Exception:
        return None


def elements_without_standards(elements: List[str], microscope_id: str, beam_energy_keV: Any,
                               meas_mode: str = dflt.measurement_mode) -> Optional[List[str]]:
    """Detectable elements that cannot be quantified (None: no standards file at all)."""
    available = quantifiable_elements(microscope_id, beam_energy_keV, meas_mode)
    if available is None:
        return None
    return [el for el in elements if el not in available and el not in _DEFAULT_UNDETECTABLE_ELS]


def common_saved_values(summaries: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Saved quantification options shared by all ``summaries`` (key: value); options that differ are omitted."""
    out: Dict[str, Any] = {}
    opts = [s.get("quant_opts") for s in summaries]
    if not opts or any(o is None for o in opts):
        return out
    for key in opts[0]:
        values = [o.get(key) for o in opts]
        if all(v == values[0] for v in values) and values[0] is not None:
            out[key] = values[0]
    return out


def sample_summary(sample_dir: str) -> Dict[str, Any]:
    """Quick summary of a sample read from its ledger JSON (no validation), for the samples table."""
    path = Path(sample_dir)
    ledger_path = path / LEDGER_NAME
    mtime = ledger_path.stat().st_mtime
    cached = _summary_cache.get(sample_dir)
    if cached and cached[0] == mtime:
        return cached[1]
    with open(ledger_path, encoding="utf-8") as fh:
        ledger = json.load(fh)
    quants = ledger.get("quantifications") or []
    active_id = ledger.get("active_quant")
    active = next((q for q in quants if q.get("quantification_id") == active_id), quants[-1] if quants else None)
    cfgs = ledger.get("configs") or {}
    elements = (active or {}).get("sample_elements") or (cfgs.get("sample_cfg") or {}).get("elements") or []
    substrate = (active or {}).get("substrate_elements")
    if substrate is None:
        substrate = (cfgs.get("sample_substrate_cfg") or {}).get("elements") or []
    spectra = ledger.get("spectra") or []
    n_quant = 0
    if active is not None:
        qid = active.get("quantification_id")
        n_quant = sum(
            any(r.get("quantification_id") == qid and r.get("composition_atomic_fractions")
                for r in (sp.get("quantification_results") or []))
            for sp in spectra
        )
    summary = {
        "sample": path.name,
        "dir": str(path),
        "date": _acquisition_date(path, ledger),
        "n_spectra": len(spectra),
        "n_quantified": n_quant,
        "n_runs": len(quants),
        "elements": ", ".join(elements),
        "substrate": ", ".join(substrate),
        "quant_opts": _saved_quant_options(cfgs, active),
        "microscope": (cfgs.get("microscope_cfg") or {}).get("ID", dflt.microscope_ID),
        "beam_energy_keV": (cfgs.get("measurement_cfg") or {}).get("beam_energy_keV"),
        "meas_mode": (cfgs.get("measurement_cfg") or {}).get("mode", dflt.measurement_mode),
        "meas_type": (cfgs.get("measurement_cfg") or {}).get("type", dflt.measurement_type),
        "energy_zero": (cfgs.get("microscope_cfg") or {}).get("energy_zero"),
        "bin_width": (cfgs.get("microscope_cfg") or {}).get("bin_width"),
    }
    _summary_cache[sample_dir] = (mtime, summary)
    return summary


def quantification_runs(info: "SampleInfo") -> List[Dict[str, Any]]:
    """One row per quantification run of a sample: settings, spectra quantified and quant flags."""
    rows = []
    for q in info.ledger.quantifications:
        flags: Dict[int, int] = {}
        n_comp = 0
        n_done = 0
        for sp in info.ledger.spectra:
            rec = next((r for r in reversed(sp.quantification_results) if r.quantification_id == q.quantification_id), None)
            if rec is None:
                continue
            n_done += 1
            if rec.composition_atomic_fractions:
                n_comp += 1
            if rec.quant_flag is not None:
                flags[int(rec.quant_flag)] = flags.get(int(rec.quant_flag), 0) + 1
        opts = q.options or {}
        rows.append({
            "id": q.quantification_id,
            "active": q.quantification_id == info.active_quant,
            "label": q.label or "",
            "elements": ", ".join(q.sample_elements),
            "substrate": ", ".join(q.substrate_elements),
            "n_processed": n_done,
            "n_quantified": n_comp,
            "n_spectra": len(info.ledger.spectra),
            "flags": flags,
            "fit_tolerance": opts.get("fit_tolerance"),
            "spectrum_lims": opts.get("spectrum_lims"),
            "min_total_counts_fraction": opts.get("min_total_counts_fraction"),
            "n_analyses": len(q.clustering_analyses),
        })
    return rows


@dataclass
class AnalysisRef:
    """One clustering analysis saved in the ledger, and the folder holding its outputs."""

    quant_id: int
    clust_index: int
    clust_id: int
    folder: str
    exists: bool
    is_active: bool
    config: ClusteringConfig
    result: Optional[Any]

    @property
    def key(self) -> str:
        return f"{self.quant_id}:{self.clust_index}"

    @property
    def label(self) -> str:
        cfg = self.config
        k = len(self.result.centroids) if self.result is not None else None
        parts = [f"quant{self.quant_id} / clust{self.clust_id}", cfg.method]
        if k is not None:
            parts.append(f"k={k}")
        parts.append(cfg.geometry)
        if cfg.ref_formulae:
            parts.append(", ".join(cfg.ref_formulae[:4]) + ("…" if len(cfg.ref_formulae) > 4 else ""))
        label = " · ".join(parts)
        if self.is_active:
            label += "  (active)"
        if not self.exists:
            label += "  [no output folder]"
        return label


@dataclass
class SampleInfo:
    sample_id: str
    sample_dir: str
    elements: List[str]
    substrate: List[str]
    n_spectra: int
    n_quantified: int
    quant_ids: List[int]
    active_quant: Optional[int]
    analyses: List[AnalysisRef]
    ledger: Any = field(repr=False)

    @property
    def active_analysis(self) -> Optional[AnalysisRef]:
        active = [a for a in self.analyses if a.is_active]
        return active[-1] if active else (self.analyses[-1] if self.analyses else None)


def _active_quant_config(ledger) -> Optional[Any]:
    return get_active_quantification_config(ledger)


def undetectable_elements(ledger) -> List[str]:
    """Elements the sample's microscope cannot detect, from its calibration module."""
    import autoemx.calibrations as calibs

    try:
        calibs.load_microscope_calibrations(
            ledger.configs.microscope_cfg.ID, ledger.configs.measurement_cfg.mode
        )
        return list(getattr(calibs, "undetectable_els", _DEFAULT_UNDETECTABLE_ELS))
    except Exception:
        return list(_DEFAULT_UNDETECTABLE_ELS)


def load_sample_info(sample_dir: str) -> SampleInfo:
    """Load the ledger of a sample and summarise it."""
    ledger = load_sample_ledger(os.path.join(sample_dir, LEDGER_NAME))
    active_q = _active_quant_config(ledger)
    analyses: List[AnalysisRef] = []
    for q in ledger.quantifications:
        is_active_quant = active_q is not None and q.quantification_id == active_q.quantification_id
        for i, ca in enumerate(q.clustering_analyses):
            folder = os.path.join(
                sample_dir, f"analysis_quant{q.quantification_id}_clust{ca.config.clustering_id}"
            )
            analyses.append(
                AnalysisRef(
                    quant_id=q.quantification_id,
                    clust_index=i,
                    clust_id=ca.config.clustering_id,
                    folder=folder,
                    exists=os.path.isdir(folder),
                    is_active=is_active_quant and q.active_clustering_analysis_index == i,
                    config=ca.config,
                    result=ca.result,
                )
            )
    elements = list(active_q.sample_elements) if active_q and active_q.sample_elements else list(ledger.configs.sample_cfg.elements)
    substrate = list(active_q.substrate_elements) if active_q and active_q.substrate_elements else list(ledger.configs.sample_substrate_cfg.elements)
    n_quantified = 0
    if active_q is not None:
        for sp in ledger.spectra:
            if any(
                r.quantification_id == active_q.quantification_id and r.composition_atomic_fractions
                for r in sp.quantification_results
            ):
                n_quantified += 1
    return SampleInfo(
        sample_id=Path(sample_dir).name,
        sample_dir=str(sample_dir),
        elements=elements,
        substrate=substrate,
        n_spectra=len(ledger.spectra),
        n_quantified=n_quantified,
        quant_ids=[q.quantification_id for q in ledger.quantifications],
        active_quant=ledger.active_quant,
        analyses=analyses,
        ledger=ledger,
    )


# =============================================================================
# Parameters
# =============================================================================
@dataclass
class ParamSpec:
    """One GUI control. ``key`` is '<section>.<name>'."""

    section: str
    name: str
    label: str
    kind: str  # bool, int, float, int_opt, float_opt, choice, text, formulae, elements, flags, pair, pair_opt
    default: Any = None
    choices: Optional[List[Any]] = None
    help: str = ""
    placeholder: str = ""
    unit: str = ""  # shown to the right of the field

    @property
    def key(self) -> str:
        return f"{self.section}.{self.name}"


# Analysis form
SECTIONS: List[Tuple[str, str]] = [
    ("filter", "Spectra filtering"),
    ("clust", "Clustering"),
    ("dbscan", "DBSCAN"),
    ("aitchison", "Aitchison geometry"),
    ("merge", "Cluster merging"),
    ("mixture", "Mixture decomposition"),
    ("plot", "Plot output"),
]

# Quantification form (batch_quantify_and_analyze)
QUANT_SECTIONS: List[Tuple[str, str]] = [
    ("quant", "Quantification settings"),
    ("qrun", "Which spectra"),
    ("qafter", "After quantification"),
]

# Single-spectrum form (fit_and_quantify_spectrum[_from_ledger])
SINGLE_SECTIONS: List[Tuple[str, str]] = [
    ("single", "Elements"),
    ("sfit", "Fit and quantification"),
]

# Acquisition form (batch_acquire_and_analyze)
ACQ_SECTIONS: List[Tuple[str, str]] = [
    ("amicro", "Microscope"),
    ("asample", "Sample & Holder"),
    ("aacq", "Acquisition"),
    ("aimg", "Images"),
    ("apowder", "Powder acquisition params"),
    ("abulk", "Bulk grid acquisition params"),
    ("aquant", "Quantification"),
]

SECTION_TITLES: Dict[str, str] = dict(SECTIONS + QUANT_SECTIONS + SINGLE_SECTIONS + ACQ_SECTIONS)

# Choices of options that are kept per sample unless overridden
SAVED, YES, NO = "saved", "yes", "no"
_TRISTATE = {SAVED: None, YES: True, NO: False}

# Sub-models whose fields are all exposed. New fields added to these models appear in the GUI automatically.
_SUBMODEL_SECTIONS: Dict[str, type] = {
    "dbscan": DBSCANParams,
    "aitchison": AitchisonParams,
    "merge": ClusterMergeParams,
    "mixture": MixtureParams,
}
# Sub-models with only some fields exposed (the others keep their default values)
_SUBMODEL_FIELDS: Dict[str, List[str]] = {
    "mixture": ["max_n_phases", "collapse_equivalent_mixtures", "equivalent_span_tol"],
}

_HELP: Dict[str, str] = {
    # Quantification
    "quant.min_bckgrnd_cnts": "Minimum background counts under the reference peaks; spectra below are flagged (quant_flag 8). Empty = keep each sample's saved value.",
    "quant.min_total_counts_fraction": "Minimum total counts, as a fraction of the target acquisition counts; spectra below are flagged (quant_flag 2). 0 disables the check. Empty = keep each sample's saved value.",
    "quant.spectrum_lims": "Lower and upper channel indices of the fitted spectral range, e.g. 14, 1100. Empty = keep each sample's saved value.",
    "quant.fit_tolerance": "Tolerance for fit convergence. Empty = keep each sample's saved value (default 1e-4).",
    "quant.use_project_specific_std_dict": "Load the P/B standards file from the results folder instead of the default calibration file. 'saved' keeps each sample's setting.",
    "quant.is_known_phase_mixture": "Sample is a mixture of known phases (its candidate phases): quantified with the P/B values of their standards. 'saved' keeps each sample's setting.",
    "quant.identify_missing_elements": "Before fitting each spectrum, look for peaks of elements missing from the known elements. Elements surely present are added (quantified if a standard exists, otherwise only fitted); doubtful ones are only reported. The elements found are listed per spectrum and per sample. Slower. 'saved' keeps each sample's setting.",
    "qrun.interrupt_fits_bad_spectra": "Stop fitting spectra as soon as they are found to give a poor quantification. Much faster.",
    "qrun.force_requantification": "Quantify all spectra again, in a new quantification run, even if a run with the same settings exists.",
    "qrun.requantify_only_unquantified_spectra": "Re-quantify only spectra without a composition (never quantified, or previously skipped/flagged).",
    "qrun.max_spectra_to_quantify": "Quantify at most this number of spectra per sample, e.g. for a quick test. Empty = all.",
    "qrun.num_CPU_cores": "CPU cores used for fitting. Empty = half of the available cores.",
    "qafter.run_analysis": "Run the clustering analysis with each sample's saved settings after quantifying it.",
    "qafter.max_analytical_error": "Maximum analytical error (w%) of the spectra used by the analysis.",
    # Single spectrum
    "single.els_sample": "Elements known to be in the sample, quantified. Required for external spectra, unless 'identify missing elements' is on.",
    "single.els_substrate": "Elements fitted but not quantified (e.g. C, O, Al for carbon tape). They are never identified as sample elements either: to quantify an element also present in the substrate (e.g. O in an oxide), list it among the known sample elements instead.",
    "single.is_standard": "The spectrum is from a standard of known composition: report the measured P/B ratios.",
    "single.std_formula": "Formula of the standard (e.g. Al2O3). Empty = composition saved in the ledger.",
    "sfit.quantify": "Quantify the spectrum after fitting it. If off, only fits it.",
    "sfit.identify_elements": "Before the fit, look for peaks of elements missing from the known elements (which may be empty). Elements surely present are added (quantified if a standard exists, otherwise only fitted); doubtful ones are only reported. Slower.",
    "sfit.is_particle": "Apply particle (rough sample) geometry corrections. Off for flat samples.",
    "sfit.fit_tol": "Tolerance for fit convergence.",
    "sfit.spectrum_lims": "Lower and upper channel indices of the fitted range. Empty = saved value (or default for external spectra).",
    "sfit.max_undetectable_w_fr": "Maximum mass fraction of undetectable elements (e.g. Li); the total of fitted elements is constrained to [1 - value, 1].",
    "sfit.force_single_iteration": "Run a single quantification iteration.",
    "sfit.interrupt_fits_bad_spectra": "Stop the fit as soon as the spectrum is found to give a poor quantification.",
    "sfit.free_area_el_lines": "Lines whose area is fitted freely, comma-separated (e.g. Fe_La, Cu_La).",
    # Filtering
    "filter.max_analytical_error_percent": "Spectra whose analytical error exceeds this value (w%) are excluded from clustering. Empty disables the check.",
    "filter.quant_flags_accepted": "Quantification flags kept for clustering.",
    "filter.min_bckgrnd_cnts": "Spectra with fewer background counts under the reference peaks are flagged (quant_flag 8). Applied without re-quantifying; empty keeps the saved value.",
    # Clustering
    "clust.method": "Clustering algorithm.",
    "clust.geometry": "'aitchison' clusters log-ratios (CLR) of compositions; 'auto' chooses per sample.",
    "clust.features": "Cluster atomic (at_fr) or mass (w_fr) fractions.",
    "clust.k_forced": "Number of clusters (k-means). Empty = determined automatically.",
    "clust.k_finding_method": "Method used to find the number of clusters automatically.",
    "clust.max_k": "Maximum number of clusters tested when k is found automatically.",
    "clust.auto_merge_clusters": "Merge k-means clusters that are pieces of one continuous population (e.g. a mixture line).",
    "clust.ref_formulae": "Candidate phases, one formula per line or comma-separated (e.g. MgO, Al2O3, MgAl2O4).",
    "clust.do_matrix_decomposition": "Decompose clusters into mixtures of candidate phases. Slow with many candidates.",
    # DBSCAN
    "dbscan.eps": "Neighbourhood radius. Empty = 0.05 (Euclidean) or 0.3 (Aitchison).",
    "dbscan.min_samples": "Minimum points in a neighbourhood to form a dense region.",
    "dbscan.metric": "Distance metric.",
    # Aitchison
    "aitchison.detection_limit_percent": "Zeros and sub-detection values of trace elements are replaced with 0.65 × this limit before the log-ratio transform.",
    "aitchison.auto_near_zero_percent": "'auto' geometry: threshold below which a major element counts as near zero.",
    "aitchison.auto_max_near_zero_fraction": "'auto' geometry: use Euclidean if at least this fraction of spectra has a near-zero major element.",
    # Merge
    "merge.dip_alpha": "Dip-test p-value above which the projection of two clusters is unimodal.",
    "merge.max_gap_ratio": "Maximum empty gap between two clusters, relative to the spread of the larger one.",
    "merge.dbscan_min_samples": "min_samples of the DBSCAN connectivity check.",
    "merge.dbscan_eps_factor": "eps of the connectivity check, as a multiple of the median k-NN distance.",
    # Mixture
    "mixture.max_n_phases": "Maximum number of candidate phases combined in a mixture (2 = binary only).",
    "mixture.max_recon_error": "Reconstruction error below which a mixture explains a cluster.",
    "mixture.max_recon_error_binary": "Binary mixtures below this error are listed (for inspection).",
    "mixture.recon_error_alpha": "Exponent of the reconstruction-error metric.",
    "mixture.conf_sigma": "Width of the error → confidence mapping.",
    "mixture.single_phase_max_rms_dist": "Clusters tighter than this RMS distance that match a candidate are single-phase.",
    "mixture.single_phase_min_ref_conf": "Minimum candidate confidence for a single-phase cluster.",
    "mixture.nmf_min_mixture_conf": "If no mixture reaches this confidence, a free NMF decomposition is also run.",
    "mixture.equivalent_recon_error_tol": "Mixtures whose errors differ less than this are equivalent.",
    "mixture.max_reported_mixtures": "Maximum mixtures listed per cluster in Clusters.csv.",
    "mixture.report_within_conf_ratio": "Mixtures within this ratio of the best confidence are always listed.",
    "mixture.min_reported_conf_ratio": "Mixtures below this ratio of the best confidence are not listed.",
    "mixture.collapse_equivalent_mixtures": "Hide mixtures spanning the same mixing line/plane as a better one.",
    "mixture.equivalent_span_tol": "Tolerance (fractions) for equivalent mixing lines/planes.",
    # Plot
    "plot.els_to_plot": "Elements forced onto the axes of the saved clustering plot (2 or 3).",
    "plot.els_excluded_clust_plot": "Elements excluded from the saved clustering plot.",
    "plot.show_unused_comps_clust": "Show discarded compositions in the saved clustering plot.",
    "plot.show_legend_clustering": "Show the legend in the saved clustering plot.",
    "plot.plot_best_mixture": "Save a plot of the best mixture of each cluster.",
    "plot.intermixing_violin_plots": "Violin plots of the intermixing of the reported mixtures. 'saved' keeps the sample's setting (by default, plotted for known phase mixtures).",
}

_CHOICES: Dict[str, List[Any]] = {
    "clust.method": list(ClusteringConfig.ALLOWED_METHODS),
    "clust.geometry": list(ClusteringConfig.ALLOWED_GEOMETRIES),
    "clust.features": [cnst.AT_FR_CL_FEAT, cnst.W_FR_CL_FEAT],
    "clust.k_finding_method": ["silhouette", "calinski_harabasz", "elbow"],
    "dbscan.metric": ["euclidean", "manhattan", "chebyshev", "cosine"],
}


def _kind_from_annotation(annotation: Any) -> str:
    origin = typing.get_origin(annotation)
    args = [a for a in typing.get_args(annotation) if a is not type(None)]
    is_optional = origin in (typing.Union, types.UnionType) and len(args) == 1
    base = args[0] if is_optional else annotation
    if base is bool:
        return "bool"
    if base is int:
        return "int_opt" if is_optional else "int"
    if base is float:
        return "float_opt" if is_optional else "float"
    return "text"


def _model_specs(section: str, model: type, names: Optional[List[str]] = None) -> List[ParamSpec]:
    specs = []
    for name, info in model.model_fields.items():
        if names is not None and name not in names:
            continue
        key = f"{section}.{name}"
        kind = "choice" if key in _CHOICES else _kind_from_annotation(info.annotation)
        default = info.get_default(call_default_factory=True)
        specs.append(
            ParamSpec(
                section=section,
                name=name,
                label=" ".join("elements" if w == "els" else w for w in name.split("_")),
                kind=kind,
                default=default,
                choices=_CHOICES.get(key),
                help=_HELP.get(key, info.description or ""),
            )
        )
    if names is not None:
        specs.sort(key=lambda s: names.index(s.name))
    return specs


def _finish_specs(specs: List[ParamSpec]) -> List[ParamSpec]:
    for s in specs:
        if not s.help:
            s.help = _HELP.get(s.key, "")
    return specs


def build_param_specs() -> List[ParamSpec]:
    """Parameters of the analysis form, in display order."""
    specs: List[ParamSpec] = [
        ParamSpec("filter", "max_analytical_error_percent", "max analytical error (w%)", "float_opt", 5.0,
                  placeholder="none"),
        ParamSpec("filter", "min_bckgrnd_cnts", "min background counts", "float_opt", 5.0, placeholder="saved"),
        ParamSpec("filter", "quant_flags_accepted", "accepted quant flags", "flags", [0, -1],
                  choices=list(QUANT_FLAG_MEANINGS)),
    ]
    specs += _model_specs(
        "clust", ClusteringConfig,
        ["method", "geometry", "features", "k_forced", "k_finding_method", "max_k",
         "auto_merge_clusters", "ref_formulae", "do_matrix_decomposition"],
    )
    for section, model in _SUBMODEL_SECTIONS.items():
        specs += _model_specs(section, model, _SUBMODEL_FIELDS.get(section))
    specs += _model_specs(
        "plot", PlotConfig,
        ["els_to_plot", "els_excluded_clust_plot", "show_unused_comps_clust",
         "show_legend_clustering", "plot_best_mixture"],
    )
    specs.append(ParamSpec("plot", "intermixing_violin_plots", "intermixing violin plots", "choice", SAVED,
                           choices=[SAVED, YES, NO]))
    for s in specs:
        if s.key in ("clust.ref_formulae",):
            s.kind = "formulae"
        elif s.key in ("plot.els_to_plot", "plot.els_excluded_clust_plot"):
            s.kind = "elements"
        if s.key == "clust.k_forced":
            s.placeholder = "auto"
        elif s.key == "dbscan.eps":
            s.placeholder = "default"
    return _finish_specs(specs)


def build_quant_specs() -> List[ParamSpec]:
    """Parameters of the quantification form. Sample/substrate elements are set per sample in the table."""
    tri = [SAVED, YES, NO]
    return _finish_specs([
        ParamSpec("quant", "min_bckgrnd_cnts", "min background counts", "float_opt", None, placeholder="saved"),
        ParamSpec("quant", "min_total_counts_fraction", "min total counts fraction", "float_opt", None,
                  placeholder="saved"),
        ParamSpec("quant", "spectrum_lims", "spectrum limits (channels)", "pair_opt", None, placeholder="saved"),
        ParamSpec("quant", "fit_tolerance", "fit tolerance", "float_opt", None, placeholder="saved"),
        ParamSpec("quant", "use_project_specific_std_dict", "project-specific standards", "choice", SAVED,
                  choices=tri),
        ParamSpec("quant", "is_known_phase_mixture", "known phase mixture", "choice", SAVED, choices=tri),
        ParamSpec("quant", "identify_missing_elements", "identify missing elements", "choice", SAVED, choices=tri),
        ParamSpec("qrun", "interrupt_fits_bad_spectra", "interrupt fits of bad spectra", "bool", True),
        ParamSpec("qrun", "force_requantification", "force requantification (new run)", "bool", False),
        ParamSpec("qrun", "requantify_only_unquantified_spectra", "re-quantify only unquantified spectra",
                  "bool", False),
        ParamSpec("qrun", "max_spectra_to_quantify", "max spectra per sample", "int_opt", None, placeholder="all"),
        ParamSpec("qrun", "num_CPU_cores", "CPU cores", "int_opt", None, placeholder="auto"),
        ParamSpec("qafter", "run_analysis", "run analysis afterwards", "bool", True),
        ParamSpec("qafter", "max_analytical_error", "max analytical error (w%)", "float", 5.0),
    ])


def build_single_specs() -> List[ParamSpec]:
    """Parameters of the single-spectrum form."""
    qdef = QuantificationOptionsConfig()
    return _finish_specs([
        ParamSpec("single", "els_sample", "known sample elements", "text", ""),
        ParamSpec("single", "els_substrate", "known substrate elements", "text", ", ".join(dflt.substrate_els)),
        ParamSpec("single", "is_standard", "standard of known composition", "bool", False),
        ParamSpec("single", "std_formula", "standard formula", "text", "", placeholder="from ledger"),
        ParamSpec("sfit", "quantify", "quantify", "bool", True),
        ParamSpec("sfit", "identify_elements", "identify missing elements", "bool", dflt.identify_missing_elements),
        ParamSpec("sfit", "is_particle", "particle geometry", "bool", True),
        ParamSpec("sfit", "fit_tol", "fit tolerance", "float", qdef.fit_tolerance),
        ParamSpec("sfit", "spectrum_lims", "spectrum limits (channels)", "pair_opt", None, placeholder="saved"),
        ParamSpec("sfit", "max_undetectable_w_fr", "max undetectable mass fraction", "float", 0.0),
        ParamSpec("sfit", "force_single_iteration", "single iteration", "bool", False),
        ParamSpec("sfit", "interrupt_fits_bad_spectra", "interrupt fit of bad spectra", "bool", False),
        ParamSpec("sfit", "free_area_el_lines", "free-area lines", "text", "", placeholder="e.g. Fe_La"),
    ])


ACQ_SAMPLE_TYPES = [cnst.S_POWDER_SAMPLE_TYPE, cnst.S_BULK_SAMPLE_TYPE, cnst.S_POWDER_CONTINUOUS_SAMPLE_TYPE,
                    cnst.S_BULK_ROUGH_SAMPLE_TYPE]
# Sample types acquired on particles (powder options) or on a grid of spots (bulk options)
ACQ_PARTICLE_TYPES = [cnst.S_POWDER_SAMPLE_TYPE]
ACQ_GRID_TYPES = [cnst.S_BULK_SAMPLE_TYPE, cnst.S_POWDER_CONTINUOUS_SAMPLE_TYPE, cnst.S_BULK_ROUGH_SAMPLE_TYPE]


def build_acq_specs() -> List[ParamSpec]:
    """Parameters of the acquisition form, with the defaults of Run_Acquisition.py."""
    from autoemx.config.runtime_configs import BulkMeasurementConfig, MeasurementConfig, PowderMeasurementConfig

    tri_help = {
        "asample.sample_type": "powder: particles on a substrate; bulk / bulk_rough: flat or rough bulk samples, "
                               "measured on a grid of spots; powder_continuous: continuous powder layer, measured on a grid.",
        "asample.sample_halfwidth": "Half-width of the sample, in mm.",
        "asample.sample_substrate_type": "Substrate on which the sample is mounted.",
        "asample.sample_substrate_shape": "Shape of the substrate (stub).",
        "asample.sample_substrate_width_mm": "Width/diameter of the substrate (e.g. Al stub), in mm.",
        "asample.is_auto_substrate_detection": "Detect the carbon tape automatically (it must appear black on a brighter stub). Only for Ctape.",
        "asample.working_distance": "Approximate working distance (mm) at which the sample is in focus. Autofocus is constrained around it.",
        "asample.working_distance_tolerance": "Maximum accepted deviation of the autofocus from the working distance, in mm.",
        "aacq.beam_energy": "Electron beam energy, in keV. The shipped P/B standards are for 15 kV.",
        "aacq.target_Xsp_counts": "Target number of counts of each spectrum.",
        "aacq.max_XSp_acquisition_time": "Maximum acquisition time per spectrum, in s. Empty = counts / 10000 × 5 s.",
        "aacq.min_n_spectra": "Number of spectra after which convergence is checked (only when quantifying during acquisition).",
        "aacq.max_n_spectra": "Maximum number of spectra collected per sample, when the clustering does not converge.",
        "aacq.is_manual_navigation": "Navigate manually: a window asks you to center each acquisition spot.",
        "aacq.auto_adjust_brightness_contrast": "Adjust brightness and contrast automatically.",
        "aacq.contrast": "Contrast, used when automatic brightness/contrast is off.",
        "aacq.brightness": "Brightness, used when automatic brightness/contrast is off.",
        "asample.els_substrate": "Substrate elements, fitted but not quantified unless also listed among the sample elements. They are never identified as sample elements either: list only elements you do not want quantified.",
        "aacq.n_spectra": "Number of spectra collected per sample.",
        "aimg.saved_images_extension": "Format of the saved SEM images (png: light; tif: lossless, larger).",
        "aimg.annotate_particle_images": "Draw the spectrum spots and scale bar on the particle images. If off, annotated copies can be made with Annotate_Particle_Images.py.",
        "aimg.save_raw_images": "Also save the non-annotated version of annotated images (frames, and particles if annotated).",
        "aacq.quantify_spectra": "Quantify spectra during acquisition, and stop when the clustering converges (between min and max spectra). Not recommended on slow microscope computers.",
        "aquant.use_project_specific_std_dict": "Load the P/B standards from the results folder.",
        "aquant.is_known_phase_mixture": "The samples are mixtures of known phases (their candidate phases): quantified with the P/B values of their standards, and their intermixing is characterized.",
        "aquant.interrupt_fits_bad_spectra": "Stop fitting spectra expected to give gross quantification errors.",
        "aquant.max_analytical_error_percent": "Maximum analytical error of the spectra used for clustering (can be changed later).",
        "aquant.min_bckgrnd_cnts": "Minimum background counts under the reference peaks for a spectrum to be used.",
        "aquant.quant_flags_accepted": "Quantification flags kept for clustering (can be changed later).",
        "aquant.max_n_clusters": "Maximum number of clusters (can be changed later).",
        "aquant.show_unused_comps_clust": "Show discarded compositions in the clustering plot.",
    }
    for key, text in tri_help.items():
        _HELP.setdefault(key, text)
    _HELP.setdefault("amicro.microscope_ID", "Microscope used for the acquisition. It sets the driver, the calibrations "
                                             "and the P/B standards used for quantification.")
    specs = [
        ParamSpec("amicro", "microscope_ID", "microscope", "choice", dflt.microscope_ID,
                  choices=available_microscopes()),
        ParamSpec("asample", "sample_type", "sample type", "choice", cnst.S_POWDER_SAMPLE_TYPE, choices=ACQ_SAMPLE_TYPES),
        ParamSpec("asample", "sample_halfwidth", "sample half-width (mm)", "float", 3.0),
        ParamSpec("asample", "sample_substrate_type", "substrate", "choice", cnst.CTAPE_SUBSTRATE_TYPE,
                  choices=[cnst.CTAPE_SUBSTRATE_TYPE, cnst.NONE_SUBSTRATE_TYPE]),
        ParamSpec("asample", "els_substrate", "known substrate elements", "text", ", ".join(dflt.substrate_els)),
        ParamSpec("asample", "sample_substrate_shape", "substrate shape", "choice", cnst.CIRCLE_SUBSTRATE_SHAPE,
                  choices=[cnst.CIRCLE_SUBSTRATE_SHAPE, cnst.SQUARE_SUBSTRATE_SHAPE]),
        ParamSpec("asample", "sample_substrate_width_mm", "substrate width (mm)", "float", 12.0),
        ParamSpec("asample", "is_auto_substrate_detection", "detect carbon tape automatically", "bool", True),
        ParamSpec("asample", "working_distance", "working distance (mm)", "float", 8.5),
        ParamSpec("asample", "working_distance_tolerance", "working distance tolerance (mm)", "float", 1.0),
        ParamSpec("aacq", "beam_energy", "beam energy (keV)", "float", 15.0),
        ParamSpec("aacq", "target_Xsp_counts", "target counts per spectrum", "int", 50000),
        ParamSpec("aacq", "max_XSp_acquisition_time", "max acquisition time (s)", "float_opt", None, placeholder="auto"),
        ParamSpec("aacq", "quantify_spectra", "quantify during acquisition", "bool", False),
        # Without quantification: a fixed number of spectra; with it: min and max (stop when converged)
        ParamSpec("aacq", "n_spectra", "number of spectra", "int", 100),
        ParamSpec("aacq", "min_n_spectra", "min spectra", "int", 50),
        ParamSpec("aacq", "max_n_spectra", "max spectra", "int", 100),
        ParamSpec("aacq", "is_manual_navigation", "manual navigation", "bool", False),
        ParamSpec("aacq", "auto_adjust_brightness_contrast", "automatic brightness/contrast", "bool", True),
        ParamSpec("aacq", "contrast", "contrast", "float_opt", None),
        ParamSpec("aacq", "brightness", "brightness", "float_opt", None),
        ParamSpec("aimg", "saved_images_extension", "image format", "choice", dflt.saved_images_extension,
                  choices=list(MeasurementConfig.ALLOWED_IMAGE_EXTENSIONS)),
        ParamSpec("aimg", "annotate_particle_images", "annotate particle images", "bool", dflt.annotate_particle_images),
        ParamSpec("aimg", "save_raw_images", "save raw copies of annotated images", "bool", dflt.save_raw_images),
    ]
    powder_choices = {
        "par_selection_mode": ["auto", "manual"],  # 'list' needs particle ids/coordinates (scripts only)
        "par_segmentation_model": list(PowderMeasurementConfig.AVAILABLE_PAR_SEGMENTATION_MODELS),
        "par_feature_selection": list(PowderMeasurementConfig.AVAILABLE_FEATURE_SELECTION),
        "par_spot_spacing": list(PowderMeasurementConfig.AVAILABLE_SPOT_SPACING_SELECTION),
    }
    for key, values in powder_choices.items():
        _CHOICES[f"apowder.{key}"] = values
    # par_spot_selection_mode 'callback' needs a Python callback (scripts only)
    powder = _model_specs("apowder", PowderMeasurementConfig,
                          [n for n in PowderMeasurementConfig.model_fields if n != "par_spot_selection_mode"])
    script_defaults = {"max_area_par": 10000.0}  # Run_Acquisition.py defaults that differ from the config model
    for spec in powder:
        spec.default = script_defaults.get(spec.name, spec.default)
        if spec.name == "par_search_frame_width_um":
            spec.placeholder = "auto"
    bulk = _model_specs("abulk", BulkMeasurementConfig)
    for spec in bulk:
        if spec.name == "image_frame_width_um":
            spec.placeholder = "10 × grid"
    # Readable labels (field names abbreviate "particle" as "par") and units shown next to the fields
    labels = {
        "par_selection_mode": ("particle selection mode", ""),
        "img_shift_tracking": ("image drift tracking", ""),
        "par_search_frame_width_um": ("particle search frame width", "µm"),
        "max_n_par_per_frame": ("max particles per frame", ""),
        "max_spectra_per_par": ("max spectra per particle", ""),
        "max_area_par": ("max particle area", "µm²"),
        "min_area_par": ("min particle area", "µm²"),
        "par_mask_margin": ("particle mask margin", "µm"),
        "xsp_spots_distance_um": ("min distance between spectrum spots", "µm"),
        "par_segmentation_model": ("particle segmentation model", ""),
        "par_brightness_thresh": ("particle brightness threshold", ""),
        "par_xy_spots_thresh": ("particle spot brightness threshold", ""),
        "par_feature_selection": ("particle spot selection", ""),
        "par_spot_spacing": ("particle spot spacing", ""),
        "grid_spot_spacing_um": ("grid spot spacing", "µm"),
        "min_xsp_spots_distance_um": ("min distance between spectrum spots", "µm"),
        "image_frame_width_um": ("image frame width", "µm"),
    }
    helps = {
        "par_selection_mode": "auto: frames are scanned and particles detected automatically; manual: a window asks you to centre each particle.",
        "img_shift_tracking": "Track and correct the image drift while acquiring spectra on a particle.",
        "par_search_frame_width_um": "Width of the frames in which particles are searched. Empty = min(20 × max particle radius, 500 µm).",
        "max_n_par_per_frame": "Maximum number of particles analysed in a frame, for a representative sampling.",
        "max_spectra_per_par": "Maximum number of spectra acquired on a particle, so that more particles are analysed.",
        "max_area_par": "Particles larger than this are ignored.",
        "min_area_par": "Particles smaller than this are ignored.",
        "par_mask_margin": "Spectra are not acquired closer than this to the particle edge.",
        "xsp_spots_distance_um": "Minimum distance between the spectrum spots on a particle.",
        "par_segmentation_model": "Model used to detect particles in the frames.",
        "par_brightness_thresh": "8-bit brightness above which pixels belong to particles (on a dark substrate).",
        "par_xy_spots_thresh": "8-bit brightness (rescaled within each particle) of the brightest, thickest regions, where spectra are acquired.",
        "par_feature_selection": "random: random spots within the bright regions; peaks: the brightest spots.",
        "par_spot_spacing": "random: unbiased spot selection; maximized: spots spread over the particle.",
        "grid_spot_spacing_um": "Distance between the spots of the acquisition grid.",
        "min_xsp_spots_distance_um": "Offset of the grid when it does not contain enough spots.",
        "image_frame_width_um": "Width of the image frames. Empty = 10 × grid spot spacing.",
        "randomize_frames": "Acquire the grid spots in random order.",
        "exclude_sample_margin": "Exclude the margin of the sample (e.g. when contaminated).",
    }
    for spec in powder + bulk:
        spec.label, spec.unit = labels.get(spec.name, (spec.label, ""))
        spec.help = spec.help or helps.get(spec.name, "")
    specs += powder + bulk
    specs += [
        ParamSpec("aquant", "use_project_specific_std_dict", "project-specific standards", "bool", False),
        ParamSpec("aquant", "is_known_phase_mixture", "known phase mixture", "bool", False),
        ParamSpec("aquant", "interrupt_fits_bad_spectra", "interrupt fits of bad spectra", "bool", True),
        ParamSpec("aquant", "max_analytical_error_percent", "max analytical error (w%)", "float", 5.0),
        ParamSpec("aquant", "min_bckgrnd_cnts", "min background counts", "float", 5.0),
        ParamSpec("aquant", "quant_flags_accepted", "accepted quant flags", "flags", [0, -1],
                  choices=list(QUANT_FLAG_MEANINGS)),
        ParamSpec("aquant", "max_n_clusters", "max clusters", "int", 6),
        ParamSpec("aquant", "show_unused_comps_clust", "show discarded compositions", "bool", True),
    ]
    return _finish_specs(specs)


PARAM_SPECS: List[ParamSpec] = build_param_specs()
QUANT_SPECS: List[ParamSpec] = build_quant_specs()
SINGLE_SPECS: List[ParamSpec] = build_single_specs()
ACQ_SPECS: List[ParamSpec] = build_acq_specs()
SPECS_BY_KEY: Dict[str, ParamSpec] = {s.key: s for s in PARAM_SPECS + QUANT_SPECS + SINGLE_SPECS + ACQ_SPECS}


def sample_param_values(info: SampleInfo, analysis: Optional[AnalysisRef] = None) -> Dict[str, Any]:
    """Analysis parameter values saved in the ledger (from ``analysis``, or the active analysis)."""
    ledger = info.ledger
    values = {s.key: s.default for s in PARAM_SPECS}
    if analysis is None:
        analysis = info.active_analysis
    cfg = analysis.config if analysis is not None else ClusteringConfig()
    values["filter.max_analytical_error_percent"] = cfg.max_analytical_error_percent
    values["filter.quant_flags_accepted"] = list(cfg.quant_flags_accepted)
    values["filter.min_bckgrnd_cnts"] = cfg.min_bckgrnd_cnts
    for name in ("method", "geometry", "features", "k_forced", "max_k",
                 "auto_merge_clusters", "ref_formulae", "do_matrix_decomposition"):
        values[f"clust.{name}"] = getattr(cfg, name)
    values["clust.k_finding_method"] = (
        cfg.k_finding_method if cfg.k_finding_method != "forced"
        else ClusteringConfig.model_fields["k_finding_method"].default
    )
    for section, attr in (("dbscan", "dbscan"), ("aitchison", "aitchison"),
                          ("merge", "cluster_merge"), ("mixture", "mixture")):
        for name, value in getattr(cfg, attr).model_dump().items():
            if f"{section}.{name}" in values:  # fields shown in the form
                values[f"{section}.{name}"] = value
    plot_cfg = ledger.configs.plot_cfg
    for name in ("els_to_plot", "els_excluded_clust_plot", "show_unused_comps_clust",
                 "show_legend_clustering", "plot_best_mixture"):
        values[f"plot.{name}"] = getattr(plot_cfg, name)
    return values


def quant_options(info: SampleInfo) -> QuantificationOptionsConfig:
    """Options of the active quantification run (defaults if never quantified)."""
    active_q = _active_quant_config(info.ledger)
    try:
        return QuantificationOptionsConfig(**(active_q.options if active_q is not None else {}))
    except Exception:
        return QuantificationOptionsConfig()


def single_param_values(info: Optional[SampleInfo]) -> Dict[str, Any]:
    """Single-spectrum parameter values: defaults, or those of the sample's active quantification."""
    values = {s.key: s.default for s in SINGLE_SPECS}
    if info is None:
        return values
    opts = quant_options(info)
    values["single.els_sample"] = ", ".join(info.elements)
    values["single.els_substrate"] = ", ".join(info.substrate)
    values["sfit.is_particle"] = bool(info.ledger.configs.sample_cfg.is_surface_rough)
    values["sfit.fit_tol"] = opts.fit_tolerance
    values["sfit.spectrum_lims"] = [int(round(v)) for v in opts.spectrum_lims]
    values["sfit.identify_elements"] = bool(getattr(opts, "identify_missing_elements", False))
    exp_cfg = info.ledger.configs.measurement_cfg.exp_stds_cfg
    values["single.is_standard"] = bool(exp_cfg is not None and getattr(exp_cfg, "is_exp_std_measurement", False))
    return values


def peak_overlaps(
    meas_type: Optional[str],
    els_sample: List[str],
    els_substrate: List[str],
    beam_energy_keV: Optional[float],
    det_ch_offset: Optional[float] = None,
    det_ch_width: Optional[float] = None,
    spectrum_lims: Optional[Any] = None,
    microscope_ID: Optional[str] = None,
    free_area_el_lines: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """Reference peaks of sample elements overlapping with peaks of free area (see autoemx.core.quantifier.peak_overlaps)."""
    from autoemx.core.quantifier.peak_overlaps import get_peak_overlaps

    lims = spectrum_lims or dflt.spectrum_lims
    energy_range = None
    if det_ch_offset is not None and det_ch_width:
        energy_range = (det_ch_offset + det_ch_width * lims[0], det_ch_offset + det_ch_width * lims[1])
    els_sample = [el for el in els_sample if el not in _DEFAULT_UNDETECTABLE_ELS]
    els_substrate = [el for el in els_substrate if el not in _DEFAULT_UNDETECTABLE_ELS and el not in els_sample]
    try:
        return get_peak_overlaps(meas_type or dflt.measurement_type, els_sample, els_substrate,
                                 float(beam_energy_keV) if beam_energy_keV else None, energy_range,
                                 microscope_ID or dflt.microscope_ID, free_area_el_lines)
    except Exception:
        return []


def sample_peak_overlaps(info: SampleInfo, els_sample: Optional[List[str]] = None,
                         els_substrate: Optional[List[str]] = None,
                         spectrum_lims: Optional[Any] = None,
                         free_area_el_lines: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """Peak overlaps for a sample, with its saved elements unless others are given."""
    cfgs = info.ledger.configs
    return peak_overlaps(
        cfgs.measurement_cfg.type,
        info.elements if els_sample is None else els_sample,
        info.substrate if els_substrate is None else els_substrate,
        cfgs.measurement_cfg.beam_energy_keV,
        cfgs.microscope_cfg.energy_zero,
        cfgs.microscope_cfg.bin_width,
        spectrum_lims or quant_options(info).spectrum_lims,
        cfgs.microscope_cfg.ID,
        free_area_el_lines,
    )


def _elements_list(text: Any) -> List[str]:
    if text is None:
        return []
    if isinstance(text, (list, tuple)):
        return [str(t) for t in text if str(t).strip()]
    from autoemx.web.pipeline import parse_elements

    return parse_elements(str(text))


def _formulae_list(text: Any) -> List[str]:
    if text is None:
        return []
    if isinstance(text, (list, tuple)):
        items = text
    else:
        items = str(text).replace("\n", ",").replace(";", ",").split(",")
    return [t.strip() for t in items if str(t).strip()]


def coerce_values(raw: Dict[str, Any], specs: Optional[List[ParamSpec]] = None) -> Dict[str, Any]:
    """Convert GUI values of a form (default: analysis) to typed values. Raises ``ValueError`` listing every invalid field."""
    out: Dict[str, Any] = {}
    errors: List[str] = []
    for spec in (specs if specs is not None else PARAM_SPECS):
        key = spec.key
        value = raw.get(key, spec.default)
        try:
            if spec.kind == "bool":
                out[key] = bool(value) if not isinstance(value, list) else bool(value)
            elif spec.kind in ("int", "int_opt"):
                if value in (None, ""):
                    if spec.kind == "int":
                        raise ValueError("a value is required")
                    out[key] = None
                else:
                    f = float(value)
                    if f != int(f):
                        raise ValueError("must be an integer")
                    out[key] = int(f)
            elif spec.kind in ("float", "float_opt"):
                if value in (None, ""):
                    if spec.kind == "float":
                        raise ValueError("a value is required")
                    out[key] = None
                else:
                    out[key] = float(value)
            elif spec.kind in ("pair", "pair_opt"):
                if spec.kind == "pair_opt" and (value is None or str(value).strip() in ("", "[]")):
                    out[key] = None
                    continue
                vals = value if isinstance(value, (list, tuple)) else str(value).replace(";", ",").split(",")
                vals = [int(float(v)) for v in vals if v is not None and str(v).strip() != ""]
                if spec.kind == "pair_opt" and not vals:
                    out[key] = None
                    continue
                if len(vals) != 2 or vals[0] >= vals[1]:
                    raise ValueError("enter both the minimum and the maximum, with min < max")
                out[key] = vals
            elif spec.kind == "flags":
                out[key] = sorted({int(v) for v in (value or [])})
            elif spec.kind == "formulae":
                out[key] = _formulae_list(value)
            elif spec.kind == "elements":
                out[key] = _elements_list(value)
            elif spec.key in ("single.els_sample", "single.els_substrate", "asample.els_substrate"):
                out[key] = _elements_list(value)
            elif spec.key == "single.std_formula":
                out[key] = _valid_formula(value) if str(value or "").strip() else None
            elif spec.key == "sfit.free_area_el_lines":
                out[key] = [t.strip() for t in str(value or "").replace(";", ",").split(",") if t.strip()] or None
            else:
                out[key] = value
        except Exception as exc:
            errors.append(f"{SECTION_TITLES[spec.section]} › {spec.label}: {exc}")
    errors += _validate_models(out)
    if errors:
        raise ValueError("; ".join(errors))
    return out


def _section_dict(values: Dict[str, Any], section: str) -> Dict[str, Any]:
    prefix = section + "."
    return {k[len(prefix):]: v for k, v in values.items() if k.startswith(prefix)}


def _validate_models(values: Dict[str, Any]) -> List[str]:
    errors = []
    for section, model in _SUBMODEL_SECTIONS.items():
        try:
            section_values = _section_dict(values, section)
            if not section_values:
                continue  # not part of this form
            if len(section_values) < len(_SUBMODEL_FIELDS.get(section) or model.model_fields):
                continue  # a field failed coercion, already reported
            model.model_validate(section_values)
        except Exception as exc:
            errors.append(f"{SECTION_TITLES[section]}: {_pydantic_msg(exc)}")
    from pymatgen.core import Composition, Element

    for formula in values.get("clust.ref_formulae", []):
        try:
            valid = all(isinstance(el, Element) for el in Composition(formula).elements)
        except Exception:
            valid = False
        if not valid:
            errors.append(f"Invalid candidate formula '{formula}'")
    if "single.els_sample" in values and not values["single.els_sample"]:
        errors.append("At least one sample element is required")
    if values.get("filter.quant_flags_accepted") == []:
        errors.append("Accept at least one quant flag")
    return errors


def _valid_formula(text: Any) -> str:
    from pymatgen.core import Composition, Element

    formula = str(text).strip()
    try:
        valid = all(isinstance(el, Element) for el in Composition(formula).elements)
    except Exception:
        valid = False
    if not valid:
        raise ValueError(f"invalid formula '{formula}'")
    return formula


def _pydantic_msg(exc: Exception) -> str:
    errs = getattr(exc, "errors", None)
    if callable(errs):
        return "; ".join(f"{'.'.join(map(str, e['loc']))}: {e['msg']}" for e in errs())
    return str(exc)


def analysis_kwargs(values: Dict[str, Any]) -> Dict[str, Any]:
    """Keyword arguments for ``analyze_sample`` from coerced GUI values."""
    k_forced = values["clust.k_forced"]
    kwargs = dict(
        ref_formulae=values["clust.ref_formulae"],
        clustering_features=values["clust.features"],
        clustering_method=values["clust.method"],
        clustering_geometry=values["clust.geometry"],
        dbscan_params=_section_dict(values, "dbscan"),
        aitchison_params=_section_dict(values, "aitchison"),
        k_forced=k_forced if k_forced is not None else False,
        k_finding_method=None if k_forced is not None else values["clust.k_finding_method"],
        max_k=values["clust.max_k"],
        auto_merge_clusters=values["clust.auto_merge_clusters"],
        cluster_merge_params=_section_dict(values, "merge"),
        do_matrix_decomposition=values["clust.do_matrix_decomposition"],
        mixture_params=_section_dict(values, "mixture"),
        max_analytical_error_percent=values["filter.max_analytical_error_percent"],
        quant_flags_accepted=values["filter.quant_flags_accepted"],
        min_bckgrnd_cnts=values["filter.min_bckgrnd_cnts"],
        els_to_plot=values["plot.els_to_plot"],
        els_excluded_clust_plot=values["plot.els_excluded_clust_plot"],
        show_unused_compositions_cluster_plot=values["plot.show_unused_comps_clust"],
        show_legend_clustering=values["plot.show_legend_clustering"],
        plot_best_mixture=values["plot.plot_best_mixture"],
        intermixing_violin_plots=_TRISTATE[values["plot.intermixing_violin_plots"]],
        show_plots=False,
        plot_custom_plots=False,
    )
    return kwargs


def quantification_kwargs(values: Dict[str, Any]) -> Dict[str, Any]:
    """Keyword arguments for ``batch_quantify_and_analyze`` from coerced quantification-form values.

    Options left empty (or 'saved') are None, i.e. each sample keeps its saved value.
    """
    lims = values["quant.spectrum_lims"]
    return dict(
        quantification_method="PB",
        min_bckgrnd_cnts=values["quant.min_bckgrnd_cnts"],
        min_total_counts_fraction=values["quant.min_total_counts_fraction"],
        spectrum_lims=tuple(lims) if lims else None,
        fit_tolerance=values["quant.fit_tolerance"],
        use_project_specific_std_dict=_TRISTATE[values["quant.use_project_specific_std_dict"]],
        is_known_phase_mixture=_TRISTATE[values["quant.is_known_phase_mixture"]],
        identify_missing_elements=_TRISTATE[values["quant.identify_missing_elements"]],
        interrupt_fits_bad_spectra=values["qrun.interrupt_fits_bad_spectra"],
        force_requantification=values["qrun.force_requantification"],
        requantify_only_unquantified_spectra=values["qrun.requantify_only_unquantified_spectra"],
        max_spectra_to_quantify=values["qrun.max_spectra_to_quantify"],
        num_CPU_cores=values["qrun.num_CPU_cores"],
        run_analysis=values["qafter.run_analysis"],
        max_analytical_error=values["qafter.max_analytical_error"],
    )


def acquisition_samples(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sample definitions for ``batch_acquire_and_analyze`` from the rows of the samples table.

    Raises ``ValueError`` listing every invalid entry.
    """
    samples, errors, ids = [], [], set()
    for i, row in enumerate(rows or [], start=1):
        sid = str(row.get("ID") or "").strip()
        if not sid and not any(str(row.get(k) or "").strip() for k in ("els", "x", "y", "cnd")):
            continue  # empty row
        name = sid or f"row {i}"
        if not sid:
            errors.append(f"{name}: missing sample ID")
        elif sid in ids:
            errors.append(f"{name}: duplicated sample ID")
        elif any(c in sid for c in '/\\:*?"<>|'):
            errors.append(f"{name}: the sample ID cannot contain / \\ : * ? \" < > |")
        ids.add(sid)
        try:
            els = _elements_list(row.get("els"))
            if not els:
                errors.append(f"{name}: no elements")
        except ValueError as exc:
            errors.append(f"{name}: {exc}")
            els = []
        try:
            pos = (float(row.get("x")), float(row.get("y")))
        except (TypeError, ValueError):
            errors.append(f"{name}: position x, y (mm) must be numbers")
            pos = (0.0, 0.0)
        cnd = _formulae_list(row.get("cnd"))
        for formula in cnd:
            try:
                _valid_formula(formula)
            except ValueError as exc:
                errors.append(f"{name}: {exc}")
        samples.append({"ID": sid, "els": els, "pos": pos, "cnd": cnd})
    if not samples and not errors:
        errors.append("Add at least one sample to acquire")
    if errors:
        raise ValueError("; ".join(errors))
    return samples


def acquisition_kwargs(values: Dict[str, Any]) -> Dict[str, Any]:
    """Keyword arguments for ``batch_acquire_and_analyze`` from coerced acquisition-form values."""
    max_time = values["aacq.max_XSp_acquisition_time"]
    if max_time is None:
        max_time = values["aacq.target_Xsp_counts"] / 10000 * 5  # as in Run_Acquisition.py
    kwargs = {name: values[f"asample.{name}"] for name in (
        "sample_type", "sample_halfwidth", "sample_substrate_type", "sample_substrate_shape",
        "sample_substrate_width_mm", "is_auto_substrate_detection", "working_distance", "working_distance_tolerance")}
    quantify = values["aacq.quantify_spectra"]
    kwargs.update(
        microscope_ID=values["amicro.microscope_ID"],
        beam_energy=values["aacq.beam_energy"],
        target_Xsp_counts=values["aacq.target_Xsp_counts"],
        max_XSp_acquisition_time=max_time,
        # Without quantification, "number of spectra" are collected: also the minimum, since the analyser
        # raises max_n_spectra to min_n_spectra (the hidden "min spectra" field would set the number)
        min_n_spectra=values["aacq.min_n_spectra"] if quantify else values["aacq.n_spectra"],
        max_n_spectra=values["aacq.max_n_spectra"] if quantify else values["aacq.n_spectra"],
        is_manual_navigation=values["aacq.is_manual_navigation"],
        auto_adjust_brightness_contrast=values["aacq.auto_adjust_brightness_contrast"],
        contrast=values["aacq.contrast"],
        brightness=values["aacq.brightness"],
        els_substrate=values["asample.els_substrate"],
        saved_images_extension=values["aimg.saved_images_extension"],
        annotate_particle_images=values["aimg.annotate_particle_images"],
        save_raw_images=values["aimg.save_raw_images"],
        powder_meas_cfg_kwargs=_section_dict(values, "apowder"),
        bulk_meas_cfg_kwargs=_section_dict(values, "abulk"),
        quantify_spectra=quantify,
        use_project_specific_std_dict=values["aquant.use_project_specific_std_dict"],
        is_known_phase_mixture=values["aquant.is_known_phase_mixture"],
        interrupt_fits_bad_spectra=values["aquant.interrupt_fits_bad_spectra"],
        max_analytical_error_percent=values["aquant.max_analytical_error_percent"],
        min_bckgrnd_cnts=values["aquant.min_bckgrnd_cnts"],
        quant_flags_accepted=values["aquant.quant_flags_accepted"],
        max_n_clusters=values["aquant.max_n_clusters"],
        show_unused_comps_clust=values["aquant.show_unused_comps_clust"],
    )
    return kwargs


def acquisition_script(samples: List[Dict[str, Any]], kwargs: Dict[str, Any], results_dir: Optional[str]) -> str:
    """A runnable Python script performing the acquisition with these settings (like Run_Acquisition.py)."""
    import pprint
    from datetime import datetime

    lines = [
        "#!/usr/bin/env python3",
        "# -*- coding: utf-8 -*-",
        '"""',
        f"Automated acquisition of X-ray spectra, exported from the AutoEMX GUI on {datetime.now():%Y-%m-%d %H:%M}.",
        "",
        "Run this script on the microscope computer to acquire the samples below with the same settings.",
        '"""',
        "",
        "samples = " + pprint.pformat(samples, sort_dicts=False, width=100),
        "",
        f"results_dir = {results_dir!r} # Folder where a sub-folder is created for each sample",
        "",
        "from autoemx.runners.batch_acquire_and_analyze import batch_acquire_and_analyze",
        "",
        "comp_analyzer = batch_acquire_and_analyze(",
        "    samples=samples,",
    ]
    for key, value in kwargs.items():
        lines.append(f"    {key}={pprint.pformat(value, sort_dicts=False, width=100)},")
    lines += ["    development_mode=False,", "    verbose=True,", "    results_dir=results_dir,", ")", ""]
    return "\n".join(lines)


def acquisition_progress(log_text: str, sample_ids: List[str], max_n: int) -> List[Dict[str, Any]]:
    """Progress of each sample of an acquisition (spectra acquired, or why it failed), parsed from its log."""
    progress = [{"sample": sid, "state": "waiting", "done": 0, "total": max_n, "error": None} for sid in sample_ids]
    positions = []
    for i, sid in enumerate(sample_ids):
        idx = log_text.find(f"Sample '{sid}'")
        if idx >= 0:
            positions.append((idx, i))
    positions.sort()
    for k, (start, i) in enumerate(positions):
        end = positions[k + 1][0] if k + 1 < len(positions) else len(log_text)
        progress[i]["state"] = "running"
        progress[i]["done"] = len(re.findall(r"Acquiring spectrum #\d+", log_text[start:end]))
        # Logged by batch_acquire_and_analyze when a sample fails (the batch goes on with the next one)
        failed = re.search(rf"Sample '{re.escape(progress[i]['sample'])}': acquisition/quantification failed: (.*)",
                           log_text[start:end])
        if failed:
            progress[i]["state"], progress[i]["error"] = "failed", failed.group(1).strip()
    for _start, i in positions[:-1]:
        if progress[i]["state"] == "running":
            progress[i]["state"] = "done"
    return progress


def microscope_status(microscope_id: str = dflt.microscope_ID) -> Tuple[bool, str]:
    """Whether the API of the microscope driver can be imported (without connecting to the microscope)."""
    import importlib.util

    api_modules = {"PhenomXL": "PyPhenom"}
    module = api_modules.get(microscope_id)
    if module is None:
        return True, f"{microscope_id} driver"
    if importlib.util.find_spec(module) is None:
        return False, f"{microscope_id}: {module} not installed, acquisition will fail on this computer"
    return True, f"{microscope_id}: {module} available"


def single_fit_kwargs(values: Dict[str, Any]) -> Dict[str, Any]:
    """Options of a single-spectrum fit from coerced single-spectrum-form values."""
    w_frs = None
    if values["single.is_standard"] and values["single.std_formula"]:
        from pymatgen.core import Composition

        from autoemx.utils import composition_as_weight_dict

        w_frs = composition_as_weight_dict(Composition(values["single.std_formula"]))
    lims = values["sfit.spectrum_lims"]
    return dict(
        els_sample=values["single.els_sample"],
        els_substrate=values["single.els_substrate"],
        is_standard=values["single.is_standard"],
        els_w_frs=w_frs,
        quantify=values["sfit.quantify"],
        identify_missing_elements=values["sfit.identify_elements"],
        is_particle=values["sfit.is_particle"],
        fit_tol=values["sfit.fit_tol"],
        spectrum_lims=tuple(lims) if lims else None,
        max_undetectable_w_fr=values["sfit.max_undetectable_w_fr"],
        force_single_iteration=values["sfit.force_single_iteration"],
        interrupt_fits_bad_spectra=values["sfit.interrupt_fits_bad_spectra"],
        free_area_el_lines=values["sfit.free_area_el_lines"],
    )


# =============================================================================
# Job runner (separate process, log captured to file)
# =============================================================================
def _job_entry(kind: str, sample_dir: str, payload: Dict[str, Any], log_path: str, result_path: str) -> None:
    """Child-process entry point."""
    if hasattr(os, "setpgrp"):
        os.setpgrp()  # own process group, so cancelling also stops worker pools
    os.environ["MPLBACKEND"] = "Agg"
    log_fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_APPEND)
    os.dup2(log_fd, 1)
    os.dup2(log_fd, 2)
    sys.stdout = os.fdopen(1, "w", buffering=1, encoding="utf-8", errors="replace", closefd=False)
    sys.stderr = os.fdopen(2, "w", buffering=1, encoding="utf-8", errors="replace", closefd=False)
    import logging

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
        force=True,
    )
    import matplotlib

    matplotlib.use("Agg")
    result: Dict[str, Any] = {"ok": False}
    try:
        result = _JOB_FUNCS[kind](sample_dir, payload)
        result.setdefault("ok", True)
    except Exception as exc:
        traceback.print_exc()
        result = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
    finally:
        sys.stdout.flush()
        tmp = result_path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(result, fh, default=_json_default)
        os.replace(tmp, result_path)


def _json_default(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return str(obj)


def _run_quant_and_analysis(sample_dir: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    if payload.get("quantify"):
        from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze

        print("=" * 60 + "\nQUANTIFICATION\n" + "=" * 60, flush=True)
        res = batch_quantify_and_analyze(
            sample_IDs=[Path(sample_dir).name],
            results_path=str(Path(sample_dir).parent),
            **payload["quant_kwargs"],
        )
        if not res:
            raise RuntimeError("Quantification failed. See the log for details.")
        out["quantified"] = True
    if payload.get("analyse", True):
        from autoemx.runners.analyze_sample import analyze_sample

        print("=" * 60 + "\nANALYSIS\n" + "=" * 60, flush=True)
        analyzer = analyze_sample(
            sample_ID=Path(sample_dir).name,
            results_path=str(Path(sample_dir).parent),
            **payload["analysis_kwargs"],
        )
        if analyzer is None:
            raise RuntimeError("Analysis failed. See the log for details.")
        out["analysis_dir"] = getattr(analyzer, "analysis_dir", None)
        out["n_clusters"] = (getattr(analyzer, "clustering_info", None) or {}).get("n_clusters")
        if out["n_clusters"] is None:
            out["warning"] = "Analysis did not produce clusters (too few valid spectra?). See the log."
    return out


def _ledger_quant_flag(quantifier: Any, els_sample: List[str], els_substrate: List[str],
                       min_bckgrnd_cnts: Optional[float]) -> Optional[int]:
    """Quant flag of a quantified spectrum, as assigned during batch quantification (see quant_flag docs)."""
    from autoemx.core.quantifier.quant_flags import check_fit_quant_validity

    if getattr(quantifier, "quant_result", None) is None and getattr(quantifier, "bad_quant_flag", None) is None:
        return None
    detectable = [el for el in els_substrate if el not in _DEFAULT_UNDETECTABLE_ELS and el not in els_sample]
    try:
        flag, _comment = check_fit_quant_validity(
            is_quant_fit_valid=getattr(quantifier, "quant_result", None) is not None,
            bad_quant_flag=getattr(quantifier, "bad_quant_flag", None),
            quantifier=quantifier,
            min_bckgrnd_ref_lines=quantifier._get_min_bckgrnd_cnts_ref_quant_lines(),
            detectable_els_substrate=detectable,
            min_bckgrnd_cnts=min_bckgrnd_cnts,
        )
        return flag
    except Exception:
        traceback.print_exc()
        return None


def _fit_payload(quantifier: Any, quant_flag: Optional[int] = None) -> Dict[str, Any]:
    """Plot series, composition and fitted peaks of a fitted (and quantified) spectrum."""
    from autoemx.web.pipeline import _peak_labels, _plot_series

    energy, counts, fitted, background = _plot_series(quantifier)
    quant = getattr(quantifier, "quant_result", None) or {}
    try:
        bckgrnd_cnts = [[line, e, v] for line, (e, v) in quantifier.get_bckgrnd_cnts_ref_lines().items()]
    except Exception:
        bckgrnd_cnts = []
    peaks = []
    ref_lines = set(getattr(quantifier, "ref_lines_for_quant", []) or [])
    for el_line, peak in (getattr(quantifier, "fitted_peaks_info", None) or {}).items():
        peaks.append({
            "line": el_line,
            "reference": el_line in ref_lines,
            "center": peak.get(cnst.PEAK_CENTER_KEY),
            "th_energy": peak.get(cnst.PEAK_TH_ENERGY_KEY),
            "fwhm": peak.get(cnst.PEAK_FWHM_KEY),
            "area": peak.get(cnst.PEAK_AREA_KEY),
            "height": peak.get(cnst.PEAK_HEIGHT_KEY),
            "pb_ratio": peak.get(cnst.PB_RATIO_KEY),
        })
    # Reference lines (used for quantification) first, then by energy
    peaks.sort(key=lambda p: (not p["reference"], p["center"] is None, p["center"] or 0))
    fit_result = getattr(quantifier, "fit_result", None)
    return {
        "energy": energy,
        "counts": counts,
        "fit": fitted,
        "background": background,
        "peak_labels": _peak_labels(quantifier),
        "bckgrnd_cnts": bckgrnd_cnts,  # [el_line, energy, background counts] under each reference line
        "peaks": peaks,
        "comp_at": dict(quant.get(cnst.COMP_AT_FR_KEY) or {}),
        "comp_w": dict(quant.get(cnst.COMP_W_FR_KEY) or {}),
        "analytical_error": quant.get(cnst.AN_ER_KEY),
        "r_squared": quant.get(cnst.R_SQ_KEY, getattr(fit_result, "rsquared", None)),
        "redchi_sq": quant.get(cnst.REDCHI_SQ_KEY, getattr(fit_result, "redchi", None)),
        "quant_flag": quant_flag,
        **_element_check_payload(quantifier),
    }


def _element_check_payload(quantifier: Any) -> Dict[str, Any]:
    """Elements added to the fit by the element identification, and those possibly present (empty if not run)."""
    if getattr(quantifier, "added_elements", None) is None:
        return {"element_check": None}
    return {"element_check": {
        "added_quantified": list(quantifier.added_elements or []),
        "added_not_quantified": list(quantifier.added_unquantified_elements or []),
        "possible": list(quantifier.possible_elements or []),
        "note": quantifier.identification_note,
    }}


def element_check_message(check: Optional[Dict[str, Any]]) -> Optional[str]:
    """One-line summary of the element identification of a fit (None if it was not run)."""
    if check is None:
        return None
    parts = []
    if check["added_quantified"]:
        parts.append(f"added and quantified: {', '.join(check['added_quantified'])}")
    if check["added_not_quantified"]:
        parts.append(f"added, fitted but not quantified (no standard): {', '.join(check['added_not_quantified'])}")
    if check["possible"]:
        parts.append(f"possibly present (not added): {', '.join(check['possible'])}")
    msg = "; ".join(parts) if parts else "no missing elements found"
    if check.get("note"):
        msg += f" ({check['note']})"
    return msg


def _run_fit(sample_dir: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Re-fit a spectrum of a sample with the settings of its active quantification (Analysis tab)."""
    from autoemx.runners.fit_and_quantify_spectrum_from_ledger import fit_and_quantify_spectrum_from_ledger

    info = load_sample_info(sample_dir)
    opts = quant_options(info)
    quantifier = fit_and_quantify_spectrum_from_ledger(
        sample_ID=info.sample_id,
        spectrum_ID=int(payload["spectrum_id"]),
        els_sample=info.elements,
        els_substrate=info.substrate,
        results_path=str(Path(sample_dir).parent),
        quantify=True,
        plot_signal=False,
        fit_tol=float(opts.fit_tolerance),
        is_particle=bool(info.ledger.configs.sample_cfg.is_surface_rough),
        use_instrument_background=bool(opts.use_instrument_background),
        identify_missing_elements=bool(getattr(opts, "identify_missing_elements", False)),
        print_results=False,
        quant_verbose=True,
        fitting_verbose=False,
        ledger=info.ledger,
    )
    if quantifier is None or getattr(quantifier, "fit_result", None) is None:
        raise RuntimeError("Fit failed. See the log for details.")
    active = info.active_analysis
    min_cnts = active.config.min_bckgrnd_cnts if active else ClusteringConfig().min_bckgrnd_cnts
    flag = _ledger_quant_flag(quantifier, info.elements, info.substrate, min_cnts)
    return {"spectrum_id": payload["spectrum_id"], **_fit_payload(quantifier, flag)}


def _run_single(sample_dir: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Fit and quantify one spectrum, from a sample ledger or an external EMSA file (Single spectrum tab)."""
    kw = dict(payload["fit_kwargs"])
    quantify = kw.pop("quantify")
    common = dict(
        quantify=quantify,
        plot_signal=False,
        print_results=False,
        quant_verbose=True,
        fitting_verbose=False,
    )
    if payload["source"] == "file":
        import autoemx.config.defaults as defaults
        from autoemx.runners.fit_and_quantify_spectrum import fit_and_quantify_spectrum
        from autoemx.web.pipeline import load_uploaded_spectrum

        path = Path(payload["path"])
        spectrum_vals, _, geometry = load_uploaded_spectrum(path)
        is_standard = kw.pop("is_standard")
        if is_standard and kw["els_w_frs"] is None:
            raise ValueError("Enter the formula of the standard to fit an external spectrum as a standard.")
        print(f"Fitting {path.name}: beam {geometry['beam_energy']} keV, emergence angle "
              f"{geometry['emergence_angle']}°, live time {geometry['sp_collection_time']} s", flush=True)
        quantifier = fit_and_quantify_spectrum(
            spectrum_vals=spectrum_vals,
            spectrum_lims=kw.pop("spectrum_lims") or defaults.spectrum_lims,
            microscope_ID=defaults.microscope_ID,
            meas_type=defaults.measurement_type,
            meas_mode=defaults.measurement_mode,
            det_ch_offset=geometry["det_ch_offset"],
            det_ch_width=geometry["det_ch_width"],
            beam_energy=geometry["beam_energy"],
            emergence_angle=geometry["emergence_angle"],
            sp_collection_time=geometry["sp_collection_time"],
            sample_ID=path.stem,
            background_vals=None,
            **kw,
            **common,
        )
        extra = {"beam_energy_kV": geometry["beam_energy"]}
    else:
        from autoemx.runners.fit_and_quantify_spectrum_from_ledger import fit_and_quantify_spectrum_from_ledger

        info = load_sample_info(sample_dir)
        quantifier = fit_and_quantify_spectrum_from_ledger(
            sample_ID=info.sample_id,
            spectrum_ID=int(payload["spectrum_id"]),
            results_path=str(Path(sample_dir).parent),
            ledger=info.ledger,
            **kw,
            **common,
        )
        extra = {"beam_energy_kV": info.ledger.configs.measurement_cfg.beam_energy_keV}
    if quantifier is None or getattr(quantifier, "fit_result", None) is None:
        raise RuntimeError("Fit failed. See the log for details.")
    min_cnts = ClusteringConfig().min_bckgrnd_cnts
    if payload["source"] != "file" and info.active_analysis is not None:
        min_cnts = info.active_analysis.config.min_bckgrnd_cnts
    flag = (_ledger_quant_flag(quantifier, kw.get("els_sample") or [], kw.get("els_substrate") or [], min_cnts)
            if quantify else None)
    return {**_fit_payload(quantifier, flag), **extra, "quantified": bool(quantify)}


def _run_acquisition(_: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Acquire (and optionally quantify) spectra of the samples (Acquisition tab)."""
    from autoemx.runners.batch_acquire_and_analyze import batch_acquire_and_analyze

    print(f"Acquiring {len(payload['samples'])} sample(s) into {payload['results_dir']}", flush=True)
    report = payload.get("report")  # run started by an external program: report each sample
    on_sample_done = None
    if report:
        from autoemx.gui import run_report

        def on_sample_done(sample_id, error):
            run_report.sample_done(report, sample_id, error)
            run_report.write_report(report)

    try:
        batch_acquire_and_analyze(
            samples=payload["samples"],
            results_dir=payload["results_dir"],
            development_mode=False,
            verbose=True,
            on_sample_done=on_sample_done,
            **payload["kwargs"],
        )
    except Exception as exc:
        if report:
            run_report.finish_report(report, "failed", f"{type(exc).__name__}: {exc}")
            run_report.write_report(report)
        raise
    if report:
        # Final here too, in case the GUI server is gone (the run goes on if its window was closed)
        run_report.finish_report(report, "finished")
        run_report.write_report(report)
    return {"samples": [s["ID"] for s in payload["samples"]]}


def inspect_spectra_folder(folder: str) -> Dict[str, Any]:
    """EMSA spectra of a folder and their acquisition settings, read from the file headers."""
    from autoemx.runners.quantify_external_spectra import read_energy_calibration
    from autoemx.utils import load_msa
    from autoemx.web.pipeline import parse_emsa_geometry

    path = Path(folder).expanduser()
    if not path.is_dir():
        raise FileNotFoundError(f"Folder not found: {folder}")
    files = sorted(p for p in path.iterdir()
                   if p.is_file() and p.suffix.lower() in cnst.EMSA_SPECTRUM_EXTENSIONS)
    beams = set()
    for f in files:
        try:
            _, _, metadata = load_msa(str(f))
            if any("BEAMKV" in k.upper() for k in metadata):
                beams.add(round(parse_emsa_geometry(metadata)["beam_energy"], 3))
        except Exception:
            pass
    calibration, calibration_error = None, None
    if files:
        try:
            calibration = read_energy_calibration(files)
        except Exception as exc:
            calibration_error = str(exc)
    return {
        "folder": str(path),
        "n_files": len(files),
        "beam_energies": sorted(beams),
        "calibration": calibration,
        "calibration_error": calibration_error,
    }


IMPORT_SAMPLE_TYPES = [cnst.S_POWDER_SAMPLE_TYPE, cnst.S_POWDER_CONTINUOUS_SAMPLE_TYPE,
                       cnst.S_BULK_SAMPLE_TYPE, cnst.S_BULK_ROUGH_SAMPLE_TYPE]


def import_kwargs(folder: str, results_dir: str, sample_id: str, elements: Any, substrate: Any,
                  sample_type: str, microscope_id: str, beam_energy: Any) -> Dict[str, Any]:
    """Validated kwargs of quantify_external_spectra to import a folder of spectra as a new sample."""
    sample_id = (sample_id or "").strip()
    if not sample_id or not re.fullmatch(r"[A-Za-z0-9_.+\-]+", sample_id):
        raise ValueError("The sample ID can only contain letters, digits and _ . + -")
    if not results_dir or not Path(results_dir).expanduser().is_dir():
        raise ValueError("Choose the results folder at the top of the page first.")
    results_dir = str(Path(results_dir).expanduser())
    if (Path(results_dir) / sample_id).exists():
        raise ValueError(f"A folder '{sample_id}' already exists in the results folder.")
    if not folder or not Path(folder).expanduser().is_dir():
        raise ValueError("Choose the folder of spectra.")
    els = _elements_list(elements)
    if not els:
        raise ValueError("Give the elements of the sample.")
    try:
        beam = float(beam_energy)
    except (TypeError, ValueError):
        raise ValueError("Give the beam energy (keV).") from None
    if not 1 <= beam <= 40:
        raise ValueError("The beam energy must be between 1 and 40 keV.")
    return {
        "samples": [{"ID": sample_id, "els": els, "spectra_dir": str(Path(folder).expanduser()),
                     "type": sample_type}],
        "microscope_ID": microscope_id,
        "beam_energy": beam,
        "els_substrate": _elements_list(substrate),
        "sample_substrate_type": "Ctape" if _elements_list(substrate) else "None",
        "sample_type": sample_type,
        "results_path": results_dir,
    }


def remove_unfinished_sample(sample_dir: Optional[str]) -> None:
    """Delete a sample folder left by an import that did not finish (one without a ledger)."""
    if sample_dir and os.path.isdir(sample_dir) and not os.path.exists(os.path.join(sample_dir, LEDGER_NAME)):
        import shutil

        shutil.rmtree(sample_dir, ignore_errors=True)


def _run_import(_: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Copy a folder of EMSA spectra into a new sample folder and build its ledger (Quantification tab)."""
    from autoemx.runners.quantify_external_spectra import quantify_external_spectra

    sample_dir = Path(payload["kwargs"]["results_path"]) / payload["kwargs"]["samples"][0]["ID"]
    if sample_dir.exists():
        raise RuntimeError(f"A folder '{sample_dir.name}' already exists in the results folder.")
    try:
        quantify_external_spectra(quantify=False, **payload["kwargs"])
    finally:
        remove_unfinished_sample(str(sample_dir))  # created by this import: no partial sample left
    if not (sample_dir / LEDGER_NAME).exists():
        raise RuntimeError("The sample could not be imported. See the log for details.")
    return {"sample_dir": str(sample_dir)}


QUANT_SAMPLE_MARKER = "@@ AutoEMX-GUI sample"


def _run_batch_quant(_: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    """Quantify several samples one after the other (Quantification tab)."""
    from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze

    samples = payload["samples"]
    results = []
    for i, sample in enumerate(samples):
        sample_dir = Path(sample["dir"])
        print(f"\n{QUANT_SAMPLE_MARKER} {i + 1}/{len(samples)}: {sample_dir.name}", flush=True)
        try:
            res = batch_quantify_and_analyze(
                sample_IDs=[sample_dir.name],
                results_path=str(sample_dir.parent),
                els_sample=sample.get("els_sample"),
                els_substrate=sample.get("els_substrate"),
                **payload["kwargs"],
            )
            ok = bool(res)
        except Exception:
            traceback.print_exc()
            ok = False
        results.append({"sample": sample_dir.name, "dir": str(sample_dir), "ok": ok})
        sys.stdout.flush()
    n_failed = sum(not r["ok"] for r in results)
    out: Dict[str, Any] = {"samples": results}
    if n_failed:
        out["warning"] = f"{n_failed} sample(s) failed. See the log."
    return out


def quant_progress(log_text: str, sample_names: List[str]) -> List[Dict[str, Any]]:
    """Progress of each sample of a batch quantification, parsed from its log."""
    progress = [{"sample": name, "state": "waiting", "done": 0, "total": None} for name in sample_names]
    blocks = re.split(re.escape(QUANT_SAMPLE_MARKER) + r" (\d+)/\d+: [^\n]*\n", log_text)
    # blocks = [preamble, idx1, text1, idx2, text2, ...]
    for idx_str, text in zip(blocks[1::2], blocks[2::2]):
        i = int(idx_str) - 1
        if i >= len(progress):
            continue
        total = re.findall(r"Starting quantification of (\d+) spectra", text)
        progress[i]["state"] = "running"
        progress[i]["total"] = int(total[-1]) if total else None
        if total:
            after = text[text.rfind("Starting quantification of"):]
            progress[i]["done"] = len(re.findall(r"^ Spectrum #\d+/\d+:", after, flags=re.M))
    started = [i for i, p in enumerate(progress) if p["state"] == "running"]
    for i in started[:-1]:
        progress[i]["state"] = "done"
    return progress


_JOB_FUNCS = {
    "analysis": _run_quant_and_analysis,
    "fit": _run_fit,
    "quant": _run_batch_quant,
    "single": _run_single,
    "acquisition": _run_acquisition,
    "import": _run_import,
}
# Jobs writing to sample ledgers: only one at a time
_LEDGER_WRITERS = {"analysis", "quant", "acquisition", "import"}



@dataclass
class Job:
    job_id: str
    kind: str
    sample_dir: str
    description: str
    process: Any
    log_path: str
    result_path: str
    started: float
    cancelled: bool = False

    def is_running(self) -> bool:
        return self.process.is_alive()

    def result(self) -> Optional[Dict[str, Any]]:
        if self.is_running() or not os.path.exists(self.result_path):
            if not self.is_running() and not os.path.exists(self.result_path):
                msg = "Cancelled." if self.cancelled else f"Process exited unexpectedly (code {self.process.exitcode})."
                return {"ok": False, "error": msg}
            return None
        with open(self.result_path, encoding="utf-8") as fh:
            return json.load(fh)

    def log_tail(self, max_chars: int = 20000) -> str:
        try:
            with open(self.log_path, "rb") as fh:
                fh.seek(0, os.SEEK_END)
                size = fh.tell()
                fh.seek(max(0, size - max_chars))
                text = fh.read().decode("utf-8", errors="replace")
        except FileNotFoundError:
            return ""
        return text if size <= max_chars else "…\n" + text.split("\n", 1)[-1]

    def elapsed(self) -> float:
        return time.time() - self.started


class JobManager:
    """Runs one job at a time per kind, each in its own process."""

    def __init__(self) -> None:
        self._jobs: Dict[str, Job] = {}
        self._lock = threading.Lock()
        self._tmpdir = tempfile.mkdtemp(prefix="autoemx_gui_")
        self._ctx = mp.get_context("spawn")

    def start(self, kind: str, sample_dir: str, payload: Dict[str, Any], description: str) -> Job:
        with self._lock:
            for job in self._jobs.values():
                if not job.is_running():
                    continue
                if job.kind == kind or (kind in _LEDGER_WRITERS and job.kind in _LEDGER_WRITERS):
                    raise RuntimeError(f"Wait for the running job to finish: {job.description}")
            job_id = uuid.uuid4().hex[:10]
            log_path = os.path.join(self._tmpdir, f"{job_id}.log")
            result_path = os.path.join(self._tmpdir, f"{job_id}.json")
            open(log_path, "w").close()
            proc = self._ctx.Process(
                target=_job_entry,
                args=(kind, sample_dir, payload, log_path, result_path),
                daemon=False,
            )
            proc.start()
            job = Job(job_id, kind, sample_dir, description, proc, log_path, result_path, time.time())
            self._jobs[job_id] = job
            return job

    def get(self, job_id: Optional[str]) -> Optional[Job]:
        return self._jobs.get(job_id) if job_id else None

    def cancel(self, job_id: str) -> None:
        job = self._jobs.get(job_id)
        if job is None or not job.is_running():
            return
        job.cancelled = True
        try:
            if hasattr(os, "killpg"):
                os.killpg(job.process.pid, signal.SIGTERM)
            else:
                # Windows: also stop the worker processes of the job (e.g. parallel quantification)
                subprocess.run(["taskkill", "/T", "/F", "/PID", str(job.process.pid)], capture_output=True,
                               timeout=30, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                job.process.terminate()
        except (OSError, subprocess.SubprocessError):  # ProcessLookupError: already ended
            pass
        job.process.join(timeout=5)
        if job.process.is_alive():
            job.process.kill()

    def shutdown(self) -> None:
        for job_id in list(self._jobs):
            self.cancel(job_id)


# =============================================================================
# Analysis results
# =============================================================================
@dataclass
class AnalysisData:
    """Everything needed to plot one clustering analysis."""

    sample_id: str
    folder: str
    quant_id: Optional[int]
    elements: List[str]
    detectable: List[str]
    features: str  # at_fr or w_fr
    comps: pd.DataFrame  # one row per spectrum, columns: spectrum, particle, cluster, status, <el> (fractions), ...
    centroids: np.ndarray  # (k, n_elements), fractions, in ``elements`` order
    stdevs: np.ndarray
    n_points: List[int]
    ref_formulae: List[str]
    ref_comps: pd.DataFrame  # index: formula, columns: elements (fractions)
    mixtures: List[List[Dict[str, Any]]]
    candidates: List[List[Tuple[str, float]]]  # per cluster: (candidate phase, confidence)
    clusters_table: Optional[pd.DataFrame]
    config: ClusteringConfig
    summary: Dict[str, Any]
    images: List[str]
    config_summary: str

    @property
    def unit(self) -> str:
        return "w%" if self.features == cnst.W_FR_CL_FEAT else "at%"


def reference_compositions(
    formulae: List[str], elements: List[str], features: str, undetectable: List[str]
) -> pd.DataFrame:
    """Candidate-phase compositions, computed as in the analyser (``_initialise_ref_phases``)."""
    from pymatgen.core import Composition

    from autoemx.utils import composition_as_weight_dict

    detectable = [el for el in elements if el not in undetectable]
    rows, names, seen = [], [], set()
    for formula in formulae:
        try:
            comp = Composition(formula)
        except Exception:
            continue
        if comp.reduced_formula in seen:
            continue
        seen.add(comp.reduced_formula)
        w_fr = composition_as_weight_dict(comp)
        if not any(el in w_fr for el in detectable):
            continue
        if features == cnst.W_FR_CL_FEAT:
            phase = {el: w_fr.get(el, 0.0) for el in detectable}
        else:
            det_w = {el: w for el, w in w_fr.items() if el in detectable}
            phase = comp.from_weight_dict(det_w).fractional_composition.as_dict()
        rows.append(phase)
        names.append(formula)
    df = pd.DataFrame(rows, columns=elements, index=names).fillna(0.0)
    return df.astype(float)


def _spectrum_key(spectrum_id: Any, index: int) -> str:
    sid = spectrum_id if spectrum_id is not None else index
    try:
        return str(int(float(sid)))
    except (TypeError, ValueError):
        return str(sid)


def _pixel(details: Any, axis: int) -> Optional[int]:
    coords = details.spot_coordinates if details is not None else None
    if coords is None or coords.pixel_coordinates is None:
        return None
    return int(coords.pixel_coordinates[axis])


def spectra_table(info: SampleInfo, quant_id: Optional[int], elements: List[str], features: str) -> pd.DataFrame:
    """One row per ledger spectrum with its composition for quantification ``quant_id``."""
    rows = []
    for idx, sp in enumerate(info.ledger.spectra):
        rec = next(
            (r for r in reversed(sp.quantification_results) if r.quantification_id == quant_id), None
        )
        details = sp.acquisition_details
        at = dict(rec.composition_atomic_fractions or {}) if rec else {}
        w = dict(rec.composition_weight_fractions or {}) if rec else {}
        fit = rec.fit_result if rec else None
        row: Dict[str, Any] = {
            "spectrum": _spectrum_key(sp.spectrum_id, idx),
            "particle": details.particle_id if details else None,
            "frame": details.frame_id if details else None,
            "total_counts": sp.total_counts,
            # Pixel position of the spot in the SEM image of its particle/frame
            "px": _pixel(details, 0),
            "py": _pixel(details, 1),
            "quant_flag": rec.quant_flag if rec else None,
            "an_err": rec.analytical_error * 100 if rec and rec.analytical_error is not None else None,
            "r_squared": fit.r_squared if fit else None,
            "redchi_sq": fit.reduced_chi_squared if fit else None,
            "comment": rec.comment if rec else None,
        }
        main = w if features == cnst.W_FR_CL_FEAT else at
        for el in elements:
            row[el] = main.get(el, 0.0) if main else np.nan
        for el in elements:
            row[f"{el} at%"] = at[el] * 100 if el in at else np.nan
        for el in elements:
            row[f"{el} w%"] = w[el] * 100 if el in w else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _cluster_labels(folder: str) -> Optional[Dict[str, Any]]:
    path = os.path.join(folder, "Compositions.csv")
    if not os.path.exists(path):
        return None
    raw = pd.read_csv(path)
    if cnst.SP_ID_DF_KEY not in raw or cnst.CL_ID_DF_KEY not in raw:
        return None
    return {
        _spectrum_key(sid, i): cl
        for i, (sid, cl) in enumerate(zip(raw[cnst.SP_ID_DF_KEY], raw[cnst.CL_ID_DF_KEY]))
    }


def load_analysis(info: SampleInfo, analysis: Optional[AnalysisRef]) -> AnalysisData:
    """Load one clustering analysis (or only the quantified compositions, if ``analysis`` is None)."""
    quant_id = analysis.quant_id if analysis is not None else info.active_quant
    quant = next((q for q in info.ledger.quantifications if q.quantification_id == quant_id), None)
    elements = list(dict.fromkeys(quant.sample_elements if quant and quant.sample_elements else info.elements))
    cfg = analysis.config if analysis is not None else ClusteringConfig()
    features = cfg.features
    undetectable = undetectable_elements(info.ledger)
    folder = analysis.folder if analysis is not None else ""

    comps = spectra_table(info, quant_id, elements, features)
    labels = _cluster_labels(folder) if folder else None
    comps["cluster"] = comps["spectrum"].map(labels) if labels is not None else np.nan
    quantified = comps[elements].notna().any(axis=1)
    if labels is None:
        comps["status"] = np.where(quantified, "not analysed", "not quantified")
    else:
        cl = comps["cluster"]
        comps["status"] = np.where(
            ~quantified, "not quantified",
            np.where(cl.isna(), "discarded", np.where(cl.fillna(-1) >= 0, "clustered", "noise")),
        )

    result = analysis.result if analysis is not None else None
    if result is not None:
        # Ledger centroids are stored in the clustering features (at_fr or w_fr), in ``elements`` order.
        cent_src, std_src = result.centroids, result.els_std_dev_per_cluster
        centroids = np.atleast_2d(np.asarray(cent_src, dtype=float))
        stdevs = np.atleast_2d(np.asarray(std_src, dtype=float))
        n_points = list(result.n_points_per_cluster)
        mixtures = [list(m or []) for m in result.clusters_assigned_mixtures]
        candidates = [_candidate_matches(row) for row in (result.refs_assigned_rows or [])]
        summary = {
            "n_clusters": len(result.centroids),
            "wcss": result.wcss,
            "silhouette": result.sil_score,
            "n_total": result.tot_n_points,
        }
    else:
        centroids = np.zeros((0, len(elements)))
        stdevs = np.zeros((0, len(elements)))
        n_points, mixtures, candidates = [], [], []
        summary = {"n_clusters": 0}
    if centroids.size == 0:
        centroids = np.zeros((0, len(elements)))
        stdevs = np.zeros((0, len(elements)))
    counts = comps["status"].value_counts().to_dict()
    summary.update({
        "n_spectra": len(comps),
        "n_clustered": int(counts.get("clustered", 0)),
        "n_noise": int(counts.get("noise", 0)),
        "n_discarded": int(counts.get("discarded", 0)),
        "n_not_quantified": int(counts.get("not quantified", 0)),
        "has_labels": labels is not None,
    })

    clusters_table = None
    clusters_path = os.path.join(folder, "Clusters.csv")
    if folder and os.path.exists(clusters_path):
        try:
            clusters_table = pd.read_csv(clusters_path)
        except Exception:
            clusters_table = None
    images = sorted(str(p) for p in Path(folder).glob("*.png")) if folder else []
    config_summary = ""
    summary_path = os.path.join(folder, "Analysis_config_summary.txt")
    if folder and os.path.exists(summary_path):
        config_summary = Path(summary_path).read_text(encoding="utf-8", errors="replace")

    return AnalysisData(
        sample_id=info.sample_id,
        folder=folder,
        elements=elements,
        detectable=[el for el in elements if el not in undetectable],
        features=features,
        comps=comps,
        centroids=centroids,
        stdevs=stdevs,
        n_points=n_points,
        ref_formulae=list(cfg.ref_formulae),
        quant_id=quant_id,
        ref_comps=reference_compositions(list(cfg.ref_formulae), elements, features, undetectable),
        mixtures=mixtures,
        candidates=candidates,
        clusters_table=clusters_table,
        config=cfg,
        summary=summary,
        images=images,
        config_summary=config_summary,
    )


def _candidate_matches(row: Dict[str, Any]) -> List[Tuple[str, float]]:
    """(candidate, confidence) pairs from a ``refs_assigned_rows`` entry ({'Cnd1': ..., 'CS_cnd1': ...})."""
    out = []
    j = 1
    while f"Cnd{j}" in row:
        name, conf = row.get(f"Cnd{j}"), row.get(f"CS_cnd{j}")
        if name:
            out.append((str(name), float(conf) if conf is not None else float("nan")))
        j += 1
    return out


def find_analysis(info: SampleInfo, key: Optional[str]) -> Optional[AnalysisRef]:
    if key:
        for a in info.analyses:
            if a.key == key:
                return a
    return info.active_analysis


def set_active_analysis(info: SampleInfo, key: Optional[str]) -> str:
    """Make the analysis selected in the GUI the active one of its sample; returns a message for the user."""
    from autoemx.runners.set_active_analysis import set_active_analysis as _set_active_analysis

    ref = find_analysis(info, key) if key else None
    if ref is None:
        raise ValueError("Select an analysis first.")
    if ref.is_active:
        return f"quant{ref.quant_id} / clust{ref.clust_id} is already the active analysis."
    _set_active_analysis(ref.quant_id, ref.clust_id, sample_path=info.sample_dir)
    return f"quant{ref.quant_id} / clust{ref.clust_id} is now the active analysis."


# =============================================================================
# Raw spectra
# =============================================================================
def load_raw_spectrum(info: SampleInfo, spectrum_id: Any) -> Dict[str, Any]:
    """Counts and energy axis of one spectrum, read from its file (no fit)."""
    ledger = info.ledger
    entry = None
    for idx, sp in enumerate(ledger.spectra):
        sid = sp.spectrum_id if sp.spectrum_id is not None else str(idx)
        try:
            match = int(float(sid)) == int(float(spectrum_id))
        except (TypeError, ValueError):
            match = str(sid) == str(spectrum_id)
        if match:
            entry = sp
            break
    if entry is None or not entry.spectrum_relpath:
        raise FileNotFoundError(f"Spectrum {spectrum_id} not found in the ledger.")
    counts = np.asarray(
        ledger._load_counts_from_pointer_file(Path(info.sample_dir) / entry.spectrum_relpath), dtype=float
    )
    mcfg = ledger.configs.microscope_cfg
    energy = mcfg.energy_zero + mcfg.bin_width * np.arange(len(counts))
    details = entry.acquisition_details
    return {
        "energy": energy,
        "counts": counts,
        "total_counts": entry.total_counts,
        "live_time": entry.live_acquisition_time,
        "particle": details.particle_id if details else None,
        "frame": details.frame_id if details else None,
        "path": str(Path(info.sample_dir) / entry.spectrum_relpath),
    }


# =============================================================================
# SEM images
# =============================================================================
SEM_IMAGE_EXTENSIONS = (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp")
# Whole frames, with the detected particles circled
_FRAME_RE = re.compile(r"_fr_?(?P<frame>.+?)_particles$")


@dataclass
class SemImage:
    """An SEM image saved during acquisition."""

    path: str
    name: str
    kind: str  # "spots" (spectrum spots of a particle or frame), "frame" (whole frame) or "other"
    particle: Optional[int] = None
    frame: Optional[str] = None

    @property
    def label(self) -> str:
        if self.kind == "spots" and self.particle is not None:
            return f"Particle {self.particle}" + (f" · frame {self.frame}" if self.frame not in (None, "None") else "")
        if self.kind in ("spots", "frame") and self.frame is not None:
            return f"Frame {self.frame}"
        return self.name


def find_sem_images(sample_dir: str) -> List[SemImage]:
    """SEM images of a sample: particle/frame images with spectrum spots first, then the others."""
    images: List[SemImage] = []
    folders = [Path(sample_dir) / cnst.IMAGES_DIR, Path(sample_dir)]
    for folder in folders:
        if not folder.is_dir():
            continue
        for path in sorted(folder.iterdir()):
            if path.suffix.lower() not in SEM_IMAGE_EXTENSIONS or path.stem.endswith("_raw"):
                continue
            parsed = parse_xsp_spots_image_name(path.stem)
            if parsed is not None:
                images.append(SemImage(str(path), path.name, "spots", parsed[0], parsed[1]))
            else:
                frame_match = _FRAME_RE.search(path.stem)
                if frame_match:
                    images.append(SemImage(str(path), path.name, "frame", frame=frame_match.group("frame")))
                else:
                    images.append(SemImage(str(path), path.name, "other"))

    def natural(text: Any) -> List[Any]:
        return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(text))]

    def order(im: SemImage):
        rank = {"spots": 0, "frame": 1, "other": 2}[im.kind]
        return (rank, im.particle if im.particle is not None else -1, natural(im.frame), natural(im.name))

    return sorted(images, key=order)


def sem_image_for_spectrum(images: List[SemImage], particle: Any, frame: Any) -> Optional[SemImage]:
    """Image showing the spot of a spectrum: the image of its particle, else of its frame."""
    spots = [im for im in images if im.kind == "spots"]
    try:
        par = int(float(particle))
    except (TypeError, ValueError):
        par = None
    if par is not None:
        matches = [im for im in spots if im.particle == par]
        if matches:
            exact = [im for im in matches if frame is not None and str(im.frame) == str(frame)]
            return (exact or matches)[0]
        return None
    if frame not in (None, "", "None"):
        matches = [im for im in spots if im.particle is None and str(im.frame) == str(frame)]
        if matches:
            return matches[0]
    return None


def spectra_on_image(comps: pd.DataFrame, image: SemImage) -> pd.DataFrame:
    """Rows of ``comps`` whose spots are in ``image`` (with pixel coordinates)."""
    if image.kind != "spots":
        return comps.iloc[0:0]
    rows = comps[comps["px"].notna() & comps["py"].notna()]
    if image.particle is not None:
        par = pd.to_numeric(rows["particle"], errors="coerce")
        return rows[par == image.particle]
    return rows[rows["frame"].astype(str) == str(image.frame)]


def read_image_png(path: str, max_size: Optional[int] = None) -> Tuple[bytes, Tuple[int, int]]:
    """PNG bytes of an image (first page of TIFFs), optionally downscaled, and its original size."""
    import io

    from PIL import Image

    with Image.open(path) as im:
        size = im.size
        im = im.convert("RGB")
        if max_size:
            im.thumbnail((max_size, max_size))
        buf = io.BytesIO()
        im.save(buf, format="PNG", optimize=False)
    return buf.getvalue(), size
