#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Load an uploaded EMSA spectrum and run AutoEMX fit/quantify."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Element

import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
from autoemx.runners.fit_and_quantify_spectrum import fit_and_quantify_spectrum
from autoemx.utils import load_msa
from autoemx.web.reader_report import SpectrumReadError, build_reader_failure_report

SUPPORTED_UPLOAD_EXTENSIONS = frozenset(cnst.EMSA_SPECTRUM_EXTENSIONS)

_DEFAULT_BEAM_ENERGY_KEV = 15.0
_DEFAULT_EMERGENCE_ANGLE_DEG = 28.5
_DEFAULT_LIVETIME_S = 10.0
# Shipped P/B standards are Phenom XL at 15 kV. Compositions are not valid otherwise.
QUANT_BEAM_KV = 15.0
_QUANT_BEAM_KV_TOLERANCE = 0.25


def beam_energy_supports_quantification(beam_energy_kV: Optional[float]) -> bool:
    """True if the spectrum beam energy matches the 15 kV quantification standards."""
    if beam_energy_kV is None:
        return False
    return abs(float(beam_energy_kV) - QUANT_BEAM_KV) <= _QUANT_BEAM_KV_TOLERANCE


@dataclass
class SpectrumFitResult:
    """Arrays and composition from one fitted/quantified spectrum."""

    filename: str
    energy_keV: np.ndarray
    counts: np.ndarray
    fit: np.ndarray
    background: np.ndarray
    composition_at: Dict[str, float]
    composition_wt: Dict[str, float]
    analytical_error: Optional[float] = None
    r_squared: Optional[float] = None
    reduced_chi_sq: Optional[float] = None
    quant_flag: Optional[int] = None
    peak_labels: List[Tuple[float, str]] = field(default_factory=list)
    els_sample: List[str] = field(default_factory=list)
    els_substrate: List[str] = field(default_factory=list)
    is_particle: bool = False
    error: Optional[str] = None
    error_stage: Optional[str] = None
    reader_report: Optional[str] = None
    beam_energy_kV: Optional[float] = None


def parse_elements(text: str) -> List[str]:
    """Parse a comma/space-separated element list and validate symbols."""
    if not text or not str(text).strip():
        return []
    tokens = [tok.strip() for tok in str(text).replace(";", ",").replace(" ", ",").split(",")]
    elements: List[str] = []
    seen = set()
    for token in tokens:
        if not token:
            continue
        symbol = token[0].upper() + token[1:].lower() if len(token) > 1 else token.upper()
        try:
            Element(symbol)
        except (ValueError, KeyError) as exc:
            raise ValueError(f"Unrecognized element symbol: '{token}'") from exc
        if symbol not in seen:
            elements.append(symbol)
            seen.add(symbol)
    return elements


def _metadata_value(metadata: Dict[str, str], *key_substrings: str) -> Optional[str]:
    for key, value in metadata.items():
        if any(sub in key.upper() for sub in key_substrings):
            return value
    return None


def parse_emsa_geometry(metadata: Dict[str, str]) -> Dict[str, float]:
    """Read detector/beam geometry from an EMSA header (keV units)."""
    units = str(metadata.get("XUNITS", "eV"))
    if "keV" in units:
        energy_scale = 1.0
    elif "eV" in units:
        energy_scale = 1000.0
    else:
        raise ValueError(
            f"Energy axis unit '{units}' unrecognized. Expected eV or keV."
        )

    offset = float(metadata.get("OFFSET", 0.0)) / energy_scale
    scale = float(metadata.get("XPERCHAN", 1.0)) / energy_scale

    beam_raw = _metadata_value(metadata, "BEAMKV")
    beam_energy = float(beam_raw) if beam_raw is not None else _DEFAULT_BEAM_ENERGY_KEV

    elev_raw = _metadata_value(metadata, "ELEVANGLE")
    emergence_angle = (
        float(elev_raw) if elev_raw is not None else _DEFAULT_EMERGENCE_ANGLE_DEG
    )

    live_raw = _metadata_value(metadata, "LIVETIME")
    livetime = float(live_raw) if live_raw is not None else _DEFAULT_LIVETIME_S

    return {
        "det_ch_offset": offset,
        "det_ch_width": scale,
        "beam_energy": beam_energy,
        "emergence_angle": emergence_angle,
        "sp_collection_time": livetime,
    }


def _plot_series(quantifier: Any) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    energy = np.asarray(quantifier.energy_vals, dtype=float)
    background = np.asarray(quantifier.background_vals, dtype=float)
    fitted = np.asarray(quantifier.fit_result.eval(), dtype=float)
    counts = np.asarray(quantifier.spectrum_vals, dtype=float)
    if not getattr(quantifier, "fit_background", True):
        counts = background + counts
        fitted = background + fitted
    return energy, counts, fitted, background


def _peak_labels(quantifier: Any) -> List[Tuple[float, str]]:
    labels: List[Tuple[float, str]] = []
    peaks = getattr(quantifier, "fitted_peaks_info", None) or {}
    ref_lines = set(getattr(quantifier, "ref_lines_for_quant", []) or [])
    for el_line, info in peaks.items():
        if ref_lines and el_line not in ref_lines:
            continue
        center = info.get(cnst.PEAK_CENTER_KEY)
        if center is None:
            continue
        labels.append((float(center), str(el_line)))
    return labels


def load_uploaded_spectrum(spectrum_path: str | Path) -> Tuple[Any, Dict[str, str], Dict[str, float]]:
    """Load counts and geometry from an EMSA file, or raise ``SpectrumReadError``."""
    path = Path(spectrum_path)
    if path.suffix.lower() not in SUPPORTED_UPLOAD_EXTENSIONS:
        exc = ValueError(
            f"Unsupported spectrum extension '{path.suffix}'. "
            f"Use one of: {', '.join(sorted(SUPPORTED_UPLOAD_EXTENSIONS))}"
        )
        raise SpectrumReadError(str(exc), build_reader_failure_report(path, exc)) from exc
    try:
        _, spectrum_vals, metadata = load_msa(str(path))
        if getattr(spectrum_vals, "size", len(spectrum_vals)) == 0:
            raise ValueError(
                "No count data was found after the EMSA header. "
                "The SPECTRUM marker or data columns may use an unsupported variant."
            )
        geometry = parse_emsa_geometry(metadata)
    except SpectrumReadError:
        raise
    except Exception as exc:
        raise SpectrumReadError(
            str(exc), build_reader_failure_report(path, exc)
        ) from exc
    return spectrum_vals, metadata, geometry


def fit_uploaded_spectrum(
    spectrum_path: str | Path,
    els_sample: Sequence[str],
    els_substrate: Sequence[str],
    is_particle: bool = True,
    spectrum_lims: Tuple[int, int] = dflt.spectrum_lims,
    fit_tol: float = 1e-3,
    microscope_ID: str = dflt.microscope_ID,
    meas_type: str = dflt.measurement_type,
    meas_mode: str = dflt.measurement_mode,
) -> SpectrumFitResult:
    """Fit and quantify one EMSA file. Does not open a matplotlib window."""
    path = Path(spectrum_path)
    if not els_sample:
        raise ValueError("Sample elements are required (the fitter does not auto-identify peaks).")

    spectrum_vals, _metadata, geometry = load_uploaded_spectrum(path)

    quantifier = fit_and_quantify_spectrum(
        spectrum_vals=spectrum_vals,
        microscope_ID=microscope_ID,
        meas_type=meas_type,
        meas_mode=meas_mode,
        det_ch_offset=geometry["det_ch_offset"],
        det_ch_width=geometry["det_ch_width"],
        beam_energy=geometry["beam_energy"],
        emergence_angle=geometry["emergence_angle"],
        sp_collection_time=geometry["sp_collection_time"],
        sample_ID=path.stem,
        background_vals=None,
        els_sample=list(els_sample),
        els_substrate=list(els_substrate),
        spectrum_lims=spectrum_lims,
        quantify_plot=True,
        plot_signal=False,
        fit_tol=fit_tol,
        is_particle=is_particle,
        print_results=False,
        quant_verbose=False,
        fitting_verbose=False,
    )
    if quantifier is None or getattr(quantifier, "fit_result", None) is None:
        raise RuntimeError(f"Fitting failed for '{path.name}'.")

    energy, counts, fitted, background = _plot_series(quantifier)
    quant_result = getattr(quantifier, "quant_result", None) or {}
    composition_at = dict(quant_result.get(cnst.COMP_AT_FR_KEY) or {})
    composition_wt = dict(quant_result.get(cnst.COMP_W_FR_KEY) or {})

    return SpectrumFitResult(
        filename=path.name,
        energy_keV=energy,
        counts=counts,
        fit=fitted,
        background=background,
        composition_at=composition_at,
        composition_wt=composition_wt,
        analytical_error=quant_result.get(cnst.AN_ER_KEY),
        r_squared=quant_result.get(cnst.R_SQ_KEY),
        reduced_chi_sq=quant_result.get(cnst.REDCHI_SQ_KEY),
        quant_flag=getattr(quantifier, "bad_quant_flag", None),
        peak_labels=_peak_labels(quantifier),
        els_sample=list(els_sample),
        els_substrate=list(els_substrate),
        is_particle=is_particle,
        beam_energy_kV=float(geometry["beam_energy"]),
    )
