#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automatic element identification
================================

Identifies the elements present in an X-ray spectrum, starting from a (possibly empty) list of
elements known to be present.

Candidates are generated from the X-ray line table and verified through full-spectrum fitting,
starting from the most intense unexplained peak:

1. Fit the spectrum with the elements already in the model (fixed, substrate and identified elements).
2. Search the fit residual for significant peaks with a top-hat filter, which suppresses any locally
   linear mismatch of the background (Statham, Anal. Chem. 49 (1977) 2149).
3. For the most intense unexplained peak, list the elements with a strong line within the detector
   resolution, plus their known overlap partners (common_peak_overlaps.json).
4. Screen the candidates with a filtered weighted least-squares fit of the templates of their line
   families to the residual. Lines of a family have fixed relative weights, so a candidate whose other
   lines are not observed is penalised. If a higher-energy family of the candidate is well excited,
   it must be observed too (confirmation lines; Newbury, Microsc. Microanal. 11 (2005) 545).
5. Verify the best candidates with a full refit. A candidate is accepted if it decreases the chi-square
   of the fit (Poisson variance plus a relative model uncertainty) by more than a threshold, and if the fitted area of its line exceeds the Currie
   detection limit L_D = 2.71 + 3.29 sqrt(B) (Currie, Anal. Chem. 40 (1968) 586).
6. Repeat until no significant peak is left; then remove identified elements that the final fit does
   not need (backward pruning).

This is the iterative qualitative-quantitative strategy of Newbury & Ritchie, Microsc. Microanal.
24 (2018) 350. Statistical tests on nested fits are not strictly valid when the null hypothesis is a
peak area of zero (Protassov et al., ApJ 571 (2002) 545), so all thresholds are empirical.

The fitting is delegated to a callable (see `make_xsp_fit_function`), so the identification loop is
independent of the fitter.
"""
from __future__ import annotations

import time
import warnings
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Element
from scipy.optimize import nnls
from scipy.signal import find_peaks

import autoemx.calibrations as calibs
from autoemx._logging import get_logger
from autoemx.core.fitter.detector_response import DetectorResponseFunction
from autoemx.core.fitter.fitter import get_reference_xray_line
from autoemx.core.fitter.peaks import OVERLAP_SEPARATION_SIGMAS
from autoemx.data.Xray_lines import get_el_xray_lines

logger = get_logger(__name__)

# Lines are fitted only if the overvoltage exceeds this value (see XSp_Fitter._define_xray_lines)
MIN_FIT_OVERVOLTAGE = 1.2
# Energy of the Si Ka line (keV), for escape peaks of Si detectors
SI_KA_ENERGY = 1.740
MIN_Z = 5  # B
MAX_Z = 94  # Pu
FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))
# Overlaps of lines below this energy (keV) are not reported (lines strongly absorbed, rarely visible)
MIN_OVERLAP_NOTE_ENERGY = 0.3


# =============================================================================
# Settings and results
# =============================================================================
@dataclass
class PeakIDSettings:
    """
    Thresholds of the element identification.

    Attributes
    ----------
    min_peak_significance : float
        Minimum significance (in standard deviations) of a residual peak, and of the candidate line family
        that explains it.
    low_energy_threshold_keV, low_energy_multiplier : float
        Thresholds are multiplied by `low_energy_multiplier` for lines below `low_energy_threshold_keV`, where
        absorption edges and chemical effects leave structure in the residual.
    element_multipliers : dict
        Threshold multipliers for specific elements (none by default).
    artifact_min_fraction : dict
        Elements also produced by the detector (Si: internal fluorescence of the Si detector, ~0.05% of the counts):
        they are kept only if the fitted area of their main line exceeds this fraction of the total counts.
    candidate_window_sigma : float
        Lines within this many detector sigmas of a residual peak are considered as candidates.
    merged_peak_shift_keV : float
        Extra width of the candidate window, accounting for the shift of the maximum of unresolved
        La+Lb and Ma+Mb peaks from the tabulated energies.
    min_line_weight : float
        Minimum weight of a line (relative to the reference line of its family) to be used to match a peak.
    confirm_overvoltage, confirm_significance : float
        A higher-energy family with overvoltage above `confirm_overvoltage` (default: the overvoltage used by
        XSp_Quantifier to select ideal reference lines, ideal_ref_line_overvoltage) must be observed with at least
        `confirm_significance` standard deviations for a candidate to be accepted.
    model_rel_uncertainty : float
        Relative uncertainty of the spectrum model (peak shapes, line weights, background), added to the
        Poisson variance as (f·model)². Systematic model errors grow with the counts while Poisson noise grows
        with their square root, so without this term every small model imperfection becomes significant in
        high-count (e.g. summed) spectra.
    min_peak_gain : float
        Local fit improvement (decrease of the chi-square within a resolved peak of a candidate) for the peak to
        support the candidate; a decrease of the same amount with the opposite sign makes the peak contradict it.
        A candidate is rejected if any of its resolved peaks predicted above the detection limit contradicts it
        (a peak predicted where there is none), whatever the total fit improvement.
    min_supported_peaks_overlap_check : int
        Minimum number of resolved peaks supporting a candidate found by the overlap check (an element absorbed
        by overlapping peaks must improve the fit at several of its peaks to be told apart from model mismatch).
    min_delta_deviance : float
        Minimum decrease of the fit statistic Σ (y − μ)² / (μ + (f·μ)²) for a candidate to be accepted. In the
        full-fit pruning, an element is also kept if its filtered (peak) gain reaches this value while the total
        statistic still decreases.
    ambiguity_margin : float
        Candidates whose decrease of the fit statistic is within this margin of the best one are reported
        as ambiguous.
    max_full_fits_per_peak : int
        Maximum number of screened candidates verified with a full fit, when full fits are needed (overlaps
        with elements in the model). Otherwise all screened candidates are tested with fast linear fits.
    max_new_elements : int
        Maximum number of identified elements.
    max_iterations : int
        Maximum number of residual peaks examined.
    z_range : (int, int)
        Atomic numbers of the elements proposed as candidates (default B to Bi: all elements above Bi are
        radioactive and essentially absent from samples; use (5, 92) to include Th and U).
    excluded_elements : tuple of str
        Elements never proposed as candidates (default: Tc and Pm, without stable isotopes, and the noble gases
        Ne, Kr, Xe, Rn, absent from solids; Ar is allowed, as it can be implanted by ion milling or sputtering).
    strong_min_filtered_gain : float
        Identified elements are classified as 'strong' (to be added to the fit and quantified) if (1) in the full
        fit without them, at least one of their resolved peaks is a clear positive peak of the filtered residual
        (min_peak_significance, with the threshold multipliers), (2) the chi-square of the top-hat filtered residual
        (insensitive to smooth background changes) improves by at least `strong_min_filtered_gain` when they are
        added, (3) the background is reliable, and (4) no other candidate explains the peak equally well.
        Otherwise they are 'possible' (only reported, for manual checking): e.g. heavy elements whose many weak
        lines fill a background or absorption mismatch, without a peak of their own in the residual.
    background_ratio_bounds, background_bands_keV
        Background health: ratio of data to model in each energy band, in channels without significant peaks
        (excluding fitted peaks and unexplained peaks of the residual). Bands are checked separately, so that a
        wrong background shape (deficit in one band, excess in another) is not averaged out. If a band is outside
        the bounds in the starting fit, the fit is repeated with the background scale K fixed (see
        `background_refit_function` of identify_elements), and the refit is kept if it reduces the worst band
        deviation without creating a background deficit that the starting fit did not have (a model above the
        data is harmless, e.g. when a missing major element changes the background composition). If the data still exceed the model by more than the upper bound in a band (background deficit,
        which any element could fill), also after adding the identified elements (a missing major element changes
        the background too), all identified elements are 'possible'.
    max_start_redchi_fraction : float or None
        If set, the identification is skipped if the reduced chi-square of the starting fit exceeds this fraction
        of the total counts (0.0016 is the criterion of XSp_Quantifier for unreliable fits, quant flag 1). Off by
        default: a missing element also makes the starting fit poor.
    linear_prefilter_fraction : float
        For candidates overlapping peaks of the current fit, the linear test underestimates the fit improvement
        (frozen shapes of the overlapping peaks have partly absorbed the candidate), so it is only used as a
        pre-filter at this fraction of min_delta_deviance; the full fit then applies the threshold.
    contradiction_rel_loss : float
        A predicted peak contradicts a candidate only if the local fit worsens by more than this fraction of the
        candidate's total fit improvement (and by more than min_peak_gain), so that small mismatches of weak lines
        (line weights) do not reject an element with strong evidence elsewhere.
    min_local_gain_fraction : float
        Minimum fraction of the fit improvement brought by a candidate that must come from peak regions (where
        the candidate or other elements have significant peaks), rather than from the background. Rejects elements
        that improve the fit by changing the background (e.g. a heavy element rescuing a fit stuck with too little
        background), rather than by explaining peaks.
    prune : bool
        Whether to remove identified elements that the final fit does not need. Uses full fits (one with all
        identified elements, one without each of them), run only if elements were identified.
    overlap_check : bool
        Whether to test the overlap partners (overlap table) of the elements in the model, even without residual
        peaks: an element whose peaks are absorbed by overlapping peaks of other elements leaves no residual
        peak, but still improves the fit. Partners are tested with linear fits (tied shifts); the best one
        passing the thresholds is confirmed with a full fit, and the check is repeated until none is confirmed.
    full_fit_overlaps : bool
        Whether candidates overlapping peaks of the current fit are tested with full fits instead of linear
        fits (in both cases, with energy shifts tied to the overlapping peaks).
    refit_after_accept : bool
        Whether to confirm each element accepted with a linear test by a full fit, so that the next residual
        peaks are searched with the fitter's peak shapes. Needed for high-count spectra, where the mismatch
        of the linear templates becomes significant; off by default for single spectra.
    final_fit : bool
        Whether the final model must come from a full fit. Off by default, since the caller refits anyway
        with the identified elements (e.g. quantification after the identification of missing elements).
    """
    min_peak_significance: float = 5.0
    low_energy_threshold_keV: float = 1.0
    low_energy_multiplier: float = 1.5
    element_multipliers: Dict[str, float] = field(default_factory=dict)
    artifact_min_fraction: Dict[str, float] = field(default_factory=lambda: {'Si': 0.002})
    candidate_window_sigma: float = 2.5
    merged_peak_shift_keV: float = 0.025
    min_line_weight: float = 0.1
    confirm_overvoltage: Optional[float] = None
    confirm_significance: float = 3.0
    model_rel_uncertainty: float = 0.01
    min_peak_gain: float = 9.0
    min_supported_peaks_overlap_check: int = 2
    min_delta_deviance: float = 100.0
    ambiguity_margin: float = 30.0
    max_full_fits_per_peak: int = 4
    max_new_elements: int = 10
    max_iterations: int = 30
    z_range: Tuple[int, int] = (5, 83)
    excluded_elements: Tuple[str, ...] = ('Tc', 'Pm', 'Ne', 'Kr', 'Xe', 'Rn')
    prune: bool = True
    max_start_redchi_fraction: Optional[float] = None
    strong_min_filtered_gain: float = 100.0
    background_ratio_bounds: Tuple[float, float] = (0.8, 1.25)
    background_bands_keV: Tuple[Tuple[float, float], ...] = ((1.0, 2.5), (2.5, 5.0), (5.0, 100.0))
    min_local_gain_fraction: float = 0.2
    contradiction_rel_loss: float = 0.05
    linear_prefilter_fraction: float = 0.25
    full_fit_overlaps: bool = False
    overlap_check: bool = True
    refit_after_accept: bool = False
    final_fit: bool = False

    def multiplier(self, el: str, energy_keV: float) -> float:
        """Threshold multiplier for a line of an element."""
        mult = self.element_multipliers.get(el, 1.0)
        if energy_keV < self.low_energy_threshold_keV:
            mult *= self.low_energy_multiplier
        return mult


@dataclass
class FitOutcome:
    """Result of a fit with a given set of elements, as returned by the fit function."""
    elements: Tuple[str, ...]
    model: np.ndarray
    deviance: float
    line_areas: Dict[str, float]  # Fitted area (counts) of each characteristic line ('Fe_Ka1': area)
    # Peaks of each line family, including escape and pile-up peaks ('Fe_Ka1': array), and the rest of
    # the model (background), so that candidates can be tested with a linear fit of fixed-shape components
    components: Dict[str, np.ndarray] = field(default_factory=dict)
    background: Optional[np.ndarray] = None
    approximate: bool = False  # True for linear refits of fixed-shape components (see _linear_refit)
    line_centers: Dict[str, float] = field(default_factory=dict)  # Fitted energy (keV) of each reference line
    redchi: Optional[float] = None  # Reduced chi-square of the (unweighted) fit, as used by the quantifier


@dataclass
class IdentifiedElement:
    """An element of the final model, with the evidence for its presence."""
    element: str
    status: str  # 'fixed', 'substrate', 'identified', 'ambiguous', 'pruned'
    matched_line: Optional[str] = None
    line_energy_keV: Optional[float] = None
    peak_significance: Optional[float] = None
    delta_deviance: Optional[float] = None
    net_area: Optional[float] = None
    detection_limit: Optional[float] = None
    ambiguous_with: List[str] = field(default_factory=list)
    overlaps_with: List[str] = field(default_factory=list)
    has_standard: Optional[bool] = None
    from_pool: bool = False  # Identified among the priority elements, without a full candidate search
    via_overlap_check: bool = False  # Found by testing the overlap partners of elements in the model
    from_peak_assignment: bool = False  # Assigned to peaks of the raw spectrum, before the starting fit
    filtered_gain: Optional[float] = None  # Full-fit improvement of the top-hat filtered residual (peaks only)
    residual_peak_sigma: Optional[float] = None  # Largest positive filtered-residual peak at its lines, fit without it
    residual_peak_ok: Optional[bool] = None  # Whether that peak is clear (see _residual_peak_check)
    note: str = ''  # Why the element was removed or only reported
    tier: str = ''  # 'strong' (add and quantify) or 'possible' (report only), for identified elements

    def describe(self) -> str:
        msg = self.element + (f" [{self.tier}]" if self.tier else '')
        if self.matched_line:
            msg += f" ({self.matched_line} {self.line_energy_keV:.3f} keV"
            if self.peak_significance is not None:
                msg += f", {self.peak_significance:.0f}σ"
            if self.delta_deviance is not None:
                msg += f", ΔD = {self.delta_deviance:.0f}"
            if self.filtered_gain is not None:
                msg += f", peak ΔD = {self.filtered_gain:.0f}"
            if self.residual_peak_sigma is not None:
                msg += f", residual peak {self.residual_peak_sigma:.0f}σ"
            msg += ")"
        if self.from_pool:
            msg += " [pool]"
        if self.via_overlap_check:
            msg += " [overlap check]"
        if self.from_peak_assignment:
            msg += " [peak assignment]"
        if self.ambiguous_with:
            msg += f" — ambiguous with {', '.join(self.ambiguous_with)}"
        if self.overlaps_with:
            msg += f" — overlaps {', '.join(self.overlaps_with)}"
        if self.has_standard is False:
            msg += " — no standard, cannot be quantified"
        if self.note:
            msg += f" — {self.note}"
        return msg


@dataclass
class UnexplainedPeak:
    """A significant residual peak that no candidate element could explain."""
    energy_keV: float
    significance: float
    net_area: float
    note: str = ''


@dataclass
class ElementIDResult:
    """Result of the element identification."""
    fixed_elements: List[str]
    substrate_elements: List[str]
    elements: List[IdentifiedElement]
    unexplained_peaks: List[UnexplainedPeak]
    energy_vals: np.ndarray
    spectrum_vals: np.ndarray
    model: Optional[np.ndarray]
    n_fits: int
    elapsed_s: float
    n_candidate_tests: int = 0
    skipped_reason: str = ''
    background_note: str = ''

    @property
    def new_elements(self) -> List[str]:
        """Identified elements (including ambiguous identifications), in order of identification."""
        return [e.element for e in self.elements if e.status in ('identified', 'ambiguous')]

    @property
    def strong_new_elements(self) -> List[str]:
        """Identified elements with strong evidence, to be added to the fit."""
        return [e.element for e in self.elements if e.status in ('identified', 'ambiguous') and e.tier == 'strong']

    @property
    def possible_new_elements(self) -> List[str]:
        """Identified elements with weak evidence, only reported as possibly present."""
        return [e.element for e in self.elements if e.status in ('identified', 'ambiguous') and e.tier == 'possible']

    @property
    def quantifiable_new_elements(self) -> List[str]:
        """Strong elements that can be quantified (unknown standards are assumed available)."""
        return [e.element for e in self.elements
                if e.status in ('identified', 'ambiguous') and e.tier == 'strong' and e.has_standard is not False]

    @property
    def unquantifiable_new_elements(self) -> List[str]:
        """Strong elements without a standard: fitted and reported as present, but not quantified."""
        return [e.element for e in self.elements
                if e.status in ('identified', 'ambiguous') and e.tier == 'strong' and e.has_standard is False]

    def summary(self) -> str:
        lines = []
        if self.skipped_reason:
            return f"Element identification not performed: {self.skipped_reason}"
        fixed = [e for e in self.elements if e.status in ('fixed', 'substrate')]
        new = [e for e in self.elements if e.status in ('identified', 'ambiguous')]
        pruned = [e for e in self.elements if e.status == 'pruned']
        lines.append(f"Elements in the starting model: {', '.join(e.element for e in fixed) or 'none'}")
        strong = [e for e in new if e.tier == 'strong']
        possible = [e for e in new if e.tier != 'strong']
        if strong:
            lines.append("Identified elements (added to the fit):")
            lines.extend(f"   • {e.describe()}" for e in strong)
        if possible:
            lines.append("Possibly present (weak evidence, not added; check manually):")
            lines.extend(f"   • {e.describe()}" for e in possible)
        if not new:
            lines.append("No new elements identified.")
        if pruned:
            lines.append(f"Removed in final check: {', '.join(e.element for e in pruned)}")
        if self.background_note:
            lines.append(f"Background: {self.background_note}")
        if self.unexplained_peaks:
            lines.append("Unexplained peaks:")
            for p in self.unexplained_peaks:
                note = f" ({p.note})" if p.note else ''
                lines.append(f"   • {p.energy_keV:.3f} keV, {p.significance:.0f}σ, ~{p.net_area:.0f} counts{note}")
        lines.append(f"{self.n_fits} full fits, {self.n_candidate_tests} candidate tests, {self.elapsed_s:.1f} s")
        return "\n".join(lines)

    def plot(self, ax=None, title: str = ''):
        """Plot spectrum, final model, residual and identified elements."""
        import matplotlib.pyplot as plt
        if ax is None:
            fig, (ax, ax_res) = plt.subplots(2, 1, sharex=True, figsize=(10, 6), height_ratios=(3, 1))
        else:
            fig, ax_res = ax.figure, None
        E, y = self.energy_vals, self.spectrum_vals
        ax.plot(E, y, color='0.4', lw=0.8, label='Spectrum')
        if self.model is not None:
            ax.plot(E, self.model, color='C0', lw=1, label='Fit')
            if ax_res is not None:
                ax_res.plot(E, (y - self.model) / np.sqrt(np.maximum(self.model, 1)), color='0.3', lw=0.6)
                ax_res.axhline(0, color='C0', lw=0.6)
                ax_res.set_ylabel('Residual (σ)')
                ax_res.set_xlabel('Energy (keV)')
        for e in self.elements:
            if e.line_energy_keV is None or e.status == 'pruned':
                continue
            color = 'C3' if e.status in ('identified', 'ambiguous') else 'C2'
            ax.axvline(e.line_energy_keV, color=color, lw=0.6, ls='--')
            ax.annotate(e.element, (e.line_energy_keV, ax.get_ylim()[1] * 0.92), color=color, ha='center')
        for p in self.unexplained_peaks:
            ax.axvline(p.energy_keV, color='C1', lw=0.6, ls=':')
        ax.set_ylabel('Counts')
        ax.set_title(title)
        ax.legend()
        return fig


# =============================================================================
# Detector resolution and line table
# =============================================================================
@lru_cache(maxsize=None)
def _detector_params(microscope_ID: str) -> Tuple[float, float, float]:
    """Detector resolution parameters (conv_eff, elec_noise, F) from the detector model of a microscope."""
    model = calibs.load_detector_model(microscope_ID)
    return model.conv_eff, model.elec_noise, model.F


def get_sigma_function(microscope_ID: str) -> Callable:
    """Detector Gaussian sigma (keV) as a function of energy (keV), for a microscope."""
    params = _detector_params(microscope_ID)
    return lambda E: DetectorResponseFunction.det_sigma(np.asarray(E, dtype=float), *params)


@lru_cache(maxsize=None)
def _xray_lines(el: str) -> Dict[str, dict]:
    """Cached get_el_xray_lines (each call builds a pandas table)."""
    return get_el_xray_lines(el)


def _family_reference_line(el: str, line: str, el_lines: Dict[str, dict]) -> str:
    """Reference line of the family of a line, as in XSp_Fitter (see get_reference_xray_line)."""
    # Imported here: the quantifier imports this module
    from autoemx.core.quantifier.quantifier import XSp_Quantifier
    el_ref_line = get_reference_xray_line(el, line, el_lines, XSp_Quantifier.xray_quant_ref_lines)
    return el_ref_line.split('_', 1)[1]


@lru_cache(maxsize=None)
def _element_families(el: str, beam_energy: float, e_min: float, e_max: float) -> Dict[str, Tuple[Tuple[str, float, float], ...]]:
    """
    Fitted line families of an element: {reference line: ((line, energy, weight), ...)}.

    Lines are included if excited by the beam and within the energy range, as in XSp_Fitter.
    """
    try:
        el_lines = _xray_lines(el)
    except (KeyError, ValueError):
        return {}
    families: Dict[str, List[Tuple[str, float, float]]] = {}
    for line, info in el_lines.items():
        en = float(info['energy (keV)'])
        if beam_energy / en <= MIN_FIT_OVERVOLTAGE or not e_min <= en <= e_max:
            continue
        ref = _family_reference_line(el, line, el_lines)
        families.setdefault(ref, []).append((line, en, float(info['weight'])))
    # Keep families whose reference line is fitted
    return {
        ref: tuple(lines) for ref, lines in families.items()
        if any(line == ref for line, _, _ in lines)
    }


@lru_cache(maxsize=None)
def _strong_lines_table(beam_energy: float, e_min: float, e_max: float, min_weight: float) -> Tuple[Tuple[float, str, str, str], ...]:
    """All strong fitted lines of elements B to Pu: ((energy, element, line, family reference line), ...)."""
    table = []
    for Z in range(MIN_Z, MAX_Z + 1):
        el = Element.from_Z(Z).symbol
        for ref, lines in _element_families(el, beam_energy, e_min, e_max).items():
            for line, en, w in lines:
                if w >= min_weight:
                    table.append((en, el, line, ref))
    table.sort()
    return tuple(table)


def element_line_families(el: str, beam_energy: float, energy_range: Tuple[float, float]) -> Dict[str, Tuple[Tuple[str, float, float], ...]]:
    """Public accessor of the fitted line families of an element (see _element_families)."""
    return _element_families(el, float(beam_energy), float(energy_range[0]), float(energy_range[1]))


# =============================================================================
# Overlap table
# =============================================================================
def _overlap_partners(
    el: str,
    ref_line: str,
    sigma_fn: Callable,
    meas_type: str = 'EDS',
) -> List[Tuple[str, str, str]]:
    """
    Known overlap partners of a line, from common_peak_overlaps.json: [(element, line, severity), ...].

    Severity is computed from the line separation in units of the RMS detector sigma, with the thresholds
    defined in the table.
    """
    # Imported here: peak_overlaps imports the quantifier, which imports the fitter package
    from autoemx.core.quantifier.peak_overlaps import _load_method_table
    table = _load_method_table(meas_type)
    if table is None:
        return []
    thresholds = table.get('severity_thresholds_sigma', {})
    partners = []
    for pair in table['pairs']:
        lines = pair['lines']
        for this, other in ((lines[0], lines[1]), (lines[1], lines[0])):
            if this['element'] == el and this['line'] == ref_line and other['element'] != el:
                e_a, e_b = this['energy_keV'], other['energy_keV']
                sep = abs(e_a - e_b) / np.sqrt((sigma_fn(e_a) ** 2 + sigma_fn(e_b) ** 2) / 2)
                severity = next((label for label, thr in thresholds.items() if sep < thr), None)
                if severity is not None:
                    partners.append((other['element'], other['line'], severity))
    return partners


# =============================================================================
# Filtering, templates, residual peaks
# =============================================================================
def top_hat_matrix(energy_vals: np.ndarray, sigma_fn: Callable) -> np.ndarray:
    """
    Top-hat filter as a matrix (filtered = F @ spectrum), with width following the detector resolution.

    Each row has a positive central lobe about one FWHM wide and two negative side lobes of half that
    width, summing to zero, so that any background that is linear across the filter is removed.
    Rows whose filter does not fit within the spectrum are zero.

    References: F. H. Schamber, in X-Ray Fluorescence Analysis of Environmental Samples, T. G. Dzubay (ed.),
    Ann Arbor Science (1977) 241; P. J. Statham, Anal. Chem. 49 (1977) 2149.
    """
    n = len(energy_vals)
    ch_width = energy_vals[1] - energy_vals[0]
    F = np.zeros((n, n))
    fwhm_ch = FWHM_PER_SIGMA * sigma_fn(energy_vals) / ch_width
    for i in range(n):
        m = max(1, int(round(fwhm_ch[i] / 2)))  # Central lobe: 2m + 1 channels
        v = m  # Each side lobe: v channels
        lo, hi = i - m - v, i + m + v
        if lo < 0 or hi >= n:
            continue
        F[i, i - m:i + m + 1] = 1 / (2 * m + 1)
        F[i, lo:i - m] = -0.5 / v
        F[i, i + m + 1:hi + 1] = -0.5 / v
    return F


def gaussian_template(energy_vals: np.ndarray, lines: Iterable[Tuple[float, float]], sigma_fn: Callable) -> np.ndarray:
    """Sum of Gaussian peaks of area `weight` (counts) at the given (energy, weight) pairs."""
    ch_width = energy_vals[1] - energy_vals[0]
    out = np.zeros_like(energy_vals, dtype=float)
    for en, w in lines:
        s = float(sigma_fn(en))
        out += w * ch_width / (s * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((energy_vals - en) / s) ** 2)
    return out


def fit_statistic(y: np.ndarray, mu: np.ndarray, rel_uncertainty: float) -> float:
    """Chi-square of a fit, with Poisson variance plus a relative model uncertainty: Σ (y − μ)² / (μ + (f·μ)²)."""
    mu = np.asarray(mu, dtype=float)
    var = np.maximum(mu, 1.0) + (rel_uncertainty * mu) ** 2
    return float(np.sum((np.asarray(y, dtype=float) - mu) ** 2 / var))


def poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    """Poisson deviance 2 Σ [y ln(y/μ) − (y − μ)]."""
    y = np.asarray(y, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1e-9, None)
    with np.errstate(divide='ignore', invalid='ignore'):
        term = np.where(y > 0, y * np.log(y / mu), 0.0)
    return float(2 * np.sum(term - (y - mu)))


@dataclass
class ResidualPeak:
    channel: int
    energy_keV: float
    significance: float
    net_area: float


class _ResidualAnalysis:
    """Filtered residual of a fit, with its covariance, for peak search and candidate screening."""

    def __init__(self, y: np.ndarray, model: Optional[np.ndarray], energy_vals: np.ndarray,
                 F: np.ndarray, sigma_fn: Callable, rel_uncertainty: float = 0.0):
        self.energy_vals = energy_vals
        self.F = F
        self.sigma_fn = sigma_fn
        if model is None:
            residual = y.astype(float)
            var = np.maximum(y, 1.0) + (rel_uncertainty * y) ** 2
        else:
            residual = y - model
            var = np.maximum(0.5 * (model + y), 1.0) + (rel_uncertainty * model) ** 2
        self.var = var
        self.r_f = F @ residual
        self.var_f = (F ** 2) @ var
        self.valid = self.var_f > 0
        self.sig = np.zeros_like(self.r_f)
        self.sig[self.valid] = self.r_f[self.valid] / np.sqrt(self.var_f[self.valid])

    def sandwich(self, B: np.ndarray) -> np.ndarray:
        """
        Bᵀ Cov B, with Cov = F diag(var) Fᵀ the covariance of the filtered residual (filtered channels are
        correlated), computed without forming the n x n covariance matrix.
        """
        G = self.F.T @ B
        return G.T @ (G * self.var[:, None])

    def find_peaks(self, settings: PeakIDSettings) -> List[ResidualPeak]:
        """Significant positive residual peaks, sorted by decreasing net area."""
        thr = np.array([settings.min_peak_significance * settings.multiplier('', en) for en in self.energy_vals])
        idx, _ = find_peaks(self.sig, height=thr, distance=3)
        peaks = []
        for i in idx:
            en = float(self.energy_vals[i])
            # Filter response at channel i to a unit-area peak at energy en
            unit = self.F[i] @ gaussian_template(self.energy_vals, [(en, 1.0)], self.sigma_fn)
            net_area = float(self.r_f[i] / unit) if unit > 0 else 0.0
            peaks.append(ResidualPeak(int(i), en, float(self.sig[i]), net_area))
        peaks.sort(key=lambda p: -p.net_area)
        return peaks


@dataclass
class _Screening:
    """Screening of a candidate element against the filtered residual."""
    element: str
    matched_line: str
    matched_ref: str
    line_energy: float
    score: float  # Decrease of the filtered chi-square
    family_amplitudes: Dict[str, float]
    family_significance: Dict[str, float]
    passed: bool
    reason: str = ''
    offsets: Dict[str, float] = field(default_factory=dict)  # Tied energy shift (keV) of each family


def _n_resolved_peaks(el: str, beam_energy: float, energy_range: Tuple[float, float], sigma_fn: Callable,
                      min_weight: float, min_energy: float) -> int:
    """
    Number of resolved peaks of an element in the fitted range, above `min_energy` (lines at low energy are
    strongly absorbed and rarely form usable peaks): lines less than 2 detector sigmas apart are merged, and
    a merged peak counts if its total weight (relative to the reference line of its family) is at least
    `min_weight`.
    """
    lines = sorted((en, w) for fam in element_line_families(el, beam_energy, energy_range).values()
                   for _, en, w in fam if en >= min_energy)
    peaks: List[List[float]] = []  # [last energy, total weight]
    for en, w in lines:
        if peaks and en - peaks[-1][0] <= 2 * float(sigma_fn(en)):
            peaks[-1][0] = en
            peaks[-1][1] += w
        else:
            peaks.append([en, w])
    return sum(1 for _, w in peaks if w >= min_weight)


@dataclass
class _PeakEvidence:
    """Local evidence for a candidate at each of its resolved peaks (see PeakIDSettings.min_peak_gain)."""
    supported: List[float] = field(default_factory=list)  # Energies (keV) of peaks supporting the candidate
    contradicted: List[float] = field(default_factory=list)  # Energies of predicted peaks absent in the data
    local_fraction: float = 1.0  # Fraction of the total fit improvement coming from peak regions (not background)

    def describe(self) -> str:
        msg = f"{len(self.supported)} supporting peaks, {self.local_fraction:.0%} of the gain within peak regions"
        if self.contradicted:
            msg += f", predicted peaks absent at {', '.join(f'{en:.2f}' for en in self.contradicted)} keV"
        return msg


def _peak_evidence(
    el: str,
    offsets: Dict[str, float],
    y: np.ndarray,
    energy_vals: np.ndarray,
    model_without: np.ndarray,
    outcome: FitOutcome,
    var: np.ndarray,
    beam_energy: float,
    sigma_fn: Callable,
    settings: PeakIDSettings,
) -> _PeakEvidence:
    """
    Compare the fit with and without a candidate within each of its resolved peaks (lines above the low-energy
    threshold, merged if less than 2 detector sigmas apart). A peak supports the candidate if the fit improves
    locally; it contradicts it if the fit worsens locally and its predicted counts (from the fitted contribution
    of the candidate) are above the detection limit.
    """
    e_range = (float(energy_vals[0]), float(energy_vals[-1]))
    families = element_line_families(el, beam_energy, e_range)
    contribution = sum((outcome.components[f"{el}_{ref}"] for ref in families if f"{el}_{ref}" in outcome.components),
                       np.zeros_like(y))
    lines = sorted((en + offsets.get(ref, 0.0), w) for ref, fam in families.items()
                   for _, en, w in fam if en >= settings.low_energy_threshold_keV)
    peaks: List[List[float]] = []  # [lowest energy, highest energy, total weight]
    for en, w in lines:
        if peaks and en - peaks[-1][1] <= 2 * float(sigma_fn(en)):
            peaks[-1][1] = en
            peaks[-1][2] += w
        else:
            peaks.append([en, en, w])
    evidence = _PeakEvidence()
    gain_per_channel = ((y - model_without) ** 2 - (y - outcome.model) ** 2) / var
    total_gain = float(np.sum(gain_per_channel))
    max_loss = max(settings.min_peak_gain, settings.contradiction_rel_loss * total_gain)
    # Peak regions, for the fraction of the gain within peaks: channels where the candidate adds counts, or where
    # the peaks of any element are significant compared with the background (the gain of an overlapping element
    # also comes from the readjustment of the peaks it overlaps). The rest is gain from the background.
    in_peaks = contribution > max(0.5, 0.02 * float(contribution.max(initial=0.0)))
    if outcome.background is not None and outcome.components:
        all_peaks = sum(outcome.components.values())
        in_peaks |= all_peaks > 0.1 * np.maximum(outcome.background, 1.0)
    for lo, hi, w in peaks:
        # Peaks are selected by their predicted counts in the fit, not by table weights (e.g. Pt Mz is weak in
        # the table, but visible in the fit)
        center = 0.5 * (lo + hi)
        half = 1.2 * FWHM_PER_SIGMA * float(sigma_fn(center))
        win = (energy_vals >= lo - half) & (energy_vals <= hi + half)
        gain = float(np.sum(gain_per_channel[win]))
        if gain >= settings.min_peak_gain:
            evidence.supported.append(center)
        elif gain <= -max_loss:
            # Only peaks predicted above the detection limit can contradict the candidate
            if contribution[win].sum() >= 2.71 + 3.29 * np.sqrt(max(model_without[win].sum(), 0.0)):
                evidence.contradicted.append(center)
    if total_gain > 0:
        evidence.local_fraction = float(np.sum(gain_per_channel[in_peaks])) / total_gain
    return evidence


def filtered_gain(y: np.ndarray, model_without: np.ndarray, model_with: np.ndarray, F: np.ndarray,
                  var: np.ndarray) -> float:
    """Chi-square improvement of the top-hat filtered residual: fit improvement of peaks, not of the background."""
    var_f = (F ** 2) @ var
    ok = var_f > 0
    r0 = (F @ (y - model_without))[ok]
    r1 = (F @ (y - model_with))[ok]
    return float(np.sum((r0 ** 2 - r1 ** 2) / var_f[ok]))


def _residual_peak_check(el: str, outcome_without: FitOutcome, y: np.ndarray, energy_vals: np.ndarray,
                         F: np.ndarray, sigma_fn: Callable, beam_energy: float,
                         settings: PeakIDSettings) -> Tuple[float, bool]:
    """
    Largest positive peak of the filtered residual of a fit without an element, within one FWHM of its lines
    (weight >= min_line_weight), and whether it reaches min_peak_significance (with the threshold multipliers of
    that line's energy) at some line: "an element is added only if it is clearly a peak in the residual".
    """
    e_range = (float(energy_vals[0]), float(energy_vals[-1]))
    res = _ResidualAnalysis(y, outcome_without.model, energy_vals, F, sigma_fn, settings.model_rel_uncertainty)
    best, ok = 0.0, False
    for lines in element_line_families(el, beam_energy, e_range).values():
        for _, en, w in lines:
            if w < settings.min_line_weight:
                continue
            win = np.abs(energy_vals - en) <= FWHM_PER_SIGMA * float(sigma_fn(en))
            if not np.any(win):
                continue
            sig = float(res.sig[win].max())
            best = max(best, sig)
            if sig >= settings.min_peak_significance * settings.multiplier(el, en):
                ok = True
    return best, ok


def background_band_ratios(y: np.ndarray, energy_vals: np.ndarray, outcome: FitOutcome,
                           bands: Sequence[Tuple[float, float]], exclude: Optional[np.ndarray] = None,
                           min_model_counts: float = 50.0) -> List[Tuple[Tuple[float, float], float]]:
    """
    Ratio of data to model in each energy band, in channels without significant peaks (excluded with `exclude`).
    Bands whose peak-free model counts are below `min_model_counts` are skipped.
    """
    if outcome.background is None:
        return []
    peaks = sum(outcome.components.values()) if outcome.components else np.zeros_like(y)
    peak_free = peaks < 0.1 * np.maximum(outcome.background, 1.0)
    if exclude is not None:
        peak_free &= ~exclude
    ratios = []
    for lo, hi in bands:
        mask = peak_free & (energy_vals >= lo) & (energy_vals < hi)
        model_sum = float(outcome.model[mask].sum())
        if model_sum >= min_model_counts:
            ratios.append(((lo, hi), float(y[mask].sum()) / model_sum))
    return ratios


def _format_band_ratios(ratios: List[Tuple[Tuple[float, float], float]]) -> str:
    return ", ".join(f"{r:.2f} ({lo:g}–{min(hi, 99):g} keV)" if hi < 99 else f"{r:.2f} (>{lo:g} keV)"
                     for (lo, hi), r in ratios)


def _tied_offsets(
    el: str,
    current: Optional[FitOutcome],
    beam_energy: float,
    energy_range: Tuple[float, float],
    sigma_fn: Callable,
) -> Dict[str, float]:
    """
    Energy shift of each line family of a candidate, tied to the overlapping peaks of the current fit as in
    Peaks_Model._fix_overlapping_ref_peaks: a reference line of the candidate overlapping fitted reference
    peaks (separation below OVERLAP_SEPARATION_SIGMAS, at tabulated energies) takes the shift of the
    largest of them (the anchor). Families without overlaps are not shifted.
    """
    offsets: Dict[str, float] = {}
    if current is None or not current.line_centers:
        return offsets
    model_refs = []
    for model_ref, center in current.line_centers.items():
        m_el, m_line = model_ref.split('_', 1)
        try:
            e_tab = float(_xray_lines(m_el)[m_line]['energy (keV)'])
        except (KeyError, ValueError):
            continue
        model_refs.append((e_tab, center - e_tab, current.line_areas.get(model_ref, 0.0)))
    for ref, lines in element_line_families(el, beam_energy, energy_range).items():
        en_ref = next(en for line, en, _ in lines if line == ref)
        group = [
            (area, offset) for e_tab, offset, area in model_refs
            if abs(en_ref - e_tab) < OVERLAP_SEPARATION_SIGMAS * np.sqrt((sigma_fn(en_ref) ** 2 + sigma_fn(e_tab) ** 2) / 2)
        ]
        if group:
            offsets[ref] = max(group)[1]
    return offsets


def _confirm_overvoltage(settings: PeakIDSettings) -> float:
    """Overvoltage above which a higher-energy line family must be observed (see PeakIDSettings)."""
    if settings.confirm_overvoltage is not None:
        return settings.confirm_overvoltage
    # Imported here: the quantifier imports this module
    from autoemx.core.quantifier.quantifier import XSp_Quantifier
    return XSp_Quantifier.ideal_ref_line_overvoltage


def _screen_candidate(
    el: str,
    matched_line: str,
    matched_ref: str,
    line_energy: float,
    res: _ResidualAnalysis,
    beam_energy: float,
    energy_range: Tuple[float, float],
    settings: PeakIDSettings,
    offsets: Optional[Dict[str, float]] = None,
) -> _Screening:
    """
    Fit the templates of the line families of a candidate to the filtered residual (NNLS).

    `offsets` are the energy shifts of the families, tied to overlapping peaks of the current fit.
    """
    offsets = offsets or {}
    families = element_line_families(el, beam_energy, energy_range)
    refs, cols = [], []
    for ref, lines in families.items():
        shift = offsets.get(ref, 0.0)
        t_f = res.F @ gaussian_template(res.energy_vals, [(en + shift, w) for _, en, w in lines], res.sigma_fn)
        if np.any(t_f[res.valid] != 0):
            refs.append(ref)
            cols.append(t_f)
    if matched_ref not in refs:
        return _Screening(el, matched_line, matched_ref, line_energy, 0.0, {}, {}, False, 'line outside fitted range', offsets)

    v = res.valid
    w = 1 / np.sqrt(res.var_f[v])
    A = np.column_stack(cols)[v] * w[:, None]
    b = res.r_f[v] * w
    a, rnorm = nnls(A, b)
    score = float(b @ b - rnorm ** 2)

    # Uncertainty of the amplitudes (sandwich estimator, accounting for correlated filtered channels)
    signif = {ref: 0.0 for ref in refs}
    active = np.flatnonzero(a > 0)
    if active.size:
        Aa = A[:, active]
        M = np.linalg.pinv(Aa.T @ Aa)
        B = np.zeros((len(res.r_f), active.size))
        B[v] = Aa * w[:, None]  # Rows of the valid channels, weighted
        cov_a = M @ res.sandwich(B) @ M
        for k, j in enumerate(active):
            sd = np.sqrt(max(cov_a[k, k], 1e-12))
            signif[refs[j]] = float(a[j] / sd)
    amps = {ref: float(a[j]) for j, ref in enumerate(refs)}

    mult = settings.multiplier(el, line_energy)
    if signif[matched_ref] < settings.min_peak_significance * mult:
        return _Screening(el, matched_line, matched_ref, line_energy, score, amps, signif, False,
                          f'{matched_ref} family at {signif[matched_ref]:.1f}σ', offsets)

    # Confirmation: well-excited families at higher energy must be observed
    ref_energy = {ref: next(en for line, en, _ in families[ref] if line == ref) for ref in refs}
    for ref in refs:
        if ref_energy[ref] > ref_energy[matched_ref] and beam_energy / ref_energy[ref] >= _confirm_overvoltage(settings):
            if signif[ref] < settings.confirm_significance:
                return _Screening(el, matched_line, matched_ref, line_energy, score, amps, signif, False,
                                  f'missing confirmation line {el} {ref}', offsets)
    return _Screening(el, matched_line, matched_ref, line_energy, score, amps, signif, True, '', offsets)


def _assign_peaks(
    y: np.ndarray,
    energy_vals: np.ndarray,
    F: np.ndarray,
    sigma_fn: Callable,
    beam_energy: float,
    energy_range: Tuple[float, float],
    line_table: Tuple[Tuple[float, str, str, str], ...],
    settings: PeakIDSettings,
    known_elements: Iterable[str],
    excluded: set,
) -> List[_Screening]:
    """
    Assign elements to the peaks of the raw spectrum, without fitting, to build the starting fit when no
    sample element is known (a fit with the substrate elements only cannot describe the spectrum).

    The line families of the known (substrate) and assigned elements are fitted to the top-hat filtered
    spectrum (NNLS, linear); the most intense peak left in the filtered residual is given the candidate with
    the best screening score (whole line families, with the confirmation lines), and the fit is repeated.
    The assigned elements are verified afterwards like any other identified element (full-fit pruning,
    residual requirement, tiers).
    """
    res = _ResidualAnalysis(y, None, energy_vals, F, sigma_fn, settings.model_rel_uncertainty)
    r_f_raw = res.r_f.copy()
    v = res.valid
    w = 1 / np.sqrt(res.var_f[v])
    table_energies = np.array([row[0] for row in line_table])
    known = list(dict.fromkeys(known_elements))
    cols: List[np.ndarray] = []

    def add_columns(el: str) -> None:
        for lines in element_line_families(el, beam_energy, energy_range).values():
            t_f = F @ gaussian_template(energy_vals, [(en, wt) for _, en, wt in lines], sigma_fn)
            if np.any(t_f[v] != 0):
                cols.append(t_f)

    def update_residual() -> None:
        r_f = r_f_raw
        if cols:
            A = np.column_stack(cols)
            a, _ = nnls(A[v] * w[:, None], r_f_raw[v] * w)
            r_f = r_f_raw - A @ a
        res.r_f = r_f
        res.sig = np.zeros_like(r_f)
        res.sig[v] = r_f[v] / np.sqrt(res.var_f[v])

    for el in known:
        add_columns(el)
    assigned: List[_Screening] = []
    unassignable: List[float] = []
    for _ in range(settings.max_iterations):
        if len(assigned) >= settings.max_new_elements:
            break
        update_residual()
        peaks = [p for p in res.find_peaks(settings)
                 if all(abs(p.energy_keV - en) > 2 * float(sigma_fn(en)) for en in unassignable)]
        if not peaks:
            break
        peak = peaks[0]
        window = settings.candidate_window_sigma * float(sigma_fn(peak.energy_keV)) + settings.merged_peak_shift_keV
        lo, hi = np.searchsorted(table_energies, [peak.energy_keV - window, peak.energy_keV + window])
        candidates: Dict[str, Tuple[str, str, float]] = {}
        for en, el, line, ref in line_table[lo:hi]:
            if el in known or el in excluded:
                continue
            if el not in candidates or abs(en - peak.energy_keV) < abs(candidates[el][2] - peak.energy_keV):
                candidates[el] = (line, ref, en)
        passed = [sc for sc in (_screen_candidate(el, line, ref, en, res, beam_energy, energy_range, settings)
                                for el, (line, ref, en) in candidates.items()) if sc.passed]
        if not passed:
            unassignable.append(peak.energy_keV)  # Left to the identification after the starting fit
            continue
        best = max(passed, key=lambda sc: sc.score)
        assigned.append(best)
        known.append(best.element)
        add_columns(best.element)
    return assigned


# =============================================================================
# Fit function using XSp_Fitter
# =============================================================================
def make_xsp_fit_function(
    spectrum_vals: np.ndarray,
    energy_vals: np.ndarray,
    spectrum_lims: Tuple[int, int],
    microscope_ID: str,
    meas_mode: str,
    det_ch_offset: float,
    det_ch_width: float,
    beam_energy: float,
    emergence_angle: float,
    tot_sp_counts: float,
    substrate_elements: Sequence[str] = (),
    is_particle: bool = False,
    sp_collection_time: Optional[float] = None,
    fit_tolerance: float = 1e-3,
    initial_par_vals: Optional[Dict[str, float]] = None,
    fixed_par_vals: Optional[Dict[str, float]] = None,
) -> Callable[[Sequence[str]], FitOutcome]:
    """
    Fit function for `identify_elements`, using XSp_Fitter.

    `spectrum_vals` and `energy_vals` are the values within the spectrum limits, while `tot_sp_counts` must be
    the total counts of the full spectrum, as in XSp_Quantifier: the fitter derives the initial background scale
    from it, and the background fit is sensitive to this starting value. The returned callable
    takes the sample elements (fitted as elements to quantify; substrate elements are added) and returns
    a FitOutcome. Results are cached by element set.

    `initial_par_vals` are passed to XSp_Fitter.fit_spectrum, e.g. the starting background scale
    {'K': XSp_Quantifier.get_starting_K_val()} for particles, as in quantification (without it, particle
    geometry parameters may compensate for the background intensity and the fit may end with too little
    background). `fixed_par_vals` are set and not varied, e.g. {'K': K} to keep the background scale at the
    value fitted above 5 keV, when the particle geometry parameters drift the background away from the data.
    """
    # Imported here to avoid a circular import (the quantifier imports this module)
    from autoemx.core.fitter.fitter import XSp_Fitter
    from autoemx.core.quantifier.quantifier import XSp_Quantifier

    y = np.asarray(spectrum_vals, dtype=float)
    E = np.asarray(energy_vals, dtype=float)
    ch_width = E[1] - E[0]
    cache: Dict[Tuple[str, ...], FitOutcome] = {}

    def fit(elements: Sequence[str]) -> FitOutcome:
        key = tuple(sorted(set(elements)))
        if key in cache:
            return cache[key]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fitter = XSp_Fitter(
                y, E, spectrum_lims, microscope_ID, meas_mode, det_ch_offset, det_ch_width,
                beam_energy, emergence_angle,
                is_particle=is_particle,
                els_to_quantify=list(key),
                els_substrate=list(substrate_elements),
                tot_sp_counts=tot_sp_counts,
                sp_collection_time=sp_collection_time,
                xray_quant_ref_lines=XSp_Quantifier.xray_quant_ref_lines,
            )
            if fixed_par_vals:
                fitter._make_spectrum_mod_pars()
                for name, val in fixed_par_vals.items():
                    if name in fitter.spectrum_pars:
                        fitter.spectrum_pars[name].set(value=val, vary=False)
                result, _ = fitter.fit_spectrum(parameters=fitter.spectrum_pars, function_tolerance=fit_tolerance,
                                                initial_par_vals=initial_par_vals)
            else:
                result, _ = fitter.fit_spectrum(function_tolerance=fit_tolerance, initial_par_vals=initial_par_vals)
        areas = {}
        for name, par in result.params.items():
            if name.endswith('_area'):
                el_line = name[:-len('_area')]
                if XSp_Fitter.escape_peaks_str in el_line or XSp_Fitter.pileup_peaks_str in el_line:
                    continue
                # Peak models are evaluated per channel, so area / channel width gives counts
                areas[el_line] = float(par.value) / ch_width
        # Group the peak components by line family (reference line); the rest of the model is background
        model = np.asarray(result.best_fit, dtype=float)
        components: Dict[str, np.ndarray] = {}
        for prefix, comp in result.eval_components(x=E).items():
            ref = fitter.el_lines_weight_refs_dict.get(prefix.rstrip('_'))
            if ref is not None:
                components[ref] = components.get(ref, 0) + np.broadcast_to(comp, E.shape)
        background = model - sum(components.values()) if components else model.copy()
        centers = {ref: float(result.params[f"{ref}_center"].value)
                   for ref in components if f"{ref}_center" in result.params}
        outcome = FitOutcome(key, model, poisson_deviance(y, model), areas, components, background,
                             line_centers=centers, redchi=float(result.redchi))
        cache[key] = outcome
        return outcome

    return fit


def _linear_variance(outcome: FitOutcome, rel_uncertainty: float) -> np.ndarray:
    """Variance of the spectrum for the linear tests: Poisson from the current model plus model uncertainty."""
    return np.maximum(outcome.model, 1.0) + (rel_uncertainty * outcome.model) ** 2


def _weighted_chi2(y: np.ndarray, model: np.ndarray, var: np.ndarray) -> float:
    return float(np.sum((y - model) ** 2 / var))


def _linear_refit(
    y: np.ndarray,
    base: FitOutcome,
    var: np.ndarray,
    add: Optional[Dict[str, np.ndarray]] = None,
    drop: Iterable[str] = (),
    elements: Tuple[str, ...] = (),
) -> FitOutcome:
    """
    Weighted non-negative linear fit of the components of a fit (fixed shapes, free amplitudes), with
    components added (unit-area templates, so that their amplitudes are areas in counts) or dropped.

    Peak shapes, energies and background shape are those of the base fit, so the result is an
    approximation of a full refit, which is good enough to compare candidates.
    """
    add = add or {}
    drop = set(drop)
    names = [k for k in base.components if k not in drop] + list(add)
    cols = [base.background] + [base.components[k] for k in names if k not in add] + [add[k] for k in add]
    A = np.column_stack(cols)
    w = 1 / np.sqrt(var)
    coef, _ = nnls(A * w[:, None], y * w)
    model = A @ coef
    components = {name: coef[i + 1] * cols[i + 1] for i, name in enumerate(names)}
    line_areas = {}
    for el_line, area in base.line_areas.items():
        ref = el_line if el_line in base.components else None
        if ref is not None and ref not in drop:
            line_areas[el_line] = area * coef[1 + names.index(ref)]
    for name in add:
        line_areas[name] = float(coef[1 + names.index(name)])
    centers = {ref: en for ref, en in base.line_centers.items() if ref not in drop}
    return FitOutcome(tuple(sorted(elements)), model, poisson_deviance(y, model), line_areas,
                      components, coef[0] * base.background, approximate=True, line_centers=centers)


# =============================================================================
# Identification loop
# =============================================================================
def _artifact_note(energy: float, outcome: Optional[FitOutcome], sigma_fn: Callable) -> str:
    """Note on a possible sum or escape peak of the major fitted lines at a given energy."""
    if outcome is None or not outcome.line_areas:
        return ''
    majors = sorted(outcome.line_areas.items(), key=lambda kv: -kv[1])[:6]
    major_lines = []
    for el_line, area in majors:
        el, line = el_line.split('_', 1)
        try:
            major_lines.append((el_line, float(_xray_lines(el)[line]['energy (keV)'])))
        except KeyError:
            continue
    tol = 2 * float(sigma_fn(energy))
    for i, (name_a, en_a) in enumerate(major_lines):
        if abs(en_a - SI_KA_ENERGY - energy) < tol:
            return f'possible Si escape peak of {name_a}'
        for name_b, en_b in major_lines[i:]:
            if abs(en_a + en_b - energy) < tol:
                return f'possible sum peak {name_a} + {name_b}'
    return ''


def identify_elements(
    spectrum_vals: np.ndarray,
    energy_vals: np.ndarray,
    fit_function: Callable[[Sequence[str]], FitOutcome],
    beam_energy: float,
    sigma_fn: Callable,
    fixed_elements: Sequence[str] = (),
    substrate_elements: Sequence[str] = (),
    settings: Optional[PeakIDSettings] = None,
    has_standard: Optional[Callable[[str], bool]] = None,
    excluded_elements: Iterable[str] = (),
    initial_outcome: Optional[FitOutcome] = None,
    priority_elements: Iterable[str] = (),
    background_refit_function: Optional[Callable[[Sequence[str]], FitOutcome]] = None,
    meas_type: str = 'EDS',
    verbose: bool = False,
) -> ElementIDResult:
    """
    Identify the elements present in a spectrum, in addition to the fixed and substrate elements.

    Parameters
    ----------
    spectrum_vals, energy_vals : array
        Spectrum counts and energies (keV), within the fitted spectrum limits.
    fit_function : callable
        Takes a list of sample elements and returns a FitOutcome of the spectrum fitted with those elements
        (substrate elements are handled by the fit function). See `make_xsp_fit_function`.
    beam_energy : float
        Beam energy (keV).
    sigma_fn : callable
        Detector sigma (keV) as a function of energy (keV). See `get_sigma_function`.
    fixed_elements : list of str
        Elements known to be present. May be empty: the peaks of the raw spectrum are then assigned first
        (see `_assign_peaks`), and the starting fit includes the assigned elements, which are verified like
        the other identified elements.
    substrate_elements : list of str
        Substrate elements, always part of the model. Their peaks cannot be attributed to the sample.
    settings : PeakIDSettings, optional
        Identification thresholds.
    has_standard : callable, optional
        Returns whether an element can be quantified. Used only for reporting.
    excluded_elements : iterable of str
        Elements never proposed (in addition to settings.excluded_elements).
    initial_outcome : FitOutcome, optional
        Fit of the spectrum with the fixed elements, if already available.
    background_refit_function : callable, optional
        Fit function with the background scale fixed (e.g. make_xsp_fit_function with fixed_par_vals={'K': K},
        K from the fit above 5 keV). Used instead of `fit_function` if the background of the starting fit is
        outside settings.background_ratio_bounds.
    priority_elements : iterable of str
        Pool of elements likely to be present (e.g. identified in other spectra of the same sample).
        For each residual peak, pool elements that pass the screening are verified first, and the other
        candidates are fitted only if no pool element explains the peak. This saves full fits, at the cost
        of a bias towards the pool for overlapping candidates.
    meas_type : str
        Measurement method, for the overlap table.

    Returns
    -------
    ElementIDResult
    """
    t_start = time.time()
    settings = settings or PeakIDSettings()
    y = np.asarray(spectrum_vals, dtype=float)
    E = np.asarray(energy_vals, dtype=float)
    e_range = (float(E[0]), float(E[-1]))
    F = top_hat_matrix(E, sigma_fn)
    line_table = _strong_lines_table(float(beam_energy), *e_range, settings.min_line_weight)
    table_energies = np.array([row[0] for row in line_table])

    fixed = list(dict.fromkeys(fixed_elements))
    substrate = [el for el in dict.fromkeys(substrate_elements) if el not in fixed]
    excluded = set(settings.excluded_elements) | set(excluded_elements)
    excluded |= {Element.from_Z(z).symbol for z in range(MIN_Z, MAX_Z + 1)
                 if not settings.z_range[0] <= z <= settings.z_range[1]}
    priority = set(priority_elements) - excluded
    accepted: List[IdentifiedElement] = []
    unexplained: List[UnexplainedPeak] = []
    n_fits = 0
    n_tests = 0

    def fit(elements):
        nonlocal n_fits
        n_fits += 1
        return fit_function(list(elements))

    def model_elements():
        return fixed + [a.element for a in accepted]

    if initial_outcome is None and not fixed:
        # No known sample element: assign the peaks of the raw spectrum first, and start from a fit with them
        for sc in _assign_peaks(y, E, F, sigma_fn, float(beam_energy), e_range, line_table, settings,
                                substrate, excluded):
            accepted.append(IdentifiedElement(
                element=sc.element, status='identified', matched_line=sc.matched_line,
                line_energy_keV=sc.line_energy, peak_significance=sc.family_significance.get(sc.matched_ref),
                from_peak_assignment=True,
            ))
        if verbose and accepted:
            logger.info(f"Peaks assigned before the starting fit: {', '.join(a.element for a in accepted)}")

    if initial_outcome is not None:
        current = initial_outcome
    elif model_elements() or substrate:
        current = fit(model_elements())
    else:
        current = None

    # Background health of the starting fit; refit with fixed background scale if needed
    background_note = ''
    background_reliable = True
    lo_b, hi_b = settings.background_ratio_bounds

    def peak_exclusion(outcome: FitOutcome) -> np.ndarray:
        """Channels within 2 FWHM of significant unexplained peaks of the residual."""
        exclude = np.zeros_like(y, dtype=bool)
        res0 = _ResidualAnalysis(y, outcome.model, E, F, sigma_fn, settings.model_rel_uncertainty)
        for p in res0.find_peaks(settings):
            exclude |= np.abs(E - p.energy_keV) <= 2 * FWHM_PER_SIGMA * float(sigma_fn(p.energy_keV))
        return exclude

    def band_ratios(outcome: FitOutcome):
        return background_band_ratios(y, E, outcome, settings.background_bands_keV, peak_exclusion(outcome))

    def worst_deviation(ratios) -> float:
        return max((abs(np.log(r)) for _, r in ratios if r > 0), default=0.0)

    def has_deficit(ratios) -> bool:
        return any(r > hi_b for _, r in ratios)

    ratios = band_ratios(current) if current is not None else []
    if any(not lo_b <= r <= hi_b for _, r in ratios):
        background_note = f"data/model in peak-free channels {_format_band_ratios(ratios)} in the starting fit"
        if background_refit_function is not None:
            refit_start = background_refit_function(model_elements())
            n_fits += 1
            refit_ratios = band_ratios(refit_start)
            creates_deficit = has_deficit(refit_ratios) and not has_deficit(ratios)
            if not creates_deficit and worst_deviation(refit_ratios) < worst_deviation(ratios):
                # The refit with fixed background scale is better (and does not create a background deficit,
                # the harmful direction): use it for all further fits
                fit_function = background_refit_function
                current, ratios = refit_start, refit_ratios
                background_note += f"; {_format_band_ratios(ratios)} after refitting with fixed background scale"
            elif creates_deficit:
                background_note += (f"; refit with fixed background scale not used (it creates a background "
                                    f"deficit: {_format_band_ratios(refit_ratios)})")
            else:
                background_note += "; refit with fixed background scale not better"
        if has_deficit(ratios):
            background_reliable = False
            background_note += " — background unreliable, identified elements only reported as possible"
        if verbose:
            logger.info(f"Background: {background_note}")

    if (settings.max_start_redchi_fraction is not None and current is not None and current.redchi is not None
            and current.redchi > settings.max_start_redchi_fraction * y.sum()):
        result = ElementIDResult(
            fixed_elements=fixed, substrate_elements=substrate,
            elements=[IdentifiedElement(el, 'fixed') for el in fixed] + [IdentifiedElement(el, 'substrate') for el in substrate],
            unexplained_peaks=[], energy_vals=E, spectrum_vals=y, model=current.model, n_fits=n_fits,
            elapsed_s=time.time() - t_start,
            skipped_reason=f"poor starting fit (reduced chi-square {current.redchi:.0f} > "
                           f"{settings.max_start_redchi_fraction * 100:.2f}% of total counts)",
        )
        if verbose:
            logger.info(result.summary())
        return result

    stats: Dict[Tuple[str, ...], float] = {}

    def stat(outcome: FitOutcome) -> float:
        """Fit statistic of an outcome (decisions use it instead of the Poisson deviance, see PeakIDSettings)."""
        if outcome.approximate:
            return fit_statistic(y, outcome.model, settings.model_rel_uncertainty)
        if outcome.elements not in stats:
            stats[outcome.elements] = fit_statistic(y, outcome.model, settings.model_rel_uncertainty)
        return stats[outcome.elements]

    def family_templates(el: str, offsets: Dict[str, float]) -> Dict[str, np.ndarray]:
        """Unit-area (reference line) templates of the line families of an element, with tied shifts."""
        return {
            f"{el}_{ref}": gaussian_template(E, [(en + offsets.get(ref, 0.0), w) for _, en, w in lines], sigma_fn)
            for ref, lines in element_line_families(el, beam_energy, e_range).items()
        }

    def _overlaps_model(s: _Screening) -> bool:
        """
        Whether a reference line of a candidate overlaps a fitted reference peak of the current model, with
        the criterion of the fitter (Peaks_Model._fix_overlapping_ref_peaks). Without fitted centers, the
        overlap table is used instead.
        """
        if current is None:
            return False
        if not current.line_centers:
            in_model = set(model_elements()) | set(substrate)
            return any(p_el in in_model for p_el, _, _ in _overlap_partners(s.element, s.matched_ref, sigma_fn, meas_type))
        ref_energies = [next(en for line, en, _ in lines if line == ref)
                        for ref, lines in element_line_families(s.element, beam_energy, e_range).items()]
        for center in current.line_centers.values():
            for en in ref_energies:
                max_sep = OVERLAP_SEPARATION_SIGMAS * np.sqrt((sigma_fn(en) ** 2 + sigma_fn(center) ** 2) / 2)
                if abs(en - center) < max_sep:
                    return True
        return False

    def _verify(to_fit: List[_Screening], in_pool: bool, in_model: set) -> list:
        """
        Test each screened candidate; returns those passing the fit-statistic and Currie criteria.

        If the current fit provides its components, a candidate is tested with a linear fit of the current
        components (fixed shapes, free amplitudes) and of the candidate's family templates, which takes
        milliseconds. Each candidate is tested with a full fit instead if there is no starting model, or if
        a candidate overlaps a peak already in the model: the fitter then ties the energy shifts of the
        overlapping peaks (Peaks_Model._fix_overlapping_ref_peaks), the fixed peak shapes of the current fit
        may have adapted to the missing element, and only a full refit can redistribute the counts between
        the overlapping peaks.
        """
        nonlocal n_tests
        n_tests += len(to_fit)
        verified = []
        overlapped = settings.full_fit_overlaps and any(_overlaps_model(s) for s in to_fit)
        linear = current is not None and current.background is not None and not overlapped
        if not linear:
            to_fit = to_fit[:settings.max_full_fits_per_peak]
        if linear:
            var = _linear_variance(current, settings.model_rel_uncertainty)
            base = _linear_refit(y, current, var, elements=current.elements)
            chi_base = _weighted_chi2(y, base.model, var)
        for s in to_fit:
            if linear:
                outcome = _linear_refit(y, current, var, add=family_templates(s.element, s.offsets),
                                        elements=tuple(model_elements()) + (s.element,))
                d_dev = chi_base - _weighted_chi2(y, outcome.model, var)
            else:
                outcome = fit(model_elements() + [s.element])
                d_dev = np.inf if current is None else stat(current) - stat(outcome)
            mult = settings.multiplier(s.element, s.line_energy)
            area = outcome.line_areas.get(f"{s.element}_{s.matched_ref}", 0.0)
            in_win = np.abs(E - s.line_energy) <= 1.2 * FWHM_PER_SIGMA * float(sigma_fn(s.line_energy))
            if current is not None:
                background = current.model[in_win].sum()
            else:
                background = max(y[in_win].sum() - area, 0.0)
            l_d = 2.71 + 3.29 * np.sqrt(max(background, 0.0))
            ok = d_dev >= settings.min_delta_deviance * mult and area >= l_d * mult
            evidence = None
            prefilter_ok = d_dev >= settings.linear_prefilter_fraction * settings.min_delta_deviance * mult
            if current is not None and (ok or (linear and prefilter_ok and area >= l_d * mult)):
                if linear:
                    evidence = _peak_evidence(s.element, s.offsets, y, E, base.model, outcome, var,
                                              beam_energy, sigma_fn, settings)
                    weak = d_dev < settings.min_delta_deviance * mult
                    if ((evidence.contradicted or weak) and _overlaps_model(s)
                            and d_dev >= settings.linear_prefilter_fraction * settings.min_delta_deviance * mult):
                        # The frozen shapes of overlapping peaks of the current fit may have adapted to the
                        # missing element, so the linear test cannot judge those peaks: decide with a full fit
                        outcome = fit(model_elements() + [s.element])
                        d_dev = stat(current) - stat(outcome)
                        area = outcome.line_areas.get(f"{s.element}_{s.matched_ref}", 0.0)
                        evidence = _peak_evidence(s.element, s.offsets, y, E, current.model, outcome,
                                                  _linear_variance(current, settings.model_rel_uncertainty),
                                                  beam_energy, sigma_fn, settings)
                        ok = d_dev >= settings.min_delta_deviance * mult and area >= l_d * mult
                else:
                    evidence = _peak_evidence(s.element, s.offsets, y, E, current.model, outcome,
                                              _linear_variance(current, settings.model_rel_uncertainty),
                                              beam_energy, sigma_fn, settings)
                ok = (ok and not evidence.contradicted
                      and evidence.local_fraction >= settings.min_local_gain_fraction)
            if verbose:
                pool_str = ' [pool]' if in_pool else ''
                ev_str = f", {evidence.describe()}" if evidence is not None else ''
                logger.info(f"   {s.element}{pool_str}: ΔD = {d_dev:.0f}, area = {area:.0f} (L_D = {l_d:.0f}){ev_str} {'✓' if ok else '✗'}")
            if ok:
                verified.append((d_dev, s, outcome, area, l_d))
        return verified

    for _ in range(settings.max_iterations):
        if len(accepted) >= settings.max_new_elements:
            break
        res = _ResidualAnalysis(y, None if current is None else current.model, E, F, sigma_fn,
                                settings.model_rel_uncertainty)
        peaks = [
            p for p in res.find_peaks(settings)
            if all(abs(p.energy_keV - u.energy_keV) > 2 * sigma_fn(u.energy_keV) for u in unexplained)
        ]
        if not peaks:
            break
        peak = peaks[0]
        in_model = set(model_elements()) | set(substrate)

        # Candidates from the line table
        window = settings.candidate_window_sigma * float(sigma_fn(peak.energy_keV)) + settings.merged_peak_shift_keV
        lo, hi = np.searchsorted(table_energies, [peak.energy_keV - window, peak.energy_keV + window])
        candidates: Dict[str, Tuple[str, str, float]] = {}
        for en, el, line, ref in line_table[lo:hi]:
            if el in in_model or el in excluded:
                continue
            if el not in candidates or abs(en - peak.energy_keV) < abs(candidates[el][2] - peak.energy_keV):
                candidates[el] = (line, ref, en)

        screened = [
            _screen_candidate(el, line, ref, en, res, beam_energy, e_range, settings,
                              _tied_offsets(el, current, beam_energy, e_range, sigma_fn))
            for el, (line, ref, en) in candidates.items()
        ]
        passed = sorted((s for s in screened if s.passed), key=lambda s: -s.score)
        if verbose:
            logger.info(f"Residual peak at {peak.energy_keV:.3f} keV ({peak.significance:.1f}σ, ~{peak.net_area:.0f} counts): "
                        + ", ".join(f"{s.element} {'✓' if s.passed else '✗ ' + s.reason}" for s in screened))
        if not passed:
            unexplained.append(UnexplainedPeak(peak.energy_keV, peak.significance, peak.net_area,
                                               _artifact_note(peak.energy_keV, current, sigma_fn)))
            continue

        # Full fits: elements of the priority pool first, then the best candidates and their overlap partners
        priority_passed = [s for s in passed if s.element in priority]
        verified = _verify(priority_passed, True, in_model) if priority_passed else []
        from_pool = bool(verified)
        if not verified:
            others = [s for s in passed if s.element not in priority]
            if others:
                # Overlap partners of the best candidate are tested first, after it
                partners = {p_el for p_el, _, _ in _overlap_partners(others[0].element, others[0].matched_ref, sigma_fn, meas_type)}
                others = others[:1] + sorted(others[1:], key=lambda s: s.element not in partners)
            verified = _verify(others, False, in_model)

        if not verified:
            unexplained.append(UnexplainedPeak(peak.energy_keV, peak.significance, peak.net_area,
                                               _artifact_note(peak.energy_keV, current, sigma_fn)))
            continue

        verified.sort(key=lambda v: -v[0])
        d_dev, best, outcome, area, l_d = verified[0]
        ambiguous = [v[1].element for v in verified[1:] if d_dev - v[0] < settings.ambiguity_margin]
        best_partners = _overlap_partners(best.element, best.matched_ref, sigma_fn, meas_type)
        overlaps = sorted({f"{p_el} {p_line} ({sev})" for p_el, p_line, sev in best_partners
                           if p_el in in_model})
        accepted.append(IdentifiedElement(
            element=best.element,
            status='ambiguous' if ambiguous else 'identified',
            matched_line=best.matched_line,
            line_energy_keV=best.line_energy,
            peak_significance=best.family_significance.get(best.matched_ref),
            delta_deviance=None if not np.isfinite(d_dev) else float(d_dev),
            net_area=area,
            detection_limit=float(l_d),
            ambiguous_with=ambiguous,
            overlaps_with=overlaps,
            from_pool=from_pool,
        ))
        if outcome.approximate and settings.refit_after_accept and current is not None:
            # Confirm with a full fit: the linear test cannot redistribute counts between overlapping peaks
            refit = fit(model_elements())
            d_full = stat(current) - stat(refit)
            if d_full < settings.min_delta_deviance * settings.multiplier(best.element, best.line_energy):
                if verbose:
                    logger.info(f"   {best.element} rejected by full fit (ΔD = {d_full:.0f})")
                accepted.pop()
                unexplained.append(UnexplainedPeak(peak.energy_keV, peak.significance, peak.net_area,
                                                   f'{best.element} rejected by full fit'))
                continue
            accepted[-1].delta_deviance = float(d_full)
            current = refit
        else:
            current = outcome

    # Overlap check: partners of the elements in the model, tested even without residual peaks
    if settings.overlap_check and current is not None and current.background is not None:
        rejected_partners: set = set()
        while len(accepted) < settings.max_new_elements:
            if current.approximate:
                current = fit(model_elements())
            in_model = set(model_elements()) | set(substrate)
            partners: Dict[str, str] = {}
            for el in sorted(in_model):
                for ref in element_line_families(el, beam_energy, e_range):
                    for p_el, _, _ in _overlap_partners(el, ref, sigma_fn, meas_type):
                        if (p_el not in in_model and p_el not in excluded and p_el not in rejected_partners
                                and settings.z_range[0] <= Element(p_el).Z <= settings.z_range[1]):
                            partners.setdefault(p_el, f"{el} {ref}")
            if not partners:
                break
            n_tests += len(partners)
            var = _linear_variance(current, settings.model_rel_uncertainty)
            base_model = _linear_refit(y, current, var, elements=current.elements).model
            chi_base = _weighted_chi2(y, base_model, var)
            res_current = _ResidualAnalysis(y, current.model, E, F, sigma_fn, settings.model_rel_uncertainty)
            tests = []
            for p_el, overlapped in partners.items():
                offsets = _tied_offsets(p_el, current, beam_energy, e_range, sigma_fn)
                add = family_templates(p_el, offsets)
                if not add:
                    continue
                outcome = _linear_refit(y, current, var, add=add, elements=tuple(model_elements()) + (p_el,))
                d_lin = chi_base - _weighted_chi2(y, outcome.model, var)
                ref_name = max(add, key=lambda k: outcome.line_areas.get(k, 0.0))
                ref = ref_name.split('_', 1)[1]
                line_e = next(en for line, en, _ in element_line_families(p_el, beam_energy, e_range)[ref] if line == ref)
                line_e += offsets.get(ref, 0.0)
                # Confirmation lines, as in the screening: with the family carrying most of the fitted area as the
                # matched one, well-excited higher-energy families must be observed in the residual
                screen = _screen_candidate(p_el, ref, ref, line_e, res_current, beam_energy, e_range, settings, offsets)
                if screen.reason.startswith('missing confirmation'):
                    if verbose:
                        logger.info(f"Overlap check candidate {p_el} (overlaps {overlapped}): {screen.reason}")
                    continue
                area = outcome.line_areas.get(ref_name, 0.0)
                in_win = np.abs(E - line_e) <= 1.2 * FWHM_PER_SIGMA * float(sigma_fn(line_e))
                l_d = 2.71 + 3.29 * np.sqrt(max(current.model[in_win].sum(), 0.0))
                mult = settings.multiplier(p_el, line_e)
                # Linear test as pre-filter only (it underestimates the gain of overlapping elements)
                if d_lin >= settings.linear_prefilter_fraction * settings.min_delta_deviance * mult and area >= l_d * mult:
                    evidence = _peak_evidence(p_el, offsets, y, E, base_model, outcome, var, beam_energy, sigma_fn, settings)
                    if verbose:
                        logger.info(f"Overlap check candidate {p_el} (overlaps {overlapped}): linear ΔD = {d_lin:.0f}, "
                                    f"{evidence.describe()}")
                    # Partners overlap peaks of the current fit by definition, so contradictions in the linear
                    # test are not conclusive (frozen shapes); the full fit below checks them
                    if len(evidence.supported) >= settings.min_supported_peaks_overlap_check:
                        tests.append((d_lin, p_el, ref, line_e, area, l_d, overlapped, mult, offsets))
            if not tests:
                break
            d_lin, p_el, ref, line_e, area, l_d, overlapped, mult, offsets = max(tests, key=lambda t: t[0])
            refit = fit(model_elements() + [p_el])
            d_full = stat(current) - stat(refit)
            evidence = _peak_evidence(p_el, offsets, y, E, current.model, refit, var, beam_energy, sigma_fn, settings)
            if verbose:
                logger.info(f"   full fit with {p_el}: ΔD = {d_full:.0f}, {evidence.describe()}")
            if (evidence.contradicted or len(evidence.supported) < settings.min_supported_peaks_overlap_check
                    or evidence.local_fraction < settings.min_local_gain_fraction):
                d_full = -np.inf  # Not confirmed peak by peak by the full fit
            if verbose:
                logger.info(f"Overlap check: {p_el} (overlaps {overlapped}), linear ΔD = {d_lin:.0f}, "
                            f"full-fit ΔD = {d_full:.0f} {'✓' if d_full >= settings.min_delta_deviance * mult else '✗'}")
            if d_full < settings.min_delta_deviance * mult:
                rejected_partners.add(p_el)
                break
            accepted.append(IdentifiedElement(
                element=p_el, status='identified', matched_line=ref, line_energy_keV=float(line_e),
                delta_deviance=float(d_full), net_area=float(area), detection_limit=float(l_d),
                overlaps_with=[overlapped], via_overlap_check=True,
            ))
            current = refit

    # Backward pruning with full fits: fit with all identified elements, then without each of them, and
    # remove those whose removal does not worsen the fit statistic beyond the threshold. Full fits are
    # needed here, since linear refits cannot redistribute counts between overlapping peaks.
    pruned: List[IdentifiedElement] = []
    if settings.prune and accepted and current is not None:
        if current.approximate:
            current = fit(model_elements())
        # Elements with a single resolved peak are tested last: they are more easily absorbed by overlapping
        # peaks of multi-peak elements, which would make them look unnecessary
        n_peaks = {a.element: _n_resolved_peaks(a.element, beam_energy, e_range, sigma_fn, settings.min_line_weight,
                                                settings.low_energy_threshold_keV)
                   for a in accepted}
        for item in sorted(accepted, key=lambda a: (n_peaks[a.element] <= 1, a.delta_deviance or 0.0)):
            others = [el for el in model_elements() if el != item.element]
            if not others and not substrate:
                continue
            without = fit(others)
            d_dev = stat(without) - stat(current)
            item.delta_deviance = float(d_dev)
            # Peak-by-peak check with the full fits, whose line weights and shapes may differ from the templates
            # of the linear tests: a predicted peak absent in the data rejects the element
            var_without = _linear_variance(without, settings.model_rel_uncertainty)
            evidence = _peak_evidence(item.element, {}, y, E, without.model, current, var_without,
                                      beam_energy, sigma_fn, settings)
            item.filtered_gain = filtered_gain(y, without.model, current.model, F, var_without)
            item.residual_peak_sigma, item.residual_peak_ok = _residual_peak_check(
                item.element, without, y, E, F, sigma_fn, beam_energy, settings)
            if not item.residual_peak_ok:
                item.note = f"no clear residual peak at its lines ({item.residual_peak_sigma:.1f}σ)"
            # Elements also produced by the detector (e.g. Si fluorescence): kept only above a fraction of the counts
            artifact_fraction = settings.artifact_min_fraction.get(item.element)
            below_artifact = False
            if artifact_fraction is not None:
                main_area = max((a for k, a in current.line_areas.items()
                                 if k.split('_', 1)[0] == item.element and k in current.components), default=0.0)
                if main_area < artifact_fraction * float(y.sum()):
                    below_artifact = True
                    item.note = (f"area {main_area:.0f} counts below {artifact_fraction:.1%} of the total counts "
                                 f"(possible detector artifact)")
            if verbose:
                logger.info(f"Pruning check {item.element}: ΔD = {d_dev:.0f}, {evidence.describe()}, "
                            f"filtered gain {item.filtered_gain:.0f}, residual peak {item.residual_peak_sigma:.1f}σ")
            # Kept if the fit improves enough, either overall or in its peaks (filtered residual, insensitive to
            # background changes) without worsening overall
            min_gain = settings.min_delta_deviance * settings.multiplier(item.element, item.line_energy_keV or 10)
            improves = d_dev >= min_gain or (item.filtered_gain >= min_gain and d_dev > 0)
            if (not improves or evidence.contradicted or evidence.local_fraction < settings.min_local_gain_fraction
                    or below_artifact):
                item.status = 'pruned'
                accepted.remove(item)
                pruned.append(item)
                current = without

    # Final full fit with the identified elements (candidates were tested with linear fits)
    if current is not None and current.approximate and settings.final_fit:
        current = fit(model_elements())

    # Overlaps of the identified elements with the elements of the final model (overlap table)
    final_elements = set(model_elements()) | set(substrate)
    for item in accepted:
        notes = set()
        for ref, lines in element_line_families(item.element, beam_energy, e_range).items():
            ref_energy = next(en for line, en, _ in lines if line == ref)
            if ref_energy < MIN_OVERLAP_NOTE_ENERGY:
                continue
            for p_el, p_line, sev in _overlap_partners(item.element, ref, sigma_fn, meas_type):
                p_energy = _xray_lines(p_el).get(p_line, {}).get('energy (keV)', 0.0)
                if p_el in final_elements and p_el != item.element and p_energy >= MIN_OVERLAP_NOTE_ENERGY:
                    notes.add(f"{p_el} {p_line} ({sev})")
        item.overlaps_with = sorted(notes)

    # Background check of the final fit: the starting fit may be off only because an element was missing
    if not background_reliable and accepted and current is not None and not current.approximate:
        final_ratios = band_ratios(current)
        if final_ratios and not has_deficit(final_ratios):
            background_reliable = True
            background_note = background_note.replace(
                " — background unreliable, identified elements only reported as possible", "")
            background_note += f"; {_format_band_ratios(final_ratios)} with the identified elements"

    for item in accepted:
        f_gain = item.filtered_gain
        strong = (background_reliable and item.status == 'identified' and bool(item.residual_peak_ok)
                  and f_gain is not None and f_gain >= settings.strong_min_filtered_gain)
        item.tier = 'strong' if strong else 'possible'

    elements = [IdentifiedElement(el, 'fixed') for el in fixed]
    elements += [IdentifiedElement(el, 'substrate') for el in substrate]
    elements += accepted + pruned
    if has_standard is not None:
        for e in elements:
            e.has_standard = bool(has_standard(e.element))

    result = ElementIDResult(
        fixed_elements=fixed,
        substrate_elements=substrate,
        elements=elements,
        unexplained_peaks=unexplained,
        energy_vals=E,
        spectrum_vals=y,
        model=None if current is None else current.model,
        n_fits=n_fits,
        elapsed_s=time.time() - t_start,
        n_candidate_tests=n_tests,
        background_note=background_note,
    )
    if verbose:
        logger.info(result.summary())
    return result
