#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quantification flags (``quant_flag``) of spectra: whether a fit or quantification is reliable and, if not, why.

Pre-fit checks (flags 1-3) are made on the measured spectrum; post-fit checks (flags 0, -1, 4-9) combine the problems
reported by the quantifier (``XSp_Quantifier.bad_quant_flag``) with checks of the fitted peaks. See the documentation
page "Quantification flags" for the meaning of each flag.

==== ================================================================================================
Flag Meaning
==== ================================================================================================
0    Quantification OK
-1   As 0, but the quantification did not converge
1    No spectral data
2    Total counts too low
3    Background counts too low (pre-fit)
4    Poor fit
5    Excessively high analytical error
6    Excessive X-ray absorption
7    Substrate peaks too intense
8    Background counts under the reference peaks too low (re-evaluated at analysis)
9    Fit interrupted for unknown reasons
==== ================================================================================================
"""

from typing import Any, Optional, Sequence

import numpy as np

import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
from autoemx.config.runtime_configs import QuantificationOptionsConfig


def format_min_total_counts_percent(fraction: float) -> str:
    """Format a total-counts fraction as a compact percentage string."""
    return f"{float(fraction) * 100:.4g}%"


def detect_prefit_spectrum_issues(
    *,
    spectrum: Optional[np.ndarray],
    energy_vals: Sequence[float],
    sp_start: int,
    sp_end: int,
    target_acquisition_counts: int,
    min_bckgrnd_cnts: Optional[float],
    min_total_counts_fraction: float = QuantificationOptionsConfig.DEFAULT_MIN_TOTAL_COUNTS_FRACTION,
) -> tuple[Optional[int], Optional[str]]:
    """Detect pre-fit spectrum quality issues that map to quant flags 1-3."""
    if spectrum is None or len(spectrum) == 0:
        return 1, "No spectral data present"

    if np.sum(spectrum) < min_total_counts_fraction * target_acquisition_counts:
        return 2, "Total counts too low"

    n_vals_considered = dflt.prefit_bckgrnd_n_points
    filter_len = dflt.prefit_bckgrnd_smoothing_channels
    en_threshold = dflt.prefit_bckgrnd_max_energy_keV

    xy_data = zip(energy_vals, spectrum[sp_start:sp_end])
    spectrum_data_to_consider = [cnts for en, cnts in xy_data if cnts > 0 and en < en_threshold]
    if len(spectrum_data_to_consider) == 0:
        return 3, "Background counts too low"

    spectrum_smooth = np.convolve(
        spectrum_data_to_consider,
        np.ones(filter_len) / filter_len,
        mode='same',
    )
    min_vals = np.sort(spectrum_smooth)[:n_vals_considered]
    if min_bckgrnd_cnts is not None and len(min_vals) > 0 and all(min_vals < min_bckgrnd_cnts):
        return 3, "Background counts too low"

    return None, None


def should_interrupt_prefit_issue(
    prefit_quant_flag: Optional[int],
    interrupt_fits_bad_spectra: bool,
) -> bool:
    """Return whether a detected pre-fit issue should skip fitting entirely."""
    if prefit_quant_flag is None:
        return False
    if prefit_quant_flag == 1:
        return True
    if prefit_quant_flag in (2, 3):
        return interrupt_fits_bad_spectra
    return False


def merge_prefit_quant_flag(
    quant_flag: int,
    comment: str,
    prefit_quant_flag: Optional[int],
    prefit_comment: Optional[str],
) -> tuple[int, str]:
    """Apply pre-fit flags 2/3 when fitting proceeded without interruption."""
    if prefit_quant_flag in (2, 3) and quant_flag == 0:
        return prefit_quant_flag, prefit_comment or comment
    return quant_flag, comment


def flag_spectrum_for_clustering(
    *,
    min_bckgrnd_ref_lines: float,
    quantifier: Any,
    detectable_els_substrate: Sequence[str],
    min_bckgrnd_cnts: Optional[float],
) -> tuple[bool, int, str]:
    """Flag 7 (substrate peaks too intense) or 8 (background under reference peaks too low) of a quantified spectrum."""
    sub_peak_int_threshold = dflt.max_substrate_peak_percent_total_counts
    sub_peak_int_thresh_cnts = quantifier.tot_sp_counts * sub_peak_int_threshold / 100

    els_substrate_intensities = {el: 0 for el in detectable_els_substrate}
    for el_line, peak_info in quantifier.fitted_peaks_info.items():
        el = el_line.split('_')[0]
        if el in detectable_els_substrate:
            els_substrate_intensities[el] += peak_info[cnst.PEAK_INTENSITY_KEY]

    for el, peak_int in els_substrate_intensities.items():
        if peak_int > sub_peak_int_thresh_cnts:
            return False, 7, f"{el} {peak_int:.0f} counts > {sub_peak_int_threshold} % of total counts"

    if min_bckgrnd_cnts is not None and min_bckgrnd_ref_lines < min_bckgrnd_cnts:
        return False, 8, f"Reference background counts too low ({min_bckgrnd_ref_lines:.1f} < {min_bckgrnd_cnts})"

    return True, 0, ""


def check_fit_quant_validity(
    *,
    is_quant_fit_valid: bool,
    bad_quant_flag: Optional[int],
    quantifier: Any,
    min_bckgrnd_ref_lines: Any,
    detectable_els_substrate: Sequence[str],
    min_bckgrnd_cnts: Optional[float],
) -> tuple[int, str]:
    """Quantification flag and comment of a spectrum after fitting and quantification (flags 0, -1, 4-9)."""
    start_str_comments = 'Fit interrupted due to ' if not is_quant_fit_valid else ''

    if bad_quant_flag == 1:
        comment = start_str_comments + 'poor fit'
        quant_flag = 4
    elif bad_quant_flag == 2:
        comment = start_str_comments + 'excessively high analytical error'
        quant_flag = 5
    elif bad_quant_flag == 3:
        comment = start_str_comments + 'excessive X-ray absorption'
        quant_flag = 6
    elif not is_quant_fit_valid:
        comment = 'Fit interrupted for unknown reasons'
        quant_flag = 9
    else:
        _, quant_flag, comment = flag_spectrum_for_clustering(
            min_bckgrnd_ref_lines=min_bckgrnd_ref_lines,
            quantifier=quantifier,
            detectable_els_substrate=detectable_els_substrate,
            min_bckgrnd_cnts=min_bckgrnd_cnts,
        )

    if bad_quant_flag == -1 and quant_flag == 0:
        comment = f"{comment} - Quantification did not converge." if comment else 'Quantification did not converge.'
        quant_flag = -1

    return quant_flag, comment
