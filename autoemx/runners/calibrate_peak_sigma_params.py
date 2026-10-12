#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runner: Calibrate the Peak Width (Sigma) Parameters of the Detector Model

The Gaussian sigma of the X-ray peaks depends on the energy E (keV) as

    sigma = conv_eff * sqrt(elec_noise**2 + E * F / conv_eff)

(see DetectorResponseFunction.det_sigma), with conv_eff the energy per electron-hole pair (keV), elec_noise the
electronic noise (electrons) and F the Fano factor (fixed, from the detector model). Equivalently,
sigma**2 = (conv_eff * elec_noise)**2 + conv_eff * F * E, the form of
J. Osán et al., X-Ray Spectrom. 30 (2001) 419–426. https://doi.org/10.1002/xrs.523.

This runner fits conv_eff and elec_noise to the peak sigmas fitted in spectra of standards (read from their ledgers)
and/or to sigmas given directly, and compares them with the values of the detector model of the microscope.

Use high-count spectra of bulk standards, ideally pure elements, with well-separated lines spanning the energy range.
Only reference lines are used (Ka1, La1, Ma1 by default): the sigmas of the other lines of an element are tied to
them, and small peaks have their sigma fixed during fitting. Since fitted sigmas are bounded to 0.93-1.1 times the
value of the detector model, sigmas at these bounds are reported: update the detector model, refit the spectra, and
calibrate again.

Created on Mon Oct 12 2026
@author: Andrea
"""

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import curve_fit

import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import get_active_quantification_config, load_sample_ledger
from autoemx.core.fitter.detector_response import DetectorResponseFunction

DEFAULT_REFERENCE_LINES = ("Ka1", "La1", "Ma1")
# Bounds of the fitted sigmas, relative to the value of the detector model (Peaks_Model._get_gaussian_sigma_param)
SIGMA_BOUNDS = (0.93, 1.1)


def _collect_sigmas(sample_path: Path, el_lines: Optional[List[str]]) -> Tuple[str, Dict[str, Tuple[float, List[float]]]]:
    """Fitted sigmas of the lines of a sample, from its active quantification (accepted spectra only)."""
    ledger = load_sample_ledger(sample_path / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
    quant_config = get_active_quantification_config(ledger)
    microscope_ID = ledger.configs.microscope_cfg.ID
    sigmas: Dict[str, Tuple[float, List[float]]] = {}
    if quant_config is None:
        return microscope_ID, sigmas
    line_pattern = re.compile(rf"^[A-Z][a-z]?_({'|'.join(DEFAULT_REFERENCE_LINES)})$")
    for spectrum in ledger.spectra:
        for result in spectrum.quantification_results:
            if result.quantification_id != quant_config.quantification_id or result.quant_flag != 0:
                continue
            if result.fit_result is None:
                continue
            for el_line, peak in result.fit_result.fitted_peaks.items():
                selected = el_line in el_lines if el_lines else bool(line_pattern.match(el_line))
                if selected and peak.sigma and peak.theoretical_energy:
                    sigmas.setdefault(el_line, (peak.theoretical_energy, []))[1].append(peak.sigma)
    return microscope_ID, sigmas


def calibrate_peak_sigma_params(
    samples: Optional[Dict[str, Optional[List[str]]]] = None,
    measured_sigmas: Optional[Dict[str, Tuple[float, float]]] = None,
    microscope_ID: Optional[str] = None,
    plot: bool = True,
) -> Dict[str, Any]:
    """
    Fit conv_eff and elec_noise of the detector model to measured peak sigmas.

    Parameters
    ----------
    samples : dict, optional
        {path of a fitted sample folder (with ledger.json): el_lines to use, or None for its Ka1, La1, Ma1 lines}.
        Sigmas are read from the active quantification of each sample, for accepted spectra (quant_flag 0).
    measured_sigmas : dict, optional
        Sigmas measured otherwise, as {label: (energy in keV, sigma in keV)}.
    microscope_ID : str, optional
        Microscope whose detector model is compared (and whose Fano factor is used). Default: that of the samples.
    plot : bool, optional
        If True (default), plots the measured sigmas with the fitted and current model.

    Returns
    -------
    dict
        conv_eff, elec_noise (fitted), current_conv_eff, current_elec_noise, F, microscope_ID, and points:
        a list of {label, energy, sigma, stdev, n, current_ratio, at_fit_bound}.
    """
    points: List[Dict[str, Any]] = []
    sample_microscopes = set()
    for sample_path, el_lines in (samples or {}).items():
        sample_microscope, sigmas = _collect_sigmas(Path(sample_path).expanduser(), el_lines)
        sample_microscopes.add(sample_microscope)
        if not sigmas:
            print(f"⚠️ No fitted sigmas found in '{sample_path}' (requested lines: {el_lines or 'Ka1, La1, Ma1'}).")
        for el_line, (energy, values) in sigmas.items():
            points.append({"label": f"{Path(sample_path).name}: {el_line}", "energy": energy, "values": values,
                           "fitted": True})
    for label, (energy, sigma) in (measured_sigmas or {}).items():
        points.append({"label": label, "energy": energy, "values": [sigma], "fitted": False})

    if microscope_ID is None:
        if len(sample_microscopes) != 1:
            raise ValueError(f"Set microscope_ID: the samples were measured with {sorted(sample_microscopes) or 'none'}.")
        microscope_ID = sample_microscopes.pop()
    model = calibs.load_detector_model(microscope_ID)
    F = model.F

    def current_sigma(E):
        return DetectorResponseFunction.det_sigma(E, model.conv_eff, model.elec_noise, F)

    # Average by line; drop fitted sigmas fixed at the model value (small peaks)
    used = []
    for point in points:
        values = np.asarray(point.pop("values"), dtype=float)
        ratios = values / current_sigma(point["energy"])
        fitted = point.pop("fitted")
        if fitted:
            values, ratios = values[~np.isclose(ratios, 1.0, rtol=0, atol=1e-9)], ratios[~np.isclose(ratios, 1.0, rtol=0, atol=1e-9)]
            if values.size == 0:
                continue
        point.update(sigma=float(values.mean()), stdev=float(values.std()), n=int(values.size))
        point["current_ratio"] = point["sigma"] / current_sigma(point["energy"])
        point["at_fit_bound"] = fitted and bool(np.any(np.isclose(ratios[:, None], SIGMA_BOUNDS, atol=1e-3)))
        used.append(point)
    points = sorted(used, key=lambda p: p["energy"])
    if len(points) < 2:
        raise ValueError(f"At least 2 lines with fitted sigmas are needed to calibrate, found {len(points)}.")

    energies = np.array([p["energy"] for p in points])
    sigmas = np.array([p["sigma"] for p in points])
    errors = np.array([p["stdev"] / np.sqrt(p["n"]) if p["n"] > 1 and p["stdev"] > 0 else np.nan for p in points])
    if np.all(np.isfinite(errors)):
        sigma_weights, absolute = errors, True
    else:
        sigma_weights, absolute = None, False

    popt, pcov = curve_fit(
        lambda E, conv_eff, elec_noise: DetectorResponseFunction.det_sigma(E, conv_eff, elec_noise, F),
        energies, sigmas, p0=[model.conv_eff, model.elec_noise], sigma=sigma_weights, absolute_sigma=absolute,
        bounds=([1e-6, 0], [np.inf, np.inf]),
    )
    conv_eff, elec_noise = (float(v) for v in popt)
    conv_eff_err, elec_noise_err = (float(v) for v in np.sqrt(np.diag(pcov)))

    results = {
        "microscope_ID": microscope_ID, "F": F,
        "conv_eff": conv_eff, "elec_noise": elec_noise, "conv_eff_err": conv_eff_err, "elec_noise_err": elec_noise_err,
        "current_conv_eff": model.conv_eff, "current_elec_noise": model.elec_noise, "points": points,
    }
    _print_results(results)
    if plot:
        _plot_results(results, current_sigma)
    return results


def _print_results(results: Dict[str, Any]) -> None:
    print("\n--- Peak sigma calibration ---")
    print(f"{'Line':40s} {'E (keV)':>8s} {'sigma (eV)':>11s} {'stdev':>7s} {'n':>3s} {'/current':>9s}")
    for p in results["points"]:
        flag = "  ⚠️ at fit bound" if p["at_fit_bound"] else ""
        print(f"{p['label']:40s} {p['energy']:8.3f} {p['sigma'] * 1000:11.2f} {p['stdev'] * 1000:7.2f} {p['n']:3d} "
              f"{p['current_ratio']:9.3f}{flag}")
    print(f"\nFano factor (fixed, from the detector model): {results['F']}")
    print(f"conv_eff   : {results['conv_eff']:.6f} ± {results['conv_eff_err']:.6f} keV "
          f"(current: {results['current_conv_eff']})")
    print(f"elec_noise : {results['elec_noise']:.3f} ± {results['elec_noise_err']:.3f} electrons "
          f"(current: {results['current_elec_noise']})")
    if any(p["at_fit_bound"] for p in results["points"]):
        print("⚠️ Some sigmas are at the bounds of the fit (0.93-1.1 times the current model): update the detector "
              "model, refit the spectra and calibrate again.")
    print(f"To use these values, update 'conv_eff' and 'elec_noise' in "
          f"autoemx/calibrations/{results['microscope_ID']}/{cnst.DETECTOR_MODEL_FILENAME}.")
    print("--- End ---\n")


def _plot_results(results: Dict[str, Any], current_sigma) -> None:
    import matplotlib.pyplot as plt

    points = results["points"]
    energies = np.array([p["energy"] for p in points])
    E = np.linspace(0, max(energies) * 1.1, 200)
    fitted = DetectorResponseFunction.det_sigma(E, results["conv_eff"], results["elec_noise"], results["F"])
    plt.figure(figsize=(8, 6))
    plt.errorbar(energies, [p["sigma"] * 1000 for p in points], yerr=[p["stdev"] * 1000 for p in points],
                 fmt="o", label="Measured")
    for p in points:
        plt.annotate(p["label"].split(": ")[-1], (p["energy"], p["sigma"] * 1000), textcoords="offset points",
                     xytext=(4, 4), fontsize=8)
    plt.plot(E, fitted * 1000, "r-", label="Fitted")
    plt.plot(E, current_sigma(E) * 1000, "k--", label="Current detector model")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Peak sigma (eV)")
    plt.title(f"Peak sigma calibration – {results['microscope_ID']}")
    plt.legend()
    plt.grid(True)
    plt.show()
