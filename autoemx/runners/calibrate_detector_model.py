#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runners: Calibrate the Peak Shape and Background Functions of the Detector Model

The detector model of a microscope (autoemx/calibrations/<microscope_ID>/detector_model.json) describes, as functions
of energy or atomic number, parameters that are fixed during fitting:

- Peak shape (incomplete charge collection, ICC): gamma, f_tail, R_e and F_loss, each a sigmoid of the X-ray energy E
  (keV): scale * (offset + amplitude / (1 + exp(-steepness * (E - center)))). gamma and f_tail apply to lines below
  1.84 keV, R_e and F_loss to lines between 1.18 and 5 keV.
- Background (modified Duncumb model): P (quadratic in the mean atomic number Z), F (1 - exp(-(Z + z_offset) /
  z_scale)) and beta (sigmoid + Gaussian in Z).

These runners refit spectra of standards with some of these parameters fitted freely (chosen by the user), collect
their fitted values, and fit the functions of the detector model to them, starting from the current coefficients and
weighting each point by its standard error. Only the coefficients chosen by the user are fitted; the others keep
their current values. Since the parameters are correlated, calibrate them one at a time, updating the detector model
before calibrating the next one (see the documentation for the recommended order).

Use high-count spectra of bulk standards, ideally pure elements, spanning the energy (peak shape) or atomic number
(background) range. The results are printed with the current values, plotted, and given as a snippet for
detector_model.json: the detector model is not changed.

Created on Mon Oct 12 2026
@author: Andrea
"""

import json
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import curve_fit

import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.runtime_configs import SampleConfig
from autoemx.core.fitter.background import Background_Model
from autoemx.core.fitter.peaks import Peaks_Model
from autoemx.data.Xray_lines import get_el_xray_lines
from autoemx.runners.fit_and_quantify_spectrum_from_ledger import fit_and_quantify_spectrum_from_ledger

SIGMOID_COEFFS = ("offset", "amplitude", "center", "steepness")


# =============================================================================
# Refitting of standards
# =============================================================================
def _refit_sample(sample_path: str, max_spectra: Optional[int], **free_params) -> Tuple[Dict[str, Any], List[Any]]:
    """Fit the spectra of a sample with the given parameters free. Returns its settings and the fitted params."""
    sample_dir = Path(sample_path).expanduser().resolve()
    ledger = load_sample_ledger(sample_dir / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}")
    configs = ledger.configs
    settings = {
        "sample": sample_dir.name,
        "microscope_ID": configs.microscope_cfg.ID,
        "meas_mode": configs.measurement_cfg.mode,
        "is_particle": configs.sample_cfg.type in SampleConfig.POWDER_SAMPLES_TYPES,
    }
    spectrum_ids = [entry.spectrum_id if entry.spectrum_id is not None else i for i, entry in enumerate(ledger.spectra)]
    fitted_params = []
    for spectrum_id in spectrum_ids[:max_spectra]:
        quantifier = fit_and_quantify_spectrum_from_ledger(
            sample_ID=ledger.sample_id, spectrum_ID=spectrum_id, ledger=ledger, results_path=str(sample_dir.parent),
            is_standard=True, quantify=False, plot_signal=False, print_results=False, quant_verbose=False,
            fitting_verbose=False, is_particle=settings["is_particle"], **free_params,
        )
        fit_result = getattr(quantifier, "fit_result", None)
        if fit_result is None:
            print(f"⚠️ Spectrum {spectrum_id} of '{sample_dir.name}' could not be fitted and is skipped.")
            continue
        fitted_params.append(fit_result.params)
    return settings, fitted_params


def _check_same_settings(all_settings: List[Dict[str, Any]]) -> Tuple[str, str]:
    instruments = {(s["microscope_ID"], s["meas_mode"]) for s in all_settings}
    if len(instruments) != 1:
        raise ValueError(f"The standards must be measured with the same microscope and mode, found {sorted(instruments)}.")
    return instruments.pop()


def _at_bound(param) -> bool:
    return any(bound is not None and np.isfinite(bound) and np.isclose(param.value, bound, rtol=1e-3, atol=0)
               for bound in (param.min, param.max))


def _summarize(values: List[float], at_bound: List[bool]) -> Dict[str, Any]:
    return {"value": float(np.mean(values)), "stdev": float(np.std(values)), "n": len(values),
            "n_at_bound": int(sum(at_bound))}


def _errors(points: List[Dict[str, Any]], scale: float = 1.0) -> np.ndarray:
    """Standard errors of the mean values of the points (0 for single values: unweighted fit)."""
    return np.array([p["stdev"] / np.sqrt(p["n"]) if p["n"] > 1 else 0.0 for p in points]) / scale


def _fit_coefficients(func: Callable, x: np.ndarray, y: np.ndarray, current: Dict[str, float],
                      coeffs_to_fit: Sequence[str], errors: Optional[np.ndarray] = None) -> Optional[Dict[str, float]]:
    """
    Fit some coefficients of func(x, **coefficients), the others keeping their current values (None: too few points).

    The fit is weighted by the inverse squared errors of the points, if all are known and positive.
    """
    coeffs_to_fit = list(coeffs_to_fit)
    if len(x) <= len(coeffs_to_fit):
        return None

    def partial(x_vals, *fitted):
        return func(x_vals, **{**current, **dict(zip(coeffs_to_fit, fitted))})

    weighted = errors is not None and bool(np.all(np.asarray(errors) > 0))
    popt, _ = curve_fit(partial, x, y, p0=[current[c] for c in coeffs_to_fit], sigma=errors if weighted else None,
                        maxfev=20000)
    return {**current, **{c: float(v) for c, v in zip(coeffs_to_fit, popt)}}


def _print_snippet(key: str, meas_mode: str, payload: Dict[str, Any], microscope_ID: str, any_fitted: bool) -> None:
    if not any_fitted:
        print("\nNo function fitted: nothing to update in the detector model.")
        return
    print(f"\nUpdated '{key}' for autoemx/calibrations/{microscope_ID}/{cnst.DETECTOR_MODEL_FILENAME}:")
    print(json.dumps({key: {meas_mode: payload}}, indent=2))


# =============================================================================
# Peak shape (ICC)
# =============================================================================
def _sigmoid(E, offset, amplitude, center, steepness):
    return offset + amplitude / (1 + np.exp(-steepness * (E - center)))


def calibrate_peak_shape_params(
    samples: Dict[str, Sequence[str]],
    params: Sequence[str] = Peaks_Model.PEAK_SHAPE_PARAMS,
    sigmoid_coeffs: Sequence[str] = SIGMOID_COEFFS,
    max_spectra_per_sample: Optional[int] = 5,
    plot: bool = True,
) -> Dict[str, Any]:
    """
    Calibrate the peak shape (ICC) functions of the detector model.

    Parameters
    ----------
    samples : dict
        {path of a standard (sample folder with ledger.json): elements whose peak shapes are fitted freely}.
    params : sequence of str, optional
        Peak shape parameters fitted freely, among 'gamma', 'f_tail', 'R_e', 'F_loss'. Default: all. The others keep
        the values of the detector model during the fits.
    sigmoid_coeffs : sequence of str, optional
        Coefficients of the sigmoids fitted to the values of each parameter, among 'offset', 'amplitude', 'center',
        'steepness'. Default: all. The others keep their current values. A sigmoid is fitted only with more points
        (lines) than coefficients.
    max_spectra_per_sample : int or None, optional
        Spectra fitted per standard (default 5; None: all).
    plot : bool, optional
        If True (default), plots the fitted values with the current and fitted functions.

    Returns
    -------
    dict
        microscope_ID, meas_mode, and for each parameter: points (line, energy, value, stdev, n, n_at_bound),
        current and fitted (None if not fitted) sigmoid coefficients.
    """
    unknown = (set(params) - set(Peaks_Model.PEAK_SHAPE_PARAMS)) | (set(sigmoid_coeffs) - set(SIGMOID_COEFFS))
    if unknown:
        raise ValueError(f"Unknown parameters or coefficients: {sorted(unknown)}.")

    values: Dict[str, Dict[str, Tuple[float, List[float], List[bool]]]] = {param: {} for param in params}
    all_settings = []
    for sample_path, els in samples.items():
        settings, fitted_params = _refit_sample(
            sample_path, max_spectra_per_sample, free_peak_shape_els=list(els), free_peak_shape_params=list(params))
        all_settings.append(settings)
        for fit_params in fitted_params:
            for name, par in fit_params.items():
                for param in params:
                    el_line = name[:-len(param) - 1]
                    if not (name.endswith(f"_{param}") and par.vary and el_line.split("_")[0] in els):
                        continue
                    el, line = el_line.split("_", 1)
                    energy = get_el_xray_lines(el)[line]["energy (keV)"]
                    entry = values[param].setdefault(f"{settings['sample']}: {el_line}", (energy, [], []))
                    entry[1].append(par.value)
                    entry[2].append(_at_bound(par))
    microscope_ID, meas_mode = _check_same_settings(all_settings)
    model = calibs.load_detector_model(microscope_ID)
    shapes = model.peak_shape[meas_mode]

    results: Dict[str, Any] = {"microscope_ID": microscope_ID, "meas_mode": meas_mode}
    for param in params:
        current_sigmoid = getattr(shapes, param)
        points = sorted(
            ({"label": label, "energy": energy, **_summarize(vals, bounds)}
             for label, (energy, vals, bounds) in values[param].items()),
            key=lambda p: p["energy"],
        )
        current = {c: getattr(current_sigmoid, c) for c in SIGMOID_COEFFS}
        fitted = None
        if points:
            E = np.array([p["energy"] for p in points])
            y = np.array([p["value"] for p in points]) / current_sigmoid.scale
            fitted = _fit_coefficients(_sigmoid, E, y, current, sigmoid_coeffs, _errors(points, current_sigmoid.scale))
        results[param] = {"points": points, "current": current, "fitted": fitted, "scale": current_sigmoid.scale}

    _print_peak_shape_results(results, params, sigmoid_coeffs, shapes)
    if plot:
        _plot_results(results, params, "Energy (keV)", "energy", (0.05, 6),
                      lambda param, r, coeffs, x: _sigmoid(x, **coeffs) * r["scale"])
    return results


def _print_peak_shape_results(results, params, sigmoid_coeffs, shapes) -> None:
    print("\n--- Peak shape (ICC) calibration ---")
    payload = shapes.model_dump()
    for param in params:
        r = results[param]
        print(f"\n{param}:")
        if not r["points"]:
            print("  No fitted values (the lines of the free elements have no such parameter at their energies).")
            continue
        for p in r["points"]:
            flag = f"  ⚠️ {p['n_at_bound']} at fit bound" if p["n_at_bound"] else ""
            print(f"  {p['label']:35s} E={p['energy']:6.3f}  {p['value']:.4g} ± {p['stdev']:.2g} (n={p['n']}){flag}")
        if r["fitted"] is None:
            print(f"  Sigmoid not fitted: {len(r['points'])} lines for {len(sigmoid_coeffs)} coefficients.")
            continue
        for c in sigmoid_coeffs:
            print(f"  {c:10s} {r['current'][c]:.6g} -> {r['fitted'][c]:.6g}")
        payload[param].update(r["fitted"])
    _print_snippet("peak_shape", results["meas_mode"], payload, results["microscope_ID"],
                   any(results[param]["fitted"] for param in params))
    print("--- End ---\n")


# =============================================================================
# Background (modified Duncumb model)
# =============================================================================
def _duncumb_P(Z, a2, a1, a0):
    return a2 * Z**2 + a1 * Z + a0


def _duncumb_F(Z, z_offset, z_scale):
    return 1 - np.exp(-(Z + z_offset) / z_scale)


def _beta(Z, sigmoid_offset, sigmoid_amplitude, sigmoid_steepness, sigmoid_center,
          gaussian_amplitude, gaussian_center, gaussian_sigma):
    return (sigmoid_offset + sigmoid_amplitude / (1 + np.exp(-sigmoid_steepness * (Z - sigmoid_center)))
            + gaussian_amplitude * np.exp(-((Z - gaussian_center) ** 2) / (2 * gaussian_sigma ** 2)))


def _background_coefficients(background) -> Dict[str, Dict[str, float]]:
    a2, a1, a0 = background.duncumb_P
    return {"P": {"a2": a2, "a1": a1, "a0": a0}, "F": background.duncumb_F.model_dump(),
            "beta": background.beta.model_dump()}


_BACKGROUND_FUNCTIONS = {"P": _duncumb_P, "F": _duncumb_F, "beta": _beta}


def calibrate_background_params(
    samples: Sequence[str],
    params: Sequence[str] = Background_Model.BACKGROUND_PARAMS,
    fit_functions: Optional[Sequence[str]] = None,
    coeffs: Optional[Dict[str, Sequence[str]]] = None,
    max_spectra_per_sample: Optional[int] = 5,
    plot: bool = True,
) -> Dict[str, Any]:
    """
    Calibrate the background functions (modified Duncumb model) of the detector model.

    Parameters
    ----------
    samples : sequence of str
        Paths of standards (sample folders with ledger.json), ideally pure elements: the values are attributed to the
        mean atomic number Z of each standard (for beta, exact only for pure elements).
    params : sequence of str, optional
        Background parameters fitted freely, among 'P', 'F', 'beta'. Default: all. The others are set by the
        functions of the detector model during the fits (also beta, otherwise fitted freely for bulk samples).
    fit_functions : sequence of str, optional
        Parameters whose functions are fitted to their values, among params. Default: all params. E.g. params
        ['P', 'F', 'beta'] with fit_functions ['P'] fits the spectra with all three free, and only the function of P.
    coeffs : dict, optional
        Coefficients fitted for each parameter, e.g. {'P': ['a0']} (P: 'a2', 'a1', 'a0', for a2*Z**2 + a1*Z + a0;
        F: 'z_offset', 'z_scale'; beta: the fields of BetaParams). Default: all. The others keep their current
        values. A function is fitted only with more points (standards) than coefficients.
    max_spectra_per_sample : int or None, optional
        Spectra fitted per standard (default 5; None: all).
    plot : bool, optional
        If True (default), plots the fitted values with the current and fitted functions.

    Returns
    -------
    dict
        microscope_ID, meas_mode, and for each parameter: points (sample, Z, value, stdev, n, n_at_bound),
        current and fitted (None if not fitted) coefficients.
    """
    unknown = set(params) - set(Background_Model.BACKGROUND_PARAMS)
    if unknown:
        raise ValueError(f"Unknown background parameters {sorted(unknown)}. Available: {Background_Model.BACKGROUND_PARAMS}")
    fit_functions = list(params if fit_functions is None else fit_functions)
    if set(fit_functions) - set(params):
        raise ValueError(f"fit_functions {sorted(set(fit_functions) - set(params))} must be among the free params {list(params)}.")

    points: Dict[str, List[Dict[str, Any]]] = {param: [] for param in params}
    all_settings = []
    for sample_path in samples:
        settings, fitted_params = _refit_sample(sample_path, max_spectra_per_sample,
                                                free_background_params=list(params))
        all_settings.append(settings)
        if not fitted_params:
            continue
        Z = float(np.mean([fit_params["Z"].value for fit_params in fitted_params]))
        for param in params:
            free = [fit_params[param] for fit_params in fitted_params if fit_params[param].vary]
            if free:
                points[param].append({"label": settings["sample"], "Z": Z,
                                      **_summarize([par.value for par in free], [_at_bound(par) for par in free])})
    microscope_ID, meas_mode = _check_same_settings(all_settings)
    model = calibs.load_detector_model(microscope_ID)
    current_coeffs = _background_coefficients(model.background[meas_mode])

    results: Dict[str, Any] = {"microscope_ID": microscope_ID, "meas_mode": meas_mode}
    for param in params:
        param_points = sorted(points[param], key=lambda p: p["Z"])
        coeffs_to_fit = (coeffs or {}).get(param, list(current_coeffs[param]))
        unknown = set(coeffs_to_fit) - set(current_coeffs[param])
        if unknown:
            raise ValueError(f"Unknown coefficients of {param}: {sorted(unknown)}. Available: {list(current_coeffs[param])}")
        fitted = None
        if param_points and param in fit_functions:
            fitted = _fit_coefficients(_BACKGROUND_FUNCTIONS[param], np.array([p["Z"] for p in param_points]),
                                       np.array([p["value"] for p in param_points]), current_coeffs[param],
                                       coeffs_to_fit, _errors(param_points))
        results[param] = {"points": param_points, "current": current_coeffs[param], "fitted": fitted,
                          "coeffs": list(coeffs_to_fit), "fit_function": param in fit_functions}

    _print_background_results(results, params, model.background[meas_mode])
    if plot:
        _plot_results(results, params, "Mean atomic number Z", "Z", (1, 92),
                      lambda param, r, coeffs, x: _BACKGROUND_FUNCTIONS[param](x, **coeffs))
    return results


def _print_background_results(results, params, background) -> None:
    print("\n--- Background calibration ---")
    payload = background.model_dump()
    for param in params:
        r = results[param]
        print(f"\n{param}:")
        if not r["points"]:
            print("  No fitted values.")
            continue
        for p in r["points"]:
            flag = f"  ⚠️ {p['n_at_bound']} at fit bound" if p["n_at_bound"] else ""
            print(f"  {p['label']:25s} Z={p['Z']:6.2f}  {p['value']:.4g} ± {p['stdev']:.2g} (n={p['n']}){flag}")
        if not r["fit_function"]:
            print("  Function not fitted (not in fit_functions).")
            continue
        if r["fitted"] is None:
            print(f"  Function not fitted: {len(r['points'])} standards for {len(r['coeffs'])} coefficients.")
            continue
        for c in r["coeffs"]:
            print(f"  {c:20s} {r['current'][c]:.6g} -> {r['fitted'][c]:.6g}")
        if param == "P":
            payload["duncumb_P"] = [r["fitted"]["a2"], r["fitted"]["a1"], r["fitted"]["a0"]]
        elif param == "F":
            payload["duncumb_F"] = r["fitted"]
        else:
            payload["beta"] = r["fitted"]
    _print_snippet("background", results["meas_mode"], payload, results["microscope_ID"],
                   any(results[param]["fitted"] for param in params))
    print("--- End ---\n")


# =============================================================================
# Plot
# =============================================================================
def _plot_results(results, params, x_label, x_key, x_range, curve) -> None:
    """Fitted values of each parameter, with its current and fitted function curve(param, result, coeffs, x)."""
    import matplotlib.pyplot as plt

    params = [param for param in params if results[param]["points"]]
    if not params:
        return
    fig, axes = plt.subplots(1, len(params), figsize=(5 * len(params), 4), squeeze=False)
    x_curve = np.linspace(*x_range, 300)
    for ax, param in zip(axes[0], params):
        r = results[param]
        ax.errorbar([p[x_key] for p in r["points"]], [p["value"] for p in r["points"]],
                    yerr=[p["stdev"] for p in r["points"]], fmt="o", label="Fitted values")
        for coeffs, style, label in ((r["current"], "k--", "Current"), (r["fitted"], "r-", "Fitted")):
            if coeffs is not None:
                ax.plot(x_curve, curve(param, r, coeffs, x_curve), style, label=label)
        ax.set_xlabel(x_label)
        ax.set_ylabel(param)
        ax.grid(True)
        ax.legend()
    fig.suptitle(f"{results['microscope_ID']} ({results['meas_mode']})")
    fig.tight_layout()
    plt.show()
