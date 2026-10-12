#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Peak shape and background parameters fitted freely, and calibration of the detector model functions.

Fits use one or two spectra of a copy of the fitted Al standard (loading a ledger from a new folder rewrites it).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from autoemx.core.fitter.background import Background_Model
from autoemx.core.fitter.peaks import Peaks_Model
from autoemx.runners import calibrate_detector_model as cdm
from autoemx.runners.fit_and_quantify_spectrum_from_ledger import fit_and_quantify_spectrum_from_ledger

AL_STANDARD = Path(__file__).resolve().parent / "inputs" / "SDD_calibration_Al"


@pytest.fixture
def al_standard(tmp_path: Path) -> Path:
    return Path(shutil.copytree(AL_STANDARD, tmp_path / "Al"))


def _fit_params(sample: Path, **free_params):
    quantifier = fit_and_quantify_spectrum_from_ledger(
        sample_ID="Al", spectrum_ID=0, results_path=str(sample.parent), is_standard=True, quantify=False,
        plot_signal=False, print_results=False, quant_verbose=False, fitting_verbose=False, is_particle=False,
        **free_params,
    )
    return quantifier.fit_result.params


def test_default_fit_has_fixed_peak_shapes_and_background_functions(al_standard: Path):
    params = _fit_params(al_standard)
    assert not any(params[f"Al_Ka1_{p}"].vary for p in Peaks_Model.PEAK_SHAPE_PARAMS)
    assert params["P"].expr and params["F"].expr
    assert params["beta"].vary  # always free for bulk samples


def test_chosen_peak_shape_params_fitted_for_chosen_elements(al_standard: Path):
    params = _fit_params(al_standard, free_peak_shape_els=["Al"], free_peak_shape_params=["R_e"])
    assert params["Al_Ka1_R_e"].vary
    assert not params["Al_Ka1_F_loss"].vary and not params["Al_Ka1_gamma"].vary
    # Other lines of Al: free parameters equal to those of the reference line, the others at their own values
    assert params["Al_Ka2_R_e"].expr == "Al_Ka1_R_e"
    assert params["Al_Ka2_F_loss"].expr is None and not params["Al_Ka2_F_loss"].vary
    # Other elements unchanged
    assert not params["O_Ka1_gamma"].vary


def test_all_peak_shape_params_fitted_by_default(al_standard: Path):
    params = _fit_params(al_standard, free_peak_shape_els=["Al"])
    assert all(params[f"Al_Ka1_{p}"].vary for p in Peaks_Model.PEAK_SHAPE_PARAMS)
    assert all(params[f"Al_Ka2_{p}"].expr == f"Al_Ka1_{p}" for p in Peaks_Model.PEAK_SHAPE_PARAMS)


def test_chosen_background_params_fitted(al_standard: Path):
    params = _fit_params(al_standard, free_background_params=["P"])
    assert params["P"].vary and params["P"].expr is None
    assert params["F"].expr and not params["F"].vary
    assert params["beta"].expr and not params["beta"].vary  # set by the detector model, also for bulk samples

    params = _fit_params(al_standard, free_background_params=["F", "beta"])
    assert params["P"].expr and params["F"].vary and params["beta"].vary


def test_unknown_params_rejected():
    with pytest.raises(ValueError, match="Unknown background parameters"):
        Background_Model(is_particle=False, free_background_params=["Q"])
    with pytest.raises(ValueError, match="Unknown parameters or coefficients"):
        cdm.calibrate_peak_shape_params({}, params=["sigma"])


def test_only_chosen_coefficients_fitted():
    current = {"offset": 1.0, "amplitude": 3.5, "center": 0.88, "steepness": 9.63}
    E = np.linspace(0.1, 3, 10)
    y = cdm._sigmoid(E, **{**current, "offset": 2.0, "center": 1.2})

    fitted = cdm._fit_coefficients(cdm._sigmoid, E, y, current, ["offset", "center"])
    assert fitted["offset"] == pytest.approx(2.0) and fitted["center"] == pytest.approx(1.2)
    assert (fitted["amplitude"], fitted["steepness"]) == (3.5, 9.63)
    assert cdm._fit_coefficients(cdm._sigmoid, E[:2], y[:2], current, ["offset", "center"]) is None  # too few points

    # Weighted by the errors of the points: an outlier with a large error barely moves the fit
    y_outlier, errors = y.copy(), np.full(len(y), 1e-3)
    y_outlier[0], errors[0] = y[0] + 1.0, 10.0
    fitted = cdm._fit_coefficients(cdm._sigmoid, E, y_outlier, current, ["offset", "center"], errors)
    assert fitted["offset"] == pytest.approx(2.0, abs=1e-3)

    Z = np.array([6.0, 13, 29, 47, 79])
    P_current = {"a2": 2.041e-5, "a1": -0.004076, "a0": 1.252}
    P_fitted = cdm._fit_coefficients(cdm._duncumb_P, Z, cdm._duncumb_P(Z, **{**P_current, "a0": 1.3}), P_current, ["a0"])
    assert P_fitted["a0"] == pytest.approx(1.3)


def test_peak_shape_calibration_runs(al_standard: Path):
    results = cdm.calibrate_peak_shape_params({str(al_standard): ["Al"]}, params=["R_e", "gamma"],
                                              max_spectra_per_sample=1, plot=False)
    assert (results["microscope_ID"], results["meas_mode"]) == ("PhenomXL", "point")
    R_e_points = results["R_e"]["points"]
    assert [p["label"] for p in R_e_points] == ["Al: Al_Ka1"]
    assert R_e_points[0]["energy"] == pytest.approx(1.4867, abs=1e-3)
    assert results["R_e"]["fitted"] is None  # 1 line for 4 coefficients
    assert results["R_e"]["scale"] == 1e-7


def test_background_calibration_runs(al_standard: Path):
    results = cdm.calibrate_background_params([str(al_standard)], params=["P", "F"], max_spectra_per_sample=1,
                                              plot=False)
    [P_point] = results["P"]["points"]
    assert P_point["Z"] == pytest.approx(13) and P_point["n"] == 1
    assert results["P"]["fitted"] is None  # 1 standard for 3 coefficients
    assert "beta" not in results

    with pytest.raises(ValueError, match="fit_functions"):
        cdm.calibrate_background_params([str(al_standard)], params=["P"], fit_functions=["F"], plot=False)
