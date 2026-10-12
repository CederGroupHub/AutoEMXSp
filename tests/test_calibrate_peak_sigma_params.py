#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibration of the peak sigma parameters (conv_eff, elec_noise) of the detector model."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.core.fitter.detector_response import DetectorResponseFunction
from autoemx.runners.calibrate_peak_sigma_params import calibrate_peak_sigma_params

AL_STANDARD = Path(__file__).resolve().parent / "inputs" / "SDD_calibration_Al"


@pytest.fixture
def al_standard(tmp_path: Path) -> Path:
    """Copy of the fitted Al standard (loading a ledger from a new folder rewrites it)."""
    return Path(shutil.copytree(AL_STANDARD, tmp_path / "Al"))


def test_known_parameters_recovered():
    model = calibs.load_detector_model("PhenomXL")
    conv_eff, elec_noise = 0.0038, 5.2
    energies = {"B_Ka1": 0.183, "O_Ka1": 0.525, "Al_Ka1": 1.487, "Ti_Ka1": 4.511, "Cu_Ka1": 8.048}
    measured = {label: (E, DetectorResponseFunction.det_sigma(E, conv_eff, elec_noise, model.F))
                for label, E in energies.items()}

    results = calibrate_peak_sigma_params(measured_sigmas=measured, microscope_ID="PhenomXL", plot=False)

    assert results["conv_eff"] == pytest.approx(conv_eff, rel=1e-6)
    assert results["elec_noise"] == pytest.approx(elec_noise, rel=1e-6)
    assert results["current_elec_noise"] == model.elec_noise
    assert [p["label"] for p in results["points"]] == list(energies)  # sorted by energy


def test_sigmas_read_from_fitted_standard(al_standard: Path):
    results = calibrate_peak_sigma_params({str(al_standard): None}, plot=False)

    points = {p["label"].split(": ")[1]: p for p in results["points"]}
    # Reference lines only (not Ka2, Kb1, pileups), and not C_Ka1, whose sigma was fixed (small peak)
    assert sorted(points) == ["Al_Ka1", "O_Ka1"]
    assert points["Al_Ka1"]["n"] == 5
    assert points["O_Ka1"]["n"] == 2  # sigma fixed in the other 3 spectra
    assert results["microscope_ID"] == "PhenomXL"
    assert not any(p["at_fit_bound"] for p in results["points"])


def test_sigma_at_fit_bound_reported(al_standard: Path):
    sample = al_standard
    ledger_path = sample / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
    ledger = json.loads(ledger_path.read_text())
    model = calibs.load_detector_model("PhenomXL")
    peak = next(result["fit_result"]["fitted_peaks"]["Al_Ka1"] for spectrum in ledger["spectra"]
                for result in spectrum["quantification_results"] if result.get("fit_result"))
    peak["sigma"] = 1.1 * DetectorResponseFunction.det_sigma(peak["theoretical_energy"], model.conv_eff,
                                                             model.elec_noise, model.F)
    ledger_path.write_text(json.dumps(ledger))

    results = calibrate_peak_sigma_params({str(sample): ["Al_Ka1", "O_Ka1"]}, plot=False)
    assert [p["at_fit_bound"] for p in results["points"]] == [False, True]


def test_too_few_lines(al_standard: Path):
    with pytest.raises(ValueError, match="At least 2 lines"):
        calibrate_peak_sigma_params({str(al_standard): ["Al_Ka1"]}, plot=False)
