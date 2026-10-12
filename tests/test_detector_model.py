#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Detector model (detector_model.json, in the package) and detector response cache (in the local AutoEMX folder).

The local AutoEMX folder is a temporary folder in tests (autoemx_home fixture in conftest.py).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import autoemx._paths as paths
import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.config.schema_models.detector_model import DetectorModel
from autoemx.core.fitter.detector_response import DetectorResponseFunction

MICROSCOPE_ID = "PhenomXL"


def _package_model() -> DetectorModel:
    return DetectorModel.from_json_file(paths.package_calibration_dir(MICROSCOPE_ID) / cnst.DETECTOR_MODEL_FILENAME)


def test_model_values_set_as_calibration_attributes():
    calibs.load_microscope_calibrations(MICROSCOPE_ID, "point")
    model = calibs.detector_model
    assert calibs.elec_noise == model.elec_noise
    assert calibs.available_meas_modes == ("point", "map")
    assert calibs.escape_peak_probability == {"point": 0.03}
    gamma, f_tail, R_e, F_loss = calibs.get_calibrated_peak_shape_params(1.4867, "point")
    assert gamma == pytest.approx(1 + 3.5 / (1 + np.exp(-9.63 * (1.4867 - 0.88))))
    assert R_e == pytest.approx((80.67 + 101.43 / (1 + np.exp(-8.6 * (1.4867 - 1.9051)))) * 1e-7)
    P, F, beta = calibs.get_calibrated_background_params("point")
    assert P == "2.041e-5*Z**2 - 0.004076*Z + 1.252"
    assert float(beta.subs("Z", 28)) == pytest.approx(0.0419 + 0.2105 / (1 + np.exp(-0.4788 * (28 - 45.9098))) + 0.2807)


def test_model_read_from_package(autoemx_home: Path):
    calibs.load_microscope_calibrations(MICROSCOPE_ID)
    assert calibs.detector_model == _package_model()
    assert not list(autoemx_home.rglob(cnst.DETECTOR_MODEL_FILENAME))


def test_invalid_mode_and_unknown_microscope(autoemx_home: Path):
    with pytest.raises(ValueError, match="Available modes"):
        calibs.load_microscope_calibrations(MICROSCOPE_ID, "line_scan")
    with pytest.raises(ValueError, match="Could not find calibrations"):
        calibs.load_microscope_calibrations("NoSuchMicroscope")
    assert not (autoemx_home / "calibrations" / "NoSuchMicroscope").exists()


def test_model_of_newer_autoemx_rejected():
    payload = _package_model().model_dump(mode="json")
    payload["schema_version"] = 99
    with pytest.raises(ValueError, match="update AutoEMX"):
        DetectorModel.model_validate(payload)


def test_mode_without_calibrated_functions():
    calibs.load_microscope_calibrations(MICROSCOPE_ID)
    with pytest.raises(ValueError, match="Missing peak shape parameters for meas_mode 'map'"):
        calibs.get_calibrated_peak_shape_params(1.0, "map")


def test_fingerprint_changes_with_any_value():
    model = calibs.load_detector_model(MICROSCOPE_ID)
    changed = model.model_copy(deep=True)
    changed.peak_shape["point"].R_e.scale = 2e-7
    assert model.fingerprint() == calibs.load_detector_model(MICROSCOPE_ID).fingerprint()
    assert changed.fingerprint() != model.fingerprint()


def test_detector_response_cached_per_model(monkeypatch: pytest.MonkeyPatch):
    calibs.load_microscope_calibrations(MICROSCOPE_ID)
    offset, width = -0.0324, 0.00994  # channel setup used only here, so not cached by other tests
    n_computed = []
    compute = DetectorResponseFunction._compute_and_cache_conv_matrices.__func__

    def counting_compute(cls, **kwargs):
        n_computed.append(1)
        return compute(cls, **kwargs)

    monkeypatch.setattr(DetectorResponseFunction, "_compute_and_cache_conv_matrices", classmethod(counting_compute))
    cache_dir = paths.user_cache_dir(MICROSCOPE_ID) / "detector_response"
    n_cached_before = len(list(cache_dir.glob("*.npz")))

    assert DetectorResponseFunction.ensure_conv_matrices_cached(offset, width, MICROSCOPE_ID, verbose=False) is True
    DetectorResponseFunction.setup_detector_response_vars(offset, width, (14, 1100), MICROSCOPE_ID, verbose=False)
    computed_op = DetectorResponseFunction.det_response_op.toarray()
    assert len(n_computed) == 1 and len(list(cache_dir.glob("*.npz"))) == n_cached_before + 1

    # Detector model changed: recomputed, in a separate file
    original_model = calibs.detector_model
    monkeypatch.setattr(calibs, "detector_model", original_model.model_copy(update={"elec_noise": 5.5}))
    assert DetectorResponseFunction.ensure_conv_matrices_cached(offset, width, MICROSCOPE_ID, verbose=False) is True
    assert len(n_computed) == 2 and len(list(cache_dir.glob("*.npz"))) == n_cached_before + 2

    # Back to the original model: loaded from the cache, identical to the computed matrices
    monkeypatch.setattr(calibs, "detector_model", original_model)
    DetectorResponseFunction.setup_detector_response_vars(offset, width, (14, 1100), MICROSCOPE_ID, verbose=False)
    assert len(n_computed) == 2
    assert np.array_equal(DetectorResponseFunction.det_response_op.toarray(), computed_op)
