#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for the AutoEMX web GUI helpers (no full quantification)."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

from autoemx.utils import load_msa
from autoemx.web.exports import (
    figure_to_png_bytes,
    fitted_spectrum_figure,
    format_composition_txt,
)
from autoemx.web.pipeline import (
    SUPPORTED_UPLOAD_EXTENSIONS,
    SpectrumFitResult,
    beam_energy_supports_quantification,
    load_uploaded_spectrum,
    parse_elements,
    parse_emsa_geometry,
)
from autoemx.web.reader_report import (
    REPORT_EMAIL,
    SpectrumReadError,
    mailto_reader_report_url,
)
from ci_paths import TESTS_DIR


def _example_msa():
    return (
        TESTS_DIR.parent
        / "autoemx"
        / "scripts"
        / "input"
        / "Example_spectrum.msa"
    )


def test_quantification_requires_15_kv():
    assert beam_energy_supports_quantification(15.0)
    assert beam_energy_supports_quantification(15.1)
    assert not beam_energy_supports_quantification(20.0)
    assert not beam_energy_supports_quantification(None)

    off_voltage = SpectrumFitResult(
        filename="20kv.msa",
        energy_keV=np.array([1.0]),
        counts=np.array([1.0]),
        fit=np.array([1.0]),
        background=np.array([0.0]),
        composition_at={"Fe": 1.0},
        composition_wt={"Fe": 1.0},
        beam_energy_kV=20.0,
    )
    text = format_composition_txt(off_voltage)
    assert "NOT valid" in text
    assert "20.000 kV" in text


def test_emsa_is_accepted_extension():
    assert ".emsa" in SUPPORTED_UPLOAD_EXTENSIONS
    assert ".msa" in SUPPORTED_UPLOAD_EXTENSIONS


def test_load_msa_accepts_emsa_alias():
    src = _example_msa()
    with tempfile.TemporaryDirectory(prefix="autoemx_emsa_") as tmp:
        dest = Path(tmp) / "spectrum_0.emsa"
        shutil.copy2(src, dest)

        energy_msa, counts_msa, meta_msa = load_msa(str(src))
        energy_emsa, counts_emsa, meta_emsa = load_msa(str(dest))

        np.testing.assert_allclose(energy_msa, energy_emsa)
        np.testing.assert_allclose(counts_msa, counts_emsa)
        assert meta_emsa.get("NPOINTS") == meta_msa.get("NPOINTS")


def test_parse_emsa_geometry_from_phenom_header():
    _, _, metadata = load_msa(str(_example_msa()))
    geometry = parse_emsa_geometry(metadata)
    assert geometry["beam_energy"] > 0
    assert geometry["det_ch_width"] > 0
    assert geometry["emergence_angle"] > 0
    assert geometry["sp_collection_time"] > 0


def test_parse_elements_validates_and_deduplicates():
    assert parse_elements("Bi, Fe, O") == ["Bi", "Fe", "O"]
    assert parse_elements("bi fe o, Fe") == ["Bi", "Fe", "O"]
    with pytest.raises(ValueError, match="Unrecognized"):
        parse_elements("Bi, Xx")


def test_composition_txt_and_png_export():
    result = SpectrumFitResult(
        filename="demo.emsa",
        energy_keV=np.linspace(0.1, 10.0, 50),
        counts=np.linspace(10.0, 100.0, 50),
        fit=np.linspace(12.0, 95.0, 50),
        background=np.full(50, 8.0),
        composition_at={"Bi": 0.2, "Fe": 0.3, "O": 0.5},
        composition_wt={"Bi": 0.55, "Fe": 0.25, "O": 0.20},
        analytical_error=0.012,
        r_squared=0.999,
        reduced_chi_sq=1.4,
        quant_flag=0,
        els_sample=["Bi", "Fe", "O"],
        els_substrate=["C", "O", "Al"],
        is_particle=True,
    )
    text = format_composition_txt(result)
    assert "Bi" in text and "At%" in text and "55.00" in text
    assert "15 kV" in text

    fig = fitted_spectrum_figure(result)
    png = figure_to_png_bytes(fig)
    fig.clf()
    assert png[:8] == b"\x89PNG\r\n\x1a\n"


def test_reader_failure_builds_report():
    with tempfile.TemporaryDirectory(prefix="autoemx_bad_msa_") as tmp:
        bad = Path(tmp) / "broken.msa"
        bad.write_text(
            "#FORMAT : not a real EMSA file\n#XUNITS : furlongs\n1, 2, 3\n",
            encoding="utf-8",
        )
        with pytest.raises(SpectrumReadError) as caught:
            load_uploaded_spectrum(bad)
        report = caught.value.report
        assert "broken.msa" in report
        assert "XUNITS" in report
        assert REPORT_EMAIL in report
        assert "Traceback" in report

    mail = mailto_reader_report_url("broken.msa", "bad units")
    assert mail.startswith(f"mailto:{REPORT_EMAIL}?")
    assert "broken.msa" in mail
