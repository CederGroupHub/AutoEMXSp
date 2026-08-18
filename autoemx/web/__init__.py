#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Browser UI for fitting and quantifying uploaded EMSA spectra."""

from autoemx.web.pipeline import (
    SUPPORTED_UPLOAD_EXTENSIONS,
    SpectrumFitResult,
    fit_uploaded_spectrum,
    load_uploaded_spectrum,
    parse_elements,
    parse_emsa_geometry,
)
from autoemx.web.reader_report import SpectrumReadError

__all__ = [
    "SUPPORTED_UPLOAD_EXTENSIONS",
    "SpectrumFitResult",
    "SpectrumReadError",
    "fit_uploaded_spectrum",
    "load_uploaded_spectrum",
    "parse_elements",
    "parse_emsa_geometry",
]
