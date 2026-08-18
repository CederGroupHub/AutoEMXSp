#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a maintainer report when EMSA/MSA reading fails."""

from __future__ import annotations

import platform
import sys
import traceback
from pathlib import Path
from typing import Optional
from urllib.parse import quote

try:
    from importlib.metadata import PackageNotFoundError, version as pkg_version
except ImportError:  # pragma: no cover
    from importlib_metadata import PackageNotFoundError, version as pkg_version

REPORT_EMAIL = "agiunto@lbl.gov"
GITHUB_NEW_ISSUE_URL = "https://github.com/CederGroupHub/AutoEMX/issues/new"

_MAX_HEADER_LINES = 80
_MAX_DATA_PREVIEW_LINES = 15


class SpectrumReadError(Exception):
    """Raised when an uploaded EMSA/MSA file cannot be parsed."""

    def __init__(self, message: str, report: str):
        super().__init__(message)
        self.report = report


def _package_version() -> str:
    try:
        return pkg_version("autoemx")
    except PackageNotFoundError:
        return "unknown (editable / unpackaged)"


def extract_emsa_file_preview(
    path: Path,
    max_header_lines: int = _MAX_HEADER_LINES,
    max_data_lines: int = _MAX_DATA_PREVIEW_LINES,
) -> str:
    """Header lines plus a few data lines — enough to extend ``load_msa``."""
    raw = path.read_bytes()
    for encoding in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            text = raw.decode(encoding)
            used = encoding
            break
        except UnicodeDecodeError:
            continue
    else:
        return (
            f"(Could not decode file as text; size={len(raw)} bytes. "
            "Please attach the original file to the email.)"
        )

    header: list[str] = []
    data: list[str] = []
    in_data = False
    for line in text.splitlines():
        stripped = line.strip()
        is_spectrum_marker = stripped.startswith("#SPECTRUM") or stripped.startswith(
            "##SPECTRUM"
        )
        if is_spectrum_marker:
            header.append(line.rstrip())
            in_data = True
            continue
        if not in_data:
            if len(header) < max_header_lines:
                header.append(line.rstrip())
            continue
        if stripped and len(data) < max_data_lines:
            data.append(line.rstrip())

    parts = [
        f"Decoded as: {used}",
        f"File size (bytes): {len(raw)}",
        "",
        "--- Header (and SPECTRUM marker) ---",
        *header,
        "",
        f"--- First {len(data)} data line(s) ---",
        *data,
    ]
    if not data and in_data:
        parts.append("(No data lines after the SPECTRUM marker.)")
    if not in_data:
        parts.append("(No #SPECTRUM / ##SPECTRUM marker found.)")
    return "\n".join(parts)


def build_reader_failure_report(
    path: Path,
    exc: BaseException,
    tb_text: Optional[str] = None,
) -> str:
    """Plain-text report for email / GitHub issue attachment."""
    if tb_text is None:
        tb_text = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    try:
        preview = extract_emsa_file_preview(path)
    except OSError as preview_exc:
        preview = f"(Could not read file preview: {preview_exc})"

    lines = [
        "AutoEMX EMSA/MSA reader failure report",
        "Please email this file to " + REPORT_EMAIL,
        "and attach the original spectrum if you can share it.",
        "GitHub: " + GITHUB_NEW_ISSUE_URL,
        "",
        f"AutoEMX version: {_package_version()}",
        f"Python: {sys.version.split()[0]} ({platform.system()} {platform.release()})",
        f"Filename: {path.name}",
        f"Suffix: {path.suffix.lower()}",
        "",
        "--- Error ---",
        f"{type(exc).__name__}: {exc}",
        "",
        "--- Traceback ---",
        tb_text.rstrip(),
        "",
        "--- File preview ---",
        preview,
        "",
    ]
    return "\n".join(lines)


def mailto_reader_report_url(filename: str, error_message: str) -> str:
    """Open the user's mail client; they still attach the downloaded report."""
    subject = f"AutoEMX EMSA reader report: {filename}"
    body = (
        "The AutoEMX EMSA/MSA reader failed on this file.\n\n"
        f"Filename: {filename}\n"
        f"Error: {error_message}\n\n"
        "I have attached the downloaded reader report (.txt). "
        "I can also attach the original spectrum if needed.\n"
    )
    return f"mailto:{REPORT_EMAIL}?subject={quote(subject)}&body={quote(body)}"
