#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Streamlit UI: upload EMSA spectra, fit/quantify, download PNG and TXT."""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# Repo root on sys.path so the app runs on Streamlit Cloud without `pip install -e .`
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import streamlit as st

from autoemx.web.exports import (
    figure_to_png_bytes,
    fitted_spectrum_figure,
    format_batch_composition_txt,
    format_composition_txt,
)
from autoemx.web.pipeline import (
    QUANT_BEAM_KV,
    SUPPORTED_UPLOAD_EXTENSIONS,
    SpectrumFitResult,
    beam_energy_supports_quantification,
    fit_uploaded_spectrum,
    parse_elements,
)
from autoemx.web.reader_report import (
    GITHUB_NEW_ISSUE_URL,
    REPORT_EMAIL,
    SpectrumReadError,
    mailto_reader_report_url,
)

_DEMO_FILE_CAP = 2
_ACCEPT = sorted(SUPPORTED_UPLOAD_EXTENSIONS)
_LICENSE_CONTACT = "IPO@lbl.gov"


def _is_hosted_demo() -> bool:
    if os.environ.get("AUTOEMX_DEMO"):
        return True
    # Streamlit Community Cloud mounts the repo here.
    return os.path.isdir("/mount/src")


def _write_upload(upload, dest_dir: Path) -> Path:
    suffix = Path(upload.name).suffix.lower()
    if suffix not in SUPPORTED_UPLOAD_EXTENSIONS:
        raise ValueError(
            f"Unsupported file '{upload.name}'. Use {', '.join(_ACCEPT)}."
        )
    dest = dest_dir / Path(upload.name).name
    dest.write_bytes(upload.getbuffer())
    return dest


def _failed_result(
    filename: str,
    els_sample,
    els_substrate,
    is_particle: bool,
    error: str,
    error_stage: str,
    reader_report: str | None = None,
) -> SpectrumFitResult:
    return SpectrumFitResult(
        filename=filename,
        energy_keV=np.array([]),
        counts=np.array([]),
        fit=np.array([]),
        background=np.array([]),
        composition_at={},
        composition_wt={},
        els_sample=list(els_sample),
        els_substrate=list(els_substrate),
        is_particle=is_particle,
        error=error,
        error_stage=error_stage,
        reader_report=reader_report,
    )


def _results_table(result: SpectrumFitResult):
    rows = []
    elements = list(result.composition_at.keys()) or list(result.composition_wt.keys())
    for el in elements:
        rows.append(
            {
                "Element": el,
                "At%": round(result.composition_at.get(el, 0.0) * 100.0, 2),
                "Wt%": round(result.composition_wt.get(el, 0.0) * 100.0, 2),
            }
        )
    return rows


def main() -> None:
    st.set_page_config(page_title="AutoEMX spectrum fit", layout="wide")
    hosted = _is_hosted_demo()

    st.title("AutoEMX spectrum fitting")
    st.caption(
        "Upload one or more EMSA spectra (``.msa``, ``.emsa``, ``.msg``), "
        "fit and quantify them, then save the overlay as PNG and the composition as TXT."
    )
    st.info(
        "**Free use for non-commercial use only.** "
        f"Contact [{_LICENSE_CONTACT}](mailto:{_LICENSE_CONTACT}) "
        "for commercial purposes."
    )
    st.warning(
        f"**Compositions are valid only for spectra collected at {QUANT_BEAM_KV:.0f} kV.** "
        "The shipped peak-to-background standards are for 15 kV. "
        "A spectrum acquired at any other accelerating voltage can still be fitted, "
        "but the reported composition must not be treated as quantitative."
    )

    if hosted:
        st.caption(
            "This is a **public demo**. Each fit typically takes 0.5–3 minutes, "
            "and the app may sleep after idle time. For production work, install "
            "AutoEMX on your machine (`pip install autoemx`) and run "
            "`python -m autoemx.web`."
        )
    else:
        st.markdown(
            "Local AutoEMX GUI. You must specify which elements are present — "
            "the engine does not auto-identify unknown peaks."
        )

    with st.sidebar:
        st.header("Sample")
        els_sample_text = st.text_input(
            "Sample elements",
            value="Bi, Fe, O",
            help="Comma-separated symbols. Required.",
        )
        els_substrate_text = st.text_input(
            "Substrate elements",
            value="C, O, Al",
            help="Typical carbon-tape substrate is C, O, Al.",
        )
        is_particle = st.checkbox(
            "Particle / powder geometry",
            value=True,
            help="Uncheck for bulk / flat samples.",
        )
        st.header("Files")
        st.caption(
            f"Collect spectra at **{QUANT_BEAM_KV:.0f} kV**. "
            "Compositions from other voltages are not valid."
        )
        uploads = st.file_uploader(
            "EMSA spectra",
            type=[ext.lstrip(".") for ext in _ACCEPT],
            accept_multiple_files=True,
        )
        run = st.button("Fit and quantify", type="primary", use_container_width=True)

    if "results" not in st.session_state:
        st.session_state.results = None

    if not run and st.session_state.results is None:
        st.markdown(
            "Provide elements, upload spectra, then click **Fit and quantify**. "
            f"Spectra **must be collected at {QUANT_BEAM_KV:.0f} kV** for the "
            "composition to be valid. An example file is "
            "`autoemx/scripts/input/Example_spectrum.msa` (Bi–Fe–O on carbon tape, 15 kV)."
        )
        return

    if not run:
        results = st.session_state.results
        _render_results(results)
        return

    try:
        els_sample = parse_elements(els_sample_text)
        els_substrate = parse_elements(els_substrate_text)
    except ValueError as exc:
        st.error(str(exc))
        return

    if not els_sample:
        st.error("Enter at least one sample element.")
        return
    if not uploads:
        st.error("Upload at least one spectrum file.")
        return

    if hosted and len(uploads) > _DEMO_FILE_CAP:
        st.error(
            f"This demo accepts at most {_DEMO_FILE_CAP} files. "
            "Install AutoEMX locally to quantify larger batches."
        )
        return

    progress = st.progress(0, text="Starting…")
    results: list[SpectrumFitResult] = []
    n_files = len(uploads)

    with tempfile.TemporaryDirectory(prefix="autoemx_web_") as tmp:
        tmp_dir = Path(tmp)
        for i, upload in enumerate(uploads):
            progress.progress(
                i / n_files,
                text=f"Fitting {upload.name} ({i + 1}/{n_files})…",
            )
            try:
                path = _write_upload(upload, tmp_dir)
                result = fit_uploaded_spectrum(
                    path,
                    els_sample=els_sample,
                    els_substrate=els_substrate,
                    is_particle=is_particle,
                )
            except SpectrumReadError as exc:
                result = _failed_result(
                    upload.name,
                    els_sample,
                    els_substrate,
                    is_particle,
                    error=str(exc),
                    error_stage="read",
                    reader_report=exc.report,
                )
            except Exception as exc:
                result = _failed_result(
                    upload.name,
                    els_sample,
                    els_substrate,
                    is_particle,
                    error=str(exc),
                    error_stage="fit",
                )
            results.append(result)
        progress.progress(1.0, text="Done.")

    st.session_state.results = results
    _render_results(results)


def _render_results(results: list[SpectrumFitResult]) -> None:
    ok = [r for r in results if not r.error]
    failed = [r for r in results if r.error]
    st.success(f"Finished {len(ok)} of {len(results)} spectrum(s).")
    if failed:
        st.warning(
            "Failed: " + ", ".join(f"{r.filename} ({r.error})" for r in failed)
        )

    if ok:
        st.download_button(
            "Download all compositions (TXT)",
            data=format_batch_composition_txt(ok),
            file_name="autoemx_compositions.txt",
            mime="text/plain",
        )

    for result in results:
        st.divider()
        st.subheader(result.filename)
        if result.error:
            st.error(result.error)
            if result.error_stage == "read":
                _render_reader_report(result)
            continue

        metrics = st.columns(3)
        metrics[0].metric(
            "R²",
            f"{result.r_squared:.5f}" if result.r_squared is not None else "—",
        )
        metrics[1].metric(
            "Reduced χ²",
            f"{result.reduced_chi_sq:.1f}" if result.reduced_chi_sq is not None else "—",
        )
        metrics[2].metric(
            "Analytical error",
            (
                f"{result.analytical_error * 100:.2f} w%"
                if result.analytical_error is not None
                else "—"
            ),
        )
        if result.beam_energy_kV is not None:
            if beam_energy_supports_quantification(result.beam_energy_kV):
                st.caption(
                    f"Header beam energy: {result.beam_energy_kV:.3g} kV "
                    f"(matches the required {QUANT_BEAM_KV:.0f} kV)."
                )
            else:
                st.error(
                    f"This spectrum was collected at **{result.beam_energy_kV:.3g} kV**, "
                    f"not {QUANT_BEAM_KV:.0f} kV. The composition below is **not valid** "
                    "with the shipped 15 kV standards."
                )

        fig = fitted_spectrum_figure(result)
        st.pyplot(fig, clear_figure=False)
        png = figure_to_png_bytes(fig)
        fig.clf()

        table = _results_table(result)
        if table:
            st.dataframe(table, hide_index=True, use_container_width=True)
        else:
            st.info("No composition values were returned for this spectrum.")

        col_png, col_txt = st.columns(2)
        stem = Path(result.filename).stem
        with col_png:
            st.download_button(
                "Save fitted spectrum (PNG)",
                data=png,
                file_name=f"{stem}_fitted.png",
                mime="image/png",
                key=f"png-{result.filename}",
            )
        with col_txt:
            st.download_button(
                "Save composition (TXT)",
                data=format_composition_txt(result),
                file_name=f"{stem}_composition.txt",
                mime="text/plain",
                key=f"txt-{result.filename}",
            )


def _render_reader_report(result: SpectrumFitResult) -> None:
    st.markdown(
        "The EMSA/MSA reader could not parse this file. Download a report and "
        f"email it to **{REPORT_EMAIL}** so the reader can be updated for this variant. "
        "Attach the original spectrum as well if you can share it."
    )
    report = result.reader_report or result.error or ""
    stem = Path(result.filename).stem
    col_dl, col_mail = st.columns(2)
    with col_dl:
        st.download_button(
            "Download reader report",
            data=report,
            file_name=f"{stem}_emsa_reader_report.txt",
            mime="text/plain",
            key=f"reader-report-{result.filename}",
        )
    with col_mail:
        st.link_button(
            "Open email to maintainer",
            mailto_reader_report_url(result.filename, result.error or ""),
        )
    st.caption(
        f"You can also open a GitHub issue at {GITHUB_NEW_ISSUE_URL} "
        "and attach the same report."
    )


main()
