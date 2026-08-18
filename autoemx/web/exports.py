#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PNG and TXT exports for fitted spectra and compositions."""

from __future__ import annotations

import io
from typing import Iterable, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from autoemx.web.pipeline import (
    QUANT_BEAM_KV,
    SpectrumFitResult,
    beam_energy_supports_quantification,
)


def format_composition_txt(result: SpectrumFitResult) -> str:
    """Human-readable composition table for download."""
    lines = [
        "# AutoEMX quantification",
        f"# File: {result.filename}",
        f"# Sample elements: {', '.join(result.els_sample) or '(none)'}",
        f"# Substrate elements: {', '.join(result.els_substrate) or '(none)'}",
        f"# Particle geometry: {result.is_particle}",
        f"# Required beam energy for a valid composition: {QUANT_BEAM_KV:.0f} kV",
    ]
    if result.beam_energy_kV is not None:
        lines.append(f"# Spectrum beam energy (header): {result.beam_energy_kV:.3f} kV")
        if not beam_energy_supports_quantification(result.beam_energy_kV):
            lines.append(
                f"# WARNING: Composition is NOT valid. Spectrum was not collected at {QUANT_BEAM_KV:.0f} kV."
            )
        else:
            lines.append(
                f"# Composition validity: spectrum beam energy matches the {QUANT_BEAM_KV:.0f} kV standards."
            )
    if result.r_squared is not None:
        lines.append(f"# R-squared: {result.r_squared:.6f}")
    if result.reduced_chi_sq is not None:
        lines.append(f"# Reduced chi-squared: {result.reduced_chi_sq:.2f}")
    if result.analytical_error is not None:
        lines.append(f"# Analytical error (w%): {result.analytical_error * 100:.2f}")
    if result.quant_flag is not None:
        lines.append(f"# Quantification flag: {result.quant_flag}")
    if result.error:
        lines.append(f"# Error: {result.error}")
    lines.append("#")
    lines.append("Element\tAt%\tWt%")

    elements = list(result.composition_at.keys()) or list(result.composition_wt.keys())
    if not elements:
        lines.append("# (no composition available)")
    for el in elements:
        at_pct = result.composition_at.get(el, 0.0) * 100.0
        wt_pct = result.composition_wt.get(el, 0.0) * 100.0
        lines.append(f"{el}\t{at_pct:.2f}\t{wt_pct:.2f}")
    lines.append("")
    return "\n".join(lines)


def format_batch_composition_txt(results: Iterable[SpectrumFitResult]) -> str:
    """One TXT with a section per spectrum."""
    blocks = [format_composition_txt(result) for result in results]
    return "\n".join(blocks)


def fitted_spectrum_figure(result: SpectrumFitResult, title: Optional[str] = None):
    """Matplotlib overlay of data, fit, and background (does not show)."""
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.plot(result.energy_keV, result.counts, "o", ms=2.5, label="Data", color="C0")
    ax.plot(result.energy_keV, result.fit, "-", lw=1.4, label="Fit", color="C1")
    ax.plot(
        result.energy_keV,
        result.background,
        "--",
        lw=1.6,
        label="Background",
        color="C3",
    )
    y_max = float(np.nanmax(result.counts)) if len(result.counts) else 1.0
    for energy, label in result.peak_labels:
        ax.text(
            energy,
            y_max * 1.02,
            label,
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title(title or f"Fitted spectrum — {result.filename}")
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig


def figure_to_png_bytes(fig) -> bytes:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=200, bbox_inches="tight")
    buffer.seek(0)
    return buffer.read()
