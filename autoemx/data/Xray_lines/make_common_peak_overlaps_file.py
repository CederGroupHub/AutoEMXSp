"""
Generate common_peak_overlaps.json: pairs of X-ray lines from different elements that lie
too close to be reliably separated by an EDS detector, so that users can be warned about
potentially inaccurate quantification.

The table is organised by measurement method. Only 'EDS' is defined here. The table stores all
candidate pairs closer than MAX_CANDIDATE_SEPARATION_EV; their severity is assigned at runtime
(autoemx.core.quantifier.peak_overlaps) from the separation in units of the detector sigma at the
line energies, computed with the calibration of the microscope used:

- severe: separation < 2 sigma_eff (Sparrow limit: two equal Gaussians merge into one peak)
- significant: separation < 3 sigma_eff (overlap criterion of the fitter, see Peaks_Model)

with sigma_eff = sqrt((sigma1^2 + sigma2^2) / 2). The resolution model and these thresholds are
specific to energy-dispersive detectors: whoever adds another measurement method to autoemx must
define its own resolution model, severity categories and candidate pairs under a new key of "methods".

Only lines whose area is fitted freely are considered, i.e. lines that are their own weight
reference in XSp_Fitter (Ka1, La1 or Ll when La1 is absent, Ma1 or Mz1). All other lines have
their area tied to a reference line of the same element, so they do not cause ambiguity.
Lines made free by the user through `free_area_el_lines` are not included, as they depend
on the fit configuration.

Line energies are taken from LineEnergies.csv. Re-run this script after modifying
LineEnergies.csv, LineWeights.csv or the reference-line logic of XSp_Fitter.
"""
import json
import os
from itertools import combinations

from pymatgen.core import Element

from autoemx.core.fitter.fitter import get_reference_xray_line
from autoemx.core.fitter.peaks import OVERLAP_SEPARATION_SIGMAS
from autoemx.core.quantifier.quantifier import XSp_Quantifier
from autoemx.data.Xray_lines import LINE_ENERGIES_DF, get_el_xray_lines

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "common_peak_overlaps.json")

# Candidate pairs are stored up to this separation. It must exceed the largest 'significant' threshold,
# i.e. 3 sigma at the highest fitted line energies (~320 eV at 25 keV for a 130 eV FWHM detector at Mn Ka)
MAX_CANDIDATE_SEPARATION_EV = 400
# Severity thresholds, in units of the RMS detector sigma of the two lines
SEVERITY_THRESHOLDS_SIGMA = {"severe": 2, "significant": OVERLAP_SEPARATION_SIGMAS}

SOURCES = {
    "Newbury2009": (
        "D. E. Newbury, 'Mistakes encountered during automatic peak identification of minor and "
        "trace constituents in electron-excited energy dispersive X-ray microanalysis', "
        "Scanning 31 (2009) 91-101, Table II. doi:10.1002/sca.20151"
    ),
}

# Newbury (2009), Table II: 'problem regions' for automatic peak identification (0.2-5 keV),
# with lines converted to the nomenclature of the line database
NEWBURY2009_REGIONS = [
    [("N", "Ka1"), ("Sc", "La1")],
    [("O", "Ka1"), ("V", "La1")],
    [("F", "Ka1"), ("Mn", "La1"), ("Fe", "La1")],
    [("Ne", "Ka1"), ("Ni", "La1")],
    [("Cu", "La1"), ("Pr", "Ma1")],
    [("Na", "Ka1"), ("Zn", "La1"), ("Pm", "Ma1")],
    [("Mg", "Ka1"), ("As", "La1"), ("Tb", "Ma1")],
    [("Al", "Ka1"), ("Br", "La1"), ("Yb", "Ma1")],
    [("Si", "Ka1"), ("Rb", "La1"), ("Sr", "La1"), ("Ta", "Ma1"), ("W", "Ma1")],
    [("P", "Ka1"), ("Zr", "La1"), ("Pt", "Ma1")],
    [("Nb", "La1"), ("Au", "Ma1"), ("Hg", "Ma1")],
    [("S", "Ka1"), ("Mo", "La1"), ("Pb", "Ma1")],
    [("Tc", "La1"), ("Bi", "Ma1")],
    [("Cl", "Ka1"), ("Rh", "La1")],
    [("Ar", "Ka1"), ("Ag", "La1"), ("Th", "Ma1")],
    [("Cd", "La1"), ("U", "Ma1")],
    [("K", "Ka1"), ("In", "La1"), ("U", "Mb")],
    [("Ca", "Ka1"), ("Sb", "La1"), ("Te", "La1")],
    [("Sc", "Ka1"), ("Xe", "La1")],
    [("Ti", "Ka1"), ("Ba", "La1")],
    [("Ti", "Kb1"), ("V", "Ka1"), ("Ce", "La1"), ("Pr", "La1")],
]


def _get_free_lines(el):
    """Return the lines of element `el` that are their own weight reference in XSp_Fitter."""
    el_xray_lines = get_el_xray_lines(el)
    free_lines = []
    for line, info in el_xray_lines.items():
        ref_el_line = get_reference_xray_line(el, line, el_xray_lines, XSp_Quantifier.xray_quant_ref_lines)
        if ref_el_line == f"{el}_{line}":
            free_lines.append({
                "element": el,
                "line": line,
                "energy_keV": round(float(info["energy (keV)"]), 4),
                "is_quant_line": line in XSp_Quantifier.xray_quant_ref_lines,
            })
    return free_lines


def _get_sources(line_a, line_b):
    a = (line_a["element"], line_a["line"])
    b = (line_b["element"], line_b["line"])
    if any(a in region and b in region for region in NEWBURY2009_REGIONS):
        return ["Newbury2009"]
    return []


def build_table():
    free_lines = []
    for Z in LINE_ENERGIES_DF.index:
        free_lines.extend(_get_free_lines(Element.from_Z(int(Z)).symbol))

    pairs = []
    for a, b in combinations(free_lines, 2):
        if a["element"] == b["element"]:
            continue
        sep_ev = abs(a["energy_keV"] - b["energy_keV"]) * 1000
        if sep_ev > MAX_CANDIDATE_SEPARATION_EV:
            continue
        lines = sorted([a, b], key=lambda line: line["energy_keV"])
        pairs.append({
            "lines": lines,
            "separation_eV": round(sep_ev, 1),
            "sources": _get_sources(a, b),
        })
    pairs.sort(key=lambda p: (p["lines"][0]["energy_keV"], p["lines"][1]["energy_keV"]))

    return {
        "description": (
            "Candidate pairs of X-ray lines from different elements that may overlap in the measured spectra, by "
            "measurement method. Only lines whose area is fitted freely (lines that are their own weight reference) "
            "are included. When both elements of a pair are fitted, the deconvolution of the two peaks may be "
            "unreliable, leading to inaccurate quantification of the element whose quantification line is "
            "involved. The resolution model and severity categories are method-specific: any new measurement "
            "method must define its own."
        ),
        "methods": {
            "EDS": {
                "energy_source": "LineEnergies.csv (generated by make_common_peak_overlaps_file.py)",
                "max_candidate_separation_eV": MAX_CANDIDATE_SEPARATION_EV,
                "resolution_model": (
                    "Detector Gaussian sigma(E) = conv_eff * sqrt(elec_noise^2 + E * F / conv_eff), with the "
                    "calibration parameters of the microscope (DetectorResponseFunction.det_sigma); "
                    "sigma_eff = sqrt((sigma1^2 + sigma2^2) / 2)."
                ),
                "severity_thresholds_sigma": SEVERITY_THRESHOLDS_SIGMA,
                "severity_definitions": {
                    "severe": (
                        "Separation < 2 sigma_eff: the two peaks merge into a single peak with no dip (Sparrow "
                        "limit), so their areas are separated only through the assumed line energies and weights."
                    ),
                    "significant": (
                        "2 sigma_eff <= separation < 3 sigma_eff: the fitter treats the peaks as overlapping and "
                        "ties their energy shifts; fitting relies heavily on accurate peak shapes and energies."
                    ),
                },
                "sources": SOURCES,
                "pairs": pairs,
            },
        },
    }

if __name__ == "__main__":
    table = build_table()
    pairs = table["methods"]["EDS"]["pairs"]
    # One pair per line, to keep the file compact and readable
    table["methods"]["EDS"]["pairs"] = "__PAIRS__"
    pairs_text = "[\n" + ",\n".join(" " * 8 + json.dumps(p, ensure_ascii=False) for p in pairs) + "\n      ]"
    text = json.dumps(table, indent=2, ensure_ascii=False).replace('"__PAIRS__"', pairs_text)
    with open(OUTPUT_FILE, "w", encoding="utf-8") as file:
        file.write(text + "\n")
    n_lit = sum(bool(p["sources"]) for p in pairs)
    print(f"Wrote {len(pairs)} candidate EDS pairs ({n_lit} from literature) to {OUTPUT_FILE}")
