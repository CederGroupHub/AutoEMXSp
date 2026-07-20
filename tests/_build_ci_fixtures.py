#!/usr/bin/env python3
"""Build minimal CI ledger fixtures under tests/inputs/."""

from __future__ import annotations

import json
import shutil
from copy import deepcopy
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
INPUTS = Path(__file__).resolve().parent / "inputs"
WULF_SRC = REPO / "examples" / "Results" / "Wulfenite_example"
K412_SRC = REPO / "examples" / "Results" / "K-412_NISTstd_example"


def build_wulfenite_mini() -> Path:
    dst = INPUTS / "Wulfenite_mini"
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "spectra").mkdir(exist_ok=True)

    wulf = json.loads((WULF_SRC / "ledger.json").read_text(encoding="utf-8"))
    keep_ids = {"0", "1"}
    spectra = []
    for src_spectrum in wulf["spectra"]:
        if str(src_spectrum["spectrum_id"]) not in keep_ids:
            continue
        spectrum = deepcopy(src_spectrum)
        spectrum["quantification_results"] = []
        rel = spectrum["spectrum_relpath"]
        src_msa = WULF_SRC / rel
        dst_msa = dst / rel
        dst_msa.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_msa, dst_msa)
        spectra.append(spectrum)

    if len(spectra) != 2:
        raise RuntimeError(f"Expected 2 Wulfenite spectra, got {len(spectra)}")

    quants = deepcopy(wulf["quantifications"])
    for quant in quants:
        for analysis in quant.get("clustering_analyses") or []:
            analysis["result"] = None

    # sample_path must be absolute for schema validation; load_sample_ledger
    # rewrites it to the directory that contains ledger.json at runtime.
    mini = {
        "sample_id": "Wulfenite_mini",
        "sample_path": str(dst.resolve()),
        "configs": wulf["configs"],
        "spectra": spectra,
        "active_quant": 0,
        "quantifications": quants,
    }
    (dst / "ledger.json").write_text(json.dumps(mini, indent=2) + "\n", encoding="utf-8")
    return dst


def build_k412_cluster_mini() -> Path:
    dst = INPUTS / "K412_cluster_mini"
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "spectra").mkdir(exist_ok=True)

    k412 = json.loads((K412_SRC / "ledger.json").read_text(encoding="utf-8"))
    good = []
    for spectrum in k412["spectra"]:
        qrs = spectrum.get("quantification_results") or []
        if not qrs:
            continue
        qr = qrs[0]
        if qr.get("quant_flag") not in (0, -1):
            continue
        an_err = qr.get("analytical_error")
        if an_err is None or abs(an_err) > 0.05:
            continue
        si = (qr.get("composition_weight_fractions") or {}).get("Si")
        if si is None:
            continue
        good.append((si, spectrum))

    good.sort(key=lambda item: item[0])
    picked = [good[i][1] for i in (0, 1, 2, -3, -2, -1)]
    selected = []
    seen = set()
    for spectrum in picked:
        sid = str(spectrum["spectrum_id"])
        if sid in seen:
            continue
        seen.add(sid)
        entry = deepcopy(spectrum)
        rel = entry.get("spectrum_relpath") or f"spectra/spectrum_{sid}.msa"
        entry["spectrum_relpath"] = rel
        entry["instrument_background_relpath"] = None
        src_msa = K412_SRC / rel
        if not src_msa.is_file():
            raise FileNotFoundError(f"Missing source MSA for K-412 spectrum {sid}: {src_msa}")
        dst_msa = dst / rel
        dst_msa.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_msa, dst_msa)
        selected.append(entry)

    if len(selected) < 5:
        raise RuntimeError(f"Need >=5 K-412 compositions, got {len(selected)}")

    quants = deepcopy(k412["quantifications"])
    for quant in quants:
        for analysis in quant.get("clustering_analyses") or []:
            analysis["result"] = None
            cfg = analysis["config"]
            cfg["k_forced"] = 2
            cfg["k_finding_method"] = "forced"
            cfg["k_resolved"] = 2

    mini = {
        "sample_id": "K412_cluster_mini",
        "sample_path": str(dst.resolve()),
        "configs": k412["configs"],
        "spectra": selected,
        "active_quant": 0,
        "quantifications": quants,
    }
    (dst / "ledger.json").write_text(json.dumps(mini, indent=2) + "\n", encoding="utf-8")
    return dst


def main() -> None:
    wulf = build_wulfenite_mini()
    k412 = build_k412_cluster_mini()
    print(f"Wrote {wulf}")
    print(f"Wrote {k412}")


if __name__ == "__main__":
    main()
