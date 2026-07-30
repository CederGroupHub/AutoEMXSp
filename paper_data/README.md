# Paper Data Reproducibility

All paper results in this directory can be reproduced by:

1. Re-quantifying all spectra.
2. Re-running the clustering analysis.

Samples in the archived data are **already quantified**. The scripts here are provided to reproduce the published results and to re-run analysis (and, optionally, re-quantification) from the paper.

Follow the tutorials in the documentation for the full workflow and parameter guidance:
- [Quantify spectra tutorial](https://cedergrouphub.github.io/AutoEMX/user/quantify_external_spectra.html)
- [EDS compositional analysis tutorial](https://cedergrouphub.github.io/AutoEMX/user/comp_analysis.html)

Practical workflow from this folder:

1. Ensure all paper data files are fully available locally (if needed, fetch Git LFS data).
2. Run the quantification workflow on all spectra using the project scripts.
3. Re-run the clustering/statistical analysis on the newly quantified results.

Repository helper scripts in this folder:
- Run_Quantification_PaperData.py
- Run_Analysis_PaperData.py

Expected runtime notes:
- Quantification takes a few minutes per spectrum.
- For multiple spectra, quantification is parallelized, so runtime does not scale linearly with the number of spectra.
- Analysis typically takes a few seconds.

---

## Folder structure

```
paper_data/
├── README.md                                      # This file
├── Employed_standards.json                        # Standards table used for paper quantification
│
├── Commercial precursors.zip                      # Git LFS — commercial powder precursors
├── Minerals&EDSstandards.zip                      # Git LFS — minerals and EDS standards
├── Synthetic samples.zip                          # Git LFS — synthetic compounds / mixtures
│
├── Run_Quantification_PaperData.py                # Batch re-quantification (+ analysis)
├── Run_Analysis_PaperData.py                      # Re-run clustering / compositional analysis
├── Fit_and_Quantify_Single_Spectrum_PaperData.py  # Inspect one spectrum fit/quantification
└── Benchmark_Composition_Accuracy.py              # Aggregate accuracy plots vs known formulae
```

After fetching LFS objects and unzipping the archives, sample folders appear under `paper_data/` (or inside category folders from each zip). Those extracted directories are **not** committed (see `.gitignore`: `paper_data/*/`); only the zip archives, scripts, and `Employed_standards.json` above are tracked.

Typical sample IDs match the names used in the scripts, for example:

| Archive | Example sample IDs |
| --- | --- |
| `Minerals&EDSstandards.zip` | `Bornite_mineral`, `Benitoite_mineral`, `Wulfenite_mineral`, … |
| `Commercial precursors.zip` | `Al2O3_precursor`, `NiO_precursor`, `LiCoPO4_precursor`, … |
| `Synthetic samples.zip` | `NaGe2(PO4)3`, `Bi2Fe4O9`, `MgCuP2O7`, … |

Approximate archive sizes (after LFS checkout): ~2.3 GB minerals/standards, ~3.0 GB commercial precursors, ~3.7 GB synthetic samples.

---

## Download with Git LFS

A normal `git clone` leaves **pointer files** (~130 bytes) instead of the real `.zip` contents. You will see text like:

```text
version https://git-lfs.github.com/spec/v1
oid sha256:...
size ...
```

To download the datasets, from the repository root:

```bash
# 1. Install Git LFS (only needed once per machine)
git lfs install

# 2. Fetch and materialize LFS files
git lfs fetch --all
git lfs checkout
```

Then extract each archive into `paper_data/`. Confirm a zip is real data by checking its size is gigabytes, not ~135 bytes.

Alternatively, download the zip files from the GitHub UI for this folder:
https://github.com/CederGroupHub/AutoEMX/tree/main/paper_data

`Employed_standards.json` is a normal Git file (not LFS) and is available immediately after clone.

---

## What each script does

Edit `sample_ID` / `sample_IDs` (and related options) at the top of each script before running. AutoEMX must be installed in the active environment.

The quantification and analysis entry points set `results_path` to the repository root, so sample folders under `paper_data/` are discovered automatically. `Benchmark_Composition_Accuracy.py` reads analysis outputs relative to `paper_data/` itself.

### `Run_Quantification_PaperData.py`

Batch re-quantification for a list of paper samples (`batch_quantify_and_analyze`).

- Re-fits and re-quantifies spectra (peak-to-background / PB method by default) to reproduce quantification from scratch.
- Continues across missing or failing samples.
- Can run clustering analysis afterward (runner default).
- Tune `force_requantification`, `min_bckgrnd_cnts`, `num_CPU_cores`, and related flags at the top of the file.

Because the archives already include quantified results, set `force_requantification=True` only when you intentionally want to overwrite existing quantification.

### `Run_Analysis_PaperData.py`

Re-runs clustering and compositional analysis for a **single** sample (`analyze_sample`), using the (already) quantified results.

- Loads quantified results and config for the selected `sample_ID`.
- Optional controls for cluster count (`k_forced` / `k_finding_method`), matrix decomposition, spectral filtering, and plot options.
- Primary script for reproducing the paper’s clustering / statistical analysis without re-quantifying.

### `Fit_and_Quantify_Single_Spectrum_PaperData.py`

Fit and quantify one spectrum from a sample ledger (`fit_and_quantify_spectrum_from_ledger`).

- Choose `sample_ID` and `spectrum_ID`.
- Useful to inspect fitting quality and quantification for a single measurement (plots optional).

### `Benchmark_Composition_Accuracy.py`

Aggregate accuracy evaluation across many samples with known reference formulae.

- Compares quantified / clustered compositions to expected formulae for minerals, commercial precursors, and synthetic samples.
- Builds histograms of absolute and relative deviations, summary statistics, and related plots under `paper_data/Quant Analysis/`.
- Select which sample sets to include (all, binaries only, synthetics only, etc.) and filtering options near the top of the file.

### `Employed_standards.json`

JSON standards library documenting the EDS standards (peak-to-background values, metadata, beam energy, etc.) used for the paper quantification. Kept here for reproducibility and reference alongside the spectral archives.
