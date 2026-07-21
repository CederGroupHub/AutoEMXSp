# Paper Data Reproducibility

All paper results in this directory can be reproduced by:

1. Re-quantifying all spectra.
2. Re-running the clustering analysis.

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
