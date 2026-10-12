=========================
Quickstart Guide
=========================

This guide introduces how to use ``AutoEMX`` for typical workflows:

- Identify how many phases you have in your sample, and measure their compositions via SEM-EDS
- Fit and optionally quantify a single EDS spectrum to evaluate the model performance
- Measure the particle size distribution in your sample via SEM
- Quantify the extent of intermixing of precursors prior a solid-state reaction

.. warning::

   If this is the first time `AutoEMX` is run on your microscope, note that there are a few steps required to set it up before `AutoEMX` can be properly run. Refer to :ref:`Maintainer <advanced_user_index>` docs.

.. warning::

   Ensure the EDS detector is periodically recalibrated for optimal EDS quantification. See :ref:`EDS Detector Calibration <advanced_sdd_calib>`.

Expected runtime on a normal desktop computer:

- Quantification typically takes a 0.5-3 minutes per spectrum.
- For multiple spectra, quantification is parallelized, so runtime does not
   scale linearly with the number of spectra. 100 spectra can be quantified in 10 minutes with 10 CPU units.
- Composition analysis (clustering/statistical analysis) typically takes a few seconds.


Workflows
----------------------------------------------------

Each workflow can be run in two ways, which use the same engine and write the same sample
folders:

- **GUI**: the AutoEMX GUI, a point-and-click interface running locally in your browser.
  Create its double-click app once with ``python -m autoemx.gui --create-launcher``
  (:ref:`Tutorial <gui_getting_started_tutorial>`).
- **Scripts**: the scripts in ``autoemx/scripts/``, which require a minimal set of user-defined
  parameters.


EDS compositional analysis for phase identification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With one click, `AutoEMX` handles the full workflow from EDS spectral acquisition and quantification, to
rule-based filtering of the quantified compositions and unsupervised machine-learning analysis to identify
the different phase compositions in your sample.

- **GUI**: Acquisition, Quantification and Analysis tabs, with interactive 3D plots of the
  results (:ref:`Tutorial <gui_comp_analysis_tutorial>`).
- **Scripts**: ``Run_Acquisition.py`` (:ref:`Tutorial <comp_analysis_tutorial>`).

See full workflow at:

    A. Giunto *et al.*, *Accurate SEM-EDS Quantification, Automation, and Machine
    Learning Enable High-Throughput Compositional Characterization of Powders*,
    *Nature Communications* **17**, 9735 (2026).
    https://doi.org/10.1038/s41467-026-76633-x


Fit and quantify EDS spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fit--and optionally quantify-- EDS spectra acquired using ``AutoEMX`` or exported by commercial
EDS software (.msa, .emsa, .msg spectra files), one at a time to evaluate the model performance,
or in batch.

- **GUI**: Single spectrum tab for one spectrum; Quantification tab for many spectra,
  including folders of exported spectra (:ref:`Tutorial <gui_fit_spectra_tutorial>`).
- **Scripts**: ``Fit_Quant_Single_AutoEMX_Spectrum.py`` (:ref:`Tutorial <fit_autoemx_spectrum_tutorial>`)
  or ``Fit_Quant_Single_MSA_Spectrum.py`` (:ref:`Tutorial <fit_msa_spectrum_tutorial>`) for one
  spectrum, which print the full process in the terminal, the employed fit parameters and their
  final values; ``quantify_external_spectra.py`` for exported spectra in batch
  (:ref:`Tutorial <quantify_external_spectra_tutorial>`).


Measure particle size distribution via SEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Have `AutoEMX` collect multiple images, detect particles, and quantify their size distribution.

- **Scripts**: ``Collect_Particle_Statistics.py`` (:ref:`Tutorial <particle_size_tutorial>`).
  Not available in the GUI yet.


Quantify the extent of intermixing in precursor powders
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use EDS to evaluate the extent of spatial intermixing of different precursor powders, known to affect
the output of solid-state reactions. `AutoEMX` offers a method to quantify the intermixing, helping
the rationalization of impurity formation in solid-state reactions. See for example Fig. 6 in:

    Chem. Mater. 2025, 37, 6807−6822 (https://pubs.acs.org/doi/10.1021/acs.chemmater.5c01573)

- **GUI**: as for the compositional analysis, with the precursor options
  (:ref:`Tutorial <gui_comp_analysis_tutorial>`, last section).
- **Scripts**: ``Run_Acquisition.py`` with ``is_known_phase_mixture = True`` (:ref:`Tutorial <precursor_mix_tutorial>`).
