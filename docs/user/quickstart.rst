=========================
Quickstart Guide
=========================

This guide introduces how to use ``AutoEMX`` for typical workflows:

- Identify how many phases you have in your sample, and measure their compositions via SEM-EDS
- Fit and optionally quantify a single EDS spectrum to evaluate the model performance
- Measure the particle size distribution in your sample via SEM
- Quantify the extent of intermixing of precursors prior a solid-state reaction

.. warning::

   If this is the first time `AutoEMX` is run on your microscope, note that there are a few steps required to set it up before `AutoEMX` can be properly run. Refer to :ref:`Advanced User <advanced_user_index>` docs.

.. warning::

   Ensure the EDS detector is periodically recalibrated for optimal EDS quantification. See :ref:`EDS Detector Calibration <advanced_sdd_calib>`.

Expected runtime on a normal desktop computer:

- Quantification typically takes a 0.5-3 minutes per spectrum.
- For multiple spectra, quantification is parallelized, so runtime does not
   scale linearly with the number of spectra. 100 spectra can be quantified in 10 minutes with 10 CPU units.
- Composition analysis (clustering/statistical analysis) typically takes a few seconds.


Workflows
----------------------------------------------------

``AutoEMX`` comes with a selection of scripts (located at autoemx/scripts/) that require a minimal set of user-defined
parameters for running:


EDS compositional analysis for phase identification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the ``run_acquisition_quant_analysis.py`` script (:ref:`Tutorial <comp_analysis_tutorial>`):

With one click, `AutoEMX` handles the full workflow from EDS spectral acquisition and quantification, to
rule-based filtering of the quantified compositions and unsupervised machine-learning analysis to identify
the different phase compositions in your sample.

See full workflow at:

    A. Giunto *et al.*, *Accurate SEM-EDS Quantification, Automation, and Machine
    Learning Enable High-Throughput Compositional Characterization of Powders*,
    Nature Communications (2026), in press.

    (https://doi.org/10.1038/s41467-026-76633-x)


Fit and quantify a single EDS spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the ``fit_quant_single_autoemx_spectrum.py`` (:ref:`Tutorial <fit_autoemx_spectrum_tutorial>`)
or the ``fit_quant_single_msa_spectrum.py`` script (:ref:`Tutorial <fit_msa_spectrum_tutorial>`):

Fit--and optionally quantify-- a single EDS spectrum acquired using ``AutoEMX`` or exported by commercial
EDS software (.msa, .emsa, .msg spectra files).

These scripts print the full process in the terminal, the employed fit parameters and their final values. 
They also show the fitted spectrum for visual evaluation of goodness of fit.

Alternatively, launch the local GUI (same engine; PNG and TXT download)::

   python -m autoemx.web

See :ref:`web_gui_tutorial`.


Measure particle size distribution via SEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the ``collect_particle_statistics.py`` script (:ref:`Tutorial <particle_size_tutorial>`).

Have `AutoEMX` collect multiple images, detect particles, and quantify their size distribution.



Quantify the extent of intermixing in precursor powders
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the ``Run_Acquisition_PrecursorMix.py`` script (:ref:`Tutorial <precursor_mix_tutorial>`):

Use EDS to evaluate the extent of spatial intermixing of different precursor powders, known to affect
the output of solid-state reactions. `AutoEMX` offers a method to quantify the intermixing, helping
the rationalization of impurity formation in solid-state reactions. See for example Fig. 6 in:

    Chem. Mater. 2025, 37, 6807−6822 (https://pubs.acs.org/doi/10.1021/acs.chemmater.5c01573)
