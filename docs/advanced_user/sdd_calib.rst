.. _advanced_sdd_calib:

EDS Detector Calibration
================================================

This guide shows how to calibrate the energy scale of the microscope Silicon Drift Detector (SDD) with the
`AutoEMX` script ``Run_SDD_Calibration.py`` (in ``autoemx/scripts/Microscope calibration/``).

Each detector channel *i* corresponds to the energy ``offset + i * scale`` (keV). These two detector channel
parameters drift over time, and shifted peak positions degrade spectral fitting and quantification.
Recalibrate them regularly, and after any service of the detector.

How it works
------------

The script measures two pure bulk standards, fits their spectra, and compares the fitted center of one reference
peak per standard with its theoretical energy. The new offset and scale follow from this two-point calibration.

1. Spectra are acquired automatically on each standard (by default 5 spectra of 50,000 counts).
2. The spectra are fitted, and the centers of the reference peaks are averaged over the spectra accepted in the fit
   (quantification flag 0).
3. The new detector channel parameters are saved as a new file
   ``<date>_detector_channel_params_calibs.json``. `AutoEMX` uses the most recent file for all following acquisitions.

Setting up the script
---------------------

Edit the variables at the top of the script:

- ``microscope_ID``: the folder of your microscope in ``autoemx/calibrations``.
- ``std_list``: exactly two standards, each defined by

  - ``ID``: a name of your choice;
  - ``formula``: the composition, for accurate spectral fitting;
  - ``ref_el`` and ``ref_peak``: the element and line used as reference (e.g. ``'Cu'``, ``'Ka1'``);
  - ``pos``: the position of the standard on the stage (mm);
  - ``sample_type``: ``'bulk'`` (preferred) or ``'powder'``;
  - ``is_manual_meas``: ``True`` to select the spots to measure manually.

- ``n_target_spectra``, ``target_Xsp_counts``: number of spectra per standard, and their counts.

.. tip::

   Choose pure elements with intense, isolated peaks that are well separated in energy, preferably above 2 keV,
   for example Cu K\ :sub:`α`\ (8.05 keV) and Al K\ :sub:`α`\ (1.49 keV), as in the default configuration.

Running the script
------------------

.. code-block:: bash

   python "autoemx/scripts/Microscope calibration/Run_SDD_Calibration.py"

The script prints the current and new offset and scale, and saves the new file.

To repeat only the fitting and calculation on spectra already acquired, set ``reuse_session`` to the name of the
session folder (e.g. ``'20260521_09h22m'``): no new spectra are acquired.

Where the calibrations are saved
--------------------------------

The detector channel parameters belong to your instrument, and are kept in its local calibration folder,
``~/.autoemx/calibrations/<microscope_ID>/Detector_channel_params_calibs/``
(on Windows, ``C:\Users\<user>\.autoemx\...``):

- ``<date>_detector_channel_params_calibs.json``: the detector channel parameters of each calibration, the most
  recent being used. A file also holds the beam current (spot size) used for each measurement mode.
- ``SDD calibrations/<session>/``: the spectra, ledgers and images of each calibration session.

The first time `AutoEMX` runs on a computer, this folder is created with the initial values shipped with the
package (``autoemx/calibrations/<microscope_ID>/``). After that, `AutoEMX` never changes the package files.
