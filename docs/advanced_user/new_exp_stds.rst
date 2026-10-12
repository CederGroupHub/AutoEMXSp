.. _advanced_new_exp_stds:

EDS Experimental Standard Collection
================================================

This guide shows how to use `AutoEMX` for automatically collecting peak-to-background (P/B) reference values
from experimental standards of known composition, required for composition extraction from the measured
EDS spectra.

`AutoEMX` quantifies spectra with the P/B method, using a single experimental standard per element: the P/B ratio
of each X-ray line measured in the sample is compared with the reference P/B ratio of the same line measured in a
standard, collected once for each instrument and beam energy (see Giunto *et al.*, Nature Communications 17, 9735
(2026), https://doi.org/10.1038/s41467-026-76633-x).

Choosing the standards
----------------------

- Use **bulk, polished** standards: pure metals when available, otherwise simple binary or ternary compounds
  (in the paper: GaN, BaSO\ :sub:`4`, MgO, CaF\ :sub:`2`, NaCl, InP, KBr, LaB\ :sub:`6`).
- Collect them at the acceleration voltage of your measurements: a standards file is required for every
  acceleration voltage used (``EDS_Stds_<voltage>keV.json``).
- For accurate reference values, average many high-count spectra: in the paper, fifteen spectra of 250,000 counts
  per standard.

Collecting the standards
------------------------

Use ``Run_Experimental_Standard_Collection.py`` (in ``autoemx/scripts/Experimental standards/``). It acquires the
spectra, fits them, and computes the reference P/B values of the standards.

Edit the variables at the top of the script:

- ``std_list``: the standards, each defined by

  - ``ID``: a name of your choice;
  - ``formula``: the composition (e.g. ``'GaN'``);
  - ``pos``: the position of the standard on the stage (mm), or ``()`` with ``'is_manual_meas': True`` to navigate
    to it manually;
  - ``sample_type``: ``'bulk'`` (preferred) or ``'powder'``;
  - optionally ``els_to_use_for_mean_PB_calc``: the elements of this standard used for the reference values
    (default ``['all']`` for bulk standards).

- ``n_target_spectra`` and ``target_Xsp_counts``: number of spectra per standard and their counts.
- ``update_std_library``: if ``True``, the fitted P/B values are written to the standards file of the microscope.
  Set it to ``False`` for trial runs.
- ``exp_stds_meas_cfg_kwargs``: which values are accepted:

  - ``min_acceptable_PB_ratio``: minimum P/B ratio of a peak to be used as reference (weak peaks are discarded);
  - ``quant_flags_accepted``: quantification flags of the spectra used (see :ref:`quant_flags`);
  - ``els_to_use_for_mean_PB_calc`` and ``generate_separate_std_dict``: see the precursor mixtures below.

The settings are saved in the ledger of each standard, so the standards can be refitted later with the same
settings.

.. note::

   Standards of the precursors of a powder mixture (``is_known_phase_mixture = True``) are saved in a separate
   standards file, and not used for the reference values of the elements: see :ref:`precursor_mix_tutorial`.

Where the reference values are saved
------------------------------------

The reference values belong to your instrument, and are kept in its local calibration folder,
``~/.autoemx/calibrations/<microscope_ID>/EDS_Stds_<voltage>keV.json``. For each X-ray line, the file lists the
P/B value of each standard, and their mean used for quantification.

A summary of the standards available (elements, lines, P/B values) is written next to the file,
``experimental_standards_details_<voltage>keV.txt``, every time the standards are updated.

The first time `AutoEMX` runs on a computer, the standards shipped with the package for the microscope are copied
to this folder. After that, `AutoEMX` only reads and updates the local file.

Managing the standards
----------------------

Other scripts in ``autoemx/scripts/Experimental standards/``:

- ``Refit_Experimental_Standards.py``: refit standards already acquired, e.g. after recalibrating the detector model
  (:ref:`advanced_new_XSp_calibs`) or changing which values are accepted, without measuring them again.
- ``Modify_EDS_Standards_File.py``: delete standards or X-ray lines, or exclude standards from the mean reference
  values. A backup of the file is saved before it is changed.
- ``Extract_Experimental_Standards_Details.py``: print the summary of the standards available.
- ``Publish_Standards.py`` (maintainers): copy the local standards file into the package, with its summary, to
  distribute it as the initial standards of the microscope. Runs only in a git checkout of `AutoEMX`.
