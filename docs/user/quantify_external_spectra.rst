.. _quantify_external_spectra_tutorial:


Tutorial: Quantify externally exported spectra in batch
=======================================================

This tutorial shows how to import, quantify, and optionally cluster externally
exported EDS spectra using the ``quantify_external_spectra.py`` script.

Use this workflow when spectra were acquired outside AutoEMX (for example with
commercial EDS software) and exported as ``.msa`` or ``.msg`` files.

For single-file inspection and fitting, use
:ref:`fit_msa_spectrum_tutorial` instead.


What this script does
---------------------

``autoemx/scripts/quantify_external_spectra.py`` performs the full pipeline:

- Copies spectra into the AutoEMX sample layout (``sample_ID/spectra``)
- Creates a complete ``ledger.json`` for each sample
- Runs quantification
- Optionally runs clustering/statistical analysis


Step 1 - Open script to edit
----------------------------

Open ``autoemx/scripts/quantify_external_spectra.py``.


Step 2 - Define samples and spectra locations
---------------------------------------------

Edit the ``samples`` list. Each sample is a dictionary with:

- ``ID`` (required): sample identifier; output folder name.
- ``els`` (required): list of sample elements to quantify.
- ``spectra_dir`` (optional): folder containing exported ``.msa``/``.msg`` files.
  If omitted, files are expected already in ``<results_dir>/<ID>/spectra/``.
- ``cnd`` (optional): candidate phases for clustering (for example ``['PbMoO4']``).

Example:

.. code-block:: python

   samples = [
       {
           'ID': 'Wulfenite_imported',
           'els': ['Pb', 'Mo', 'O'],
           'spectra_dir': '/path/to/wulfenite_spectra',
           'cnd': ['PbMoO4'],
       }
   ]


Step 3 - Set instrument and substrate parameters
------------------------------------------------

Set parameters so they match the acquisition conditions of the exported spectra:

- ``microscope_ID`` and ``measurement_mode``
- ``beam_energy`` (must match export acquisition voltage)
- ``els_substrate``, ``sample_substrate_type``, ``sample_substrate_shape``
- ``sample_type`` (``powder``, ``bulk``, ``bulk_rough``, ``powder_continuous``)

These values affect calibration lookup and quantification corrections.


Step 4 - Configure quantification and clustering
------------------------------------------------

Common parameters to adjust:

- ``interrupt_fits_bad_spectra``: speed-focused early interruption of poor fits
- ``min_bckgrnd_cnts``: minimum background counts under reference peaks
- ``max_analytical_error``: filtering threshold used for clustering
- ``quant_flags_accepted``: accepted quantification flags in clustering
- ``max_n_clusters``: upper bound for cluster search
- ``use_project_specific_std_dict``: load standards from project folder

If you only want quantification now, set ``run_analysis=False`` in the function
call block.


Step 5 - Set output location
----------------------------

Set ``results_dir``:

- ``None``: uses ``<current working directory>/Results``
- absolute path: writes to your chosen project folder

Each sample output is saved in ``<results_dir>/<sample_ID>/``.


Step 6 - Run the script
-----------------------

Run from your AutoEMX environment:

.. code-block:: bash

   python autoemx/scripts/quantify_external_spectra.py


Output
------

For each successfully processed sample, AutoEMX writes:

- ``ledger.json`` with spectra registry and quantification configuration
- copied/standardized spectra files under ``spectra/``
- quantification outputs in the sample folder
- clustering outputs when ``run_analysis=True``


Notes and troubleshooting
-------------------------

- Supported input extensions are ``.msa``, ``.msg``, and ``.emsa``.
- If ``spectra_dir`` contains no valid files, the sample is skipped.
- If a ``ledger.json`` already exists, ingestion is skipped by default.
  Set ``overwrite_existing=True`` to rebuild from current script settings.
- For re-analysis with different clustering filters (without re-ingestion), use
  ``run_analysis.py``.
