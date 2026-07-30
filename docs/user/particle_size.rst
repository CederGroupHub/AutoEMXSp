.. _particle_size_tutorial:


Tutorial: Measure particle size distribution via SEM
=======================================================================

This tutorial explains how to collect particle-size statistics with
``AutoEMX`` using the script:

``autoemx/scripts/Particle size Statistics measurements/Collect_Particle_Statistics.py``

It also explains how to post-process results and, for advanced users, how to
add a separate particle-segmentation model.


What this script does
---------------------

``Collect_Particle_Statistics.py``:

- Acquires SEM frames across the selected sample region.
- Segments particles in each frame.
- Filters particles by size.
- Continues until the requested number of particles is reached (or frames are exhausted).
- Saves per-particle measurements and summary statistics.
- Saves a particle-size histogram.


Step 1 - Open script to edit
----------------------------

Open:

``autoemx/scripts/Particle size Statistics measurements/Collect_Particle_Statistics.py``

When you run it, a GUI appears to collect the main inputs.


Step 2 - Define required inputs in the GUI
------------------------------------------

Fill in the following fields:

- ``Sample ID``
	Name used to create the output folder and output files.

- ``Sample Center X, Y``
	Stage coordinates (mm) of the center of the sample.

- ``Number of Particles to Analyse``
	Target number of particles to analyse.

- ``Min & Max Particle Diameter (um)``
	Accepted particle-diameter range (used to compute area filters).
	Ranges spanning more than one order of magnitude are split automatically into
	decade bands (largest first). Each band uses a search FOV matched to its maximum
	particle size. Coarse frames use normal labels (for example ``A0``); finer bands
	scan sub-frames inside visited parent FOVs with hierarchical ids such as
	``A0_a0``. Particle geometry is stored in ``ledger.json`` as ``ParticleInfo``
	records (area, equivalent diameter, absolute stage coordinates in mm, and
	``frame_id``), and the usual Par_sizes / Par_size_stats CSV outputs are still
	written from that ledger.

- ``Carbon Tape Diameter (mm)``
	Used to estimate the effective sample area to scan over the carbon tape.

- ``Working Distance (mm)``
	Nominal working distance used by autofocus to prevent catastrophic focus drift.

- ``Segmentation Model``
	Particle segmentation model. Default is ``threshold_bright``.
	This default model works best when there is clear contrast between carbon
	tape and particles. For this reason, electron backscatter detector mode is
	recommended to maximize Z-contrast.
	Particles should also be well separated; touching/overlapping particles can
	be segmented as a single particle.
	For more complex samples, machine-learning segmentation models are recommended
	(for example ``Rettenberger2024``).
	You can also add your own segmentation model (see the advanced section below).

- ``Manual Navigation``
	If enabled, you manually navigate to regions/particles and determine the frame size.

- ``Auto Detect Carbon Tape``
	Enables automatic substrate detection for carbon tape workflows.

- ``Auto Adjust Brightness and Contrast``
	If set to ``No``, you must also provide:

	- ``Brightness``
	- ``Contrast``

After submitting, choose ``results_dir`` (save directory) in the folder picker.


Step 3 - Run the script
-----------------------

Note that the measurement runs using the microscope's current settings, including detector
type, beam current, and acceleration voltage.

Run from your AutoEMX environment, for example:

.. code-block:: bash

	 python "autoemx/scripts/Particle size Statistics measurements/Collect_Particle_Statistics.py"


Output
------

In the sample output directory, the workflow writes:

- ``ledger.json``
	Sample ledger whose ``particles`` list holds ``ParticleInfo`` records used for
	the size-distribution analysis (including absolute stage coordinates in mm and
	hierarchical ``frame_id`` values when multiple size bands are scanned).

- ``<SampleID>_Par_sizes.csv``
	Per-particle table with particle ID, frame ID, area, and equivalent diameter.

- ``<SampleID>_Par_size_stats.csv``
	Summary statistics (mean, stdev, median, min, max, D10, D25, D75, D90).

- ``<SampleID>_Par_size_distribution_hist.png``
	Histogram plot of equivalent particle diameters.

- SEM images collected during the scan (depending on script options).


Optional: Reprocess particle statistics after filtering
-------------------------------------------------------

If some particles or full frames were incorrectly segmented and you want to
exclude them before recomputing the size distribution, use:

``autoemx/scripts/Particle size Statistics measurements/Process_Particle_Stats_Files.py``

Main editable fields in that script are:

- ``sample_ID``: sample folder to process
- ``input_dir``: parent path where sample folder is stored
- ``particles_IDs_to_filter``: list of particle IDs to ignore
- ``frame_IDs_to_filter``: list of frame IDs to ignore (if a full frame is bad, all its particles are ignored)

The script generates processed outputs with suffix ``_processed``.


Optional (Advanced): Add a separate segmentation model
------------------------------------------------------

Most users do not need this section. Use it only if the built-in segmentation
models are not suitable for your images.

``AutoEMX`` automatically discovers custom segmentation models from:

``autoemx/core/em_runtime/particle_segmentation_models/``

To add a new model:

1. Create a new Python file in that folder, for example:
	 ``my_segmentation_model.py``

2. Use
	 ``autoemx/core/em_runtime/particle_segmentation_models/segmentation_model_template.py``
	 as the template.

3. Implement a function with this signature:

.. code-block:: python

	 def segment_particles(frame_image, powder_meas_config=None, save_image=False, EM=None):
			 ...
			 return par_mask

4. Return either:

	 - a binary mask (background 0, particles > 0), or
	 - a labeled image (background 0, each particle with a unique positive label).

5. Restart the particle-statistics script. Your model name should appear in the
	 ``Segmentation Model`` dropdown automatically.

6. Select it in the GUI (or set ``par_segmentation_model`` directly in
	 ``powder_meas_cfg_kwargs``).

Notes:

- If your model requires external files (for example ONNX weights), keep paths
	stable and preferably relative to the model file location.
- If a model name is invalid, AutoEMX falls back to ``threshold_bright``.