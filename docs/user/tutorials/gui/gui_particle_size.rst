.. _gui_particle_size_tutorial:

Tutorial: particle size statistics in the GUI
=============================================

This tutorial shows how to measure the particle size distribution of powders in the GUI. How
the statistics are computed is described in :ref:`particle_size_statistics`.


Step 1 - Acquire
----------------

In the **Acquisition** tab, choose the **Measurement**:

- **Particle sizes only**: particle size statistics, without EDS.
- **EDS + particle sizes**: EDS compositional analysis, then particle size statistics on the
  same samples. The particles in the frames searched by EDS count towards the target, so only
  the missing frames are acquired.
- **EDS composition**: EDS only.

In **Particle size statistics**, set the diameter range of the particles to measure
(**min/max particle diameter**), the **particles per size band**, and the **particle
segmentation model** (see Step 3). Then start the acquisition as in
:ref:`gui_comp_analysis_tutorial`.


Step 2 - View the statistics
----------------------------

In the **Analysis** tab, select the sample and open **Particle sizes**. Samples without EDS
spectra open directly on it.

- **This sample**: the frames with the particles found (see Step 4), then the size distribution
  and its statistics. What each statistic means: :ref:`psd_reported_statistics`.
- **Compare samples**: the distributions of several samples of the folder, with **Export CSV**.
- **Compare segmentations**: the results of this sample (see Step 3).


.. _gui_segment_again:

Step 3 - Segment again with another model
-----------------------------------------

Every frame image is saved, so the frames can be segmented again: in the **Segment again** section,
choose the **model** and click **Segment again**. Each segmentation gives its own result, kept
next to the others: choose one in the results menu and click **Set as active** to use it by
default.

Models:

- **threshold_bright** (default): fast; set the **threshold** (8-bit brightness of the
  particles). Touching particles are measured as one.
- **SAM2** and **microSAM**: segment each particle on its own, also when touching (~30 s per
  frame for SAM2 on a CPU). They are optional: until installed they are listed as *not
  installed*, and selecting one shows how to install it. Install them in your AutoEMX
  environment, then restart the GUI:

  .. code-block:: bash

     pip install "autoemx[segmentation]"                       # SAM2
     python -m autoemx.utils.install_segmentation micro-sam    # microSAM

  For microSAM, ignore pip's warnings that micro-sam requires PyQt6, segment-anything-py
  and trackastra: they are left out on purpose.

**options** (optional) changes model settings, e.g. ``{"pred_iou_thresh": 0.7}`` for SAM2 to
find more particles.


Step 4 - Leave out frames and particles
---------------------------------------

Below the statistics, the **Frames** browser shows each frame with the outlines of the particles
found by the shown result: counted (green), left out (red), or not counted because they touch the
frame edge or are outside the size range (amber, dashed). Click the dot at the centre of a particle
to leave it out (click again to count it), or tick **Leave this frame out**, then click
**Recalculate and save**. This saves a new result, made the active one; the earlier ones are kept.
