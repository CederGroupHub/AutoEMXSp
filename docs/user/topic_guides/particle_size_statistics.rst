.. _particle_size_statistics:

Particle size statistics
========================

How AutoEMX turns the particles measured in SEM frames into statistics of the sample.
To run a measurement, see the :ref:`particle size tutorial <particle_size_tutorial>`.

.. contents:: On this page
   :local:
   :depth: 1


Size bands
----------

Particle size is the area-equivalent diameter, :math:`d = 2\sqrt{A/\pi}`. A frame
(10 × the largest diameter) only measures particles covering at least
``min_par_pixels`` pixels (default 80, ~10 pixels across). Wider size ranges are
split into **size bands**, each with frames matched to its sizes, and
``n_par_target`` particles are collected in each band.


Corrections
-----------

Each particle gets a weight, so that the statistics describe the sample and not the
way it was imaged:

- **Frame edge.** Particles touching the frame edge are left out, which would
  under-count large particles. Each measured particle (bounding box
  :math:`w \times h` in a :math:`W \times H` frame) is counted
  :math:`c = WH / ((W-w-7)(H-h-7))` times, the Miles-Lantuéjoul correction
  [Miles1974]_ [Lantuejoul1978]_.
- **Bands.** Each band counts in proportion to its particles per unit area scanned,
  :math:`\sum c_i / S_b`, not to the number collected in it.
- **EDS particles** of the same sample, chosen automatically, are included in the
  band of their size (``include_eds_auto_particles``), without changing its weight.

The unweighted statistics are kept in ``ledger.json`` for comparison.


Reported statistics
-------------------

- **Number-weighted:** mean, standard deviation, D10, D25, D50 (median), D75, D90.
- **Volume-weighted**, for comparison with laser diffraction (each particle counted
  as a sphere, :math:`\propto d^3`): ``D43`` (:math:`D[4,3]` [ISO9276]_), ``Dv10``,
  ``Dv50``, ``Dv90``. They are dominated by the largest particles and need many more
  of them to be precise.
- **95 % confidence intervals**, by bootstrap [Efron1993]_: the frames of each band
  are resampled, with their particles and area. Given once each band has 5 frames.

Weighted percentiles use type 5 of Hyndman and Fan [Hyndman1996]_.


Stopping on precision
---------------------

Set ``size_stats_ci_tolerance`` (e.g. ``0.05``) to collect particles until the
confidence interval of ``size_stats_ci_quantities`` (default ``["D50"]``) is within
± that fraction, up to ``size_stats_max_n_frames`` frames (default 500).
``n_par_target`` is then the minimum per band.


Output
------

``ledger.json`` (particles and ``particle_stats``), ``<sample>_Par_sizes.csv``
(with each particle's ``Weight``), ``<sample>_Par_size_stats.csv`` (with
``*_CI95_low/high``), and a plot of the weighted distribution
(``size_stats_plot_format``). An interrupted measurement resumes when run again
with the same settings.


References
----------

.. [Miles1974] R. E. Miles, "On the elimination of edge effects in planar sampling",
   in *Stochastic Geometry* (Harding and Kendall, eds.), Wiley, 1974, pp. 228–247.
.. [Lantuejoul1978] C. Lantuéjoul, in *Geometrical Probability and Biological
   Structures*, Lecture Notes in Biomathematics 23, Springer, 1978, pp. 323–329.
   https://doi.org/10.1007/978-3-642-93089-8_27
.. [ISO9276] ISO 9276-2, *Representation of results of particle size analysis —
   Part 2*.
.. [Efron1993] B. Efron and R. J. Tibshirani, *An Introduction to the Bootstrap*,
   Chapman & Hall, 1993.
.. [Hyndman1996] R. J. Hyndman and Y. Fan, *The American Statistician* **50**,
   361–365 (1996). https://doi.org/10.1080/00031305.1996.10473566
