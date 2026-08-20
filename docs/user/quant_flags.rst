.. _quant_flags:

Quantification flags (``quant_flag``)
=====================================

Every fitted or quantified spectrum is assigned an integer ``quant_flag``.
The flag records whether AutoEMX considers the result reliable and, if not,
why. Flags are stored on each spectrum's quantification record (and appear as
``Quant_flag`` in legacy ``Data.csv`` outputs).

At analysis time, ``quant_flags_accepted`` (in the clustering configuration) decides which of
those flags are kept for compositional clustering. Changing
``quant_flags_accepted`` does **not** re-run quantification; it only changes
which already-quantified compositions enter the clustering step.

Default for analysis workflows is ``quant_flags_accepted = [0, -1]``.


Flag meanings
-------------

======= ========================================================================
Flag    Meaning
======= ========================================================================
``0``   Quantification is OK (may still carry large analytical error).
``-1``  Same as ``0``, but quantification did not converge within the iteration
        limit.
``1``   Error during EDS acquisition / no spectral data. No fit executed.
``2``   Total counts below 90% of the target acquisition counts (often bad
        segmentation or a truncated acquisition). Fitting is skipped when
        ``interrupt_fits_bad_spectra=True``.
``3``   Too little low-energy signal (background under ~2 keV too low), which
        compromises P/B ratios. Fitting is skipped when
        ``interrupt_fits_bad_spectra=True``.
``4``   Poor fit (reduced chi-squared too high relative to total counts). Fit
        may be interrupted when ``interrupt_fits_bad_spectra=True``.
``5``   Excessively high analytical error (> 50 wt%). Often caused by a
        missing element or another major model error. Fit may be interrupted
        when ``interrupt_fits_bad_spectra=True``.
``6``   Excessive X-ray absorption (particle geometry). Fit may be
        interrupted when ``interrupt_fits_bad_spectra=True``.
``7``   Excessive substrate contamination (a substrate-element peak intensity
        above 10% of total counts).
``8``   Background counts under a reference peak are too low (i.e. below ``min_bckgrnd_cnts``).
``9``   Fit interrupted for an unknown reason.
``10``  Experimental-standards measurements only: a required reference peak is
        missing.
======= ========================================================================


How to use them
---------------

- Prefer keeping ``[0, -1]`` for phase identification unless you intentionally want
  to inspect poorer spectra.
- To include borderline spectra in clustering, add their flags to
  ``quant_flags_accepted`` (for example, you may consider including ``8``, which is often quite accurate anyways (see Fig. 4 in the paper)).
- Prefit issues (``1``, ``2``, ``3``) can stop fitting early when
  ``interrupt_fits_bad_spectra=True``. Flags ``2`` and ``3`` can still be
  written on a completed quantification if fitting was allowed to proceed.

Related parameters: ``interrupt_fits_bad_spectra``, ``min_bckgrnd_cnts``, and
``max_analytical_error_percent``. See :ref:`comp_analysis_tutorial` and
:ref:`quantify_external_spectra_tutorial`.
