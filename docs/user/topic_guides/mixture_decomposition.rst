.. _mixture_decomposition:

Mixture decomposition
=====================

A cluster does not always come from a single phase. When the electron beam probes
several phases at once (fine intermixed powders, unreacted precursors, partly
reacted samples), AutoEMX describes the cluster as a **mixture of candidate
phases** (``ref_formulae``): which phases, and in which proportions.

.. admonition:: Defaults: change them only if unsatisfied with the results

   The decomposition runs by default and only needs the candidate phases,
   ``ref_formulae``. A phase that is not a candidate cannot be found.

.. contents:: On this page
   :local:
   :depth: 1


The idea
--------

Each spot measures a weighted average of the phases it probes, so the spots of a
mixture fill the region between its phases:

- a **binary** mixture spreads along the **line** between two phases;
- a **ternary** mixture spreads over the **triangle** between three phases.

.. list-table::
   :class: borderless

   * - .. figure:: /_static/Example_mixture_plot_binary.png
          :alt: Binary mixture in a ternary diagram
          :width: 100%

          **Binary.** CaCO₃ + 2 Ta₂O₅ after reaction: the spots lie along the
          Ta₂O₅–Ca₄Ta₂O₉ line.
     - .. figure:: /_static/Example_mixture_plot_ternary.png
          :alt: Ternary mixture in a ternary diagram
          :width: 100%

          **Ternary.** 2 SrCO₃ + Ta₂O₅, partly reacted: the spots fill the
          triangle between the precursors and the product Sr₄Ta₂O₉.

AutoEMX uses **the fewest phases that reproduce the cluster**: pairs first, then
three or four phases only if no pair fits. Clusters that are compact and match a
single candidate phase are not decomposed.


Reading the results
-------------------

Each mixture in ``Clusters.csv`` gives its phases (``Mix``), their molar fractions
(``Mol_Ratio``, ``X_means``) and a confidence score, ``CS_mix``, between 0 and 1:

- **above ~0.7**: the phases reproduce the cluster within measurement noise;
- **below ~0.5**: they do not. A phase is probably missing from ``ref_formulae``,
  and AutoEMX also reports two phases of unknown composition, with approximate
  formulas.

Several mixtures can be listed per cluster, best first. Only the most plausible
are shown; the ``Mix_more`` row says how many more are saved in ``ledger.json``
and why they were left out. All columns are described in :doc:`../tutorials/scripts/comp_analysis`.


Plots
-----

``Mixture_plot_cl<i>.png`` (and ``_zoomed``) shows the best mixture of cluster
``i`` with its phases (blue stars), the other candidate phases (grey stars) and
the spots. Check that the spots actually lie along the line, or inside the
triangle, of the reported phases.

- **3 elements**: ternary diagram of the elements.
- **2 phases, more elements**: the two elements that differ most between the
  phases (or ``els_to_plot``).
- **3 phases, 4 or more elements**: ternary diagram of the phase fractions.

.. figure:: /_static/Example_mixture_plot_2d.png
   :alt: Binary mixture of a six-element sample in 2D
   :width: 45%
   :align: center

   Mixture of nepheline (NaAlSiO₄) and LiMn₁.₅Ni₀.₅O₄: Mn and Si are the elements
   that differ most between the two phases. Each spot is linked to its position on
   the mixture line (red crosses), which gives its phase fractions.

For known phase mixtures (``is_known_phase_mixture = True``), a violin plot
of the molar fractions is also drawn. Set ``intermixing_violin_plots`` to ``True``
or ``False`` to draw them or not regardless.


Things to keep in mind
----------------------

- **Phases on the same line cannot be told apart.** If a candidate lies on the
  line between two others (e.g. CuAl₂O₄ between CuO and Al₂O₃), different pairs
  fit equally well. AutoEMX reports the pair closest to the spots and hides the
  equivalent ones.
- **Elements that are not quantified are ignored** (e.g. H, Li, and C when not in
  the sample): LiNb₂O₅ and Nb₂O₅ look the same in atomic fraction (``at_fr``) space.
  They will look different in mass fraction (``wt_fr``) space, as the Li mass will be
  taken into account
- **With 3 elements, three phases surrounding the cluster always fit it.** A
  ternary result means that no pair fits, but which three phases are named
  depends on the candidates. Candidates close to a pure element (e.g. ``CO``,
  pure O when C is not quantified) can enclose almost anything.


Changing the defaults
---------------------

The options are in ``MixtureParams``, passed through
:func:`analyze_sample <autoemx.runners.analyze_sample.analyze_sample>`. The most
useful one:

.. code-block:: python

   analyze_sample(sample_ID, mixture_params={"max_n_phases": 2})   # binary mixtures only

Use it when a mixture you know is binary is reported with three phases, and the
third phase is not plausible in the mixture plot. The other options are listed
in the ``MixtureParams`` docstring.
