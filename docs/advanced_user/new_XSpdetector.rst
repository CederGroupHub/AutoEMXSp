.. _advanced_new_XSp_calibs:

==============================================================
Calibrate a new Detector for EDS fit and quantification
==============================================================

This guide explains how to calibrate the detector model used by `AutoEMX` to fit and quantify EDS spectra, for a new
detector or microscope.

`AutoEMX` models each spectrum with peak shapes and a background whose parameters depend on the detector. These
parameters are fixed during fitting, as functions of the X-ray energy or of the atomic number, calibrated once per
detector on spectra of standards. The calibration of the Phenom XL is shipped with `AutoEMX` and described in the
Supplementary Information of Giunto *et al.*, Nature Communications 17, 9735 (2026),
https://doi.org/10.1038/s41467-026-76633-x (Supplementary Notes *Parameter calibration for peak shape modeling* and
*Calibration of P, F and β background model parameters*).

.. note::

   Calibration is required for optimal accuracy, not to use `AutoEMX`: the parameters of the Phenom XL generalize
   reasonably well to other detectors.

The energy calibration of the detector (offset and gain) is a separate, regular calibration: see
:ref:`advanced_sdd_calib`.

The detector model
------------------

The detector model of each microscope is saved in the package, in ``autoemx/calibrations/<microscope_ID>/``:

- ``detector_model.json``: the values of the model, validated when loaded (see
  ``autoemx/config/schema_models/detector_model.py`` for the meaning of each value). It contains:

  - detector characteristics: number of channels, emergence (take-off) angle, Si density, undetectable elements,
    available measurement modes;
  - **peak width**: Fano factor ``F``, electronic noise ``elec_noise``, conversion efficiency ``conv_eff``;
  - **peak shape** (``peak_shape``): parameters of the low-energy tailing and incomplete charge collection (ICC),
    ``gamma``, ``f_tail``, ``R_e`` and ``F_loss``, each a sigmoid of the line energy;
  - **background** (``background``): functions of the atomic number ``P(Z)``, ``F(Z)`` and ``β(Z)`` of the modified
    Duncumb model of the generated background;
  - fine-tuning factors: upper bounds of escape and pileup peaks, Ll/Kα ratio, background and strobe peak scaling.

- ``detector_efficiency.msa``: the detector efficiency as a function of energy, as provided by the manufacturer.

To add a new microscope, create a folder named with its microscope ID, starting from a copy of the files of a
similar microscope, and set the detector characteristics. The microscope also needs a driver (see
:ref:`advanced_new_EM_driver`): by default, the driver named as the microscope ID. For a second instrument of a
supported model, with its own calibrated detector, set ``"driver"`` in its ``detector_model.json`` to the driver of
the model (e.g. ``"driver": "PhenomXL"``).

.. note::

   The detector model is read from the package. To calibrate it, work in a git checkout of `AutoEMX` installed in
   editable mode (``pip install -e .``), edit ``detector_model.json`` there, and consider contributing the
   calibration of your microscope to `AutoEMX`.

Calibration procedure
---------------------

The parameters of the model are correlated: the fitted value of each depends on the values of the others.
**Calibrate one parameter at a time**, updating ``detector_model.json`` after each step, in the following order,
and repeat the sequence until the values no longer change significantly.

1. **Energy calibration.** Calibrate the detector channel parameters (:ref:`advanced_sdd_calib`).

2. **Spectra of standards.** Acquire high-count spectra of bulk, polished standards, ideally pure elements, spanning
   the energy range of the X-ray lines (peak shape) and the atomic numbers of your samples (background), e.g. with
   ``Run_Experimental_Standard_Collection.py`` and ``update_std_library = False`` (:ref:`advanced_new_exp_stds`).
   The paper used ten to fifteen spectra of 250,000 counts per standard. Acquire them at the acceleration voltage of
   your measurements: the background functions, β in particular, are expected to depend on it.

3. **Peak width**: ``Calibrate_Peak_Sigma_Params.py`` (below).

4. **Peak shape**: ``Calibrate_Peak_Shape_Params.py``, one parameter at a time (below).

5. **Background**: ``Calibrate_Background_Params.py``, in the order P, β, F (below).

6. **Reference P/B values.** Once the model is calibrated, refit the standards with
   ``Refit_Experimental_Standards.py`` and ``update_std_library = True`` (or collect new ones), since the reference
   values depend on the spectral fit.

All calibration scripts are in ``autoemx/scripts/Microscope calibration/``. They print the fitted values with the
current ones, plot them, and print the new values to paste into ``detector_model.json``: they never change the model
themselves. After changing the model, refit the standards before calibrating the next parameter. The detector
response matrices computed from the model are recomputed automatically (about 10 s, then cached in
``~/.autoemx/cache``).

.. tip::

   Fitted values at the bounds allowed during fitting are flagged (⚠️ at fit bound). They mean that the current
   model is far from the measured value: update the model, refit, and calibrate again.

Peak width
^^^^^^^^^^

The Gaussian width of the peaks depends on the line energy *E* (keV) as

.. math::

   \sigma = \mathrm{conv\_eff} \sqrt{\mathrm{elec\_noise}^2 + E \cdot F / \mathrm{conv\_eff}}

with the Fano factor *F* fixed (0.12). ``Calibrate_Peak_Sigma_Params.py`` reads the peak widths fitted in the
spectra of the standards from their ledgers (Kα\ :sub:`1`, Lα\ :sub:`1` and Mα\ :sub:`1` lines, or the lines you
list), and fits ``conv_eff`` and ``elec_noise``.

.. code-block:: python

   samples = {
       '/path/to/Std_measurements/Al': None,                 # Ka1, La1, Ma1 lines
       '/path/to/Std_measurements/Cu': ['Cu_Ka1', 'Cu_La1'],
   }

Widths can only vary between 0.93 and 1.1 times the current model during fitting: if many are at these bounds,
update the model, refit the standards, and calibrate again.

Peak shape
^^^^^^^^^^

The peaks of lines below 1.84 keV have a low-energy tail (``gamma``: skewness, ``f_tail``: tail fraction), and those
between 1.18 and 5 keV are broadened by incomplete charge collection (``R_e``: charge cloud radius, ``F_loss``:
fraction of X-rays lost at the detector surface); between 1.18 and 1.84 keV, both apply. Each parameter is a sigmoid
of the line energy *E*, ``c + L / (1 + exp(-k (E - x0)))``, with coefficients ``offset`` (*c*), ``amplitude`` (*L*),
``center`` (*x0*) and ``steepness`` (*k*).

``Calibrate_Peak_Shape_Params.py`` refits the spectra of the standards with the chosen parameters fitted freely for
the chosen elements, averages the fitted values of each line, and fits the chosen sigmoid coefficients:

.. code-block:: python

   samples = {
       '/path/to/Std_measurements/C': ['C'],     # elements whose peak shapes are fitted
       '/path/to/Std_measurements/Al': ['Al'],
   }
   params = ['gamma']                            # parameters fitted freely
   sigmoid_coeffs = ['offset', 'amplitude', 'center', 'steepness']

A sigmoid is fitted only with more lines than coefficients. With few standards, fit fewer coefficients (e.g. only
``offset``, to shift the current function).

Calibrate the parameters one at a time, for example ``gamma`` and ``f_tail`` first, on lines below 1.18 keV where
only the tail applies, then ``R_e`` and ``F_loss`` on lines between 1.84 and 5 keV where only the ICC applies. Then
check the lines between 1.18 and 1.84 keV, where both apply.

Background
^^^^^^^^^^

The generated background (modified Duncumb model) depends on three functions of the atomic number *Z*: the shape
factor ``P(Z)`` (a quadratic polynomial), which sets the curvature of the continuum, the attenuation factor ``β(Z)``
(a sigmoid plus a Gaussian), and the intensity factor ``F(Z) = 1 - exp(-(Z + z_offset) / z_scale)``. For compounds,
P and F are evaluated at the mean atomic number of the sample, and β is the mass-weighted average of the elemental
values: calibrate with **pure elements**, spanning the atomic numbers of your samples (C to Bi in the paper).

``Calibrate_Background_Params.py`` refits the spectra with the chosen parameters fitted freely (``params``), the
others being set by the detector model, and fits the functions of the chosen parameters (``fit_functions``). Follow
the order of the paper:

1. ``params = ['P', 'F', 'beta']``, ``fit_functions = ['P']``: all three free, fit ``P(Z)``. Update the model.
2. ``params = ['beta']``, ``fit_functions = ['beta']``: with P fixed, fit ``β(Z)``. Update the model.
3. ``params = ['F']``, ``fit_functions = ['F']``: with P and β fixed, fit ``F(Z)``. Update the model.

Optionally, fit only some coefficients of a function with ``coeffs``, e.g. ``{'P': ['a0']}``. A function is fitted
only with more standards than coefficients (P: 3, F: 2, β: 7).

Fine-tuning factors
^^^^^^^^^^^^^^^^^^^

The remaining factors (``escape_peak_probability``, ``pileup_peak_probability``, ``weight_Ll_ref_Ka1``,
``gen_background_time_scaling_factor``, ``strobe_peak_int_factor``, ``zero_strobe_peak_sigma``) are upper bounds or
initial guesses for parameters fitted in each spectrum, and are not critical: reasonable estimates suffice. To check
them, compare with the values fitted in a few high-count spectra of standards.
