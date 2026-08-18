.. _web_gui_tutorial:

Tutorial: Fit and quantify spectra with the AutoEMX GUI
=======================================================

The AutoEMX GUI wraps the same fitting and quantification engine used by
``Fit_Quant_Single_MSA_Spectrum.py``. Use it to upload one or more EMSA
spectra, inspect the fitted overlay and compositions, and save them as PNG
and TXT.

This is the **local, production** interface. A public hosted page, if
available, is a demo only (slow, may sleep, limited number of files).

.. warning::

   Spectra **must be collected at 15 kV**. The shipped peak-to-background
   standards are for 15 kV only. A fit may still run at another voltage, but
   the reported composition is **not valid**.

Launch
------

::

   pip install autoemx
   python -m autoemx.web

Then open the URL Streamlit prints (typically http://localhost:8501).

What to enter
-------------

- **Sample elements** (required): comma-separated symbols, e.g. ``Bi, Fe, O``.
  The fitter does not identify unknown elements.
- **Substrate elements**: e.g. ``C, O, Al`` for carbon tape.
- **Particle / powder geometry**: leave checked for particles; uncheck for
  bulk / flat samples.
- **Files**: ``.msa``, ``.emsa``, or ``.msg``, collected at **15 kV**.

An example spectrum is ``autoemx/scripts/input/Example_spectrum.msa``
(Bi–Fe–O on carbon tape).

Runtime
-------

Each spectrum typically takes 0.5–3 minutes. Default calibrations are
Phenom XL, **15 kV**, point mode. Compositions are valid only at that voltage.

If the file cannot be read
--------------------------

If the EMSA/MSA reader fails, the GUI shows **Download reader report** and
**Open email to maintainer**. The report includes the error, traceback, EMSA
header, and the first data lines — enough to extend ``load_msa``. Attach the
original spectrum to the email if you can share it.

Hugging Face demo (optional)
----------------------------

The GUI can be published later as a Hugging Face Space pointing at this
repository (``autoemx/web/app.py``). Set ``AUTOEMX_DEMO=1`` or rely on the
Space-provided ``SPACE_ID`` environment variable to show the demo banner
and cap uploads at two files. A Hugging Face account is not required to
run the GUI locally. Never commit Hugging Face tokens to the repository.
