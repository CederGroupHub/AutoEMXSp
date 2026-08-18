.. AutoEMX documentation master file, created by
   sphinx-quickstart on Mon Dec 15 16:11:19 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


AutoEMX Documentation
==========================

Welcome to ``AutoEMX``, a Python package for **Automated Electron Microscopy
X-ray Spectroscopy**:

https://github.com/CederGroupHub/AutoEMX

`AutoEMX` is a framework for running **automated acquisition and analysis
routines** on electron microscopes (EM), offering both X-ray spectroscopy
and image acquisition and analysis workflows.

`AutoEMX` currently supports **Energy-Dispersive X-ray Spectroscopy (EDS)**
in **Scanning Electron Microscopy (SEM)**.

The package was primarily conceived for **automated EDS compositional analysis**,
but it also includes scripts for **automated particle size distribution
measurements** based on SEM imaging.

This work is described in:
    A. Giunto *et al.*, *Accurate SEM‑EDS Quantification, Automation, and
    Machine Learning Enable High‑Throughput Compositional Characterization
    of Powders*, Nature Communications (2026), in press.
    DOI: https://doi.org/10.1038/s41467-026-76633-x
    
Please cite this work if you use ``AutoEMX``.


Key Features
------------

- **Fully automated SEM-EDS compositional analysis of samples**, integrating:
    - **Live SEM control** for particle identification and EDS spectral acquisition
    - **EDS spectra quantification** using the peak-to-background method
    - **Rule-based filtering** of compositions to discard poorly quantified spectra
    - **Unsupervised machine learning–based compositional analysis** to identify the
      compositions of individual phases in the sample
      
- Manual single/multiple EDS spectra quantification

- **Automated experimental standard collection** scripts

- **Automated particle size distribution measurement** scripts

- **Extensible architecture**, adaptable to other techniques such as:
  - Wavelength Dispersive Spectroscopy (WDS)
  - Scanning Transmission Electron Microscopy (STEM) with EDS

- **Extensible hardware support**, including a driver for the
  ThermoFisher Phenom Desktop SEM series, and adaptable to any electron
  microscope exposing a Python API
  

Supported Use Cases
-------------------

- Powders and rough samples (e.g. rough films or pellets)
- Scanning Electron Microscopy (SEM) with Energy-Dispersive Spectroscopy (EDS)


Demo
---------------------------

Watch ``AutoEMX`` in action on a desktop SEM–EDS system:
https://youtu.be/Bym58gNxlj0

For a step-by-step guide to run this workflow, see:
:ref:`Tutorial: EDS compositional analysis for phase identification <comp_analysis_tutorial>`.

Expected runtime on a normal desktop computer:

- Quantification typically takes a few minutes per spectrum.
- For multiple spectra, quantification is parallelized and therefore does not
  scale linearly with the number of spectra.
- Composition analysis (clustering/statistical analysis) typically takes a few
  seconds.


Performance
-----------

- **Benchmarked** on 74 single-phase samples spanning **38 elements**
  (from nitrogen to bismuth), achieving **<5–10% relative deviation** from
  expected values
  
- **Machine learning–based compositional analysis** detects individual phase
  compositions in **multi-phase samples**, including minor phases

- **Intermixed phases** can also be resolved

See https://doi.org/10.1038/s41467-026-76633-x for more details


Requirements
------------

- Python 3.11 or above

- Software dependencies:

  - cvxpy >=1.5.3, <=1.7.3
  - ecos >=2.0.12
  - lmfit >=1.3.2, <=1.3.4
  - matplotlib >=3.5.2, <=3.10.6
  - numpy >=1.26.4, <=2.2.6
  - pandas >=2.2.2, <=2.3.3
  - pydantic >=2.7.0, <=2.11.3
  - pymatgen >=2024.11.13, <=2025.10.7
  - scikit-learn >=1.5.1, <=1.7.2
  - scikit-image >=0.24.0, <=0.25.2
  - scipy >=1.14.1, <=1.16.2
  - seaborn ==0.13.2
  - sympy >=1.13.2, <=1.14.0
  - yellowbrick ==1.5
  - opencv-python >=4.10.0, <=4.12.0.88
  - pillow >=10.0.1

- Electron Microscope provided with an API.
  ``AutoEMX`` comes with a driver for Thermofisher PyPhenom. For different microscopes, 
  the EM_driver must be adapted (see :ref:`EM Driver Set Up <advanced_new_EM_driver>`).
  
- ``AutoEMX`` comes with EDS calibrations for the Thermofisher Phenom XL series. Different microscopes or detectors
  require recalibration (see :ref:`Calibrating EDS <advanced_new_XSp_calibs>`).



Scope of the Documentation
--------------------------

This documentation is intended for **both standard and advanced users** of the
AutoEMX package.

- **Standard users**

  You run predefined scripts without any prior knowledge of the internal code
  structure. The documentation provides **step-by-step instructions** to help
  you get started quickly.


- **Advanced users**

  You interact with AutoEMX beyond simple script execution (workflows setup,
  EDS detector calibration, etc...).
  This documentation guides you through the **initial configuration and setup**
  required to deploy AutoEMX on a **new microscope** and adapt it to new
  experimental workflows.

.. toctree::
   :maxdepth: 2
   :caption: User Documentation:

   user/index
   
   
.. toctree::
   :maxdepth: 2
   :caption: Advanced User Documentation:

   advanced_user/index