# AutoXSp  

[![PyPI version](https://badge.fury.io/py/autoxsp.svg)](https://pypi.org/project/autoxsp/)
[![Python Version](https://img.shields.io/pypi/pyversions/autoxsp.svg)](https://pypi.org/project/autoxsp/)
[![License](https://img.shields.io/github/license/<your-github-username>/AutoXSp.svg)](LICENSE)
[![Downloads](https://pepy.tech/badge/autoxsp)](https://pepy.tech/project/autoxsp)

**Automated Compositional Characterization via X-Ray Spectroscopy at Electron Microscopes**

AutoXSp is a **fully automated framework** for SEM-EDS workflows — from spectral acquisition and quantification to data filtering and compositional analysis — all in **one click**.

### ✨ Key Features
- **Automated acquisition & quantification** of X-ray spectra. Single spectrum quantification script available
- **Rule-based filtering** to automatically discard poorly quantified spectra from the analysis
- **Machine learning–based compositional analysis** to identify the compositions of individual phases in the sample  
- **Automated experimental standard collection** scripts included
- **Extensible architecture** — adaptable to other techniques such as  
  - Wavelength Dispersive Spectroscopy (WDS)  
  - Scanning Transmission Electron Microscopy (STEM) with EDS  
- **Extensible hardware support** — includes driver for ThermoFisher Phenom Desktop SEM series, and can be extended to any electron microscope with a Python API  

### 📊 Performance
- **Benchmarked** on 74 single-phase samples with compositions spanning **38 elements** (from nitrogen to bismuth), it achieved **<5–10% relative deviation** from expected values  
    *(See publication: TO ADD)* 
- **Machine learning** compositional analysis detects individual phase composition in **multi-phase samples**, including minor phases
- **Intermixed phases** can also be resolved

### 🧪 Supported Use Cases
- Powder, bulk, and rough samples  
- Scanning Electron Microscopy (SEM) with Energy-Dispersive Spectroscopy (EDS)  

### ⚙️ Requirements
- Works on all major platforms  
- Quick installation  
- Requires some calibration for use with different electron microscopes  

---

## 🚀 Installation

You can install **AutoXSp** in just one command.

### Using pip (Python Package Index)
pip install autoxsp

### Or directly from GitHub:
pip install git+https://github.com/CederGroupHub/AutoXSp

---

## 📑 Table of Contents
- [✨ Features](#-features)
- [🎥 Demo](#-demo)
- [🆕 Coming Soon](#-coming-soon)
- [🚀 Installation](#-installation)
- [🖥 Quick Start](#-quick-start)
- [📦 Requirements](#-requirements)
- [📂 Project Structure](#-project-structure)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📬 Contact](#-contact)

---

## ✨ Features
- Automatic processing of X-ray spectral data
- Open-source and customizable for non-commercial-use

---

## 🎥 Demo
- Watch Auto-XSp in action on a desktop SEM-EDS system at https://youtu.be/Bym58gNxlj0

---

## 🆕 Coming Soon
Here’s what’s planned for future releases of **AutoXSp**:
- ⚡ GPU acceleration for faster data processing
- 🐍 Upgrade to **Python 3.12** for improved performance, modern syntax features, and better compatibility with the latest scientific libraries
- 🚀 Integration of a **forked `lmfit`** version accepting `Model.fit(data, fit_kws={'full_output': False})` to avoid covariance computations and speed up computations

---

## 🚀 Installation

You can install **AutoXSp** in just one command.

### Using pip
pip install autoxsp

### Or directly from GitHub:
pip install git+https://github.com/CederGroupHub/AutoXSp

### Using conda
conda install -c conda-forge autoxsp

---

## 🖥 Quick Start

AutoXSp supports two main automated workflows:

1. **Experimental Standard Collection** — acquire and fit EDS/WDS spectra from known-composition samples to generate reference peak-to-background ratios.
2. **Sample Acquisition & Analysis** — acquire spectra from unknown samples, quantify them, and perform compositional phase analysis.

---

### 1️⃣ Acquire Experimental Standards
from autoxsp.runners import batch_acquire_experimental_stds

# Define standards(s) to analyse (additional options available):
# - 'ID': unique standard identifier
# - 'formula': standard composition
# - 'pos': stage position (x, y) in mm
# - 'sample_type': bulk or powder
# - 'is_manual_meas': Manually select spots if standard is not bulk, nor powder

std_list = [
    {
        'id': 'Al_std',
        'formula': 'Al',
        'pos': (0, 0),
        'sample_type': 'bulk',
        'is_manual_meas': False
    },
]

# Run experimental standard acquisition at the microscope computer
batch_acquire_experimental_stds(stds=std_list)


### 2️⃣ Acquire & Analyse Samples
from autoxsp.runners import batch_acquire_and_analyze

# Define sample(s) to analyse (additional options available):
# - 'ID': unique sample identifier
# - 'els': list of possible elements in the sample
# - 'pos': stage position (x, y) in mm
# - 'cnd' (optional): list of candidate phases/formulas

samples = [
    {
        'id': 'Anorthite_mineral',
        'els': ['Ca', 'Al', 'Si', 'O'],
        'pos': (-37.5, -37.5),
        'cnd': ['CaAl2Si2O8']
    },
]

# Run acquisition and analysis at the microscope computer
batch_acquire_and_analyze(samples)


---

## 📦 Requirements

Python 3.8 or newer
Dependencies are installed automatically with pip or conda

Electron Microscope driver developed for the Thermofisher Phenom Desktop SEMs. Will work with any microscope of the series, equipped with PPI (Phenom Programming Interface)

---

## 📂 Project Structure

python
Copy code
AutoXSp/
├── autoxsp/           # Main package source code
├── scripts/           # Helper scripts
├── tests/             # Unit tests
├── LICENSE
├── README.md
└── pyproject.toml

---

## 🤝 Contributing

Contributions are welcome!

Open to collaborations to extend this code to different microscopes or to different types of samples, for example thin films

---

## 📄 License

This project is licensed under an academic, nonprofit, internal, research & development, NON-COMMERCIAL USE ONLY,
LICENSE — see the LICENSE file for details.

---

## 📬 Contact

For questions or issues, please open an issue on GitHub.


---

