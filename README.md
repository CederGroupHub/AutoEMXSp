<div align="center">

# AutoEMXSp 

[![PyPI version](https://badge.fury.io/py/autoemxsp.svg)](https://pypi.org/project/autoemxsp/)
[![Python Version](https://img.shields.io/pypi/pyversions/autoemxsp.svg)](https://pypi.org/project/autoemxsp/)
[![License: Custom Non-Commercial](https://img.shields.io/badge/license-Custom%20Non--Commercial-blue.svg)](https://github.com/CederGroupHub/AutoEMXSp/blob/main/LICENSE.txt)

**Automated Electron Microscopy X-Ray Spectroscopy for Compositional Characterization of Materials**

</div>

AutoEMXSp is a **fully automated framework** for SEM-EDS workflows — from spectral acquisition and quantification to data filtering and compositional analysis — all in **one click**.

📖 This work is described in:  
A. Giunto *et al.*, *Harnessing Automated SEM-EDS and Machine Learning to Unlock High-Throughput Compositional Characterization of Powder Materials*, 2025.  
DOI: [https://doi.org/10.21203/rs.3.rs-7837297/v1](https://doi.org/10.21203/rs.3.rs-7837297/v1)

### ✨ Key Features
- **Automated acquisition & quantification** of X-ray spectra using the peak-to-background method. Single spectrum quantification also available
- **Automated rule-based filtering** of compositions to discard poorly quantified spectra from the analysis
- **Automated machine learning–based compositional analysis** to identify the compositions of individual phases in the sample  
- **Automated experimental standard collection** scripts included
- **Extensible architecture** — adaptable to other techniques such as  
  - Wavelength Dispersive Spectroscopy (WDS)  
  - Scanning Transmission Electron Microscopy (STEM) with EDS  
- **Extensible hardware support** — includes driver for ThermoFisher Phenom Desktop SEM series, and can be extended to any electron microscope with a Python API  

### 📊 Performance
- **Benchmarked** on 74 single-phase samples with compositions spanning **38 elements** (from nitrogen to bismuth), it achieved **<5–10% relative deviation** from expected values  
- **Machine learning** compositional analysis detects individual phase composition in **multi-phase samples**, including minor phases
- **Intermixed phases** can also be resolved

### 🧪 Supported Use Cases
- Powder, bulk, and rough samples  
- Scanning Electron Microscopy (SEM) with Energy-Dispersive Spectroscopy (EDS)  

### ⚙️ Requirements
- Cross-platform: runs on **Linux, macOS, and Windows**
- Quick installation  
- Requires calibration for use with different electron microscopes  

---

## 📑 Table of Contents
- [🎥 Demo](#-demo)
- [🚀 Installation](#-installation)
- [🖥 Quick Start](#-quick-start)
- [📦 Requirements](#-requirements)
- [🆕 Coming Soon](#-coming-soon)
- [📂 Project Structure](#-project-structure)
- [📁 Scripts](#-scripts)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📖 Citation](#-citation)
- [📂 Paper Data](#-paper-data)
- [📬 Contact](#-contact)

---

## 🎥 Demo
- Watch AutoEMXSp in action on a desktop SEM-EDS system at https://youtu.be/Bym58gNxlj0

---

## 🚀 Installation

You can install **AutoEMXSp** in just one command.

### Using pip
```bash
pip install autoemxsp
```

### Or directly from GitHub:
```bash
pip install git+https://github.com/CederGroupHub/AutoEMXSp
```

### Using conda
```bash
conda install -c conda-forge autoemxsp
```

---

## 🖥 Quick Start

AutoEMXSp supports three main automated workflows:

1. **Experimental Standard Collection** — acquire and fit X-ray spectra from known-composition samples to generate reference peak-to-background ratios.
2. **Sample Acquisition & Analysis** — acquire spectra from unknown samples, quantify them, and perform compositional phase analysis.
3. **Particle Size Statistical Analysis** - control EM to search for particles and collect statistics on their size distribution.

---

### 1️⃣ Acquire Experimental Standards
See **Run_Experimental_Standard_Collection.py** script

```python
from autoemxsp.runners import batch_acquire_experimental_stds

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
```

### 2️⃣ Acquire & Analyse Samples
See **Run_Acquisition_Quant_Analysis.py** script

```python
from autoemxsp.runners import batch_acquire_and_analyze

# Define sample(s) to analyse (additional options available):
# - 'id': unique sample identifier
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
```


### 3️⃣ Particle Size Statistical Analysis
See **Collect_Particle_Statistics.py** script

```python
from autoemxsp.runners import collect_particle_statistics

# Define sample(s) to analyse (additional options available):
# - 'id': unique sample identifier
# - 'pos': stage position (x, y) in mm

samples = [
    {
        'id': 'Anorthite_mineral',
        'pos': (-37.5, -37.5),
    },
]

# Run acquisition and analysis at the microscope computer
collect_particle_statistics(samples)
```

---

## 📦 Requirements

- **Python 3.11 or newer**  
- All dependencies are installed automatically via `pip` or `conda`.  
- Tested versions of dependencies are specified in `pyproject.toml`.  
  > The package **may work** with more recent versions, but these have **not been tested**.
  
---

### Electron Microscope Support
- ✅ Developed and tested for **Thermo Fisher Phenom Desktop SEMs**.  
- ✅ Compatible with any Phenom microscope equipped with **PPI (Phenom Programming Interface)**.  
- ⚠️ For other microscope models, the driver must be adapted to the appropriate API commands.  


---

## 🆕 Coming Soon
Here’s what’s planned for future releases of **AutoEMXSp**:
- ⚡ GPU acceleration for faster spectral fitting
- 🐍 Verify with the latest **Python** version for improved compatibility with current scientific libraries
- 📏 New scripts for **spectral parameter calibration** to extend the `XSp_calibs` library to your own instrument.
- 🤖 Integration of ML models for particle segmentation and improved size distribution analysis

---

## 📂 Project Structure

The repository is organized as follows:

```text
AutoEMXSp/
├── autoemxsp/                 # Main package source code
│   ├── core/                   # Core objects and source code
│   ├── runners/                # Runner functions calling on core objects
│   ├── lib/                    # Libraries of X-ray data
│   ├── tools/                  # Miscellaneous helper functions
│   │   └── custom_fnctns.py        # Customizable clustering plot function
│   ├── EM_driver/              # Electron Microscope driver (⚠️ adapt to your own instrument)
│   ├── XSp_calibs/             # X-ray spectral calibrations (⚠️ adapt to your own instrument)
│   ├── scripts/                # Helper scripts (see Scripts below)
│   └── Results/                # Example acquired data (used for unit tests)
│
├── examples/                  # Example scripts for fitting, quantification and compositional analysis of example data
├── tests/                     # Unit tests for fitting, quantification, compositional analysis and image processing
│                               # (Acquisition tests require proper EM drivers & calibration)
├── paper_data/                # Raw paper data uploaded on Git LFS (Dowload instructions in Paper Data section below)
│
├── LICENSE.txt
├── README.md
└── pyproject.toml
```

---

## 📁 Scripts

This repository includes a collection of scripts that streamline the use of **AutoEMXSp**.  
Each script is tailored for a specific task in spectral acquisition, calibration, quantification, or analysis.

### 🔬 Acquisition, Quantification & Analysis
- **Run_Acquisition_Quant_Analysis.py** — Acquire X-ray spectra and optionally perform quantification and composition analysis.  
- **Run_Quantification_Analysis.py** — Quantify acquired spectra (single or multiple samples) and perform machine-learning analysis.  
- **Run_Analysis.py** — Launch customized machine-learning analysis on previously quantified data.  

### 🛠️ Miscellaneous
- **Collect_Particle_Statistics.py** - Analyse sample collecting particle size statistics and distribution.
- **Fit_Quant_Single_Spectrum.py** — Fit and optionally quantify a single spectrum. Prints fitting parameters and plots fitted spectrum for detailed inspection of model performance.  
- **Run_Experimental_Standard_Collection.py** — Acquire and fit experimental standards.  
- **Run_SDD_Calibration.py** — Perform calibration of the SDD detector.

### ⚗️ Characterize Extent of Intermixing in Known Powder Mixtures  
*(see [Chem. Mater. 2015, 27, 20, 7084–7094](https://pubs.acs.org/doi/10.1021/acs.chemmater.5c01573) for example)*  
- **Run_Acquisition_PrecursorMix.py** — Acquire spectra for powder precursor mixtures.  
- **Run_Quantification_PrecursorMix.py** — Quantify spectra for one or multiple powder mixtures and run machine-learning analysis.
- Customized analysis can be performed using the **Run_Analysis.py** script

👉 All scripts can be executed directly from the command line or imported into a Python environment, making them accessible from anywhere on your system.  

---

## 🤝 Contributing

Contributions are welcome!

Open to collaborations to extend this package to different tools or to different types of samples, for example thin films.
Please contact me at agiunto@lbl.gov

---

## 📄 License

This project is licensed under an academic, nonprofit, internal, research & development, NON-COMMERCIAL USE ONLY,
LICENSE — see the LICENSE file for details.

---

## 📖 Citation

If you use **AutoEMXSp** in your research, please cite the following publication:

> A. Giunto, Y. Fei, P. Nevatia, B. Rendy, N. Szymanski and G. Ceder
> *Harnessing Automated SEM-EDS and Machine Learning to Unlock High-Throughput Compositional Characterization of Powder Materials*, 2025.  
> DOI: [https://doi.org/10.21203/rs.3.rs-7837297/v1](https://doi.org/10.21203/rs.3.rs-7837297/v1)

### BibTeX
```bibtex
@article{Giunto2025AutoEMXSp,
  author  = {Giunto, Andrea and Fei, Yuxing and Nevatia, Pragnay and Rendy, Bernardus and Szymanski, Nathan and Ceder, Gerbrand},
  title   = {Harnessing Automated SEM-EDS and Machine Learning to Unlock High-Throughput Compositional Characterization of Powder Materials},
  year    = {2025},
  doi     = {10.21203/rs.3.rs-7837297/v1},
  url     = {https://doi.org/10.21203/rs.3.rs-7837297/v1}
}
```

---

## 📂 Paper Data

The raw data used in the associated publication is stored in the `paper_data/` directory.  
These files are tracked with **Git LFS** (Large File Storage).

### 🔽 Download with Git LFS
If you cloned the repository without Git LFS, you may only see placeholder text files instead of the actual datasets.  
To download the full data, on the terminal go to the repo directory and:

```bash
# 1. Install Git LFS (only needed once per machine)
git lfs install

# 2. Fetch the data files
git lfs pull
```

Alternatively, download manually from the github repo Download button.

After downloading, move the raw paper data into the Results/ folder to analyze it with AutoEMXSp, or add the folder's path to 'results_dir' within the provided analysis and quantification scripts.

---

## 📬 Contact

For questions or issues, please open an issue on GitHub.


---

