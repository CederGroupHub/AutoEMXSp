from typing import Tuple

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoEMX default values to define microscope environment.

Defined Variables
------------------
microscope_ID : str
    Identifier for the microscope hardware.
    Must correspond to a calibration folder in `./XSp_calibs/Microscopes/<ID>` (MicroscopeConfig.ID).
    Default is `'PhenomXL'`.
microscope_type : str
    Type of microscope. Allowed: `'SEM'` (implemented), `'STEM'` (not implemented).
    Default is `'SEM'` (MicroscopeConfig.type).
measurement_type : str
    Measurement type. Allowed: `'EDS'` (implemented), `'WDS'` (not implemented).
    Default is `'EDS'` (MeasurementConfig.type).
measurement_mode : str
    Acquisition mode (e.g., `'point'`, `'map'`), defining beam/detector calibration settings.
    Default is `'point'` (MeasurementConfig.mode).
detector_type : str
    Detector type to employ when navigating SEM and collecting images.
quantification_method : str
    Quantification method. Currently only `'PB'` (Phi-Rho-Z) is implemented.
    Default is `'PB'` (QuantificationOptionsConfig.method).
spectrum_lims : tuple of float
    Lower and upper energy limits for spectrum fitting in eV.
    Default is `(14, 1100)` (QuantificationOptionsConfig.spectrum_lims).
use_instrument_background : bool
    Whether to use instrument background files during fitting.
    If False, background is computed during fitting.
    Default is `False` (QuantificationOptionsConfig.use_instrument_background).
min_total_counts_fraction : float
    Minimum accepted total spectrum counts as a fraction of ``target_acquisition_counts``.
    Default is `0.9` (QuantificationOptionsConfig.min_total_counts_fraction).
identify_missing_elements : bool
    Whether to look for elements missing from the fit of each spectrum (added_elements) before quantifying it.
    Default is `True` (QuantificationOptionsConfig.identify_missing_elements). Quantification runs of ledgers saved
    before this option existed load with it off.
identification_seed_spectra : int
    Number of spectra checked without prior knowledge before the elements found in them are tested first
    in the remaining spectra. Default is `10` (QuantificationOptionsConfig.identification_seed_spectra).
identification_z_range : tuple of int
    Atomic numbers of the elements the identification of missing elements may propose. Default is `(5, 83)` (B to Bi).
identification_excluded_elements : tuple of str
    Elements the identification of missing elements never proposes. Default: Tc and Pm (no stable isotopes) and the noble
    gases Ne, Kr, Xe, Rn.
RAW_SPECTRUM_EXT : str
    Extension of spectrum files when no microscope driver is loaded (e.g. when AutoEMX writes spectrum files itself).
    During acquisition, the extension defined by the microscope driver (RAW_SPECTRUM_EXT, the format used by the EDS
    software to export spectra) is used instead. Default is `'.msa'` (EMSA/MAS format).
saved_images_extension : str
    Image extension used when saving SEM frames.
    Default is `'png'`.
save_raw_images : bool
    Whether to also save the non-annotated version of annotated SEM images (frame images, and
    particle images when ``annotate_particle_images`` is True). Default is `False`.
annotate_particle_images : bool
    Whether to draw the X-ray spectrum spots (dot + spectrum ID) and a scale bar on the particle
    images saved during acquisition. If False, raw particle images are saved; annotated copies can
    be made afterwards with ``Annotate_Particle_Images.py``. Default is `False`.

Created on Sun Dec 21 18:59:50 2025

@author: Andrea
"""

microscope_ID: str = 'PhenomXL'

microscope_type: str = 'SEM'

measurement_type: str = 'EDS'

measurement_mode: str = 'point'

detector_type: str = 'BSD'

quantification_method: str = 'PB'

spectrum_lims: Tuple[int, int] = (14, 1100)

substrate_els = ['C']  # O and Al (e.g. Al stub) are found by the identification of missing elements if present

substrate_els_without_identification = ['C', 'O', 'Al']  # Used when the identification of missing elements is off (e.g. experimental standards), since O and Al (e.g. Al stub) would otherwise be missed

use_instrument_background: bool = False

min_total_counts_fraction: float = 0.9

identify_missing_elements: bool = True

identification_seed_spectra: int = 10

identification_z_range: tuple = (5, 83)

identification_excluded_elements: tuple = ('Tc', 'Pm', 'Ne', 'Kr', 'Xe', 'Rn')

RAW_SPECTRUM_EXT: str = '.msa'

saved_images_extension: str = 'png'

save_raw_images: bool = False

annotate_particle_images: bool = False

# These values are used as a general estimate but they should be properly defined in the microscope calibration module XS_calibrations.py
escape_peak_probability = 0.03
pileup_peak_probability = 0.003
weight_Ll_ref_Ka1 = 0.05