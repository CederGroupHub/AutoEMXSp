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
RAW_SPECTRUM_EXT : str
    File extension used when writing per-spectrum raw data pointer files.
    Default is `'.msa'` (EMSA/MAS format).
saved_images_extension : str
    Image extension used when saving SEM frames.
    Default is `'png'`.
save_raw_images : bool
    Whether to persist the non-annotated SEM image alongside the annotated one.
    Default is `False`.

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

substrate_els = ['C', 'O', 'Al']

use_instrument_background: bool = False

RAW_SPECTRUM_EXT: str = '.msa'

saved_images_extension: str = 'png'

save_raw_images: bool = False

# These values are used as a general estimate but they should be properly defined in the microscope calibration module XS_calibrations.py
escape_peak_probability = 0.03
pileup_peak_probability = 0.003
weight_Ll_ref_Ka1 = 0.05