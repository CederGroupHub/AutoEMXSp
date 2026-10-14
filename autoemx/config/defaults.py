from typing import Tuple

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoEMX default values to define microscope environment.

Defined Variables
------------------
microscope_ID : str
    Identifier for the microscope hardware.
    Must correspond to a calibration folder in `autoemx/calibrations/<ID>` (MicroscopeConfig.ID).
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
min_bckgrnd_cnts : float
    Minimum background counts under the reference peaks: spectra below are flagged (quant_flag 8). A filter of the
    analysis, re-evaluated from the stored background counts of each spectrum. Default is `5.0`
    (ClusteringConfig.min_bckgrnd_cnts).

Quality checks of spectra and quantifications (quant_flag values, see the user documentation page "Quantification
flags"). Fixed thresholds, not settings.

prefit_bckgrnd_max_energy_keV : float
    quant_flag 3: the background is checked on the spectrum below this energy. Default is `2` keV.
prefit_min_bckgrnd_cnts : float
    quant_flag 3: minimum background counts of the spectrum below prefit_bckgrnd_max_energy_keV. Default is `5.0`.
prefit_bckgrnd_n_points : int
    quant_flag 3: number of lowest (smoothed) spectrum points that must all be below prefit_min_bckgrnd_cnts.
    Default is `20`.
prefit_bckgrnd_smoothing_channels : int
    quant_flag 3: width of the moving average applied to the spectrum before the check, in channels. Default is `3`.
max_redchi_percent_total_counts : float
    quant_flag 4 (poor fit): maximum reduced chi-squared of the fit, in % of the total spectrum counts. Default is `0.16`.
max_quant_analytical_error : float
    quant_flag 5: maximum analytical error (deviation of the total mass fraction from 1). Default is `0.5` (50 w%).
max_absorption_increase : float
    quant_flag 6: maximum increase of the X-ray absorption around the reference peaks of particles, relative to a
    bulk sample. Default is `0.7` (absorption up to 170% of the bulk one).
max_substrate_peak_percent_total_counts : float
    quant_flag 7: maximum intensity of a substrate element peak, in % of the total spectrum counts. Default is `10`.
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
    Whether to also save the non-annotated version of annotated SEM images (when
    ``annotate_particle_images`` is True). Default is `False`.
annotate_particle_images : bool
    Whether to annotate the images saved during acquisition: X-ray spectrum spots on particle images,
    particles found on frame images, and a scale bar. If False, raw images are saved; annotated copies
    can be made afterwards with ``Annotate_Particle_Images.py``. Default is `False`.
Default values of runner parameters (also used by the configuration classes with the same default).

sample_halfwidth_mm : float
    Half-width of the sample area, in mm (SampleConfig.half_width_mm). Default is `2.9`.
substrate_width_mm : float
    Lateral dimension of the substrate holder (stub), in mm (SampleSubstrateConfig.stub_w_mm). Default is `12`.
working_distance_mm : float
    Approximate working distance at which the sample is in focus, in mm. Default is `5`.
working_distance_tolerance_mm : float
    Allowed deviation from the working distance during autofocus, in mm
    (MeasurementConfig.working_distance_tolerance). Default is `1`.
beam_energy_keV : float
    Electron beam energy, in keV (MeasurementConfig.beam_energy_keV). Default is `15.0`.
target_Xsp_counts : int
    Target counts of each spectrum (MeasurementConfig.target_acquisition_counts). Default is `50000`.
exp_std_target_Xsp_counts : int
    Target counts of each spectrum of experimental standards. Default is `250000`.
min_n_spectra : int
    Number of spectra after which acquisition runners check the convergence of the analysis. Default is `50`.
max_n_spectra : int
    Maximum number of spectra collected per sample (MeasurementConfig.max_n_spectra). Default is `100`.
max_analytical_error_percent : float
    Maximum analytical error, in w%, of the compositions used for clustering
    (ClusteringConfig.max_analytical_error_percent). Default is `5.0`.
max_n_clusters : int
    Maximum number of clusters (ClusteringConfig.max_k). Default is `6`.
fit_tolerance : float
    lmfit tolerance for fit convergence (QuantificationOptionsConfig.fit_tolerance). Default is `1e-4`.
max_undetectable_w_fr : float
    Maximum total mass fraction of undetectable elements, in single-spectrum fits. Default is `0`.
particle_stats_contrast, particle_stats_brightness : float
    SEM contrast and brightness used for particle statistics when automatic adjustment is off (PhenomXL values).
    Defaults are `4.3877` and `0.4504`.
peak_weights_sum_cv_threshold : float
    Peak-weight calibration: maximum coefficient of variation of the summed line weights. Default is `0.12`.
peak_weights_n_sig_figs : int
    Peak-weight calibration: significant figures of the calibrated weights. Default is `5`.
homogeneous_cluster_max_rms_dist : float
    Clusters whose RMS distance of the spectra from the centroid, in atomic fractions, is below this value are
    reported as likely homogeneous in Clusters.csv ("Homogeneous": "Likely"). Default is `0.03` (3 at%).
inhomogeneous_cluster_min_rms_dist : float
    Clusters whose RMS distance, in atomic fractions, is above this value are reported as unlikely to be
    homogeneous ("Unlikely"); between the two thresholds, "Unsure". Default is `0.05` (5 at%).
homogeneity_min_n_points : int
    Clusters with fewer spectra than this are too small to judge homogeneity ("Not Applicable"). Default is `5`.

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

min_bckgrnd_cnts: float = 5.0

# Quality checks of spectra and quantifications (fixed thresholds, see the module docstring)
prefit_bckgrnd_max_energy_keV: float = 2
prefit_min_bckgrnd_cnts: float = 5.0
prefit_bckgrnd_n_points: int = 20
prefit_bckgrnd_smoothing_channels: int = 3
max_redchi_percent_total_counts: float = 0.16
max_quant_analytical_error: float = 0.5
max_absorption_increase: float = 0.7
max_substrate_peak_percent_total_counts: float = 10

identify_missing_elements: bool = True

identification_seed_spectra: int = 10

identification_z_range: tuple = (5, 83)

identification_excluded_elements: tuple = ('Tc', 'Pm', 'Ne', 'Kr', 'Xe', 'Rn')

RAW_SPECTRUM_EXT: str = '.msa'

saved_images_extension: str = 'png'

save_raw_images: bool = False

annotate_particle_images: bool = False

# Default values of runner parameters (see the module docstring)
sample_halfwidth_mm: float = 2.9
substrate_width_mm: float = 12
working_distance_mm: float = 5
working_distance_tolerance_mm: float = 1
beam_energy_keV: float = 15.0
target_Xsp_counts: int = 50000
exp_std_target_Xsp_counts: int = 250000
min_n_spectra: int = 50
max_n_spectra: int = 100
max_analytical_error_percent: float = 5.0
max_n_clusters: int = 6
fit_tolerance: float = 1e-4
max_undetectable_w_fr: float = 0
particle_stats_contrast: float = 4.3877
particle_stats_brightness: float = 0.4504
peak_weights_sum_cv_threshold: float = 0.12
peak_weights_n_sig_figs: int = 5

homogeneous_cluster_max_rms_dist: float = 0.03

inhomogeneous_cluster_min_rms_dist: float = 0.05

homogeneity_min_n_points: int = 5

# These values are used as a general estimate but they should be properly defined in the detector model of the microscope (detector_model.json)
escape_peak_probability = 0.03
pileup_peak_probability = 0.003
weight_Ll_ref_Ka1 = 0.05