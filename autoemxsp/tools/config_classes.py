#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoEMXSp configuration dataclasses.

Created on Mon Jul 28 10:33:35 2025

@author: Andrea

This module provides configuration dataclasses for all stages of an automated X-ray spectroscopy workflow,
including microscope setup, sample and substrate definition, measurement and acquisition settings,
spectrum fitting, quantification and filtering, powder measurement, clustering, and plotting.

Configurations:

- MicroscopeConfig: Settings for microscope hardware, calibration, and imaging parameters.
- SampleConfig: Defines the sample’s identity, elements, and spatial properties.
- SampleSubstrateConfig: Specifies the substrate composition and geometry supporting the sample.
- MeasurementConfig: Controls measurement type, beam parameters, and acquisition settings.
- QuantConfig: Options for spectral fitting and quantification.
- PowderMeasurementConfig: Settings for analyzing powder samples and particle selection.
- ClusteringConfig: Configures clustering algorithms and filtering of X-ray spectra.
- PlotConfig: Options for saving, displaying, and customizing plots.

Each dataclass includes attribute documentation and input validation.
"""
import re
import numpy as np
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple, Dict

from pymatgen.core.periodic_table import Element
from pymatgen.core import Composition

import autoemxsp.tools.constants as cnst
import autoemxsp.core.particle_segmentation_models as par_seg_models

@dataclass
class MicroscopeConfig:
    """
    Configuration for the microscope hardware.

    Attributes:
        ID (str): Identifier for the microscope, defining instrument calibrations at ./XSp_calibs/Microscopes/ID.
        type (str): Type of microscope. Allowed: 'SEM' (implemented), 'STEM' (not implemented).
        is_auto_BC (bool): If True, brightness/contrast are set automatically.
        brightness (Optional[float]): Manual brightness value; required if is_auto_BC is False.
        contrast (Optional[float]): Manual contrast value; required if is_auto_BC is False.
        energy_zero (Optional[float]): Set from detector calibration files during spectral collection.
        bin_width (Optional[float]): Set from detector calibration files during spectral collection.

    Notes:
        - If `is_auto_BC` is False, both `brightness` and `contrast` must be provided.
        - STEM mode is not implemented and will raise NotImplementedError.
        - The microscope ID must correspond to a folder at ./XSp_calibs/Microscopes/ID containing all necessary calibration files.
    """

    ID: str
    type: str
    detector_type: str = 'BSD'
    is_auto_BC: bool = True
    brightness: Optional[float] = None
    contrast: Optional[float] = None
    energy_zero: Optional[float] = None
    bin_width: Optional[float] = None

    ALLOWED_TYPES = ("SEM", "STEM")
    ALLOWED_DETECTOR_TYPES = ("BSD")

    def __post_init__(self) -> None:
        import os
        from pathlib import Path

        # Check that the calibration folder for this microscope exists
        parent_dir = str(Path(__file__).resolve().parent.parent)
        calib_path = os.path.join(parent_dir, cnst.CALIBS_DIR, cnst.MICROSCOPES_CALIBS_DIR, self.ID)
        if not os.path.isdir(calib_path):
            raise FileNotFoundError(
                f"Calibration folder for microscope ID '{self.ID}' not found at '{calib_path}'.\n"
                "Please add all necessary calibration files and ensure the folder is named with the same ID."
            )

        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"Microscope type must be one of {self.ALLOWED_TYPES}, got '{self.type}'.")

        if self.type == "STEM":
            # STEM mode is not supported yet.
            raise NotImplementedError("STEM mode is not implemented yet.")
            
        if self.detector_type not in self.ALLOWED_DETECTOR_TYPES:
            raise ValueError(f"Detector type must be one of {self.ALLOWED_DETECTOR_TYPES}, got '{self.detector_type}'.")
                
        if not self.is_auto_BC:
            if self.brightness is None or self.contrast is None:
                raise ValueError(
                    "If is_auto_BC is False, both brightness and contrast must be provided."
                )


@dataclass
class SampleConfig:
    """
    Configuration for the sample.

    Attributes:
        ID (str): Identifier for the sample.
        elements (List[str]): List of elemental symbols (e.g., ['Fe', 'O']).
        type (str): Sample type. Allowed: 'powder' (implemented), 'bulk', 'film' (not implemented).
        w_frs (Dict[str,float]): Dict of elemental mass fractions to be kept fixed (e.g., {'Fe': 0.4, 'O': 0.6}).
            Normally not used
        center_pos (Tuple[float, float]): (x, y) center position of the sample on the stage, in mm.
        half_width_mm (float): Half-width of the sample in millimeters.

    Notes:
        - Only 'powder' and 'bulk' type are implemented. 'film' will raise NotImplementedError.
        - Element symbols are validated. An error is raised if any symbol is unrecognized.
    """

    ID: str
    elements: List[str]
    type: str = cnst.S_POWDER_SAMPLE_TYPE
    w_frs: Dict[str, float] = None
    center_pos: Tuple[float, float] = (0.0, 0.0) # in mm
    half_width_mm: float = 2.9  # in mm

    ALLOWED_TYPES = (cnst.S_POWDER_SAMPLE_TYPE, cnst.S_BULK_SAMPLE_TYPE , cnst.S_FILM_SAMPLE_TYPE)

    def __post_init__(self) -> None:
        # Validate and clean sample ID
        self.ID = self._clean_ID(self.ID)
        
        # Validate sample type
        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"Sample type must be one of {self.ALLOWED_TYPES}, got '{self.type}'.")

        if self.type in (cnst.S_FILM_SAMPLE_TYPE):
            raise NotImplementedError(f"Sample type '{self.type}' is not implemented yet.")

        # Validate element symbols using pymatgen
        for symbol in self.elements:
            try:
                Element(symbol)
            except Exception:
                raise ValueError(f"Element symbol '{symbol}' is not a recognized element.")
    
    @staticmethod
    def _clean_ID(ID: str) -> str:
        """Remove trailing whitespace and invisible characters from the ID to avoid errors in output saving."""
        cleaned_ID = ID.rstrip()
        if cleaned_ID != ID:
            print(f"Warning: ID '{ID}' contained trailing whitespace or invisible characters. Using cleaned ID: '{cleaned_ID}'")
        # Remove leading and trailing invisible characters
        cleaned_ID = re.sub(r'^\W+|\W+$', '', cleaned_ID)
        return cleaned_ID


@dataclass
class SampleSubstrateConfig:
    """
    Configuration for the sample substrate.

    Attributes:
        elements (List[str]): List of element symbols present in the sample substrate.
        type (str): Type of the sample substrate. Allowed values: 'Ctape'.
        shape (str): Shape of the sample substrate. Allowed values: 'circle', 'rectangle'.
        auto_detection (bool): Whether to attempt automatic detection of substrate. (implemented only for type = Ctape & shape = 'circle')
        stub_w_mm (float): Lateral dimension of substrate holder in mm, used for determining image size for auto_detection.

    Notes:
        - Element symbols are validated. An error is raised if any symbol is unrecognized.
    """
    elements: List[str] = field(default_factory=lambda: ['C', 'O', 'Al'])
    type: str = cnst.CTAPE_SUBSTRATE_TYPE
    shape: str = cnst.CIRCLE_SUBSTRATE_SHAPE
    auto_detection: bool = True
    stub_w_mm: float = 12

    ALLOWED_TYPES = (cnst.CTAPE_SUBSTRATE_TYPE, cnst.NONE_SUBSTRATE_TYPE)
    ALLOWED_SHAPES = (cnst.CIRCLE_SUBSTRATE_SHAPE, cnst.SQUARE_SUBSTRATE_SHAPE)
    ALLOWED_AUTO_DETECTION_TYPES = (cnst.CTAPE_SUBSTRATE_TYPE)

    def __post_init__(self):
        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"SampleSubstrate type must be one of {self.ALLOWED_TYPES}")
        if self.shape not in self.ALLOWED_SHAPES:
            raise ValueError(f"SampleSubstrate shape must be one of {self.ALLOWED_SHAPES}")
        if self.auto_detection and self.type != cnst.CTAPE_SUBSTRATE_TYPE:
            raise NotImplementedError(f"auto_detection is only implemented for types {self.ALLOWED_AUTO_DETECTION_TYPES}.")
        # Validate element symbols using pymatgen
        for symbol in self.elements:
            try:
                Element(symbol)
            except Exception:
                raise ValueError(f"Element symbol '{symbol}' is not a recognized element.")


@dataclass
class MeasurementConfig:
    """
    Configuration for the measurement/acquisition.

    Attributes:
        type (str): Measurement type. Allowed: 'EDS' (implemented), 'WDS' (not implemented).
        mode (str): Measurement mode (e.g., 'point'). Defines set of measurement parameters (i.e., beam current), determining detector calibration parameters
        working_distance (Optional[float]): Working distance to use for current measurement, in mm. Takes it from EM_driver if left unspecified. 
        beam_energy_keV (float): Electron beam energy in keV.
        beam_current (Optional[float]): Beam current; must be provided at initialization or via detector channel calibration file.
        emergence_angle (Optional[float]): Emergence angle; updated from microscope driver file if not provided.
        is_manual_navigation (bool): If True, instrument navigation is performed manually.
        max_acquisition_time (float): Maximum X-ray spectral acquisition time in seconds.
        target_acquisition_counts (int): Target number of counts for acquisition of X-ray spectrum.
        min_n_spectra (int): Minimum number of spectra to acquire.
        max_n_spectra (int): Maximum number of spectra to acquire.

    Notes:
        - Only 'EDS' type is implemented. 'WDS' will raise NotImplementedError.
        - If `beam_current` or `emergence_angle` are not provided, they should be set via calibration or microscope driver.
    """

    type: str = "EDS"
    mode: str = "point"
    working_distance: float = None # mm
    beam_energy_keV: float = 15.0  # in keV
    beam_current: Optional[float] = None  # Provide at initialization or via calibration
    emergence_angle: Optional[float] = None  # Updated from microscope driver if not provided
    is_manual_navigation: bool = False
    max_acquisition_time: float = 30.0  # seconds
    target_acquisition_counts: int = 50000
    min_n_spectra: int = 30
    max_n_spectra: int = 100
    
    PARTICLE_STATS_MEAS_TYPE_KEY = "particle_stats"
    ALLOWED_TYPES = ("EDS", "WDS", PARTICLE_STATS_MEAS_TYPE_KEY)

    def __post_init__(self) -> None:
        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"Measurement type must be one of {self.ALLOWED_TYPES}, got '{self.type}'.")

        if self.type == "WDS":
            raise NotImplementedError("WDS measurement type is not implemented yet.")
            
            
@dataclass
class QuantConfig:
    """
    Configuration for X-ray spectrum fitting and quantification.

    Attributes:
        method (str): Method to use for quantification. Currently only accepts 'PB'
        spectrum_lims (Tuple[float, float]): Lower and upper spectral index limits.
        fit_tolerance (float): lmfit tolerance for fit convergence
        use_instrument_background (bool): Whether to use the instrument background in the fit (Default: False).
            If False, AutoEMXSp computes the background while fitting.
        interrupt_fits_bad_spectra (bool): If True, fitting will stop early for spectra identified as poor quality.
        min_bckgrnd_cnts (Optional[int]): Minimum background counts required for spectrum not to be filtered out. Can be None.
    """
    method: str = 'PB'
    spectrum_lims: Tuple[float, float] = (14, 1100)
    fit_tolerance: float = 1e-3
    use_instrument_background: bool = False
    interrupt_fits_bad_spectra: bool = True
    min_bckgrnd_cnts: Optional[int] = 5  # Can be None
    
    ALLOWED_METHODS = ['PB']
    
    def __post_init__(self) -> None:
        if self.method not in self.ALLOWED_METHODS:
            raise ValueError(f"Quantification method must be one of {self.ALLOWED_METHODS}, got '{self.method}'."
                             "Currently no other method is implemented.")


@dataclass
class PowderMeasurementConfig:
    """
    Configuration for powder measurement.

    Attributes:
        is_manual_particle_selection (bool): Whether to manually navigate sample to select particles to analyse (Default = False).
        is_known_powder_mixture_meas (bool): Whether sample is a known binary mixture of powders. Used to characterize precursor extent of intermixing (Default = False).
        par_search_frame_width_um (float, optional): Frame width used when searching for particles, in um.
            Default: min(20*max_par_radius, 500 um)
        max_n_par_per_frame (int): Maximum number of particles analyzed in a single frame. 
            Used to ensure spatial representation of the analyzed sample.
        max_spectra_per_par (int): Maximum number of spot X-ray spectra collected in a single particle.
            Limiting this ensures more particles are analyzed.
        max_area_par (float): Maximum area (in µm²) for a particle to be considered.
        min_area_par (float): Minimum area (in µm²) for a particle to be considered.
        par_mask_margin (float): Margin (in µm) from particle edge where X-ray spectra should not be collected.
        xsp_spots_distance_um (float): Min distance between X-ray spectrum acquisition points
        par_segmentation_model (str) : Model to use for particle segmentation. Default: "threshold_bright"
        par_brightness_thresh (int): Intensity threshold in 8-bit image that defines a particle over a dark background.
        par_xy_spots_thresh (int): Intensity threshold in 8-bit image that defines bright (i.e., thickest) regions in particles.
            X-ray spectra are acquired only from these regions.
            Particle pixel intensities are scaled to 8-bit prior threhsolding, i.e., darkest pixel will be set to 0, and brightest to 255.
        par_feature_selection (str): 'random' for random selection of points within bright regions, 'peaks' for brightest peak spots (default: 'random').
        par_spot_spacing (str): 'random' for unbiased spot selecton, 'maximized' for maximized spot spacing over particle (default: 'random').
    """
    DEFAULT_PAR_SEGMENTATION_MODEL = "threshold_bright"

    is_manual_particle_selection: bool = False
    is_known_powder_mixture_meas: bool = False
    par_search_frame_width_um: float = None 
    max_n_par_per_frame: int = 30
    max_spectra_per_par: int = 3
    max_area_par: float = 300.0    # µm²
    min_area_par: float = 10.0     # µm²
    par_mask_margin: float = 1.0   # µm
    xsp_spots_distance_um: float = 1.0 # µm
    par_segmentation_model : str =  DEFAULT_PAR_SEGMENTATION_MODEL
    par_brightness_thresh: int = 100 # in 8-bit image
    par_xy_spots_thresh: int = 100  # considering particle pixel intensities are scaled to 8-bit image
    par_feature_selection: str = 'random'
    par_spot_spacing: str = 'random'
    
    AVAILABLE_PAR_SEGMENTATION_MODELS = [DEFAULT_PAR_SEGMENTATION_MODEL] + par_seg_models.AVAILABLE_SEGMENTATION_MODELS
    AVAILABLE_FEATURE_SELECTION = ('random', 'peaks')
    AVAILABLE_SPOT_SPACING_SELECTION = ('random', 'maximized')
    
    def __post_init__(self):
        # --- 0. Check validity of passed variables
        if self.par_segmentation_model not in self.AVAILABLE_PAR_SEGMENTATION_MODELS:
            raise ValueError(
                f'Value of "par_segmentation_model" set to {self.par_segmentation_model} is invalid. '
                f'Must be one of {self.AVAILABLE_PAR_SEGMENTATION_MODELS}.'
            )
        if self.par_feature_selection not in self.AVAILABLE_FEATURE_SELECTION:
            raise ValueError(
                f'Value of "par_feature_selection" set to {self.par_feature_selection} is invalid. '
                f'Must be one of {self.available_feature_selection}.'
            )
        if self.par_spot_spacing not in self.AVAILABLE_SPOT_SPACING_SELECTION:
            raise ValueError(
                f'Value of "par_spot_spacing" set to {self.par_spot_spacing} is invalid. '
                f'Must be one of {self.available_spot_spacing}.'
            )
        # Additional checks can be added here (e.g., for numeric bounds)
        if self.min_area_par < 0 or self.max_area_par < 0:
            raise ValueError("Particle area thresholds must be non-negative.")
        if self.max_area_par < self.min_area_par:
            raise ValueError("max_area_par must be greater than or equal to min_area_par.")
        if self.max_n_par_per_frame <= 0:
            raise ValueError("max_n_par_per_frame must be positive.")
        if self.max_spectra_per_par <= 0:
            raise ValueError("max_spectra_per_par must be positive.")
        if self.par_mask_margin < 0:
            raise ValueError("par_mask_margin must be non-negative.")
        if not (0 <= self.par_brightness_thresh <= 255):
            raise ValueError("par_brightness_thresh must be in 0..255.")
        if not (0 <= self.par_xy_spots_thresh <= 255):
            raise ValueError("par_xy_spots_thresh must be in 0..255.")
            
        # --- 1. Define default par_search_frame_width_um if None
        if self.par_search_frame_width_um is None:
            max_par_radius = np.sqrt(self.max_area_par / np.pi)  # in µm
            self.par_search_frame_width_um = min(20 * max_par_radius, 500.0)  # µm
 

@dataclass
class BulkMeasurementConfig:
    """
    Configuration for characterization or bulk-like samples.

    Attributes
    ----------
    grid_spot_spacing_um : float
        Distance between grid points to measure, in micrometers (µm).
    min_xsp_spots_distance_um : float
        Offset distance for acquisition spot grid if the original grid
        does not contain enough spots to measure the required number
        of spectra, in micrometers (µm).
    image_frame_width_um : float, optional
        Width of the image frame in micrometers (µm). If not specified,
        defaults to 10 × grid_spot_spacing_um.
    randomize_frames : bool
        Whether to randomize the order of spectra acquisition in the constructed grid.
    exclude_sample_margin : bool
        Whether to exclude the margin of the sample (useful if contaminated).
    """
    grid_spot_spacing_um: float = 100.0  # µm
    min_xsp_spots_distance_um: float = 5.0  # µm
    image_frame_width_um: float = None # µm
    randomize_frames: bool = False
    exclude_sample_margin: bool = False

    def __post_init__(self):
        # Validate grid spot spacing
        if not (self.grid_spot_spacing_um > 0):
            raise ValueError("grid_spot_spacing_um must be positive.")

        # Validate minimum spot distance
        if not (self.min_xsp_spots_distance_um > 0):
            raise ValueError("min_xsp_spots_distance_um must be positive.")

        if self.min_xsp_spots_distance_um > self.grid_spot_spacing_um:
            raise ValueError(
                "min_xsp_spots_distance_um should not exceed grid_spot_spacing_um."
            )

        # Set default image frame width if unspecified
        if self.image_frame_width_um is None:
            self.image_frame_width_um = 10 * self.grid_spot_spacing_um
 

@dataclass
class ExpStandardsConfig:
    """
    Configuration for the collection of experimental standards.

    Attributes:
        is_exp_std_measurement (bool): 
            Whether the configuration corresponds to the measurement of an experimental standard (Default = False)
            If True, a valid `formula` must be provided and weight fractions will be automatically calculated.

        formula (str): 
            Chemical formula of the experimental standard. Required if `is_exp_std_measurement` is True.
            Must be parseable by `pymatgen.core.Composition`.
        
        use_for_mean_PB_calc (bool):
            Whether the acquired experimental standards should be used to calculate the average PB, which is the 
            reference standard value employed generally during spectral quantification (Default = True).
            This should be set to False when collecting powder standards for quantifying the extent of intermixing
            in powder standards.
        
        min_acceptable_PB_ratio (float): 
            Minimum PB ratio required for a peak to be accepted as a standard. in cnts/cnts*keV^-1 (Deafult = 10).

        quant_flags_accepted (List[int]): 
            List of quantification flags considered acceptable. Other spectra are filtered out before clustering.
            Quantification flags indicate whether the quantification or the fit of each spectrum is likely to be 
            affected by large errors:
                0  : Quantification is ok, although it may be affected by large analytical error.
               -1  : As above, but quantification did not converge within 30 steps.
                1  : Error during EDS acquisition. No fit executed.
                2  : Total counts < 95% of target counts, likely due to wrong segmentation. No fit executed.
                3  : Too little low-energy signal, causing poor quantification in that region. No fit executed.
                4  : Poor fit. Fit interrupted if interrupt_fits_bad_spectra=True.
                5  : High analytical error (>50%), possibly due to missing element or other major error. Fit interrupted if interrupt_fits_bad_spectra=True.
                6  : Excessive X-ray absorption. Fit interrupted if interrupt_fits_bad_spectra=True.
                7  : Excessive contamination from substrate.
                8  : Too few background counts below reference peak, likely leading to large quantification errors.
                9  : Unknown fitting error.
                10 : (Only for measurement of experimental standards) Reference peak missing.

        w_frs (Optional[Dict[str, float]]): 
            Dictionary of element symbols and their corresponding weight fractions (computed via pymatgen)
            if `is_exp_std_measurement` is True and `formula` is valid; otherwise None.

    Raises:
        ValueError: 
            If `is_exp_std_measurement` is True but `formula` is missing or invalid.
    """

    
    is_exp_std_measurement: bool = False
    formula: str = ''
    use_for_mean_PB_calc: bool = True
    min_acceptable_PB_ratio: float = 10
    quant_flags_accepted: List[int] = field(default_factory=lambda: [0])
    w_frs: Optional[Dict[str, float]] = None  # Will hold calculated weight fractions

    def __post_init__(self) -> None:
        if self.is_exp_std_measurement:
            if not self.formula:
                raise ValueError("Formula must be provided when is_exp_std_measurement is True.")
            try:
                comp = Composition(self.formula)
                # Convert FloatWithUnit to plain float
                self.w_frs = {el: float(w) for el, w in comp.to_weight_dict.items()}
            except Exception as e:
                raise ValueError(f"Invalid chemical formula '{self.formula}': {e}")
        
        

@dataclass
class ClusteringConfig:
    """
    Configuration for clustering of compositions and their filtering.

    Attributes:
        method (str): Clustering algorithm to use. Allowed: 'kmeans' (implemented), 'dbscan' (not implemented).
        features (List[Any]): Feature set to use for clustering.
        k (Optional[int]): If provided, defines a fixed number of clusters.
        k_finding_method (str): Method to determine the number of clusters. Set to "forced" if a value of 'k' is specified manually.
            Allowed methods are "silhouette", "calinski_harabasz", "elbow".
        max_k (int): Maximum allowed number of clusters.
        ref_formulae (List[str]): List of possible phases present in the sample, as chemical formula strings.
        max_analytical_error_percent (float): Maximum analytical error acceptable for composition to be considered in phase determination, expressed as w%. Can be None.
        quant_flags_accepted (List[int]): List of quantification flags considered acceptable, others are filtered out prior clustering.
            Quantification flags indicate whether the quantification or the fit of each spectrum is likely to be affected by large errors:
                0: Quantification is ok, although it may be affected by large analytical error
               -1: As above, but quantification did not converge within 30 steps
                1: Error during EDS acquisition. No fit executed
                2: Total number of counts is lower than 95% of target counts, likely due to wrong segmentation. No fit executed
                3: Spectrum has too low signal in its low-energy portion, leading to poor quantification in this region. No fit executed
                4: Poor fit. Fit interrupted if interrupt_fits_bad_spectra=True
                5: Too high analytical error (>50%) indicating a missing element or other major sources of error. Fit interrupted if interrupt_fits_bad_spectra=True
                6: Excessive X-ray absorption. Fit interrupted if interrupt_fits_bad_spectra=True
                7: Excessive signal contamination from substrate
                8: Too few background counts below reference peak, likely leading to large quantification errors
                9: Unknown fitting error
    """
    method: str = 'kmeans'
    features: List[Any] = field(default_factory=lambda: cnst.AT_FR_CL_FEAT)
    k: Optional[int] = None
    DEFAULT_K_FINDING_METHOD = 'silhouette'
    k_finding_method: str = DEFAULT_K_FINDING_METHOD
    max_k: int = 6
    ref_formulae: List[str] = field(default_factory=list)
    max_analytical_error_percent: float = 5  # w%, Can be None
    quant_flags_accepted: List[int] = field(default_factory=lambda: [0, -1])
    
    FORCED_K_METHOD_KEY = 'forced'
    ALLOWED_METHODS = ("kmeans", "dbscan")
    ALLOWED_FEATURE_SETS = (cnst.W_FR_CL_FEAT, cnst.AT_FR_CL_FEAT)
    ALLOWED_K_FINDING_METHODS = ("silhouette", "calinski_harabasz", "elbow", FORCED_K_METHOD_KEY)

    def __post_init__(self):    
        if self.method not in self.ALLOWED_METHODS:
            raise ValueError(f"Clustering method must be one of {self.ALLOWED_METHODS}, got '{self.method}'.")
        if self.method == "dbscan":
            raise NotImplementedError("DBSCAN clustering is not implemented yet.")
        if not any(self.features == allowed for allowed in self.ALLOWED_FEATURE_SETS):
            raise ValueError(
                f"Invalid value for features: {self.features}. "
                f"Expected one of: {self.ALLOWED_FEATURE_SETS}."
            )
        if self.k_finding_method not in self.ALLOWED_K_FINDING_METHODS:
            raise ValueError(
                f"k_finding_method must be one of {self.ALLOWED_K_FINDING_METHODS}, "
                f"got '{self.k_finding_method}'."
            )
        elif isinstance(self.k, int):
            self.k_finding_method = self.FORCED_K_METHOD_KEY
        elif self.k_finding_method == self.FORCED_K_METHOD_KEY:
            raise ValueError(
                f"'k_finding_method' should not be set to {self.FORCED_K_METHOD_KEY} if "
                f"'k' is left unspecified. Setting 'k_finding_method = {self.DEFAULT_K_FINDING_METHOD}'."
            )
            
    
    
@dataclass
class PlotConfig:
    """
    Configuration for plotting.

    Attributes:
        show_unused_comps_clust (bool): Whether to plot unused data points in clustering plot.
        els_excluded_clust_plot (List[str]): Elements to exclude in cluster plot when more than 3 elements are present.
        save_plots (bool): Whether to save plots to disk.
        show_plots (bool): Whether to display plots interactively.
        use_custom_plots (bool): Whether to use custom plotting routines.
    """
    show_unused_comps_clust: bool = True
    els_excluded_clust_plot: List[str] = field(default_factory=list)
    save_plots: bool = True
    show_plots: bool = False
    use_custom_plots: bool = False
    
    
# Dictionary of all dataclasses. Loaded for data import
config_classes_dict = {
    cnst.SAMPLE_CFG_KEY: SampleConfig,
    cnst.MICROSCOPE_CFG_KEY: MicroscopeConfig,
    cnst.MEASUREMENT_CFG_KEY: MeasurementConfig,
    cnst.SAMPLESUBSTRATE_CFG_KEY: SampleSubstrateConfig,
    cnst.QUANTIFICATION_CFG_KEY: QuantConfig,
    cnst.CLUSTERING_CFG_KEY: ClusteringConfig,
    cnst.PLOT_CFG_KEY: PlotConfig,
    cnst.POWDER_MEASUREMENT_CFG_KEY: PowderMeasurementConfig,
    cnst.BULK_MEASUREMENT_CFG_KEY: BulkMeasurementConfig,
    cnst.EXP_STD_MEASUREMENT_CFG_KEY: ExpStandardsConfig,
}