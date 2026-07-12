#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AutoEMX configuration dataclasses.

Created on Mon Jul 28 10:33:35 2025

@author: Andrea

This module provides configuration dataclasses for all stages of an automated X-ray spectroscopy workflow,
including microscope setup, sample and substrate definition, measurement and acquisition settings,
spectrum fitting, quantification and filtering, powder measurement, and plotting.

Configurations:

- MicroscopeConfig: Settings for microscope hardware, calibration, and imaging parameters.
- SampleConfig: Defines the sample’s composition/type and spatial properties.
- SampleSubstrateConfig: Specifies the substrate composition and geometry supporting the sample.
- MeasurementConfig: Controls measurement type, beam parameters, and acquisition settings.
- QuantificationOptionsConfig: Runtime options for spectral fitting and quantification.
- PowderMeasurementConfig: Settings for analyzing powder samples and particle selection.
- BulkMeasurementConfig: Settings for analyzing non-powder samples.
- PlotConfig: Options for saving, displaying, and customizing plots.

Each dataclass includes attribute documentation and input validation.
"""
import numpy as np
import multiprocessing
from typing import Any, ClassVar, Dict, List, Literal, Optional, Tuple

from pydantic import BaseModel, ConfigDict, Field, model_validator

from pymatgen.core.periodic_table import Element # type: ignore
from pymatgen.core import Composition # type: ignore

import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
import autoemx.core.em_runtime.particle_segmentation_models as par_seg_models

class MicroscopeConfig(BaseModel):
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

    ID: str = dflt.microscope_ID
    type: str = dflt.microscope_type
    detector_type: str = dflt.detector_type
    is_auto_BC: bool = True
    brightness: Optional[float] = None
    contrast: Optional[float] = None
    energy_zero: Optional[float] = None
    bin_width: Optional[float] = None

    ALLOWED_TYPES: ClassVar[Tuple[str, ...]] = ("SEM", "STEM")
    ALLOWED_DETECTOR_TYPES: ClassVar[Tuple[str, ...]] = ("BSD",)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _validate(self) -> "MicroscopeConfig":
        import os
        from pathlib import Path

        parent_dir = str(Path(__file__).resolve().parent.parent)
        calib_path = os.path.join(parent_dir, cnst.CALIBS_DIR, self.ID)
        if not os.path.isdir(calib_path):
            raise FileNotFoundError(
                f"Calibration folder for microscope ID '{self.ID}' not found at '{calib_path}'.\n"
                "Please add all necessary calibration files and ensure the folder is named with the same ID."
            )

        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"Microscope type must be one of {self.ALLOWED_TYPES}, got '{self.type}'.")

        if self.type == "STEM":
            raise NotImplementedError("STEM mode is not implemented yet.")

        if self.detector_type not in self.ALLOWED_DETECTOR_TYPES:
            raise ValueError(f"Detector type must be one of {self.ALLOWED_DETECTOR_TYPES}, got '{self.detector_type}'.")

        if self.is_auto_BC:
            # Do not persist manual BC values when auto BC is enabled.
            self.brightness = None
            self.contrast = None
        else:
            if self.brightness is None or self.contrast is None:
                raise ValueError(
                    "If is_auto_BC is False, both brightness and contrast must be provided."
                )

        return self


class SampleConfig(BaseModel):
    """
    Configuration for the sample.

    Attributes:
        elements (List[str]): List of elemental symbols (e.g., ['Fe', 'O']).
        type (str): Sample type. Allowed types:
            - powder: Expects particles, and uses geometrical correction factors during spectral fits.
            - powder_continuous: Expects a quasi-continuous powder mixture, sampled in a grid (NO PARTICLE DETECTION APPLIED).
                    Applies geometrical correction factors during spectral fits.
            - bulk: Expects a flat, continuous surface, sampled in a grid. Does not apply geometrical factors.
            - bulk_rough: Expects a continuous surface, sampled in a grid. Applies geometrical factors.
            - film: NOT IMPLEMENTED YET
        w_frs (Dict[str,float]): Dict of elemental mass fractions to be kept fixed (e.g., {'Fe': 0.4, 'O': 0.6}).
            Normally not used
        center_pos (Tuple[float, float]): (x, y) center position of the sample on the stage, in mm.
        half_width_mm (float): Half-width of the sample in millimeters.

    Notes:
        - Only 'powder' and 'bulk' type are implemented. 'film' will raise NotImplementedError.
        - Element symbols are validated. An error is raised if any symbol is unrecognized.
    """

    elements: List[str]
    type: str = cnst.S_POWDER_SAMPLE_TYPE
    w_frs: Optional[Dict[str, float]] = None
    center_pos: Tuple[float, float] = (0.0, 0.0)  # in mm
    half_width_mm: float = 2.9  # in mm

    ALLOWED_TYPES: ClassVar[Tuple[str, ...]] = (
        cnst.S_POWDER_SAMPLE_TYPE,
        cnst.S_POWDER_CONTINUOUS_SAMPLE_TYPE,
        cnst.S_BULK_SAMPLE_TYPE,
        cnst.S_BULK_ROUGH_SAMPLE_TYPE,
        cnst.S_FILM_SAMPLE_TYPE,
    )

    POWDER_SAMPLES_TYPES: ClassVar[List[str]] = [cnst.S_POWDER_SAMPLE_TYPE, cnst.S_POWDER_CONTINUOUS_SAMPLE_TYPE]
    is_powder_sample: bool = False  # overwritten based on type

    ROUGH_SURFACE_TYPES: ClassVar[List[str]] = POWDER_SAMPLES_TYPES + [cnst.S_BULK_ROUGH_SAMPLE_TYPE]
    is_surface_rough: bool = False  # overwritten based on type

    GRID_ACQUISITION_TYPES: ClassVar[Tuple[str, ...]] = (
        cnst.S_BULK_SAMPLE_TYPE,
        cnst.S_POWDER_CONTINUOUS_SAMPLE_TYPE,
        cnst.S_BULK_ROUGH_SAMPLE_TYPE,
        cnst.S_FILM_SAMPLE_TYPE,
    )
    is_grid_acquisition: bool = False  # overwritten based on type

    PARTICLE_ACQUISITION_TYPES: ClassVar[str] = cnst.S_POWDER_SAMPLE_TYPE
    is_particle_acquisition: bool = False  # overwritten based on type

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _validate(self) -> "SampleConfig":
        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"Sample type must be one of {self.ALLOWED_TYPES}, got '{self.type}'.")

        if self.type in (cnst.S_FILM_SAMPLE_TYPE):
            raise NotImplementedError(f"Sample type '{self.type}' is not implemented yet.")

        for symbol in self.elements:
            try:
                Element(symbol)
            except Exception:
                raise ValueError(f"Element symbol '{symbol}' is not a recognized element.")

        self.is_powder_sample = self.type in self.POWDER_SAMPLES_TYPES
        self.is_surface_rough = self.type in self.ROUGH_SURFACE_TYPES
        self.is_grid_acquisition = self.type in self.GRID_ACQUISITION_TYPES
        self.is_particle_acquisition = self.type in self.PARTICLE_ACQUISITION_TYPES

        return self


class SampleSubstrateConfig(BaseModel):
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

    elements: List[str] = Field(default_factory=lambda: ['C', 'O', 'Al'])
    type: str = cnst.CTAPE_SUBSTRATE_TYPE
    shape: str = cnst.CIRCLE_SUBSTRATE_SHAPE
    auto_detection: bool = True
    stub_w_mm: float = 12

    ALLOWED_TYPES: ClassVar[Tuple[str, ...]] = (cnst.CTAPE_SUBSTRATE_TYPE, cnst.NONE_SUBSTRATE_TYPE)
    ALLOWED_SHAPES: ClassVar[Tuple[str, ...]] = (cnst.CIRCLE_SUBSTRATE_SHAPE, cnst.SQUARE_SUBSTRATE_SHAPE)
    ALLOWED_AUTO_DETECTION_TYPES: ClassVar[str] = cnst.CTAPE_SUBSTRATE_TYPE

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _validate(self) -> "SampleSubstrateConfig":
        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"SampleSubstrate type must be one of {self.ALLOWED_TYPES}")
        if self.shape not in self.ALLOWED_SHAPES:
            raise ValueError(f"SampleSubstrate shape must be one of {self.ALLOWED_SHAPES}")
        if self.auto_detection and self.type != cnst.CTAPE_SUBSTRATE_TYPE:
            raise NotImplementedError(f"auto_detection is only implemented for types {self.ALLOWED_AUTO_DETECTION_TYPES}.")
        for symbol in self.elements:
            try:
                Element(symbol)
            except Exception:
                raise ValueError(f"Element symbol '{symbol}' is not a recognized element.")
        return self


class MeasurementConfig(BaseModel):
    """
    Configuration for the measurement/acquisition session.

    Attributes:
        type (str): Measurement type. Allowed: 'EDS' (implemented), 'WDS' (not implemented).
        mode (str): Measurement mode (e.g., 'point'). Defines set of measurement parameters (i.e., beam current), determining detector calibration parameters
        working_distance (Optional[float]): Working distance to use for current measurement, in mm. Takes it from EM_driver if left unspecified.
        working_distance_tolerance (Optional[float]): Defines maximum accepted deviation of working distance from its typical value, in mm.
            Used to prevent gross mistakes from EM autofocus. Default: 1 mm.
        beam_energy_keV (float): Electron beam energy in keV.
        beam_current (Optional[float]): Beam current; must be provided at initialization or via detector channel calibration file.
        emergence_angle (Optional[float]): Emergence angle; updated from microscope driver file if not provided.
        is_manual_navigation (bool): If True, instrument navigation is performed manually.
        max_acquisition_time (float): Maximum X-ray spectral acquisition time in seconds.
        target_acquisition_counts (int): Target number of counts for acquisition of X-ray spectrum.
        min_n_spectra (int): Minimum number of spectra to acquire.
        max_n_spectra (int): Maximum number of spectra to acquire.
        powder_meas_cfg (Optional[PowderMeasurementConfig]): Powder acquisition settings.
        bulk_meas_cfg (Optional[BulkMeasurementConfig]): Bulk/grid acquisition settings.
        exp_stds_cfg (Optional[ExpStandardsConfig]): Experimental standards settings.
        saved_images_extension (str): Extension used when saving SEM frame images.
        save_raw_images (bool): Whether to save the non-annotated SEM image.
    """

    type: str = dflt.measurement_type
    mode: str = dflt.measurement_mode
    working_distance: Optional[float] = None  # mm
    working_distance_tolerance: float = 1  # mm
    beam_energy_keV: float = 15.0  # in keV
    beam_current: Optional[float] = None
    emergence_angle: Optional[float] = None
    is_manual_navigation: bool = False
    max_acquisition_time: float = 30.0  # seconds
    target_acquisition_counts: int = 50000
    min_n_spectra: int = 30
    max_n_spectra: int = 100
    powder_meas_cfg: Optional["PowderMeasurementConfig"] = Field(default_factory=lambda: PowderMeasurementConfig())
    bulk_meas_cfg: Optional["BulkMeasurementConfig"] = Field(default_factory=lambda: BulkMeasurementConfig())
    exp_stds_cfg: Optional["ExpStandardsConfig"] = None
    saved_images_extension: str = dflt.saved_images_extension
    save_raw_images: bool = dflt.save_raw_images
    ALLOWED_IMAGE_EXTENSIONS: ClassVar[Tuple[str, ...]] = (
        "tif", "tiff", "png", "jpg", "jpeg", "webp", "bmp",
    )

    PARTICLE_STATS_MEAS_TYPE_KEY: ClassVar[str] = "particle_stats"
    ALLOWED_TYPES: ClassVar[Tuple[str, ...]] = ("EDS", "WDS", "particle_stats")

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _validate(self) -> "MeasurementConfig":
        if self.type not in self.ALLOWED_TYPES:
            raise ValueError(f"Measurement type must be one of {self.ALLOWED_TYPES}, got '{self.type}'.")
        if self.type == "WDS":
            raise NotImplementedError("WDS measurement type is not implemented yet.")
        ext = str(self.saved_images_extension).strip().lower().lstrip(".")
        if ext not in self.ALLOWED_IMAGE_EXTENSIONS:
            raise ValueError(
                f"saved_images_extension must be one of {self.ALLOWED_IMAGE_EXTENSIONS}, got '{self.saved_images_extension}'."
            )
        self.saved_images_extension = ext
        return self


class QuantificationOptionsConfig(BaseModel):
    """
    Configuration for X-ray spectrum fitting and quantification.

    Attributes:
        method (str): Method to use for quantification. Currently only accepts 'PB'
        spectrum_lims (Tuple[float, float]): Lower and upper spectral index limits.
        fit_tolerance (float): lmfit tolerance for fit convergence
        use_instrument_background (bool): Whether to use the instrument background in the fit (Default: False).
            If False, AutoEMX computes the background while fitting.
        use_project_specific_std_dict (bool): If True, tries to load the dictionary of reference standards from the project folder.
            If not found, uses the default file "EDS_Stds_{beamenergy}keV.json" at XSp_calibs/Microscopes/your_microscope.
    """
    method: str = dflt.quantification_method
    spectrum_lims: Tuple[float, float] = dflt.spectrum_lims
    fit_tolerance: float = 1e-4
    use_instrument_background: bool = dflt.use_instrument_background
    use_project_specific_std_dict: bool = False

    ALLOWED_METHODS: ClassVar[List[str]] = ['PB']

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def _normalize_legacy_fields(cls, data: Any) -> Any:
        if isinstance(data, dict):
            cleaned = dict(data)
            cleaned.pop("interrupt_fits_bad_spectra", None)
            cleaned.pop("min_bckgrnd_cnts", None)
            cleaned.pop("is_particle", None)
            cleaned.pop("beam_energy_keV", None)
            cleaned.pop("emergence_angle", None)
            cleaned.pop("det_ch_offset", None)
            cleaned.pop("det_ch_width", None)
            return cleaned
        return data

    @model_validator(mode="after")
    def _validate(self) -> "QuantificationOptionsConfig":
        if self.method not in self.ALLOWED_METHODS:
            raise ValueError(
                f"Quantification method must be one of {self.ALLOWED_METHODS}, got '{self.method}'."
                "Currently no other method is implemented."
            )
        return self


class PowderMeasurementConfig(BaseModel):
    """
    Configuration for powder measurement.

    Attributes:
        par_selection_mode (str): 'auto' for automatic particle navigation, 'manual' to prompt user to center each particle (default: 'auto').
        is_known_powder_mixture_meas (bool): Whether sample is a known binary mixture of powders. Used to characterize precursor extent of intermixing (Default = False).
        img_shift_tracking (bool): Whether to use image shift tracking during acquisition (Default = True).
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
        par_spot_selection_mode (str): 'auto' for built-in spot selection, 'callback' to supply spots via xsp_spot_selector (default: 'auto').
    """
    DEFAULT_PAR_SEGMENTATION_MODEL: ClassVar[str] = "threshold_bright"
    AVAILABLE_SPOT_SELECTION_MODES: ClassVar[Tuple[str, ...]] = ('auto', 'callback')
    AVAILABLE_PAR_SELECTION: ClassVar[Tuple[str, ...]] = ('auto', 'manual')

    par_selection_mode: Literal['auto', 'manual'] = 'auto'
    is_known_powder_mixture_meas: bool = False
    img_shift_tracking: bool = True
    par_search_frame_width_um: Optional[float] = None
    max_n_par_per_frame: int = 30
    max_spectra_per_par: int = 3
    max_area_par: float = 300.0    # µm²
    min_area_par: float = 10.0     # µm²
    par_mask_margin: float = 1.0   # µm
    xsp_spots_distance_um: float = 1.0  # µm
    par_spot_selection_mode: Literal['auto', 'callback'] = 'auto'
    par_segmentation_model: str = DEFAULT_PAR_SEGMENTATION_MODEL    # "threshold_bright"    
    par_brightness_thresh: int = 100  # in 8-bit image
    par_xy_spots_thresh: int = 100
    par_feature_selection: str = 'random'
    par_spot_spacing: str = 'random'

    AVAILABLE_PAR_SEGMENTATION_MODELS: ClassVar[List[str]] = ["threshold_bright"] + par_seg_models.AVAILABLE_SEGMENTATION_MODELS
    AVAILABLE_FEATURE_SELECTION: ClassVar[Tuple[str, ...]] = ('random', 'peaks')
    AVAILABLE_SPOT_SPACING_SELECTION: ClassVar[Tuple[str, ...]] = ('random', 'maximized')

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def _normalize_legacy_fields(cls, data: Any) -> Any:
        if isinstance(data, dict):
            cleaned = dict(data)
            legacy_manual = cleaned.pop("is_manual_particle_selection", None)
            if legacy_manual is not None and "par_selection_mode" not in cleaned:
                cleaned["par_selection_mode"] = "manual" if legacy_manual else "auto"
            return cleaned
        return data

    @model_validator(mode="after")
    def _validate(self) -> "PowderMeasurementConfig":
        if self.par_segmentation_model not in self.AVAILABLE_PAR_SEGMENTATION_MODELS:
            raise ValueError(
                f'Value of "par_segmentation_model" set to {self.par_segmentation_model} is invalid. '
                f'Must be one of {self.AVAILABLE_PAR_SEGMENTATION_MODELS}.'
            )
        if self.par_feature_selection not in self.AVAILABLE_FEATURE_SELECTION:
            raise ValueError(
                f'Value of "par_feature_selection" set to {self.par_feature_selection} is invalid. '
                f'Must be one of {self.AVAILABLE_FEATURE_SELECTION}.'
            )
        if self.par_spot_spacing not in self.AVAILABLE_SPOT_SPACING_SELECTION:
            raise ValueError(
                f'Value of "par_spot_spacing" set to {self.par_spot_spacing} is invalid. '
                f'Must be one of {self.AVAILABLE_SPOT_SPACING_SELECTION}.'
            )
        if self.par_spot_selection_mode not in self.AVAILABLE_SPOT_SELECTION_MODES:
            raise ValueError(
                f'Value of "par_spot_selection_mode" set to {self.par_spot_selection_mode} is invalid. '
                f'Must be one of {self.AVAILABLE_SPOT_SELECTION_MODES}.'
            )
        if self.par_selection_mode not in self.AVAILABLE_PAR_SELECTION:
            raise ValueError(
                f'Value of "par_selection_mode" set to {self.par_selection_mode} is invalid. '
                f'Must be one of {self.AVAILABLE_PAR_SELECTION}.'
            )
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
        if self.par_search_frame_width_um is None:
            max_par_radius = np.sqrt(self.max_area_par / np.pi)  # in µm
            self.par_search_frame_width_um = min(20 * max_par_radius, 500.0)  # µm
        return self


class BulkMeasurementConfig(BaseModel):
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
    image_frame_width_um: Optional[float] = None  # µm
    randomize_frames: bool = False
    exclude_sample_margin: bool = False

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="after")
    def _validate(self) -> "BulkMeasurementConfig":
        if not (self.grid_spot_spacing_um > 0):
            raise ValueError("grid_spot_spacing_um must be positive.")
        if not (self.min_xsp_spots_distance_um > 0):
            raise ValueError("min_xsp_spots_distance_um must be positive.")
        if self.min_xsp_spots_distance_um > self.grid_spot_spacing_um:
            raise ValueError("min_xsp_spots_distance_um should not exceed grid_spot_spacing_um.")
        if self.image_frame_width_um is None:
            self.image_frame_width_um = 10 * self.grid_spot_spacing_um
        return self


class ExpStandardsConfig(BaseModel):
    """
    Configuration for the collection of experimental standards.

    Attributes:
        is_exp_std_measurement (bool):
            Whether the configuration corresponds to the measurement of an experimental standard (Default = False)
            If True, a valid `formula` must be provided and weight fractions will be automatically calculated.

        formula (str):
            Chemical formula of the experimental standard. Required if `is_exp_std_measurement` is True.
            Must be parseable by `pymatgen.core.Composition`.

        els_to_use_for_mean_PB_calc (Optional[List[str]]):
            List of element symbols (e.g., ["Fe", "O"]) to use for mean PB calculation.
            Special values:
                ["all"] or ["All"]: use all elements/lines as standards (default)
                None, [], ["none"], or ["None"]: use no elements/lines

        generate_separate_std_dict (bool):
            Whether the acquired reference standard values are added to the current reference dictionary. If True,
            copies the current standard dictionary to the project folder and updates it. If such .json file is already
            present in the project folder, then it updates it. This is generally used when measuring the extent of
            powder precursor intermixing (i.e., powder_meas_cfg_kwargs["is_known_powder_mixture_meas"] = True).

        min_acceptable_PB_ratio (float):
            Minimum PB ratio required for a peak to be accepted as a standard. in cnts/cnts*keV^-1 (Default = 10).

        quant_flags_accepted (List[int]):
            List of quantification flags considered acceptable. Other spectra are filtered out before clustering.
            Quantification flags indicate whether the quantification or the fit of each spectrum is likely to be
            affected by large errors:
                - 0  : Quantification is ok, although it may be affected by large analytical error.
                - -1  : As above, but quantification did not converge within 30 steps.
                - 1  : Error during EDS acquisition. No fit executed.
            2  : Total counts < 90% of target counts, likely due to wrong segmentation. Fit interrupted if interrupt_fits_bad_spectra=True.
                - 3  : Too little low-energy signal, causing poor quantification in that region. Fit interrupted if interrupt_fits_bad_spectra=True.
                - 4  : Poor fit. Fit interrupted if interrupt_fits_bad_spectra=True.
                - 5  : High analytical error (>50%), possibly due to missing element or other major error. Fit interrupted if interrupt_fits_bad_spectra=True.
                - 6  : Excessive X-ray absorption. Fit interrupted if interrupt_fits_bad_spectra=True.
                - 7  : Excessive contamination from substrate.
                - 8  : Too few background counts below reference peak, likely leading to large quantification errors.
                - 9  : Unknown fitting error.
                - 10 : (Only for measurement of experimental standards) Reference peak missing.

        w_frs (Optional[Dict[str, float]]):
            Dictionary of element symbols and their corresponding weight fractions (computed via pymatgen)
            if `is_exp_std_measurement` is True and `formula` is valid; otherwise None.

    Raises:
        ValueError:
            If `is_exp_std_measurement` is True but `formula` is missing or invalid.
    """

    is_exp_std_measurement: bool = False
    formula: str = ''
    els_to_use_for_mean_PB_calc: Optional[List[str]] = ["all"]
    generate_separate_std_dict: bool = False
    min_acceptable_PB_ratio: float = 10
    quant_flags_accepted: List[int] = Field(default_factory=lambda: [0])
    w_frs: Optional[Dict[str, float]] = None

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def _normalize_legacy_els_to_use(cls, data: Any) -> Any:
        # Legacy compatibility: convert bool to list or None
        if isinstance(data, dict):
            val = data.get("use_for_mean_PB_calc", None)
            if isinstance(val, bool):
                if val:
                    data["els_to_use_for_mean_PB_calc"] = ["all"]
                else:
                    data["els_to_use_for_mean_PB_calc"] = None
            # Remove legacy key if present
            data.pop("use_for_mean_PB_calc", None)
        return data

    @model_validator(mode="after")
    def _validate(self) -> "ExpStandardsConfig":
        # Normalize els_to_use_for_mean_PB_calc
        v = self.els_to_use_for_mean_PB_calc
        if v is None or v == [] or v == ["none"] or v == ["None"]:
            self.els_to_use_for_mean_PB_calc = None
        elif v == ["all"] or v == ["All"]:
            self.els_to_use_for_mean_PB_calc = ["all"]
        elif not isinstance(v, list):
            raise ValueError("els_to_use_for_mean_PB_calc must be a list of element symbols, ['all'], or None.")
        # Validate elements if not special values
        if self.els_to_use_for_mean_PB_calc not in (None, ["all"]):
            for el in self.els_to_use_for_mean_PB_calc:
                try:
                    Element(el)
                except Exception:
                    raise ValueError(f"Element symbol '{el}' in els_to_use_for_mean_PB_calc is not a recognized element.")
        if self.is_exp_std_measurement:
            if not self.formula:
                raise ValueError("Formula must be provided when is_exp_std_measurement is True.")
            try:
                comp = Composition(self.formula)
                try:
                    self.w_frs = {el: float(w) for el, w in comp.as_weight_dict.items()} # type: ignore
                except Exception:
                    import warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=FutureWarning)
                        self.w_frs = {el: float(w) for el, w in comp.to_weight_dict.items()}
            except Exception as e:
                raise ValueError(f"Invalid chemical formula '{self.formula}': {e}")
        return self




class PlotConfig(BaseModel):
    """
    Configuration for plotting.

    Attributes:
        show_unused_comps_clust (bool): Whether to plot unused data points in clustering plot.
        els_excluded_clust_plot (List[str]): Elements to exclude in cluster plot when more than 3 elements are present.
        show_legend_clustering bool : Whether to show the legend in the clustering plot. Default: True
        save_plots (bool): Whether to save plots to disk.
        show_plots (bool): Whether to display plots interactively.
        use_custom_plots (bool): Whether to use custom plotting routines.
        custom_plot_file (Optional[str]): Path to a user-editable custom plotting file.
    """

    show_unused_comps_clust: bool = True
    els_excluded_clust_plot: List[str] = Field(default_factory=list)
    show_legend_clustering: bool = True
    save_plots: bool = True
    show_plots: bool = False
    use_custom_plots: bool = False
    custom_plot_file: Optional[str] = None

    model_config = ConfigDict(extra="forbid")


# Dictionary of all config models. Loaded for data import
config_classes_dict = {
    cnst.SAMPLE_CFG_KEY: SampleConfig,
    cnst.MICROSCOPE_CFG_KEY: MicroscopeConfig,
    cnst.MEASUREMENT_CFG_KEY: MeasurementConfig,
    cnst.SAMPLESUBSTRATE_CFG_KEY: SampleSubstrateConfig,
    cnst.QUANTIFICATION_CFG_KEY: QuantificationOptionsConfig,
    cnst.PLOT_CFG_KEY: PlotConfig,
    cnst.POWDER_MEASUREMENT_CFG_KEY: PowderMeasurementConfig,
    cnst.BULK_MEASUREMENT_CFG_KEY: BulkMeasurementConfig,
    cnst.EXP_STD_MEASUREMENT_CFG_KEY: ExpStandardsConfig,
}