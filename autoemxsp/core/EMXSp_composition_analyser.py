#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EMXSp_Composition_Analyzer

Main class for automated compositional analysis of electron microscopy X-ray spectroscopy (EMXSp) data.

Can be run from Run_Acquisition_Quant_Analysis.py

Features:
- Structured configuration for microscope, sample, measurement, and analysis parameters.
- Automated acquisition and quantification of X-ray spectra at electron microscope.
- Filtering and clustering of compositional data.
- Phase identification, mixture analysis, and comprehensive results export.
- Utilities for plotting, saving, and reporting analysis results.

Example Usage
-------------
    # Create analyzer instance
    >>> analyzer = EMXSp_Composition_Analyzer(
            microscope_cfg=microscope_cfg,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            quant_cfg=quant_cfg,
            clustering_cfg=clustering_cfg,
            powder_meas_cfg=powder_meas_cfg,
            plot_cfg=plot_cfg,
            is_acquisition=True,
            development_mode=False,
            output_filename_suffix='',
            verbose=True,
        )

    # Acquire and quantify spectra, and analyse compositions
    >>> analyzer.run_collection_and_quantification(quantify=True)

    # Alternatively, acquire only, then quantify:
    >>> analyzer.run_collection_and_quantification(quantify=False)
    >>> quantify and analyse on another machine using Run_Quantification.py


@author: Andrea
Created on Mon Jul 22 17:43:35 2024
"""

# Standard library imports
import os
import json
import time
import shutil
import itertools
import warnings
import traceback
from datetime import datetime
from dataclasses import asdict
from typing import Any, Optional, Tuple, List, Dict, Iterable, Union

# Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import seaborn as sns
from pymatgen.core.composition import Composition
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score, silhouette_samples
import cvxpy as cp
from yellowbrick.cluster import KElbowVisualizer, SilhouetteVisualizer

# Project-specific imports
from autoemxsp.core.XSp_quantifier import XSp_Quantifier, Quant_Corrections
from autoemxsp.core.EM_controller import EM_Controller, EM_Sample_Finder
import autoemxsp.XSp_calibs as calibs
import autoemxsp.tools.constants as cnst
import autoemxsp.tools.custom_fnctns as custom_fnctns
from autoemxsp.tools.utils import (
    print_single_separator,
    print_double_separator,
    to_latex_formula,
    make_unique_path,
    weight_to_atomic_fr
)
from autoemxsp.tools.config_classes import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    QuantConfig,
    ClusteringConfig,
    PowderMeasurementConfig,
    BulkMeasurementConfig,
    ExpStandardsConfig,
    PlotConfig,
)

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#%% EMXSp_Composition_Analyzer class
class EMXSp_Composition_Analyzer:
    """
    Main class for electron microscopy X-ray spectroscopy (EMXSp) composition analysis.

    This class orchestrates the acquisition, quantification, clustering, and plotting
    of X-ray spectra and composition data, using structured configuration objects for
    all instrument and analysis settings.

    Parameters
    ----------
    microscope_cfg : MicroscopeConfig
        Configuration for the microscope hardware.
    sample_cfg : SampleConfig
        Configuration for the sample.
    measurement_cfg : MeasurementConfig
        Configuration for the measurement/acquisition.
    sample_substrate_cfg : SampleSubstrateConfig
        Configuration for the sample substrate.
    quant_cfg : QuantConfig
        Configuration for spectrum fitting and quantification.
    clustering_cfg : ClusteringConfig
        Configuration for clustering of spectra/compositions.
    powder_meas_cfg : PowderMeasurementConfig
        Configuration for powder measurement.
    bulk_meas_cfg : BulkMeasurementConfig
        Configuration for measurements of bulk or bulk-like samples.
    exp_stds_cfg : ExpStandardsConfig
        Configuration for measurements of experimental standards.
    plot_cfg : PlotConfig
        Configuration for plotting.
    is_acquisition : bool, optional
        If True, indicates class is being used for automated acquisition (default: False).
    standards_dict : dict, optional
        Dictionary of reference PB values from experimental standards. Default : None.
        If None, dictionary of standards is loaded from the XSp_calibs/Your_Microscope_ID directory.
        Provide standards_dict only when providing different standards from those normally used for quantification.
    development_mode : bool, optional
        If True, enables development/debug features (default: False).
    output_filename_suffix : str, optional
        String to append to saved filenames (default: '').
    verbose : bool, optional
        If True, enables verbose output (default: True).
    results_dir : Optional[str], optional
        Directory to save results (default: None). If None, uses default directory, created inside package folder
            - Results, for sample analysis
            - Std_measurements, for experimental standard measurements

    Attributes
    ----------
    All configuration objects are stored as attributes with the same name as their parameter.
    """

    def __init__(
        self,
        microscope_cfg: MicroscopeConfig,
        sample_cfg: SampleConfig,
        measurement_cfg: MeasurementConfig,
        sample_substrate_cfg: SampleSubstrateConfig,
        quant_cfg: QuantConfig = QuantConfig(),
        clustering_cfg: ClusteringConfig = ClusteringConfig(),
        powder_meas_cfg: PowderMeasurementConfig = PowderMeasurementConfig(),
        bulk_meas_cfg: BulkMeasurementConfig = BulkMeasurementConfig(),
        exp_stds_cfg: ExpStandardsConfig = ExpStandardsConfig(),
        plot_cfg: PlotConfig = PlotConfig(),
        is_acquisition: bool = False,
        standards_dict: dict = None,
        development_mode: bool = False,
        output_filename_suffix: str = '',
        verbose: bool = True,
        results_dir: Optional[str] = None,
    ):
        """
        Initialize the EMXSp_Composition_Analyzer with all configuration objects.

        See class docstring for parameter documentation.
        """
        # --- Record process time
        self.start_process_time = time.time()
        if verbose:
            print_double_separator()
            print(f"Starting compositional analysis of sample {sample_cfg.ID}")
            
            
        # --- Define use of class instance
        self.is_acquisition = is_acquisition
        is_XSp_measurement = measurement_cfg.type != measurement_cfg.PARTICLE_STATS_MEAS_TYPE_KEY
        self.development_mode = development_mode
        
        
        # --- System characteristics
        self.microscope_cfg = microscope_cfg
        
        if is_XSp_measurement:
            # Load microscope calibrations for this instrument and mode
            calibs.load_microscope_calibrations(microscope_cfg.ID, measurement_cfg.mode, load_detector_channel_params=is_acquisition)
            if not measurement_cfg.emergence_angle:
                measurement_cfg.emergence_angle = calibs.emergence_angle # Fixed by instrument geometry
        
        
        # --- Measurement configurations
        self.measurement_cfg = measurement_cfg
        self.powder_meas_cfg = powder_meas_cfg
        self.bulk_meas_cfg = bulk_meas_cfg
        self.exp_stds_cfg = exp_stds_cfg
        
        if is_XSp_measurement:
            if is_acquisition:
                # Loaded latest detector calibration values
                meas_modes_calibs = calibs.detector_channel_params
                energy_zero = meas_modes_calibs[measurement_cfg.mode][cnst.OFFSET_KEY]
                bin_width = meas_modes_calibs[measurement_cfg.mode][cnst.SCALE_KEY]
                beam_current = meas_modes_calibs[measurement_cfg.mode][cnst.BEAM_CURRENT_KEY]
                # Store, because needed to call XS_Quantifier
                self.microscope_cfg.energy_zero = energy_zero
                self.microscope_cfg.bin_width = bin_width
                if not measurement_cfg.beam_current:
                    self.measurement_cfg.beam_current = beam_current
                    
                # --- Type checking ---
                for var_name, var_value in [
                    ('energy_zero', energy_zero),
                    ('bin_width', bin_width),
                    ('beam_current', beam_current)
                ]:
                    if not isinstance(var_value, float):
                        raise TypeError(f"{var_name} must be a float, got {type(var_value).__name__}: {var_value}")
            else:
                energy_zero = microscope_cfg.energy_zero
                bin_width = microscope_cfg.bin_width
                beam_current = measurement_cfg.beam_current
                
            self.det_ch_offset = energy_zero
            self.det_ch_width = bin_width

            # Max and min number of EDS spectra to be collected
            self.min_n_spectra = measurement_cfg.min_n_spectra
            if measurement_cfg.max_n_spectra < self.min_n_spectra:
                self.max_n_spectra = self.min_n_spectra
            else:
                self.max_n_spectra = measurement_cfg.max_n_spectra
            
            
        # --- Sample characteristics
        self.sample_cfg = sample_cfg         # Elements possibly present in the sample
        self.sample_substrate_cfg = sample_substrate_cfg
        if is_XSp_measurement:
            # Elements possibly present in the sample
            self.all_els_sample = list(dict.fromkeys(sample_cfg.elements)) #remove any eventual duplicate, keeping original order
            # Detectable elements possibly present in the sample 
            self.detectable_els_sample = [el for el in self.all_els_sample if el not in calibs.undetectable_els]
            # Elements present in the substrate, which have to be subtracted if not present in the sample
            self.all_els_substrate = list(dict.fromkeys(sample_substrate_cfg.elements)) #remove any eventual duplicate, keeping original order
            detectable_els_substrate = [el for el in self.all_els_substrate if el not in calibs.undetectable_els] # remove undetectable elements
            self.detectable_els_substrate = [el for el in detectable_els_substrate if el not in self.detectable_els_sample] #remove any eventual duplicate
            self.is_particle = True if sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE else False
        
        
        # --- Fitting and Quantification
        self.quant_cfg = quant_cfg
        self.standards_dict = standards_dict
        if is_XSp_measurement:
            # Set EDS detector channels to include in the quantification
            self.sp_start, self.sp_end = quant_cfg.spectrum_lims
            # Compute values of energies corresponding to detector channels
            if energy_zero and bin_width:
                self.energy_vals = np.array([energy_zero + bin_width * i for i in range(self.sp_start, self.sp_end)])
            elif is_acquisition and is_XSp_measurement:
                raise ValueError("Missing detector calibration values.\n Please add detector calibration file at {calibs.calibration_files_dir}")
            # Set a threshold value below which counts are considered to be too low
            # Used to filter "bad" spectra out from clustering analysis. All spectra having less counts than this threshold are filtered out
            # Used also to avoid fitting spectra with excessive absorption, which inevitably lead to large quantification errors
            if not quant_cfg.min_bckgrnd_cnts:
                min_bckgrnd_cnts = measurement_cfg.target_acquisition_counts/(2*10**4) # empirical value
                self.quant_cfg.min_bckgrnd_cnts = min_bckgrnd_cnts

        # --- Clustering
        self.clustering_cfg = clustering_cfg
        if is_XSp_measurement:
            self.ref_formulae = list(dict.fromkeys(clustering_cfg.ref_formulae)) # Remove duplicates
            self._calc_reference_phases_df() # Calculate dataframe with reference compositions, and any possible analyitical error deriving from undetectable elements
        
        
        # --- Plotting
        self.plot_cfg = plot_cfg
        
        
        # --- Output
        # Create a new directory if acquiring
        if is_acquisition:
            if results_dir is None:
                if self.exp_stds_cfg.is_exp_std_measurement:
                    results_folder = cnst.STDS_DIR
                else:
                    results_folder = cnst.RESULTS_DIR
                results_dir = make_unique_path(os.path.join(parent_dir, results_folder), sample_cfg.ID)
            else:
                results_dir = make_unique_path(results_dir, sample_cfg.ID)
            os.makedirs(results_dir)

        self.sample_result_dir = results_dir
        self.output_filename_suffix = output_filename_suffix
        self.verbose = verbose

        if is_XSp_measurement:
            # --- Variable initialization
            self.XSp_std_dict = None
            self.sp_coords = [] # List containing particle number + relative coordinates on the image to retrieve exact position where spectra were collected
            self.particle_cntr = -1 # Counter to save particle number
            
            # Initialise lists containing spectral data and comments that will be saved with the quantification data
            self.spectral_data = {key : [] for key in cnst.LIST_SPECTRAL_DATA_KEYS}
            
            # List containing the quantification results of each of the collected spectra (composition and analytical error)
            self.spectra_quant = []
            
        
        # --- Save configurations
        # Save spectrum collection info when class is used to collect spectra
        if is_acquisition:
            self._save_experimental_config(is_XSp_measurement)
        
        
        # --- Initialisations
        # Initialise microscope and XSp analyser
        if is_acquisition:
            if microscope_cfg.type == 'SEM':
                self._initialise_SEM()
            
            if is_XSp_measurement:
                self._initialise_Xsp_analyzer()
        

    #%% Instrument initializations
    # =============================================================================
    def _initialise_SEM(self) -> None:
        """
        Initialize the SEM (Scanning Electron Microscope) and related analysis tools.
    
        Sets up the instrument controller, directories, and, if applicable,
        initializes the particle finder for automated powder sample analysis.
        For circular sample substrates, it automatically detects the C-tape position.
    
        Raises
        ------
        FileNotFoundError
            If the sample result directory does not exist and cannot be created.
        NotImplementedError
            If the sample type is not 'powder'.
        """
        # Determine collection and detection modes based on user-configured settings
        is_manual_navigation = self.measurement_cfg.is_manual_navigation
        is_auto_substrate_detection = self.sample_substrate_cfg.auto_detection
    
        # If using automated collection with a circular sample substrate, detect C-tape and update sample coordinates (center, radius)
        if not is_manual_navigation and is_auto_substrate_detection:
            sample_finder = EM_Sample_Finder(
                microscope_ID=self.microscope_cfg.ID,
                center_pos=self.sample_cfg.center_pos,
                sample_half_width_mm=self.sample_cfg.half_width_mm,
                substrate_width_mm=self.sample_substrate_cfg.stub_w_mm,
                results_dir=self.sample_result_dir,
                verbose=self.verbose
            )
            if self.sample_substrate_cfg.type == cnst.CTAPE_SUBSTRATE_TYPE:
                Ctape_coords = sample_finder.detect_Ctape()
                if Ctape_coords:
                    center_pos, C_tape_r = Ctape_coords
                    # Update detected center position and half-width
                    self.sample_cfg.center_pos = center_pos
                    self.sample_cfg.half_width_mm = C_tape_r
            else:
                warnings.warn(f"Automatic detection is only implemented for {cnst.ALLOWED_AUTO_DETECTION_TYPES}")
        
        # Set up image directory for this sample
        EM_images_dir = os.path.join(self.sample_result_dir, cnst.IMAGES_DIR)
        if not os.path.exists(EM_images_dir):
            try:
                os.makedirs(EM_images_dir)
            except Exception as e:
                raise FileNotFoundError(f"Could not create results directory: {EM_images_dir}") from e
    
        # Initialise instrument controller
        self.EM_controller = EM_Controller(
            self.microscope_cfg,
            self.sample_cfg,
            self.measurement_cfg,
            self.sample_substrate_cfg,
            self.powder_meas_cfg,
            self.bulk_meas_cfg,
            results_dir=EM_images_dir,
            verbose=self.verbose
        )
        self.EM_controller.initialise_SEM()
        self.EM_controller.initialise_sample_navigator(self.EM_controller, exclude_sample_margin=True)
        
        # Update employed working distance
        self.measurement_cfg.working_distance = self.EM_controller.measurement_cfg.working_distance
            
            
    def _initialise_Xsp_analyzer(self):
        """
        Initialize the X-ray spectroscopy analyzer according to the measurement configuration.
    
        If the measurement type is 'EDS', this initializes the EDS (Energy Dispersive X-ray Spectroscopy)
        analyzer via the associated EM_controller. For any other measurement type, a NotImplementedError
        is raised.
    
        Raises
        ------
        NotImplementedError
            If the measurement type is not 'EDS'.
        """
        # Only EDS is supported at present
        if self.measurement_cfg.type == 'EDS':
            self.EM_controller.initialise_XS_analyzer()
        elif self.measurement_cfg.type not in self.measurement_cfg.ALLOWED_TYPES:
            raise NotImplementedError(
                f"X-ray spectroscopy analyzer initialization for measurement type '{self.measurement_cfg.type}' is not currently implemented."
            )
        
    #%% Other initializations
    # =============================================================================
    def _make_analysis_dir(self) -> None:
        """
        Create a unique directory for saving analysis results using `make_unique_dir`.
    
        The directory name is based on `cnst.ANALYSIS_DIR` and `self.output_filename_suffix`.
        If a directory with the target name exists, a counter is appended to ensure uniqueness.
    
        The resulting directory path is stored in `self.analysis_dir`.
    
        Raises
        ------
        FileExistsError
            If unable to create a unique analysis directory due to file system errors.
        """
        # Compose base directory name for analysis results
        base_name = cnst.ANALYSIS_DIR + self.output_filename_suffix
        try:
            analysis_dir = make_unique_path(self.sample_result_dir, base_name)
            os.makedirs(analysis_dir)
        except Exception as e:
            raise FileExistsError(f"Could not create analysis directory in '{self.sample_result_dir}' with base name '{base_name}'.") from e
    
        self.analysis_dir = analysis_dir
    
    
    def _initialise_std_dict(self) -> None:
        """
        Initialise the dictionary of X-ray standards for quantification.
    
        This method determines how the `XSp_std_dict` attribute is initialised
        based on the sample configuration and measurement type:
    
        - If the sample is a powder sample (`S_POWDER_SAMPLE_TYPE`) and the
          measurement is of a known powder mixture, the standards dictionary
          is compiled from reference data using `_compile_standards_from_references()`.
    
        - Otherwise, the standards dictionary is expected to be loaded directly
          within the `XSp_Quantifier` and is set to `None` here.
    
        Returns
        -------
        None
            This method modifies the `self.XSp_std_dict` attribute in place.
        """
        is_powder_sample = self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE
        is_known_mixture = getattr(self.powder_meas_cfg, "is_known_powder_mixture_meas", False)
    
        if is_powder_sample and is_known_mixture:
            self.XSp_std_dict = self._compile_standards_from_references()
        else:
            # Standards dictionary will be loaded directly within the `XSp_Quantifier`
            self.XSp_std_dict = None


    def _calc_reference_phases_df(self) -> None:
        """
        Calculate the compositions of reference phases and store them in a pd.DataFrame.
    
        For each reference formula in `self.ref_formulae`, this method:
          - Computes the composition using pymatgen's Composition class.
          - Computes either mass or atomic fractions, depending on clustering configuration.
          - Accounts for undetectable elements and calculates the maximum analytical error due to their presence.
          - Stores the resulting phase compositions in `self.ref_phases_df` and the weights in `self.ref_weights_in_mixture`.
    
        If no reference formulae are provided, the function exits without error.
    
        Warnings
        --------
        Issues a warning if a formula cannot be parsed or if no detectable elements are found in a formula.
    
        Raises
        ------
        ValueError
            If an unknown clustering feature set is specified.
        """
        import warnings
    
        undetectable_an_err = 0
        ref_phases = []
        ref_weights_in_mixture = []
        
        # Check if self.ref_formulae is set to None
        if not self.ref_formulae:
            # No reference formulae provided; nothing to do
            self.ref_phases_df = pd.DataFrame(columns=self.all_els_sample)
            self.ref_weights_in_mixture = []
            self.undetectable_an_er = 0
            return
        
        for formula in self.ref_formulae:
            # Use pymatgen class Composition
            try:
                comp = Composition(formula)
            except Exception as e:
                warnings.warn(f"Invalid chemical formula '{formula}': {e}")
                continue
    
            # Get mass fractions as dictionary el: w_fr
            w_fr_dict = comp.to_weight_dict
    
            # Check for detectable elements at the beginning
            detectable_in_formula = [el for el in self.detectable_els_sample if el in w_fr_dict]
            if not detectable_in_formula:
                warnings.warn(f"No detectable elements found in formula '{formula}'.")
                continue
    
            # Calculate analytical error due to undetectable elements
            for el, w_fr in w_fr_dict.items():
                if el in calibs.undetectable_els:
                    undetectable_an_err = max(undetectable_an_err, w_fr)
    
            if self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
                # Mass fractions are not normalised, so a negative analytical error is possible when undetectable elements are present
                # Calculate reference dictionary considering only quantified elements (e.g. Li is ignored)
                phase = {el: w_fr_dict.get(el, 0) for el in self.detectable_els_sample}
                # Store weight of reference in an eventual mixture, which is simply equal to the compound molar weight
                ref_weights_in_mixture.append(comp.weight)
    
            elif self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
                # Atomic fractions are normalised, so for the purpose of reference phases we should calculate it normalising the
                # mass fractions, after discarding the undetectable elements
                detectable_w_frs = {el: w_fr for el, w_fr in w_fr_dict.items() if el in self.detectable_els_sample}
                # Transform to Composition class
                detectable_comp = comp.from_weight_dict(detectable_w_frs)
                # Get dictionary of el : at_fr
                phase = detectable_comp.fractional_composition.as_dict()
                # Store weight of reference in an eventual mixture, which is equal to the number of atoms in reference formula, without undetectable elements
                ref_weight = sum(at_n for el, at_n in comp.get_el_amt_dict().items() if el in self.detectable_els_sample)
                ref_weights_in_mixture.append(ref_weight)
            else:
                raise ValueError(f"Unknown clustering feature set: {self.clustering_cfg.features}")
    
            ref_phases.append(phase)
    
        # Convert to pd.DataFrame and store it
        ref_phases_df = pd.DataFrame(ref_phases, columns=self.all_els_sample).fillna(0)
        self.ref_phases_df = ref_phases_df
    
        # Store values of reference weights used to calculate molar fractions from mixtures
        self.ref_weights_in_mixture = ref_weights_in_mixture
    
        # Calculate negative analytical error accepted to compensate for elements undetectable by EDS (H, He, Li, Be)
        self.undetectable_an_er = undetectable_an_err
        
        
    #%% Single spectrum operations
    # =============================================================================            
    def _acquire_spectrum(self, x: float, y: float) -> Tuple:
        """
        Acquire an X-ray spectrum at the specified stage position and store the results.
    
        Parameters
        ----------
        x, y : float
            X, Y coordinatesfor the spectrum acquisition.
            Coordinate System
            ----------------
            The coordinates are expressed in a normalized, aspect-ratio-correct system centered at the image center:
    
                - The origin (0, 0) is at the image center.
                - The x-axis is horizontal, increasing to the right, ranging from -0.5 (left) to +0.5 (right).
                - The y-axis is vertical, increasing downward, and scaled by the aspect ratio (height/width):
                    * Top edge:    y = -0.5 × (height / width)
                    * Bottom edge: y = +0.5 × (height / width)
                
                |        (-0.5, -0.5*height/width)         (0.5, -0.5*height/width)
                |                       +-------------------------+
                |                       |                         |
                |                       |                         |
                |                       |           +(0,0)        |-----> +x
                |                       |                         |
                |                       |                         |
                v  +y                   +-------------------------+
                        (-0.5,  0.5*height/width)         (0.5, 0.5*height/width)
    
            This ensures the coordinate system is always centered and aspect-ratio-correct, regardless of image size.
    
        Returns
        -------
        spectrum_data : np.array
            The acquired spectrum data.
        background_data : np.array, None
            The acquired background data.
        real_time : float
            Real acquisition time used.
        live_time : float
            Live acquisition time used.
    
        Notes
        -----
        - Results are appended to self.spectral_data using the keys defined in `cnst`.
        """
        # Get spectral data from the EM controller
        spectrum_data, background_data, real_time, live_time = self.EM_controller.acquire_XS_spot_spectrum(
            x, y,
            self.measurement_cfg.max_acquisition_time,
            self.measurement_cfg.target_acquisition_counts
        )
    
        # Store results in the spectral_data dictionary
        self.spectral_data[cnst.SPECTRUM_DF_KEY].append(spectrum_data)
        self.spectral_data[cnst.BACKGROUND_DF_KEY].append(background_data)
        self.spectral_data[cnst.REAL_TIME_DF_KEY].append(real_time)
        self.spectral_data[cnst.LIVE_TIME_DF_KEY].append(live_time)
    
        return spectrum_data, background_data, real_time, live_time
    
    
    def _fit_exp_std_spectrum(
        self,
        spectrum: Iterable,
        background: Optional[Iterable] = None,
        sp_collection_time: float = None,
        els_w_frs: Optional[Dict[str,float]] = None,
        sp_id: str = '',
        verbose: bool = True
    ) -> Optional[Dict]:
        """
        Quantify a single X-ray spectrum.
    
        This method checks if the spectrum is valid for fitting, runs the quantification,
        flags the result as necessary, and appends comments and quantification flags to
        the spectral data attributes.
    
        Parameters
        ----------
        spectrum : Iterable
            The spectrum data to be quantified.
        background : Iterable, optional
            The background data associated with the spectrum.
        sp_collection_time : float, optional
            The collection time for the spectrum.
        sp_id: str, optional
            The spectrum ID, used as label for printing
        verbose : bool, optional
            If True, enables verbose output (default: True).
    
        Returns
        -------
        fit_result : Dict or None
            Dictionary returned by XSp_Quantifier, containing calculated composition in atomic fractions and
            analytical error, or None if the spectrum is not suitable for quantification or fitting fails.
    
        Notes
        -----
        - Filtering flags are appended through function _check_fit_quant_validity().
        """
        if verbose:
            if sp_id != '':
                sp_id_str = " #" + sp_id
            else:
                sp_id_str = '...'
            print_single_separator()
            print('Fitting spectrum' + sp_id_str)
            start_quant_time = time.time()
                            
        # Check if spectrum is worth fitting
        if not self._is_spectrum_valid_for_fitting(spectrum, background):
            return None
        
        # Initialize class to quantify spectrum
        quantifier = XSp_Quantifier(
            spectrum_vals=spectrum,
            spectrum_lims=(self.sp_start, self.sp_end),
            microscope_ID=self.microscope_cfg.ID,
            meas_type=self.measurement_cfg.type,
            meas_mode=self.measurement_cfg.mode,
            det_ch_offset=self.det_ch_offset,
            det_ch_width=self.det_ch_width,
            beam_e=self.measurement_cfg.beam_energy_keV,
            emergence_angle=self.measurement_cfg.emergence_angle,
            energy_vals=None,
            background_vals=background,
            els_sample=self.all_els_sample,
            els_substrate=self.detectable_els_substrate,
            els_w_fr=self.exp_stds_cfg.w_frs,
            is_particle=self.is_particle,
            sp_collection_time=sp_collection_time,
            max_undetectable_w_fr=self.undetectable_an_er,
            fit_tol=self.quant_cfg.fit_tolerance,
            standards_dict=self.XSp_std_dict,
            verbose=False,
            fitting_verbose=False
        )
        
        try:
            bad_quant_flag = quantifier.initialize_and_fit_spectrum(print_results=self.verbose)
            is_fit_valid = True
            min_bckgrnd_ref_lines = quantifier._get_min_bckgrnd_cnts_ref_quant_lines()
        except Exception as e:
            is_fit_valid = False
            print(f"{type(e).__name__}: {e}")
            traceback.print_exc()  # prints the full traceback
            self._check_fit_quant_validity(is_fit_valid, None, None, None)
            return None
        
        fit_results_dict, are_all_ref_peaks_present = self._assemble_fit_info(quantifier)
        
        if are_all_ref_peaks_present:
            self._check_fit_quant_validity(is_fit_valid, bad_quant_flag, quantifier, min_bckgrnd_ref_lines)
        else:
            self.spectral_data[cnst.COMMENTS_DF_KEY].append("Reference peak missing")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(10)
        
        if verbose:
            fit_time = time.time() - start_quant_time
            print(f"Fitting took {fit_time:.2f} s")
    
        return fit_results_dict
    
    
    def _assemble_fit_info(self, quantifier):
        are_all_ref_peaks_present = True
        
        # Get fit result data to retrieve PB ratio 
        fit_data = quantifier.fitted_peaks_info
        
        reduced_chi_squared = quantifier.fit_result.redchi
        r_squared = 1 - quantifier.fit_result.residual.var() / np.var(quantifier.spectrum_vals)
        
        # Initialise variables
        PB_ratios_d = {} # Dictionary used to store the PB ratios of each line fitted in the spectrum
        
        # Store PB ratios from fitted peaks
        el_lines = [el_line for el_line in fit_data.keys() if 'esc' not in el_line and 'pileup' not in el_line]
        for el_line in el_lines:
            el, line = el_line.split('_')[:2]
            
            if el not in self.detectable_els_sample:
                continue # Do not store PB ratios for substrate elements
            
            meas_PB_ratio = fit_data[el_line][cnst.PB_RATIO_KEY]

            # Assign a nan value if PB ratio is too low, to later filter only the significant peaks
            if meas_PB_ratio < self.exp_stds_cfg.min_acceptable_PB_ratio:
                meas_PB_ratio = np.nan
            
            # Store PB ratio information
            if line in quantifier.xray_quant_ref_lines:
                # Store PB-ratio value, only for reference peaks
                PB_ratios_d[el_line] = meas_PB_ratio
                
                # Store theoretical energy values for fitted peaks
                self._th_peak_energies[el_line] = fit_data[el_line][cnst.PEAK_TH_ENERGY_KEY] 
                
                other_xray_ref_lines = [l for l in quantifier.xray_quant_ref_lines if l != line]

                # Elements of the standard must be properly fitted, and possess a background with enough counts
                if el in self.detectable_els_sample and all(el + '_' + l not in fit_data.keys() for l in other_xray_ref_lines):
                    # Check if peak is present
                    if not meas_PB_ratio > 0:
                        are_all_ref_peaks_present = False
                        # Reference peak not present
                        if self.verbose:
                            print(f"{el_line} reference peak missing.")
                    
        # Create dictionary of fit results
        fit_results_dict = {**PB_ratios_d, cnst.R_SQ_KEY : r_squared, cnst.REDCHI_SQ_KEY : reduced_chi_squared}

        # Append to list of results
        return fit_results_dict, are_all_ref_peaks_present
        
    
    def _fit_quantify_spectrum(
        self,
        spectrum: Iterable,
        background: Optional[Iterable] = None,
        sp_collection_time: float = None,
        sp_id: str = '',
        verbose: bool = True
    ) -> Optional[Dict]:
        """
        Quantify a single X-ray spectrum.
    
        This method checks if the spectrum is valid for fitting, runs the quantification,
        flags the result as necessary, and appends comments and quantification flags to
        the spectral data attributes.
    
        Parameters
        ----------
        spectrum : Iterable
            The spectrum data to be quantified.
        background : Iterable, optional
            The background data associated with the spectrum.
        sp_collection_time : float, optional
            The collection time for the spectrum.
        sp_id: str, optional
            The spectrum ID, used as label for printing
        verbose : bool, optional
            If True, enables verbose output (default: True).
    
        Returns
        -------
        quant_result : Dict or None
            Dictionary returned by XSp_Quantifier, containing calculated composition in atomic fractions and
            analytical error, or None if the spectrum is not suitable for quantification or fitting fails.
    
        Notes
        -----
        - Filtering flags are appended through function _check_fit_quant_validity().
        """
        if verbose:
            if sp_id != '':
                sp_id_str = " #" + sp_id
            else:
                sp_id_str = '...'
            print_single_separator()
            print('Quantifying spectrum' + sp_id_str)
            start_quant_time = time.time()
                            
        # Check if spectrum is worth fitting
        if not self._is_spectrum_valid_for_fitting(spectrum, background):
            return None
        
        # Initialize class to quantify spectrum
        quantifier = XSp_Quantifier(
            spectrum_vals=spectrum,
            spectrum_lims=(self.sp_start, self.sp_end),
            microscope_ID=self.microscope_cfg.ID,
            meas_type=self.measurement_cfg.type,
            meas_mode=self.measurement_cfg.mode,
            det_ch_offset=self.det_ch_offset,
            det_ch_width=self.det_ch_width,
            beam_e=self.measurement_cfg.beam_energy_keV,
            emergence_angle=self.measurement_cfg.emergence_angle,
            energy_vals=None,
            background_vals=background,
            els_sample=self.all_els_sample,
            els_substrate=self.detectable_els_substrate,
            els_w_fr=self.sample_cfg.w_frs,
            is_particle=self.is_particle,
            sp_collection_time=sp_collection_time,
            max_undetectable_w_fr=self.undetectable_an_er,
            fit_tol=self.quant_cfg.fit_tolerance,
            standards_dict=self.XSp_std_dict,
            verbose=False,
            fitting_verbose=False
        )
        
        try:
            # Returns dictionary containing calculated composition in atomic fractions + analytical error
            quant_result, min_bckgrnd_ref_lines, bad_quant_flag = quantifier.quantify_spectrum(
                print_result=False,
                interrupt_fits_bad_spectra=self.quant_cfg.interrupt_fits_bad_spectra
            )
            is_quant_fit_valid = True
        except Exception as e:
            print(f"{type(e).__name__}: {e}")
            is_quant_fit_valid = False
            self._check_fit_quant_validity(is_quant_fit_valid, None, None, None)
            return None
        else:
            self._check_fit_quant_validity(is_quant_fit_valid, bad_quant_flag, quantifier, min_bckgrnd_ref_lines)
        
        if verbose and quant_result:
            quantification_time = time.time() - start_quant_time
            for el in quant_result[cnst.COMP_AT_FR_KEY].keys():
                print(f"{el} at%: {quant_result[cnst.COMP_AT_FR_KEY][el]*100:.2f}%")
            print(f"An. er.: {quant_result[cnst.AN_ER_KEY]*100:.2f}%")
            print(f"Quantification took {quantification_time:.2f} s")
    
        return quant_result


    def _check_fit_quant_validity(
        self,
        is_quant_fit_valid: bool,
        bad_quant_flag: int,
        quantifier: Any,
        min_bckgrnd_ref_lines: Any
    ) -> None:
        """
        Append comments and flags to spectral data based on fit and quantification outcome.
    
        Parameters
        ----------
        is_quant_fit_valid : bool,
            Whether an error occurred during fit or quantification
        bad_quant_flag : int
            Flag returned by fit and quantification indicating the type of issue detected.
            (1: poor fit, 2: high analytical error, 3: excessive absorption, -1: non-converged fit)
        quantifier : object
            The quantifier instance used for this spectrum (may be used for further checks).
        min_bckgrnd_ref_lines : Any
            Reference value for background lines, used for further spectrum checks.
    
        Returns
        -------
        None
        """
        # Prefix for comments if fit was interrupted
        start_str_comments = 'Fit interrupted due to ' if not is_quant_fit_valid else ''
    
        if bad_quant_flag == 1:
            if self.verbose and is_quant_fit_valid:
                print("Flagged for poor fit")
            self.spectral_data[cnst.COMMENTS_DF_KEY].append(start_str_comments + "poor fit")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(4)
        elif bad_quant_flag == 2:
            if self.verbose and is_quant_fit_valid:
                print("Flagged for excessively high analytical error")
            self.spectral_data[cnst.COMMENTS_DF_KEY].append(start_str_comments + "excessively high analytical error")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(5)
        elif bad_quant_flag == 3:
            if self.verbose and is_quant_fit_valid:
                print("Flagged for excessive X-ray absorption")
            self.spectral_data[cnst.COMMENTS_DF_KEY].append(start_str_comments + "excessive X-ray absorption")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(6)
        elif not is_quant_fit_valid:
            self.spectral_data[cnst.COMMENTS_DF_KEY].append("Fit interrupted for unknown reasons")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(9)
        else:
            # Fit completed with no apparent issue; check for low background counts, etc.
            self._flag_spectrum_for_clustering(min_bckgrnd_ref_lines, quantifier)
    
        # If fit did not converge, append convergence comment and flag
        if bad_quant_flag == -1:
            self.spectral_data[cnst.COMMENTS_DF_KEY][-1] += " - Quantification did not converge."
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY][-1] = -1  # Signal non-convergence
        
        
    def _is_spectrum_valid_for_fitting(
        self, 
        spectrum: np.ndarray, 
        background: np.ndarray = None
    ) -> bool:
        """
        Check if a spectrum is valid for quantification fitting.
    
        This method applies several criteria to determine if a spectrum should be processed:
          - No spectrum data present.
          - Total counts are too low.
          - Too many low-count channels in the low-energy range.
    
        For each failure, a comment and quantification flag are appended to `self.spectral_data`, and
        a message is printed if `self.verbose` is True.
    
        Parameters
        ----------
        spectrum : np.ndarray
            The spectrum data to be validated.
        background : np.ndarray, optional
            The background data (not used in this method).
    
        Returns
        -------
        is_spectrum_valid : bool
            True if the spectrum is valid for fitting, False otherwise.
    
        Notes
        -----
        - Assumes all class attributes and keys are correctly initialized.
        - Uses constants from `cnst` for comment and flag keys.
        """
        is_spectrum_valid = True
    
        if spectrum is None:
            # Check if spectrum data is present
            is_spectrum_valid = False
            self.spectral_data[cnst.COMMENTS_DF_KEY].append("No spectral data present")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(1)
            if self.verbose:
                print("Error during spectrum collection. No quantification was done.")
        elif np.sum(spectrum) < 0.9 * self.measurement_cfg.target_acquisition_counts:
            # Skip quantification of spectrum when counts are too low
            is_spectrum_valid = False
            self.spectral_data[cnst.COMMENTS_DF_KEY].append("Total counts too low")
            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(2)
            if self.verbose:
                print(f"Quantification skipped due to spectrum counts lower than 90% of the target counts of {self.measurement_cfg.target_acquisition_counts}")
        else:
            # Skip quantification if too many low values, which leads to errors due to imprecise fitting
            n_vals_considered = 20  # Number of data channels that must be low for spectrum to be excluded
            filter_len = 3
            en_threshold = 2  # keV
    
            # Prepare (energy, counts) pairs for the relevant region
            xy_data = zip(self.energy_vals, spectrum[self.sp_start: self.sp_end])
            # Consider only data with counts > 0 and energy < threshold
            spectrum_data_to_consider = [cnts for en, cnts in xy_data if cnts > 0 and en < en_threshold]
            # Smoothen spectrum to reduce noise
            spectrum_smooth = np.convolve(spectrum_data_to_consider, np.ones(filter_len)/filter_len, mode='same')
            # Get the n lowest values in the smoothed spectrum
            min_vals = np.sort(spectrum_smooth)[:n_vals_considered]
            if all(min_vals < self.quant_cfg.min_bckgrnd_cnts):
                is_spectrum_valid = False
                self.spectral_data[cnst.COMMENTS_DF_KEY].append("Background counts too low")
                self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(3)
                if self.verbose:
                    print(f"Quantification skipped due to at least {n_vals_considered} spectrum points with E < {en_threshold} keV having a count lower than {self.quant_cfg.min_bckgrnd_cnts}")
                    print("This generally indicates an excessive absorption of X-rays before they reach the detector, which compromises accurate measurements of PB ratios.")
    
        return is_spectrum_valid
    
    
    def _flag_spectrum_for_clustering(
        self,
        min_bckgrnd_ref_lines: float,
        quantifier: Any,
    ) -> bool:
        """
        Check spectrum validity for clustering based on substrate peak intensities and background counts.
    
        This method:
          - Flags spectra where any substrate element has a peak intensity larger than a set percentage
            of total counts.
          - Flags spectra where the minimum background counts under reference peaks are too low.
          - Appends comments and quantification flags to `self.spectral_data` using keys from `cnst`.
          - Prints warnings if `self.verbose` is True.
    
        Parameters
        ----------
        min_bckgrnd_ref_lines : float
            Minimum average counts under reference peaks in the spectrum.
        quantifier : Any
            The quantification object containing fitting information.
    
        Returns
        -------
        is_spectrum_valid : bool
            True if the spectrum passes all checks, False otherwise.
    
        Notes
        -----
        - Assumes all class attributes and keys are correctly initialized.
        """
        is_spectrum_valid = True
    
        # Check that substrate signal is not too high
        sub_peak_int_threshold = 10  # % of total counts
        sub_peak_int_thresh_cnts = quantifier.tot_sp_counts * sub_peak_int_threshold / 100
    
        # Sum intensities from substrate peaks
        els_substrate_intensities = {el: 0 for el in self.detectable_els_substrate}  # initialise dictionary of peak intensities
        for el_line, peak_info in quantifier.fitted_peaks_info.items():
            el = el_line.split('_')[0]
            if el in self.detectable_els_substrate:
                els_substrate_intensities[el] += peak_info[cnst.PEAK_INTENSITY_KEY]
    
        # Check that no substrate element has too high intensity
        for el, peak_int in els_substrate_intensities.items():
            if peak_int > sub_peak_int_thresh_cnts:
                is_spectrum_valid = False
                self.spectral_data[cnst.COMMENTS_DF_KEY].append(
                    f"{el} {peak_int:.0f} counts > {sub_peak_int_threshold} % of total counts"
                )
                self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(7)
                if self.verbose:
                    print(f"Intensity of substrate element {el} is {peak_int:.0f} cnts, larger than {sub_peak_int_threshold}% of total counts")
                    print("This is likely to lead to large quantification errors.")
                break  # Stop if one element has too high intensity
    
        # Check that background intensity is high enough
        if is_spectrum_valid:
            self.spectral_data[cnst.COMMENTS_DF_KEY].append(f"{min_bckgrnd_ref_lines:.1f} min. ref. bckgrnd counts")
            # Spectrum is not valid if any of the reference peaks has average counts lower than self.quant_cfg.min_bckgrnd_cnts
            if min_bckgrnd_ref_lines < self.quant_cfg.min_bckgrnd_cnts:
                is_spectrum_valid = False
                self.spectral_data[cnst.COMMENTS_DF_KEY][-1] += ', too low'
                self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(8)
                if self.verbose:
                    print(f"Counts below a reference peak are on average < {self.quant_cfg.min_bckgrnd_cnts}")
                    print("This is likely to lead to large quantification errors.")
            else:
                self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(0)  # Quantification is ok
    
        return is_spectrum_valid  # Not used, but returned for completeness
    
    
    #%% Spectra acquisition and quantification routines
    # ============================================================================= 
    def _collect_spectra(
        self,
        n_spectra_to_collect: int,
        n_tot_sp_collected: int = 0,
        quantify: bool = True
    ) -> Tuple[int, bool]:
        """
        Acquire and optionally quantify spectra from particles.
    
        This method supports two operational modes:
          - Collection and quantification (default): For each spot, acquire and immediately quantify the spectrum.
          - Collection only: Only acquire spectra and update coordinates; quantification is deferred.
                              Useful when quantifying spectra separately
    
        Parameters
        ----------
        n_spectra_to_collect : int
            Number of new spectra to collect.
        n_tot_sp_collected : int, optional
            The running total of spectra already collected (default: 0).
        quantify : bool, optional
            If True, perform spectra quantification (default: True).
    
        Returns
        -------
        n_tot_sp_collected : int
            The updated total number of spectra collected after this session.
        success : bool
            False if collection was interrupted by user, or if no more particles could be found. True otherwise.
    
        Notes
        -----
        - If `quantify` is True, quantification occurs immediately after each collection.
        """
        success = False
    
        n_spectra_collected = 0
        n_spectra_init = n_tot_sp_collected

        while n_spectra_collected < n_spectra_to_collect:
            success, spots_xy_list, particle_cntr = self.EM_controller.get_XSp_coords(n_tot_sp_collected)
            
            if not success:
                break
            
            self.particle_cntr = particle_cntr
            frame_ID = self.EM_controller.frame_labels[self.EM_controller._frame_cntr-1]
            
            latest_spot_id = None # For image annotations
            for i, (x, y) in enumerate(spots_xy_list):
                latest_spot_id = i
                value_map = {
                    cnst.SP_ID_DF_KEY: n_tot_sp_collected,
                    cnst.FRAME_ID_DF_KEY : frame_ID,
                    cnst.SP_X_COORD_DF_KEY: f'{x:.3f}',
                    cnst.SP_Y_COORD_DF_KEY: f'{y:.3f}'
                }
                # Add particle ID only if not None
                if self.particle_cntr is not None:
                    value_map[cnst.PAR_ID_DF_KEY] = self.particle_cntr
                    
                self.sp_coords.append({
                    key: value_map[key]
                    for key in cnst.LIST_SPECTRUM_COORDINATES_KEYS
                    if key in value_map
                }) # Ensures any modification of keys is done at the level of LIST_SPECTRUM_COORDINATES_KEYS
                    # This allows correct loading when quantifying or analysing spectra after acquisition

                if self.verbose:
                    print_single_separator()
                    print(f'Acquiring spectrum #{n_tot_sp_collected}...')

                n_tot_sp_collected += 1
                spectrum_data, background_data, collection_time, live_time = self._acquire_spectrum(x, y)
                
                if self.verbose:
                    print(f"Acquisition took {collection_time:.2f} s")
                
                # Contamination check: skip quantification if counts are too low (only at first measurement spot)
                if i==0 and self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE:
                    if sum(spectrum_data) < 0.95 * self.measurement_cfg.target_acquisition_counts:
                        if quantify:
                            self.spectra_quant.append(None)
                            self.spectral_data[cnst.COMMENTS_DF_KEY].append("Total counts too low")
                            self.spectral_data[cnst.QUANT_FLAG_DF_KEY].append(2)
                        if self.verbose:
                                print('Current particle is unlikely to be part of the sample.\nSkipping to the next particle.')
                                print('Increase measurement_cfg.max_acquisition_time if this behavior is undesired.')
                        break
            
            # Save image of particle, with ID of acquired XSp spots
            if latest_spot_id is not None:
                # Prepare save path
                par_cntr_str = f"_par{self.particle_cntr}" if self.particle_cntr is not None else ''
                filename = f"{self.sample_cfg.ID}{par_cntr_str}_fr{frame_ID}_xyspots"
                # Construct annotation dictionary
                im_annotations = []
                for i, xy_coords in enumerate(spots_xy_list):
                    # Skip if latest_spot_id is None or i is out of range
                    if latest_spot_id is None or i > latest_spot_id:
                        break
                
                    xy_center = self.EM_controller.convert_XS_coords_to_pixels(xy_coords)
                    if xy_center is None:
                        continue
                    
                    im_annotations.append({
                        self.EM_controller.an_text_key: (
                            str(n_tot_sp_collected - 1 - latest_spot_id + i),
                            (xy_center[0] - 30, xy_center[1] - 15)
                        ),
                        self.EM_controller.an_circle_key: (10, xy_center, -1)
                    })
                # Save image with annotations
                self.EM_controller.save_frame_image(filename, im_annotations = im_annotations)
                
            if quantify:
                self._fit_and_quantify_spectra()

            n_spectra_collected = n_tot_sp_collected - n_spectra_init
    
        return n_tot_sp_collected, success
    
    
    def _fit_and_quantify_spectra(self, quantify: bool = True) -> None:
        """
        Fit and )optionally) quantify all collected spectra that have not yet been processed.
        
        Used when spectral acquisition and quantification are performed separately,
         or when measuring standards, for which quantify msut be set to False.
         
        Parameters
        ----------
        quantify: bool
            If False, does not perform quantification. Used for fitting of standards (Default = True).
        
        This method:
          - Iterates over all unquantified spectra in self.spectral_data.
          - Retrieves the spectrum, (optionally) background, and collection time for each.
          - If quantify == True, performs quantification using self._fit_quantify_spectrum,
               otherwise only fits the spectra.
          - Prints results and timing information if self.verbose is True.
          - Appends each quantification result to self.spectra_quant.
        """
        # Get position of first spectrum to be quantified
        quant_sp_cntr = len(self.spectra_quant)
    
        # Quantify all spectra that have not been quantified yet
        tot_spectra_collected = len(self.spectral_data[cnst.SPECTRUM_DF_KEY])
        n_spectra_to_quant = tot_spectra_collected - quant_sp_cntr
        if self.verbose and n_spectra_to_quant>0:
            print_single_separator()
            quant_str = "quantification" if quantify else "fitting"
            print(f"Starting {quant_str} of {n_spectra_to_quant} spectra")
        for i in range(quant_sp_cntr, tot_spectra_collected):
    
            # Get spectral (and background) data
            spectrum = self.spectral_data[cnst.SPECTRUM_DF_KEY][i]
            if self.quant_cfg.use_instrument_background:
                background = self.spectral_data[cnst.BACKGROUND_DF_KEY][i]
            else:
                background = None
    
            sp_collection_time = self.spectral_data[cnst.LIVE_TIME_DF_KEY][i]
    
            # Quantify spectrum. Returns dictionary containing calculated composition in atomic fractions + analytical error
            if quantify:
                result = self._fit_quantify_spectrum(
                                        spectrum,
                                        background,
                                        sp_collection_time,
                                        sp_id = f'{i}/{tot_spectra_collected - 1}'
                            )
            else:
                result = self._fit_exp_std_spectrum(
                                        spectrum,
                                        background,
                                        sp_collection_time,
                                        sp_id = f'{i}/{tot_spectra_collected - 1}'
                            )
    
            # Store result
            self.spectra_quant.append(result)
    
    #%% Find number of clusters in kmeans
    # ============================================================================= 
    def _find_optimal_k(self, compositions_df, k, compute_k_only_once = False):
        """
        Determine the optimal number of clusters for k-means.
    
        Returns
        -------
        k : int
            Optimal number of clusters.
        """
        if not k:
            # Check if there is only one single cluster, or no clusters
            is_single_cluster = EMXSp_Composition_Analyzer._is_single_cluster(compositions_df, verbose=self.verbose)
            if is_single_cluster or self.clustering_cfg.max_k <= 1:
                k = 1
            elif compute_k_only_once:
                # Get number of clusters (k) and optionally save the plot
                results_dir = self.analysis_dir if self.plot_cfg.save_plots else None
                k = EMXSp_Composition_Analyzer._get_k(
                    compositions_df, self.clustering_cfg.max_k, self.clustering_cfg.k_finding_method,
                    show_plot=self.plot_cfg.show_plots, results_dir=results_dir
                )
            else:
                # Calculate most frequent number of clusters (k) with elbow method. Does not save the plot
                k = EMXSp_Composition_Analyzer._get_most_freq_k(
                    compositions_df, self.clustering_cfg.max_k, self.clustering_cfg.k_finding_method, verbose=self.verbose
                )
        elif self.verbose:
            print_single_separator()
            print(f"Number of clusters was forced to be {k}")
        return k
    
    
    @staticmethod
    def _get_most_freq_k(
        compositions_df: 'pd.DataFrame',
        max_k: int,
        k_finding_method: str,
        verbose: bool = False,
        show_plot: bool = False,
        results_dir: str = None
    ) -> int:
        """
        Determine the most frequent optimal number of clusters (k) for the given compositions.
    
        This method repeatedly runs the k-finding algorithm and selects the most robust k value.
        It loops until it finds a value of k that is at least twice as frequent as the second most frequent value,
        or until a maximum number of iterations is reached.
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame containing the compositions to cluster.
        max_k : int
            Maximum number of clusters to test.
        k_finding_method : str
            Method used to determine the optimal k (passed to _get_k).
        verbose : bool, optional
            If True, print progress and summary information.
        show_plot : bool, optional
            If True, display plots for each k-finding run.
        results_dir : str, optional
            Directory to save plots/results (if applicable).
    
        Returns
        -------
        k : int
            The most robustly determined number of clusters.
    
        Raises
        ------
        ValueError
            If there are not enough data points to determine clusters.
    
        Notes
        -----
        - The function tries up to 5 times to find a dominant k value.
        - If a tie or ambiguity remains, it picks the smallest k with frequency ≥ half of the most frequent.
    
        """
        if len(compositions_df) < 2:
            raise ValueError("Not enough data points to determine clusters (need at least 2).")
    
        if verbose:
            print_single_separator()
            print("Computing number of clusters k...")
    
        k_found = []
        k = None
        max_iter = 5
        i = 0
    
        max_allowed_k = len(compositions_df) - 1  # KElbowVisualizer throws error if max_k > n_samples - 1
        if max_k > max_allowed_k:
            max_k = max_allowed_k
            warnings.warn(
                f"Maximum number of clusters reduced to {max_allowed_k} because number of clustered points is {max_allowed_k + 1}.",
                UserWarning
            )
    
        while k is None:
            i += 1
            for n in range(20):
                k_val = EMXSp_Composition_Analyzer._get_k(
                    compositions_df, max_k, k_finding_method,
                    show_plot=show_plot, results_dir=results_dir
                )
                if not isinstance(k_val, int) or k_val < 1:
                    continue  # skip invalid k values
                k_found.append(k_val)
    
            if not k_found:
                raise ValueError("No valid cluster counts were found.")
    
            counts = np.bincount(k_found)
            total = counts.sum()
            sorted_k = np.argsort(-counts)  # descending
            first_k = sorted_k[0]
            first_count = counts[first_k]
    
            if len(sorted_k) > 1:
                second_k = sorted_k[1]
                second_count = counts[second_k]
            else:
                second_k = None
                second_count = 0
    
            # Check if first is at least twice as common as second
            if second_count == 0 or first_count >= 2 * second_count:
                k = first_k
            elif i >= max_iter:  # max_iter reached
                # Pick smallest of all k values whose frequency is ≥ half of the most frequent
                threshold = first_count / 2
                k = min([k_val for k_val, count in enumerate(counts) if count >= threshold])
            else:
                k = None
    
        if verbose:
            print(f"Most frequent k: {first_k} (count = {first_count}, frequency = {first_count / total:.2%})")
            if second_k is not None and second_k != 0:
                print(f"Second most frequent k: {second_k} (count = {second_count}, frequency = {second_count / total:.2%})")
            if len(np.where(counts == first_count)[0]) > 1:
                print(f"Tie detected among: {np.where(counts == first_count)[0].tolist()} (choosing {k})")
    
        return int(k)
    
    
    @staticmethod
    def _get_k(
        compositions_df: 'pd.DataFrame',
        max_k: int = 6,
        method: str = 'silhouette',
        model: 'KMeans' = None,
        results_dir: str = None,
        show_plot: bool = False
    ) -> int:
        """
        Determine the optimal number of clusters for the data using visualizer methods.
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame containing the compositions to cluster.
        max_k : int, optional
            Maximum number of clusters to test (default: 6).
        method : str, optional
            Method for evaluating the number of clusters. One of 'elbow', 'silhouette', or 'calinski_harabasz' (default: 'silhouette').
        model : KMeans or compatible, optional
            Clustering model to use (default: KMeans(n_init='auto')).
        results_dir : str, optional
            Directory to save the plot (if provided).
        show_plot : bool, optional
            If True, show the plot interactively (default: False).
    
        Returns
        -------
        optimal_k : int
            The optimal number of clusters found.
    
        Raises
        ------
        ValueError
            If an unsupported method is provided.
    
        Notes
        -----
        - Uses yellowbrick's KElbowVisualizer.
        - For 'elbow', finds the inflection point; for 'silhouette' or 'calinski_harabasz', finds the k with the highest score.
        - If cluster finding fails, returns 1 and prints a warning.
        - If `show_plot` is True, the plot is shown interactively.
        - If `results_dir` is provided, the plot is saved as 'Elbow_plot.png'.
        """
        if model is None:
            from sklearn.cluster import KMeans
            model = KMeans(n_init='auto')
    
        # Map 'elbow' to 'distortion' for yellowbrick, but keep original for logic
        user_method = method
        if method == 'elbow':
            yb_method = 'distortion'
        elif method in ['silhouette', 'calinski_harabasz']:
            yb_method = method
        else:
            raise ValueError(f"Unsupported method '{method}' for evaluating number of clusters.")
    
        plt.figure(figsize=(10, 8))
        visualizer = KElbowVisualizer(model, k=max_k, metric=yb_method, timings=True, show=False)
    
        try:
            visualizer.fit(compositions_df)
        except ValueError as er:
            warnings.warn(f"Number of clusters could not be identified due to the following error:\n{er}\nForcing k = 1.", UserWarning)
            return 1
    
        # Get optimal number of clusters
        if user_method == 'elbow':
            optimal_k = visualizer.elbow_value_
        elif user_method in ['silhouette', 'calinski_harabasz']:
            # For silhouette and calinski_harabasz, k_scores_ is indexed from k=2
            optimal_k = np.argmax(visualizer.k_scores_) + 2
            visualizer.elbow_value_ = optimal_k  # For correct plotting
    
        # Add labels
        ax1, ax2 = visualizer.axes
        ax1.set_ylabel(f'{user_method} score')
        ax1.set_xlabel('k: number of clusters')
        ax2.set_ylabel('Fit time (sec)')
        ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    
        if show_plot:
            plt.ion()
            visualizer.show()
            plt.pause(0.001)
    
        if results_dir:
            fig = visualizer.fig
            fig.savefig(os.path.join(results_dir, 'Elbow_plot.png'))
    
        if not show_plot:
            plt.close(visualizer.fig)
    
        # Set k to 1 if elbow method was unsuccessful
        if optimal_k is None:
            optimal_k = 1
    
        return int(optimal_k)
    
    
    @staticmethod
    def _is_single_cluster(
        compositions_df: 'pd.DataFrame',
        verbose: bool = False
    ) -> bool:
        """
        Determine if the data effectively forms a single cluster using k-means and silhouette analysis.
    
        This method:
          - Fits k-means with k=1 and calculates the RMS distance from the centroid.
          - Fits k-means with k=2 multiple times, keeping the best silhouette score and inertia.
          - Uses empirically determined thresholds on silhouette score, centroid distance, and inertia ratio
            to decide if the data forms a single cluster or multiple clusters.
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to analyze.
        verbose : bool, optional
            If True, print detailed output of the clustering metrics.
    
        Returns
        -------
        is_single_cluster : bool
            True if the data is best described as a single cluster, False otherwise.
    
        Notes
        -----
        - Uses silhouette score and inertia ratio as main criteria.
        - Empirical thresholds: mean centroid distance < 0.025, silhouette < 0.5, or inertia ratio < 1.5.
        """
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        import numpy as np
    
        if verbose:
            print_single_separator()
            print('Checking if more than 1 cluster is present...')
    
        # Fit k-means for k=1
        kmeans_1 = KMeans(n_clusters=1, random_state=0, n_init='auto')
        kmeans_1.fit(compositions_df)
        inertia_1 = kmeans_1.inertia_
        rms_distance_1 = np.sqrt(inertia_1 / len(compositions_df))
    
        # Fit k-means for k=2, keep best silhouette score and inertia
        best_silhouette_score_2 = -1
        best_inertia_2 = None
        for _ in range(10):
            kmeans_2 = KMeans(n_clusters=2, random_state=None, n_init='auto')
            labels_2 = kmeans_2.fit_predict(compositions_df)
            inertia_2 = kmeans_2.inertia_
            sil_score = silhouette_score(compositions_df, labels_2)
            if sil_score > best_silhouette_score_2:
                best_silhouette_score_2 = sil_score
                best_inertia_2 = inertia_2
    
        ratio_inertias = inertia_1 / best_inertia_2 if best_inertia_2 else float('inf')
    
        if verbose:
            print(f"RMS distance for k=1: {rms_distance_1*100:.1f}%")
            print(f"Inertia for k=1: {inertia_1:.3f}")
            print(f"Inertia for k=2: {best_inertia_2:.3f}")
            print(f"Ratio of inertia for k=1 over k=2: {ratio_inertias:.2f}")
            print(f"Silhouette Score for k=2: {best_silhouette_score_2:.2f}")
    
        # Empirical decision logic
        if rms_distance_1 < 0.03:
            is_single_cluster = True
            reason_str = 'd_rms < 3%'
        elif best_silhouette_score_2 < 0.5:
            is_single_cluster = True
            reason_str = 's < 0.5'
        elif best_silhouette_score_2 > 0.6:
            is_single_cluster = False
            reason_str = 's > 0.6'
        elif ratio_inertias < 1.5:
            is_single_cluster = True
            reason_str = 'ratio of inertias < 1.5'
        else:
            is_single_cluster = False
            reason_str = 'ratio of inertias > 1.5 and 0.5 < s < 0.6'
    
        if verbose:
            if is_single_cluster:
                print(reason_str + ": The data effectively forms a single cluster.")
            else:
                print(reason_str + ": The data forms multiple clusters.") 
    
        return is_single_cluster
    #%% Clustering operations
    # ============================================================================= 
    def _run_kmeans_clustering(self, k, compositions_df):
        """
        Run k-means clustering multiple times and select the best solution by silhouette score.
    
        Returns
        -------
        kmeans : KMeans
            The best fitted KMeans instance.
        labels : np.ndarray
            Cluster labels for each composition.
        sil_score : float
            Best silhouette score obtained.
        """
        if k > 1:
            n_clustering_eval = 20  # Number of clustering evaluations to run
            best_sil_score = -np.inf  # Initialise best silhouette score
            for _ in range(n_clustering_eval):
                # K-means is not ideal for clusters with varying sizes/densities. Consider alternatives (e.g., GMM).
                is_clustering_ok = False
                max_loops_nonneg_silh = 50  # Max loops to find clustering solutions with no negative silhouette values
                n_loop = 0
                while not is_clustering_ok and n_loop < max_loops_nonneg_silh:
                    n_loop += 1
                    kmeans, labels = self._get_clustering_kmeans(k, compositions_df)
                    silhouette_vals = silhouette_samples(compositions_df, labels)
                    if np.all(silhouette_vals > 0):
                        # Clustering is accepted only if all silhouette values are positive (no wrong clustering)
                        is_clustering_ok = True
                    sil_score = silhouette_score(compositions_df, labels)
                if sil_score > best_sil_score:
                    best_kmeans = kmeans
                    best_labels = labels
                    best_sil_score = sil_score
            return best_kmeans, best_labels, best_sil_score
        else:
            # Clustering with k = 1 is trivial, and has no silhouette score
            kmeans, labels = self._get_clustering_kmeans(k, compositions_df)
            return kmeans, labels, np.nan
    

    def _prepare_composition_dataframes(self, compositions_list_at, compositions_list_w):
        """
        Convert lists of compositions to DataFrames for clustering.
    
        Returns
        -------
        compositions_df : pd.DataFrame
            DataFrame of compositions for clustering (feature set selected).
        compositions_df_other_fr : pd.DataFrame
            DataFrame of compositions in the alternate fraction representation.
        """
        # Substitute nan with 0
        if self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
            compositions_df = pd.DataFrame(compositions_list_at).fillna(0)
            compositions_df_other_fr = (pd.DataFrame(compositions_list_w)).fillna(0) 
        elif self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
            compositions_df = pd.DataFrame(compositions_list_w).fillna(0)
            compositions_df_other_fr = (pd.DataFrame(compositions_list_at)).fillna(0)
            
        return compositions_df, compositions_df_other_fr
    
    
    def _get_clustering_kmeans(
        self,
        k: int,
        compositions_df: 'pd.DataFrame'
    ) -> Tuple['KMeans', 'np.ndarray']:
        """
        Perform k-means clustering on the given compositions.
    
        Parameters
        ----------
        k : int
            The number of clusters to find.
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to cluster.
    
        Returns
        -------
        kmeans : KMeans
            The fitted KMeans object.
        labels : np.ndarray
            Array of cluster (phase) labels for each composition point.
    
        Raises
        ------
        ValueError
            If clustering is unsuccessful due to invalid data or parameters.
    
        Notes
        -----
        - Uses k-means++ initialization and scikit-learn's default settings.
        - n_init='auto' requires scikit-learn >= 1.2.0.
        """
        try:
            kmeans = KMeans(n_clusters=k, init='k-means++', n_init='auto')
            # Perform clustering. Returns labels (array of cluster (= phase) ID each composition point belongs to)
            labels = kmeans.fit_predict(compositions_df)
        except Exception as e:
            raise ValueError(f"Clustering unsuccessful due to the following error:\n{e}")
    
        return kmeans, labels

    
    def _get_clustering_dbscan(
        self,
        compositions_df: 'pd.DataFrame'
    ) -> Tuple['np.ndarray', int]:
        """
        Perform DBSCAN clustering on the given compositions.
        CURRENTLY NOT SUPPORTED
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to cluster.
    
        Returns
        -------
        labels : np.ndarray
            Array of cluster labels for each composition point. Noise points are labeled as -1.
        num_labels : int
            Number of unique clusters found (excluding noise points).
    
        Raises
        ------
        ValueError
            If clustering is unsuccessful due to invalid data or parameters.
    
        Notes
        -----
        - Uses eps=0.1 and min_samples=1 as DBSCAN parameters by default.
        - The number of clusters excludes noise points (label -1).
        """
        try:
            dbscan = DBSCAN(eps=0.1, min_samples=1)
            labels = dbscan.fit_predict(compositions_df)
        except Exception as e:
            raise ValueError(f"Clustering unsuccessful due to the following error:\n{e}")
    
        # Get the number of unique labels, excluding noise (-1)
        num_labels = len(set(labels)) - (1 if -1 in labels else 0)
    
        return labels, num_labels


    def _compute_cluster_statistics(self, compositions_df, compositions_df_other_fr, centroids, labels):
        """
        Compute statistics for each cluster, including WCSS, standard deviations, and centroids
        in terms of both atomic and mass fractions.
    
        Returns
        -------
        wcss_per_cluster : list
            Within-Cluster Sum of Squares for each cluster.
        rms_dist_cluster : list
            Standard deviation of distances to centroid for each cluster.
        rms_dist_cluster_other_fr : list
            Standard deviation of distances in alternate fraction representation.
        n_points_per_cluster : list
            Number of points in each cluster.
        els_std_dev_per_cluster : list
            Elemental standard deviations within each cluster.
        els_std_dev_per_cluster_other_fr : list
            Elemental standard deviations in alternate fraction representation.
        centroids_other_fr : list
            Centroids in alternate fraction representation.
        max_cl_rmsdist : float
            Maximum standard deviation across all clusters.
        """
        wcss_per_cluster = []
        rms_dist_cluster = []
        rms_dist_cluster_other_fr = []
        n_points_per_cluster = []
        els_std_dev_per_cluster = []
        els_std_dev_per_cluster_other_fr = []
        centroids_other_fr = []
        for i, centroid in enumerate(centroids):
            # Save data using the elemental fraction employed as feature
            cluster_points = compositions_df[labels == i].to_numpy()
            n_points_per_cluster.append(len(cluster_points))
            els_std_dev_per_cluster.append(np.std(cluster_points, axis=0, ddof=1))
            distances = np.linalg.norm(cluster_points - centroid, axis=1)
            wcss_per_cluster.append(np.sum(distances ** 2))
            rms_dist_cluster.append(np.sqrt(np.mean(distances ** 2)))
    
            # Save data also for the other elemental fraction type
            cluster_points_other = compositions_df_other_fr[labels == i].to_numpy()
            centroid_other_fr = np.mean(cluster_points_other, axis=0)
            centroids_other_fr.append(centroid_other_fr)
            els_std_dev_per_cluster_other_fr.append(np.std(cluster_points_other, axis=0, ddof=1))
            distances_other_fr = np.linalg.norm(cluster_points_other - centroid_other_fr, axis=1)
            rms_dist_cluster_other_fr.append(np.sqrt(np.mean(distances_other_fr ** 2)))
            
        max_cl_rmsdist = max(rms_dist_cluster)
        
        return (wcss_per_cluster, rms_dist_cluster, rms_dist_cluster_other_fr, n_points_per_cluster,
                els_std_dev_per_cluster, els_std_dev_per_cluster_other_fr, centroids_other_fr, max_cl_rmsdist)
    
    
    #%% Data compositional analysis
    # =============================================================================     
    def analyse_data(self, max_analytical_error_percent, k=None, compute_k_only_once=False):
        """
        Analyse quantified spectra, perform clustering, assign reference phases and mixtures, and save results.
    
        This function orchestrates the workflow:
          1. Selects good compositions for clustering.
          2. Prepares DataFrames for clustering.
          3. Determines the optimal number of clusters (k).
          4. Runs clustering and computes cluster statistics.
          5. Assigns reference phases and detects mixtures.
          6. Saves results and related plots.
    
        Parameters
        ----------
        max_analytical_error : float or None
            Maximum allowed analytical error for a composition to be considered valid, expressed as w%.
        k : int, optional
            Number of clusters to use (if not provided, determined automatically).
        compute_k_only_once : bool, optional
            If True, compute k only once; otherwise, use the most frequent k.
    
        Returns
        -------
        success : bool
            True if analysis was successful, False otherwise.
        max_cl_rmsdist : float
            Maximum standard deviation across clusters.
        min_conf : float or None
            Minimum confidence among assigned reference phases.
        """
        # 1. Select compositions to use for clustering
        max_analytical_error = max_analytical_error_percent / 100
        (compositions_list_at, compositions_list_w, unused_compositions_list,
         df_indices, n_datapts) = self._select_good_compositions(max_analytical_error)
        n_datapts_used = len(compositions_list_at)
    
        if n_datapts_used < 5:
            print_single_separator()
            print(f"Only {n_datapts_used} spectra were considered 'good', but a minimum of 5 data points are required for clustering.")
            # Print additional messages with how many spectra were discarded for which reason
            self._report_n_discarded_spectra(n_datapts, max_analytical_error)
            return False, 0, 0  # zeroes are placeholders
    
        if self.verbose:
            print_single_separator()
            print('Spectra selection:')
            print(f"{n_datapts_used} data points are used, out of {n_datapts} collected spectra.")
            self._report_n_discarded_spectra(n_datapts, max_analytical_error)
    
        # 2. Make analysis directory to save results
        self._make_analysis_dir()
    
        # 3. Prepare DataFrames for clustering
        compositions_df, compositions_df_other_fr = self._prepare_composition_dataframes(compositions_list_at, compositions_list_w)

        # 4. Perform clustering
        if k is None:
            # Extract k value from configurations (None if not provided)
            k = self.clustering_cfg.k
        if self.clustering_cfg.method == 'kmeans':
            k = self._find_optimal_k(compositions_df, k, compute_k_only_once)
            kmeans, labels, sil_score = self._run_kmeans_clustering(k, compositions_df)
            centroids = kmeans.cluster_centers_
            wcss = kmeans.inertia_
        elif self.clustering_cfg.method == 'dbscan':
            # labels, num_labels = self._get_clustering_dbscan(compositions_df)
            print('Clustering via DBSCAN is not implemented yet')
            return False, 0, 0  # zeroes are placeholders
    
        # 5. Compute cluster statistics
        (wcss_per_cluster, rms_dist_cluster, rms_dist_cluster_other_fr, 
         n_points_per_cluster, els_std_dev_per_cluster, els_std_dev_per_cluster_other_fr,
         centroids_other_fr, max_cl_rmsdist) = self._compute_cluster_statistics(
            compositions_df, compositions_df_other_fr, centroids, labels
        )
    
        # 6. Assign reference phases
        min_conf, max_raw_confs, refs_assigned_df = self._assign_reference_phases(centroids, rms_dist_cluster)
    
        # 7. Assign mixtures
        clusters_assigned_mixtures = self._assign_mixtures(
            k, labels, compositions_df, rms_dist_cluster, max_raw_confs, n_points_per_cluster
        )
    
        # 8. Save and store results
        if self.is_acquisition:
            # When collecting, save collected spectra, their quantification, and to which cluster they are assigned
            self._save_collected_data(labels, df_indices, backup_previous_data=True, include_spectral_data=True)
        else:
            # During analysis of Data.csv, save the compositions, together with their assigned phases, in the Analysis folder
            self._save_collected_data(labels, df_indices, backup_previous_data=True, include_spectral_data=False)
    
        self._save_result_and_stats(
            centroids, els_std_dev_per_cluster, centroids_other_fr, els_std_dev_per_cluster_other_fr,
            n_points_per_cluster, wcss_per_cluster, rms_dist_cluster, rms_dist_cluster_other_fr,
            refs_assigned_df, wcss, sil_score, n_datapts, max_analytical_error, clusters_assigned_mixtures
        )
    
        # 9. Save plots
        if self.plot_cfg.save_plots:
            self._save_plots(kmeans, compositions_df, centroids, labels, els_std_dev_per_cluster, unused_compositions_list)
    
        return True, max_cl_rmsdist, min_conf
    
    
    def _select_good_compositions(self, max_analytical_error):
        """
        Select compositions for clustering, filtering out those with high analytical error or bad quantification flags.
    
        Returns
        -------
        compositions_list_at : list
            List of atomic fractions for good spectra.
        compositions_list_w : list
            List of mass fractions for good spectra.
        unused_compositions_list : list
            List of compositions not used for clustering (for plotting).
        df_indices : list
            Indices of rows used for phase identification.
        n_datapts : int
            Total number of spectra considered.
        """
        # Initialise counters for spectra filtered out
        self.n_sp_too_low_counts = 0
        self.n_sp_too_high_an_err = 0
        self.n_sp_bad_quant = 0
    
        compositions_list_at = []
        compositions_list_w = []
        unused_compositions_list = []
        df_indices = []
        n_datapts = len(self.spectra_quant)
    
        for i in range(n_datapts):
            if self.spectra_quant[i] is not None:
                is_comp_ok = True
                spectrum_quant_result_at = self.spectra_quant[i][cnst.COMP_AT_FR_KEY]
                spectrum_quant_result_w = self.spectra_quant[i][cnst.COMP_W_FR_KEY]
                analytical_error = self.spectra_quant[i][cnst.AN_ER_KEY]
                quant_flag = self.spectral_data[cnst.QUANT_FLAG_DF_KEY][i]
    
                # Check if composition was flagged as bad during quantification
                if quant_flag not in self.clustering_cfg.quant_flags_accepted:
                    is_comp_ok = False
                    self.n_sp_bad_quant += 1
    
                elif max_analytical_error is None:
                    # Analytical error check is disabled
                    pass
    
                # Check if analytical error is too high
                elif analytical_error < - (max_analytical_error + self.undetectable_an_er) or analytical_error > max_analytical_error:
                    is_comp_ok = False
                    self.n_sp_too_high_an_err += 1
    
                # Append composition to list of used or unused datapoints
                if is_comp_ok:
                    df_indices.append(i)
                    # Construct dictionary that includes all elements that are supposed to be in the sample.
                    comp_at = {el: spectrum_quant_result_at.get(el, 0) for el in self.all_els_sample}
                    compositions_list_at.append(comp_at)
                    comp_w = {el: spectrum_quant_result_w.get(el, 0) for el in self.all_els_sample}
                    compositions_list_w.append(comp_w)
                else:
                    # Collect unused data points to show them in the clustering plot
                    if self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
                        comp = [spectrum_quant_result_at.get(el, 0) for el in self.all_els_sample]
                    elif self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
                        comp = [spectrum_quant_result_w.get(el, 0) for el in self.all_els_sample]
                    unused_compositions_list.append(comp)
            else:
                self.n_sp_too_low_counts += 1
    
        return compositions_list_at, compositions_list_w, unused_compositions_list, df_indices, n_datapts


    def _correlate_centroids_to_refs(
        self,
        centroids: 'np.ndarray',
        cluster_radii: 'np.ndarray',
        ref_phases_df: 'pd.DataFrame'
    ) -> Tuple[List[float], 'pd.DataFrame']:
        """
        Correlate each cluster centroid to reference phases and compute confidence scores.
    
        For each centroid, selects all reference phases within a hypersphere of radius
        max(0.1, 5 * cluster_radius) around the centroid, and computes a confidence score
        for each reference. The highest confidence for each cluster is stored.
    
        Parameters
        ----------
        centroids : np.ndarray, shape (n_clusters, n_features)
            Array of cluster centroids (in elemental fraction space).
        cluster_radii : np.ndarray, shape (n_clusters,)
            Array of standard deviations (radii) for each cluster.
        ref_phases_df : pd.DataFrame
            DataFrame where each row is a reference phase (elemental fractions).
    
        Returns
        -------
        max_raw_confs : list of float
            List of maximum confidence scores for each cluster.
        refs_assigned_df : pd.DataFrame
            DataFrame with reference names and their confidences for each cluster.
    
        Notes
        -----
        - Only reference phases within 5 times the cluster radius (or at least 0.1) from the centroid are considered.
        - Confidence is computed using EMXSp_Composition_Analyzer._get_ref_confidences.
        """
        # Get all reference phase compositions as a numpy array
        all_ref_phases = ref_phases_df.to_numpy()
        refs_dict = []        # For DataFrame of references and their confidences
        max_raw_confs = []    # For convergence checks
    
        # For each cluster, assign to reference(s) if present and calculate confidence
        for centroid, radius in zip(centroids, cluster_radii):
            # Calculate distances from centroid to each reference phase
            distances = np.linalg.norm(all_ref_phases - centroid, axis=1)
            # Select all reference phases within 5*radius (min 0.1) of centroid
            indices = np.where(distances < max(0.1, 5 * radius))[0]
            # Get chemical formulae and compositions of selected reference phases
            ref_names = [self.ref_formulae[i] for i in indices]
            ref_phases = [all_ref_phases[i] for i in indices]
            # Calculate confidences based on distance between centroid and reference
            max_raw_conf, refs_dict_row = EMXSp_Composition_Analyzer._get_ref_confidences(
                centroid, ref_phases, ref_names
            )
            # Store maximum confidence for this cluster
            max_raw_confs.append(max_raw_conf)
            # Store dictionary of reference names and confidences for this cluster
            refs_dict.append(refs_dict_row)
    
        # Create DataFrame with information on reference phases assigned to clusters
        refs_assigned_df = pd.DataFrame(refs_dict)
    
        return max_raw_confs, refs_assigned_df


    def _assign_reference_phases(self, centroids, rms_dist_cluster):
        """
        Assign reference phases to clusters if reference formulae are provided.
    
        Returns
        -------
        min_conf : float or None
            Minimum confidence among all clusters assigned to a reference.
        max_raw_confs : list or None
            Maximum raw confidence scores for each cluster.
        refs_assigned_df : pd.DataFrame or None
            DataFrame of reference assignments.
        """
        min_conf = None
        max_raw_confs = None
        refs_assigned_df = None
        if self.ref_formulae is not None:
            # Correlate calculated centroids to the reference phases
            max_raw_confs, refs_assigned_df = self._correlate_centroids_to_refs(
                centroids, rms_dist_cluster, self.ref_phases_df
            )
            # Get lowest value among the highest confidences assigned to each cluster, used for convergence
            if len(max_raw_confs) > 0:
                max_confs_num = [conf for conf in max_raw_confs if conf is not None]
                if len(max_confs_num) > 0:
                    min_conf = min(max_confs_num)
        return min_conf, max_raw_confs, refs_assigned_df


    @staticmethod
    def _get_ref_confidences(
        centroid: 'np.ndarray',
        ref_phases: 'np.ndarray',
        ref_names: List[str]
    ) -> Tuple[Optional[float], Dict]:
        """
        Compute confidence scores for reference phases near a cluster centroid.
    
        For each reference phase within a cluster's neighborhood, this function:
          - Computes the Euclidean distance to the centroid.
          - Assigns a confidence score using a Gaussian function of the distance.
          - Reduces confidences if multiple references are nearby, to account for ambiguity.
          - Returns a dictionary of references and their confidences (sorted by confidence), and the maximum raw confidence.
    
        Parameters
        ----------
        centroid : np.ndarray
            Cluster centroid in feature space (shape: n_features,).
        ref_phases : np.ndarray
            Array of reference phase compositions (shape: n_refs, n_features).
        ref_names : list of str
            Names of reference phases.
    
        Returns
        -------
        max_raw_conf : float or None
            Maximum confidence score among the references, or None if no reference is close.
        refs_dict : Dict
            Dictionary of reference names and their confidences, sorted by confidence.
            Keys: 'Ref1', 'Conf1', 'Ref2', 'Conf2', etc.
    
        Notes
        -----
        - Only confidences above 1% are included in the output dictionary.
        - The confidence spread (sigma) is set to 0.03 for the Gaussian.
        - Nearby references reduce each other's confidence using a secondary Gaussian weighting.
        """
        if ref_phases == [] or len(ref_phases) == 0:
            # No reference phase is close enough to the centroid
            refs_dict = {'Ref1': np.nan, 'Conf1': np.nan}
            max_raw_conf = None
        else:
            # Calculate distances from centroid to each reference phase
            distances = np.linalg.norm(ref_phases - centroid, axis=1)
    
            # Assign confidence using a Gaussian function (sigma = 0.03)
            confidences = np.exp(-distances**2 / (2 * 0.03**2))
    
            # Reduce confidences for ambiguity if multiple references are close
            weights_conf = np.exp(-(1 - confidences)**2 / (2 * 0.3**2))
            weights_conf /= np.sum(weights_conf)  # Normalize
    
            # Adjust confidences by their weights
            confidences *= weights_conf
    
            # Get maximum raw confidence
            max_raw_conf = float(np.max(confidences))
    
            # Sort references by confidence, descending
            sorted_indices = np.argsort(-confidences)
            sorted_ref_names = np.array(ref_names)[sorted_indices]
            sorted_confidences = confidences[sorted_indices]
    
            # Build dictionary of references and confidences (only those > 1%)
            refs_dict = {}
            for i, (ref_name, conf) in enumerate(zip(sorted_ref_names, sorted_confidences)):
                if conf > 0.01:
                    refs_dict[f'Ref{i+1}'] = ref_name
                    refs_dict[f'Conf{i+1}'] = np.round(conf, 2)
    
        return max_raw_conf, refs_dict


    #%% Binary cluster decomposition
    # =============================================================================       
    def _assign_mixtures(self, k, labels, compositions_df, rms_dist_cluster, max_raw_confs, n_points_per_cluster):
        """
        Determine if clusters are mixtures or single phases, using reference phases and NMF if needed.
    
        Returns
        -------
        clusters_assigned_mixtures : list
            List of mixture assignments for each cluster.
            
        Potential improvements
        ----------------------
        Instead of using the cluster standard deviation, use covariance of elemental fractions
        to discern clusters that may originate from binary phase mixtures or solid solutions.
        """
        clusters_assigned_mixtures = []
        for i in range(k):
            # Get compositions of data points included in cluster as np.array (only detectable elements)
            cluster_data = compositions_df[self.detectable_els_sample].iloc[labels == i].values
            max_mix_conf = 0
            mixtures_dicts = []
            
            # # Use log-ratio transformations, which map the data from the simplex to real Euclidean space
            # if len(cluster_data) > 1:
            #     # Suppose X is your n × m array of normalized compositions
            #     X_clr = clr(cluster_data + 1e-6) # to avoid zeroes
            #     # Compute covariance matrix on CLR-transformed data
            #     cov_clr = np.cov(X_clr.T, rowvar=True)
            #     print(cov_clr)
                
            #     # 5. Compute correlation matrix
            #     std_dev = np.sqrt(np.diag(cov_clr))
            #     corr_matrix = cov_clr / np.outer(std_dev, std_dev)
                
            #     # 6. Plot covariance heatmap
            #     plt.figure(figsize=(10, 4))
                
            #     plt.subplot(1, 2, 1)
            #     sns.heatmap(cov_clr, xticklabels=self.detectable_els_sample, yticklabels=self.detectable_els_sample,
            #                 cmap='coolwarm', center=0, annot=True, fmt=".3f")
            #     plt.title('Covariance matrix (CLR space)')
                
            #     # 7. Plot correlation heatmap
            #     plt.subplot(1, 2, 2)
            #     sns.heatmap(corr_matrix, xticklabels=self.detectable_els_sample, yticklabels=self.detectable_els_sample,
            #                 cmap='coolwarm', center=0, annot=True, fmt=".3f")
            #     plt.title('Correlation matrix (CLR space)')
                
            #     plt.tight_layout()
            #     plt.show()
            
            max_rmsdist_single_cluster = 0.03
            if rms_dist_cluster[i] < max_rmsdist_single_cluster:
                if max_raw_confs is None or len(max_raw_confs) < 1:
                    is_cluster_single_phase = n_points_per_cluster[i] > 3
                elif max_raw_confs[i] is not None and max_raw_confs[i] > 0.5:
                    is_cluster_single_phase = True
                else:
                    is_cluster_single_phase = False
            else:
                is_cluster_single_phase = False
    
            if is_cluster_single_phase:
                # Cluster determined to stem from a single phase
                pass
            elif len(self.ref_formulae) > 1:
                max_mix_raw_conf, mixtures_dicts = self._identify_mixture_from_refs(cluster_data)
                max_mix_conf = max(max_mix_conf, max_mix_raw_conf)
            if not is_cluster_single_phase and max_mix_conf < 0.5:
                mix_nmf_conf, mixture_dict = self._identify_mixture_nmf(cluster_data)
                if mixture_dict is not None:
                    mixtures_dicts.append(mixture_dict)
                max_mix_conf = max(max_mix_conf, mix_nmf_conf)
            clusters_assigned_mixtures.append(mixtures_dicts)
        return clusters_assigned_mixtures
    
    
    def _identify_mixture_from_refs(self, X: 'np.ndarray') -> Tuple[float, List[Dict]]:
        """
        Identify mixtures within a cluster by testing all pairs of reference phases using constrained optimization.
    
        For each possible pair of reference phases, tests if the cluster compositions (X)
        can be well described by a linear combination of the two reference phases, using
        non-negative matrix factorization (NMF) with fixed bases.
    
        Parameters
        ----------
        X : np.ndarray
            Cluster data (compositions), shape (n_samples, n_features).
    
        Returns
        -------
        max_confidence : float
            The highest confidence score among all tested mixtures.
        mixtures_dicts : list of Dict
            List of mixture descriptions for all successful reference pairs.
    
        Notes
        -----
        - Each mixture is described by a dictionary, as returned by _get_mixture_dict_with_conf.
        - Only pairs with acceptable reconstruction error are included.
        - The confidence metric and acceptance criteria are defined in _get_mixture_dict_with_conf.
        """
        # Generate all possible pairs of reference phases
        ref_pair_combinations = list(itertools.combinations(range(len(self.ref_phases_df)), 2))
    
        mixtures_dicts = []
        max_confidence = 0
    
        for ref_comb in ref_pair_combinations:
            # Get the names of the reference phases in this pair
            ref_names = [self.ref_formulae[ref_i] for ref_i in ref_comb]
    
            # Ratio of weights of references, for molar concentrations of parent phases
            ref_w_r = self.ref_weights_in_mixture[ref_comb[0]] / self.ref_weights_in_mixture[ref_comb[1]]
    
            # Get matrix of basis vectors (H) for the two reference phases
            H = np.array([
                self.ref_phases_df[self.detectable_els_sample].iloc[ref_i].values
                for ref_i in ref_comb
            ])
            
            # Perform NMF with fixed H to fit the cluster data as a mixture of the two reference phases
            W, _ = self._nmf_with_constraints(X, n_components=2, fixed_H=H)
    
            # Compute reconstruction error for the fit
            recon_er = self._calc_reconstruction_error(X, W, H)
    
            # If the pair yields an acceptable reconstruction error, store the result
            pair_dict, conf = self._get_mixture_dict_with_conf(W, ref_w_r, recon_er, ref_names)
            if pair_dict is not None:
                mixtures_dicts.append(pair_dict)
                max_confidence = max(max_confidence, conf)
    
        return max_confidence, mixtures_dicts


    def _calc_reconstruction_error(
        self,
        X: 'np.ndarray',
        W: 'np.ndarray',
        H: 'np.ndarray'
    ) -> float:
        """
        Calculate the reconstruction error for a matrix factorization X ≈ W @ H.
    
        The error metric is an exponential penalty (with parameter alpha) applied to the
        absolute difference between X and its reconstruction W @ H, normalized by the
        number of elements in X. This penalizes large deviations more strongly.
    
        Parameters
        ----------
        X : np.ndarray
            Original data matrix of shape (m, n).
        W : np.ndarray
            Weight matrix of shape (m, k).
        H : np.ndarray
            Basis matrix of shape (k, n).
    
        Returns
        -------
        normalized_norm : float
            The normalized exponential reconstruction error.
    
        Notes
        -----
        - The penalty parameter alpha is set to 15 by default.
        """
        # Compute the approximation WH
        WH = np.dot(W, H)
    
        # Compute the Frobenius norm of the difference (X - WH), using an exponential form to penalize deviations more strongly
        alpha = 15
        norm = np.sum(np.exp(alpha * np.abs(X - WH)) - 1)
    
        # Get dimensions of the matrix X
        m, n = X.shape
    
        # Normalize the error by the number of entries
        normalized_norm = norm / (m * n)
    
        return normalized_norm
    
    
    def _get_mixture_dict_with_conf(
        self,
        W: 'np.ndarray',
        ref_w_r: float,
        reconstruction_error: float,
        ref_names: List[str]
    ) -> Tuple[Optional[Dict], float]:
        """
        Evaluate if a cluster is a mixture of two reference phases, and compute a confidence score.
    
        If the reconstruction error is below a set threshold, computes a confidence score and
        transforms the NMF coefficients into molar fractions. Returns a dictionary describing
        the mixture and the confidence score.
    
        Parameters
        ----------
        W : np.ndarray
            Matrix of NMF coefficients for each point in the cluster (shape: n_points, 2).
        ref_w_r : float
            Ratio of weights of the two reference phases (for molar concentration conversion).
        reconstruction_error : float
            Reconstruction error for the mixture fit.
        ref_names : list of str
            Names of the two reference phases.
    
        Returns
        -------
        mixture_dict : Dict or None
            Dictionary with mixture information if fit is acceptable, else None.
            Keys: 'refs', 'conf_score', 'mean', 'stddev'.
        conf : float
            Confidence score for this mixture (0 if not acceptable).
    
        Notes
        -----
        - The minimum acceptable reconstruction error is set empirically to 2 (a.u.).
        - Confidence is calculated as a Gaussian function of the reconstruction error (sigma=0.5).
        - Molar fractions are derived from NMF coefficients and normalized.
        - Only the first component's mean and stddev are returned in the dictionary.
        """
        # Set a minimum reconstruction error threshold for accepting mixtures
        min_acceptable_recon_error = 2  # Empirically determined
        
        save_violin_plot = self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE and self.powder_meas_cfg.is_known_powder_mixture_meas
        
        if reconstruction_error < min_acceptable_recon_error or save_violin_plot:
            # Calculate confidence score: 0.66 when error is 0.5 (empirical)
            gauss_sigma = 0.5
            conf = np.exp(-reconstruction_error**2 / (2 * gauss_sigma**2))
    
            # Transform NMF coefficients into molar fractions (see documentation for derivation)
            W_mol_frs = []
            for c1, c2 in W:
                # x1, x2 are the molar fractions
                x2 = c2 * ref_w_r / (1 - c2 * (1 - ref_w_r))
                x1 = c1 * (1 + x2 * (1 / ref_w_r - 1))
                W_mol_frs.append([x1, x2])
            W_mol_frs = np.array(W_mol_frs)
    
            # Calculate mean and standard deviation of molar fractions (not normalized)
            mol_frs_norm_means = np.mean(W_mol_frs, axis=0)
            mol_frs_norm_stddevs = np.std(W_mol_frs, axis=0)
            
            if save_violin_plot:
                self._save_violin_plot_powder_mixture(W_mol_frs, ref_names)
            
            # Store mixture information
            mixture_dict = {
                cnst.REF_NAME_KEY: ref_names,
                cnst.CONF_SCORE_KEY: conf,
                cnst.MOLAR_FR_MEAN_KEY: mol_frs_norm_means[0],
                cnst.MOLAR_FR_STDEV_KEY: mol_frs_norm_stddevs[0]
            }
        else:
            mixture_dict = None
            conf = 0
    
        return mixture_dict, conf
    
    
    def _nmf_with_constraints(
        self,
        X: 'np.ndarray',
        n_components: int,
        fixed_H: 'np.ndarray' = None
    ) -> Tuple['np.ndarray', 'np.ndarray']:
        """
        Perform Non-negative Matrix Factorization (NMF) with optional constraints on the factor matrices.
    
        This function alternates between optimizing two non-negative matrices W and H, such that X ≈ W @ H:
          - If H is fixed (provided via fixed_H), only W is updated.
          - If H is not fixed, both W and H are updated via alternating minimization.
    
        Constraints:
          - Both W and H are non-negative.
          - The rows of both W (sum of coefficients) and H (sum of elemental fractions) sum to 1.
          - Sparsity regularization (L1) is applied to H when it is updated, to favor bases with limited elements.
    
        Parameters
        ----------
        X : np.ndarray, shape (m, n)
            The input matrix to be factorized, where m is the number of samples and n is the number of elements.
        n_components : int
            The number of latent components (rank of the decomposition). For binary mixtures, n_components=2.
        fixed_H : np.ndarray or None, optional
            If provided, a fixed basis matrix of shape (n_components, n). If None, H is updated during optimization.
    
        Returns
        -------
        W : np.ndarray, shape (m, n_components)
            The non-negative coefficient matrix learned during the factorization.
        H : np.ndarray, shape (n_components, n)
            The non-negative basis matrix learned during the factorization.
            If fixed_H is provided, this matrix is not modified.
    
        Notes
        -----
        - Uses alternating minimization, solving for one matrix while keeping the other fixed.
        - Convergence is based on the Frobenius norm of the change in W and H between iterations.
        - Stops when the change is smaller than the specified tolerance (convergence_tol = 1e-3), or max_iter is reached.
        - Regularization is applied to H (if it is updated) to encourage sparsity and avoid all elements being present in both parent phases.
        """
        max_iter = 1000
        convergence_tol = 1e-3  # Algorithm converges when change in coefficients or el_fr is less than 0.1%
        lambda_H = 0  # Regularization parameter for sparsity in H. Set >0 to favor sparse basis matrix.
    
        # Initialize W and H with non-negative random values if not provided
        W = np.random.rand(X.shape[0], n_components)
        if fixed_H is None:
            H = np.random.rand(n_components, X.shape[1])  # Initialize H if not provided
        else:
            H = fixed_H
    
        prev_W, prev_H = np.inf, np.inf
        convergence = np.inf
        i = 0
    
        while convergence > convergence_tol and i < max_iter:
            # Solve for W with H fixed (or fixed_H provided)
            W_var = cp.Variable((X.shape[0], n_components), nonneg=True)
            objective_W = cp.Minimize(cp.norm(X - W_var @ H, 'fro'))
            constraints_W = [cp.sum(W_var, axis=1) == 1]  # Constrain rows of W to sum to 1
            problem_W = cp.Problem(objective_W, constraints_W)
            problem_W.solve()
            W = W_var.value  # Update W
    
            # If H is not fixed, solve for H as well (alternating minimization)
            if fixed_H is None:
                H_var = cp.Variable((n_components, X.shape[1]), nonneg=True)
                objective_H = cp.Minimize(cp.norm(X - W @ H_var, 'fro') + lambda_H * cp.norm(H_var, 1))
                constraints_H = [cp.sum(H_var, axis=1) == 1]  # Constrain rows of H to sum to 1
                problem_H = cp.Problem(objective_H, constraints_H)
                problem_H.solve()
                H = H_var.value  # Update H
    
            # Compute convergence based on the changes in W and H
            convergence_W = np.linalg.norm(W - prev_W, 'fro')
            convergence_H = np.linalg.norm(H - prev_H, 'fro') if fixed_H is not None else 0
            convergence = max(convergence_W, convergence_H)
    
            prev_W, prev_H = W, H
            i += 1
    
        return W, H
        
    
    def _identify_mixture_nmf(
        self,
        X: 'np.ndarray',
        n_components: int = 2
    ) -> Tuple[float, Optional[Dict]]:
        """
        Identify a mixture within a cluster using unconstrained NMF (Non-negative Matrix Factorization).
    
        This method fits the cluster data X to n_components using NMF with constraints (rows of W and H sum to 1),
        evaluates the reconstruction error, and if acceptable, returns a dictionary describing the mixture and a confidence score.
    
        Parameters
        ----------
        X : np.ndarray
            Cluster data (compositions), shape (n_samples, n_features).
        n_components : int, optional
            Number of components (phases) to fit (default: 2).
    
        Returns
        -------
        conf : float
            Confidence score for the mixture (0 if not acceptable).
        mixture_dict : Dict or None
            Dictionary describing the mixture if reconstruction is acceptable, else None.
    
        Notes
        -----
        - The confidence and mixture dictionary are computed using _get_mixture_dict_with_conf.
        - Elemental fractions lower than 0.5% are set to 0 using _get_pretty_formulas_nmf.
        
        Potential improvements
        ----------------------
        Consider re-running algorithms with fixed zeroes after purifying compositions
        """
        mixture_dict = None
        conf = 0
    
        # Run NMF, constraining coefficients and values of bases to add up to 1
        W, H = self._nmf_with_constraints(X, n_components)
    
        # Compute the reconstruction error
        recon_er = self._calc_reconstruction_error(X, W, H)
    
        # Get human-readable formulas and weights for the NMF bases
        ref_names, ref_weights = self._get_pretty_formulas_nmf(H, n_components)
    
        # Calculate ratio of reference weights, needed to compute molar fractions of parent phases
        ref_w_r = ref_weights[0] / ref_weights[1]

        # If pair of bases yields an acceptable reconstruction error, store the mixture info
        mixture_dict, conf = self._get_mixture_dict_with_conf(W, ref_w_r, recon_er, ref_names)  # Returns (None, 0) if error is too high

        return conf, mixture_dict
    
    
    def _get_pretty_formulas_nmf(
        self,
        phases: 'np.ndarray',
        n_components: int
    ) -> Tuple[List[str], List[float]]:
        """
        Generate human-readable (pretty) formulas from NMF bases, accounting for data noise.
    
        For each component, filters out small fractions, constructs a composition dictionary,
        and returns a formula string and a weight or atom count, depending on the clustering feature.
    
        Parameters
        ----------
        phases : np.ndarray
            Array of shape (n_components, n_elements), each row is a basis vector from NMF.
        n_components : int
            Number of components (phases) to process.
    
        Returns
        -------
        ref_names : list of str
            List of pretty formula strings for each phase.
        ref_weights : list of float
            List of weights (for mass fractions) or atom counts (for atomic fractions), for each phase.
    
        Notes
        -----
        - Fractions below 0.5% are set to zero for formula construction.
        - Uses `Composition` to generate formulas and weights.
        - For mass fractions, the weight of the phase is used.
        - For atomic fractions, the total atom count in the formula is used.
        """
        ref_names = []
        ref_weights = []
    
        for i in range(n_components):
            # Filter out too small fractions (set <0.5% to zero)
            frs = phases[i, :].copy()
            frs[frs < 0.005] = 0
    
            # Build dictionary for the parent phase composition
            fr_dict = {el: fr for el, fr in zip(self.detectable_els_sample, frs)}
    
            # Generate a Composition object based on the selected feature type
            if self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
                comp = Composition().from_weight_dict(fr_dict)
            elif self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
                comp = Composition(fr_dict)
    
            # Get integer formula and construct a pretty formula
            formula = comp.get_integer_formula_and_factor()[0]
            ref_integer_comp = Composition(formula)
            min_at_n = min(ref_integer_comp.get_el_amt_dict().values())
            pretty_at_frs = {el: round(n / min_at_n, 1) for el, n in ref_integer_comp.get_el_amt_dict().items()}
            pretty_comp = Composition(pretty_at_frs)
            pretty_formula = pretty_comp.formula
            ref_names.append(pretty_formula)
    
            # Store weight or atom count for the phase
            if self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
                ref_weights.append(pretty_comp.weight)
            elif self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
                n_atoms_in_formula = sum(pretty_comp.get_el_amt_dict().values())
                ref_weights.append(n_atoms_in_formula)
    
        return ref_names, ref_weights
    

    def _build_mixtures_df(
        self,
        clusters_assigned_mixtures: List[List[Dict]]
    ) -> 'pd.DataFrame':
        """
        Build a DataFrame summarizing mixture assignments for each cluster.
    
        For each cluster, sorts mixture dictionaries by confidence score and extracts:
          - Reference phase names (as a comma-separated string)
          - Confidence score
          - Molar ratio (mean / (1 - mean))
          - Mean and standard deviation of the main component's molar fraction
    
        Parameters
        ----------
        clusters_assigned_mixtures : list of list of Dict
            Outer list: clusters; inner list: mixture dictionaries for each cluster.
    
        Returns
        -------
        mixtures_df : pd.DataFrame
            DataFrame summarizing mixture assignments for each cluster.
            Columns: Mix1, ConfMix1, Mol_Ratio1, Icomp_Mol_Fr_Mean1, Stddev1, etc.
    
        Notes
        -----
        - If no mixtures are assigned for a cluster, the entry will be an empty dictionary.
        - The DataFrame is intended for addition to Clusters.csv or similar summary files.
        """
        mixtures_strings_dict = []
        for mixtures_dict in clusters_assigned_mixtures:
            if mixtures_dict:
                # Sort mixture dictionaries by decreasing confidence score
                sorted_mixtures = sorted(mixtures_dict, key=lambda x: -x[cnst.CONF_SCORE_KEY])
                cluster_mix_dict = {}
                for i, mixture_dict in enumerate(sorted_mixtures, start=1):
                    cluster_mix_dict[f'Mix{i}'] = ', '.join(mixture_dict[cnst.REF_NAME_KEY])
                    cluster_mix_dict[f'ConfMix{i}'] = np.round(mixture_dict[cnst.CONF_SCORE_KEY], 2)
                    cluster_mix_dict[f'Mol_Ratio{i}'] = np.round(mixture_dict[cnst.MOLAR_FR_MEAN_KEY] / (1 - mixture_dict[cnst.MOLAR_FR_MEAN_KEY]), 2)
                    cluster_mix_dict[f'Icomp_Mol_Fr_Mean{i}'] = np.round(mixture_dict[cnst.MOLAR_FR_MEAN_KEY], 2)
                    cluster_mix_dict[f'Stddev{i}'] = np.round(mixture_dict[cnst.MOLAR_FR_STDEV_KEY], 2)
                mixtures_strings_dict.append(cluster_mix_dict)
            else:
                mixtures_strings_dict.append({})
    
        # Create DataFrame to add to Clusters.csv
        mixtures_df = pd.DataFrame(mixtures_strings_dict)
    
        return mixtures_df
    
    
    #%% Run algorithms
    # =============================================================================     
    def run_collection_and_quantification(
        self,
        quantify: bool = True,
    ) -> Tuple[bool, bool]:
        """
        Perform iterative collection (and optional quantification) of spectra, followed by phase analysis and convergence check.
    
        This method:
          - Iteratively collects spectra in batches (and quantifies them if `quantify` is True).
          - After each batch, saves collected data and (if quantification is enabled) performs phase analysis (clustering).
          - Checks for convergence based on clustering statistics and confidence.
          - Stops early if convergence is achieved and minimum spectra is reached, or if no more particles are available.
    
        Parameters
        ----------
        quantify : bool, optional
            If True (default), spectra are quantified after each batch and clustering is performed.
            If False, only spectra collection is performed; quantification and clustering are skipped.
    
        Returns
        -------
        is_analysis_successful : bool
            if quantify == True: True if analysis was successful, False otherwise.
            if quantify == False: True if collection of target number of spectra was successful, False otherwise.
        is_converged : bool
            True if phase identification converged to acceptable errors, False otherwise.
    
        Notes
        -----
        - During experimental standard collection, "quantify" in fact determines whether spectra are "fitted" in-situ
        - Saves data after each batch to prevent data loss.
        - Prints a summary and processing times at the end.
        """
        tot_n_spectra = 0  # Total number of collected spectra
        max_n_sp_per_iter = 10  # Max spectra to collect per iteration (for saving in between)
        tot_spectra_to_collect = self.max_n_spectra
        n_spectra_to_collect = min(max_n_sp_per_iter, tot_spectra_to_collect, self.min_n_spectra)
        is_converged = False
        is_analysis_successful = False
        is_acquisition_successful = True
        is_exp_std_measurement = self.exp_stds_cfg.is_exp_std_measurement
        self._initialise_std_dict() # Initialise dictionary of standards to (optionally) pass onto XSp_Quantifier. Only used with known powder mixtures
        
        if quantify:
            if is_exp_std_measurement:
                quant_str = ' and fitting'
            else:
                quant_str = ' and quantification'
        else:
            quant_str = ''

        
        if self.verbose:
            print_double_separator()
            print(f"Starting collection{quant_str} of {tot_spectra_to_collect} spectra.")
        
        while tot_n_spectra < tot_spectra_to_collect:
            if self.verbose:
                print_double_separator()
                print(f"Collecting{quant_str} {n_spectra_to_collect} spectra...")
    
            # Collect the next batch of spectra (and quantify if requested)
            tot_n_spectra, is_acquisition_successful = self._collect_spectra(
                n_spectra_to_collect,
                n_tot_sp_collected=tot_n_spectra,
                quantify=(quantify and not is_exp_std_measurement)
            )
    
            # Save temporary data file to avoid data loss
            if self.is_acquisition:
                self._save_collected_data(None, None, backup_previous_data=False, include_spectral_data=True)
                
            if self.verbose:
                print_single_separator()
                print(f"{tot_n_spectra}/{tot_spectra_to_collect} spectra collected and saved.")
                
            # Collect additional spectra in next iteration
            n_spectra_to_collect = min(
                max_n_sp_per_iter,
                tot_spectra_to_collect - tot_n_spectra,
                self.min_n_spectra
            )
            
            if quantify and tot_n_spectra > 0:
                if is_exp_std_measurement:
                    # Fit spectra and check if target number of good spectra has been collected
                    is_analysis_successful, is_converged = self._evaluate_exp_std_fit(tot_n_spectra)
                else:
                    # Perform clustering analysis and check for convergence
                    is_analysis_successful, is_converged = self._evaluate_clustering_convergence(tot_n_spectra, n_spectra_to_collect)
                    
                if is_converged:
                    break
                
                
            # Stop if no more particles are available on the sample
            if not is_acquisition_successful:
                if self.verbose:
                    print("Acquisition interrupted.")
                    if self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE:
                        print(f'Not enough particles were found on the sample to collect all {tot_spectra_to_collect} spectra.')
                    elif self.sample_cfg.type == cnst.S_BULK_SAMPLE_TYPE:
                        print(f'The specified spectrum spacing did not allow to collect all {tot_spectra_to_collect} spectra.\n'
                              "Change spacing in bulk_meas_cfg to collect more spectra.")
                break
    
        print_double_separator()
        print('Sample ID: %s' % self.sample_cfg.ID)
        par_str = f' over {self.particle_cntr} particles' if self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE else ''
        print(f'{tot_n_spectra} spectra were collected{par_str}.')
        process_time = (time.time() - self.start_process_time) / 60
        print(f'Total compositional analysis time: {process_time:.1f} min')
        print_single_separator()
    
        if quantify and not is_exp_std_measurement:
            if is_analysis_successful:
                if is_converged:
                    print('Clustering converged to small errors. All phases identified with confidence higher than 0.8.')
                else:
                    print('Phases could not be identified with confidence higher than 0.8.')
    
                self.print_results()
    
            elif not is_acquisition_successful:
                print('This did not allow to determine which phases are present in the sample.')
            else:
                print(f'Phases could not be identified with the allowed maximum of {self.max_n_spectra} collected spectra.')
                is_analysis_successful = False
                is_converged = False
        else:
            is_analysis_successful = is_acquisition_successful
    
        return is_analysis_successful, is_converged
    
    
    def _evaluate_clustering_convergence(
        self,
        tot_n_spectra: int,
        n_spectra_to_collect: int
    ) -> Tuple[bool, bool]:
        """
        Evaluate whether compositional clustering analysis has converged.
    
        This method checks the results of the clustering analysis after a given number of spectra 
        have been collected. It determines whether the analysis has converged based on the 
        clustering standard deviation and minimum confidence, and whether additional spectra 
        should be collected.
    
        Parameters
        ----------
        tot_n_spectra : int
            Total number of spectra collected so far.
        n_spectra_to_collect : int
            Total number of spectra to be collected.
    
        Returns
        -------
        Tuple[bool, bool]
            A tuple containing:
            - is_analysis_successful (bool): Whether the clustering analysis ran successfully.
            - is_converged (bool): Whether the compositional analysis has converged.
    
        Raises
        ------
        RuntimeError
            If `analyse_data` returns unexpected results.
        """
    
        if self.verbose:
            print_double_separator()
            print(f"Analysing phases after collection of {tot_n_spectra} spectra...")
    
        try:
            is_analysis_successful, max_cl_rmsdist, min_conf = self.analyse_data(
                self.clustering_cfg.max_analytical_error_percent,
                k=self.clustering_cfg.k
            )
        except Exception as e:
            raise RuntimeError(f"Error during clustering analysis: {e}") from e
    
        # Default value in case convergence check is skipped
        is_converged = False
    
        if is_analysis_successful:
            if self.verbose:
                print_double_separator()
                print("Clustering analysis performed")
    
            # Check whether phase identification converged
            try:
                is_converged = self._is_comp_analysis_converged(max_cl_rmsdist, min_conf)
            except Exception as e:
                raise RuntimeError(f"Error while checking convergence: {e}") from e
    
            if tot_n_spectra >= self.min_n_spectra:
                if is_converged:
                    return is_analysis_successful, is_converged
                elif self.verbose and n_spectra_to_collect > 0:
                    print("Compositional analysis did not converge, more spectra will be collected.")
            elif tot_n_spectra >= self.max_n_spectra:
                print(f"Maximum allowed number of {self.max_n_spectra} was acquired.")
            else:
                if self.verbose:
                    print(f"Collecting additional spectra to reach minimum number of {self.min_n_spectra}.")
    
        elif self.verbose:
            print("Clustering analysis unsuccessful.")
            if n_spectra_to_collect > 0:
                print(", more spectra will be collected.")
    
        return is_analysis_successful, is_converged
    
    
    def _is_comp_analysis_converged(
        self,
        stdev: float,
        min_conf: Optional[float]
    ) -> bool:
        """
        Determine if the clustering analysis has converged based on cluster statistics.
        Used when collecting and quantifying spectra in real time.
    
        Convergence criteria:
          - If no reference phases are present or assigned (min_conf is None), require cluster standard deviation < 2.5%.
          - If reference phases are assigned, require minimum confidence > 0.8 and cluster standard deviation < 3%.
    
        Parameters
        ----------
        stdev : float
            Maximum standard deviation among clusters (fractional units, e.g., 0.025 for 2.5%).
        min_conf : float or None
            Minimum confidence among all clusters assigned to reference phases. If None, no references are assigned.
    
        Returns
        -------
        is_converged : bool
            True if convergence criteria are met, False otherwise.
    
        Notes
        -----
        - The thresholds are empirically determined for robust phase identification.
        """
        if min_conf is None:
            # No reference phases present or assigned; require tighter cluster homogeneity
            is_converged = stdev < 0.025
        else:
            # Require high confidence and allow slightly larger within-cluster spread
            is_converged = (min_conf > 0.8) and (stdev < 0.03)
    
        return is_converged
    

    def run_quantification(self) -> None:
        """
        Perform quantification of all collected spectra and save the results.
    
        This method quantifies the spectra using self._fit_and_quantify_spectra(),
        then saves the quantification results to file using self._save_collected_data().
    
        Notes
        -----
        - The arguments (None, None) to _save_collected_data indicate that all spectra are to be saved.
        - Assumes that spectra have been correctly saved in self.quant_results
        """
        self._initialise_std_dict() # Initialise dictionary of standards to (optionally) pass onto XSp_Quantifier. Only used with known powder mixtures
        self._fit_and_quantify_spectra()
        self._save_collected_data(None, None, backup_previous_data=True, include_spectral_data=True)
        
    
    def run_exp_std_collection(
        self,
        fit_during_collection: bool,
        update_std_library: bool
    ) -> None:
        """
        Collect, fit, and optionally update the library of experimental standards.
    
        This method automates the acquisition and fitting of spectra from experimental 
        standards, ensuring that all required elemental fractions are defined before 
        proceeding.
    
        Parameters
        ----------
        fit_during_collection : bool
            If True, spectra will be fitted in real-time during collection.
            If False, fitting must be performed after collection.
        update_std_library : bool
            If True, the experimental standard library will be updated with the 
            newly fitted PB ratios.
    
        Raises
        ------
        ValueError
            If `self.exp_stds_cfg.is_exp_std_measurement` is not set to True.
        KeyError
            If any element in `self.sample_cfg.elements` is missing from 
            `self.exp_stds_cfg.w_frs`.
        """
    
        if not self.exp_stds_cfg.is_exp_std_measurement:
            raise ValueError(
                "Experimental standard collection mode is not active. "
                "Set `self.exp_stds_cfg.is_exp_std_measurement = True` before running."
            )
        
        # Ensure all elemental fractions are defined in the experimental standard configuration
        missing = [el for el in self.sample_cfg.elements if el not in self.exp_stds_cfg.w_frs]
        if missing:
            raise KeyError(
                f"The following elements are missing from `exp_stds_cfg.formula`: "
                f"{', '.join(str(m) for m in missing)}. "
                f"Ensure the formula contains all elements defined in `self.sample_cfg.elements`."
            )
        
        if self.verbose:
            print_double_separator()
            print(f"Experimental standard acquisition of {self.sample_cfg.ID}")
        
        # Run collection and quantification (fitting optionally performed during collection)
        self._th_peak_energies = {} # Initialise
        self.run_collection_and_quantification(quantify=fit_during_collection)
        
        # Fit standards and save results
        std_ref_lines, results_df, Z_sample = self._fit_stds_and_save_results(backup_previous_data=False)
        
        # Optionally update the standards library with the new results
        if update_std_library and std_ref_lines is not None and std_ref_lines is not {}: 
            self._update_standard_library(std_ref_lines, results_df, Z_sample)
        
    #%% Save Plots
    def _save_plots(
        self,
        kmeans: 'KMeans',
        compositions_df: 'pd.DataFrame',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list
    ) -> None:
        """
        Generate and save clustering and silhouette plots for the clustering analysis.
    
        This method:
          - Saves a silhouette plot if more than one cluster is present.
          - Determines which elements to include in the clustering plot (max 3 for 3D).
          - Excludes elements as specified in plot configuration.
          - Warns if only one element is available for plotting.
          - Saves a clustering plot (2D or 3D) using either a custom or default plotting function.
    
        Parameters
        ----------
        kmeans : KMeans
            Fitted KMeans clustering object.
        compositions_df : pd.DataFrame
            DataFrame of sample compositions used for clustering.
        centroids : np.ndarray
            Array of cluster centroids (shape: n_clusters, n_features).
        labels : np.ndarray
            Cluster labels for each sample.
        els_std_dev_per_cluster : list
            List of standard deviations for each element in each cluster.
        unused_compositions_list : list
            List of compositions excluded from clustering.
    
        Returns
        -------
        None
    
        Notes
        -----
        - Only up to 3 elements can be plotted; others are excluded.
        - If only one element is left after exclusions, no plot is generated.
        - Uses either a custom or default plotting function as configured.
        """
        # Silhouette plot (only if more than one cluster)
        if len(centroids) > 1:
            EMXSp_Composition_Analyzer._save_silhouette_plot(
                kmeans, compositions_df, self.analysis_dir, show_plot=self.plot_cfg.show_plots
            )
    
        # Determine which elements can be used for plotting
        can_plot_clustering = True
    
        # Elements for plotting (excluding those set for exclusion)
        els_for_plot = list(set(self.detectable_els_sample) - set(self.plot_cfg.els_excluded_clust_plot))
        els_excluded_clust_plot = list(set(self.all_els_sample) - set(els_for_plot))
        n_els = len(els_for_plot)
    
        if n_els == 1:
            # Cannot plot with only 1 detectable element
            can_plot_clustering = False
            print_single_separator()
            warnings.warn("Cannot generate clustering plot with a single element.", UserWarning)
            if len(self.detectable_els_sample) > 1:
                print('Too many elements were excluded from the clustering plot via the use of "els_excluded_clust_plot".')
                print(f'Consider removing one or more among the list: {self.plot_cfg.els_excluded_clust_plot}')
        elif n_els > 3:
            # Only 3 elements can be plotted at once (for 3D)
            els_excluded_clust_plot += els_for_plot[3:]
            els_for_plot = els_for_plot[:3]
    
        # Determine indices to remove for excluded elements
        indices_to_remove = [self.all_els_sample.index(el) for el in els_excluded_clust_plot]
        # Update values to exclude the selected elements
        els_for_plot = [el for i, el in enumerate(self.all_els_sample) if i not in indices_to_remove]
        centroids = np.array([
            [coord for i, coord in enumerate(row) if i not in indices_to_remove]
            for row in centroids
        ])
        els_std_dev_per_cluster = [
            [stddev for i, stddev in enumerate(row) if i not in indices_to_remove]
            for row in els_std_dev_per_cluster
        ]
        unused_compositions_list = [
            [fr for i, fr in enumerate(row) if i not in indices_to_remove]
            for row in unused_compositions_list
        ]
    
        # Generate and save the clustering plot if possible
        if can_plot_clustering:
            # List of lists, where each list is populated with the atomic fractions of one element in all data points
            els_comps_list = compositions_df[els_for_plot].to_numpy().T
    
            # Use custom or default plotting function
            if self.plot_cfg.use_custom_plots:
                custom_fnctns._save_clustering_plot_custom_3D(
                    els_for_plot, els_comps_list, centroids, labels,
                    els_std_dev_per_cluster, unused_compositions_list,
                    self.clustering_cfg.features, self.ref_phases_df,
                    self.ref_formulae, self.plot_cfg.show_plots, self.sample_cfg.ID
                )
            else:
                self._save_clustering_plot(
                    els_for_plot, els_comps_list, centroids, labels,
                    els_std_dev_per_cluster, unused_compositions_list
                )
        elif self.verbose:
            print('Clusters were not plotted because only one detectable element was present in the sample.')
            print(f"Elements {calibs.undetectable_els} cannot be detected at the employed instrument.")
            
            
    def _save_clustering_plot(
        self,
        elements: List[str],
        els_comps_list: 'np.ndarray',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list
    ) -> None:
        """
        Generate and save a 2D or 3D clustering plot with centroids, standard deviation ellipses/ellipsoids,
        unused compositions, and reference phases.
    
        Parameters
        ----------
        elements : list of str
            List of element symbols to plot (max 3).
        els_comps_list : np.ndarray
            Array of shape (n_elements, n_samples) with elemental fractions for each sample (used for clustering).
        centroids : np.ndarray
            Array of cluster centroids (shape: n_clusters, n_elements).
        labels : np.ndarray
            Cluster labels for each sample.
        els_std_dev_per_cluster : list
            List of standard deviations for each element in each cluster.
        unused_compositions_list : list
            List of compositions excluded from clustering.
        
        Returns
        -------
        None
    
        Notes
        -----
        - The plot is saved as 'Clustering_plot.png' in the analysis directory.
        - Uses matplotlib for plotting (2D or 3D based on the number of elements).
        - Reference phases and centroids are annotated; standard deviation is shown as ellipses (2D) or ellipsoids (3D).
        """
        # Set font parameters
        plt.rcParams['font.family'] = 'Arial'
        fontsize = 14
        labelpad = 12
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.titlesize'] = fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize
    
        # Define axis label suffix
        axis_label_add = ' (w%)' if self.clustering_cfg.features == cnst.W_FR_CL_FEAT else ' (at%)'
        ticks = np.arange(0, 1, 0.1)
        ticks_labels = [f"{x*100:.0f}" for x in ticks]
    
        # Create figure and axes
        fig = plt.figure(figsize=(6, 6))
        if len(elements) == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.set_zlabel(elements[2] + axis_label_add, labelpad=labelpad * 0.95)
            ax.set_zlim(0, 1)
            ax.set_zticks(ticks)
            ax.set_zticklabels(ticks_labels)
        else:
            ax = fig.add_subplot(111)
    
        # Plot clustered datapoints
        ax.scatter(*els_comps_list, c=labels, cmap='viridis', marker='o')
    
        # Plot centroids
        ax.scatter(*centroids.T, c='red', marker='x', s=100, label='Centroids')
    
        # Plot standard deviation ellipses or ellipsoids
        first_ellipse = True
        for centroid, stdevs in zip(centroids, els_std_dev_per_cluster):
            if ~np.any(np.isnan(stdevs)):
                if len(elements) == 3:  # 3D plot
                    x_c, y_c, z_c = centroid
                    rx, ry, rz = stdevs
    
                    # Create the ellipsoid
                    u = np.linspace(0, 2 * np.pi, 100)
                    v = np.linspace(0, np.pi, 100)
                    x = x_c + rx * np.outer(np.cos(u), np.sin(v))
                    y = y_c + ry * np.outer(np.sin(u), np.sin(v))
                    z = z_c + rz * np.outer(np.ones_like(u), np.cos(v))
    
                    # Plot the surface with transparency
                    ax.plot_surface(x, y, z, color='red', alpha=0.1, edgecolor='none')
    
                    if first_ellipse:
                        first_ellipse = False
                        ax.plot([], [], [], color='red', alpha=0.1, label='Stddev')
                else:  # 2D plot
                    x_c, y_c = centroid
                    rx, ry = stdevs
    
                    # Plot the ellipse with transparency
                    ellipse = patches.Ellipse((x_c, y_c), rx, ry, edgecolor='red', facecolor='red', linestyle='--', alpha=0.2)
                    if first_ellipse:
                        ellipse.set_label('Stddev')
                        first_ellipse = False
                    ax.add_patch(ellipse)
    
        # Plot unused compositions (discarded from clustering)
        if unused_compositions_list and self.plot_cfg.show_unused_comps_clust:
            ax.scatter(*np.array(unused_compositions_list).T, c='grey', marker='^', label='Discarded comps.')
    
        # Plot reference phases
        if self.ref_formulae is not None:
            first_ref = True
            ref_phases_df = self.ref_phases_df[elements]
            for index, row in ref_phases_df.iterrows():
                label = 'Candidate phases' if first_ref else None
                ax.scatter(*row.values, c='blue', marker='*', s=100, label=label)
                ref_label = to_latex_formula(self.ref_formulae[index])
                ax.text(*row.values, ref_label, color='black', fontsize=fontsize, ha='left', va='bottom')
                first_ref = False
    
        # Annotate centroids with their cluster labels
        for i, centroid in enumerate(centroids):
            ax.text(*centroid, str(i), color='black', fontsize=fontsize, ha='right', va='bottom')
    
        # Set axis labels and limits
        ax.set_xlabel(elements[0] + axis_label_add, labelpad=labelpad)
        ax.set_ylabel(elements[1] + axis_label_add, labelpad=labelpad)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks_labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks_labels)
        ax.set_title(f'{self.clustering_cfg.method} clustering {self.sample_cfg.ID}')

        ax.legend(fontsize = fontsize)
        
        # plt.ion()
        # plt.show()
        # plt.pause(0.001)
        fig.savefig(os.path.join(self.analysis_dir, cnst.CLUSTERING_PLOT_FILENAME + cnst.CLUSTERING_PLOT_FILEEXT))
        # plt.close(fig)

    
    def _save_violin_plot_powder_mixture(
        self,
        W_mol_frs,
        ref_names
    ) -> None:
        """
        Generate and save a violin plot visualizing the distribution of precursor molar fractions in a binary powder mixture.
    
        The plot displays:
          - The kernel density estimate (KDE) of the measured molar fractions
          - Individual measured values
          - The mean and standard deviation of the distribution
    
        Note:
            Prior to plotting, the precursor molar fractions are normalized so their sum equals 1.
            As a result, the standard deviation is identical for both precursors in the mixture.
    
        Parameters
        ----------
        W_mol_frs : list(float)
            Measured molar fractions of the precursors for the current cluster, represented as a binary mixture of two powders.
        ref_names : list(str)
            Chemical formulas (or names) of the two parent phases forming the mixture.
        """
    
        # --- Plot styling ---
        plt.rcParams['font.family'] = 'Arial'
        fontsize = 17
        labelpad = 0
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.titlesize'] = fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize
        purple_cmap = cm.get_cmap('Purples')
        yellow_cmap = cm.get_cmap('autumn')
    
        # Extract coordinates from W
        y_vals = W_mol_frs[:, 0]
    
        fig, ax_left = plt.subplots(figsize=(4, 4))
        
        mean = np.mean(y_vals)
        std = np.std(y_vals)
        
        # Violin plot (default zorder is 1)
        ax_left = sns.violinplot(data=y_vals, inner=None, color=purple_cmap(0.3),
                                 linewidth=1.5, density_norm='area', width=1, zorder=1)
        
        # Swarmplot (zorder 2)
        sns.swarmplot(data=y_vals, color=purple_cmap(0.8),
                      edgecolor=purple_cmap(1.0), linewidth=2, size=5, label='data', zorder=2)
        
        # Error bars (zorder 3 and 4)
        ax_left.errorbar(0, mean, yerr=std / 2, fmt='none', color=yellow_cmap(0.9),
                         label='Mean ±1 Std Dev', capsize=5, elinewidth=1,
                         zorder=4, markerfacecolor=yellow_cmap(0.9),
                         markeredgecolor='black', markeredgewidth=1,
                         marker='o', linestyle='none')
        ax_left.errorbar(0, mean, yerr=std / 2, fmt='none', color='none',
                         label='_nolegend_', capsize=6, elinewidth=2,
                         zorder=3, markerfacecolor='none',
                         markeredgecolor='black', markeredgewidth=2,
                         marker='o', linestyle='none', ecolor='black')
        
        # Mean point (highest zorder, plotted last)
        ax_left.scatter(0, mean, color=yellow_cmap(0.9), marker='o', s=50,
                        edgecolors='k', linewidths=1, label='Mean', zorder=10)
    
        ax_left = plt.gca()
        ax_left.set_xticks([])
        ax_left.set_yticks([0, 1])  # Show ticks at 0 and 1 on left y-axis
        ax_left.set_frame_on(True)
        for spine in ax_left.spines.values():
            spine.set_color('black')
            spine.set_linewidth(0.5)
        plt.grid(False)
    
        plt.xlim(-0.5, 0.5)
        ylim_bottom = 0
        ylim_top = 1
        ax_left.set_ylim(ylim_bottom, ylim_top)
    
        # Left y-axis label (0→1)
        left_formula = to_latex_formula(ref_names[0], include_dollar_signs=False)
        ax_left.set_ylabel(rf"$x_{{\mathrm{{{left_formula}}}}}$", labelpad=labelpad)
        
        # Right y-axis (inverted 1→0)
        ax_right = ax_left.twinx()
        ax_right.set_ylim(ylim_top, ylim_bottom)
        ax_right.set_yticks([1, 0])  # Inverted ticks on right y-axis
        right_formula = to_latex_formula(ref_names[1], include_dollar_signs=False)
        ax_right.set_ylabel(rf"$x_{{\mathrm{{{right_formula}}}}}$", labelpad=labelpad)
    
        # Add standard deviation inside the plot
        ax_left.text(0.03, 0.03, rf"$\sigma_x = {std*100:.1f}$%", fontsize=fontsize,
                     color='black', ha='left', va='bottom', transform=ax_left.transAxes)
    
        ax_left.set_title(f'Violin plot {self.sample_cfg.ID}')
    
        # Save figure
        fig.savefig(
            os.path.join(self.analysis_dir,
                         cnst.POWDER_MIXTURE_PLOT_FILENAME + f"_{ref_names[0]}_{ref_names[1]}" + cnst.CLUSTERING_PLOT_FILEEXT),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0
        )

            
    @staticmethod
    def _save_silhouette_plot(
        model: 'KMeans',
        compositions_df: 'pd.DataFrame',
        results_dir: str,
        show_plot: bool
    ) -> None:
        """
        Generate and save a silhouette plot for the clustering results.
    
        Parameters
        ----------
        model : KMeans
            Fitted clustering model.
        compositions_df : pd.DataFrame
            DataFrame of sample compositions used for clustering.
        results_dir : str
            Directory where the plot will be saved.
        show_plot : bool
            If True, the plot will be displayed interactively.
    
        Returns
        -------
        None
    
        Notes
        -----
        - Uses Yellowbrick's SilhouetteVisualizer for plotting.
        - Suppresses harmless sklearn warnings during visualization.
        - The plot is saved as 'Silhouette_plot.png' in the results directory.
        """
        plt.figure(figsize=(10, 8))
        sil_visualizer = SilhouetteVisualizer(model, colors='yellowbrick')
        with warnings.catch_warnings():
            # Suppress harmless sklearn warnings
            warnings.simplefilter("ignore", UserWarning)
            sil_visualizer.fit(compositions_df)  # Fit the data to the visualizer
    
        plt.ylabel('Cluster label')
        plt.xlabel('Silhouette coefficient values')
        plt.legend(loc='upper right', frameon=True)
    
        if show_plot:
            plt.ion()
            sil_visualizer.show()
            plt.pause(0.001)
    
        fig = sil_visualizer.fig
        fig.savefig(os.path.join(results_dir, 'Silhouette_plot.png'))
    
        # Close the figure if not displaying
        if not show_plot:
            plt.close(fig)
    
    
    #%% Save Data
    # =============================================================================
    def _save_result_and_stats(
        self,
        centroids: 'np.ndarray',
        els_std_dev_per_cluster: list,
        centroids_other_fr: 'np.ndarray',
        els_std_dev_per_cluster_other_fr: list,
        n_points_per_cluster: list,
        wcss_per_cluster: list,
        rms_dist_cluster: list,
        rms_dist_cluster_other_fr: list,
        refs_assigned_df: 'pd.DataFrame',
        wcss: float,
        sil_score: float,
        tot_n_points: int,
        max_analytical_error: float,
        clusters_assigned_mixtures: list
    ) -> None:
        """
        Save and store clustering results and statistics, including centroids, standard deviations, reference assignments, mixture assignments, and summary statistics.
    
        This method:
          - Constructs a DataFrame of cluster statistics and assignments.
          - Adds reference phase and mixture assignments if available.
          - Saves the DataFrame to CSV and stores it as an attribute.
          - Saves general clustering information to a JSON file and stores it as an attribute.
    
        Parameters
        ----------
        centroids : np.ndarray
            Array of cluster centroids (shape: n_clusters, n_features).
        els_std_dev_per_cluster : list
            Standard deviations for each elemental fraction in each cluster (same shape as centroids).
        centroids_other_fr : np.ndarray
            Centroids in the alternate fraction representation.
        els_std_dev_per_cluster_other_fr : list
            Standard deviations of elemental fractions in the alternate fraction representation.
        n_points_per_cluster : list
            Number of points in each cluster.
        wcss_per_cluster : list
            Within-cluster sum of squares for each cluster.
        rms_dist_cluster : list
            Standard deviation of distances to centroid for each cluster.
        rms_dist_cluster_other_fr : list
            Standard deviation of distances in the alternate fraction representation.
        refs_assigned_df : pd.DataFrame
            DataFrame of reference assignments for each cluster.
        wcss : float
            Total within-cluster sum of squares.
        sil_score : float
            Silhouette score for the clustering.
        tot_n_points : int
            Total number of spectra considered.
        max_analytical_error : float
            Maximum allowed analytical error for spectra used in clustering.
        clusters_assigned_mixtures : list
            List of mixture assignments for each cluster.
    
        Returns
        -------
        None
    
        Notes
        -----
        - The cluster DataFrame is saved as '<cnst.CLUSTERS_FILENAME>.csv' in the analysis directory.
        - General clustering info is saved as '<cnst.CLUSTERING_INFO_FILENAME>.json'.
        - Both are also stored as attributes for later use.
    
        Raises
        ------
        OSError
            If the analysis directory cannot be created or files cannot be written.
    
        Suggestions
        -----------
        - Consider using more explicit type hints for lists, e.g., List[float] or List[int], and for DataFrames, use pd.DataFrame directly.
        - If you expect large data, consider saving DataFrames in a binary format (e.g., Parquet) for efficiency.
        """
    
        # Prepare cluster statistics as dictionaries for DataFrame construction
        els_fr = np.transpose(centroids)
        els_other_fr = np.transpose(centroids_other_fr)
        els_stdevs = np.transpose(els_std_dev_per_cluster)
        els_stdevs_other_fr = np.transpose(els_std_dev_per_cluster_other_fr)
    
        # Select keys for fraction and standard deviation columns based on configuration
        if self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
            fr_key = cnst.AT_FR_DF_KEY
            other_fr_key = cnst.W_FR_DF_KEY
        elif self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
            fr_key = cnst.W_FR_DF_KEY
            other_fr_key = cnst.AT_FR_DF_KEY
        else:
            # Suggestion: handle unexpected feature settings
            raise ValueError(f"Unknown clustering feature: {self.clustering_cfg.features}")
    
        stddev_key = cnst.STDEV_DF_KEY + fr_key
        other_stddev_key = cnst.STDEV_DF_KEY + other_fr_key
    
        # Prepare dictionaries for DataFrame columns
        els_fr_dict = {el + fr_key: np.round(el_comps * 100, 2) for el, el_comps in zip(self.all_els_sample, els_fr)}
        els_fr_stdevs_dict = {el + stddev_key: np.round(el_stddev * 100, 2) for el, el_stddev in zip(self.all_els_sample, els_stdevs)}
        els_other_fr_dict = {el + other_fr_key: np.round(el_comps * 100, 2) for el, el_comps in zip(self.all_els_sample, els_other_fr)}
        els_other_fr_stdevs_dict = {el + other_stddev_key: np.round(el_stddev * 100, 2) for el, el_stddev in zip(self.all_els_sample, els_stdevs_other_fr)}
    
        # Compose main cluster DataFrame
        clusters_dict = {
            cnst.N_PTS_DF_KEY: n_points_per_cluster,
            **els_fr_dict,
            **els_fr_stdevs_dict,
            **els_other_fr_dict,
            **els_other_fr_stdevs_dict,
            cnst.RMS_DIST_DF_KEY + fr_key: (np.array(rms_dist_cluster) * 100).round(2),
            cnst.RMS_DIST_DF_KEY + other_fr_key: (np.array(rms_dist_cluster_other_fr) * 100).round(2),
            cnst.WCSS_DF_KEY + fr_key: (np.array(wcss_per_cluster) * 10000).round(2)
        }
        clusters_df = pd.DataFrame(clusters_dict)
    
        # Add reference assignments if available
        if self.ref_formulae:
            clusters_df = pd.concat([clusters_df.reset_index(drop=True), refs_assigned_df.reset_index(drop=True)], axis=1)
    
        # Add mixture assignments if any
        mixtures_df = self._build_mixtures_df(clusters_assigned_mixtures)
        clusters_df = pd.concat([clusters_df.reset_index(drop=True), mixtures_df.reset_index(drop=True)], axis=1)
    
        # Ensure the analysis directory exists
        try:
            os.makedirs(self.analysis_dir, exist_ok=True)
        except Exception as e:
            raise OSError(f"Could not create analysis directory '{self.analysis_dir}': {e}")
    
        # Save and store DataFrame
        clusters_csv_path = os.path.join(self.analysis_dir, cnst.CLUSTERS_FILENAME + '.csv')
        try:
            clusters_df.to_csv(clusters_csv_path, index=True, header=True)
        except Exception as e:
            raise OSError(f"Could not write clusters DataFrame to '{clusters_csv_path}': {e}")
    
        self.clusters_df = clusters_df
    
        # Save general clustering info and store for printing
        now = datetime.now()
    
        # Gather configuration dataclasses as dictionaries
        cfg_dataclasses = {
            cnst.QUANTIFICATION_CFG_KEY: asdict(self.quant_cfg),
            cnst.CLUSTERING_CFG_KEY: asdict(self.clustering_cfg),
            cnst.PLOT_CFG_KEY: asdict(self.plot_cfg),
        }
    
        clustering_info = {
            cnst.DATETIME_KEY: now.strftime("%Y-%m-%d %H:%M:%S"),
            cnst.N_SP_ACQUIRED_KEY: tot_n_points,
            cnst.N_SP_USED_KEY: sum(n_points_per_cluster),
            cnst.N_CLUST_KEY: len(centroids),
            cnst.WCSS_KEY: wcss,
            cnst.SIL_SCORE_KEY: sil_score,
            **cfg_dataclasses,
        }
        clustering_json_path = os.path.join(self.analysis_dir, cnst.CLUSTERING_INFO_FILENAME + '.json')
        try:
            with open(clustering_json_path, 'w', encoding='utf-8') as file:
                json.dump(clustering_info, file, indent=2, ensure_ascii=False)
        except Exception as e:
            raise OSError(f"Could not write clustering info to '{clustering_json_path}': {e}")
    
        self.clustering_info = clustering_info
        
        
    def _save_collected_data(
        self,
        labels: List,
        df_indices: List,
        backup_previous_data: Optional[bool] = True,
        include_spectral_data: Optional[bool] = True,
    ) -> None:
        """
        Save the collected spectra, their quantification (optionally), and their cluster assignments.
    
        This method builds a DataFrame with quantification and clustering information for each spectrum,
        along with spectral data if requested. It ensures unique output files and backs up existing data.
        Spectra with insufficient counts are handled gracefully.
    
        Parameters
        ----------
        labels : List
            List of cluster labels assigned to each spectrum (by index in df_indices).
        df_indices : List
            List of indices mapping labels to DataFrame rows.
        backup_previous_data : bool, optional
            Backs up previous data file if present (Default = True).
        include_spectral_data : bool, optional
            If True, includes raw spectral and background data in the output (default: True).
            
        Returns
        -------
        data_df: pd.Dataframe
            Dataframe object containing the saved data. Used only when measuring experimental standards
            None if no spectrum to save was acquired
    
        Notes
        -----
        - If a file with the intended name exists and contains quantification data, a counter is appended to the filename.
        - If include_spectral_data is False, only compositions are saved (to the analysis directory).
        - Existing main data files are backed up before being overwritten.
        - Uses make_unique_dir for unique directory creation if necessary.
    
        Raises
        ------
        OSError
            If the output directory cannot be created or files cannot be written.
        """
    
        # Get list with quantification data
        quant_result_list = self.spectra_quant
        is_standards_measurements = self.exp_stds_cfg.is_exp_std_measurement
    
        # Determine the number of spectra to process
        if include_spectral_data:
            n_spectra = len(self.spectral_data[cnst.SPECTRUM_DF_KEY])
        else:
            n_spectra = len(quant_result_list)
    
        if n_spectra == 0:
            return None
        
        data_list = []
        for i in range(n_spectra):
            # Check if spectrum has been quantified
            if i < len(quant_result_list) and quant_result_list[i] is not None:
                quant_dict = quant_result_list[i]
                
                if is_standards_measurements:
                    exp_std_comp_d = self.exp_stds_cfg.w_frs
                    std_els = list(exp_std_comp_d.keys())
                    std_w_frs = list(exp_std_comp_d.values())
                    std_at_frs = weight_to_atomic_fr(std_w_frs,std_els,verbose= False)
                    atomic_comp = {el + cnst.AT_FR_DF_KEY: round(fr *100, 2) for el, fr in zip(std_els, std_at_frs)}
                    weight_comp = {el + cnst.W_FR_DF_KEY: round(fr *100, 2) for el, fr in exp_std_comp_d.items()}
                    meas_data = {**atomic_comp, **weight_comp, **quant_dict}
                else:
                    # Unpack spectral quantification results and convert from elemental fraction to % for readability
                    atomic_comp = {el + cnst.AT_FR_DF_KEY: round(fr *100, 2) for el, fr in quant_dict[cnst.COMP_AT_FR_KEY].items()}
                    weight_comp = {el + cnst.W_FR_DF_KEY: round(fr *100, 2) for el, fr in quant_dict[cnst.COMP_W_FR_KEY].items()}
                    analytical_er = {cnst.AN_ER_DF_KEY: round(quant_dict[cnst.AN_ER_KEY] *100, 2)}
        
                    # Extract cluster label if assigned
                    try:
                        label_index = df_indices.index(i)
                        cluster_n = labels[label_index]
                    except ValueError:
                        cluster_n = pd.NA
                    except AttributeError:
                        cluster_n = pd.NA
                    
                    # Compose row of data to be saved
                    meas_data = {
                        cnst.CL_ID_DF_KEY: cluster_n,
                        **atomic_comp,
                        **analytical_er,
                        **weight_comp,
                        cnst.R_SQ_KEY: quant_dict[cnst.R_SQ_KEY],
                        cnst.REDCHI_SQ_KEY: quant_dict[cnst.REDCHI_SQ_KEY]
                    } 
                    
                # Compose row of data to be saved
                data_row = {
                    **self.sp_coords[i],
                    **meas_data
                }
            else:
                # Counts in this spectrum were too low
                data_row = self.sp_coords[i]
            
            # Add comment and quantification flag columns, if available
            try:
                data_row[cnst.COMMENTS_DF_KEY] = self.spectral_data[cnst.COMMENTS_DF_KEY][i]
                data_row[cnst.QUANT_FLAG_DF_KEY] = self.spectral_data[cnst.QUANT_FLAG_DF_KEY][i]
            except Exception:
                pass
    
            # Add spectral data
            if include_spectral_data:
                try:  # If background present
                    background_entry = '[' + ','.join(map(str, self.spectral_data[cnst.BACKGROUND_DF_KEY][i])) + ']'
                except Exception:
                    background_entry = None  # Use None so pandas will recognize as missing
    
                real_time = self.spectral_data[cnst.REAL_TIME_DF_KEY][i]
                live_time = self.spectral_data[cnst.LIVE_TIME_DF_KEY][i]
    
                if real_time is not None:
                    real_time = round(real_time, 2)
                if live_time is not None:
                    live_time = round(live_time, 2)
    
                # Format strings to avoid truncation when saving dataframe into csv
                data_row = {
                    **data_row,
                    cnst.REAL_TIME_DF_KEY: real_time,
                    cnst.LIVE_TIME_DF_KEY: live_time,
                    cnst.SPECTRUM_DF_KEY: '[' + ','.join(map(str, self.spectral_data[cnst.SPECTRUM_DF_KEY][i])) + ']',
                    cnst.BACKGROUND_DF_KEY: background_entry
                }
    
            data_list.append(data_row)
    
        # Convert list of dictionaries to DataFrame
        data_df = pd.DataFrame(data_list)
    
        # Remove Cluster ID column if no value has been assigned
        if cnst.CL_ID_DF_KEY in data_df.columns:
            if data_df[cnst.CL_ID_DF_KEY].isna().all():
                data_df.pop(cnst.CL_ID_DF_KEY)
            else:
                # Convert to nullable integer Int64 dtype
                data_df[cnst.CL_ID_DF_KEY] = data_df[cnst.CL_ID_DF_KEY].astype('Int64')
    
        # Remove background column if it is entirely None or NaN
        if cnst.BACKGROUND_DF_KEY in data_df.columns:
            if data_df[cnst.BACKGROUND_DF_KEY].isna().all():
                data_df.pop(cnst.BACKGROUND_DF_KEY)
    
        # Reorder columns to ensure spectral data is at the end
        columns = data_df.columns.tolist()
        last_columns = [
            cnst.R_SQ_KEY, cnst.REDCHI_SQ_KEY,
            cnst.QUANT_FLAG_DF_KEY, cnst.COMMENTS_DF_KEY,
            cnst.REAL_TIME_DF_KEY, cnst.LIVE_TIME_DF_KEY,
            cnst.SPECTRUM_DF_KEY, cnst.BACKGROUND_DF_KEY
        ]
        remaining_columns = [col for col in columns if col not in last_columns]
        new_column_order = remaining_columns + [col for col in last_columns if col in columns]
        data_df = data_df[new_column_order]
    
        # Save dataframe
        if data_df is not None and include_spectral_data:
            # Data.csv must always be saved, as it's used for analysis
            base_name = f'{cnst.DATA_FILENAME}' if not is_standards_measurements else f'{cnst.STDS_MEAS_FILENAME}'
            extension = f'{cnst.DATA_FILEEXT}'
            data_path = os.path.join(self.sample_result_dir, base_name + extension)
            
            # Unless it's saving a temporary file, check if Data.csv already exists and backup old version
            if os.path.exists(data_path) and backup_previous_data:
                # Avoid backing up acquisition files that do not contain any quantification or fitting, simply replace them
                # Checks if current data file has Quant_flag in the headers
                if cnst.QUANT_FLAG_DF_KEY in pd.read_csv(data_path, nrows=0).columns:
                    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                    backup_path = make_unique_path(self.sample_result_dir, f'{base_name}_old_{timestamp}', extension)
                    try: 
                        shutil.copyfile(data_path, backup_path)
                    except Exception as e:
                        backup_successful = False
                        raise OSError(f"Could not backup previous {data_path} file. This file was not overwritten."
                                      f"Ensure you ovewrite it with a copy of {base_name}{extension} prior analysis:"
                                      f"{e}")
                    else:
                        backup_successful = True
                else:
                    backup_successful = True
            else:
                backup_successful = True
            
            # Save new version
            try:
                if backup_successful:
                    data_df.to_csv(data_path, index=False, header=True)
                    if isinstance(self.output_filename_suffix, str) and self.output_filename_suffix != '':
                        if backup_previous_data:
                            data_path_with_suffix = make_unique_path(self.sample_result_dir, base_name + self.output_filename_suffix, extension)
                        else:
                            data_path_with_suffix = os.path.join(self.sample_result_dir, base_name + self.output_filename_suffix + extension)
                        data_df.to_csv(data_path_with_suffix, index=False, header=True)
                else:
                    data_df.to_csv(data_path_with_suffix, index=False, header=True)
                    
            except Exception as e:
                raise OSError(f"Could not write '{data_path}': {e}")
    
        else:
            # Save only compositions in Compositions.csv.
            # Used when cluster analysis is performed and cluster_IDs are assigned to each spectrum
            comp_path = os.path.join(self.analysis_dir, cnst.COMPOSITIONS_FILENAME + '.csv')
            try:
                data_df.to_csv(comp_path, index=False, header=True)
            except Exception as e:
                raise OSError(f"Could not write compositions data to '{comp_path}': {e}")
            
        return data_df
    
    
    def _save_experimental_config(self, is_XSp_measurement) -> None:
        """
        Save all relevant configuration dataclasses and metadata related to the
        current spectrum collection/acquisition to a JSON file.
    
        The saved file includes:
            - Timestamp of saving
            - All configuration dataclasses
    
        This function is intended to be called after acquisition to ensure
        reproducibility and traceability of the experimental configuration.
    
        Raises
        ------
        OSError
            If the output directory cannot be created or file cannot be written.
        """
    
        now = datetime.now()
    
        # Gather configuration dataclasses as dictionaries
        cfg_dataclasses = {
            cnst.SAMPLE_CFG_KEY: asdict(self.sample_cfg),
            cnst.MICROSCOPE_CFG_KEY: asdict(self.microscope_cfg),
            cnst.MEASUREMENT_CFG_KEY: asdict(self.measurement_cfg),
            cnst.SAMPLESUBSTRATE_CFG_KEY: asdict(self.sample_substrate_cfg),
        }
        
        if is_XSp_measurement:
            cfg_dataclasses[cnst.QUANTIFICATION_CFG_KEY] = asdict(self.quant_cfg)
    
        # Include dataclass corresponding to sample type
        if self.sample_cfg.type == cnst.S_POWDER_SAMPLE_TYPE:
            cfg_dataclasses[cnst.POWDER_MEASUREMENT_CFG_KEY] = asdict(self.powder_meas_cfg)
        elif self.sample_cfg.type == cnst.S_BULK_SAMPLE_TYPE:
            cfg_dataclasses[cnst.BULK_MEASUREMENT_CFG_KEY] = asdict(self.bulk_meas_cfg)
        
        # Include dataclasses corresponding to measurement type
        if self.exp_stds_cfg.is_exp_std_measurement:
            cfg_dataclasses[cnst.EXP_STD_MEASUREMENT_CFG_KEY] = asdict(self.exp_stds_cfg)
        elif is_XSp_measurement:
            cfg_dataclasses[cnst.CLUSTERING_CFG_KEY] = asdict(self.clustering_cfg)
            cfg_dataclasses[cnst.PLOT_CFG_KEY] = asdict(self.plot_cfg)


    
        # Compose the metadata dictionary for saving
        spectrum_collection_info = {
            cnst.DATETIME_KEY: now.strftime("%Y-%m-%d %H:%M:%S"),
            **cfg_dataclasses
        }
    
        # Ensure the output directory exists before saving
        try:
            os.makedirs(self.sample_result_dir, exist_ok=True)
        except Exception as e:
            raise OSError(f"Could not create output directory '{self.sample_result_dir}': {e}")
    
        # Save data to JSON file in a human-readable format
        output_path = os.path.join(self.sample_result_dir, f"{cnst.ACQUISITION_INFO_FILENAME}.json")
        try:
            with open(output_path, 'w', encoding='utf-8') as file:
                json.dump(spectrum_collection_info, file, indent=2, ensure_ascii=False)
        except Exception as e:
            raise OSError(f"Could not write spectrum collection info to '{output_path}': {e}")
            
    #%% Print Results
    # =============================================================================            
    def _report_n_discarded_spectra(
        self,
        n_datapts: int,
        max_analytical_error: float
    ) -> None:
        """
        Print a summary report of discarded spectra and the reasons for their exclusion.
    
        Parameters
        ----------
        n_datapts : int
            Total number of spectra considered.
        max_analytical_error : float
            Maximum allowed analytical error (as a fraction, e.g., 0.05 for 5%). If None, this check is skipped.
    
        Notes
        -----
        - Prints detailed reasons for spectrum exclusion if self.verbose is True.
        - Advises user to check comments in the exported CSV files for more information.
        - Quantification flags indicate whether the quantification or the fit of each spectrum is likely to be affected by large errors:
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
        is_any_spectrum_discarded = self.n_sp_too_low_counts + self.n_sp_bad_quant + self.n_sp_too_high_an_err > 0
        if not (self.verbose and n_datapts > 0 and is_any_spectrum_discarded):
            return
    
        print_single_separator()
        print("Summary of Discarded Spectra")
        print("  → For details, see the 'Comments' column in Compositions.csv.")
    
        # Discarded due to low counts, insufficient background, or acquisition/fitting errors
        if self.n_sp_too_low_counts > 0:
            print(
                f"  • {self.n_sp_too_low_counts} spectra were discarded due to insufficient total counts, "
                f"background counts below the threshold ({self.quant_cfg.min_bckgrnd_cnts}), "
                "or errors during spectrum collection/fitting."
            )
    
        # Discarded due to quantification flags
        if self.n_sp_bad_quant > 0:
            print(
                f"  • {self.n_sp_bad_quant} spectra were discarded because they were flagged during quantification."
            )
    
        # Warning if more than half of the spectra were flagged
        if self.n_sp_bad_quant / n_datapts > 0.5:
            print_single_separator()
            print("  Warning: More than 50% of spectra were flagged during quantification!")
            print(
                "  Common causes for poor fits (quant_flag = 4) include missing elements in the fit.\n"
                "  Ensure that all elements present in your sample have been specified in the 'elements' argument "
                "  when initializing EMXSp_Composition_Analyzer."
            )
            
        # Discarded due to high analytical error
        if self.n_sp_too_high_an_err > 0:
            print(
                f"  • {self.n_sp_too_high_an_err} spectra were discarded because their analytical error "
                f"exceeded the maximum allowed value of {max_analytical_error*100:.1f}%."
            )
    
        print_single_separator()
            
        
    def print_results(self) -> None:
        """
        Print a summary of clustering results, including clustering configuration, metrics,
        and a table of identified phases with elemental fractions, standard deviations,
        and reference/mixture assignments if present.
    
        The method:
          - Prints main clustering configuration and metrics.
          - Prints a table of phases, each with number of points, elemental fractions (with stddev),
            cluster stddev, WCSS, reference assignments, and mixture information if available.
    
        Raises
        ------
        AttributeError
            If required attributes (clustering_info, clusters_df, etc.) are missing.
        KeyError
            If expected keys are missing from clustering_info or clusters_df.
        """
        # Print clustering info
        print_double_separator()
        print(f"Compositional analysis results for sample {self.sample_cfg.ID}:")
        print_single_separator()
        try:
            print('Clustering method: %s' % self.clustering_cfg.method)
            print('Clustering features: %s' % self.clustering_cfg.features)
            print('k finding method: %s' % self.clustering_cfg.k_finding_method)
            print('Number of clusters: %d' % self.clustering_info[cnst.N_CLUST_KEY])
            print('WCSS (%%): %.2f' % (self.clustering_info[cnst.WCSS_KEY] * 10000))
            print('Silhouette score: %.2f' % self.clustering_info[cnst.SIL_SCORE_KEY])
        except KeyError as e:
            raise KeyError(f"Missing key in clustering_info: {e}")
        except AttributeError as e:
            raise AttributeError(f"Missing attribute: {e}")
    
        # Print details on identified phases
        print_single_separator()
        print('Identified phases:')
        # Print stddev in-column for ease of visualization
        try:
            clusters_df = self.clusters_df
            el_fr_feature_key = cnst.AT_FR_DF_KEY if self.clustering_cfg.features == cnst.AT_FR_CL_FEAT else cnst.W_FR_DF_KEY
            fr_labels = [el + el_fr_feature_key for el in self.all_els_sample]
            stddev_labels = [el + cnst.STDEV_DF_KEY + el_fr_feature_key for el in self.all_els_sample]
            df_mod_to_print = []
            for index, row in clusters_df.iterrows():
                els_dict = {}
                df_mod_to_print.append({cnst.N_PTS_DF_KEY: row[cnst.N_PTS_DF_KEY]})
                # Add conversion to atomic fraction when mass fractions are used as features
                if self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
                    at_fr_dict = {}
                    for el in self.all_els_sample:
                        label = el + cnst.AT_FR_DF_KEY
                        at_fr_dict[label] = row[label]
                    df_mod_to_print[-1].update(at_fr_dict)
                # Get elemental fractions (at_fr or w_fr)
                for element, fr_l, stddev_l in zip(self.all_els_sample, fr_labels, stddev_labels):
                    els_dict[element + el_fr_feature_key] = f"{row[fr_l]:.1f} ± {row[stddev_l]:.1f}"
                # Add elemental fractions + cluster stddev and wcss
                df_mod_to_print[-1].update({
                    **els_dict,
                })
                
                # Add cluster-level std dev entries ---
                cluster_stddev_entries = {
                    col: f"{row[col]:.1f}" 
                    for col in clusters_df.columns 
                    if col.startswith(cnst.RMS_DIST_DF_KEY)
                }
                df_mod_to_print[-1].update(cluster_stddev_entries)
                
                # Add references if present
                if self.ref_formulae:
                    ref_keys_to_print = ['Conf1', 'Conf2', 'Ref1', 'Ref2']
                    ref_dict = {key: value for key, value in row.items() if key in ref_keys_to_print}
                    df_mod_to_print[-1].update(ref_dict)
                # Add mixtures to the printed report
                mix_keys_to_print = ['Mix1', 'ConfMix1', 'Mol_Ratio1', 'Mix2', 'ConfMix2', 'Mol_Ratio2']
                mix_dict = {key: value for key, value in row.items() if key in mix_keys_to_print}
                df_mod_to_print[-1].update(mix_dict)
            # Set display options for float precision
            with pd.option_context('display.float_format', '{:,.2f}'.format):
                pd.set_option('display.max_columns', None)  # Display all columns
                print(pd.DataFrame(df_mod_to_print))
        except Exception as e:
            raise RuntimeError(f"Error printing phase results: {e}")
    

    #%% Experimental standard PB ratios
    # =============================================================================
    def _compile_standards_from_references(self) -> dict:
        """
        Compile a standards dictionary for the current sample by using the input
        reference phases, if present in the list of standards.
    
        This function loads the standards library, iterates over all elements in
        the current sample, and for each X-ray reference line:
          - Verifies if the reference phase compositions are present in the standards.
          - If no reference phases are found, a warning is issued and existing
            standards are used.
          - If references are found, the function computes the mean of the
            corrected PB values and substitutes them into the standards dictionary.
    
        Returns:
            dict: The updated standards dictionary, where each entry for a given
                  element-line combination contains either existing standards or
                  a single mean standard to be fed to XSp_Quantifier
    
        Warns:
            UserWarning: If none of the input reference phase compositions are
                    present for a given reference line in the standards file.
                    
        Note
        ----
        Currently only used when analysing mixtures of known powder precursors
        (i.e., powder_meas_cfg.is_known_powder_mixture_meas = True)
        """
        std_dict_all_modes, std_dir = self._load_xsp_standards()
        std_dict_all_lines = std_dict_all_modes[self.measurement_cfg.mode]
        ref_lines = XSp_Quantifier.xray_quant_ref_lines
        ref_formulae = self.clustering_cfg.ref_formulae
        
        filtered_std_dict = {}
        for el in self.detectable_els_sample:
            for line in ref_lines:
                el_line = f"{el}_{line}"
                if el_line not in std_dict_all_lines:
                    continue

                # Gather matching reference entries by comparing chemical formulas
                ref_entries = []
                for i, std_dict in enumerate(std_dict_all_lines[el_line]):
                    if std_dict[cnst.STD_ID_KEY] != cnst.STD_MEAN_ID_KEY:
                        try:
                            std_comp = Composition(std_dict[cnst.STD_FORMULA_KEY])
                            for ref_formula in ref_formulae:
                                if std_comp.reduced_formula == Composition(ref_formula).reduced_formula:
                                    ref_entries += [i]
                        except Exception as e:
                            warnings.warn(
                                f"Could not parse formula '{std_dict[cnst.STD_FORMULA_KEY]}' "
                                f"or compare with reference formulas {ref_formulae}. Error: {e}"
                            )

                if len(ref_entries) < 1 and not self.exp_stds_cfg.is_exp_std_measurement:
                    text_line = "provided standards" if std_dir == "" else f"standards file at: {std_dir}"
                    warnings.warn(
                        f"None of the input reference phases {ref_formulae} "
                        f"is present for line {el_line} in the {text_line}. "
                        "Using other available standards."
                    )
                    continue
    
                # Compute mean PB value from all available references
                new_std_ref_list = [std_d for i, std_d in enumerate(std_dict_all_lines[el_line]) if i in ref_entries]
                list_PB = [ref_line[cnst.COR_PB_DF_KEY] for ref_line in new_std_ref_list]
    
                std_dict_mean = {
                    cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY,
                    cnst.COR_PB_DF_KEY: float(np.mean(list_PB)),
                }
                filtered_std_dict[el_line] = [std_dict_mean]
                
        return filtered_std_dict

    
    def _fit_stds_and_save_results(self, backup_previous_data: bool = True) -> Union[Tuple, None]:
        """
        Fit spectra collected from experimental standards, process results, and save them to disk.
    
        Parameters
        ----------
        backup_previous_data : bool, optional
            Backs up previous data file if present (Default = True).
    
        Returns
        -------
        Tuple or None
            If data was successfully processed:
                - std_ref_lines : Any
                    Data structure containing assembled standard PB data.
                - results_df : pandas.DataFrame
                    DataFrame containing averaged PB ratios and corrected values.
                - Z_sample : Any
                    Sample average atomic number computed with different methods.
            If no measurement data was available:
                Returns `(None, None, None)`.
    
        Raises
        ------
        RuntimeError
            If fitting, saving, or PB correction fails unexpectedly.
        """
        
        # Initialize return variables
        std_ref_lines = None
        results_df = None
        Z_sample = None
        
        try:
            # Fit spectra and assemble results
            self._fit_and_quantify_spectra(quantify=False)
        except Exception as e:
            raise RuntimeError(f"Error during fitting and quantification: {e}") from e
        
        try:
            # Save per-spectrum measurement results
            data_df = self._save_collected_data(
                None,
                None,
                backup_previous_data=backup_previous_data,
                include_spectral_data=True
            )
        except Exception as e:
            raise RuntimeError(f"Error while saving collected data: {e}") from e
        
        if data_df is not None and not data_df.empty:
            try:
                # Assemble PB data and calculate corrections
                std_ref_lines = self._assemble_std_PB_data(data_df)
                if std_ref_lines != {}:
                    PB_corrected, Z_sample = self._calc_corrected_PB(std_ref_lines)
                    
                    # Save averaged PB results
                    results_df = self._save_std_results(std_ref_lines, PB_corrected)
                else:
                    if self.verbose:
                        print("No valid standard measurement acquired.")
            except Exception as e:
                raise RuntimeError(f"Error while processing standard results: {e}") from e
    
            return std_ref_lines, results_df, Z_sample
        
        # No data available to process
        return None, None, None

    
    def _evaluate_exp_std_fit(self, tot_n_spectra: int) -> Tuple[bool, bool]:
        """
        Evaluate the experimental standard fitting results after collecting a given number of spectra.
    
        This method attempts to fit the experimental standards using the currently collected spectra.
        It determines whether the fit was successful and whether the minimum required number of valid
        spectra has been reached. The method also provides verbose feedback if enabled.
    
        Parameters
        ----------
        tot_n_spectra : int
            Total number of spectra collected so far.
    
        Returns
        -------
        Tuple[bool, bool]
            A Tuple containing:
            - is_fit_successful (bool): Whether the fitting process produced valid results.
            - is_converged (bool): Whether the minimum number of valid spectra was reached.
        """
        is_fit_successful = False
        is_converged = False
    
        try:
            if self.verbose:
                print_double_separator()
                print(f"Fitting after collection of {tot_n_spectra} spectra...")
    
            _, results_df, _ = self._fit_stds_and_save_results(backup_previous_data=False)
    
            if results_df is not None and not results_df.empty:
                is_fit_successful = True
    
                # Retrieve the minimum number of valid spectra from the results
                try:
                    num_valid_spectra = int(np.min(results_df[cnst.N_SP_USED_KEY]))
                except (KeyError, ValueError, TypeError) as e:
                    raise RuntimeError(f"Results DataFrame missing or invalid '{cnst.N_SP_USED_KEY}' column.") from e
    
                is_converged = num_valid_spectra >= self.min_n_spectra
    
                if self.verbose:
                    print_double_separator()
                    print("Fitting performed.")
                    print(f"{num_valid_spectra} valid spectra were collected.")
                    if is_converged:
                        print(f"Target number of {self.min_n_spectra} was reached.")
            else:
                if self.verbose:
                    print_double_separator()
                    print("No valid spectrum collected.")
    
            # If not converged, provide feedback
            if not is_converged and self.verbose:
                if tot_n_spectra >= self.max_n_spectra:
                    print(f"Maximum allowed number of {self.max_n_spectra} spectra was acquired, "
                          f"but target number of {self.min_n_spectra} was not reached.")
                else:
                    print(f"More spectra will be collected to reach target number of {self.min_n_spectra}.")
    
        except Exception as e:
            raise RuntimeError("An error occurred while evaluating the experimental standard fit.") from e
    
        return is_fit_successful, is_converged
    
        
    def _assemble_std_PB_data(
        self,
        data_df: "pd.DataFrame"
    ) -> Dict[str, Dict[str, Union[float, List[float]]]]:
        """
        Assemble Peak-to-Background (PB) ratio data for the experimental standard references to use during quantification.
    
        This method processes the provided DataFrame of spectral data to:
        1. Remove any X-ray peaks whose PB ratio is absent or below the acceptable threshold for all spectra.
        2. Exclude spectra that do not meet the accepted quantification flags.
        3. Compile PB ratio statistics (mean, std. dev) and corresponding theoretical energies for each relevant element line.
    
        Parameters
        ----------
        data_df : pd.DataFrame
            DataFrame containing PB ratio measurements and associated metadata.
            Must contain:
            - `cnst.QUANT_FLAG_DF_KEY` column for quantification flags.
            - Columns for PB ratios of element lines (e.g., 'Fe_Ka', 'Cu_Ka', etc.).
            - Possibly NaN values where peaks are absent.
    
        Returns
        -------
        Dict[str, Dict[str, float | List[float]]]
            Dictionary mapping each fitted standard element line to a sub-dictionary containing:
            - cnst.PB_RATIO_KEY: List of measured PB ratios.
            - cnst.MEAN_PB_KEY: mean PB ratio (ignoring NaN).
            - cnst.STDEV_PB_DF_KEY: standard deviation of PB ratios (ignoring NaN).
            - cnst.PEAK_TH_ENERGY_KEY: theoretical peak energy for that element line.
    
        Notes
        -----
        - Assumes that `self._th_peak_energies` is a dictionary mapping element lines to their theoretical energies.
        """
        # Filter out X-ray peaks whose PB ratio is absent for all spectra
        data_df = data_df.dropna(axis=1, how="all")
        # Filter out rows corresponding to spectra that should be discarded
        try:
            data_filtered_df = data_df[data_df[cnst.QUANT_FLAG_DF_KEY].isin(self.exp_stds_cfg.quant_flags_accepted)]
        except KeyError as e:
            raise RuntimeError(f"Missing required column '{cnst.QUANT_FLAG_DF_KEY}' in input DataFrame.") from e
            
        # Get fitted element lines for elements in the standard
        all_fitted_el_lines = [
            el_line for el_line in self._th_peak_energies.keys()
            if el_line in data_filtered_df.columns
        ]
        fitted_std_el_lines = [
            el_line for el_line in all_fitted_el_lines
            if el_line.split("_")[0] in self.detectable_els_sample
        ]

        # Update lists of measured PB ratios, their means, stddev, and corresponding theoretical energies
        std_ref_lines = {}
        for el_line in fitted_std_el_lines:
            meas_PB_ratios = data_filtered_df[el_line].tolist()
            if len(meas_PB_ratios) > 0:
                std_ref_lines[el_line] = {
                    cnst.PB_RATIO_KEY: meas_PB_ratios,
                    cnst.MEAN_PB_KEY: float(np.nanmean(meas_PB_ratios)),
                    cnst.STDEV_PB_DF_KEY: float(np.nanstd(meas_PB_ratios)),
                    cnst.PEAK_TH_ENERGY_KEY: self._th_peak_energies[el_line]
                }

        return std_ref_lines
        
        
    def _calc_corrected_PB(
        self,
        std_ref_lines: Dict[str, Dict[str, Union[float, List[float]]]]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate ZAF-corrected Peak-to-Background (PB) ratios for the experimental standard.
    
        This method applies ZAF (atomic number, bnackscattering, absorption) corrections
        to the measured PB ratios of the standard's element lines. The corrected PB ratios are normalized
        by the mass fraction of each element to obtain the pure element PB ratios.
    
        Parameters
        ----------
        std_ref_lines : Dict[str, Dict[str, float | List[float]]]
            Dictionary mapping each element line (e.g., 'Fe_Ka') to its PB ratio statistics and theoretical peak energy.
            Must contain:
            - cnst.PEAK_TH_ENERGY_KEY: float, theoretical peak energy.
            - cnst.MEAN_PB_KEY: float, mean measured PB ratio.
            - Corresponding element's mass fraction in `self.exp_stds_cfg.w_frs`.
    
        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            - PB_corrected (np.ndarray): ZAF-corrected PB ratios, normalized by mass fractions.
            - Z_sample (np.ndarray): ZAF correction factors for the sample.
    
        Notes
        -----
        - Expects std_ref_lines to be non-empty
        - Relies on `self.exp_stds_cfg.w_frs` for element mass fractions.
        - Uses `Quant_Corrections.get_ZAF_mult_f_pb` for ZAF factor computation.
        - Assumes `self.detectable_els_sample` contains the list of detectable element symbols.
        """
        peak_energies_dict: Dict[str, float] = {}
        means_PB: List[float] = []
        w_frs: List[float] = []
    
        # Extract peak energies, mean PB ratios, and corresponding mass fractions
        for el_line, el_line_dict in std_ref_lines.items():
            try:
                peak_energies_dict[el_line] = el_line_dict[cnst.PEAK_TH_ENERGY_KEY]
                means_PB.append(float(el_line_dict[cnst.MEAN_PB_KEY]))
            except KeyError as e:
                raise RuntimeError(f"Missing expected key in std_ref_lines for element line '{el_line}'.") from e
    
            el = el_line.split('_')[0]
            try:
                w_frs.append(self.exp_stds_cfg.w_frs[el])
            except KeyError as e:
                raise RuntimeError(f"Mass fraction for element '{el}' not found in exp_stds_cfg.w_frs.") from e
    
        # Initialize ZAF correction calculator (second-order corrections for PB method)
        ZAF_calculator = Quant_Corrections(
            elements=self.detectable_els_sample,
            beam_energy=self.measurement_cfg.beam_energy_keV,
            emergence_angle=self.measurement_cfg.emergence_angle,
            meas_mode=self.measurement_cfg.mode,
            verbose=False
        )
    
        # Get nominal mass fractions for all detectable elements
        missing_elements = [el for el in self.detectable_els_sample if el not in self.exp_stds_cfg.w_frs]
        if missing_elements:
            raise RuntimeError(
                f"Missing mass fraction(s) for detectable element(s): {', '.join(missing_elements)}."
            )
        
        nominal_w_frs = [self.exp_stds_cfg.w_frs[el] for el in self.detectable_els_sample]

        # Calculate ZAF corrections
        # Returns arrays with dimensions corresponding to w_frs
        ZAF_pb, Z_sample = ZAF_calculator.get_ZAF_mult_f_pb(nominal_w_frs, peak_energies_dict)
    
        # Apply ZAF correction and normalize by mass fractions to get pure element PB ratios
        PB_corrected = ZAF_pb * np.array(means_PB) / np.array(w_frs)
    
        return PB_corrected, Z_sample
    
    
    def _save_std_results(
        self,
        std_ref_lines: Dict[str, Dict[str, Any]],
        PB_corrected: List[float]
    ) -> Optional[pd.DataFrame]:
        """
        Save and return a summary table of standard reference line results.
    
        Constructs a table summarizing the mean, standard deviation, relative error,
        and number of spectra used for each reference line, then saves it as a CSV file.
    
        Parameters
        ----------
        std_ref_lines : Dict
            Dictionary mapping line identifiers (str) to dictionaries containing
            statistical results for each line. Each inner dictionary must contain keys
            for mean, standard deviation, and a list or Series of PB ratios.
        PB_corrected : list of float
            List of corrected PB values, one per reference line (order must match keys of std_ref_lines).
    
        Returns
        -------
        results_df : pandas.DataFrame or None
            DataFrame with summary statistics for each reference line, or None if no lines are provided.
    
        Raises
        ------
        ValueError
            If the number of PB_corrected values does not match the number of reference lines.
    
        Notes
        -----
        The DataFrame is saved as a CSV file in Std_measurements, with a filename
        including the measurement mode and output filename suffix.
        """
        if not std_ref_lines:
            return None
    
        # Extract statistics for each reference line
        means_PB = []
        stdevs_PB = []
        n_spectra_per_line = []
        line_keys = list(std_ref_lines.keys())
    
        if len(PB_corrected) != len(line_keys):
            raise ValueError("Length of PB_corrected does not match number of reference lines.")
    
        for el_line in line_keys:
            el_line_dict = std_ref_lines[el_line]
            means_PB.append(el_line_dict[cnst.MEAN_PB_KEY])
            stdevs_PB.append(el_line_dict[cnst.STDEV_PB_DF_KEY])
            pb_ratios = el_line_dict[cnst.PB_RATIO_KEY]
            n_spectra_used = sum(
                (x is not None) and (not (isinstance(x, float) and np.isnan(x)))
                for x in pb_ratios
            )
            n_spectra_per_line.append(n_spectra_used)
    
        # Construct the results DataFrame
        results_df = pd.DataFrame({
            cnst.MEAS_PB_DF_KEY: means_PB,
            cnst.STDEV_PB_DF_KEY: stdevs_PB,
            cnst.COR_PB_DF_KEY: PB_corrected,
            cnst.REL_ER_PERCENT_PB_DF_KEY: np.array(stdevs_PB) / np.array(means_PB) * 100,
            cnst.N_SP_USED_KEY: n_spectra_per_line
        }, index=line_keys)
    
        # Save the DataFrame as CSV
        filename = f"{cnst.STDS_RESULT_FILENAME}_{self.measurement_cfg.mode}" + self.output_filename_suffix
        results_path = os.path.join(self.sample_result_dir, filename + '.csv')
        results_df.to_csv(results_path, index=True, header=True)
    
        return results_df
    
    
    def _load_xsp_standards(self) -> Tuple[dict, str]:
        """
        Load the X-ray Spectroscopy standards library for the current measurement configuration.
        
        This function attempts to load an existing standards library based on the
        measurement type and beam energy defined in the measurement configuration.
        If the library cannot be found, a new empty dictionary is created for the
        current measurement mode. If loading fails due to an unexpected error, a
        RuntimeError is raised with the original exception preserved.
        
        Returns:
            tuple[dict, str]: 
                A tuple containing:
                - standards (dict): The standards library, indexed by measurement mode.
                - std_dir (str): The directory path where the standards are stored.
        
        Raises:
            RuntimeError: If loading the standards library fails due to an 
                          unexpected error (not just missing files).
        """
        meas_mode = self.measurement_cfg.mode
        # Load or create standards dictionary
        if self.standards_dict is None:
            try:
                standards, std_dir = calibs.load_standards(self.measurement_cfg.type, self.measurement_cfg.beam_energy_keV)
            except FileNotFoundError:
                std_dir = calibs.standards_dir
                standards = {meas_mode: {}}
            except Exception as e:
                raise RuntimeError("Failed to load standards library.") from e
        else:
            standards = self.standards_dict
            std_dir = ''

        # Ensure measurement mode exists in the standards dictionary, otherwise create it
        if meas_mode not in standards:
            standards[meas_mode] = {}
            
        return standards, std_dir
    
    
    def _update_standard_library(
        self,
        std_ref_lines: Dict[str, Dict[str, Union[float, List[float]]]],
        results_df: pd.DataFrame,
        Z_sample: np.ndarray
    ) -> None:
        """
        Update the standards library with new Peak-to-Background (PB) ratio measurements.
    
        This method:
        1. Loads the current standards library from disk (or creates a new one if missing).
        2. Removes any previous entries for the current standard.
        3. Appends the new measurements for each element line.
        4. Recalculates the 'Mean' reference PB ratio and associated uncertainty for each element line.
        5. Saves the updated standards library back to disk.
    
        Parameters
        ----------
        std_ref_lines : Dict[str, Dict[str, float | List[float]]]
            Dictionary mapping element lines (e.g., 'Fe_Ka') to PB ratio data and metadata.
        results_df : pd.DataFrame
            DataFrame containing measured PB ratios, corrected PB ratios, standard deviations,
            and relative errors for each element line.
        Z_sample : np.ndarray
            Mean sampel atomic number, computed using different methods.
    
        Raises
        ------
        RuntimeError
            If the standards library cannot be loaded or updated due to missing keys or invalid data.
        """
        meas_mode = self.measurement_cfg.mode
        
        # Load standards
        if self.standards_dict is not None:
            warnings.warn("The 'standards_dict' provided when initializing EMXSp_Composition_Analyzer will be ignored."
                          f"Loading standards dictionary from XSp_calibs/{self.microscope_cfg.ID}", UserWarning())
            self.standards_dict = None
        standards, std_dir = self._load_xsp_standards()
        
        std_lib = standards[meas_mode]
        
        # Remove all previous entries measured from this standard
        was_standard_already_measured = False
        for el_line, stds_list in list(std_lib.items()):
            for i, std_dict in enumerate(list(stds_list)):
                if std_dict.get(cnst.STD_ID_KEY) == self.sample_cfg.ID:
                    std_lib[el_line].pop(i)
                    was_standard_already_measured = True
                    break
        if was_standard_already_measured and self.verbose:
            print_single_separator()
            print(f"Previously measured values for standard '{self.sample_cfg.ID}' were found and removed.")
    
        # Add new standards
        now = datetime.now()
        for el_line in std_ref_lines.keys():
            # Validate presence of required result_df fields
            for key in [
                cnst.COR_PB_DF_KEY,
                cnst.MEAS_PB_DF_KEY,
                cnst.STDEV_PB_DF_KEY,
                cnst.REL_ER_PERCENT_PB_DF_KEY
            ]:
                if key not in results_df.columns:
                    raise RuntimeError(f"Missing required column '{key}' in results_df.")
    
            std_dict_new = {
                cnst.STD_ID_KEY: self.sample_cfg.ID,
                cnst.STD_FORMULA_KEY: self.exp_stds_cfg.formula,
                cnst.STD_TYPE_KEY: self.sample_cfg.type,
                cnst.DATETIME_KEY: now.strftime("%Y-%m-%d %H:%M:%S"),
                cnst.COR_PB_DF_KEY: results_df.at[el_line, cnst.COR_PB_DF_KEY],
                cnst.MEAS_PB_DF_KEY: results_df.at[el_line, cnst.MEAS_PB_DF_KEY],
                cnst.STDEV_PB_DF_KEY: results_df.at[el_line, cnst.STDEV_PB_DF_KEY],
                cnst.REL_ER_PERCENT_PB_DF_KEY: results_df.at[el_line, cnst.REL_ER_PERCENT_PB_DF_KEY],
                cnst.STD_USE_FOR_MEAN_KEY : self.exp_stds_cfg.use_for_mean_PB_calc,
                cnst.STD_Z_KEY: Z_sample
            }
    
            # Add or append standard measurement
            if el_line in std_lib:
                if self.verbose:
                    print(f"Added the measured standard PB value for {el_line} to the current list.")
                std_lib[el_line].append(std_dict_new)
            else:
                if self.verbose:
                    print(f"Created a new list for the {el_line} line PB standard values.")
                std_lib[el_line] = [std_dict_new]
    
            # Recalculate mean of standards (excluding previous mean entry)
            std_el_line_entries = [
                std for std in std_lib[el_line]
                if std.get(cnst.STD_ID_KEY) != cnst.STD_MEAN_ID_KEY
            ]
            # Select corrected PB ratios that should be used for calculating the mean (i.e., PB ratios computed from the mean)
            list_PB_for_mean = [std[cnst.COR_PB_DF_KEY] for std in std_el_line_entries if std[cnst.STD_USE_FOR_MEAN_KEY]]
            if len(list_PB_for_mean) > 0: 
                mean_PB = float(np.mean(list_PB_for_mean)) if list_PB_for_mean else float("nan")
                stddev_mean_PB = float(np.std(list_PB_for_mean)) if list_PB_for_mean else float("nan")
                error_mean_PB = (stddev_mean_PB / mean_PB * 100) if mean_PB else float("nan")
        
                std_dict_mean = {
                    cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY,
                    cnst.DATETIME_KEY: now.strftime("%Y-%m-%d %H:%M:%S"),
                    cnst.COR_PB_DF_KEY: mean_PB,
                    cnst.STDEV_PB_DF_KEY: stddev_mean_PB,
                    cnst.REL_ER_PERCENT_PB_DF_KEY: error_mean_PB
                }
                std_el_line_entries.append(std_dict_mean)
            std_lib[el_line] = std_el_line_entries
    
        # Save updated file with standards
        try:
            with open(std_dir, "w") as file:
                json.dump(standards, file, indent=2)
        except Exception as e:
            raise RuntimeError(f"Failed to save updated standards to {std_dir}.") from e
 