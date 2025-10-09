#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated X-Ray Experimental Standard Acquisition and Analysis

This module configures and runs automated collection and fitting
of EDS/WDS spectra from experimental standards (i.e., sampels of known composition)
to generate reference values of peak-to-background ratios.

Import this module in your own Python code and call the
`batch_acquire_experimental_stds()` function, passing your desired configuration
and sample list as arguments. This enables integration into larger
automation workflows or pipelines.

Workflow includes:
    - Configuration of microscope, measurement, substrate, and fitting parameters
    - Sample setup (elements, position, reference formulae)
    - Calibration and quantification options
    - Automated or manual navigation and acquisition modes

Requirements:
    - Proper instrument calibration files and instrument driver for the selected microscope

Typical usage:
    - Edit the 'std_list' list to define your standards
    - Adjust configuration parameters as needed, either directly in the script or by passing them to `batch_acquire_experimental_stds`
    - Run the script, or import and call `batch_acquire_experimental_stds()` to collect experimental standards for one or multiple samples at a time

Parameters
----------
stds : list of dict
    List of experimental standard definitions.  
    Each dictionary must contain:
        - 'ID' (str): Identifier for the standard sample (SampleConfig.ID).
        - 'formula' (str): Chemical formula of the standard (ExpStandardsConfig.formula).
        - 'pos' (tuple of float): (x, y) stage coordinates in mm (SampleConfig.center_pos).
        - 'sample_type' (str): Sample type ('powder', 'bulk', etc.) (SampleConfig.type).
        - 'is_manual_meas' (bool): If True, navigation to positions is manual (MeasurementConfig.is_manual_navigation).
microscope_ID : str, optional
    Identifier for the microscope hardware.  
    Must correspond to a calibration folder in `./XSp_calibs/Microscopes/<ID>` (MicroscopeConfig.ID).  
    Default is `'PhenomXL'`.
microscope_type : str, optional
    Type of microscope. Allowed: `'SEM'` (implemented), `'STEM'` (not implemented).  
    Default is `'SEM'` (MicroscopeConfig.type).
measurement_type : str, optional
    Measurement type. Allowed: `'EDS'` (implemented), `'WDS'` (not implemented).  
    Default is `'EDS'` (MeasurementConfig.type).
measurement_mode : str, optional
    Acquisition mode (e.g., `'point'`, `'map'`), defining beam/detector calibration settings.  
    Default is `'point'` (MeasurementConfig.mode).
sample_halfwidth : float, optional
    Half-width of the sample area in mm for mapping/acquisition.  
    Default is `3.0` (SampleConfig.half_width_mm).
sample_substrate_type : str, optional
    Type of sample substrate. Allowed: `'Ctape'`, `'None'`.  
    Default is `'Ctape'` (SampleSubstrateConfig.type).
sample_substrate_shape : str, optional
    Shape of the substrate. Allowed: `'circle'`, `'square'`.  
    Default is `'circle'` (SampleSubstrateConfig.shape).
working_distance : float, optional
    Working distance in mm for acquisition. If None, taken from microscope driver.  
    Default is `5.0` (MeasurementConfig.working_distance).
beam_energy : float, optional
    Electron beam energy in keV.  
    Default is `15.0` (MeasurementConfig.beam_energy_keV).
spectrum_lims : tuple of float, optional
    Lower and upper energy limits for spectrum fitting in eV.  
    Default is `(14, 1100)` (QuantConfig.spectrum_lims).
use_instrument_background : bool, optional
    Whether to use instrument background files during fitting.  
    If False, background is computed during fitting.  
    Default is `False` (QuantConfig.use_instrument_background).
min_bckgrnd_cnts : float, optional
    Minimum background counts required for a spectrum not to be filtered out.  
    Default is `5` (QuantConfig.min_bckgrnd_cnts).
fit_during_collection : bool, optional
    If True, fit spectra during acquisition; otherwise fit later.  
    Default is `True`.
update_std_library : bool, optional
    If True, update the stored experimental standards library after acquisition.  
    Default is `False`.
is_auto_substrate_detection : bool, optional
    If True, substrate elements are detected automatically.  
    Implemented only for `'Ctape'` substrates.  
    Default is `False` (SampleSubstrateConfig.auto_detection).
auto_adjust_brightness_contrast : bool, optional
    If True, brightness/contrast are set automatically.  
    Default is `True` (MicroscopeConfig.is_auto_BC).
contrast : float, optional
    Manual contrast setting (required if auto_adjust_brightness_contrast is False).  
    Default is `4.3877` (MicroscopeConfig.contrast).
brightness : float, optional
    Manual brightness setting (required if auto_adjust_brightness_contrast is False).  
    Default is `0.4504` (MicroscopeConfig.brightness).
min_n_spectra : int, optional
    Minimum number of spectra to acquire.  
    Default is `50` (MeasurementConfig.min_n_spectra).
max_n_spectra : int, optional
    Maximum number of spectra to acquire.  
    Default is `100` (MeasurementConfig.max_n_spectra).
target_Xsp_counts : int, optional
    Target counts for spectrum acquisition.  
    Default is `250000` (MeasurementConfig.target_acquisition_counts).
max_XSp_acquisition_time : float, optional
    Maximum acquisition time in seconds. If None, estimated from target counts.  
    Default is `None` (MeasurementConfig.max_acquisition_time).
els_substrate : list of str, optional
    List of substrate element symbols.  
    Default is `['C', 'O', 'Al']` (SampleSubstrateConfig.elements).
powder_meas_cfg_kwargs : dict, optional
    Additional keyword arguments for PowderMeasurementConfig.
bulk_meas_cfg_kwargs : dict, optional
    Additional keyword arguments for BulkMeasurementConfig.
exp_stds_meas_cfg_kwargs : dict, optional
    Additional keyword arguments for ExpStandardsConfig.  
    Used to customize experimental standard acquisition and PB ratio filtering.
output_filename_suffix : str, optional
    String appended to output filenames.  
    Default is `''`.
development_mode : bool, optional
    If True, enables development/debug features.  
    Default is `False`.
verbose : bool, optional
    If True, print verbose output.  
    Default is `True`.
exp_std_dir : str, optional
    Directory where experimental standard results are saved.  
    Default is `None`.
    
    
Created on Fri Aug 20 09:34:34 2025

@author: Andrea
"""

import logging
from typing import List, Dict, Tuple, Any

from autoemxsp.core.EMXSp_composition_analyser import EMXSp_Composition_Analyzer
import autoemxsp.XSp_calibs as calibs
from autoemxsp.tools.utils import print_double_separator
from autoemxsp.tools.config_classes import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    QuantConfig,
    ClusteringConfig,
    PowderMeasurementConfig,
    BulkMeasurementConfig,
    ExpStandardsConfig
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def batch_acquire_experimental_stds(
    stds: List[Dict[str, Any]],
    microscope_ID: str = 'PhenomXL',
    microscope_type: str = 'SEM',
    measurement_type: str = 'EDS',
    measurement_mode: str = 'point',
    sample_halfwidth: float = 3.0,
    sample_substrate_type: str = 'Ctape',
    sample_substrate_shape: str = 'circle',
    working_distance: float = 5, #mm
    beam_energy: float = 15.0,
    spectrum_lims: Tuple[float, float] = (14, 1100),
    use_instrument_background: bool = False,
    min_bckgrnd_cnts: float = 5,
    fit_during_collection= True,
    update_std_library = False,
    is_auto_substrate_detection: bool = False,
    auto_adjust_brightness_contrast: bool = True,
    contrast: float = 4.3877,
    brightness: float = 0.4504,
    min_n_spectra: int = 50,
    max_n_spectra: int = 100,
    target_Xsp_counts: int = 250000,
    max_XSp_acquisition_time: float = None,
    els_substrate: List[str] = None,
    powder_meas_cfg_kwargs: Dict[str, Any] = None,
    bulk_meas_cfg_kwargs: Dict[str, Any] = None,
    exp_stds_meas_cfg_kwargs: Dict[str, Any] = None,
    output_filename_suffix: str = '',
    development_mode: bool = False,
    verbose: bool = True,
    exp_std_dir: str = None,
) -> None:
    """
    Batch acquisition (and optional quantification) of X-ray spectra for a list of powder samples.

    Parameters
    ----------
    stds : list of dict
        List of experimental standard definitions.  
        Each dictionary must contain:
            - 'ID' (str): Identifier for the standard sample (SampleConfig.ID).
            - 'formula' (str): Chemical formula of the standard (ExpStandardsConfig.formula).
            - 'pos' (tuple of float): (x, y) stage coordinates in mm (SampleConfig.center_pos).
            - 'sample_type' (str): Sample type ('powder', 'bulk', etc.) (SampleConfig.type).
            - 'is_manual_meas' (bool): If True, navigation to positions is manual (MeasurementConfig.is_manual_navigation).
    microscope_ID : str, optional
        Identifier for the microscope hardware.  
        Must correspond to a calibration folder in `./XSp_calibs/Microscopes/<ID>` (MicroscopeConfig.ID).  
        Default is `'PhenomXL'`.
    microscope_type : str, optional
        Type of microscope. Allowed: `'SEM'` (implemented), `'STEM'` (not implemented).  
        Default is `'SEM'` (MicroscopeConfig.type).
    measurement_type : str, optional
        Measurement type. Allowed: `'EDS'` (implemented), `'WDS'` (not implemented).  
        Default is `'EDS'` (MeasurementConfig.type).
    measurement_mode : str, optional
        Acquisition mode (e.g., `'point'`, `'map'`), defining beam/detector calibration settings.  
        Default is `'point'` (MeasurementConfig.mode).
    sample_halfwidth : float, optional
        Half-width of the sample area in mm for mapping/acquisition.  
        Default is `3.0` (SampleConfig.half_width_mm).
    sample_substrate_type : str, optional
        Type of sample substrate. Allowed: `'Ctape'`, `'None'`.  
        Default is `'Ctape'` (SampleSubstrateConfig.type).
    sample_substrate_shape : str, optional
        Shape of the substrate. Allowed: `'circle'`, `'square'`.  
        Default is `'circle'` (SampleSubstrateConfig.shape).
    working_distance : float, optional
        Working distance in mm for acquisition. If None, taken from microscope driver.  
        Default is `5.0` (MeasurementConfig.working_distance).
    beam_energy : float, optional
        Electron beam energy in keV.  
        Default is `15.0` (MeasurementConfig.beam_energy_keV).
    spectrum_lims : tuple of float, optional
        Lower and upper energy limits for spectrum fitting in eV.  
        Default is `(14, 1100)` (QuantConfig.spectrum_lims).
    use_instrument_background : bool, optional
        Whether to use instrument background files during fitting.  
        If False, background is computed during fitting.  
        Default is `False` (QuantConfig.use_instrument_background).
    min_bckgrnd_cnts : float, optional
        Minimum background counts required for a spectrum not to be filtered out.  
        Default is `5` (QuantConfig.min_bckgrnd_cnts).
    fit_during_collection : bool, optional
        If True, fit spectra during acquisition; otherwise fit later.  
        Default is `True`.
    update_std_library : bool, optional
        If True, update the stored experimental standards library after acquisition.  
        Default is `False`.
    is_auto_substrate_detection : bool, optional
        If True, substrate elements are detected automatically.  
        Implemented only for `'Ctape'` substrates.  
        Default is `False` (SampleSubstrateConfig.auto_detection).
    auto_adjust_brightness_contrast : bool, optional
        If True, brightness/contrast are set automatically.  
        Default is `True` (MicroscopeConfig.is_auto_BC).
    contrast : float, optional
        Manual contrast setting (required if auto_adjust_brightness_contrast is False).  
        Default is `4.3877` (MicroscopeConfig.contrast).
    brightness : float, optional
        Manual brightness setting (required if auto_adjust_brightness_contrast is False).  
        Default is `0.4504` (MicroscopeConfig.brightness).
    min_n_spectra : int, optional
        Minimum number of spectra to acquire.  
        Default is `50` (MeasurementConfig.min_n_spectra).
    max_n_spectra : int, optional
        Maximum number of spectra to acquire.  
        Default is `100` (MeasurementConfig.max_n_spectra).
    target_Xsp_counts : int, optional
        Target counts for spectrum acquisition.  
        Default is `250000` (MeasurementConfig.target_acquisition_counts).
    max_XSp_acquisition_time : float, optional
        Maximum acquisition time in seconds. If None, estimated from target counts.  
        Default is `None` (MeasurementConfig.max_acquisition_time).
    els_substrate : list of str, optional
        List of substrate element symbols.  
        Default is `['C', 'O', 'Al']` (SampleSubstrateConfig.elements).
    powder_meas_cfg_kwargs : dict, optional
        Additional keyword arguments for PowderMeasurementConfig.
    bulk_meas_cfg_kwargs : dict, optional
        Additional keyword arguments for BulkMeasurementConfig.
    exp_stds_meas_cfg_kwargs : dict, optional
        Additional keyword arguments for ExpStandardsConfig.  
        Used to customize experimental standard acquisition and PB ratio filtering.
    output_filename_suffix : str, optional
        String appended to output filenames.  
        Default is `''`.
    development_mode : bool, optional
        If True, enables development/debug features.  
        Default is `False`.
    verbose : bool, optional
        If True, print verbose output.  
        Default is `True`.
    exp_std_dir : str, optional
        Directory where experimental standard results are saved.  
        Default is `None`.
            
    Returns
    -------
    results : list(EMXSp_Composition_Analyzer)
        A list of the composition analysis objects (one per sample) containing the results and methods for further analysis.
    """
    if max_XSp_acquisition_time is None:
        max_XSp_acquisition_time = target_Xsp_counts / 10000 * 5
    if els_substrate is None:
        els_substrate = ['C', 'O', 'Al']
    
    # --- Configuration objects
    microscope_cfg = MicroscopeConfig(
        ID=microscope_ID,
        type=microscope_type,
        is_auto_BC=auto_adjust_brightness_contrast,
        brightness=brightness,
        contrast=contrast
    )
    # Load microscope calibrations for this instrument and mode
    calibs.load_microscope_calibrations(microscope_ID, measurement_mode)


    quant_cfg = QuantConfig(
        spectrum_lims=spectrum_lims,
        use_instrument_background=use_instrument_background,
        min_bckgrnd_cnts=min_bckgrnd_cnts
    )

    if powder_meas_cfg_kwargs:
        powder_meas_cfg = PowderMeasurementConfig(**powder_meas_cfg_kwargs)
    else:
        powder_meas_cfg = PowderMeasurementConfig()
        
    if bulk_meas_cfg_kwargs:
        bulk_meas_cfg = BulkMeasurementConfig(**bulk_meas_cfg_kwargs)
    else:
        bulk_meas_cfg = BulkMeasurementConfig()
        

    sample_substrate_cfg = SampleSubstrateConfig(
        elements=els_substrate,
        type=sample_substrate_type,
        shape=sample_substrate_shape,
        auto_detection=is_auto_substrate_detection
    )
    
    results = []
    
    for std_sample in stds:
        # --- Sample configuration
        sample_ID = std_sample['ID']
        formula = std_sample['formula']
        center_pos = std_sample['pos']
        sample_type = std_sample['sample_type']
        is_manual_meas = std_sample['is_manual_meas']
        
        print_double_separator()
        logging.info(f"Sample '{sample_ID}'")
        
        if exp_stds_meas_cfg_kwargs:
            exp_stds_cfg = ExpStandardsConfig(is_exp_std_measurement = True, formula = formula, **exp_stds_meas_cfg_kwargs)
        else:
            exp_stds_cfg = ExpStandardsConfig(is_exp_std_measurement = True, formula = formula)
        
        measurement_cfg = MeasurementConfig(
            type=measurement_type,
            mode=measurement_mode,
            working_distance = working_distance,
            beam_energy_keV=beam_energy,
            is_manual_navigation=is_manual_meas,
            max_acquisition_time=max_XSp_acquisition_time,
            target_acquisition_counts=target_Xsp_counts,
            min_n_spectra=min_n_spectra,
            max_n_spectra=max_n_spectra
        )
        
        
        elements = list(exp_stds_cfg.w_frs.keys())
        sample_cfg = SampleConfig(
            ID=sample_ID,
            elements=elements,
            type=sample_type,
            center_pos=center_pos,
            half_width_mm=sample_halfwidth
        )
        
        # # --- Template: Customizing Parameters Per Sample
        
        # # For any parameter you wish to override on a per-sample basis, add a key:value
        # # pair to the sample's dictionary (e.g., 'n_spectra', 'beam_energy', etc.).
        # # In the main loop, check if the key is present; if not, use the default value.
        # #
        # # Example: To override 'n_spectra' for specific samples,
        # # add 'n_spectra': <value> in the sample dictionary.
        # # Then, after loading the sample parameters, use the following pattern:
        
        # # -- Inside your main sample loop:
        # if sample.get('max_n_spectra') is not None:
        #     max_n_spectra_val = sample['max_n_spectra']
        # else:
        #     max_n_spectra_val = max_n_spectra  # or your chosen default
        
        # # -- Repeat similarly for other parameters you wish to customize per sample:
        # if sample.get('min_n_spectra') is not None:
        #     min_n_spectra_val = sample['min_n_spectra']
        # else:
        #     min_n_spectra_val = min_n_spectra  # or your chosen default
        
        # # -- Make sure to place the configuration object instantiation *after*
        # #    these assignments so that each sample's settings are correctly applied.
        # #    For example:
        # measurement_cfg = MeasurementConfig(
        #     min_n_spectra=min_n_spectra,
        #     max_n_spectra=max_n_spectra
        #     # ... other parameters ...
        # )


        # --- Run Composition Analyzer
        comp_analyzer = EMXSp_Composition_Analyzer(
            microscope_cfg=microscope_cfg,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            quant_cfg=quant_cfg,
            clustering_cfg=ClusteringConfig(),
            powder_meas_cfg=powder_meas_cfg,
            bulk_meas_cfg=bulk_meas_cfg,
            exp_stds_cfg=exp_stds_cfg,
            is_acquisition=True,
            development_mode=development_mode,
            output_filename_suffix=output_filename_suffix,
            verbose=verbose,
            results_dir=exp_std_dir
        )
        
        try:
            comp_analyzer.run_exp_std_collection(fit_during_collection= fit_during_collection, update_std_library = update_std_library)
            results.append(comp_analyzer)
        except Exception as e:
            results.append(None)
            logging.exception(f"Sample '{sample_ID}': acquisition/fitting failed: {e}")
            continue
        
    
    # Put microscope in standby after completion
    if not development_mode and len(stds) > 1:
        try:
            comp_analyzer.EM_controller.standby()
        except Exception as e:
            logging.warning(f"Could not put microscope in standby: {e}")
    
    return results
    
    
    