#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 10:14:12 2025

@author: Andrea
"""

import logging
from typing import List, Dict, Any

from autoemxsp.core.EMXSp_composition_analyser import EMXSp_Composition_Analyzer
from autoemxsp.tools.utils import print_double_separator
from autoemxsp.tools.config_classes import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    PowderMeasurementConfig,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def collect_particle_statistics(
    samples: List[Dict[str, Any]],
    n_par_target: int,
    microscope_ID: str = 'PhenomXL',
    microscope_type: str = 'SEM',
    detector_type: str = 'BSD',
    sample_halfwidth: float = 3.0,
    sample_substrate_type: str = 'Ctape',
    sample_substrate_shape: str = 'circle',
    sample_substrate_width_mm: float = 12,
    working_distance: float = 5, #mm
    is_manual_navigation: bool = False,
    is_auto_substrate_detection: bool = False,
    auto_adjust_brightness_contrast: bool = True,
    contrast: float = 4.3877,
    brightness: float = 0.4504,
    powder_meas_cfg_kwargs: Dict[str, Any] = None,
    output_filename_suffix: str = '',
    development_mode: bool = False,
    verbose: bool = True,
    results_dir: str = None
) -> None:
    """
    Batch acquisition (and optional quantification) of X-ray spectra for a list of powder samples.

    Parameters
    ----------
    samples : list of dict
        List of sample definitions. Each dictionary must contain:
            - 'ID' (str): Sample identifier (SampleConfig.ID).
            - 'pos' (tuple of float): (x, y) stage coordinates in mm (SampleConfig.center_pos).
    n_par_target : int
        Target amount of particles to analyse
    microscope_ID : str, optional
        Identifier for the microscope hardware.
        Must correspond to a calibration folder in `./XSp_calibs/Microscopes/<ID>` (MicroscopeConfig.ID).
        Default is `'PhenomXL'`.
    microscope_type : str, optional
        Type of microscope. Allowed: `'SEM'` (implemented), `'STEM'` (not implemented).
        Default is `'SEM'` (MicroscopeConfig.type).
    detector_type: str, optional
        Default : BSD
    sample_halfwidth : float, optional
        Half-width of the sample area in mm for mapping/acquisition.
        Default is `3.0` (SampleConfig.half_width_mm).
    sample_substrate_type : str, optional
        Type of sample substrate. Allowed: `'Ctape'`, `'None'`.
        Default is `'Ctape'` (SampleSubstrateConfig.type).
    sample_substrate_shape : str, optional
        Shape of the substrate. Allowed: `'circle'`, `'square'`.
        Default is `'circle'` (SampleSubstrateConfig.shape).
    sample_substrate_width_mm : float, optional
        Lateral dimension of substrate holder in mm (SampleSubstrateConfig.stub_w_mm).
    working_distance : float, optional
        Working distance in mm for acquisition. If None, taken from microscope driver.
        Default is `5.0` (MeasurementConfig.working_distance).
    is_manual_navigation : bool, optional
        If True, navigation to sample positions is manual.
        Default is `False` (MeasurementConfig.is_manual_navigation).
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
    powder_meas_cfg_kwargs : dict, optional
        Additional keyword arguments for PowderMeasurementConfig.
    output_filename_suffix : str, optional
        String appended to output filenames.
        Default is `''`.
    development_mode : bool, optional
        If True, enables development/debug features.
        Default is `False`.
    verbose : bool, optional
        If True, print verbose output.
        Default is `True`.
    results_dir : str, optional
        Directory where results are saved.
        Default is `None`.
    
            
    Returns
    -------
    comp_analyzer : EMXSp_Composition_Analyzer
        The composition analysis object containing the results and methods for further analysis.
    """
    
    # --- Configuration objects
    microscope_cfg = MicroscopeConfig(
        ID=microscope_ID,
        type=microscope_type,
        detector_type=detector_type,
        is_auto_BC=auto_adjust_brightness_contrast,
        brightness=brightness,
        contrast=contrast
    )

    measurement_cfg = MeasurementConfig(
        type=MeasurementConfig.PARTICLE_STATS_MEAS_TYPE_KEY,
        working_distance = working_distance,
        is_manual_navigation=is_manual_navigation,
    )

    if powder_meas_cfg_kwargs:
        powder_meas_cfg = PowderMeasurementConfig(**powder_meas_cfg_kwargs)
    else:
        powder_meas_cfg = PowderMeasurementConfig()
        
    sample_substrate_cfg = SampleSubstrateConfig(
        type=sample_substrate_type,
        shape=sample_substrate_shape,
        auto_detection=is_auto_substrate_detection,
        stub_w_mm=sample_substrate_width_mm
    )

    for sample in samples:
        # --- Sample configuration
        sample_ID = sample['ID']
        center_pos = sample['pos']
        
        print_double_separator()
        logging.info(f"Sample '{sample_ID}'")
        
        sample_cfg = SampleConfig(
            ID=sample_ID,
            elements=[],
            type='powder',
            center_pos=center_pos,
            half_width_mm=sample_halfwidth
        )
        
        # --- Run Composition Analyzer
        EM_analyzer = EMXSp_Composition_Analyzer(
            microscope_cfg=microscope_cfg,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            powder_meas_cfg=powder_meas_cfg,
            is_acquisition = True,
            development_mode=development_mode,
            output_filename_suffix = output_filename_suffix,
            verbose=verbose,
            results_dir=results_dir
        )
        
        
        try:
            EM_analyzer.EM_controller.particle_finder.get_particle_stats(n_par_target)
        except Exception as e:
            logging.exception(f"Sample '{sample_ID}': particle stats collection failed: {e}")
            continue
    
    # Put microscope in standby after completion
    if not development_mode and len(samples) > 1:
        try:
            EM_analyzer.EM_controller.standby()
        except Exception as e:
            logging.warning(f"Could not put microscope in standby: {e}")