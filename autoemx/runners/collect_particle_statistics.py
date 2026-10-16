#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 10:14:12 2025

@author: Andrea
"""

import logging
from typing import List, Dict, Any, Optional

from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
from autoemx.utils import print_double_separator
import autoemx.config.defaults as dflt
from autoemx.config import (
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

__all__ = ["collect_particle_statistics"]

def collect_particle_statistics(
    samples: List[Dict[str, Any]],
    n_par_target: int,
    microscope_ID: str = dflt.microscope_ID,
    microscope_type: str = dflt.microscope_type,
    detector_type: str = dflt.detector_type,
    sample_halfwidth: float = dflt.sample_halfwidth_mm,
    sample_substrate_type: str = 'Ctape',
    sample_substrate_shape: str = 'circle',
    sample_substrate_width_mm: float = dflt.substrate_width_mm,
    working_distance: float = dflt.working_distance_mm, #mm
    is_manual_navigation: bool = False,
    is_auto_substrate_detection: bool = False,
    auto_adjust_brightness_contrast: bool = True,
    contrast: float = dflt.particle_stats_contrast,
    brightness: float = dflt.particle_stats_brightness,
    saved_images_extension: str = dflt.saved_images_extension,
    save_raw_images: bool = dflt.save_raw_images,
    powder_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    output_filename_suffix: str = '',
    development_mode: bool = False,
    verbose: bool = True,
    results_dir: Optional[str] = None
) -> None:
    """
    Particle size statistics of a list of powder samples, measured at the electron microscope.

    For each sample, particles are detected in frames scanned over the sample until ``n_par_target`` particles are
    found in each size band (see ``EM_Particle_Finder.get_particle_stats``). The statistics are stored in the
    sample ledger (``ledger.json``) and saved as particle sizes, statistics and histogram files in the sample folder.

    Parameters
    ----------
    samples : list of dict
        List of sample definitions. Each dictionary must contain:
            - 'ID' (str): Sample identifier (SampleConfig.ID).
            - 'pos' (tuple of float): (x, y) stage coordinates in mm (SampleConfig.center_pos).
    n_par_target : int
        Target number of particles to analyse, in each size band.
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
    saved_images_extension : str, optional
        File extension for saved SEM images. Allowed: `'png'`, `'tif'`, `'jpg'`, `'webp'`.
        `'png'` and `'tif'` are both lossless (8-bit), `'png'` smaller; avoid `'jpg'` for images to reanalyse.
        Default is taken from defaults (MeasurementConfig.saved_images_extension).
    save_raw_images : bool, optional
        If True, also save raw (unannotated) SEM images alongside annotated frames.
        Default is `False` (saves only annotated images, lightweight output).
        Default is taken from defaults (MeasurementConfig.save_raw_images).
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
    None
        Results are saved in the folder of each sample.
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

    if powder_meas_cfg_kwargs:
        powder_meas_cfg = PowderMeasurementConfig(**powder_meas_cfg_kwargs)
    else:
        powder_meas_cfg = PowderMeasurementConfig()
    
    # Build MeasurementConfig with image extension and raw image options
    measurement_cfg = MeasurementConfig(
        type=MeasurementConfig.PARTICLE_STATS_MEAS_TYPE_KEY,
        working_distance=working_distance,
        is_manual_navigation=is_manual_navigation,
        powder_meas_cfg=powder_meas_cfg,
        saved_images_extension=saved_images_extension,
        save_raw_images=save_raw_images,
    )
        
    sample_substrate_cfg = SampleSubstrateConfig(
        type=sample_substrate_type,
        shape=sample_substrate_shape,
        auto_detection=is_auto_substrate_detection,
        stub_w_mm=sample_substrate_width_mm
    )
    
    EM_analyzer = None

    for sample in samples:
        # --- Sample configuration
        sample_ID = sample['ID']
        center_pos = sample['pos']
        
        print_double_separator(f"Sample '{sample_ID}'")
        
        sample_cfg = SampleConfig(
            elements=[],
            type='powder',
            center_pos=center_pos,
            half_width_mm=sample_halfwidth
        )
        
        # --- Run Composition Analyzer
        EM_analyzer = EMXSp_Composition_Analyzer(
            microscope_cfg=microscope_cfg,
            sample_id=sample_ID,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            is_acquisition = True,
            development_mode=development_mode,
            output_filename_suffix = output_filename_suffix,
            verbose=verbose,
            results_dir=results_dir
        )
        
        
        try:
            if EM_analyzer and EM_analyzer.EM_controller and EM_analyzer.EM_controller.particle_finder:
                EM_analyzer.EM_controller.particle_finder.get_particle_stats(n_par_target)
        except Exception as e:
            logging.exception(f"Sample '{sample_ID}': particle stats collection failed: {e}")
            continue
    
    # Put microscope in standby after completion
    if not development_mode and EM_analyzer and len(samples) > 1:
        try:
            EM_analyzer.EM_controller.standby()
        except Exception as e:
            logging.warning(f"Could not put microscope in standby: {e}")