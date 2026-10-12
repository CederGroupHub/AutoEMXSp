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

Created on Fri Aug 20 09:34:34 2025

@author: Andrea
"""

import logging
from pathlib import Path
from typing import List, Dict, Tuple, Any, Optional

from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
import autoemx.calibrations as calibs
import autoemx.config.defaults as dflts
import autoemx.utils.constants as cnst
from autoemx.utils import print_double_separator
from autoemx.config import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    QuantificationOptionsConfig,
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

__all__ = ["batch_acquire_experimental_stds"]

def batch_acquire_experimental_stds(
    stds: List[Dict[str, Any]],
    microscope_ID: str = dflts.microscope_ID,
    microscope_type: str = dflts.microscope_type,
    measurement_type: str = dflts.measurement_type,
    measurement_mode: str = dflts.measurement_mode,
    sample_halfwidth: float = dflts.sample_halfwidth_mm,
    sample_substrate_type: str = 'Ctape',
    sample_substrate_shape: str = 'circle',
    working_distance: float = dflts.working_distance_mm, #mm
    working_distance_tolerance: float = dflts.working_distance_tolerance_mm, #mm
    beam_energy: float = dflts.beam_energy_keV,
    spectrum_lims: Tuple[float, float] = dflts.spectrum_lims,
    use_instrument_background: bool = dflts.use_instrument_background,
    min_bckgrnd_cnts: Optional[float] = dflts.min_bckgrnd_cnts,
    fit_during_collection= True,
    update_std_library = True,
    is_auto_substrate_detection: bool = False,
    auto_adjust_brightness_contrast: bool = True,
    contrast: Optional[float] = None,
    brightness: Optional[float] = None,
    min_n_spectra: int = dflts.min_n_spectra,
    max_n_spectra: int = dflts.max_n_spectra,
    target_Xsp_counts: int = dflts.exp_std_target_Xsp_counts,
    max_XSp_acquisition_time: Optional[float] = None,
    els_substrate: Optional[List[str]] = None,
    powder_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    bulk_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    exp_stds_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    output_filename_suffix: str = '',
    saved_images_extension: str = dflts.saved_images_extension,
    save_raw_images: bool = dflts.save_raw_images,
    annotate_particle_images: bool = dflts.annotate_particle_images,
    development_mode: bool = False,
    verbose: bool = True,
    exp_std_dir: Optional[str] = None,
    is_manual_meas: bool = False,
) -> List[Optional[EMXSp_Composition_Analyzer]]:
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
    working_distance_tolerance : float, optional
        Defines maximum accepted deviation of working distance from its typical value, in mm.
            Used to prevent gross mistakes from EM autofocus. Default: 1 mm. 
    beam_energy : float, optional
        Electron beam energy in keV.  
        Default is `15.0` (MeasurementConfig.beam_energy_keV).
    spectrum_lims : tuple of float, optional
        Lower and upper energy limits for spectrum fitting in eV.  
        Default is `(14, 1100)` (QuantificationOptionsConfig.spectrum_lims).
    use_instrument_background : bool, optional
        Whether to use instrument background files during fitting.  
        If False, background is computed during fitting.  
        Default is `False` (QuantificationOptionsConfig.use_instrument_background).
    min_bckgrnd_cnts : float, optional
        Minimum background counts required for a spectrum not to be filtered out.  
        Default is `5` (ClusteringConfig.min_bckgrnd_cnts).
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
        Default is `['C', 'O', 'Al']` (defaults.substrate_els_without_identification: carbon tape on an Al stub; standards are fitted without identification of missing elements).
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
    saved_images_extension : str, optional
        Extension used for SEM image files saved during acquisition.
        Default is `'png'`.
    save_raw_images : bool, optional
        Whether to also save the non-annotated version of annotated SEM images (frame images, and
        particle images when ``annotate_particle_images`` is True). Default is `False`.
    annotate_particle_images : bool, optional
        Whether to draw the spectrum spots (dot + spectrum ID) and a scale bar on the saved particle images.
        If False, raw particle images are saved; annotated copies can be made afterwards with
        ``Annotate_Particle_Images.py``. Default is `False`.
    development_mode : bool, optional
        If True, enables development/debug features.  
        Default is `False`.
    verbose : bool, optional
        If True, print verbose output.  
        Default is `True`.
    exp_std_dir : str, optional
        Directory where experimental standard results are saved.  
        Default is `None`.
    is_manual_meas : bool, optional
        If True, all samples are measured in manual navigation mode unless overridden by 'is_manual_meas' in the sample dict. Default is False.

    Returns
    -------
    results : list(EMXSp_Composition_Analyzer)
        A list of the composition analysis objects (one per sample) containing the results and methods for further analysis.
    """
    if max_XSp_acquisition_time is None:
        max_XSp_acquisition_time = target_Xsp_counts / 10000 * 5
    if els_substrate is None:
        els_substrate = list(dflts.substrate_els_without_identification)
    
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

    quant_cfg = QuantificationOptionsConfig(
        spectrum_lims=spectrum_lims,
        use_instrument_background=use_instrument_background,
    )

    clustering_cfg = ClusteringConfig(min_bckgrnd_cnts=min_bckgrnd_cnts)

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
    
    results: List[Optional[EMXSp_Composition_Analyzer]] = []
    comp_analyzer: Optional[EMXSp_Composition_Analyzer] = None
    
    for std_sample in stds:
        # --- Sample configuration
        sample_ID = std_sample['ID']
        formula = std_sample['formula']
        center_pos = std_sample['pos']
        sample_type = std_sample['sample_type']
        is_manual_meas = std_sample.get('is_manual_meas', is_manual_meas)
        els_to_use_for_mean_PB_calc = std_sample.get('els_to_use_for_mean_PB_calc', None)

        print("\n \n")
        print_double_separator()
        print_double_separator(f"Sample '{sample_ID}'")

        exp_stds_kwargs = dict(exp_stds_meas_cfg_kwargs) if exp_stds_meas_cfg_kwargs else {}
        if els_to_use_for_mean_PB_calc is not None:
            exp_stds_kwargs['els_to_use_for_mean_PB_calc'] = els_to_use_for_mean_PB_calc
        exp_stds_cfg = ExpStandardsConfig(is_exp_std_measurement=True, formula=formula, **exp_stds_kwargs)
        measurement_cfg = MeasurementConfig(
            type=measurement_type,
            mode=measurement_mode,
            working_distance = working_distance,
            working_distance_tolerance = working_distance_tolerance,
            beam_energy_keV=beam_energy,
            is_manual_navigation=is_manual_meas,
            max_acquisition_time=max_XSp_acquisition_time,
            target_acquisition_counts=target_Xsp_counts,
            min_n_spectra=min_n_spectra,
            max_n_spectra=max_n_spectra,
            saved_images_extension=saved_images_extension,
            save_raw_images=save_raw_images,
            annotate_particle_images=annotate_particle_images,
        )
        
        
        # The composition of the standard, from its formula, is kept fixed during fits
        w_frs = exp_stds_cfg.formula_w_frs()
        if w_frs is None:
            raise ValueError(
                f"Experimental standards configuration for '{sample_ID}' is missing the formula of the standard."
            )
        elements = list(w_frs.keys())
        sample_cfg = SampleConfig(
            elements=elements,
            w_frs=w_frs,
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
            sample_id=sample_ID,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            quant_cfg=quant_cfg,
            initial_clustering_cfg=clustering_cfg,
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
    if comp_analyzer is not None and not development_mode and len(stds) > 1:
        try:
            comp_analyzer.EM_controller.standby()
        except Exception as e:
            logging.warning(f"Could not put microscope in standby: {e}")

    # Generate a standards-coverage report only when standards library was updated.
    if update_std_library:
        try:
            from autoemx.runners.extract_experimental_standards_details import extract_experimental_standards_details  # type: ignore

            if exp_std_dir is not None:
                standards_dir = Path(exp_std_dir).expanduser().resolve()
            else:
                standards_dir = Path(__file__).resolve().parents[1] / cnst.CALIBS_DIR / microscope_ID

            standards_json_filename = f"{measurement_type}_{cnst.STD_FILENAME}_{int(beam_energy):d}keV.json"
            standards_json_path = standards_dir / standards_json_filename

            report_kwargs = {
                "microscope_ID": microscope_ID,
                "voltage": beam_energy,
                "report_output_dir": str(standards_dir),
                "print_report": True,
                "print_per_peak": False,
            }
            if standards_json_path.exists():
                report_kwargs["standards_json_path"] = str(standards_json_path)

            report_path = extract_experimental_standards_details(**report_kwargs)
            logging.info("Standards details report generated at: %s", report_path)
        except Exception as e:
            logging.warning("Could not generate standards details report after collection: %s", e)
    
    return results
    
    
    