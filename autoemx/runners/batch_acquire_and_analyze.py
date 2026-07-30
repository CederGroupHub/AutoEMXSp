#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Automated X-Ray Spectral Acquisition and Analysis

This module configures and runs automated collection and (optionally) quantification
of EDS/WDS spectra for powder samples using an electron microscope (EM) with 
a specified substrate and calibration.

Import this module in your own Python code and call the
`batch_acquire_and_analyze()` function, passing your desired configuration
and sample list as arguments. This enables integration into larger
automation workflows or pipelines.

Workflow includes:
    - Configuration of microscope, measurement, substrate, and fitting parameters
    - Sample setup (elements, position, reference formulae)
    - Calibration and quantification options
    - Optional clustering for compositional analysis
    - Automated or manual navigation and acquisition modes

Requirements:
    - Proper instrument calibration files and instrument driver for the selected microscope

Typical usage:
    - Edit the 'samples' list to define your standards or unknowns
    - Adjust configuration parameters as needed, either directly in the script or by passing them to `batch_acquire_and_analyze`
    - Run the script, or import and call `batch_acquire_and_analyze()` to perform spectrum collection and (optionally) quantification for one or multiple samples at a time

Created on Fri Jul 26 09:34:34 2024

@author: Andrea
"""

import logging
from typing import List, Dict, Tuple, Any, Optional, Sequence

from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
from autoemx.core.em_runtime.xsp_spot_selection import XSpSpotSelectorCallback
import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.utils import print_double_separator
import autoemx.config.defaults as dflt
from autoemx.config import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    QuantificationOptionsConfig,
    ClusteringConfig,
    PowderMeasurementConfig,
    BulkMeasurementConfig,
    PlotConfig,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["batch_acquire_and_analyze"]

def batch_acquire_and_analyze(
    samples: List[Dict[str, Any]],
    microscope_ID: str = dflt.microscope_ID,
    microscope_type: str = dflt.microscope_type,
    measurement_type: str = dflt.measurement_type,
    measurement_mode: str = dflt.measurement_mode,
    quantification_method: str = dflt.quantification_method,
    sample_type: str = 'powder',
    sample_halfwidth: float = 3.0,
    sample_substrate_type: str = 'Ctape',
    sample_substrate_shape: str = 'circle',
    sample_substrate_width_mm: float = 12,
    working_distance: float = 5, #mm
    working_distance_tolerance: float = 1, #mm
    beam_energy: float = 15.0,
    spectrum_lims: Tuple[float, float] = dflt.spectrum_lims,
    use_instrument_background: bool = dflt.use_instrument_background,
    use_project_specific_std_dict: bool = False,
    interrupt_fits_bad_spectra: bool = True,
    max_analytical_error_percent: float = 5,
    min_bckgrnd_cnts: float = 5,
    quant_flags_accepted: Optional[List[int]] = None,
    max_n_clusters: int = 6,
    show_unused_comps_clust: bool = True,
    is_manual_navigation: bool = False,
    is_auto_substrate_detection: bool = False,
    auto_adjust_brightness_contrast: bool = True,
    contrast: Optional[float] = None,
    brightness: Optional[float] = None,
    quantify_spectra: bool = False,
    min_n_spectra: int = 50,
    max_n_spectra: int = 100,
    target_Xsp_counts: int = 50000,
    max_XSp_acquisition_time: Optional[float] = None,
    els_substrate: Optional[List[str]] = None,
    powder_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    bulk_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    standards_dict: Optional[Dict[str, float]] = None,
    output_filename_suffix: str = '',
    saved_images_extension: str = dflt.saved_images_extension,
    save_raw_images: bool = dflt.save_raw_images,
    development_mode: bool = False,
    verbose: bool = True,
    results_dir: Optional[str] = None,
    xsp_spot_selector: Optional[XSpSpotSelectorCallback] = None,
    particle_ids: Optional[Sequence[int]] = None,
    particle_coords_mm: Optional[Sequence[Sequence[float]]] = None,
) -> EMXSp_Composition_Analyzer:
    """
    Batch acquisition (and optional quantification) of X-ray spectra for a list of powder samples.

    Parameters
    ----------
    samples : list of dict
        List of sample definitions. Each dictionary must contain:
            - 'ID' (str): Sample identifier (SampleConfig.ID).
            - 'els' (list of str): List of expected element symbols (SampleConfig.elements).
            - 'pos' (tuple of float): (x, y) stage coordinates in mm (SampleConfig.center_pos).
            - 'cnd' (list of str, optional): Reference chemical formulae for known phases (ClusteringConfig.ref_formulae).
            - 'particle_ids' (list of int, optional): Per-sample override for list-mode
              ledger particle ids (requires ``par_selection_mode='list'``).
            - 'particle_coords_mm' (list of (x, y), optional): Per-sample override for
              list-mode absolute stage coordinates in mm.
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
    quantification_method : str, optional
        Quantification method. Currently only `'PB'` (Phi-Rho-Z) is implemented.
        Default is `'PB'` (QuantificationOptionsConfig.method).
    sample_type : str, optional
        Sample type. Allowed: `'powder'` (implemented), `'bulk'`, `'film'` (not implemented).
        Default is `'powder'` (SampleConfig.type).
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
    use_project_specific_std_dict : bool, optional
        If True, loads standards from project folder (i.e. results_dir) during quantification.
        Default: False
    interrupt_fits_bad_spectra : bool, optional
        If True, fitting stops early for poor-quality spectra.
        Default is `True`.
    max_analytical_error_percent : float, optional
        Maximum allowed analytical error (%) for compositions to be included in clustering.
        Default is `5` (ClusteringConfig.max_analytical_error_percent).
    min_bckgrnd_cnts : float, optional
        Minimum background counts required for a spectrum not to be filtered out.
        Default is `5` (ClusteringConfig.min_bckgrnd_cnts).
    quant_flags_accepted : list of int, optional
        List of acceptable quantification flags; others are filtered out before clustering.
        Default is `[0, -1]` (ClusteringConfig.quant_flags_accepted).
    max_n_clusters : int, optional
        Maximum number of clusters allowed in compositional clustering.
        Default is `6` (ClusteringConfig.max_k).
    show_unused_comps_clust : bool, optional
        Whether to display unused compositions in clustering plots.
        Default is `True` (PlotConfig.show_unused_comps_clust).
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
    quantify_spectra : bool, optional
        If True, perform quantification after acquisition.
        Default is `False`.
    min_n_spectra : int, optional
        Minimum number of spectra to acquire.
        Default is `50` (MeasurementConfig.min_n_spectra).
    max_n_spectra : int, optional
        Maximum number of spectra to acquire.
        Default is `100` (MeasurementConfig.max_n_spectra).
    target_Xsp_counts : int, optional
        Target counts for spectrum acquisition.
        Default is `50000` (MeasurementConfig.target_acquisition_counts).
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
    standards_dict : dict, optional
        Dictionary of reference PB values from experimental standards. Default : None.
        If None, dictionary of standards is loaded from the XSp_calibs/Your_Microscope_ID directory.
        Provide standards_dict only when providing different standards from those normally used for quantification.
    output_filename_suffix : str, optional
        String appended to output filenames.
        Default is `''`.
    saved_images_extension : str, optional
        Extension used for SEM image files saved during acquisition.
        Default is `'tif'`.
    save_raw_images : bool, optional
        Whether to save non-annotated raw SEM images in addition to annotated images.
        Default is `True`.
    development_mode : bool, optional
        If True, enables development/debug features.
        Default is `False`.
    verbose : bool, optional
        If True, print verbose output.
        Default is `True`.
    results_dir : str, optional
        Directory where results are saved.
        Default is `None`.
    xsp_spot_selector : callable, optional
        Callback invoked on each particle when
        ``powder_meas_cfg.par_spot_selection_mode='callback'``. Called repeatedly on
        the same particle until it returns ``[]``. Signature::

            spots = xsp_spot_selector(context: XSpSpotSelectionContext)

        Return pixel ``(x, y)`` coordinates in the current frame. Return ``[]``
        to finish with the current particle (including to skip it). Spots are
        validated against image bounds only. Stage navigation to the next
        particle remains automatic; spots are measured via beam positioning
        with drift correction like automated mode.
    particle_ids : sequence of int, optional
        Global list of ledger particle ids to visit when
        ``powder_meas_cfg.par_selection_mode='list'``. Mutually exclusive with
        ``particle_coords_mm``. Can be overridden per sample via
        ``sample['particle_ids']``. Missing ids raise before acquisition starts.
    particle_coords_mm : sequence of (x, y), optional
        Global list of absolute stage coordinates in mm to visit when
        ``par_selection_mode='list'``. Mutually exclusive with ``particle_ids``.
        Can be overridden per sample via ``sample['particle_coords_mm']``.
    
            
    Returns
    -------
    comp_analyzer : EMXSp_Composition_Analyzer
        The composition analysis object containing the results and methods for further analysis.
    """
    if max_XSp_acquisition_time is None:
        max_XSp_acquisition_time = target_Xsp_counts / 10000 * 5
    if els_substrate is None:
        els_substrate = ['C', 'O', 'Al']
    if quant_flags_accepted is None:
        quant_flags_accepted = [0, -1]

    resolved_spot_selector = xsp_spot_selector
    if powder_meas_cfg_kwargs and powder_meas_cfg_kwargs.get('par_spot_selection_mode') == 'callback':
        if resolved_spot_selector is None:
            raise ValueError(
                "xsp_spot_selector is required when powder_meas_cfg_kwargs sets "
                "par_spot_selection_mode='callback'."
            )
    if particle_ids is not None and particle_coords_mm is not None:
        raise ValueError(
            "Provide exactly one of particle_ids or particle_coords_mm "
            "(or pass them per-sample)."
        )
    
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

    measurement_cfg = MeasurementConfig(
        type=measurement_type,
        mode=measurement_mode,
        working_distance = working_distance,
        working_distance_tolerance = working_distance_tolerance,
        beam_energy_keV=beam_energy,
        is_manual_navigation=is_manual_navigation,
        max_acquisition_time=max_XSp_acquisition_time,
        target_acquisition_counts=target_Xsp_counts,
        min_n_spectra=min_n_spectra,
        max_n_spectra=max_n_spectra,
        saved_images_extension=saved_images_extension,
        save_raw_images=save_raw_images
    )

    quant_cfg = QuantificationOptionsConfig(
        method = quantification_method,
        spectrum_lims=spectrum_lims,
        use_instrument_background=use_instrument_background,
        use_project_specific_std_dict = use_project_specific_std_dict
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
        auto_detection=is_auto_substrate_detection,
        stub_w_mm=sample_substrate_width_mm
    )

    comp_analyzer: Optional[EMXSp_Composition_Analyzer] = None

    for sample in samples:
        # --- Sample configuration
        sample_ID = sample['ID']
        elements = sample['els']
        center_pos = sample['pos']
        ref_formulae = sample.get('cnd', [])
        smpl_type = sample.get('type', sample_type)
        sample_spot_selector = sample.get('xsp_spot_selector', resolved_spot_selector)
        if powder_meas_cfg.par_spot_selection_mode == 'callback' and sample_spot_selector is None:
            raise ValueError(
                f"Sample '{sample_ID}': xsp_spot_selector is required when "
                "par_spot_selection_mode='callback'."
            )

        sample_particle_ids = sample.get('particle_ids', particle_ids)
        sample_particle_coords_mm = sample.get('particle_coords_mm', particle_coords_mm)
        if sample_particle_ids is not None and sample_particle_coords_mm is not None:
            raise ValueError(
                f"Sample '{sample_ID}': provide exactly one of particle_ids or "
                "particle_coords_mm."
            )
        if powder_meas_cfg.par_selection_mode == 'list':
            if sample_particle_ids is None and sample_particle_coords_mm is None:
                raise ValueError(
                    f"Sample '{sample_ID}': par_selection_mode='list' requires "
                    "particle_ids or particle_coords_mm."
                )
        elif sample_particle_ids is not None or sample_particle_coords_mm is not None:
            raise ValueError(
                f"Sample '{sample_ID}': particle_ids / particle_coords_mm require "
                "powder_meas_cfg.par_selection_mode='list'."
            )
        
        print_double_separator()
        logging.info(f"Sample '{sample_ID}'")
        
        sample_cfg = SampleConfig(
            elements=elements,
            type=smpl_type,
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

        # --- Clustering configuration
        if any(el in getattr(calibs, 'undetectable_els', []) for el in elements) or len(elements) == 2:
            clustering_features = cnst.W_FR_CL_FEAT
        else:
            clustering_features = cnst.AT_FR_CL_FEAT
        k_finding_method = 'silhouette' if len(elements) > 2 else 'calinski_harabasz'

        clustering_cfg = ClusteringConfig(
            method='kmeans',
            features=clustering_features,
            k_finding_method=k_finding_method,
            max_k=max_n_clusters,
            ref_formulae=ref_formulae,
            max_analytical_error_percent=max_analytical_error_percent,
            min_bckgrnd_cnts=min_bckgrnd_cnts,
            quant_flags_accepted=quant_flags_accepted
        )

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
            plot_cfg=PlotConfig(),
            is_acquisition=True,
            development_mode=development_mode,
            standards_dict=standards_dict,
            output_filename_suffix=output_filename_suffix,
            verbose=verbose,
            results_dir=results_dir,
            xsp_spot_selector=sample_spot_selector,
            particle_ids=sample_particle_ids,
            particle_coords_mm=sample_particle_coords_mm,
        )
        
        try:
            comp_analyzer.run_collection_and_quantification(quantify=quantify_spectra, interrupt_fits_bad_spectra=interrupt_fits_bad_spectra)
        except Exception as e:
            logging.exception(f"Sample '{sample_ID}': acquisition/quantification failed: {e}")
            continue
    
    # Put microscope in standby after completion
    if (
        comp_analyzer is not None
        and getattr(comp_analyzer, "EM_controller", None) is not None
        and not development_mode
        and len(samples) > 1
    ):
        try:
            comp_analyzer.EM_controller.standby()
        except Exception as e:
            logging.warning(f"Could not put microscope in standby: {e}")
    elif comp_analyzer is None:
        logging.warning("No samples were processed; skipping microscope standby.")
    
    if comp_analyzer is None:
        raise ValueError("No samples were successfully processed. comp_analyzer was not initialized.")
    
    return comp_analyzer