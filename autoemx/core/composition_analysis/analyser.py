#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EMXSp_Composition_Analyzer

Main class for automated compositional analysis of electron microscopy X-ray spectroscopy (EMXSp) data.

Can be run from run_acquisition_quant_analysis.py

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
            sample_id='sample_001',
            microscope_cfg=microscope_cfg,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            quant_cfg=quant_cfg,
            initial_clustering_cfg=initial_clustering_cfg,
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
import re
import time
import shutil
import itertools
import warnings
import traceback
from pathlib import Path
from datetime import datetime
from dataclasses import asdict
from typing import Any, Optional, Tuple, List, Dict, Iterable, Sequence, Union
from joblib import Parallel, delayed

#TODO remove in future versions
warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    module="cvxpy.*",
)

# Third-party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.core.composition import Composition
from sklearn.cluster import KMeans
import cvxpy as cp
 

# Project-specific imports
from autoemx.core.quantifier import XSp_Quantifier
from autoemx.core.em_runtime.controller import EM_Controller
from autoemx.core.em_runtime.sample_finder import EM_Sample_Finder
import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.utils.helper import (
    print_single_separator,
    print_double_separator,
    make_unique_path,
)
from autoemx.config import (
    MicroscopeConfig,
    SampleConfig,
    MeasurementConfig,
    SampleSubstrateConfig,
    QuantificationOptionsConfig,
    PowderMeasurementConfig,
    BulkMeasurementConfig,
    ExpStandardsConfig,
    PlotConfig,
)
from autoemx.config.ledger_schemas import (

    AcquisitionDetails,
    ClusteringConfig as LedgerClusteringConfig,
    Coordinate2D,
    LedgerConfigs,
    QuantificationConfig,
    QuantificationDiagnostics,
    QuantificationResult,
    SampleLedger,
    SpotCoordinates,
    SpectrumEntry,
)
from autoemx.config.ledger_io import ingest_spectra
from autoemx.config.schema_models.clustering import ClusteringAnalysis, ClusteringResult
from autoemx.config.schema_models.standards import EDSStandardsFile
from .clustering import ClusteringModule
from autoemx.utils.legacy.ledger_bootstrap import (
    load_or_create_ledger_with_legacy_data_csv,
    resolve_legacy_bootstrap_configs,
)
from autoemx.utils.legacy.legacy_config_loader import has_legacy_data_csv
from .plotting import PlottingModule
from .reference_matching import ReferenceMatchingModule
from autoemx.utils.legacy.spectrum_pointer_writer import load_vendor_msa_template_lines, write_spectrum_pointer_file
from .standards import StandardsModule

from autoemx._logging import get_logger
logger = get_logger(__name__)

parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def _resolve_parallel_mode(requested_cores: int) -> tuple[int, Optional[str], bool]:
    """Resolve execution mode shared by quantify and fitting-only workflows.

    Returns
    -------
    n_cores : int
        Effective number of cores to use.
    backend : Optional[str]
        Joblib backend to use when parallel is enabled.
    use_parallel : bool
        Whether parallel execution is enabled.
    """
    n_cores = max(1, int(requested_cores))
    backend = 'loky'
    use_parallel = n_cores > 1
    if not use_parallel:
        n_cores = 1
        backend = None
    return n_cores, backend, use_parallel


def _worker_is_spectrum_valid_for_fitting(
    *,
    spectrum: np.ndarray,
    energy_vals: Sequence[float],
    sp_start: int,
    sp_end: int,
    target_acquisition_counts: int,
    min_bckgrnd_cnts: Optional[float],
) -> tuple[bool, Optional[int], Optional[str]]:
    """Pure worker-safe version of spectrum pre-fit validation."""
    if spectrum is None:
        return False, 1, "No spectral data present"

    if np.sum(spectrum) < 0.9 * target_acquisition_counts:
        return False, 2, "Total counts too low"

    n_vals_considered = 20
    filter_len = 3
    en_threshold = 2

    xy_data = zip(energy_vals, spectrum[sp_start:sp_end])
    spectrum_data_to_consider = [cnts for en, cnts in xy_data if cnts > 0 and en < en_threshold]
    if len(spectrum_data_to_consider) == 0:
        return False, 3, "Background counts too low"

    spectrum_smooth = np.convolve(
        spectrum_data_to_consider,
        np.ones(filter_len) / filter_len,
        mode='same',
    )
    min_vals = np.sort(spectrum_smooth)[:n_vals_considered]
    if min_bckgrnd_cnts is not None and len(min_vals) > 0 and all(min_vals < min_bckgrnd_cnts):
        return False, 3, "Background counts too low"

    return True, None, None


def _worker_flag_spectrum_for_clustering(
    *,
    min_bckgrnd_ref_lines: float,
    quantifier: Any,
    detectable_els_substrate: Sequence[str],
    min_bckgrnd_cnts: Optional[float],
) -> tuple[bool, int, str]:
    """Pure worker-safe version of clustering-validity flagging."""
    is_spectrum_valid = True
    sub_peak_int_threshold = 10
    sub_peak_int_thresh_cnts = quantifier.tot_sp_counts * sub_peak_int_threshold / 100

    els_substrate_intensities = {el: 0 for el in detectable_els_substrate}
    for el_line, peak_info in quantifier.fitted_peaks_info.items():
        el = el_line.split('_')[0]
        if el in detectable_els_substrate:
            els_substrate_intensities[el] += peak_info[cnst.PEAK_INTENSITY_KEY]

    for el, peak_int in els_substrate_intensities.items():
        if peak_int > sub_peak_int_thresh_cnts:
            return False, 7, f"{el} {peak_int:.0f} counts > {sub_peak_int_threshold} % of total counts"

    if min_bckgrnd_cnts is not None and min_bckgrnd_ref_lines < min_bckgrnd_cnts:
        return False, 8, f"Reference background counts too low ({min_bckgrnd_ref_lines:.1f} < {min_bckgrnd_cnts})"

    return True, 0, ""


def _worker_check_fit_quant_validity(
    *,
    is_quant_fit_valid: bool,
    bad_quant_flag: Optional[int],
    quantifier: Any,
    min_bckgrnd_ref_lines: Any,
    detectable_els_substrate: Sequence[str],
    min_bckgrnd_cnts: Optional[float],
) -> tuple[int, str]:
    """Pure worker-safe version of post-fit quantification validation."""
    start_str_comments = 'Fit interrupted due to ' if not is_quant_fit_valid else ''

    if bad_quant_flag == 1:
        comment = start_str_comments + 'poor fit'
        quant_flag = 4
    elif bad_quant_flag == 2:
        comment = start_str_comments + 'excessively high analytical error'
        quant_flag = 5
    elif bad_quant_flag == 3:
        comment = start_str_comments + 'excessive X-ray absorption'
        quant_flag = 6
    elif not is_quant_fit_valid:
        comment = 'Fit interrupted for unknown reasons'
        quant_flag = 9
    else:
        _, quant_flag, comment = _worker_flag_spectrum_for_clustering(
            min_bckgrnd_ref_lines=min_bckgrnd_ref_lines,
            quantifier=quantifier,
            detectable_els_substrate=detectable_els_substrate,
            min_bckgrnd_cnts=min_bckgrnd_cnts,
        )

    if bad_quant_flag == -1 and quant_flag == 0:
        comment = f"{comment} - Quantification did not converge." if comment else 'Quantification did not converge.'
        quant_flag = -1

    return quant_flag, comment


def _quantify_spectrum_worker(worker_payload: Dict[str, Any]) -> tuple[int, Optional[Dict[str, Any]], QuantificationResult, Optional[float]]:
    """Process-safe worker for one spectrum quantification or fitting-only."""
    spectrum_index = int(worker_payload['spectrum_index'])
    pointer_abs = Path(worker_payload['pointer_abs'])
    counts = SampleLedger._load_counts_from_pointer_file(pointer_abs)
    spectrum = np.asarray(counts, dtype=float)

    background = None
    background_abs = worker_payload.get('background_abs')
    if background_abs:
        background = EMXSp_Composition_Analyzer._load_background_vector_from_file(Path(background_abs))

    start_quant_time = time.time()
    is_sp_valid_for_fitting, quant_flag, comment = _worker_is_spectrum_valid_for_fitting(
        spectrum=spectrum,
        energy_vals=worker_payload['energy_vals'],
        sp_start=int(worker_payload['sp_start']),
        sp_end=int(worker_payload['sp_end']),
        target_acquisition_counts=int(worker_payload['target_acquisition_counts']),
        min_bckgrnd_cnts=worker_payload['min_bckgrnd_cnts'],
    )
    if not is_sp_valid_for_fitting:
        quant_record = QuantificationResult(
            quantification_id=int(worker_payload['quantification_id']),
            quant_flag=quant_flag,
            comment=comment,
            diagnostics=QuantificationDiagnostics(
                converged=False,
                interrupted=True,
            ),
        )
        return spectrum_index, None, quant_record, time.time() - start_quant_time

    quantifier = XSp_Quantifier(
        spectrum_vals=spectrum,
        spectrum_lims=(int(worker_payload['sp_start']), int(worker_payload['sp_end'])),
        microscope_ID=worker_payload['microscope_id'],
        meas_type=worker_payload['measurement_type'],
        meas_mode=worker_payload['measurement_mode'],
        det_ch_offset=float(worker_payload['det_ch_offset']),
        det_ch_width=float(worker_payload['det_ch_width']),
        beam_e=float(worker_payload['beam_energy_keV']),
        emergence_angle=float(worker_payload['emergence_angle']),
        energy_vals=None,
        background_vals=background,
        els_sample=worker_payload['all_els_sample'],
        els_substrate=worker_payload['detectable_els_substrate'],
        els_w_fr=worker_payload.get('sample_w_frs'),  # None for fitting-only
        is_particle=bool(worker_payload['apply_geom_factors']),
        sp_collection_time=float(worker_payload['sp_collection_time']),
        max_undetectable_w_fr=float(worker_payload['max_undetectable_w_fr']),
        fit_tol=float(worker_payload['fit_tolerance']),
        standards_dict=worker_payload['standards_dict'],
        verbose=False,
        fitting_verbose=False,
    )

    try:
        quantify = bool(worker_payload.get('quantify', True))
        if quantify:
            quant_result, min_bckgrnd_ref_lines, bad_quant_flag = quantifier.quantify_spectrum(
                print_result=False,
                interrupt_fits_bad_spectra=bool(worker_payload['interrupt_fits_bad_spectra']),
            )
            is_quant_fit_valid = quant_result is not None
        else:
            # Fitting-only path (experimental standards)
            bad_quant_flag = quantifier.initialize_and_fit_spectrum(print_results=False)
            quant_result = None
            is_quant_fit_valid = True
            min_bckgrnd_ref_lines = quantifier._get_min_bckgrnd_cnts_ref_quant_lines()
    except Exception as exc:
        try:
            logger.error(f"❌ {type(exc).__name__}: {exc}")
        except BaseException:
            pass
        quant_flag, comment = _worker_check_fit_quant_validity(
            is_quant_fit_valid=False,
            bad_quant_flag=None,
            quantifier=None,
            min_bckgrnd_ref_lines=None,
            detectable_els_substrate=worker_payload['detectable_els_substrate'],
            min_bckgrnd_cnts=worker_payload['min_bckgrnd_cnts'],
        )
        quant_record = quantifier.export_quantification_result(
            quantification_id=int(worker_payload['quantification_id']),
            quant_result=None,
            quant_flag=quant_flag,
            comment=comment,
        )
        return spectrum_index, None, quant_record, time.time() - start_quant_time

    quant_flag, comment = _worker_check_fit_quant_validity(
        is_quant_fit_valid=is_quant_fit_valid,
        bad_quant_flag=bad_quant_flag,
        quantifier=quantifier,
        min_bckgrnd_ref_lines=min_bckgrnd_ref_lines,
        detectable_els_substrate=worker_payload['detectable_els_substrate'],
        min_bckgrnd_cnts=worker_payload['min_bckgrnd_cnts'],
    )
    quant_record = quantifier.export_quantification_result(
        quantification_id=int(worker_payload['quantification_id']),
        quant_result=quant_result,
        quant_flag=quant_flag,
        comment=comment,
    )
    return spectrum_index, quant_result, quant_record, time.time() - start_quant_time

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
        Configuration for sample composition and geometry.
    measurement_cfg : MeasurementConfig
        Configuration for the measurement/acquisition.
    sample_substrate_cfg : SampleSubstrateConfig
        Configuration for the sample substrate.
    quant_cfg : QuantificationOptionsConfig
        Configuration for spectrum fitting and quantification.
    initial_clustering_cfg : autoemx.config.ledger_schemas.ClusteringConfig
        Initial clustering settings used when no active quantification config is yet available in the ledger.
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
    standards_dict : autoemx.config.schema_models.EDSStandardsFile | dict, optional
        Standards payload for reference PB values from experimental standards.
        If a dict is provided, it must conform to the standards Pydantic schema (or legacy
        shape that can be normalized by the schema compatibility adapter).
        If None, standards are loaded from the calibration directory for the active microscope.
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
    sample_id : Optional[str], optional
        Canonical identifier of the analyzed sample/run. Used for persistence
        (e.g., ledger.sample_id), output directory naming, and file prefixes.
        If omitted, it is derived from `results_dir` basename.

    Attributes
    ----------
    TO COMPLETE
    """
    #TODO
    @staticmethod
    def _coerce_optional_finite_float(value: Any) -> Optional[float]:
        """Return a finite float or None when the input is missing/non-finite."""
        if value is None:
            return None
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            return None
        if not np.isfinite(numeric_value):
            return None
        return numeric_value

    @staticmethod
    def _to_finite_float(value: Any, default: float = 0.0) -> float:
        """Return a finite float, falling back to default for NaN/inf/unparseable values."""
        try:
            numeric_value = float(value)
        except (TypeError, ValueError):
            return float(default)
        if not np.isfinite(numeric_value):
            return float(default)
        return numeric_value

    @staticmethod
    def _to_finite_float_list(values: Any) -> List[float]:
        """Convert a 1D numeric sequence to finite floats (NaN/inf -> 0.0)."""
        arr = np.asarray(values, dtype=float)
        return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0).tolist()

    @staticmethod
    def _to_finite_float_matrix(values: Any) -> List[List[float]]:
        """Convert a 2D numeric sequence to finite floats (NaN/inf -> 0.0)."""
        arr = np.asarray(values, dtype=float)
        return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0).tolist()

    @classmethod
    def _sanitize_json_compatible(cls, value: Any) -> Any:
        """Recursively replace non-finite numeric values with None in free-form payloads."""
        if isinstance(value, dict):
            return {k: cls._sanitize_json_compatible(v) for k, v in value.items()}
        if isinstance(value, list):
            return [cls._sanitize_json_compatible(v) for v in value]
        if isinstance(value, tuple):
            return [cls._sanitize_json_compatible(v) for v in value]
        if isinstance(value, np.generic):
            return cls._sanitize_json_compatible(value.item())
        if isinstance(value, (float, np.floating)):
            return None if not np.isfinite(float(value)) else float(value)
        return value

    @staticmethod
    def _clean_sample_id(sample_id: str) -> str:
        """Normalize sample_id for folder/file-safe persistence."""
        cleaned_id = str(sample_id).rstrip()
        cleaned_id = re.sub(r'^\W+|\W+$', '', cleaned_id)
        if not cleaned_id:
            raise ValueError("sample_id must contain at least one alphanumeric character.")
        return cleaned_id

    @classmethod
    def _resolve_sample_id(
        cls,
        sample_id: Optional[str],
        results_dir: Optional[str],
        is_acquisition: bool,
    ) -> str:
        """Resolve canonical sample_id from explicit input or a sample folder path."""
        if sample_id is not None and str(sample_id).strip():
            return cls._clean_sample_id(sample_id)

        if results_dir:
            derived_id = os.path.basename(os.path.normpath(results_dir))
            if derived_id:
                return cls._clean_sample_id(derived_id)

        mode_hint = "acquisition" if is_acquisition else "analysis"
        raise ValueError(
            f"sample_id is required for {mode_hint} unless it can be derived from results_dir."
        )

    def __init__(
        self,
        microscope_cfg: MicroscopeConfig,
        sample_cfg: SampleConfig,
        measurement_cfg: MeasurementConfig,
        sample_substrate_cfg: SampleSubstrateConfig,
        quant_cfg: Optional[QuantificationOptionsConfig] = QuantificationOptionsConfig(),
        initial_clustering_cfg: Optional[Any] = None,
        powder_meas_cfg: Optional[PowderMeasurementConfig] = PowderMeasurementConfig(),
        bulk_meas_cfg: Optional[BulkMeasurementConfig] = BulkMeasurementConfig(),
        exp_stds_cfg: Optional[ExpStandardsConfig] = ExpStandardsConfig(),
        plot_cfg: Optional[PlotConfig] = PlotConfig(),
        is_acquisition: bool = False,
        standards_dict: Optional[Union[EDSStandardsFile, Dict[str, Any]]] = None,
        development_mode: bool = False,
        output_filename_suffix: str = '',
        verbose: bool = True,
        results_dir: Optional[str] = None,
        sample_id: Optional[str] = None,
    ):
        """
        Initialize the EMXSp_Composition_Analyzer with all configuration objects.

        See class docstring for parameter documentation.
        """
        self.sample_id = self._resolve_sample_id(
            sample_id=sample_id,
            results_dir=results_dir,
            is_acquisition=is_acquisition,
        )

        # --- Record process time
        self.start_process_time = time.time()
        if verbose:
            print_double_separator()
            logger.info(f"▶️ Starting compositional analysis of sample {self.sample_id}")
            
            
        # --- Define use of class instance
        self.is_acquisition = is_acquisition
        self.development_mode = development_mode

        measurement_cfg = measurement_cfg.model_copy(
            update={
                "powder_meas_cfg": powder_meas_cfg if powder_meas_cfg is not None else measurement_cfg.powder_meas_cfg,
                "bulk_meas_cfg": bulk_meas_cfg if bulk_meas_cfg is not None else measurement_cfg.bulk_meas_cfg,
                "exp_stds_cfg": exp_stds_cfg if exp_stds_cfg is not None else measurement_cfg.exp_stds_cfg,
            }
        )
        is_XSp_measurement = measurement_cfg.type != measurement_cfg.PARTICLE_STATS_MEAS_TYPE_KEY
        
        
        # --- System characteristics
        self.microscope_cfg = microscope_cfg
        
        if is_XSp_measurement:
            # Load microscope calibrations for this instrument and mode
            calibs.load_microscope_calibrations(microscope_cfg.ID, measurement_cfg.mode, load_detector_channel_params=is_acquisition)
            if not measurement_cfg.emergence_angle:
                measurement_cfg.emergence_angle = calibs.emergence_angle # Fixed by instrument geometry
        
        
        # --- Measurement configurations
        self.measurement_cfg = measurement_cfg
        self.powder_meas_cfg = self.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig()
        self.bulk_meas_cfg = self.measurement_cfg.bulk_meas_cfg or BulkMeasurementConfig()
        self.exp_stds_cfg = self.measurement_cfg.exp_stds_cfg or ExpStandardsConfig()
        
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
            self._apply_geom_factors = True if sample_cfg.is_surface_rough else False
        
        
        # --- Fitting and Quantification
        self.quant_cfg = quant_cfg
        if standards_dict is None:
            self.standards = None
        elif isinstance(standards_dict, EDSStandardsFile):
            self.standards = standards_dict
        elif isinstance(standards_dict, dict):
            self.standards = EDSStandardsFile.from_payload(
                standards_dict,
                meas_type=self.measurement_cfg.type,
                beam_energy_keV=int(self.measurement_cfg.beam_energy_keV),
            )
        else:
            raise TypeError(
                "standards_dict must be None, EDSStandardsFile, or a dictionary payload"
            )

        if initial_clustering_cfg is None:
            self.clustering_cfg = LedgerClusteringConfig(features=cnst.AT_FR_CL_FEAT)
        elif hasattr(initial_clustering_cfg, "model_dump"):
            clustering_payload = initial_clustering_cfg.model_dump(mode="json")
            clustering_payload.pop("clustering_id", None)
            self.clustering_cfg = LedgerClusteringConfig.model_validate(clustering_payload)
        elif isinstance(initial_clustering_cfg, dict):
            clustering_payload = dict(initial_clustering_cfg)
            clustering_payload.pop("clustering_id", None)
            self.clustering_cfg = LedgerClusteringConfig.model_validate(clustering_payload)
        else:
            raise TypeError(
                "initial_clustering_cfg must be None, a mapping, or a model with model_dump()"
            )

        if is_XSp_measurement:
            # Set EDS detector channels to include in the quantification
            self.sp_start, self.sp_end = (
                int(round(quant_cfg.spectrum_lims[0])),
                int(round(quant_cfg.spectrum_lims[1])),
            )
            # Compute values of energies corresponding to detector channels
            if energy_zero is not None and bin_width is not None and float(bin_width) != 0.0:
                self.energy_vals = np.array([energy_zero + bin_width * i for i in range(self.sp_start, self.sp_end)])
            elif is_acquisition:
                raise ValueError("Missing detector calibration values.\n Please add detector calibration file at {calibs.calibration_files_dir}")
            # Set a threshold value below which counts are considered to be too low
            # Used to filter "bad" spectra out from clustering analysis. All spectra having less counts than this threshold are filtered out
            # Used also to avoid fitting spectra with excessive absorption, which inevitably lead to large quantification errors
            if self.clustering_cfg.min_bckgrnd_cnts is None:
                min_bckgrnd_cnts = measurement_cfg.target_acquisition_counts / (2 * 10**4)  # empirical value
                self.clustering_cfg.min_bckgrnd_cnts = int(round(min_bckgrnd_cnts))

        # --- Clustering
        if is_XSp_measurement:
            self.ref_formulae = list(dict.fromkeys(self.clustering_cfg.ref_formulae)) # Remove duplicates
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
                results_dir = make_unique_path(os.path.join(os.getcwd(), results_folder), self.sample_id)
                os.makedirs(results_dir)
            else:
                # When a project directory is provided, save under a per-sample folder.
                # If the provided path already points to a sample folder (resume), reuse it.
                normalized_results_dir = os.path.normpath(results_dir)
                is_sample_dir = os.path.basename(normalized_results_dir) == self.sample_id
                has_sample_artifacts = (
                    os.path.exists(os.path.join(results_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT))
                    or os.path.exists(os.path.join(results_dir, cnst.IMAGES_DIR))
                    or os.path.exists(os.path.join(results_dir, cnst.SPECTRA_DIR))
                )

                if is_sample_dir or has_sample_artifacts:
                    sample_result_dir = results_dir
                else:
                    sample_result_dir = os.path.join(results_dir, self.sample_id)

                os.makedirs(sample_result_dir, exist_ok=True)
                results_dir = sample_result_dir

        self.sample_result_dir = results_dir
        self.output_filename_suffix = output_filename_suffix
        self.verbose = verbose

        if is_XSp_measurement:
            # --- Variable initialization
            self.standards_dict = None
            self.particle_cntr = -1 # Counter to save particle number

            # List containing the quantification result records for each collected spectrum.
            # spectra_quant_records[i] is a QuantificationResult (or None when not yet quantified).
            # All composition/error/fit data needed for analysis is read directly from these records.
            self.spectra_quant_records = []
            self.current_quant_config: Optional[QuantificationConfig] = None
            self.current_quantification_id: Optional[int] = None
            
        
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
                self._initialise_acquisition_ledger()
        

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
    
        # Initialise instrument controller — configs are read from the ledger bundle
        self.EM_controller = EM_Controller(
            self._build_ledger_configs(),
            sample_id=self.sample_id,
            results_dir=EM_images_dir,
            verbose=self.verbose,
        )
        self.EM_controller.initialise_SEM()
        self.EM_controller.initialise_sample_navigator(exclude_sample_margin=True)
        
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

    def _initialise_acquisition_ledger(self) -> None:
        """Initialize acquisition ledger and reconcile any existing spectra pointers.

        This guarantees that if a ledger is missing but spectra files already exist
        under ``spectra/``, those spectra are ingested into the newly created ledger
        before acquisition restarts.
        """
        ledger = self._load_or_create_ledger()
        if (
            self.verbose
            and getattr(self, "_last_acq_ledger_created", False)
            and getattr(self, "_last_acq_ledger_ingested_spectra_count", 0) > 0
            and ledger is not None
        ):
            logger.info(
                "ℹ️ No ledger was found, but %d pre-existing spectrum file/s were detected in spectra/. "
                "These spectra were ingested into a newly created ledger, assuming they were acquired "
                "with the same configurations as the current run.",
                getattr(self, "_last_acq_ledger_ingested_spectra_count", 0),
            )

    #%% Other initializations
    # =============================================================================
    def _make_analysis_dir(self) -> None:
        """
        Create a deterministic directory for saving analysis results.

        The directory name is based on the active quantification and clustering config ids:
        ``analysis_Quant{quantification_id}_Clust{clustering_id}``.
        If the directory already exists (same active config pair), it is reused as-is.
        Existing files are preserved unless overwritten by newly generated outputs.
        A different active config pair results in a different directory name.
    
        The resulting directory path is stored in `self.analysis_dir`.
    
        Raises
        ------
        RuntimeError
            If active quantification/clustering ids cannot be resolved.
        OSError
            If directory creation fails.
        """
        quantification_id, clustering_id = self._resolve_active_analysis_config_ids()
        base_name = f"analysis_quant{quantification_id}_clust{clustering_id}"
        analysis_dir = os.path.join(self.sample_result_dir, base_name)

        try:
            os.makedirs(analysis_dir, exist_ok=True)
        except Exception as e:
            raise OSError(f"Could not create analysis directory '{analysis_dir}': {e}") from e

        self.analysis_dir = analysis_dir


    def _save_analysis_config_summary(self) -> None:
        """Write a human-readable plain-text summary of the quantification and
        clustering configurations used to produce this analysis folder.

        The file is written to ``self.analysis_dir/Analysis_config_summary.txt``.
        It is always overwritten so it reflects the configs that generated the
        current files in the folder.
        """
        lines: list[str] = []

        # Quantification section
        q  = self.quant_cfg             # QuantificationOptionsConfig
        qc = self.current_quant_config  # QuantificationConfig from ledger (may be None)
        if qc is None and self.sample_result_dir is not None:
            try:
                _ledger = self._load_or_create_ledger()
                qc = self._get_active_quantification_config(_ledger)
            except Exception:
                pass

        # Prefer options from the active QuantificationConfig stored in ledger,
        # because those are the exact scientific inputs used for this analysis.
        active_options = qc.options if (qc is not None and isinstance(qc.options, dict)) else {}
        method_used = str(active_options.get("method", q.method))
        fit_tolerance_used = float(active_options.get("fit_tolerance", q.fit_tolerance))
        use_instr_background_used = bool(
            active_options.get("use_instrument_background", q.use_instrument_background)
        )
        use_proj_std_dict_used = bool(
            active_options.get("use_project_specific_std_dict", q.use_project_specific_std_dict)
        )
        is_particle_used = bool(
            active_options.get("is_particle", getattr(self, "_apply_geom_factors", False))
        )
        beam_energy_used = float(
            active_options.get("beam_energy_keV", self.measurement_cfg.beam_energy_keV)
        )
        emergence_angle_used = active_options.get("emergence_angle", self.measurement_cfg.emergence_angle)
        emergence_angle_used = (
            float(emergence_angle_used)
            if emergence_angle_used is not None
            else None
        )
        det_ch_offset_used = active_options.get("det_ch_offset", self.det_ch_offset)
        det_ch_offset_used = (
            float(det_ch_offset_used)
            if det_ch_offset_used is not None
            else None
        )
        det_ch_width_used = active_options.get("det_ch_width", self.det_ch_width)
        det_ch_width_used = (
            float(det_ch_width_used)
            if det_ch_width_used is not None
            else None
        )
        spectrum_lims_used = active_options.get("spectrum_lims", q.spectrum_lims)
        if isinstance(spectrum_lims_used, (list, tuple)) and len(spectrum_lims_used) == 2:
            sp_low = int(float(spectrum_lims_used[0]))
            sp_high = int(float(spectrum_lims_used[1]))
        else:
            sp_low = int(float(q.spectrum_lims[0]))
            sp_high = int(float(q.spectrum_lims[1]))

        lines.append("=" * 60)
        lines.append("QUANTIFICATION CONFIGURATIONS")
        lines.append("=" * 60)
        if qc is not None:
            lines.append(f"  Quantification ID    : {qc.quantification_id}")
            if qc.label:
                lines.append(f"  Label                : {qc.label}")
            if qc.sample_elements:
                lines.append(f"  Sample elements      : {', '.join(qc.sample_elements)}")
            if qc.substrate_elements:
                lines.append(f"  Substrate elements   : {', '.join(qc.substrate_elements)}")
        lines.append(f"  Method               : {method_used}")
        lines.append(f"  Spectrum limits      : {sp_low} - {sp_high} channels")
        lines.append(f"  Fit tolerance        : {fit_tolerance_used:.2e}")
        lines.append(f"  Instrument background: {'yes' if use_instr_background_used else 'no'}")
        lines.append(f"  Project std dict     : {'yes' if use_proj_std_dict_used else 'no'}")
        lines.append(f"  Particle corrections : {'yes' if is_particle_used else 'no'}")
        lines.append(f"  Beam energy (keV)    : {beam_energy_used:.4f}")
        lines.append(
            f"  Emergence angle (deg): {emergence_angle_used:.4f}"
            if emergence_angle_used is not None
            else "  Emergence angle (deg): unavailable"
        )
        lines.append(
            f"  Detector offset      : {det_ch_offset_used:.6f}"
            if det_ch_offset_used is not None
            else "  Detector offset      : unavailable"
        )
        lines.append(
            f"  Detector width       : {det_ch_width_used:.6f}"
            if det_ch_width_used is not None
            else "  Detector width       : unavailable"
        )
        if qc is not None and qc.reference_lines_by_element:
            sample_elements = set(qc.sample_elements or [])
            reference_lines_by_element = [
                el_line
                for el, el_line in qc.reference_lines_by_element.items()
                if el in sample_elements
            ]
            ref_lines_str = ', '.join(reference_lines_by_element)
            if ref_lines_str:
                lines.append(f"  Reference lines      : {ref_lines_str}")

        # Clustering section
        cc = self.clustering_cfg

        lines.append("")
        lines.append("=" * 60)
        lines.append("CLUSTERING CONFIGURATIONS")
        lines.append("=" * 60)
        if hasattr(cc, "clustering_id"):
            lines.append(f"  Clustering ID        : {cc.clustering_id}")
        lines.append(f"  Method               : {cc.method}")
        lines.append(f"  Features             : {cc.features}")
        if cc.k_forced is not None:
            lines.append(f"  k (forced)           : {cc.k_forced}")
        else:
            lines.append(f"  k finding method     : {cc.k_finding_method}")
            lines.append(f"  Max k                : {cc.max_k}")
        if getattr(cc, "k_resolved", None) is not None:
            lines.append(f"  k (resolved)         : {cc.k_resolved}")
        lines.append(f"  Max analytical error : {cc.max_analytical_error_percent} w%")
        lines.append(
            f"  Min background counts: {cc.min_bckgrnd_cnts if cc.min_bckgrnd_cnts is not None else 'disabled'}"
        )
        lines.append(f"  Accepted quant flags : {cc.quant_flags_accepted}")
        lines.append(f"  Matrix decomposition : {'yes' if cc.do_matrix_decomposition else 'no'}")
        if cc.ref_formulae:
            lines.append("  Reference formulae   :")
            for formula in cc.ref_formulae:
                lines.append(f"    - {formula}")
        else:
            lines.append("  Reference formulae   : (none)")
        lines.append("")

        summary_path = os.path.join(
            self.analysis_dir,
            cnst.ANALYSIS_CONFIG_SUMMARY_FILENAME + ".txt",
        )
        try:
            with open(summary_path, "w", encoding="utf-8") as fh:
                fh.write("\n".join(lines))
        except OSError as e:
            logging.warning(f"Could not write analysis config summary: {e}")


    def _resolve_active_analysis_config_ids(self) -> Tuple[int, int]:
        """Resolve active quantification and clustering config ids for analysis folder naming."""
        quant_config = self.current_quant_config

        if quant_config is None:
            ledger = self._load_or_create_ledger()
            quant_config = self._get_active_quantification_config(ledger)

        if quant_config is None:
            raise RuntimeError("No active quantification config is available to name analysis directory")

        active_clustering_analysis = quant_config.get_active_clustering_analysis()
        if active_clustering_analysis is None:
            raise RuntimeError("No active clustering config is available to name analysis directory")

        return int(quant_config.quantification_id), int(active_clustering_analysis.config.clustering_id)


    @staticmethod
    def _clear_directory_contents(dir_path: str) -> None:
        """Remove all files/subdirectories from a directory while keeping the directory itself."""
        try:
            for entry in os.scandir(dir_path):
                if entry.is_dir(follow_symlinks=False):
                    shutil.rmtree(entry.path)
                else:
                    os.remove(entry.path)
        except Exception as e:
            raise OSError(f"Could not clear existing analysis directory '{dir_path}': {e}") from e
    
    
    def _initialise_std_dict(self) -> None:
        """
        Initialise the dictionary of X-ray standards for quantification.
    
        This method determines how the `standards_dict` attribute is initialised
        based on the sample configuration and measurement type:
    
        - If the measurement is of a known powder mixture, the standards dictionary
          is compiled from reference data using `_compile_standards_from_references()`.
    
        - Otherwise, the standards dictionary is expected to be loaded directly
          within the `XSp_Quantifier` and is set to `None` here.
    
        Returns
        -------
        None
            This method modifies the `self.standards_dict` attribute in place.
        """
        is_known_mixture = getattr(self.powder_meas_cfg, "is_known_powder_mixture_meas", False)
        
        if is_known_mixture:
            self.standards_dict = StandardsModule._compile_standards_from_references(self)
        elif self.quant_cfg.use_project_specific_std_dict:
            standards, _ = StandardsModule._load_standards(self)
            self.standards_dict = standards.to_standards_dict().get(self.measurement_cfg.mode, {})
        else:
            # Standards dictionary will be loaded directly within the `XSp_Quantifier`
            self.standards_dict = None


    def _calc_reference_phases_df(self) -> None:
        """
        Calculate the compositions of candidate phases and store them in a pd.DataFrame.
    
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
        
        valid_formulae = []
        valid_compositions = set()  # store normalized Composition keys, to check for duplicates

        for formula in self.ref_formulae:
            # Use pymatgen class Composition
            try:
                comp = Composition(formula)
            except Exception as e:
                warnings.warn(f"Invalid chemical formula '{formula}': {e}")
                continue
            
            # Normalize composition to a string key to check duplicates
            comp_key = comp.reduced_formula  # or str(comp) if you want exact
            if comp_key in valid_compositions:
                continue  # skip duplicate compositions
            
            valid_compositions.add(comp_key)
            valid_formulae.append(formula)
            
            # Get mass fractions as dictionary el: w_fr
            try:
                w_fr_dict = comp.as_weight_dict()
            except AttributeError:
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
                # Atomic fractions are normalised, so for the purpose of candidate phases we should calculate it normalising the
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
        
        # Copy all valid formulae back onto self.ref_formulae attribute
        self.ref_formulae = valid_formulae
        
        # Convert to pd.DataFrame and store it
        ref_phases_df = pd.DataFrame(ref_phases, columns=self.all_els_sample).fillna(0)
        self.ref_phases_df = ref_phases_df
    
        # Store values of reference weights used to calculate molar fractions from mixtures
        self.ref_weights_in_mixture = ref_weights_in_mixture
    
        # Calculate negative analytical error accepted to compensate for elements undetectable by EDS (H, He, Li, Be)
        self.undetectable_an_er = undetectable_an_err
        
        
    #%% Single spectrum operations
    # =============================================================================            
    def _acquire_spectrum(
        self,
        x: float,
        y: float,
        spectrum_id: str,
        msa_file_path: str,
    ) -> int:
        """
        Acquire an X-ray spectrum at the specified stage position and store the results.

        Parameters
        ----------
        x, y : float
            The x, y machine coordinates for the spectrum acquisition.
        spectrum_id : str
            Unique identifier for the spectrum.
        msa_file_path : str
            Path to save the acquired spectrum (.msa file).

        Returns
        -------
        total_counts : int
            Total counts in the acquired spectrum, derived from persisted data.

        Notes
        -----
        All spectrum data is saved to and loaded from .msa files. No in-memory spectrum caching is used.
        """
        # Acquire at the instrument and rely on persisted pointer files as source of truth.
        background_elements = None
        if self.quant_cfg.use_instrument_background:
            background_elements = list(getattr(self, "detectable_els_sample", []) or [])

        os.makedirs(os.path.dirname(msa_file_path), exist_ok=True)

        # Measure real acquisition time
        _acq_start = time.time()
        spectrum_data, background_data = self.EM_controller.acquire_XS_spot_spectrum(
            x, y,
            self.measurement_cfg.max_acquisition_time,
            self.measurement_cfg.target_acquisition_counts,
            elements=background_elements,
            msa_file_path=msa_file_path,
        )
        _acq_end = time.time()
        measured_real_time = _acq_end - _acq_start

        counts_arr = np.asarray(spectrum_data, dtype=float)

        pointer_path = Path(msa_file_path)
        if not pointer_path.exists():
            from autoemx.utils.legacy.spectrum_pointer_writer import write_spectrum_pointer_file
            write_spectrum_pointer_file(
                str(pointer_path),
                list(counts_arr),
                getattr(self, 'energy_vals', []),
                real_time=measured_real_time
            )
        loaded_counts = SampleLedger._load_counts_from_pointer_file(pointer_path)
        counts_arr = np.asarray(loaded_counts, dtype=float)

        # For user info, print the acquisition time
        if self.verbose:
            logger.info(f"✅ Acquisition took {measured_real_time:.2f} s")

        if self.quant_cfg.use_instrument_background and background_data is not None:
            self._write_manufacturer_background_vector(
                spectrum_id=spectrum_id,
                background_vals=list(map(float, background_data)),
            )
        elif self.quant_cfg.use_instrument_background and background_data is None:
            warnings.warn(
                "Instrument background retrieval failed during acquisition; "
                "falling back to automatic background subtraction for the "
                "remaining spectra in this run."
            )
            self.quant_cfg.use_instrument_background = False

        return int(np.round(np.sum(counts_arr)))
    
    
    def _fit_exp_std_spectrum(
        self,
        spectrum: Iterable,
        background: Optional[Iterable] = None,
        sp_collection_time: float = None,
        els_w_frs: Optional[Dict[str,float]] = None,
        sp_id: str = '',
        verbose: bool = True
    ) -> QuantificationResult:
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
        QuantificationResult
            Schema result for this spectrum including quant flag/comment diagnostics and
            fit metrics under ``fit_result``.
    
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
            logger.info('🔬 Fitting spectrum' + sp_id_str)
            start_quant_time = time.time()
                            
        # Check if spectrum is worth fitting
        is_sp_valid_for_fitting, quant_flag, comment = self._is_spectrum_valid_for_fitting(spectrum, background)
        if not is_sp_valid_for_fitting:
            quant_record = QuantificationResult(
                quantification_id=self.current_quantification_id,
                quant_flag=quant_flag,
                comment=comment,
                diagnostics=QuantificationDiagnostics(
                    converged=False,
                    interrupted=True,
                ),
            )
            return quant_record
        
        try:
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
                is_particle=self._apply_geom_factors,
                sp_collection_time=sp_collection_time,
                max_undetectable_w_fr=self.undetectable_an_er,
                fit_tol=self.quant_cfg.fit_tolerance,
                standards_dict=self.standards_dict,
                verbose=False,
                fitting_verbose=False
            )
            bad_quant_flag = quantifier.initialize_and_fit_spectrum(print_results=self.verbose)
            is_fit_valid = True
            min_bckgrnd_ref_lines = quantifier._get_min_bckgrnd_cnts_ref_quant_lines()
        except Exception as e:
            is_fit_valid = False
            logger.error(f"❌ {type(e).__name__}: {e}")
            traceback.print_exc()
            quant_flag, comment = self._check_fit_quant_validity(is_fit_valid, None, None, None)
            quant_record = QuantificationResult(
                quantification_id=self.current_quantification_id,
                quant_flag=quant_flag,
                comment=comment,
                diagnostics=QuantificationDiagnostics(
                    converged=False,
                    interrupted=True,
                ),
            )
            return quant_record
        
        _, are_all_ref_peaks_present = self._assemble_fit_info(quantifier)
        
        if are_all_ref_peaks_present:
            quant_flag, comment = self._check_fit_quant_validity(is_fit_valid, bad_quant_flag, quantifier, min_bckgrnd_ref_lines)
        else:
            comment = "Reference peak missing"
            quant_flag = 10
        
        if verbose:
            fit_time = time.time() - start_quant_time
            logger.info(f"✅ Fitting took {fit_time:.2f} s")

        quant_record = QuantificationResult(
            quantification_id=self.current_quantification_id,
            quant_flag=quant_flag,
            comment=comment,
            fit_result=quantifier.export_fit_result(),
            diagnostics=QuantificationDiagnostics(
                interrupted=False,
                min_background_ref_lines=min_bckgrnd_ref_lines,
            ),
        )

        return quant_record
    
    
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
                            logger.warning(f"⚠️ {el_line} reference peak missing.")
                    
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
        spectrum_index: Optional[int] = None,
        interrupt_fits_bad_spectra: bool = True,
        verbose: bool = True
    ) -> Tuple[Optional[Dict[str, Any]], QuantificationResult, Optional[float]]:
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
        Tuple[Optional[Dict[str, Any]], QuantificationResult, Optional[float]]
            Runtime quantifier payload (or None), persisted schema result, and elapsed quantification time.
    
        Notes
        -----
        - Filtering flags are appended through function _check_fit_quant_validity().
        """
        quantification_time = None
        if verbose:
            if sp_id != '':
                sp_id_str = " #" + sp_id
            else:
                sp_id_str = '...'
            print_single_separator()
            logger.info('🔬 Quantifying spectrum' + sp_id_str)
            start_quant_time = time.time()
                            
        is_sp_valid_for_fitting, quant_flag, comment = self._is_spectrum_valid_for_fitting(spectrum, background)
        if not is_sp_valid_for_fitting:
            quant_record = QuantificationResult(
                quantification_id=self.current_quantification_id,
                quant_flag=quant_flag,
                comment=comment,
                diagnostics=QuantificationDiagnostics(
                    converged=False,
                    interrupted=True,
                ),
            )
            if verbose:
                quantification_time = time.time() - start_quant_time
            return None, quant_record, quantification_time

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
            is_particle=self._apply_geom_factors,
            sp_collection_time=sp_collection_time,
            max_undetectable_w_fr=self.undetectable_an_er,
            fit_tol=self.quant_cfg.fit_tolerance,
            standards_dict=self.standards_dict,
            verbose=False,
            fitting_verbose=False
        )
        
        try:
            # Returns dictionary containing calculated composition in atomic fractions + analytical error
            quant_result, min_bckgrnd_ref_lines, bad_quant_flag = quantifier.quantify_spectrum(
                print_result=False,
                interrupt_fits_bad_spectra=interrupt_fits_bad_spectra
            )
            is_quant_fit_valid = True if quant_result is not None else False
        except Exception as e:
            logger.error(f"❌ {type(e).__name__}: {e}")
            is_quant_fit_valid = False
            quant_flag, comment = self._check_fit_quant_validity(
                is_quant_fit_valid,
                None,
                None,
                None,
            )
            quant_record = quantifier.export_quantification_result(
                quantification_id=self.current_quantification_id,
                quant_result=None,
                quant_flag=quant_flag,
                comment=comment,
            )
            if verbose:
                quantification_time = time.time() - start_quant_time
            return None, quant_record, quantification_time
        else:
            quant_flag, comment = self._check_fit_quant_validity(is_quant_fit_valid, bad_quant_flag, quantifier, min_bckgrnd_ref_lines)
            quant_record = quantifier.export_quantification_result(
                quantification_id=self.current_quantification_id,
                quant_result=quant_result,
                quant_flag=quant_flag,
                comment=comment,
            )

        if verbose:
            quantification_time = time.time() - start_quant_time

        return quant_result, quant_record, quantification_time


    def _get_ledger_path(self) -> str:
        """Return the ledger path for the current sample result directory."""
        return os.path.join(self.sample_result_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)


    def _resolve_or_create_spectrum_pointer(
        self,
        spectrum_id: str,
        spectrum_vals: List[float],
        *,
        live_time: Optional[float] = None,
        real_time: Optional[float] = None,
    ) -> str:
        """Return a relative spectrum pointer path and create the spectrum file when missing.

        This is used by ledger reconstruction/update code paths (including legacy
        Data.csv backfill). If no vendor template file is available in the sample
        directory, the spectrum is written with the minimal EMSA fallback format.
        """
        spectrum_relpath = self._build_spectrum_relpath(spectrum_id)
        spectrum_pointer_abs_path = os.path.join(self.sample_result_dir, spectrum_relpath)

        if os.path.exists(spectrum_pointer_abs_path):
            return spectrum_relpath

        if not hasattr(self, "_vendor_msa_template_lines"):
            setattr(
                self,
                "_vendor_msa_template_lines",
                load_vendor_msa_template_lines(self.sample_result_dir, cnst.MSA_SP_FILENAME),
            )

        write_spectrum_pointer_file(
            spectrum_pointer_abs_path,
            spectrum_vals,
            self.energy_vals,
            template_lines=getattr(self, "_vendor_msa_template_lines"),
            live_time=live_time,
            real_time=real_time,
        )
        return spectrum_relpath


    def _build_spectrum_relpath(self, spectrum_id: str) -> str:
        """Build the relative path for one raw spectrum pointer file."""
        filename = f"{cnst.SPECTRUM_FILENAME_PREFIX}{spectrum_id}{dflt.RAW_SPECTRUM_EXT}"
        return os.path.join(cnst.SPECTRA_DIR, filename)


    def _build_background_relpath(self, spectrum_id: str) -> str:
        """Build relative path for one manufacturer background vector file."""
        filename = (
            f"{cnst.SPECTRUM_FILENAME_PREFIX}{spectrum_id}"
            f"{cnst.SPECTRUM_MAN_BACKGROUND_SUFFIX}{cnst.VECTOR_FILEEXT}"
        )
        return os.path.join(cnst.SPECTRA_DIR, filename)


    def _write_manufacturer_background_vector(self, spectrum_id: str, background_vals: Optional[List[float]]) -> Optional[str]:
        """Persist manufacturer background counts as a companion vector file when enabled."""
        if background_vals is None or not self.quant_cfg.use_instrument_background:
            return None

        background_relpath = self._build_background_relpath(spectrum_id)
        background_abs_path = os.path.join(self.sample_result_dir, background_relpath)
        os.makedirs(os.path.dirname(background_abs_path), exist_ok=True)
        np.save(background_abs_path, np.asarray(list(map(float, background_vals)), dtype=float))
        return background_relpath


    @staticmethod
    def _load_background_vector_from_file(background_path: Path) -> Optional[np.ndarray]:
        """Load an instrument background vector from a sidecar file when available."""
        if not background_path.exists() or background_path.suffix.lower() != cnst.VECTOR_FILEEXT:
            return None
        try:
            background_vals = np.load(background_path, allow_pickle=False)
        except Exception:
            return None
        if background_vals is None:
            return None
        return np.asarray(background_vals, dtype=float)


    @staticmethod
    def _load_realtime_from_pointer_file(pointer_path: Path) -> Optional[float]:
        """Read REALTIME from an EMSA-like header when available."""
        if pointer_path.suffix.lower() not in {".msa", ".msg"}:
            return None

        try:
            with pointer_path.open("r", encoding="utf-8") as f:
                for raw_line in f:
                    line = raw_line.strip()
                    if not line.startswith("#") or ":" not in line:
                        continue
                    if line.upper().startswith("#SPECTRUM"):
                        break
                    key, value = line[1:].split(":", maxsplit=1)
                    key_norm = key.strip().replace("_", "").replace(" ", "").upper()
                    if key_norm == "REALTIME":
                        return float(value.strip())
        except Exception:
            return None

        return None


    def _load_existing_ledger(self) -> Optional[SampleLedger]:
        """Load an existing ledger if present and valid."""
        ledger_path = self._get_ledger_path()
        if not os.path.exists(ledger_path):
            return None
        try:
            return SampleLedger.from_json_file(ledger_path)
        except Exception:
            return None


    def _get_ingested_spectrum_ids(self) -> set[str]:
        """Return the spectrum ids already ingested for the current sample directory."""
        cached_sample_dir = getattr(self, "_ingested_spectrum_ids_sample_result_dir", None)
        ingested_spectrum_ids = getattr(self, "_ingested_spectrum_ids", None)

        if ingested_spectrum_ids is None or cached_sample_dir != self.sample_result_dir:
            ingested_spectrum_ids = set()
            ledger = self._load_existing_ledger()
            if ledger is not None:
                ingested_spectrum_ids.update(
                    str(entry.spectrum_id)
                    for entry in ledger.spectra
                    if entry.spectrum_id not in (None, "")
                )
            self._ingested_spectrum_ids = ingested_spectrum_ids
            self._ingested_spectrum_ids_sample_result_dir = self.sample_result_dir

        return ingested_spectrum_ids


    def _get_spectra_dir(self) -> str:
        """Return the absolute path to the spectra pointer directory."""
        return os.path.join(self.sample_result_dir, cnst.SPECTRA_DIR)


    def _list_pointer_files_in_spectra_dir(self) -> List[Path]:
        """List pointer files currently present in the spectra directory."""
        spectra_dir = Path(self._get_spectra_dir())
        if not spectra_dir.exists():
            return []

        allowed_ext = {".msa", ".msg", ".json"}
        ext_priority = {".msa": 0, ".msg": 1, ".json": 2}
        files_by_spectrum_id: Dict[str, Path] = {}
        for path in spectra_dir.iterdir():
            if not path.is_file() or path.suffix.lower() not in allowed_ext:
                continue
            stem = path.stem
            if not stem.startswith(cnst.SPECTRUM_FILENAME_PREFIX):
                continue
            if stem.endswith(cnst.SPECTRUM_MAN_BACKGROUND_SUFFIX):
                continue
            spectrum_id = stem[len(cnst.SPECTRUM_FILENAME_PREFIX):]
            existing = files_by_spectrum_id.get(spectrum_id)
            if existing is None:
                files_by_spectrum_id[spectrum_id] = path
                continue

            existing_priority = ext_priority.get(existing.suffix.lower(), 99)
            current_priority = ext_priority.get(path.suffix.lower(), 99)
            if current_priority < existing_priority:
                files_by_spectrum_id[spectrum_id] = path

        def sort_key(path: Path) -> Tuple[int, Union[int, str], str]:
            stem = path.stem
            spectrum_id = stem[len(cnst.SPECTRUM_FILENAME_PREFIX):] if stem.startswith(cnst.SPECTRUM_FILENAME_PREFIX) else stem
            if spectrum_id.isdigit():
                return (0, int(spectrum_id), path.name)
            return (1, spectrum_id.lower(), path.name)

        return sorted(files_by_spectrum_id.values(), key=sort_key)


    def _infer_max_particle_id_from_saved_images(self) -> Optional[int]:
        """Infer the maximum particle id from saved acquisition image filenames."""
        images_dir = Path(self.sample_result_dir, cnst.IMAGES_DIR)
        if not images_dir.exists():
            return None

        max_particle_id = None
        particle_pattern = re.compile(r"_par(\d+)_")
        for image_path in images_dir.iterdir():
            if not image_path.is_file():
                continue
            match = particle_pattern.search(image_path.name)
            if match is None:
                continue
            try:
                particle_id = int(match.group(1))
            except (TypeError, ValueError):
                continue
            if max_particle_id is None or particle_id > max_particle_id:
                max_particle_id = particle_id

        return max_particle_id


    @staticmethod
    def _parse_optional_int(value: Any) -> Optional[int]:
        if value is None or pd.isna(value) or value == "":
            return None
        try:
            return int(float(value))
        except Exception:
            return None


    @staticmethod
    def _parse_optional_float(value: Any) -> Optional[float]:
        if value is None or pd.isna(value) or value == "":
            return None
        try:
            return float(value)
        except Exception:
            return None


    def _build_spot_coordinates(
        self,
        machine_x: Any = None,
        machine_y: Any = None,
        pixel_x: Any = None,
        pixel_y: Any = None,
    ) -> Optional[SpotCoordinates]:
        machine_coords = None
        pixel_coords = None

        x_machine = self._parse_optional_float(machine_x)
        y_machine = self._parse_optional_float(machine_y)
        if x_machine is not None and y_machine is not None:
            machine_coords = Coordinate2D(x=x_machine, y=y_machine)

        x_pixel = self._parse_optional_float(pixel_x)
        y_pixel = self._parse_optional_float(pixel_y)
        if x_pixel is not None and y_pixel is not None:
            pixel_coords = (int(round(x_pixel)), int(round(y_pixel)))

        # If machine coordinates are missing, derive them from pixel coordinates
        # using the active frame context from EM_Controller.
        if machine_coords is None and pixel_coords is not None and hasattr(self, "EM_controller"):
            try:
                pos_abs_mm = self.EM_controller.convert_pixel_pos_to_mm(np.array(pixel_coords, dtype=float))
                machine_coords = Coordinate2D(
                    x=float(pos_abs_mm[0]),
                    y=float(pos_abs_mm[1]),
                )
            except Exception:
                machine_coords = None

        if machine_coords is None and pixel_coords is None:
            return None

        return SpotCoordinates(
            machine_coordinates=machine_coords,
            pixel_coordinates=pixel_coords,
        )


    def _build_spectrum_entry_from_pointer_file(
        self,
        pointer_file: Path,
        existing_results: Optional[List[QuantificationResult]] = None,
        acquisition_details_by_id: Optional[Dict[str, AcquisitionDetails]] = None,
        quantification_results_by_id: Optional[Dict[str, List[QuantificationResult]]] = None,
    ) -> SpectrumEntry:
        """Build a SpectrumEntry by inspecting one file under sample_path/spectra."""
        sample_root = Path(self.sample_result_dir).resolve()
        pointer_abs = pointer_file.resolve()
        pointer_rel = str(pointer_abs.relative_to(sample_root))
        stem = pointer_file.stem
        if stem.startswith(cnst.SPECTRUM_FILENAME_PREFIX):
            spectrum_id = stem[len(cnst.SPECTRUM_FILENAME_PREFIX):]
        else:
            spectrum_id = stem

        try:
            counts = SampleLedger._load_counts_from_pointer_file(pointer_abs)
            total_counts = int(round(float(np.sum(counts))))
        except Exception:
            total_counts = 0

        acquisition_details = AcquisitionDetails(frame_id=None, particle_id=None, spot_coordinates=None)
        if acquisition_details_by_id is not None:
            acquisition_details = acquisition_details_by_id.get(spectrum_id, acquisition_details)

        realtime_from_header = self._coerce_optional_finite_float(self._load_realtime_from_pointer_file(pointer_abs))
        background_relpath = None
        candidate_background = Path(
            self.sample_result_dir,
            self._build_background_relpath(spectrum_id),
        )
        if candidate_background.exists():
            background_relpath = str(candidate_background.resolve().relative_to(sample_root))

        entry_results = list(existing_results or [])
        if not entry_results and quantification_results_by_id is not None:
            entry_results = list(quantification_results_by_id.get(spectrum_id, []))

        return SpectrumEntry(
            live_acquisition_time=realtime_from_header if realtime_from_header is not None else 1.0,
            total_counts=total_counts,
            spectrum_id=spectrum_id,
            spectrum_relpath=pointer_rel,
            instrument_background_relpath=background_relpath,
            acquisition_details=acquisition_details,
            quantification_results=entry_results,
        )


    def _load_or_create_ledger(self) -> SampleLedger:
        """Load a ledger, using legacy Data.csv bootstrap only for non-acquisition workflows."""
        if self.is_acquisition:
            self._last_acq_ledger_created = False
            self._last_acq_ledger_ingested_spectra_count = 0

            spectra_dir = self._get_spectra_dir()
            os.makedirs(spectra_dir, exist_ok=True)

            ledger = self._load_existing_ledger()
            ledger_changed = False
            if ledger is None:
                ledger = SampleLedger(
                    sample_id=self.sample_id,
                    sample_path=os.path.abspath(self.sample_result_dir),
                    configs=self._build_ledger_configs(),
                    spectra=[],
                    quantifications=[],
                    active_quant=None,
                )
                self._last_acq_ledger_created = True
                ledger_changed = True

            n_ingested, n_removed = ingest_spectra(ledger)
            self._last_acq_ledger_ingested_spectra_count = n_ingested
            if n_ingested > 0 or n_removed > 0:
                ledger_changed = True

            if ledger.sample_path != os.path.abspath(self.sample_result_dir):
                ledger.sample_path = os.path.abspath(self.sample_result_dir)
                ledger_changed = True

            if ledger_changed:
                ledger.to_json_file(self._get_ledger_path())

            return ledger

        ledger = self._load_existing_ledger()
        if ledger is not None:
            return ledger

        if has_legacy_data_csv(self.sample_result_dir):
            default_ledger_configs = self._build_ledger_configs()
            default_ledger_configs, legacy_configs = resolve_legacy_bootstrap_configs(
                sample_result_dir=self.sample_result_dir,
                default_ledger_configs=default_ledger_configs,
            )

            # Keep runtime analyzer configs aligned with recovered legacy values.
            if legacy_configs is not None:
                self.microscope_cfg = default_ledger_configs.microscope_cfg
                self.sample_cfg = default_ledger_configs.sample_cfg
                self.measurement_cfg = default_ledger_configs.measurement_cfg
                self.sample_substrate_cfg = default_ledger_configs.sample_substrate_cfg
                self.powder_meas_cfg = self.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig()
                self.bulk_meas_cfg = self.measurement_cfg.bulk_meas_cfg or BulkMeasurementConfig()
                self.exp_stds_cfg = self.measurement_cfg.exp_stds_cfg or ExpStandardsConfig()
                self.plot_cfg = default_ledger_configs.plot_cfg
                if legacy_configs.get(cnst.QUANTIFICATION_CFG_KEY) is not None:
                    self.quant_cfg = legacy_configs[cnst.QUANTIFICATION_CFG_KEY]
                if legacy_configs.get(cnst.CLUSTERING_CFG_KEY) is not None:
                    self.clustering_cfg = legacy_configs[cnst.CLUSTERING_CFG_KEY]

                # Recompute energy axis using recovered legacy calibration values.
                legacy_energy_zero = self._coerce_optional_finite_float(getattr(self.microscope_cfg, "energy_zero", None))
                legacy_bin_width = self._coerce_optional_finite_float(getattr(self.microscope_cfg, "bin_width", None))
                self.sp_start, self.sp_end = (
                    int(round(self.quant_cfg.spectrum_lims[0])),
                    int(round(self.quant_cfg.spectrum_lims[1])),
                )
                if legacy_energy_zero is not None and legacy_bin_width is not None and float(legacy_bin_width) != 0.0:
                    self.energy_vals = np.array(
                        [legacy_energy_zero + legacy_bin_width * i for i in range(self.sp_start, self.sp_end)]
                    )
                elif not hasattr(self, "energy_vals"):
                    self.energy_vals = np.array([], dtype=float)

            return load_or_create_ledger_with_legacy_data_csv(
                sample_result_dir=self.sample_result_dir,
                sample_id=self.sample_id,
                microscope_id=self.microscope_cfg.ID,
                use_instrument_background=bool(self.quant_cfg.use_instrument_background),
                default_ledger_configs=default_ledger_configs,
                resolve_or_create_spectrum_pointer=self._resolve_or_create_spectrum_pointer,
                write_background_pointer=self._write_manufacturer_background_vector,
            )

        raise FileNotFoundError(
            f"No ledger.json or Data.csv found for sample '{self.sample_id}' "
            f"in '{self.sample_result_dir}'. Please ensure the sample directory contains "
            "either ledger.json or Data.csv."
        )


    def _extract_spectrum_info(self, spectrum: 'SpectrumEntry', index: int) -> dict:
        """Extract spectrum metadata from a ledger SpectrumEntry for use in CSV export."""
        par_id = frame_id = ""
        acq = spectrum.acquisition_details
        if acq is not None:
            par_id = str(acq.particle_id or "")
            frame_id = str(acq.frame_id or "")
        return {
            cnst.SP_ID_DF_KEY: str(spectrum.spectrum_id) if spectrum.spectrum_id is not None else str(index),
            cnst.PAR_ID_DF_KEY: par_id,
            cnst.FRAME_ID_DF_KEY: frame_id,
        }


    def _build_ledger_configs(self) -> LedgerConfigs:
        """Build inline ledger configs from current analyzer configuration objects."""
        self.measurement_cfg = self.measurement_cfg.model_copy(
            update={
                "powder_meas_cfg": self.powder_meas_cfg,
                "bulk_meas_cfg": self.bulk_meas_cfg,
                "exp_stds_cfg": self.exp_stds_cfg if self.exp_stds_cfg.is_exp_std_measurement else None,
            }
        )
        return LedgerConfigs(
            microscope_cfg=self.microscope_cfg,
            sample_cfg=self.sample_cfg,
            measurement_cfg=self.measurement_cfg,
            sample_substrate_cfg=self.sample_substrate_cfg,
            plot_cfg=self.plot_cfg,
        )


    def _ensure_quant_tracking_length(self, total_spectra: int) -> None:
        """Ensure in-memory quantification tracking lists are indexable for all spectra."""
        if len(self.spectra_quant_records) < total_spectra:
            self.spectra_quant_records.extend([None] * (total_spectra - len(self.spectra_quant_records)))


    def _ensure_current_quantification_run(self, force_new: bool = False) -> None:
        """Resolve the active quantification run for this launch.

        Parameters
        ----------
        force_new : bool
            If True, always create a new quantification config id.
        """
        existing_ledger = self._load_or_create_ledger()
        active_quant_config = self._get_active_quantification_config(existing_ledger)
        active_quant_has_results = (
            active_quant_config is not None
            and self._ledger_has_quantification_results(
                existing_ledger,
                int(active_quant_config.quantification_id),
            )
        )

        current_sample_elements = list(self.all_els_sample)
        current_substrate_elements = list(self.all_els_substrate)
        current_options = self._build_quantification_options()
        current_reference_values = self._get_reference_values_by_el_line(active_quant_config=active_quant_config)

        candidate_id = (
            active_quant_config.quantification_id
            if active_quant_config is not None
            else self._next_quantification_id(existing_ledger)
        )
        candidate_config = self._build_quantification_config(
            quantification_id=candidate_id,
            sample_elements=current_sample_elements,
            substrate_elements=current_substrate_elements,
            options=current_options,
            reference_values_by_el_line=current_reference_values,
        )

        if (
            not force_new
            and active_quant_config is not None
            and self._quantification_configs_match(active_quant_config, candidate_config)
        ):
            self.current_quant_config = active_quant_config
            self.current_quantification_id = active_quant_config.quantification_id
            # Even when reusing quantification settings, clustering inputs (e.g., ref_formulae
            # passed via runners as sample['cnd']) may have changed for this run.
            self._ensure_current_clustering_run(
                self.current_quant_config,
            )
            self._apply_active_clustering_config(self.current_quant_config)
            return

        if not force_new and active_quant_config is not None:
            changes = active_quant_config.fingerprint_differences(candidate_config)
            if changes:
                changed_summary = self._format_quantification_config_changes(changes)
                if not active_quant_has_results:
                    warnings.warn(
                        "Quantification scientific inputs changed; replacing unused active quantification config. "
                        f"Changed fields: {changed_summary}",
                        UserWarning,
                    )
                    self.current_quant_config = candidate_config
                    self.current_quantification_id = candidate_config.quantification_id
                    self._ensure_current_clustering_run(
                        self.current_quant_config,
                    )
                    self._apply_active_clustering_config(self.current_quant_config)
                    return
                warnings.warn(
                    "Quantification scientific inputs changed; creating a new quantification config. "
                    f"Changed fields: {changed_summary}",
                    UserWarning,
                )

        quantification_id = self._next_quantification_id(existing_ledger)
        self.current_quant_config = self._build_quantification_config(
            quantification_id=quantification_id,
            sample_elements=current_sample_elements,
            substrate_elements=current_substrate_elements,
            options=current_options,
            reference_values_by_el_line=current_reference_values,
        )
        self.current_quantification_id = quantification_id
        self._ensure_current_clustering_run(self.current_quant_config)
        self._apply_active_clustering_config(self.current_quant_config)


    def _build_quantification_options(self) -> Dict[str, Any]:
        """Build the subset of quantification options that defines result reuse."""
        return {
            "method": self.quant_cfg.method,
            "spectrum_lims": [
                float(self.quant_cfg.spectrum_lims[0]),
                float(self.quant_cfg.spectrum_lims[1]),
            ],
            "fit_tolerance": float(self.quant_cfg.fit_tolerance),
            "use_instrument_background": bool(self.quant_cfg.use_instrument_background),
            # Effective quantifier state values that materially affect quantification output.
            "is_particle": bool(self._apply_geom_factors),
            "beam_energy_keV": float(self.measurement_cfg.beam_energy_keV),
            "emergence_angle": float(self.measurement_cfg.emergence_angle),
            "det_ch_offset": float(self.det_ch_offset),
            "det_ch_width": float(self.det_ch_width),
        }


    def _get_reference_values_by_el_line(
        self,
        active_quant_config: Optional[QuantificationConfig] = None,
    ) -> Dict[str, Any]:
        """Return method-dependent reference values, reusing persisted config values when possible."""
        if active_quant_config is not None and active_quant_config.reference_values_by_el_line:
            cached_reference_values = dict(sorted(active_quant_config.reference_values_by_el_line.items()))
            current_reference_values = self._extract_reference_values_from_standards(load_if_missing=False)
            if current_reference_values and self._reference_values_changed(
                current_reference_values,
                cached_reference_values,
                decimals=2,
            ):
                warnings.warn(
                    "Reference values in standards differ from active quantification config "
                    "(beyond 2-decimal tolerance); "
                    "a new quantification config will be created.",
                    UserWarning,
                )
                return current_reference_values
            return cached_reference_values

        return self._extract_reference_values_from_standards(load_if_missing=True)


    @staticmethod
    def _rounded_reference_values(
        reference_values_by_el_line: Dict[str, Any],
        decimals: int,
    ) -> Dict[str, float]:
        """Return deterministic rounded reference values keyed by element-line."""
        rounded: Dict[str, float] = {}
        for el_line, value in reference_values_by_el_line.items():
            rounded[str(el_line)] = round(float(value), decimals)
        return dict(sorted(rounded.items()))


    @classmethod
    def _reference_values_changed(
        cls,
        current_reference_values: Dict[str, Any],
        cached_reference_values: Dict[str, Any],
        decimals: int,
    ) -> bool:
        """Return True when reference values differ beyond configured decimal tolerance."""
        rounded_current = cls._rounded_reference_values(current_reference_values, decimals)
        rounded_cached = cls._rounded_reference_values(cached_reference_values, decimals)
        return rounded_current != rounded_cached


    def _extract_reference_values_from_standards(self, load_if_missing: bool) -> Dict[str, Any]:
        """Extract per-reference-line values from standards for the configured quantification method."""
        standards_by_line = None

        if self.quant_cfg.use_project_specific_std_dict:
            if self.standards_dict is not None and not load_if_missing:
                standards_by_line = self.standards_dict
            else:
                standards, _ = StandardsModule._load_standards(self)
                standards_by_line = standards.standards_by_mode.get(self.measurement_cfg.mode, {})
        else:
            standards_by_line = self.standards_dict
            if standards_by_line is None and self.standards is not None:
                if isinstance(self.standards, EDSStandardsFile):
                    standards_by_line = self.standards.standards_by_mode.get(
                        self.measurement_cfg.mode,
                        {},
                    )
                elif self.measurement_cfg.mode in self.standards:
                    standards_by_line = self.standards[self.measurement_cfg.mode]
                else:
                    standards_by_line = self.standards
            if standards_by_line is None and load_if_missing:
                standards, _ = StandardsModule._load_standards(self)
                standards_by_line = standards.standards_by_mode.get(self.measurement_cfg.mode, {})

        if standards_by_line is None:
            return {}

        relevant_elements = set(self.detectable_els_sample) | set(self.detectable_els_substrate)
        reference_values_by_el_line: Dict[str, Any] = {}
        if self.quant_cfg.method == "PB":
            for el_line, std_values in standards_by_line.items():
                element = el_line.split("_", maxsplit=1)[0]
                if element not in relevant_elements:
                    continue

                if hasattr(std_values, "reference_mean"):
                    reference_mean = getattr(std_values, "reference_mean", None)
                    if reference_mean is None:
                        continue
                    reference_values_by_el_line[el_line] = float(reference_mean.corrected_pb)
                    continue

                if isinstance(std_values, dict):
                    reference_mean = std_values.get("reference_mean")
                    if isinstance(reference_mean, dict) and cnst.COR_PB_DF_KEY in reference_mean:
                        reference_values_by_el_line[el_line] = float(reference_mean[cnst.COR_PB_DF_KEY])
                        continue
                    std_values = std_values.get("entries", [])

                mean_std = next(
                    (std for std in std_values if std.get(cnst.STD_ID_KEY) == cnst.STD_MEAN_ID_KEY),
                    None,
                )
                if mean_std is None or cnst.COR_PB_DF_KEY not in mean_std:
                    continue
                reference_values_by_el_line[el_line] = float(mean_std[cnst.COR_PB_DF_KEY])

        return dict(sorted(reference_values_by_el_line.items()))


    @staticmethod
    def _quantification_configs_match(left: QuantificationConfig, right: QuantificationConfig) -> bool:
        """Compare quantification configs by deterministic scientific-input fingerprint."""
        return left.fingerprint() == right.fingerprint()

    @staticmethod
    def _clustering_configs_match(left: LedgerClusteringConfig, right: LedgerClusteringConfig) -> bool:
        """Compare clustering configs by deterministic scientific-input fingerprint."""
        return left.fingerprint() == right.fingerprint()

    @staticmethod
    def _format_quantification_config_changes(
        changes: Dict[str, Dict[str, Any]],
        max_items: int = 8,
        max_value_length: int = 80,
    ) -> str:
        """Return a concise deterministic summary of config changes with old->new values."""
        if not changes:
            return "none"

        sorted_items = sorted(changes.items(), key=lambda item: item[0])
        formatted_entries: List[str] = []
        for path, change in sorted_items[:max_items]:
            old_val = EMXSp_Composition_Analyzer._short_repr(change.get("old"), max_value_length)
            new_val = EMXSp_Composition_Analyzer._short_repr(change.get("new"), max_value_length)
            formatted_entries.append(f"{path}: {old_val} -> {new_val}")

        remaining = len(sorted_items) - len(formatted_entries)
        if remaining > 0:
            formatted_entries.append(f"... and {remaining} more")

        return "; ".join(formatted_entries)

    @staticmethod
    def _short_repr(value: Any, max_length: int) -> str:
        """Return a stable shortened representation suitable for warning messages."""
        text = repr(value)
        if len(text) <= max_length:
            return text
        return text[: max_length - 3] + "..."


    @staticmethod
    def _runtime_clustering_cfg_payload(clustering_cfg: LedgerClusteringConfig) -> Dict[str, Any]:
        """Convert a ledger clustering config to the runtime config payload shape."""
        payload = clustering_cfg.model_dump(mode="json")
        payload.pop("clustering_id", None)
        return payload


    def _apply_active_clustering_config(self, quant_config: QuantificationConfig) -> None:
        """Load runtime clustering options from the active clustering config of a quantification config."""
        active_clustering_analysis = quant_config.get_active_clustering_analysis()
        if active_clustering_analysis is None:
            raise ValueError(
                "Active quantification config does not contain an active clustering config"
            )

        self.clustering_cfg = active_clustering_analysis.config.model_copy(deep=True)
        self.ref_formulae = list(dict.fromkeys(self.clustering_cfg.ref_formulae))
        self._calc_reference_phases_df()


    def _build_clustering_config_descriptor(self, clustering_id: int) -> LedgerClusteringConfig:
        """Build a persisted clustering config descriptor for the current analysis settings."""
        return LedgerClusteringConfig(
            clustering_id=clustering_id,
            method=self.clustering_cfg.method,
            features=self.clustering_cfg.features,
            k_forced=self.clustering_cfg.k_forced,
            k_resolved=self.clustering_cfg.k_resolved,
            k_finding_method=self.clustering_cfg.k_finding_method,
            max_k=self.clustering_cfg.max_k,
            ref_formulae=list(self.clustering_cfg.ref_formulae),
            do_matrix_decomposition=self.clustering_cfg.do_matrix_decomposition,
            max_analytical_error_percent=self.clustering_cfg.max_analytical_error_percent,
            min_bckgrnd_cnts=self.clustering_cfg.min_bckgrnd_cnts,
            quant_flags_accepted=list(self.clustering_cfg.quant_flags_accepted),
        )


    def _ensure_current_clustering_run(
        self,
        quant_config: QuantificationConfig,
    ) -> None:
        """Ensure active quantification config tracks clustering run history and active config."""
        candidate_clustering_config = self._build_clustering_config_descriptor(
            clustering_id=self._next_clustering_config_id(quant_config)
        )
        matching_index: Optional[int] = None
        runtime_clustering_id = getattr(self.clustering_cfg, "clustering_id", None)
        for idx, analysis in enumerate(quant_config.clustering_analyses):
            if not self._clustering_configs_match(analysis.config, candidate_clustering_config):
                continue
            if runtime_clustering_id is not None and analysis.config.clustering_id == runtime_clustering_id:
                matching_index = idx
                break
            if matching_index is None:
                matching_index = idx

        if matching_index is not None:
            quant_config.active_clustering_analysis_index = matching_index
            return

        active_clustering_analysis = quant_config.get_active_clustering_analysis()
        active_clustering_config = (
            active_clustering_analysis.config if active_clustering_analysis is not None else None
        )
        if active_clustering_config is not None:
            changes = active_clustering_config.fingerprint_differences(candidate_clustering_config)
            if changes:
                changed_summary = self._format_quantification_config_changes(changes)
                warnings.warn(
                    "Clustering scientific inputs changed; appending a new clustering config to the active "
                    "quantification config. "
                    f"Changed fields: {changed_summary}",
                    UserWarning,
                )

        quant_config.clustering_analyses.append(
            ClusteringAnalysis(config=candidate_clustering_config, result=None)
        )
        quant_config.active_clustering_analysis_index = len(quant_config.clustering_analyses) - 1


    @staticmethod
    def _next_clustering_config_id(quant_config: QuantificationConfig) -> int:
        """Return the next clustering config id local to the provided quantification config."""
        if not quant_config.clustering_analyses:
            return 0
        return max(analysis.config.clustering_id for analysis in quant_config.clustering_analyses) + 1


    @staticmethod
    def _get_active_quantification_config(existing_ledger: Optional[SampleLedger]) -> Optional[QuantificationConfig]:
        """Return the active quantification config from the ledger when available."""
        if existing_ledger is None or not existing_ledger.quantifications:
            return None

        if existing_ledger.active_quant is not None:
            for quant_config in existing_ledger.quantifications:
                if quant_config.quantification_id == existing_ledger.active_quant:
                    return quant_config

        return existing_ledger.quantifications[-1]


    @staticmethod
    def _next_quantification_id(existing_ledger: Optional[SampleLedger]) -> int:
        """Return the next ledger-local integer quantification identifier."""
        if existing_ledger is None or not existing_ledger.quantifications:
            return 0
        return max(config.quantification_id for config in existing_ledger.quantifications) + 1


    @staticmethod
    def _ledger_has_quantification_results(
        existing_ledger: Optional[SampleLedger],
        quantification_id: int,
    ) -> bool:
        """Return True when at least one spectrum contains results for the given quantification id."""
        if existing_ledger is None:
            return False
        for spectrum in existing_ledger.spectra:
            for result in spectrum.quantification_results:
                if result.quantification_id == quantification_id:
                    return True
        return False


    def _build_quantification_config(
        self,
        quantification_id: int,
        sample_elements: List[str],
        substrate_elements: List[str],
        options: Dict[str, Any],
        reference_values_by_el_line: Dict[str, Any],
    ) -> QuantificationConfig:
        """Build the persisted descriptor for the current quantification run."""
        label = f"Quant {quantification_id} {self.output_filename_suffix}"
        if isinstance(self.output_filename_suffix, str) and self.output_filename_suffix:
            label += self.output_filename_suffix
        reference_lines_by_element = QuantificationConfig.derive_reference_lines_by_element(
            reference_values_by_el_line=reference_values_by_el_line,
            preferred_lines=XSp_Quantifier.xray_quant_ref_lines,
        )
        return QuantificationConfig(
            quantification_id=quantification_id,
            label=label,
            sample_elements=sample_elements,
            substrate_elements=substrate_elements,
            els_w_fr=self._get_forced_mass_fractions(),
            options=options,
            reference_values_by_el_line=reference_values_by_el_line,
            reference_lines_by_element=reference_lines_by_element,
            clustering_analyses=[],
            active_clustering_analysis_index=None,
        )


    def _get_forced_mass_fractions(self) -> Optional[Dict[str, float]]:
        """Return forced elemental mass fractions for the active run, if any."""
        if getattr(self.exp_stds_cfg, "is_exp_std_measurement", False):
            forced = getattr(self.exp_stds_cfg, "w_frs", None)
            if isinstance(forced, dict) and forced:
                return {str(el): float(fr) for el, fr in forced.items()}

        forced = getattr(self.sample_cfg, "w_frs", None)
        if isinstance(forced, dict) and forced:
            return {str(el): float(fr) for el, fr in forced.items()}

        return None


    def _upsert_current_quantification_config_on_ledger(self, ledger: SampleLedger) -> None:
        """Insert or replace the current quantification config inside the provided ledger, checking for true scientific equivalence using fingerprint."""
        if self.current_quant_config is None:
            return

        # Check if an equivalent quantification config already exists in the ledger
        for idx, config in enumerate(ledger.quantifications):
            if self._quantification_configs_match(config, self.current_quant_config):
                # Equivalent quantification inputs found: preserve the existing quantification_id
                # while updating clustering history and active clustering index from current runtime state.
                merged_config = self.current_quant_config.model_copy(
                    update={"quantification_id": config.quantification_id}
                )
                ledger.quantifications[idx] = merged_config
                self.current_quant_config = merged_config
                self.current_quantification_id = merged_config.quantification_id
                ledger.active_quant = merged_config.quantification_id
                return

        # Otherwise, upsert by quantification_id as before
        existing_index = next(
            (
                i
                for i, config in enumerate(ledger.quantifications)
                if config.quantification_id == self.current_quant_config.quantification_id
            ),
            None,
        )

        if existing_index is None:
            ledger.append_quantification_config(self.current_quant_config)
        else:
            ledger.quantifications[existing_index] = self.current_quant_config
        ledger.active_quant = self.current_quant_config.quantification_id


    def _persist_current_quantification_config(self) -> None:
        """Persist the currently selected quantification config as the ledger active run."""
        if self.current_quant_config is None or self.current_quantification_id is None:
            return

        ledger = self._load_or_create_ledger()
        self._upsert_current_quantification_config_on_ledger(ledger)

        # Ensure the full current config payload (including clustering history/index)
        # is written for the active quantification id.
        for idx, config in enumerate(ledger.quantifications):
            if config.quantification_id == self.current_quantification_id:
                ledger.quantifications[idx] = self.current_quant_config.model_copy(deep=True)
                break

        ledger.active_quant = self.current_quantification_id
        ledger.to_json_file(self._get_ledger_path())


    def _sync_existing_quantification_from_ledger(self) -> None:
        """Populate in-memory quantification slots from ledger results with the current id."""
        if self.current_quantification_id is None:
            return
        existing_ledger = self._load_or_create_ledger()
        if existing_ledger is None:
            return

        total_spectra = len(existing_ledger.spectra)
        self._ensure_quant_tracking_length(total_spectra)

        for index, spectrum in enumerate(existing_ledger.spectra):
            matching = next(
                (
                    result for result in spectrum.quantification_results
                    if result.quantification_id == self.current_quantification_id
                ),
                None,
            )
            if matching is None:
                continue

            self.spectra_quant_records[index] = matching


    def _persist_quantification_record(self, spectrum_index: int, quant_record: QuantificationResult, overwrite: bool = False) -> None:
        """Write one quantification result to the ledger immediately for interruption-safe progress.

        Parameters
        ----------
        overwrite : bool
            If True and a result with the same quantification_id already exists for this spectrum,
            replace it (used when requantify_only_unquantified_spectra=True).
        """
        existing_ledger = self._load_or_create_ledger()
        ledger = existing_ledger

        if self.current_quant_config is None:
            raise ValueError("Current quantification config is not initialized")

        self._upsert_current_quantification_config_on_ledger(ledger)

        ledger.active_quant = self.current_quantification_id

        existing_ids = {
            existing.quantification_id
            for existing in ledger.spectra[spectrum_index].quantification_results
        }
        if quant_record.quantification_id not in existing_ids:
            ledger.append_quantification_result(spectrum_index, quant_record)
        elif overwrite:
            ledger.spectra[spectrum_index].quantification_results = [
                quant_record if r.quantification_id == quant_record.quantification_id else r
                for r in ledger.spectra[spectrum_index].quantification_results
            ]

        ledger.to_json_file(self._get_ledger_path())


    def _check_fit_quant_validity(
        self,
        is_quant_fit_valid: bool,
        bad_quant_flag: int,
        quantifier: Any,
        min_bckgrnd_ref_lines: Any
    ) -> tuple[int, str]:
        """
        Determine the quantification flag and comment for a spectrum based on fit outcomes.
    
        Parameters
        ----------
        is_quant_fit_valid : bool
            Whether the spectrum fit and quantification succeeded without errors.
        bad_quant_flag : int
            Indicator of the type of issue detected during fitting:
            - 1: poor fit
            - 2: excessively high analytical error
            - 3: excessive absorption
            - -1: non-converged fit
        quantifier : object
            The quantifier instance used for this spectrum; may be used for additional checks.
        min_bckgrnd_ref_lines : Any
            Reference value for background lines, used for further spectrum checks.
    
        Returns
        -------
        quant_flag : int
            Numerical flag representing the spectrum quality after fit/quantification.
        comment : str
            Human-readable comment describing the outcome or issue detected.
        """
        # Prefix for comments if fit was interrupted
        start_str_comments = 'Fit interrupted due to ' if not is_quant_fit_valid else ''
    
        if bad_quant_flag == 1:
            if self.verbose and is_quant_fit_valid:
                logger.warning("⚠️ Flagged for poor fit")
            comment = start_str_comments + "poor fit"
            quant_flag = 4
        elif bad_quant_flag == 2:
            if self.verbose and is_quant_fit_valid:
                logger.warning("⚠️ Flagged for excessively high analytical error")
            comment = start_str_comments + "excessively high analytical error"
            quant_flag = 5
        elif bad_quant_flag == 3:
            if self.verbose and is_quant_fit_valid:
                logger.warning("⚠️ Flagged for excessive X-ray absorption")
            comment = start_str_comments + "excessive X-ray absorption"
            quant_flag = 6
        elif not is_quant_fit_valid:
            comment = "Fit interrupted for unknown reasons"
            quant_flag = 9
        else:
            # Fit completed with no apparent issue; check for low background counts, etc.
            _, quant_flag, comment = self._flag_spectrum_for_clustering(min_bckgrnd_ref_lines, quantifier)
    
        # If fit was good but did not converge, annotate comment safely and set flag
        if bad_quant_flag == -1 and quant_flag == 0:
            if comment:
                comment += " - Quantification did not converge."
            else:
                comment = "Quantification did not converge."
            quant_flag = -1  # Signal non-convergence
        
        return quant_flag, comment
        
    def _is_spectrum_valid_for_fitting(
        self, 
        spectrum: np.ndarray, 
        background: Optional[np.ndarray] = None,
    ) -> tuple[bool, int, str]:
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
        background : Optional[np.ndarray], optional
            Reserved for interface compatibility with call sites and potential
            future checks. Currently not used.
    
        Returns
        -------
        is_spectrum_valid : bool
            True if the spectrum is valid for fitting, False otherwise.
        quant_flag : int
            Numerical flag representing the spectrum quality after fit/quantification.
        comment : str
            Human-readable comment describing the outcome or issue detected.
    
        Notes
        -----
        - Assumes all class attributes and keys are correctly initialized.
        - Uses constants from `cnst` for comment and flag keys.
        """
        is_spectrum_valid = True
        quant_flag = None
        comment = None
        
        if spectrum is None:
            # Check if spectrum data is present
            is_spectrum_valid = False
            comment = "No spectral data present"
            quant_flag = 1
            if self.verbose:
                logger.error("❌ Error during spectrum collection. No quantification was done.")
        elif np.sum(spectrum) < 0.9 * self.measurement_cfg.target_acquisition_counts:
            # Skip quantification of spectrum when counts are too low
            is_spectrum_valid = False
            comment = "Total counts too low"
            quant_flag = 2
            if self.verbose:
                logger.info(f"⏭️ Quantification skipped due to spectrum counts lower than 90% of the target counts of {self.measurement_cfg.target_acquisition_counts}")
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
            min_background_threshold = self.clustering_cfg.min_bckgrnd_cnts
            if min_background_threshold is not None and all(min_vals < min_background_threshold):
                is_spectrum_valid = False
                comment = "Background counts too low"
                quant_flag = 3
                if self.verbose:
                    logger.info(f"⏭️ Quantification skipped due to at least {n_vals_considered} spectrum points with E < {en_threshold} keV having a count lower than {min_background_threshold}")
                    logger.warning("⚠️ This generally indicates an excessive absorption of X-rays before they reach the detector, which compromises accurate measurements of PB ratios.")
    
        return is_spectrum_valid, quant_flag, comment
    
    
    def _flag_spectrum_for_clustering(
        self,
        min_bckgrnd_ref_lines: float,
        quantifier: Any,
    ) -> tuple[bool, int, str]:
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
        quant_flag : int
            Numerical flag representing the spectrum quality after fit/quantification.
        comment : str
            Human-readable comment describing the outcome or issue detected.
    
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
                comment = f"{el} {peak_int:.0f} counts > {sub_peak_int_threshold} % of total counts"
                quant_flag = 7
                if self.verbose:
                    logger.warning(f"⚠️ Intensity of substrate element {el} is {peak_int:.0f} cnts, larger than {sub_peak_int_threshold}% of total counts")
                    logger.warning("⚠️ This is likely to lead to large quantification errors.")
                break  # Stop if one element has too high intensity
    
        # Check that background intensity is high enough
        if is_spectrum_valid:
            comment = ""
            min_background_threshold = self.clustering_cfg.min_bckgrnd_cnts
            # Spectrum is not valid if any reference peak has average counts lower than the configured threshold.
            if min_background_threshold is not None and min_bckgrnd_ref_lines < min_background_threshold:
                is_spectrum_valid = False
                comment = (
                    f"Reference background counts too low "
                    f"({min_bckgrnd_ref_lines:.1f} < {min_background_threshold})"
                )
                quant_flag = 8
                if self.verbose:
                    logger.warning(f"⚠️ Counts below a reference peak are on average < {min_background_threshold}")
                    logger.warning("⚠️ This is likely to lead to large quantification errors.")
            else:
                quant_flag = 0  # Quantification is ok
    
        return is_spectrum_valid, quant_flag, comment  # Not used, but returned for completeness
    
    
    #%% Spectra acquisition and quantification routines
    # ============================================================================= 
    def _collect_spectra(
        self,
        n_spectra_to_collect: int,
        n_tot_sp_collected: Optional[int] = None,
        quantify: bool = True,
        interrupt_fits_bad_spectra: bool = True,
        particle_id_offset: int = 0,
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
            Next spectrum ID to use. If None, it is inferred as max(existing_id) + 1.
        quantify : bool, optional
            If True, perform spectra quantification (default: True).
    
        Returns
        -------
        n_tot_sp_collected : int
            Updated next spectrum ID to use.
        success : bool
            False if collection was interrupted by user, or if no more particles could be found. True otherwise.
    
        Notes
        -----
        - If `quantify` is True, quantification occurs immediately after each collection.
        """
        success = False
        ingested_spectrum_ids = self._get_ingested_spectrum_ids()

        # Auto-detect next available spectrum ID if not provided
        if n_tot_sp_collected is None:
            pointer_files = self._list_pointer_files_in_spectra_dir()
            max_id = -1
            for pf in pointer_files:
                stem = pf.stem
                if stem.startswith(cnst.SPECTRUM_FILENAME_PREFIX):
                    spectrum_id = stem[len(cnst.SPECTRUM_FILENAME_PREFIX):]
                    if spectrum_id.isdigit():
                        max_id = max(max_id, int(spectrum_id))
            n_tot_sp_collected = max_id + 1

        next_spectrum_id = n_tot_sp_collected
        n_spectra_collected = 0

        while n_spectra_collected < n_spectra_to_collect:
            success, spots_xy_list, particle_cntr = self.EM_controller.get_XSp_coords(next_spectrum_id)
            
            if not success:
                break
            
            if particle_cntr:
                self.particle_cntr = particle_cntr + particle_id_offset
            else:
                self.particle_cntr = None
            frame_ID = self.EM_controller.current_frame_label
            
            latest_spot_id = None # For image annotations
            for i, (x, y) in enumerate(spots_xy_list):
                latest_spot_id = i
                xy_center = (int(x), int(y))
                current_spectrum_id = str(next_spectrum_id)

                if current_spectrum_id in ingested_spectrum_ids:
                    next_spectrum_id += 1
                    continue

                if self.verbose:
                    print_single_separator()
                    logger.info(f'🔬 Acquiring spectrum #{next_spectrum_id}...')

                spectrum_relpath = self._build_spectrum_relpath(current_spectrum_id)
                manufacturer_msa_path = os.path.join(self.sample_result_dir, spectrum_relpath)
                os.makedirs(os.path.dirname(manufacturer_msa_path), exist_ok=True)

                next_spectrum_id += 1
                total_counts = self._acquire_spectrum(
                    x,
                    y,
                    spectrum_id=current_spectrum_id,
                    msa_file_path=manufacturer_msa_path,
                )

                # Immediately update ledger after each successful acquisition
                try:
                    from autoemx.config.ledger_schemas import SpectrumEntry, SampleLedger
                except ImportError:
                    pass  # Already imported at top
                # Build SpectrumEntry directly from the pointer file and current acquisition metadata
                acq_details = AcquisitionDetails(
                    frame_id=str(frame_ID).strip() if frame_ID is not None and str(frame_ID).strip() else None,
                    particle_id=self.particle_cntr if self.particle_cntr is not None else None,
                    spot_coordinates=self._build_spot_coordinates(
                        pixel_x=xy_center[0],
                        pixel_y=xy_center[1],
                    ),
                )
                spectrum_entry = self._build_spectrum_entry_from_pointer_file(
                    Path(manufacturer_msa_path),
                    acquisition_details_by_id={current_spectrum_id: acq_details},
                )
                # Load or create ledger
                ledger_path = self._get_ledger_path()
                try:
                    ledger = self._load_existing_ledger()
                except Exception:
                    ledger = None
                if ledger is None:
                    ledger = SampleLedger(
                        sample_id=self.sample_id,
                        sample_path=os.path.abspath(self.sample_result_dir),
                        configs=self._build_ledger_configs(),
                        spectra=[],
                        quantifications=[],
                        active_quant=None,
                    )
                ledger.spectra.append(spectrum_entry)
                ledger.to_json_file(ledger_path)
                ingested_spectrum_ids.add(str(spectrum_entry.spectrum_id))
                n_spectra_collected += 1

                # Contamination check: skip quantification if counts are too low (only at first measurement spot)
                if i==0 and self.sample_cfg.is_particle_acquisition:
                    if total_counts < 0.95 * self.measurement_cfg.target_acquisition_counts:
                        if quantify:
                            self.spectra_quant_records.append(None)
                        if self.verbose:
                            logger.warning('⚠️ Current particle is unlikely to be part of the sample.\nSkipping to the next particle.')
                            logger.info('ℹ️ Increase measurement_cfg.max_acquisition_time if this behavior is undesired.')
                        break
            
            # Save image of particle, with ID of acquired XSp spots
            if latest_spot_id is not None:
                # Prepare save path
                par_cntr_str = f"_par{self.particle_cntr}" if self.particle_cntr is not None else ''
                filename = f"{self.sample_id}{par_cntr_str}_fr{frame_ID}_xyspots"
                # Construct annotation dictionary
                im_annotations = []
                for i, xy_coords in enumerate(spots_xy_list):
                    # Skip if latest_spot_id is None or i is out of range
                    if latest_spot_id is None or i > latest_spot_id:
                        break

                    xy_center = (int(xy_coords[0]), int(xy_coords[1]))
                    if xy_center is None:
                        continue
                    
                    im_annotations.append({
                        cnst.ANNOTATION_TEXT_KEY: (
                            str(next_spectrum_id - 1 - latest_spot_id + i),
                            (xy_center[0] - 30, xy_center[1] - 15)
                        ),
                        cnst.ANNOTATION_CIRCLE_KEY: (10, xy_center, -1)
                    })
                # Save image with annotations
                self.EM_controller.save_frame_image(filename, im_annotations = im_annotations)
                
            if quantify:
                self._fit_and_quantify_spectra(interrupt_fits_bad_spectra = interrupt_fits_bad_spectra)
    
        return next_spectrum_id, success
    
    
    def _fit_and_quantify_spectra(
            self,
            quantify: bool = True,
            force_requantification: bool = False,
            requantify_only_unquantified_spectra: bool = False,
            interrupt_fits_bad_spectra: bool = True,
            num_CPU_cores: Optional[int] = None,
        ) -> None:
            """
            Fit and (optionally) quantify all collected spectra.

            Parameters
            ----------
            quantify : bool
                If False, only fits spectra (used for experimental standards). Default True.
            force_requantification : bool
                Create a new quantification run and reprocess every spectrum, even when
                an identical run already exists.
            requantify_only_unquantified_spectra : bool
                Reuse the latest matching run but reprocess only spectra that have no
                composition result (never quantified, or previously skipped/flagged).
                Overwrites the prior skipped record in the ledger. Ignored when
                force_requantification=True.
            interrupt_fits_bad_spectra : bool
                Controls early-exit behaviour during iterative spectral fitting.

                If ``True`` (default), the fit is aborted mid-iteration when poor fit quality,
                excessive analytical error, or excessive X-ray absorption is detected.
                The spectrum is stored with ``QuantificationDiagnostics.interrupted=True``
                and no composition is saved.

                If ``False``, early-exit is disabled.  Any spectrum from the active
                quantification run whose record has ``interrupted=True`` is re-quantified
                and its ledger record is overwritten with the new result.
            num_CPU_cores : Optional[int]
                Number of CPU cores for parallel fitting.
                None uses half of available cores.
            """
            
            _ledger = self._load_or_create_ledger()
            tot_spectra_collected = len(_ledger.spectra) if _ledger is not None else 0
            _ledger_spectra = list(_ledger.spectra) if _ledger is not None else []

            available_cores = os.cpu_count() or 1
            default_cores = max(1, available_cores // 2)
            requested_cores = min(
                (num_CPU_cores if num_CPU_cores is not None else default_cores),
                available_cores,
            )
            _n_cores, parallel_backend, use_parallel = _resolve_parallel_mode(
                requested_cores=requested_cores,
            )

            try:
                if use_parallel:
                    logger.info(
                        "⚙️ Execution mode: parallel (backend=%s, cores=%d)",
                        parallel_backend,
                        _n_cores,
                    )
                else:
                    logger.info("⚙️ Execution mode: serial (cores=1)")
            except BaseException:
                pass

            quant_worker_payloads: List[Dict[str, Any]] = []

            if quantify:
                # Always bootstrap/sync ledger before quantification cycles.
                self._ensure_current_quantification_run(force_new=force_requantification)
                self._persist_current_quantification_config()
                self._ensure_quant_tracking_length(tot_spectra_collected)
                if force_requantification:
                    indices_to_process = list(range(tot_spectra_collected))
                elif requantify_only_unquantified_spectra:
                    self._sync_existing_quantification_from_ledger()
                    indices_to_process = [
                        i for i in range(tot_spectra_collected)
                        if self.spectra_quant_records[i] is None
                        or self.spectra_quant_records[i].composition_atomic_fractions is None
                    ]
                else:
                    self._sync_existing_quantification_from_ledger()
                    indices_to_process = [
                        i for i in range(tot_spectra_collected)
                        if self.spectra_quant_records[i] is None
                        or (
                            not interrupt_fits_bad_spectra and getattr(self.spectra_quant_records[i].diagnostics, 'interrupted', False)
                        )
                    ]
                    # Always skip spectra that have been tried (quant_record exists), regardless of interrupted, unless interrupt_fits_bad_spectra is False
                    if interrupt_fits_bad_spectra:
                        indices_to_process = [
                            i for i in indices_to_process
                            if self.spectra_quant_records[i] is None
                        ]
            else:
                # Fitting-only (experimental standards) path
                self._ensure_current_quantification_run(force_new=force_requantification)
                self._persist_current_quantification_config()
                self._ensure_quant_tracking_length(tot_spectra_collected)
                self._sync_existing_quantification_from_ledger()
                if force_requantification:
                    indices_to_process = list(range(tot_spectra_collected))
                else:
                    indices_to_process = [
                        i for i in range(tot_spectra_collected)
                        if self.spectra_quant_records[i] is None
                    ]
                
                # Build fitting payloads (same structure as quantification payloads)
                for i in indices_to_process:
                    spectrum_entry = _ledger_spectra[i]
                    background_abs = None
                    if self.quant_cfg.use_instrument_background and spectrum_entry.instrument_background_relpath:
                        background_abs = str(Path(self.sample_result_dir, spectrum_entry.instrument_background_relpath))
                    quant_worker_payloads.append({
                        'spectrum_index': i,
                        'pointer_abs': str(Path(self.sample_result_dir, spectrum_entry.spectrum_relpath)),
                        'background_abs': background_abs,
                        'sp_collection_time': (
                            float(spectrum_entry.live_acquisition_time)
                            if spectrum_entry.live_acquisition_time is not None else 1.0
                        ),
                        'quantification_id': int(self.current_quantification_id),
                        'energy_vals': list(map(float, self.energy_vals)),
                        'sp_start': int(self.sp_start),
                        'sp_end': int(self.sp_end),
                        'target_acquisition_counts': int(self.measurement_cfg.target_acquisition_counts),
                        'min_bckgrnd_cnts': self.clustering_cfg.min_bckgrnd_cnts,
                        'microscope_id': self.microscope_cfg.ID,
                        'measurement_type': self.measurement_cfg.type,
                        'measurement_mode': self.measurement_cfg.mode,
                        'det_ch_offset': float(self.det_ch_offset),
                        'det_ch_width': float(self.det_ch_width),
                        'beam_energy_keV': float(self.measurement_cfg.beam_energy_keV),
                        'emergence_angle': float(self.measurement_cfg.emergence_angle),
                        'all_els_sample': list(self.all_els_sample),
                        'detectable_els_substrate': list(self.detectable_els_substrate),
                        'sample_w_frs': None,  # Not used for fitting-only
                        'apply_geom_factors': bool(self._apply_geom_factors),
                        'max_undetectable_w_fr': float(self.undetectable_an_er),
                        'fit_tolerance': float(self.quant_cfg.fit_tolerance),
                        'standards_dict': self.standards_dict,
                        'interrupt_fits_bad_spectra': bool(interrupt_fits_bad_spectra),
                        'quantify': False,
                    })

            n_spectra_to_quant = len(indices_to_process)
        
            if self.verbose and n_spectra_to_quant > 0:
                print_single_separator()
                quant_str = "quantification" if quantify else "fitting"
                logger.info(f"▶️ Starting {quant_str} of {n_spectra_to_quant} spectra on up to {_n_cores} cores")
        
            if quantify:
                indices_to_process_set = set(indices_to_process)
                for i in range(tot_spectra_collected):
                    if i not in indices_to_process_set:
                        quant_record = self.spectra_quant_records[i]
                        if quant_record is not None and not force_requantification:
                            if self.verbose:
                                print_single_separator()
                                logger.info(
                                    "⏭️ Spectrum #%d/%d already quantified. Skipping...",
                                    i + 1,
                                    tot_spectra_collected,
                                )
                        continue

                    spectrum_entry = _ledger_spectra[i]
                    background_abs = None
                    if self.quant_cfg.use_instrument_background and spectrum_entry.instrument_background_relpath:
                        background_abs = str(Path(self.sample_result_dir, spectrum_entry.instrument_background_relpath))
                    quant_worker_payloads.append({
                        'spectrum_index': i,
                        'pointer_abs': str(Path(self.sample_result_dir, spectrum_entry.spectrum_relpath)),
                        'background_abs': background_abs,
                        'sp_collection_time': (
                            float(spectrum_entry.live_acquisition_time)
                            if spectrum_entry.live_acquisition_time is not None else 1.0
                        ),
                        'quantification_id': int(self.current_quantification_id),
                        'energy_vals': list(map(float, self.energy_vals)),
                        'sp_start': int(self.sp_start),
                        'sp_end': int(self.sp_end),
                        'target_acquisition_counts': int(self.measurement_cfg.target_acquisition_counts),
                        'min_bckgrnd_cnts': self.clustering_cfg.min_bckgrnd_cnts,
                        'microscope_id': self.microscope_cfg.ID,
                        'measurement_type': self.measurement_cfg.type,
                        'measurement_mode': self.measurement_cfg.mode,
                        'det_ch_offset': float(self.det_ch_offset),
                        'det_ch_width': float(self.det_ch_width),
                        'beam_energy_keV': float(self.measurement_cfg.beam_energy_keV),
                        'emergence_angle': float(self.measurement_cfg.emergence_angle),
                        'all_els_sample': list(self.all_els_sample),
                        'detectable_els_substrate': list(self.detectable_els_substrate),
                        'sample_w_frs': (
                            dict(self.sample_cfg.w_frs)
                            if self.sample_cfg.w_frs is not None else None
                        ),
                        'apply_geom_factors': bool(self._apply_geom_factors),
                        'max_undetectable_w_fr': float(self.undetectable_an_er),
                        'fit_tolerance': float(self.quant_cfg.fit_tolerance),
                        'standards_dict': self.standards_dict,
                        'interrupt_fits_bad_spectra': bool(interrupt_fits_bad_spectra),
                        'quantify': True,
                    })
        
            def _finalize_quant_result(idx, result, quant_record, quantification_time):
                quant_flag = quant_record.quant_flag
                comment = quant_record.comment
                fit_result = getattr(quant_record, 'fit_result', None)
                fit_was_not_run = (not quantify and fit_result is None)
                had_prior_record = self.spectra_quant_records[idx] is not None
                was_interrupted = (
                    had_prior_record
                    and getattr(self.spectra_quant_records[idx].diagnostics, 'interrupted', False)
                )

                self.spectra_quant_records[idx] = quant_record
                if self.verbose:
                    lines = [f" Spectrum #{idx}/{tot_spectra_collected - 1}:"]
                    if fit_was_not_run:
                        if comment:
                            lines.append(f"  {comment}")
                    if quantify and result is None:
                        if comment:
                            lines.append(f"  {comment}")
                    elif not fit_was_not_run:
                        if quantify:
                            for el, at_fr in result[cnst.COMP_AT_FR_KEY].items():
                                lines.append(f"  {el} at%: {at_fr * 100:.2f}%")
                            lines.append(f"  Analytical error: {result[cnst.AN_ER_KEY] * 100:.2f}%")
                        else:
                            fitted_peaks = getattr(fit_result, 'fitted_peaks', None)
                            standard_elements = set(getattr(getattr(self, 'exp_stds_cfg', None), 'w_frs', {}) or {})
                            ref_line_priority = {
                                line: rank for rank, line in enumerate(XSp_Quantifier.xray_quant_ref_lines)
                            }
                            pb_by_peak: Dict[str, tuple[int, str, float]] = {}

                            if fitted_peaks:
                                for _, peak in fitted_peaks.items():
                                    el = getattr(peak, 'element', None)
                                    line = getattr(peak, 'line', None)
                                    pb_ratio = getattr(peak, 'pb_ratio', None)

                                    if el is None or line is None or pb_ratio is None:
                                        continue
                                    if line not in ref_line_priority:
                                        continue
                                    if standard_elements and el not in standard_elements:
                                        continue

                                    peak_label = f"{el}_{line}"
                                    pb_by_peak[peak_label] = (
                                        ref_line_priority[line],
                                        peak_label,
                                        float(pb_ratio),
                                    )

                            if pb_by_peak:
                                lines.append("  Measured PB:")
                                sorted_peaks = sorted(pb_by_peak.values(), key=lambda x: (x[0], x[1]))
                                for _, peak_label, pb_ratio in sorted_peaks:
                                    lines.append(f"  {peak_label}: {pb_ratio:.1f}")
                            else:
                                lines.append("  Fit completed")
                        if quant_flag not in (None, 0) and comment:
                            lines.append(f"  {comment}")
                        if quantification_time is not None:
                            lines.append(f"  Quantification took {quantification_time:.2f} s")
                    lines.append("")  # Add an extra line for spacing
                    print_single_separator("\n".join(lines))

                if quant_record is not None:
                    overwrite = (requantify_only_unquantified_spectra and had_prior_record) or was_interrupted
                    self._persist_quantification_record(idx, quant_record, overwrite=overwrite)
            
            try:
                if not use_parallel:
                    for payload in quant_worker_payloads:
                        idx, result, quant_record, quantification_time = _quantify_spectrum_worker(payload)
                        _finalize_quant_result(idx, result, quant_record, quantification_time)
                else:
                    # Stream completed tasks back to the main process so progress is visible in real time.
                    try:
                        completed = Parallel(
                            n_jobs=_n_cores,
                            backend=parallel_backend,
                            return_as='generator_unordered',
                        )(
                            delayed(_quantify_spectrum_worker)(payload) for payload in quant_worker_payloads
                        )
                    except TypeError:
                        # Compatibility fallback for joblib versions without return_as.
                        completed = Parallel(n_jobs=_n_cores, backend=parallel_backend)(
                            delayed(_quantify_spectrum_worker)(payload) for payload in quant_worker_payloads
                        )

                    for idx, result, quant_record, quantification_time in completed:
                        _finalize_quant_result(idx, result, quant_record, quantification_time)
            except Exception as e:
                logger.warning(
                    f"⚠️ Parallel spectrum processing failed ({type(e).__name__}: {e}), "
                    "falling back to sequential execution."
                )
                for payload in quant_worker_payloads:
                    idx, result, quant_record, quantification_time = _quantify_spectrum_worker(payload)
                    _finalize_quant_result(idx, result, quant_record, quantification_time)


    def _persist_resolved_k_on_active_clustering_config(self, resolved_k: Optional[int]) -> None:
        """Persist resolved k for the active clustering config of the active quantification run."""
        if resolved_k is None:
            return

        resolved_k_int = int(resolved_k)
        self.clustering_cfg.k_resolved = resolved_k_int

        if self.current_quant_config is not None:
            active_clustering_analysis = self.current_quant_config.get_active_clustering_analysis()
            active_clustering_config = (
                active_clustering_analysis.config if active_clustering_analysis is not None else None
            )
            if active_clustering_config is not None and active_clustering_config.k_resolved != resolved_k_int:
                active_clustering_config.k_resolved = resolved_k_int
                self._persist_current_quantification_config()
            return

        ledger = self._load_or_create_ledger()
        active_quant_config = self._get_active_quantification_config(ledger)
        if active_quant_config is None:
            return

        active_clustering_analysis = active_quant_config.get_active_clustering_analysis()
        active_clustering_config = (
            active_clustering_analysis.config if active_clustering_analysis is not None else None
        )
        if active_clustering_config is None or active_clustering_config.k_resolved == resolved_k_int:
            return

        active_clustering_config.k_resolved = resolved_k_int
        ledger.to_json_file(self._get_ledger_path())


    def _persist_clustering_result_on_active_analysis(self, clustering_result: ClusteringResult) -> None:
        """Persist the latest clustering result on the active clustering analysis entry."""
        ledger = self._load_or_create_ledger()
        active_quant_config = self._get_active_quantification_config(ledger)
        if active_quant_config is None:
            return

        # Ensure the current runtime clustering settings map to a concrete analysis entry
        # on the same ledger object being persisted.
        self._ensure_current_clustering_run(active_quant_config)

        active_quant_config.set_active_clustering_result(clustering_result.model_copy(deep=True))

        if (
            self.current_quant_config is not None
            and self.current_quant_config.quantification_id == active_quant_config.quantification_id
        ):
            self.current_quant_config.set_active_clustering_result(clustering_result.model_copy(deep=True))

        ledger.to_json_file(self._get_ledger_path())


    #%% Data compositional analysis
    # =============================================================================     
    def analyse_data(self, max_analytical_error_percent, k=None, compute_k_only_once=False):
        """
        Analyse quantified spectra, perform clustering, assign candidate phases and mixtures, and save results.
    
        This function orchestrates the workflow:
          1. Selects good compositions for clustering.
          2. Prepares DataFrames for clustering.
          3. Determines the optimal number of clusters (k).
          4. Runs clustering and computes cluster statistics.
          5. Assigns candidate phases and detects mixtures.
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
            Minimum confidence among assigned candidate phases.
        """
        # Ensure in-memory state is fully synced from the ledger before analysis.
        # This is essential when analyse_data is called without a preceding
        # run_quantification (i.e. analysis-only path): we need the typed
        # QuantificationResult records populated in memory.
        if self.sample_result_dir is not None:
            existing_ledger = self._load_or_create_ledger()
            # Resolve current_quantification_id from the ledger when not already set
            # (analysis-only callers skip run_quantification so it is never assigned).
            if self.current_quantification_id is None and existing_ledger is not None:
                if existing_ledger.active_quant is not None:
                    self.current_quantification_id = existing_ledger.active_quant

            # Analysis-only runs must also register the current clustering inputs on the
            # active quantification config so clustering history/index stays in sync.
            active_quant_config = self._get_active_quantification_config(existing_ledger)
            if active_quant_config is not None:
                self.current_quant_config = active_quant_config
                self.current_quantification_id = active_quant_config.quantification_id
                self._ensure_current_clustering_run(self.current_quant_config)
                self._persist_current_quantification_config()
                self._apply_active_clustering_config(self.current_quant_config)

            tot = len(existing_ledger.spectra) if existing_ledger is not None else 0
            if tot > 0:
                self._ensure_quant_tracking_length(tot)
                self._sync_existing_quantification_from_ledger()

        # 1. Make analysis directory to save results
        self._make_analysis_dir()
    
        self._save_analysis_config_summary()

        # 2. Select compositions to use for clustering
        if max_analytical_error_percent is not None:
            max_analytical_error = max_analytical_error_percent / 100
        else:
            max_analytical_error = max_analytical_error_percent
        (compositions_list_at, compositions_list_w, unused_compositions_list,
         df_indices, n_datapts) = ClusteringModule._select_good_compositions(self, max_analytical_error)
        n_datapts_used = len(compositions_list_at)
    
        if n_datapts_used < 5:
            print_single_separator()
            logger.warning(f"⚠️ Only {n_datapts_used} spectra were considered 'good', but a minimum of 5 data points are required for clustering.")
            # Print additional messages with how many spectra were discarded for which reason
            self._report_n_discarded_spectra(n_datapts, max_analytical_error)
            # Save Composition.csv file anyways
            self._save_analysis_summary(None, None)
            return False, 0, 0  # zeroes are placeholders
    
        if self.verbose:
            print_single_separator()
            logger.info('ℹ️ Spectra selection:')
            logger.info(f"ℹ️ {n_datapts_used} data points are used, out of {n_datapts} collected spectra.")
            self._report_n_discarded_spectra(n_datapts, max_analytical_error)

        # 3. Prepare DataFrames for clustering
        compositions_df, compositions_df_other_fr = ClusteringModule._prepare_composition_dataframes(self, compositions_list_at, compositions_list_w)

        # 4. Perform clustering
        if k is None:
            if self.clustering_cfg.k_finding_method == "forced":
                k = self.clustering_cfg.k_forced
                if k is None:
                    raise ValueError("k_finding_method='forced' requires k_forced to be set")
            else:
                # Recompute k for non-forced methods on each analysis run.
                k = None
        if self.clustering_cfg.method == 'kmeans':
            k = ClusteringModule._find_optimal_k(self, compositions_df, k, compute_k_only_once)
            self._persist_resolved_k_on_active_clustering_config(k)
            kmeans, labels, sil_score = ClusteringModule._run_kmeans_clustering(self, k, compositions_df)
            centroids = kmeans.cluster_centers_
            wcss = kmeans.inertia_
        elif self.clustering_cfg.method == 'dbscan':
            # labels, num_labels = self._get_clustering_dbscan(compositions_df)
            logger.warning('⚠️ Clustering via DBSCAN is not implemented yet')
            return False, 0, 0  # zeroes are placeholders
    
        # 5. Compute cluster statistics
        (wcss_per_cluster, rms_dist_cluster, rms_dist_cluster_other_fr, 
         n_points_per_cluster, els_std_dev_per_cluster, els_std_dev_per_cluster_other_fr,
         centroids_other_fr, max_cl_rmsdist) = ClusteringModule._compute_cluster_statistics(
            self, compositions_df, compositions_df_other_fr, centroids, labels
        )
    
        # 6. Assign candidate phases
        min_conf, max_raw_confs, refs_assigned_df = ReferenceMatchingModule._assign_reference_phases(self, centroids, rms_dist_cluster)
    
        # 7. Assign mixtures
        if self.clustering_cfg.do_matrix_decomposition:
            clusters_assigned_mixtures = ReferenceMatchingModule._assign_mixtures(
                self, k, labels, compositions_df, rms_dist_cluster, max_raw_confs, n_points_per_cluster
            )
        else:
            clusters_assigned_mixtures = []
    
        # 8. Save and store results
        self._save_analysis_summary(labels, df_indices)
    
        clustering_artifacts = ClusteringResult(
            centroids=self._to_finite_float_matrix(centroids),
            els_std_dev_per_cluster=self._to_finite_float_matrix(els_std_dev_per_cluster),
            centroids_other_fr=self._to_finite_float_matrix(centroids_other_fr),
            els_std_dev_per_cluster_other_fr=self._to_finite_float_matrix(els_std_dev_per_cluster_other_fr),
            n_points_per_cluster=n_points_per_cluster,
            wcss_per_cluster=self._to_finite_float_list(wcss_per_cluster),
            rms_dist_cluster=self._to_finite_float_list(rms_dist_cluster),
            rms_dist_cluster_other_fr=self._to_finite_float_list(rms_dist_cluster_other_fr),
            refs_assigned_rows=self._sanitize_json_compatible(refs_assigned_df.to_dict(orient='records')),
            wcss=self._to_finite_float(wcss),
            sil_score=self._to_finite_float(sil_score),
            tot_n_points=n_datapts,
            clusters_assigned_mixtures=self._sanitize_json_compatible(clusters_assigned_mixtures),
        )
        self._save_result_and_stats(clustering_artifacts)
        self._persist_clustering_result_on_active_analysis(clustering_artifacts)
    
        # 9. Save plots
        if self.plot_cfg.save_plots:
            PlottingModule._save_plots(self, kmeans, compositions_df, centroids, labels, els_std_dev_per_cluster, unused_compositions_list)
    
        return True, max_cl_rmsdist, min_conf
    
    

    #%% Run algorithms
    # =============================================================================     
    def run_collection_and_quantification(
        self,
        quantify: bool = True,
        interrupt_fits_bad_spectra: bool = True,
    ) -> Tuple[bool, bool]:
        """
        Perform iterative collection (and optional quantification) of spectra, followed by phase analysis and convergence check.
    
                This method:
                    - Iteratively collects spectra in batches (and quantifies them if `quantify` is True).
                    - After each batch, persists acquisitions to the ledger and (if quantification is enabled) performs phase analysis (clustering).
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
        - Acquisition and quantification records are persisted incrementally to the ledger to prevent data loss.
        - Prints a summary and processing times at the end.
        """
        # Resume from any previously acquired spectra so file names and particle IDs
        # are always monotonically increasing across restarted acquisitions.
        _had_ledger_before = self._load_existing_ledger() is not None
        _existing_ledger = self._load_or_create_ledger()
        if _existing_ledger is not None and _existing_ledger.spectra:
            tot_n_spectra = len(_existing_ledger.spectra)
            _numeric_ids = [
                int(e.spectrum_id)
                for e in _existing_ledger.spectra
                if e.spectrum_id is not None and str(e.spectrum_id).isdigit()
            ]
            next_spectrum_id = (max(_numeric_ids) + 1) if _numeric_ids else 0
            _particle_ids = [
                e.acquisition_details.particle_id
                for e in _existing_ledger.spectra
                if e.acquisition_details is not None
                and e.acquisition_details.particle_id is not None
            ]
            _particle_id_offset = 0
            self.particle_cntr = -1
            try:
                if _particle_ids:
                    # Controller particle counter restarts from 1 on each run.
                    # Keep offset equal to the last persisted particle id so resume starts at last+1.
                    _particle_id_offset = max(int(p) for p in _particle_ids)
                    self.particle_cntr = _particle_id_offset
                else:
                    inferred_particle_id = self._infer_max_particle_id_from_saved_images()
                    if inferred_particle_id is not None:
                        _particle_id_offset = inferred_particle_id
                        self.particle_cntr = inferred_particle_id
            except (ValueError, TypeError):
                _particle_id_offset = 0
                self.particle_cntr = -1
            if self.verbose:
                if (
                    not _had_ledger_before
                    and getattr(self, "_last_acq_ledger_created", False)
                    and getattr(self, "_last_acq_ledger_ingested_spectra_count", 0) > 0
                ):
                    logger.info(
                        "ℹ️ Restart detected %d pre-existing spectrum file/s and ingested them into a new ledger "
                        "(assuming the same acquisition configurations).",
                        getattr(self, "_last_acq_ledger_ingested_spectra_count", 0),
                    )
                elif _had_ledger_before:
                    logger.info(
                        "ℹ️ Existing ledger detected; acquisition will resume from the next available spectrum and particle IDs."
                    )
            logger.info(
                "ℹ️ %d previously acquired spectrum/spectra detected. "
                "New spectra will be appended starting from index %d. "
                "Delete the sample folder before starting if a fresh acquisition is desired: %s",
                tot_n_spectra,
                next_spectrum_id,
                self.sample_result_dir,
            )
        else:
            tot_n_spectra = 0
            next_spectrum_id = 0
            _particle_id_offset = 0

        next_particle_id = _particle_id_offset + 1 if self.sample_cfg.is_particle_acquisition else "n/a"
        if next_spectrum_id > 0 and self.verbose:
            logger.info(
                "ℹ️ Resume point resolved: next spectrum ID=%d, next particle ID=%s.",
                next_spectrum_id,
                next_particle_id,
            )

        max_n_sp_per_iter = 10  # Max spectra to collect per iteration (for saving in between)
        n_spectra_collected_this_session = 0  # Track new spectra collected in this session
        is_converged = False
        is_analysis_successful = False
        is_acquisition_successful = True
        is_exp_std_measurement = self.exp_stds_cfg.is_exp_std_measurement
        is_spectral_quant = quantify and not is_exp_std_measurement
        tot_spectra_to_collect = self.max_n_spectra if not is_exp_std_measurement else self.min_n_spectra
        n_spectra_to_collect = min(max_n_sp_per_iter, max(0, tot_spectra_to_collect - tot_n_spectra), self.min_n_spectra)

        if is_spectral_quant:
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
            logger.info(f"▶️ Starting acquisition{quant_str} of {tot_spectra_to_collect} spectra.")
        
        while tot_n_spectra < self.max_n_spectra:
            if self.verbose:
                print_double_separator()
                logger.info(f"🔬 Collecting{quant_str} {n_spectra_to_collect} spectra...")
    
            # Collect the next batch of spectra (and quantify if requested)
            n_spectra_before = tot_n_spectra
            next_spectrum_id, is_acquisition_successful = self._collect_spectra(
                n_spectra_to_collect,
                n_tot_sp_collected=next_spectrum_id,
                quantify=is_spectral_quant,
                particle_id_offset=_particle_id_offset,
                interrupt_fits_bad_spectra=interrupt_fits_bad_spectra
            )
            refreshed_ledger = self._load_or_create_ledger()
            tot_n_spectra = len(refreshed_ledger.spectra) if refreshed_ledger is not None else 0
            n_spectra_collected_this_session += (tot_n_spectra - n_spectra_before)
    
            if self.verbose:
                print_single_separator()
                logger.info(f"✅ {tot_n_spectra} total spectra collected.")
            
            if quantify and tot_n_spectra > 0:
                if is_exp_std_measurement:
                    # Fit spectra and check if target number of good spectra has been collected
                    n_valid_spectra_collected, is_analysis_successful, is_converged = StandardsModule._evaluate_exp_std_fit(self, tot_n_spectra)
                else:
                    # Perform clustering analysis and check for convergence
                    n_valid_spectra_collected, is_analysis_successful, is_converged = self._evaluate_clustering_convergence(tot_n_spectra, n_spectra_to_collect)
                    
                if is_converged:
                    break
            else:
                n_valid_spectra_collected = tot_n_spectra
                
            # Stop if no more particles are available on the sample
            if not is_acquisition_successful:
                if self.verbose:
                    logger.warning("⚠️ Acquisition interrupted.")
                    if self.sample_cfg.is_particle_acquisition:
                        logger.warning(f'⚠️ Not enough particles were found on the sample to collect all {tot_spectra_to_collect} spectra.')
                    elif self.sample_cfg.is_grid_acquisition:
                        logger.warning(f'⚠️ The specified spectrum spacing did not allow to collect all {tot_spectra_to_collect} spectra.\n'
                              "Change spacing in bulk_meas_cfg to collect more spectra.")
                break

                # Collect additional spectra in next iteration
            n_spectra_to_collect = min(
                max_n_sp_per_iter,
                tot_spectra_to_collect - n_valid_spectra_collected,
                self.min_n_spectra
            )
    
        print_double_separator()
        logger.info('ℹ️ Sample ID: %s', self.sample_id)
        if self.sample_cfg.is_particle_acquisition:
            if n_spectra_collected_this_session > 0:
                par_str = f' over {self.particle_cntr} particles'
            else:
                par_str = ' (no new spectra collected in this session)'
        else:
            par_str = ''
        if n_spectra_collected_this_session > 0:
            logger.info(f'✅ {n_spectra_collected_this_session} new spectra collected (total: {tot_n_spectra}){par_str}.')
        else:
            logger.info(f'✅ {tot_n_spectra} total spectra available{par_str}.')
        process_time = (time.time() - self.start_process_time) / 60
        logger.info(f'✅ Total compositional analysis time: {process_time:.1f} min')
        print_single_separator()
    
        if is_spectral_quant:
            if is_analysis_successful:
                if is_converged:
                    logger.info('✅ Clustering converged to small errors. All phases identified with confidence higher than 0.8.')
                else:
                    logger.warning('⚠️ Phases could not be identified with confidence higher than 0.8.')
    
                self.print_results()
    
            elif not is_acquisition_successful:
                logger.warning('⚠️ This did not allow to determine which phases are present in the sample.')
            else:
                logger.warning(f'⚠️ Phases could not be identified with the allowed maximum of {self.max_n_spectra} collected spectra.')
                is_analysis_successful = False
                is_converged = False
        else:
            is_analysis_successful = is_acquisition_successful
    
        return is_analysis_successful, is_converged
    
    
    def _evaluate_clustering_convergence(
        self,
        tot_n_spectra: int,
        n_spectra_to_collect: int
    ) -> Tuple[int, bool, bool]:
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
            logger.info(f"📊 Analysing phases after collection of {tot_n_spectra} spectra...")
    
        try:
            is_analysis_successful, max_cl_rmsdist, min_conf = self.analyse_data(
                self.clustering_cfg.max_analytical_error_percent,
                k=self.clustering_cfg.k_forced if self.clustering_cfg.k_finding_method == "forced" else None
            )
        except Exception as e:
            raise RuntimeError(f"Error during clustering analysis: {e}") from e
    
        # Default value in case convergence check is skipped
        is_converged = False
    
        if is_analysis_successful:
            if self.verbose:
                print_double_separator()
                logger.info("✅ Clustering analysis performed")
    
            # Check whether phase identification converged
            try:
                is_converged = self._is_comp_analysis_converged(max_cl_rmsdist, min_conf)
            except Exception as e:
                raise RuntimeError(f"Error while checking convergence: {e}") from e
    
            if tot_n_spectra >= self.min_n_spectra:
                if is_converged:
                    return tot_n_spectra, is_analysis_successful, is_converged
                elif self.verbose and n_spectra_to_collect > 0:
                    logger.warning("⚠️ Compositional analysis did not converge, more spectra will be collected.")
            elif tot_n_spectra >= self.max_n_spectra:
                logger.info(f"ℹ️ Maximum allowed number of {self.max_n_spectra} was acquired.")
            else:
                if self.verbose:
                    logger.info(f"ℹ️ Collecting additional spectra to reach minimum number of {self.min_n_spectra}.")
    
        elif self.verbose:
            logger.error("❌ Clustering analysis unsuccessful.")
            if n_spectra_to_collect > 0:
                logger.info("ℹ️ More spectra will be collected.")
    
        return tot_n_spectra, is_analysis_successful, is_converged
    
    
    def _is_comp_analysis_converged(
        self,
        rms_dist: float,
        min_conf: Optional[float]
    ) -> bool:
        """
        Determine if the clustering analysis has converged based on cluster statistics.
        Used when collecting and quantifying spectra in real time.
    
        Convergence criteria:
          - If no candidate phases are present or assigned (min_conf is None), require cluster RMS point-to-centroid distance to be  < 2.5%.
          - If candidate phases are assigned, require minimum confidence > 0.8 and cluster standard deviation < 3%.
    
        Parameters
        ----------
        rms_dist : float
            Maximum RMS point-to-centroid distance among clusters (fractional units, e.g., 0.025 for 2.5%).
        min_conf : float or None
            Minimum confidence among all clusters assigned to candidate phases. If None, no references are assigned.
    
        Returns
        -------
        is_converged : bool
            True if convergence criteria are met, False otherwise.
    
        Notes
        -----
        - The thresholds are empirically determined for robust phase identification.
        """
        if min_conf is None:
            # No candidate phases present or assigned; require tighter cluster homogeneity
            is_converged = rms_dist < 0.025
        else:
            # Require high confidence and allow slightly larger within-cluster spread
            is_converged = (min_conf > 0.8) and (rms_dist < 0.03)
    
        return is_converged
    

    def run_quantification(
        self,
        force_requantification: bool = False,
        requantify_only_unquantified_spectra: bool = False,
        interrupt_fits_bad_spectra: bool = True,
        num_CPU_cores: Optional[int] = None,
    ) -> None:
        """
        Perform quantification of all collected spectra and save the results.
    
        Parameters
        ----------
        force_requantification : bool, optional
            If True, quantifies all spectra again even when the same quantification
            settings were already used before (creates a new run).
        requantify_only_unquantified_spectra : bool, optional
            Reuse the latest matching run but reprocess only spectra with no composition
            result (never quantified or previously skipped/flagged). Overwrites prior
            skipped records. Ignored when force_requantification=True.
        interrupt_fits_bad_spectra : bool, optional
            Controls early-exit behaviour during iterative spectral fitting.

            If ``True`` (default), the fit is aborted mid-iteration when poor fit quality,
            excessive analytical error, or excessive X-ray absorption is detected.
            The spectrum is stored with ``QuantificationDiagnostics.interrupted=True``
            and no composition is saved.

            If ``False``, early-exit is disabled.  Any spectrum from the active
            quantification run whose record has ``interrupted=True`` is re-quantified
            and its ledger record is overwritten with the new result.
        num_CPU_cores : Optional[int], optional
            Number of CPU cores for parallel fitting (non-quantify path only).
            None uses half of available cores.
        """
        self._initialise_std_dict()
        self._fit_and_quantify_spectra(
            force_requantification=force_requantification,
            requantify_only_unquantified_spectra=requantify_only_unquantified_spectra,
            interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
            num_CPU_cores=num_CPU_cores,
        )
        # Quantification-only workflows also export Compositions.csv, so ensure
        # a deterministic analysis directory exists instead of falling back to
        # the sample root.
        self._make_analysis_dir()
        self._save_analysis_summary(None, None)
        
    
    def run_exp_std_collection(
        self,
        fit_during_collection: bool = True,
        update_std_library: bool = True
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
            logger.info(f"🔬 Experimental standard acquisition of {self.sample_id}")
        
        # Run collection and quantification (fitting optionally performed during collection)
        self._th_peak_energies = {} # Initialise
        self.run_collection_and_quantification(quantify=fit_during_collection)
        
        self.run_exp_std_fit(run_fitting=not fit_during_collection, update_std_library=update_std_library)


    def run_exp_std_fit(self, run_fitting: bool = True, force_refitting: bool = False, update_std_library: bool = True) -> None:
        """
        Fit experimental standards and optionally update the standards library.
        If run_fitting is False, only aggregate/save from existing records.
        """
        fit_results = StandardsModule._fit_stds_and_save_results(
            self,
            run_fitting=run_fitting,
            force_refitting=force_refitting
        )

        if fit_results is not None and fit_results.lines:
            StandardsModule._print_pb_summary(self, fit_results, sample_id=self.sample_id)
        else:
            logger.info("⚠️ No valid fitted peaks were produced for sample '%s'.", self.sample_id)

        if update_std_library and fit_results is not None and fit_results.lines:
            StandardsModule._update_standard_library(self, fit_results)
        
    #%% Save Data
    # =============================================================================
    def _save_result_and_stats(
        self,
        clustering_artifacts: ClusteringResult,
    ) -> None:
        """
        Save and store clustering results and statistics, including centroids, standard deviations, reference assignments, mixture assignments, and summary statistics.
    
        This method:
          - Constructs a DataFrame of cluster statistics and assignments.
          - Adds candidate phase and mixture assignments if available.
          - Saves the DataFrame to CSV and stores it as an attribute.
          - Saves general clustering information to a JSON file and stores it as an attribute.
    
        Parameters
        ----------
        clustering_artifacts : ClusteringResult
            Schema model bundling all clustering outputs required for serialization.
    
        Returns
        -------
        None
    
        Notes
        -----
                - The cluster DataFrame is saved as a transposed '<cnst.CLUSTERS_FILENAME>.csv' in the analysis directory
                    for easier manual reading.
        - General clustering info is stored as an attribute for later use and persisted in the ledger.
    
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
        centroids = clustering_artifacts.centroids
        els_std_dev_per_cluster = clustering_artifacts.els_std_dev_per_cluster
        centroids_other_fr = clustering_artifacts.centroids_other_fr
        els_std_dev_per_cluster_other_fr = clustering_artifacts.els_std_dev_per_cluster_other_fr
        n_points_per_cluster = clustering_artifacts.n_points_per_cluster
        wcss_per_cluster = clustering_artifacts.wcss_per_cluster
        rms_dist_cluster = clustering_artifacts.rms_dist_cluster
        rms_dist_cluster_other_fr = clustering_artifacts.rms_dist_cluster_other_fr
        refs_assigned_df = pd.DataFrame(clustering_artifacts.refs_assigned_rows)
        wcss = clustering_artifacts.wcss
        sil_score = clustering_artifacts.sil_score
        tot_n_points = clustering_artifacts.tot_n_points
        clusters_assigned_mixtures = clustering_artifacts.clusters_assigned_mixtures

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
        mixtures_df = ReferenceMatchingModule._build_mixtures_df(self, clusters_assigned_mixtures)
        clusters_df = pd.concat([clusters_df.reset_index(drop=True), mixtures_df.reset_index(drop=True)], axis=1)
    
        # Ensure the analysis directory exists
        try:
            os.makedirs(self.analysis_dir, exist_ok=True)
        except Exception as e:
            raise OSError(f"Could not create analysis directory '{self.analysis_dir}': {e}")
    
        # Save and store DataFrame
        clusters_csv_path = os.path.join(self.analysis_dir, cnst.CLUSTERS_FILENAME + '.csv')
        try:
            clusters_csv_df = clusters_df.transpose().copy()
            clusters_csv_df.columns = [f"Cluster_{i}" for i in range(len(clusters_df))]
            clusters_csv_df.to_csv(clusters_csv_path, index=True, header=True, index_label="Metric")
        except Exception as e:
            raise OSError(f"Could not write clusters DataFrame to '{clusters_csv_path}': {e}")
    
        self.clusters_df = clusters_df
    
        # Save general clustering info and store for printing
        now = datetime.now()

        def _serialize_config(config_obj: Any) -> Dict[str, Any]:
            if hasattr(config_obj, "model_dump"):
                return config_obj.model_dump(mode="json")
            return asdict(config_obj)
    
        # Gather configuration dataclasses as dictionaries
        cfg_dataclasses = {
            cnst.QUANTIFICATION_CFG_KEY: _serialize_config(self.quant_cfg),
            cnst.CLUSTERING_CFG_KEY: _serialize_config(self.clustering_cfg),
            cnst.PLOT_CFG_KEY: _serialize_config(self.plot_cfg),
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
        self.clustering_info = clustering_info
        
        
    def _save_analysis_summary(
        self,
        labels: Optional[List],
        df_indices: Optional[List],
    ) -> None:
        """
        Save quantification outputs to Compositions.csv.
    
        Parameters
        ----------
        labels : List
            List of cluster labels assigned to each spectrum (by index in df_indices).
        df_indices : List
            List of indices mapping labels to DataFrame rows.
        Returns
        -------
        data_df: pd.DataFrame | None
            DataFrame object containing the exported composition data.
            None if no spectrum is available for export.
    
        Notes
        -----
        - Only the Compositions.csv export path is retained.
        - Legacy Data.csv handling is intentionally removed.
    
        Raises
        ------
        OSError
            If the output directory cannot be created or files cannot be written.
        """
    
        ledger = self._load_or_create_ledger()
        ledger_spectra = list(ledger.spectra) if ledger is not None else []
        n_spectra = len(ledger_spectra)
    
        if n_spectra == 0:
            return None
        
        data_list = []
        for i in range(n_spectra):
            # Retrieve the typed QuantificationResult for this spectrum (None when not quantified)
            record = self.spectra_quant_records[i] if i < len(self.spectra_quant_records) else None
            coords = self._extract_spectrum_info(ledger_spectra[i], i)
            has_composition_result = record is not None and record.composition_atomic_fractions is not None
            has_result = has_composition_result

            if has_result:
                # Unpack spectral quantification results and convert from elemental fraction to % for readability
                atomic_comp = {el + cnst.AT_FR_DF_KEY: round(fr * 100, 2) for el, fr in record.composition_atomic_fractions.items()}
                weight_comp = {el + cnst.W_FR_DF_KEY: round(fr * 100, 2) for el, fr in record.composition_weight_fractions.items()}
                analytical_er = {cnst.AN_ER_DF_KEY: round(float(record.analytical_error) * 100, 2)}

                # Fit quality metrics (present only when a full fit was performed)
                r_squared = None
                redchi_sq = None
                if record.fit_result is not None:
                    r_squared = float(f"{record.fit_result.r_squared:.5f}") if record.fit_result.r_squared is not None else None
                    redchi_sq = float(f"{record.fit_result.reduced_chi_squared:.1f}") if record.fit_result.reduced_chi_squared is not None else None

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
                    cnst.R_SQ_KEY: r_squared,
                    cnst.REDCHI_SQ_KEY: redchi_sq,
                }
                    
                # Compose row of data to be saved
                data_row = {
                    **coords,
                    **meas_data
                }
            else:
                # Counts in this spectrum were too low or quantification was interrupted
                data_row = coords
            
            # Add comment and quantification flag columns when a quantification record exists.
            if record is not None:
                data_row[cnst.COMMENTS_DF_KEY] = record.comment
                data_row[cnst.QUANT_FLAG_DF_KEY] = record.quant_flag
    
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
    
        # Reorder columns to keep fit quality and flags toward the end.
        columns = data_df.columns.tolist()
        last_columns = [
            cnst.R_SQ_KEY, cnst.REDCHI_SQ_KEY,
            cnst.QUANT_FLAG_DF_KEY, cnst.COMMENTS_DF_KEY,
        ]
        remaining_columns = [col for col in columns if col not in last_columns]
        new_column_order = remaining_columns + [col for col in last_columns if col in columns]
        data_df = data_df[new_column_order]

        # Save compositions in Compositions.csv under the analysis directory only.
        if not hasattr(self, 'analysis_dir') or self.analysis_dir is None:
            self._make_analysis_dir()
        os.makedirs(self.analysis_dir, exist_ok=True)
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

        def _serialize_config(config_obj: Any) -> Dict[str, Any]:
            if hasattr(config_obj, "model_dump"):
                return config_obj.model_dump(mode="json")
            return asdict(config_obj)
    
        # Gather configuration dataclasses as dictionaries
        self.measurement_cfg = self.measurement_cfg.model_copy(
            update={
                "powder_meas_cfg": self.powder_meas_cfg,
                "bulk_meas_cfg": self.bulk_meas_cfg,
                "exp_stds_cfg": self.exp_stds_cfg if self.exp_stds_cfg.is_exp_std_measurement else None,
            }
        )
        cfg_dataclasses = {
            cnst.SAMPLE_CFG_KEY: _serialize_config(self.sample_cfg),
            cnst.MICROSCOPE_CFG_KEY: _serialize_config(self.microscope_cfg),
            cnst.MEASUREMENT_CFG_KEY: _serialize_config(self.measurement_cfg),
            cnst.SAMPLESUBSTRATE_CFG_KEY: _serialize_config(self.sample_substrate_cfg),
        }
        
        if is_XSp_measurement:
            cfg_dataclasses[cnst.QUANTIFICATION_CFG_KEY] = _serialize_config(self.quant_cfg)
    
        # Include dataclasses corresponding to measurement type
        if is_XSp_measurement and not self.exp_stds_cfg.is_exp_std_measurement:
            cfg_dataclasses[cnst.CLUSTERING_CFG_KEY] = self._runtime_clustering_cfg_payload(self.clustering_cfg)
            cfg_dataclasses[cnst.PLOT_CFG_KEY] = _serialize_config(self.plot_cfg)


    
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
    
        try:
            ledger = self._load_or_create_ledger()
            ledger.configs = LedgerConfigs.model_validate({
                cnst.MICROSCOPE_CFG_KEY: spectrum_collection_info[cnst.MICROSCOPE_CFG_KEY],
                cnst.SAMPLE_CFG_KEY: spectrum_collection_info[cnst.SAMPLE_CFG_KEY],
                cnst.MEASUREMENT_CFG_KEY: spectrum_collection_info[cnst.MEASUREMENT_CFG_KEY],
                cnst.SAMPLESUBSTRATE_CFG_KEY: spectrum_collection_info[cnst.SAMPLESUBSTRATE_CFG_KEY],
                cnst.PLOT_CFG_KEY: spectrum_collection_info.get(cnst.PLOT_CFG_KEY, self.plot_cfg.model_dump(mode="json")),
            })

            # Persist an initial quantification+clustering descriptor at acquisition start
            # so user-provided clustering inputs (e.g., sample['cnd'] -> ref_formulae)
            # are saved to ledger even before any quantification run is executed.
            if is_XSp_measurement and not self.exp_stds_cfg.is_exp_std_measurement and not ledger.quantifications:
                quantification_id = self._next_quantification_id(ledger)
                self.current_quant_config = self._build_quantification_config(
                    quantification_id=quantification_id,
                    sample_elements=list(self.all_els_sample),
                    substrate_elements=list(self.all_els_substrate),
                    options=self._build_quantification_options(),
                    reference_values_by_el_line=self._get_reference_values_by_el_line(active_quant_config=None),
                )
                self.current_quantification_id = quantification_id
                self._ensure_current_clustering_run(self.current_quant_config)
                ledger.append_quantification_config(self.current_quant_config)
                ledger.active_quant = quantification_id

            ledger.to_json_file(self._get_ledger_path())
        except Exception as e:
            raise OSError(f"Could not persist configurations to ledger '{self._get_ledger_path()}': {e}")
            
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
        logger.info("📊 Summary of Discarded Spectra")
        print("  → For details, see the 'Comments' column in Compositions.csv.")
    
        # Discarded due to low counts, insufficient background, or acquisition/fitting errors
        if self.n_sp_too_low_counts > 0:
            print(
                f"  • {self.n_sp_too_low_counts} spectra were discarded due to insufficient total counts, "
                f"background counts below the threshold ({self.clustering_cfg.min_bckgrnd_cnts}), "
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
            logger.warning("  ⚠️ Warning: More than 50% of spectra were flagged during quantification!")
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
            
        
    def print_results(self, n_cnd_to_print = 2, n_mix_to_print = 2) -> None:
        """
        Print a summary of clustering results, including clustering configuration, metrics,
        and a table of identified phases with elemental fractions, standard deviations,
        and reference/mixture assignments if present.
    
        The method:
          - Prints main clustering configuration and metrics.
          - Prints a table of phases, each with number of points, elemental fractions (with stddev),
            cluster stddev, WCSS, reference assignments, and mixture information if available.
        
        Parameters
        ----------
        n_cnd_to_print : int
            Max number of candidate phases and relative confidence scores to show. Candidates with scores
            close to 0 are not shown.
        n_mix_to_print : int
            Max number of candidate mixtures and relative confidence scores to show. Mixtures with scores
            close to 0 are not shown.
        
        Raises
        ------
        AttributeError
            If required attributes (clustering_info, clusters_df, etc.) are missing.
        KeyError
            If expected keys are missing from clustering_info or clusters_df.
        """
        # Print clustering info
        print_double_separator()
        logger.info(f"📊 Compositional analysis results for sample {self.sample_id}:")
        print_single_separator()
        try:
            logger.info('  Clustering method: %s', self.clustering_cfg.method)
            logger.info('  Clustering features: %s', self.clustering_cfg.features)
            logger.info('  k finding method: %s', self.clustering_cfg.k_finding_method)
            logger.info('  Number of clusters: %d', self.clustering_info[cnst.N_CLUST_KEY])
            logger.info('  WCSS (%%): %.2f', self.clustering_info[cnst.WCSS_KEY] * 10000)
            logger.info('  Silhouette score: %.2f', self.clustering_info[cnst.SIL_SCORE_KEY])
        except KeyError as e:
            raise KeyError(f"Missing key in clustering_info: {e}")
        except AttributeError as e:
            raise AttributeError(f"Missing attribute: {e}")
    
        # Print details on identified phases
        print_single_separator()
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
                    ref_keys_to_print = [key for i in range(1, n_cnd_to_print + 1) for key in [f'{cnst.CS_CND_DF_KEY}{i}', f'{cnst.CND_DF_KEY}{i}']]
                    ref_dict = {key: value for key, value in row.items() if key in ref_keys_to_print}
                    df_mod_to_print[-1].update(ref_dict)
                # Add mixtures to the printed report
                mix_keys_to_print = [key for i in range(1, n_mix_to_print + 1) for key in [f'{cnst.MIX_DF_KEY}{i}', f'{cnst.CS_MIX_DF_KEY}{i}']]
                mix_dict = {key: value for key, value in row.items() if key in mix_keys_to_print}
                df_mod_to_print[-1].update(mix_dict)
            # Set display options for float precision
            with pd.option_context('display.float_format', '{:,.2f}'.format):
                pd.set_option('display.max_columns', None)  # Display all columns
                phase_table = pd.DataFrame(df_mod_to_print).to_string()
                logger.info("Identified phases:\n%s", phase_table)
        except Exception as e:
            raise RuntimeError(f"Error printing phase results: {e}")
 