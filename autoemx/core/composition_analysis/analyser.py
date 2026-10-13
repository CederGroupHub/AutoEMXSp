#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EMXSp_Composition_Analyzer

Main class for automated compositional analysis of electron microscopy X-ray spectroscopy (EMXSp) data.

Can be run from Run_Acquisition.py, Run_Quantification.py and Run_Analysis.py.

Features:
- Structured configuration for microscope, sample, measurement, and analysis parameters.
- Automated acquisition and quantification of X-ray spectra at electron microscope.
- Filtering and clustering of compositional data.
- Phase identification, mixture analysis, and comprehensive results export.
- Utilities for plotting, saving, and reporting analysis results.

Example Usage
-------------
::

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

    # Analyzer of a sample already acquired, from its ledger (folder, ledger.json or loaded ledger)
    >>> analyzer = EMXSp_Composition_Analyzer.from_ledger('Results/sample_001')
    >>> analyzer.run_quantification()
    >>> analyzer.analyse_data(max_analytical_error_percent=5)


@author: Andrea
Created on Mon Jul 22 17:43:35 2024
"""

# Standard library imports
import os
import re
import time
import warnings
from pathlib import Path
from datetime import datetime
from dataclasses import asdict
from typing import Any, Optional, Tuple, List, Dict, Union
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
from sklearn.metrics import silhouette_score
 

# Project-specific imports
from autoemx.core.quantifier import XSp_Quantifier
from autoemx.core.quantifier.peak_overlaps import warn_peak_overlaps
from autoemx.core.quantifier.quant_flags import (
    check_fit_quant_validity,
    detect_prefit_spectrum_issues,
    format_min_total_counts_percent,
    merge_prefit_quant_flag,
    should_interrupt_prefit_issue,
)
from autoemx.core.fitter import DetectorResponseFunction
from autoemx.core.em_runtime.controller import EM_Controller
from autoemx.core.em_runtime.image_utilities import xsp_spot_annotation
from autoemx.core.em_runtime.xsp_spot_selection import XSpSpotSelectorCallback
from autoemx.core.em_runtime.sample_finder import EM_Sample_Finder
import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflt
from autoemx.utils.helper import (
    print_single_separator,
    print_double_separator,
    make_unique_path,
    composition_as_weight_dict,
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
    ParticleInfo,
    QuantificationConfig,
    QuantificationDiagnostics,
    QuantificationResult,
    SampleLedger,
)
from autoemx.config.ledger_io import (
    LedgerRuntimeConfigs,
    ledger_runtime_configs,
    load_sample_ledger,
    _build_spectrum_entry,
    _list_spectrum_pointer_files,
)
from autoemx.config.schema_models.clustering import ClusteringResult
from autoemx.config.schema_models.standards import EDSStandardsFile
from .clustering import ClusteringModule, count_results_without_bckgrnd_cnts, effective_quant_flag
from .ledger import LedgerModule
from .plotting import PlottingModule
from .reference_matching import ReferenceMatchingModule
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


def _identification_settings(worker_payload: Dict[str, Any]):
    """Settings of the identification of missing elements from a worker payload (None: defaults)."""
    if not worker_payload.get('identify_missing_elements', False):
        return None
    from autoemx.core.fitter.peak_identification import PeakIDSettings
    return PeakIDSettings(
        z_range=tuple(worker_payload.get('identification_z_range', dflt.identification_z_range)),
        excluded_elements=tuple(worker_payload.get('identification_excluded_elements', dflt.identification_excluded_elements)),
    )


def _spectrum_file_read_error(pointer_abs: Path) -> Optional[str]:
    """Return why a spectrum file cannot be read, or None if it can."""
    if not pointer_abs.exists():
        return "file not found"
    try:
        SampleLedger._load_counts_from_pointer_file(pointer_abs)
    except Exception as e:
        return f"{type(e).__name__}: {e}"
    return None


def _quantify_spectrum_worker(worker_payload: Dict[str, Any]) -> tuple[int, Optional[Dict[str, Any]], QuantificationResult, Optional[float]]:
    """Process-safe worker for one spectrum quantification or fitting-only."""
    spectrum_index = int(worker_payload['spectrum_index'])
    pointer_abs = Path(worker_payload['pointer_abs'])
    start_quant_time = time.time()

    # Empty / unreadable pointer files must not abort the parallel batch: map them to
    # quant_flag 1 (no spectral data) and continue with the remaining spectra.
    try:
        counts = SampleLedger._load_counts_from_pointer_file(pointer_abs)
        spectrum = np.asarray(counts, dtype=float)
    except Exception:
        spectrum = None

    background = None
    background_abs = worker_payload.get('background_abs')
    if background_abs:
        try:
            background = LedgerModule._load_background_vector_from_file(Path(background_abs))
        except Exception:
            background = None

    interrupt_fits_bad_spectra = bool(worker_payload['interrupt_fits_bad_spectra'])
    prefit_quant_flag, prefit_comment = detect_prefit_spectrum_issues(
        spectrum=spectrum,
        energy_vals=worker_payload['energy_vals'],
        sp_start=int(worker_payload['sp_start']),
        sp_end=int(worker_payload['sp_end']),
        target_acquisition_counts=int(worker_payload['target_acquisition_counts']),
        min_bckgrnd_cnts=dflt.prefit_min_bckgrnd_cnts,
        min_total_counts_fraction=float(
            worker_payload.get(
                'min_total_counts_fraction',
                QuantificationOptionsConfig.DEFAULT_MIN_TOTAL_COUNTS_FRACTION,
            )
        ),
    )
    if should_interrupt_prefit_issue(prefit_quant_flag, interrupt_fits_bad_spectra):
        quant_record = QuantificationResult(
            quantification_id=int(worker_payload['quantification_id']),
            quant_flag=prefit_quant_flag,
            comment=prefit_comment,
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
        els_w_fr=worker_payload.get('sample_w_frs'),
        is_particle=bool(worker_payload['apply_geom_factors']),
        sp_collection_time=float(worker_payload['sp_collection_time']),
        max_undetectable_w_fr=float(worker_payload['max_undetectable_w_fr']),
        fit_tol=float(worker_payload['fit_tolerance']),
        standards_dict=worker_payload['standards_dict'],
        verbose=False,
        fitting_verbose=False,
        identify_missing_elements=bool(worker_payload.get('identify_missing_elements', False)),
        identification_pool=worker_payload.get('identification_pool'),
        identification_settings=_identification_settings(worker_payload),
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
        quant_flag, comment = check_fit_quant_validity(
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

    quant_flag, comment = check_fit_quant_validity(
        is_quant_fit_valid=is_quant_fit_valid,
        bad_quant_flag=bad_quant_flag,
        quantifier=quantifier,
        min_bckgrnd_ref_lines=min_bckgrnd_ref_lines,
        detectable_els_substrate=worker_payload['detectable_els_substrate'],
        min_bckgrnd_cnts=worker_payload['min_bckgrnd_cnts'],
    )
    quant_flag, comment = merge_prefit_quant_flag(
        quant_flag,
        comment,
        prefit_quant_flag,
        prefit_comment,
    )
    quant_record = quantifier.export_quantification_result(
        quantification_id=int(worker_payload['quantification_id']),
        quant_result=quant_result,
        quant_flag=quant_flag,
        comment=comment,
    )
    return spectrum_index, quant_result, quant_record, time.time() - start_quant_time

# Columns of Compositions.csv written when the identification of missing elements is enabled
ADDED_ELEMENTS_DF_KEY = 'Added elements'
POSSIBLE_ELEMENTS_DF_KEY = 'Possible elements'


#%% EMXSp_Composition_Analyzer class
class EMXSp_Composition_Analyzer:
    """
    Main class for electron microscopy X-ray spectroscopy (EMXSp) composition analysis.

    This class orchestrates the acquisition, quantification, clustering, and plotting
    of X-ray spectra and composition data, using structured configuration objects for
    all instrument and analysis settings. To quantify or analyse a sample already acquired,
    create the analyzer from its ledger with ``EMXSp_Composition_Analyzer.from_ledger``.

    The sample ledger (ledger.json) is read and written by ``LedgerModule``
    (composition_analysis/ledger.py); legacy samples (Data.csv) are imported by
    ``autoemx.utils.legacy``.

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
        xsp_spot_selector: Optional[XSpSpotSelectorCallback] = None,
        particle_ids: Optional[List[int]] = None,
        particle_coords_mm: Optional[List[Tuple[float, float]]] = None,
    ):
        """
        Initialize the EMXSp_Composition_Analyzer with all configuration objects.

        See class docstring for parameter documentation.

        particle_ids / particle_coords_mm
            Used when ``powder_meas_cfg.par_selection_mode='list'``. Provide exactly
            one: existing ledger particle ids, or absolute stage coordinates in mm.
        """
        self.sample_id = self._resolve_sample_id(
            sample_id=sample_id,
            results_dir=results_dir,
            is_acquisition=is_acquisition,
        )

        # --- Record process time
        self.start_process_time = time.time()
        if verbose:
            print_double_separator(f"▶️ Starting compositional analysis of sample {self.sample_id}")
            
            
        # --- Define use of class instance
        self.is_acquisition = is_acquisition
        self.development_mode = development_mode
        self.xsp_spot_selector = xsp_spot_selector
        self._particle_ids = list(particle_ids) if particle_ids is not None else None
        self._particle_coords_mm = (
            list(particle_coords_mm) if particle_coords_mm is not None else None
        )

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
                raise ValueError(f"Missing detector calibration values.\n Please add detector calibration file at {calibs.calibration_files_dir}")
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
            LedgerModule._save_experimental_config(self, is_XSp_measurement)
        
        
        # --- Initialisations
        # Reconcile ledger / existing spectra before touching the microscope so resume
        # runs with sufficient data skip instrument initialization entirely.
        if is_acquisition:
            if self.measurement_cfg.is_manual_navigation and self.powder_meas_cfg.par_selection_mode == "list":
                raise ValueError(
                    "par_selection_mode='list' cannot be used with is_manual_navigation=True. "
                    "List mode moves the stage automatically to the listed particles (by ledger id or stage "
                    "coordinates), while manual navigation never moves the stage and leaves all positioning "
                    "to the user, so the list would be ignored. Set is_manual_navigation=False to measure the "
                    "listed particles, or use par_selection_mode='auto' with manual navigation."
                )

            spectrum_acquisition_needed = True
            is_particle_stats = (
                measurement_cfg.type == measurement_cfg.PARTICLE_STATS_MEAS_TYPE_KEY
            )
            if is_XSp_measurement:
                LedgerModule._initialise_acquisition_ledger(self)
                spectrum_acquisition_needed = self._spectrum_acquisition_needed

            if microscope_cfg.type == 'SEM' and spectrum_acquisition_needed:
                self._initialise_SEM()
            
            if is_XSp_measurement and spectrum_acquisition_needed:
                self._initialise_Xsp_analyzer()

            if is_particle_stats:
                LedgerModule._attach_particle_stats_ledger(self)

            if (
                is_XSp_measurement
                and spectrum_acquisition_needed
                and self.powder_meas_cfg.par_selection_mode == "list"
            ):
                LedgerModule._configure_list_particle_selection(self)
        

    @classmethod
    def from_ledger(
        cls,
        ledger: Union[str, Path, SampleLedger],
        configs: Optional[LedgerRuntimeConfigs] = None,
        *,
        exp_stds_cfg: Optional[ExpStandardsConfig] = None,
        sample_id: Optional[str] = None,
        **kwargs: Any,
    ) -> "EMXSp_Composition_Analyzer":
        """
        Create an analyzer of a sample already acquired, from its ledger (to quantify or analyse it again).

        Parameters
        ----------
        ledger : str, Path or SampleLedger
            Sample folder, path of its ledger.json, or its loaded ledger.
        configs : LedgerRuntimeConfigs, optional
            Configurations of the sample, e.g. ledger_runtime_configs(ledger) changed by the caller. Default: those of
            the ledger, with the quantification options and clustering configuration of its active quantification
            run (see ledger_runtime_configs).
        exp_stds_cfg : ExpStandardsConfig, optional
            Settings of the experimental standard, to fit the sample as a standard (as when it was acquired). Default:
            the sample is quantified and analysed as an ordinary sample, even if it is a standard.
        sample_id : str, optional
            Sample ID. Default: that of the ledger.
        **kwargs
            Other arguments of the analyzer, e.g. standards_dict, output_filename_suffix, verbose.

        Returns
        -------
        EMXSp_Composition_Analyzer
            Analyzer of the sample (not in acquisition mode).
        """
        if not isinstance(ledger, SampleLedger):
            ledger_path = Path(ledger)
            if ledger_path.is_dir():
                ledger_path = ledger_path / f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
            ledger = load_sample_ledger(ledger_path)
        if configs is None:
            configs = ledger_runtime_configs(ledger)
        return cls(
            microscope_cfg=configs.microscope_cfg,
            sample_cfg=configs.sample_cfg,
            measurement_cfg=configs.measurement_cfg,
            sample_substrate_cfg=configs.sample_substrate_cfg,
            quant_cfg=configs.quant_cfg,
            initial_clustering_cfg=configs.clustering_cfg,
            powder_meas_cfg=configs.powder_meas_cfg,
            bulk_meas_cfg=configs.bulk_meas_cfg,
            exp_stds_cfg=exp_stds_cfg if exp_stds_cfg is not None else ExpStandardsConfig(),
            plot_cfg=configs.plot_cfg,
            is_acquisition=False,
            sample_id=sample_id if sample_id is not None else ledger.sample_id,
            results_dir=ledger.sample_path,
            **kwargs,
        )

 


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
            LedgerModule._build_ledger_configs(self),
            sample_id=self.sample_id,
            results_dir=EM_images_dir,
            verbose=self.verbose,
            xsp_spot_selector=self.xsp_spot_selector,
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
        quantification_id, clustering_id = LedgerModule._resolve_active_analysis_config_ids(self)
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
                _ledger = LedgerModule._load_or_create_ledger(self)
                qc = LedgerModule._get_active_quantification_config(_ledger)
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
        min_total_counts_fraction_used = float(
            active_options.get(
                "min_total_counts_fraction",
                getattr(
                    q,
                    "min_total_counts_fraction",
                    QuantificationOptionsConfig.DEFAULT_MIN_TOTAL_COUNTS_FRACTION,
                ),
            )
        )
        lines.append(
            "  Min total counts     : "
            f"{format_min_total_counts_percent(min_total_counts_fraction_used)} of target"
        )
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
        if getattr(cc, "method", "kmeans") == "kmeans":
            lines.append(f"  Merge clusters       : {'yes' if getattr(cc, 'auto_merge_clusters', False) else 'no'}")
        lines.append(f"  Matrix decomposition : {'yes' if cc.do_matrix_decomposition else 'no'}")
        if cc.do_matrix_decomposition:
            lines.append(f"  Max mixture phases   : {cc.mixture.max_n_phases} (max recon. error {cc.mixture.max_recon_error})")
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
            logger.warning(f"Could not write analysis config summary: {e}")


    def _initialise_std_dict(self) -> None:
        """
        Initialise the dictionary of X-ray standards for quantification.
    
        This method determines how the `standards_dict` attribute is initialised
        based on the sample configuration and measurement type:
    
        - If the sample is a known phase mixture (quant_cfg.is_known_phase_mixture), the standards dictionary
          is compiled from reference data using `_compile_standards_from_references()`.
    
        - If standards were passed to the analyzer, or project-specific standards are used,
          the standards dictionary is taken from them.

        - Otherwise, the standards dictionary is expected to be loaded directly
          within the `XSp_Quantifier` and is set to `None` here.

        Returns
        -------
        None
            This method modifies the `self.standards_dict` attribute in place.
        """
        if getattr(self.quant_cfg, "is_known_phase_mixture", False):
            self.standards_dict = StandardsModule._compile_standards_from_references(self)
        elif self.quant_cfg.use_project_specific_std_dict or self.standards is not None:
            # _load_standards returns the standards passed to the analyzer when present
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
            w_fr_dict = composition_as_weight_dict(comp)
    
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
            LedgerModule._write_manufacturer_background_vector(self, 
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
    
    
    def _min_total_counts_fraction(self) -> float:
        """Return the configured minimum total-counts fraction for pre-fit flag 2."""
        return float(
            getattr(
                self.quant_cfg,
                "min_total_counts_fraction",
                QuantificationOptionsConfig.DEFAULT_MIN_TOTAL_COUNTS_FRACTION,
            )
        )

    def _include_quantified_added_elements(self) -> None:
        """
        Add the elements found and quantified by the identification of missing elements to the elements of the
        analysis (clustering, reference phases, plots). The original element list is restored by run_quantification, so that it does
        not change the quantification run.
        """
        if not hasattr(self, '_els_sample_before_added'):
            self._els_sample_before_added = list(self.all_els_sample)
        added_elements = [
            el for el, found in self.element_identification_summary().items()
            if found['quantified'] and el not in self._els_sample_before_added
        ]
        self.all_els_sample = list(self._els_sample_before_added) + added_elements
        self.detectable_els_sample = [el for el in self.all_els_sample if el not in calibs.undetectable_els]
        self.detectable_els_substrate = [el for el in self.detectable_els_substrate if el not in self.detectable_els_sample]
        if added_elements:
            logger.info(f"🔎 Added elements quantified and included in the analysis: {', '.join(added_elements)}")

    def _restore_elements_without_added(self) -> None:
        """Undo _include_quantified_added_elements (the quantification run is defined by the original elements)."""
        original = getattr(self, '_els_sample_before_added', None)
        if original is not None:
            self.all_els_sample = list(original)
            self.detectable_els_sample = [el for el in self.all_els_sample if el not in calibs.undetectable_els]
            del self._els_sample_before_added

    def element_identification_summary(self) -> Dict[str, Dict[str, List[int]]]:
        """
        Elements found by the identification of missing elements in the spectra of the current quantification run.

        Returns
        -------
        dict
            {element: {'quantified': [spectrum indices], 'not quantified': [...], 'possible': [...]}}, where
            'not quantified' elements were added to the fit but have no standard, and 'possible' ones were not
            confirmed by the fit (to be checked manually).
        """
        summary: Dict[str, Dict[str, List[int]]] = {}
        for i, record in enumerate(getattr(self, 'spectra_quant_records', None) or []):
            diag = getattr(record, 'diagnostics', None) if record is not None else None
            if diag is None:
                continue
            for key, els in (('quantified', diag.added_elements), ('not quantified', diag.added_unquantified_elements),
                             ('possible', diag.possible_elements)):
                for el in els or []:
                    summary.setdefault(el, {'quantified': [], 'not quantified': [], 'possible': []})[key].append(i)
        return summary

    @staticmethod
    def format_element_identification_summary(summary: Dict[str, Dict[str, List[int]]], sample_ID: str = '') -> str:
        """Text report of element_identification_summary()."""
        title = f"Element identification{f' — {sample_ID}' if sample_ID else ''}"
        if not summary:
            return f"{title}: no missing elements found."
        lines = [f"{title}:"]
        labels = {'quantified': 'added and quantified', 'not quantified': 'added, not quantified (no standard)',
                  'possible': 'possibly present, check manually'}
        for el in sorted(summary, key=lambda e: -sum(len(v) for v in summary[e].values())):
            for key, label in labels.items():
                spectra = summary[el][key]
                if spectra:
                    lines.append(f"   • {el}: {label} in {len(spectra)} spectra (#{', #'.join(map(str, spectra))})")
        return "\n".join(lines)

    def report_added_elements(self) -> Dict[str, Dict[str, List[int]]]:
        """Log the elements found in the sample and write them to Element_identification_report.txt in the analysis folder."""
        summary = self.element_identification_summary()
        text = self.format_element_identification_summary(summary, getattr(self.sample_cfg, 'ID', ''))
        print_double_separator(text)
        try:
            if getattr(self, 'analysis_dir', None):
                with open(os.path.join(self.analysis_dir, 'Element_identification_report.txt'), 'w', encoding='utf-8') as file:
                    file.write(text + "\n")
        except OSError as e:
            logger.warning(f"⚠️ Could not write the element identification report: {e}")
        return summary

    def _identify_missing_elements(self) -> bool:
        """Whether the identification of missing elements (elements missing from the fit) is enabled."""
        return bool(getattr(self.quant_cfg, "identify_missing_elements", dflt.identify_missing_elements))

    def _identification_seed_spectra(self) -> int:
        """Number of spectra checked before the elements found in them are tested first in the others."""
        return int(getattr(self.quant_cfg, "identification_seed_spectra", dflt.identification_seed_spectra))

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
        ingested_spectrum_ids = LedgerModule._get_ingested_spectrum_ids(self)

        # Auto-detect next available spectrum ID if not provided
        if n_tot_sp_collected is None:
            pointer_files = _list_spectrum_pointer_files(Path(self.sample_result_dir, cnst.SPECTRA_DIR))
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
            
            particle_finder = None
            if getattr(self, "EM_controller", None) is not None:
                particle_finder = getattr(self.EM_controller, "particle_finder", None)

            if particle_cntr is not None:
                if particle_finder is not None and getattr(
                    particle_finder, "uses_absolute_particle_ids", False
                ):
                    # List-mode particle ids are already absolute ledger ids.
                    self.particle_cntr = int(particle_cntr)
                else:
                    self.particle_cntr = int(particle_cntr) + particle_id_offset
            else:
                self.particle_cntr = None
            frame_ID = self.EM_controller.current_frame_label

            particle_area_um: Optional[float] = None
            particle_coordinates: Optional[Coordinate2D] = None
            if (
                self.sample_cfg.is_particle_acquisition
                and self.particle_cntr is not None
                and particle_finder is not None
            ):
                particle_area_um = getattr(particle_finder, "current_par_area_um2", None)
                par_coords = getattr(particle_finder, "current_par_coordinates_mm", None)
                if par_coords is not None:
                    particle_coordinates = Coordinate2D(
                        x=float(par_coords[0]),
                        y=float(par_coords[1]),
                    )
            
            latest_spot_id = None # For image annotations
            for i, (x, y) in enumerate(spots_xy_list):
                latest_spot_id = i
                xy_center = (int(x), int(y))
                current_spectrum_id = str(next_spectrum_id)

                if current_spectrum_id in ingested_spectrum_ids:
                    next_spectrum_id += 1
                    continue

                if self.verbose:
                    print_single_separator(f'🔬 Acquiring spectrum #{next_spectrum_id}...')

                spectrum_relpath = LedgerModule._build_spectrum_relpath(self, current_spectrum_id)
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
                # Build SpectrumEntry directly from the pointer file and current acquisition metadata
                frame_id_str = (
                    str(frame_ID).strip()
                    if frame_ID is not None and str(frame_ID).strip()
                    else None
                )
                acq_details = AcquisitionDetails(
                    frame_id=frame_id_str,
                    particle_id=self.particle_cntr if self.particle_cntr is not None else None,
                    spot_coordinates=LedgerModule._build_spot_coordinates(self, 
                        pixel_x=xy_center[0],
                        pixel_y=xy_center[1],
                    ),
                )
                spectrum_entry = _build_spectrum_entry(
                    Path(self.sample_result_dir), Path(manufacturer_msa_path), acquisition_details=acq_details
                )
                is_particle = self.sample_cfg.is_particle_acquisition and self.particle_cntr is not None
                LedgerModule._record_acquired_spectrum(
                    self,
                    spectrum_entry,
                    particle_id=int(self.particle_cntr) if is_particle else None,
                    particle_area_um=particle_area_um,
                    frame_id=frame_id_str,
                    particle_coordinates=particle_coordinates,
                    # Pixel size of the particle image, saved with the spectrum spots
                    image_pixel_size_um=getattr(self.EM_controller, "pixel_size_um", None),
                )
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
            
            # Save image of particle. Spot pixel coordinates are stored in the ledger, so raw images
            # can be annotated afterwards with Annotate_Particle_Images.py.
            if latest_spot_id is not None:
                # Prepare save path
                par_cntr_str = f"_par{self.particle_cntr}" if self.particle_cntr is not None else ''
                filename = f"{self.sample_id}{par_cntr_str}_fr{frame_ID}_xyspots"
                if not self.measurement_cfg.annotate_particle_images:
                    self.EM_controller.save_frame_image(filename, im_annotations=None, scalebar=False)
                else:
                    # Construct annotation dictionary with the ID of the acquired XSp spots
                    im_annotations = []
                    for i, xy_coords in enumerate(spots_xy_list):
                        # Skip if latest_spot_id is None or i is out of range
                        if latest_spot_id is None or i > latest_spot_id:
                            break

                        xy_center = (int(xy_coords[0]), int(xy_coords[1]))
                        if xy_center is None:
                            continue
                        
                        im_annotations.append(
                            xsp_spot_annotation(next_spectrum_id - 1 - latest_spot_id + i, xy_center)
                        )
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
            max_spectra_to_quantify: Optional[int] = None,
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
            max_spectra_to_quantify : Optional[int]
                Process at most this number of spectra, among those that need processing.
                None processes all of them.
            """
            
            _ledger = LedgerModule._load_or_create_ledger(self)
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
            # Spectra whose file cannot be read: (spectrum index, relative path, reason)
            unreadable_spectra: List[Tuple[int, str, str]] = []

            if quantify:
                # Always bootstrap/sync ledger before quantification cycles.
                LedgerModule._ensure_current_quantification_run(self, force_new=force_requantification)
                LedgerModule._persist_current_quantification_config(self)
                LedgerModule._ensure_quant_tracking_length(self, tot_spectra_collected)
                if force_requantification:
                    indices_to_process = list(range(tot_spectra_collected))
                elif requantify_only_unquantified_spectra:
                    LedgerModule._sync_existing_quantification_from_ledger(self)
                    indices_to_process = [
                        i for i in range(tot_spectra_collected)
                        if self.spectra_quant_records[i] is None
                        or self.spectra_quant_records[i].composition_atomic_fractions is None
                    ]
                else:
                    LedgerModule._sync_existing_quantification_from_ledger(self)
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
                if max_spectra_to_quantify is not None:
                    indices_to_process = indices_to_process[:max(0, int(max_spectra_to_quantify))]
            else:
                # Fitting-only (experimental standards) path
                LedgerModule._ensure_current_quantification_run(self, force_new=force_requantification)
                LedgerModule._persist_current_quantification_config(self)
                LedgerModule._ensure_quant_tracking_length(self, tot_spectra_collected)
                LedgerModule._sync_existing_quantification_from_ledger(self)
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
                    read_error = _spectrum_file_read_error(Path(self.sample_result_dir, spectrum_entry.spectrum_relpath))
                    if read_error is not None:
                        unreadable_spectra.append((i, spectrum_entry.spectrum_relpath, read_error))
                        continue
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
                        'min_total_counts_fraction': self._min_total_counts_fraction(),
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
                            if self.sample_cfg.w_frs is not None
                            else None
                        ),
                        'apply_geom_factors': bool(self._apply_geom_factors),
                        'max_undetectable_w_fr': float(self.undetectable_an_er),
                        'fit_tolerance': float(self.quant_cfg.fit_tolerance),
                        'standards_dict': self.standards_dict,
                        'interrupt_fits_bad_spectra': bool(interrupt_fits_bad_spectra),
                        'quantify': False,
                        'identify_missing_elements': self._identify_missing_elements(),
                        'identification_pool': [],
                        'identification_z_range': list(getattr(self.quant_cfg, 'identification_z_range', dflt.identification_z_range)),
                        'identification_excluded_elements': list(getattr(self.quant_cfg, 'identification_excluded_elements', dflt.identification_excluded_elements)),
                    })

            n_spectra_to_quant = len(indices_to_process)
        
            if self.verbose and n_spectra_to_quant > 0:
                print_single_separator()
                quant_str = "quantification" if quantify else "fitting"
                logger.info(f"▶️ Starting {quant_str} of {n_spectra_to_quant} spectra on up to {_n_cores} cores.\nThis may take a while...")
        
            if quantify:
                indices_to_process_set = set(indices_to_process)
                n_already_quantified = sum(
                    1 for i in range(tot_spectra_collected)
                    if i not in indices_to_process_set and self.spectra_quant_records[i] is not None
                )
                if self.verbose and n_already_quantified > 0 and not force_requantification:
                    print_single_separator(f"⏭️ {n_already_quantified} spectra already quantified. Skipping them...")
                for i in range(tot_spectra_collected):
                    if i not in indices_to_process_set:
                        continue

                    spectrum_entry = _ledger_spectra[i]
                    read_error = _spectrum_file_read_error(Path(self.sample_result_dir, spectrum_entry.spectrum_relpath))
                    if read_error is not None:
                        unreadable_spectra.append((i, spectrum_entry.spectrum_relpath, read_error))
                        continue
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
                        'min_total_counts_fraction': self._min_total_counts_fraction(),
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
                        'identify_missing_elements': self._identify_missing_elements(),
                        'identification_pool': [],
                        'identification_z_range': list(getattr(self.quant_cfg, 'identification_z_range', dflt.identification_z_range)),
                        'identification_excluded_elements': list(getattr(self.quant_cfg, 'identification_excluded_elements', dflt.identification_excluded_elements)),
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
                            fitted_peaks = getattr(fit_result, 'fitted_peaks', None) or {}
                            pb_lines = XSp_Quantifier.format_measured_pb(
                                (
                                    (getattr(peak, 'element', None), getattr(peak, 'line', None), getattr(peak, 'pb_ratio', None))
                                    for peak in fitted_peaks.values()
                                ),
                                standard_elements=self.sample_cfg.w_frs,
                            )
                            if pb_lines:
                                lines.append("  Measured PB:")
                                lines.extend(f"  {pb_line}" for pb_line in pb_lines)
                            else:
                                lines.append("  Fit completed")
                        if quant_flag not in (None, 0) and comment:
                            lines.append(f"  {comment}")
                        if quantification_time is not None:
                            lines.append(f"  Quantification took {quantification_time:.2f} s")
                    lines.append("")  # Add an extra line for spacing
                    print_single_separator("\n".join(lines))

                if quant_record is not None:
                    overwrite = (
                        force_requantification
                        or (requantify_only_unquantified_spectra and had_prior_record)
                        or was_interrupted
                    )
                    LedgerModule._persist_quantification_record(self, idx, quant_record, overwrite=overwrite)
            
            def _execute(payloads):
                try:
                    if not use_parallel:
                        for payload in payloads:
                            idx, result, quant_record, quantification_time = _quantify_spectrum_worker(payload)
                            _finalize_quant_result(idx, result, quant_record, quantification_time)
                    else:
                        # Compute/load convolution matrices once on the main process so workers
                        # never race on a cold detector-cache miss.
                        if payloads:
                            DetectorResponseFunction.ensure_conv_matrices_cached(
                                float(self.det_ch_offset),
                                float(self.det_ch_width),
                                self.microscope_cfg.ID,
                                verbose=bool(self.verbose),
                            )

                        # Stream completed tasks back to the main process so progress is visible in real time.
                        try:
                            completed = Parallel(
                                n_jobs=_n_cores,
                                backend=parallel_backend,
                                return_as='generator_unordered',
                            )(
                                delayed(_quantify_spectrum_worker)(payload) for payload in payloads
                            )
                        except TypeError:
                            # Compatibility fallback for joblib versions without return_as.
                            completed = Parallel(n_jobs=_n_cores, backend=parallel_backend)(
                                delayed(_quantify_spectrum_worker)(payload) for payload in payloads
                            )

                        for idx, result, quant_record, quantification_time in completed:
                            _finalize_quant_result(idx, result, quant_record, quantification_time)
                except Exception as e:
                    logger.warning(
                        f"⚠️ Parallel spectrum processing failed ({type(e).__name__}: {e}), "
                        "falling back to sequential execution."
                    )
                    for payload in payloads:
                        idx, result, quant_record, quantification_time = _quantify_spectrum_worker(payload)
                        _finalize_quant_result(idx, result, quant_record, quantification_time)

            # Unreadable spectrum files. A spectrum already quantified in this run keeps its result, since
            # deleting it after quantification is assumed to be a mistake. Otherwise, it is flagged as having no data.
            for idx, relpath, read_error in unreadable_spectra:
                has_result = any(
                    result.quantification_id == self.current_quantification_id
                    for result in _ledger_spectra[idx].quantification_results
                )
                if has_result:
                    logger.warning(f"⚠️ Could not read spectrum file {relpath} ({read_error}). Keeping its existing result.")
                    continue
                logger.warning(f"⚠️ Could not read spectrum file {relpath} ({read_error}).")
                quant_record = QuantificationResult(
                    quantification_id=int(self.current_quantification_id),
                    quant_flag=1,
                    comment=f"Could not read spectrum file ({read_error})",
                    diagnostics=QuantificationDiagnostics(converged=False, interrupted=True),
                )
                _finalize_quant_result(idx, None, quant_record, None)

            def _build_identification_pool(records):
                pool = []
                for record in records:
                    diag = getattr(record, 'diagnostics', None) if record is not None else None
                    for el in (getattr(diag, 'added_elements', None) or []) + (getattr(diag, 'added_unquantified_elements', None) or []):
                        if el not in pool:
                            pool.append(el)
                return pool

            # Element pool: elements confirmed in the first spectra are tested first in the others.
            # Spectra already processed in this run (e.g. earlier particles during acquisition) count as seed.
            n_seed = self._identification_seed_spectra()
            payload_indices = {int(payload['spectrum_index']) for payload in quant_worker_payloads}
            prior_records = [
                record for i, record in enumerate(self.spectra_quant_records[:tot_spectra_collected])
                if i not in payload_indices and record is not None
            ]
            n_seed_missing = n_seed - len(prior_records)
            if self._identify_missing_elements() and n_seed > 0 and n_seed_missing <= len(quant_worker_payloads):
                n_seed_missing = max(0, n_seed_missing)
                seed_payloads, other_payloads = quant_worker_payloads[:n_seed_missing], quant_worker_payloads[n_seed_missing:]
                _execute(seed_payloads)
                seed_records = [self.spectra_quant_records[int(payload['spectrum_index'])] for payload in seed_payloads]
                pool = _build_identification_pool(prior_records + seed_records)
                if pool and seed_payloads:
                    logger.info(f"🔎 Element identification: elements found in the first {n_seed} spectra, "
                                f"tested first in the others: {', '.join(pool)}")
                for payload in other_payloads:
                    payload['identification_pool'] = list(pool)
                _execute(other_payloads)
            else:
                _execute(quant_worker_payloads)


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
        # The clustering config must record the error threshold actually used, so that a new value
        # creates a new clustering analysis in the ledger instead of matching the existing one
        self.clustering_cfg.max_analytical_error_percent = max_analytical_error_percent

        # Ensure in-memory state is fully synced from the ledger before analysis
        LedgerModule._sync_active_quantification_from_ledger(self)

        # A changed min_bckgrnd_cnts is applied from the background counts stored with each result
        n_without_bckgrnd_cnts = count_results_without_bckgrnd_cnts(self)
        if n_without_bckgrnd_cnts > 0:
            logger.warning(
                f"⚠️ min_bckgrnd_cnts = {self.clustering_cfg.min_bckgrnd_cnts} cannot be applied to "
                f"{n_without_bckgrnd_cnts} spectra quantified by an older AutoEMX version, which did not store their "
                "background counts: they keep their quant_flag. You should requantify first to apply this parameter."
            )

        # Elements found and quantified by the identification of missing elements become additional elements of the
        # analysis (zero in the spectra without them), so that particles containing them can form their own clusters
        self._include_quantified_added_elements()

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
            print_single_separator('ℹ️ Spectra selection:')
            logger.info(f"ℹ️ {n_datapts_used} data points are used, out of {n_datapts} collected spectra.")
            self._report_n_discarded_spectra(n_datapts, max_analytical_error)

        # 3. Prepare DataFrames for clustering
        compositions_df, compositions_df_other_fr = ClusteringModule._prepare_composition_dataframes(self, compositions_list_at, compositions_list_w)
        # Features clustering operates on (CLR coordinates in Aitchison geometry). Centroids and
        # cluster statistics are always reported in fraction space (compositions_df).
        # With geometry='auto', the geometry is selected here from the measured compositions.
        self.resolved_clustering_geometry = ClusteringModule._resolve_geometry(self, compositions_df)
        clustering_df = ClusteringModule._get_clustering_features(
            self, compositions_df, self.resolved_clustering_geometry
        )
        is_aitchison = self.resolved_clustering_geometry == 'aitchison'

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
            k = ClusteringModule._find_optimal_k(
                self, compositions_df, k, compute_k_only_once, clustering_df=clustering_df,
                max_analytical_error=max_analytical_error,
            )
            LedgerModule._persist_resolved_k_on_active_clustering_config(self, k)
            kmeans, labels, sil_score = ClusteringModule._run_kmeans_clustering(self, k, clustering_df)
            if (
                k > 1
                and getattr(self.clustering_cfg, 'auto_merge_clusters', False)
                and self.clustering_cfg.k_finding_method != 'forced'
            ):
                labels, k_merged = ClusteringModule._merge_connected_clusters(self, clustering_df, labels)
                if k_merged < k:
                    if self.verbose:
                        print_single_separator(f"ℹ️ Clusters merged: {k} -> {k_merged} (auto_merge_clusters).")
                    k = k_merged
                    LedgerModule._persist_resolved_k_on_active_clustering_config(self, k)
                    # The k-means model no longer describes the clusters
                    kmeans = None
                    sil_score = silhouette_score(clustering_df, labels) if k > 1 else np.nan
            if kmeans is None:
                centroids, wcss = ClusteringModule._centroids_and_wcss_from_labels(
                    compositions_df.to_numpy(), labels, k
                )
            elif is_aitchison:
                # k-means centers live in CLR space; report arithmetic means in fraction space.
                centroids, wcss = ClusteringModule._centroids_and_wcss_from_labels(
                    compositions_df.to_numpy(), labels, k
                )
            else:
                centroids = kmeans.cluster_centers_
                wcss = kmeans.inertia_
        elif self.clustering_cfg.method == 'dbscan':
            # DBSCAN determines the number of clusters from data density; k_forced/k_finding
            # settings are ignored. No KMeans model exists, so plotting must handle kmeans=None.
            kmeans = None
            labels, centroids, k, sil_score, wcss = ClusteringModule._run_dbscan_clustering(
                self, compositions_df, clustering_df=clustering_df, geometry=self.resolved_clustering_geometry
            )
            LedgerModule._persist_resolved_k_on_active_clustering_config(self, k)
            if k < 1:
                print_single_separator()
                logger.warning(
                    "⚠️ DBSCAN found no clusters (all points classified as noise) with "
                    f"eps={self.clustering_cfg.dbscan.resolved_eps(self.resolved_clustering_geometry):g} "
                    f"({self.resolved_clustering_geometry} geometry). Try increasing 'eps' or lowering 'min_samples'."
                )
                self._save_analysis_summary(None, None)
                return False, 0, 0  # zeroes are placeholders
        else:
            raise ValueError(f"Unsupported clustering method: {self.clustering_cfg.method}")
    
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
        LedgerModule._persist_clustering_result_on_active_analysis(self, clustering_artifacts)

        # 8b. Back-fill particle cluster/composition summaries when particles exist
        LedgerModule._backfill_particles_from_clustering(self, labels, df_indices, clustering_artifacts)
    
        # 9. Save plots
        if self.plot_cfg.save_plots:
            PlottingModule._save_plots(
                self, kmeans, compositions_df, centroids, labels, els_std_dev_per_cluster, unused_compositions_list,
                silhouette_df=clustering_df,
            )
            if getattr(self.plot_cfg, 'plot_best_mixture', True) and clusters_assigned_mixtures:
                PlottingModule._save_best_mixture_plots(self, compositions_df, labels, clusters_assigned_mixtures)

        # Show all open figures (elbow, silhouette, clustering, mixture) together, blocking until they are closed
        if self.plot_cfg.show_plots and plt.get_fignums():
            plt.show(block=True)
    
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
        resume_state = getattr(self, "_acquisition_resume_state", None)
        if resume_state is None:
            resume_state = LedgerModule._prepare_acquisition_resume_state(self)
        else:
            resume_state = LedgerModule._prepare_acquisition_resume_state(self, log_resume_info=False)
            self._acquisition_resume_state = resume_state
        tot_n_spectra = resume_state["tot_n_spectra"]
        next_spectrum_id = resume_state["next_spectrum_id"]
        _particle_id_offset = resume_state["particle_id_offset"]

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
            self._initialise_std_dict() # Initialise dictionary of standards to (optionally) pass onto XSp_Quantifier. Only used with known phase mixtures
        
        if quantify:
            if is_exp_std_measurement:
                quant_str, quant_verb_str = ' and fitting', ' and fitting'
            else:
                quant_str, quant_verb_str = ' and quantification', ' and quantifying'
        else:
            quant_str, quant_verb_str = '', ''

        
        if self.verbose:
            print_double_separator(f"▶️ Starting acquisition{quant_str} of {tot_spectra_to_collect} spectra.")
        if is_spectral_quant:
            self._warn_peak_overlaps()
        
        while tot_n_spectra < self.max_n_spectra:
            if self.verbose:
                print_double_separator(f"🔬 Collecting{quant_verb_str} {n_spectra_to_collect} spectra...")
    
            # Collect the next batch of spectra (and quantify if requested)
            n_spectra_before = tot_n_spectra
            next_spectrum_id, is_acquisition_successful = self._collect_spectra(
                n_spectra_to_collect,
                n_tot_sp_collected=next_spectrum_id,
                quantify=is_spectral_quant,
                particle_id_offset=_particle_id_offset,
                interrupt_fits_bad_spectra=interrupt_fits_bad_spectra
            )
            refreshed_ledger = LedgerModule._load_or_create_ledger(self)
            tot_n_spectra = len(refreshed_ledger.spectra) if refreshed_ledger is not None else 0
            n_spectra_collected_this_session += (tot_n_spectra - n_spectra_before)
    
            if self.verbose:
                print_single_separator(f"✅ {tot_n_spectra} total spectra collected.")
            
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
                final_ledger = LedgerModule._load_or_create_ledger(self)
                n_particles = len({
                    s.acquisition_details.particle_id
                    for s in (final_ledger.spectra if final_ledger is not None else [])
                    if s.acquisition_details is not None and s.acquisition_details.particle_id is not None
                })
                par_str = f' over {n_particles} particles'
            else:
                par_str = ' (no new spectra collected in this session)'
        else:
            par_str = ''
        if n_spectra_collected_this_session > 0:
            logger.info(f'✅ {n_spectra_collected_this_session} new spectra collected (total: {tot_n_spectra}{par_str}).')
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
    
                self.print_results(warn_peak_overlaps=False)
    
            elif not is_acquisition_successful:
                logger.warning('⚠️ This did not allow to determine which phases are present in the sample.')
            else:
                logger.warning(f'⚠️ Phases could not be identified with the allowed maximum of {self.max_n_spectra} collected spectra.')
                is_analysis_successful = False
                is_converged = False
        else:
            is_analysis_successful = is_acquisition_successful

        if is_spectral_quant:
            self._warn_peak_overlaps()

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
            print_double_separator(f"📊 Analysing phases after collection of {tot_n_spectra} spectra...")
    
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
                print_double_separator("✅ Clustering analysis performed")
    
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
        max_spectra_to_quantify: Optional[int] = None,
        warn_peak_overlaps_at_end: bool = True,
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
        max_spectra_to_quantify : Optional[int], optional
            Quantify at most this number of spectra, among those that need quantification
            (in acquisition order), e.g. for quick tests. None quantifies all of them.
        warn_peak_overlaps_at_end : bool, optional
            Peak overlaps that may compromise the quantification are warned about at the start of the
            quantification and, if True, also at its end. Set to False when the analysis follows, as
            print_results warns again.
        """
        self._restore_elements_without_added()
        self._initialise_std_dict()
        self._warn_peak_overlaps()
        self._fit_and_quantify_spectra(
            force_requantification=force_requantification,
            requantify_only_unquantified_spectra=requantify_only_unquantified_spectra,
            interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
            num_CPU_cores=num_CPU_cores,
            max_spectra_to_quantify=max_spectra_to_quantify,
        )
        # Quantification-only workflows also export Compositions.csv, so ensure
        # a deterministic analysis directory exists instead of falling back to
        # the sample root.
        self._make_analysis_dir()
        self._save_analysis_summary(None, None)
        if self._identify_missing_elements():
            self.report_added_elements()
        if warn_peak_overlaps_at_end:
            self._warn_peak_overlaps()
        
    
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
            `self.sample_cfg.w_frs` (the composition of the standard).
        """
    
        if not self.exp_stds_cfg.is_exp_std_measurement:
            raise ValueError(
                "Experimental standard collection mode is not active. "
                "Set `self.exp_stds_cfg.is_exp_std_measurement = True` before running."
            )
        
        # Ensure all elemental fractions of the standard are defined
        missing = [el for el in self.sample_cfg.elements if el not in (self.sample_cfg.w_frs or {})]
        if missing:
            raise KeyError(
                f"The following elements are missing from `exp_stds_cfg.formula`: "
                f"{', '.join(str(m) for m in missing)}. "
                f"Ensure the formula contains all elements defined in `self.sample_cfg.elements`."
            )
        
        if self.verbose:
            print_double_separator(f"🔬 Experimental standard acquisition of {self.sample_id}")
        
        # Run collection and quantification (fitting optionally performed during collection)
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
    
        # Homogeneity of each cluster, from the spread of its spectra in atomic fractions
        rms_dist_at_fr = rms_dist_cluster if fr_key == cnst.AT_FR_DF_KEY else rms_dist_cluster_other_fr
        homogeneity = []
        for n_points, rms_dist in zip(n_points_per_cluster, rms_dist_at_fr):
            if n_points < dflt.homogeneity_min_n_points:
                homogeneity.append("Not Applicable")
            elif rms_dist < dflt.homogeneous_cluster_max_rms_dist:
                homogeneity.append("Likely")
            elif rms_dist <= dflt.inhomogeneous_cluster_min_rms_dist:
                homogeneity.append("Unsure")
            else:
                homogeneity.append("Unlikely")

        # Compose main cluster DataFrame
        clusters_dict = {
            cnst.N_PTS_DF_KEY: n_points_per_cluster,
            **els_fr_dict,
            **els_fr_stdevs_dict,
            **els_other_fr_dict,
            **els_other_fr_stdevs_dict,
            cnst.RMS_DIST_DF_KEY + fr_key: (np.array(rms_dist_cluster) * 100).round(2),
            cnst.RMS_DIST_DF_KEY + other_fr_key: (np.array(rms_dist_cluster_other_fr) * 100).round(2),
            cnst.HOMOGENEOUS_DF_KEY: homogeneity,
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
    
        Raises
        ------
        OSError
            If the output directory cannot be created or files cannot be written.
        """
    
        ledger = LedgerModule._load_or_create_ledger(self)
        ledger_spectra = list(ledger.spectra) if ledger is not None else []
        n_spectra = len(ledger_spectra)
    
        if n_spectra == 0:
            return None

        particles_by_id: Optional[Dict[int, ParticleInfo]] = None
        if ledger.particles:
            particles_by_id = {particle.id: particle for particle in ledger.particles}
        
        data_list = []
        for i in range(n_spectra):
            # Retrieve the typed QuantificationResult for this spectrum (None when not quantified)
            record = self.spectra_quant_records[i] if i < len(self.spectra_quant_records) else None
            coords = LedgerModule._extract_spectrum_info(self, 
                ledger_spectra[i], i, particles_by_id=particles_by_id
            )
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
                # Flag for the min_bckgrnd_cnts of this analysis (quant_flag 8 re-evaluated if it changed)
                quant_flag, comment = effective_quant_flag(record, self)
                data_row[cnst.COMMENTS_DF_KEY] = comment
                data_row[cnst.QUANT_FLAG_DF_KEY] = quant_flag
                # Results of the identification of missing elements (only when the check was run)
                diag = getattr(record, 'diagnostics', None)
                if diag is not None and getattr(diag, 'added_elements', None) is not None:
                    found = list(diag.added_elements or []) + [
                        f"{el} (not quantified)" for el in (diag.added_unquantified_elements or [])
                    ]
                    data_row[ADDED_ELEMENTS_DF_KEY] = ', '.join(found)
                    data_row[POSSIBLE_ELEMENTS_DF_KEY] = ', '.join(diag.possible_elements or [])
    
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
            ADDED_ELEMENTS_DF_KEY, POSSIBLE_ELEMENTS_DF_KEY,
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
            2: Total number of counts is lower than min_total_counts_fraction of target counts (default 90%), likely due to wrong segmentation. Fit interrupted if interrupt_fits_bad_spectra=True
            3: Spectrum has too low signal in its low-energy portion, leading to poor quantification in this region. Fit interrupted if interrupt_fits_bad_spectra=True
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
    
        print_single_separator("📊 Summary of Discarded Spectra")
        logger.info("  → For details, see the 'Comments' column in Compositions.csv.")
    
        # Discarded due to low counts, insufficient background, or acquisition/fitting errors
        if self.n_sp_too_low_counts > 0:
            logger.info(
                f"  • {self.n_sp_too_low_counts} spectra were discarded due to insufficient total counts, "
                f"background counts below the threshold ({self.clustering_cfg.min_bckgrnd_cnts}), "
                "or errors during spectrum collection/fitting."
            )
    
        # Discarded due to quantification flags
        if self.n_sp_bad_quant > 0:
            logger.info(
                f"  • {self.n_sp_bad_quant} spectra were discarded because they were flagged during quantification."
            )
    
        # Warning if more than half of the spectra were flagged
        if self.n_sp_bad_quant / n_datapts > 0.5:
            print_single_separator()
            logger.warning("  ⚠️ Warning: More than 50% of spectra were flagged during quantification!")
            logger.info(
                "  Common causes for poor fits (quant_flag = 4) include missing elements in the fit.\n"
                "  Ensure that all elements present in your sample have been specified in the 'elements' argument "
                "  when initializing EMXSp_Composition_Analyzer."
            )
            
        # Discarded due to high analytical error
        if self.n_sp_too_high_an_err > 0:
            logger.info(
                f"  • {self.n_sp_too_high_an_err} spectra were discarded because their analytical error "
                f"exceeded the maximum allowed value of {max_analytical_error*100:.1f}%."
            )
    
        print_single_separator()
            
        
    def print_results(self, n_cnd_to_print = 2, n_mix_to_print = 2, warn_peak_overlaps = True) -> None:
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
        warn_peak_overlaps : bool
            If True, warns about peak overlaps that may compromise the quantification, after the results.
        
        Raises
        ------
        AttributeError
            If required attributes (clustering_info, clusters_df, etc.) are missing.
        KeyError
            If expected keys are missing from clustering_info or clusters_df.
        """
        # Print clustering info
        print_double_separator(f"📊 Compositional analysis results for sample {self.sample_id}:")
        print_single_separator()
        try:
            logger.info('  Clustering method: %s', self.clustering_cfg.method)
            geometry = self.clustering_cfg.geometry
            resolved = getattr(self, 'resolved_clustering_geometry', None)
            if geometry == 'auto' and resolved is not None:
                geometry = f'auto -> {resolved}'
            logger.info('  Clustering geometry: %s', geometry)
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

        if warn_peak_overlaps:
            self._warn_peak_overlaps()


    def _warn_peak_overlaps(self) -> None:
        """Warn about reference peaks of sample elements overlapping with other peaks whose area is fitted freely."""
        energy_vals = getattr(self, 'energy_vals', None)
        energy_range = (energy_vals[0], energy_vals[-1]) if energy_vals is not None and len(energy_vals) else None
        warn_peak_overlaps(
            self.measurement_cfg.type,
            getattr(self, 'detectable_els_sample', []),
            getattr(self, 'detectable_els_substrate', []),
            beam_energy_keV=self.measurement_cfg.beam_energy_keV,
            energy_range_keV=energy_range,
            microscope_ID=self.microscope_cfg.ID,
        )


