#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch fitting of experimental standards spectra (PB ratios) from existing ledgers.

This script is analogous to batch_acquire_experimental_stds.py, but instead of running collection and quantification,
it directly fits and updates standards from existing data using run_exp_std_fit.

Usage:
    - Edit the 'std_list' list to define your standards (same as in batch_acquire_experimental_stds.py)
    - Adjust configuration parameters as needed
    - Run the script, or import and call batch_fit_experimental_stds() to fit/update standards for one or multiple samples

Created on May 20, 2026
@author: Andrea
"""

import logging
import os
import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, cast

from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
from autoemx.utils.helper import print_double_separator, get_sample_dir
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import ClusteringConfig  # type: ignore
import autoemx.utils.constants as cnst
import autoemx.config.defaults as dflts
from autoemx.config.schema_models.standards import EDSStandardsFile
import autoemx.config as config_module

config_classes_dict = cast(Dict[str, Any], getattr(config_module, "config_classes_dict"))


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["batch_fit_experimental_stds"]


def batch_fit_experimental_stds(
    std_ids: List[str],
    spectrum_lims: Tuple[float, float] = dflts.spectrum_lims,
    use_instrument_background: bool = dflts.use_instrument_background,
    min_bckgrnd_cnts: float = 5,
    update_std_library: bool = True,
    force_refitting: bool = False,
    els_substrate: Optional[List[str]] = None,
    output_filename_suffix: str = '',
    verbose: bool = True,
    stds_path: Optional["str | Path"] = None,
    standards_dict_path: Optional["str | Path"] = None,
) -> List[Optional[EMXSp_Composition_Analyzer]]:
    """
    Batch fitting of experimental standard spectra from existing ledgers.

    Parameters
    ----------
    std_ids : list of str
        List of identifiers for the experimental standards.
    spectrum_lims : tuple of float, optional
        Lower and upper energy limits for spectrum fitting in eV.
    use_instrument_background : bool, optional
        Whether to use instrument background files during fitting.
    min_bckgrnd_cnts : float, optional
        Minimum background counts required for a spectrum not to be filtered out.
    update_std_library : bool, optional
        If True, the experimental standard library will be updated with the newly fitted PB ratios.
    force_refitting : bool, optional
        If True, force refitting even if results with same fitting configurations already exist for a standard.
        Default is False.
    els_substrate : Optional[List[str]], optional
        List of elements in the substrate. If None, defaults to all elements in the sample.
    output_filename_suffix : str, optional
        Suffix for the output filename.
    verbose : bool, optional
        If True, print verbose output.
    stds_path : Optional[str or pathlib.Path], optional
        Directory containing experimental standard data.
    standards_dict_path : Optional[str or pathlib.Path], optional
        Path to a standards dictionary file.

    Returns
    -------
    results : list(EMXSp_Composition_Analyzer)
        A list of the composition analysis objects (one per sample).
    """
    if stds_path is None:
        stds_path = os.path.join(os.getcwd(), cnst.STDS_DIR)

    standards_dict = None
    if standards_dict_path is not None:
        standards_dict_path = str(standards_dict_path)
        try:
            eds_file = EDSStandardsFile.from_json_file(
                standards_dict_path, meas_type=None, beam_energy_keV=None
            )
            standards_dict = eds_file.to_standards_dict()
        except Exception as e:
            raise ValueError(
                f"Could not load standards dictionary from {standards_dict_path}: {e}"
            )

    if els_substrate is None:
        els_substrate = dflts.substrate_els

    results: List[Optional[EMXSp_Composition_Analyzer]] = []
    comp_analyzer: Optional[EMXSp_Composition_Analyzer] = None

    # Track values from the last successfully loaded ledger so we can generate
    # a standards-coverage report at the end without re-opening files.
    last_microscope_ID: Optional[str] = None
    last_measurement_type: Optional[str] = None
    last_beam_energy: Optional[float] = None

    for std_id in std_ids:
        try:
            std_path = get_sample_dir(stds_path, std_id)
        except Exception as e:
            logging.warning("Failed to get sample directory for %s: %s", std_id, e)
            continue

        ledger_path = os.path.join(
            std_path, f"{cnst.LEDGER_FILENAME}{cnst.LEDGER_FILEEXT}"
        )

        print("\n \n")
        print_double_separator()
        print_double_separator()
        logging.info(f"Standard '{std_id}'")

        try:
            ledger = load_sample_ledger(ledger_path)
            configs = {
                cnst.MICROSCOPE_CFG_KEY: ledger.configs.microscope_cfg,
                cnst.SAMPLE_CFG_KEY: ledger.configs.sample_cfg,
                cnst.MEASUREMENT_CFG_KEY: ledger.configs.measurement_cfg,
                cnst.SAMPLESUBSTRATE_CFG_KEY: ledger.configs.sample_substrate_cfg,
                cnst.PLOT_CFG_KEY: ledger.configs.plot_cfg,
            }
            if ledger.quantifications:
                active_quant_id = ledger.active_quant
                active_quant_config = next(
                    (
                        quant_config
                        for quant_config in ledger.quantifications
                        if quant_config.quantification_id == active_quant_id
                    ),
                    ledger.quantifications[-1],
                )
                configs[cnst.QUANTIFICATION_CFG_KEY] = config_classes_dict[
                    cnst.QUANTIFICATION_CFG_KEY
                ](**active_quant_config.options)
                active_clustering_analysis = (
                    active_quant_config.get_active_clustering_analysis()
                )
                active_clustering_config = (
                    active_clustering_analysis.config
                    if active_clustering_analysis is not None
                    else None
                )
                if active_clustering_config is not None:
                    configs[cnst.CLUSTERING_CFG_KEY] = active_clustering_config
            else:
                configs[cnst.QUANTIFICATION_CFG_KEY] = config_classes_dict[
                    cnst.QUANTIFICATION_CFG_KEY
                ]()
            if ledger.configs.measurement_cfg.powder_meas_cfg is not None:
                configs[cnst.POWDER_MEASUREMENT_CFG_KEY] = (
                    ledger.configs.measurement_cfg.powder_meas_cfg
                )
            if ledger.configs.measurement_cfg.bulk_meas_cfg is not None:
                configs[cnst.BULK_MEASUREMENT_CFG_KEY] = (
                    ledger.configs.measurement_cfg.bulk_meas_cfg
                )
        except FileNotFoundError:
            logging.warning(
                f"Could not find {ledger_path}. Skipping sample '{std_id}'."
            )
            continue
        except Exception as e:
            logging.warning(
                f"Error loading {ledger_path}. Skipping sample '{std_id}': {e}"
            )
            continue

        sample_processing_time_start = time.time()

        # Retrieve configuration objects for this sample
        try:
            microscope_cfg       = configs[cnst.MICROSCOPE_CFG_KEY]
            sample_cfg           = configs[cnst.SAMPLE_CFG_KEY]
            measurement_cfg      = configs[cnst.MEASUREMENT_CFG_KEY]
            sample_substrate_cfg = configs[cnst.SAMPLESUBSTRATE_CFG_KEY]
            quant_cfg            = configs[cnst.QUANTIFICATION_CFG_KEY]
            clustering_cfg       = configs.get(cnst.CLUSTERING_CFG_KEY)
            plot_cfg             = configs[cnst.PLOT_CFG_KEY]
            powder_meas_cfg      = configs.get(cnst.POWDER_MEASUREMENT_CFG_KEY, None)
            bulk_meas_cfg        = configs.get(cnst.BULK_MEASUREMENT_CFG_KEY, None)
            exp_stds_cfg         = configs.get(cnst.EXP_STD_MEASUREMENT_CFG_KEY, None)
        except KeyError as e:
            logging.warning(
                f"Missing configuration '{e.args[0]}' in {ledger_path}. "
                f"Skipping sample '{std_id}'."
            )
            continue

        if clustering_cfg is None:
            clustering_cfg = ClusteringConfig()
        if quant_cfg is None:
            quant_cfg = config_classes_dict[cnst.QUANTIFICATION_CFG_KEY]()
        exp_stds_cfg_from_ledger = getattr(
            ledger.configs.measurement_cfg, "exp_stds_cfg", None
        )
        if exp_stds_cfg_from_ledger is None:
            logging.warning(
                f"No ExpStandardsConfig found in ledger for '{std_id}'. "
                f"Skipping — cannot fit a standard without a reference formula."
            )
            continue
        configs[cnst.EXP_STD_MEASUREMENT_CFG_KEY] = exp_stds_cfg_from_ledger

        if min_bckgrnd_cnts is not None:
            clustering_cfg.min_bckgrnd_cnts = min_bckgrnd_cnts
        if use_instrument_background is not None:
            quant_cfg.use_instrument_background = use_instrument_background
        if spectrum_lims is not None:
            quant_cfg.spectrum_lims = spectrum_lims

        # Track info for the post-loop standards report
        try:
            last_microscope_ID = getattr(microscope_cfg, "ID", last_microscope_ID)
            last_measurement_type = getattr(measurement_cfg, "type", last_measurement_type)
            last_beam_energy = getattr(measurement_cfg, "beam_energy_keV", last_beam_energy)
        except Exception:
            pass

    # --- Run Composition Analyzer (single instantiation, fit mode)
        comp_analyzer = EMXSp_Composition_Analyzer(
            microscope_cfg=microscope_cfg,
            sample_id=std_id,
            sample_cfg=sample_cfg,
            measurement_cfg=measurement_cfg,
            sample_substrate_cfg=sample_substrate_cfg,
            quant_cfg=quant_cfg,
            initial_clustering_cfg=clustering_cfg,
            powder_meas_cfg=powder_meas_cfg,
            bulk_meas_cfg=bulk_meas_cfg,
            plot_cfg=plot_cfg,
            exp_stds_cfg=exp_stds_cfg,
            is_acquisition=False,
            development_mode=False,
            standards_dict=standards_dict,
            output_filename_suffix=output_filename_suffix,
            verbose=verbose,
            results_dir=std_path,
        )

        try:
            # Always run fitting (run_fitting=True) since this runner is for
            # batch fitting from existing ledgers
            comp_analyzer.run_exp_std_fit(
                run_fitting=True,
                force_refitting = force_refitting,
                update_std_library=update_std_library,
            )
            results.append(comp_analyzer)
        except Exception as e:
            results.append(None)
            logging.exception(
                f"Sample '{std_id}': fitting failed: {e}"
            )
            continue

        total_process_time = (time.time() - sample_processing_time_start) / 60
        print_double_separator()
        logging.info(
            f"Sample '{std_id}' successfully fitted in {total_process_time:.1f} min."
        )

    # Generate a standards-coverage report only when standards library was updated.
    if update_std_library and last_microscope_ID is not None \
            and last_measurement_type is not None and last_beam_energy is not None:
        try:
            from autoemx.runners.extract_experimental_standards_details import (  # type: ignore
                extract_experimental_standards_details,
            )

            if stds_path is not None:
                standards_dir = Path(stds_path).expanduser().resolve()
            else:
                standards_dir = (
                    Path(__file__).resolve().parents[1]
                    / cnst.CALIBS_DIR
                    / last_microscope_ID
                )

            standards_json_filename = (
                f"{last_measurement_type}_{cnst.STD_FILENAME}_"
                f"{int(last_beam_energy):d}keV.json"
            )
            standards_json_path = standards_dir / standards_json_filename

            report_kwargs: Dict[str, Any] = {
                "microscope_ID": last_microscope_ID,
                "voltage": last_beam_energy,
                "report_output_dir": str(standards_dir),
                "print_report": True,
                "print_per_peak": False,
            }
            if standards_json_path.exists():
                report_kwargs["standards_json_path"] = str(standards_json_path)

            report_path = extract_experimental_standards_details(**report_kwargs)
            logging.info(
                "Standards details report generated at: %s", report_path
            )
        except Exception as e:
            logging.warning(
                "Could not generate standards details report after fitting: %s", e
            )

    return results