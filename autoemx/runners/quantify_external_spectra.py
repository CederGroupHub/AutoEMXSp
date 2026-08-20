#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import and quantify externally-acquired X-ray spectra.

This module copies ``.msa`` / ``.msg`` / ``.emsa`` spectrum files from an arbitrary directory
into the canonical ``sample_id/spectra/`` subfolder under ``results_path``, builds
a complete ``ledger.json`` (identical to what the acquisition runner writes), and
then hands off to :func:`batch_quantify_and_analyze` for quantification and
optional clustering analysis.

Typical workflow
----------------
1. Run ``quantify_external_spectra.py`` → copies spectra, writes ``ledger.json``,
   runs quantification, and (optionally) clustering analysis.
2. Re-run ``run_analysis.py`` at any time to repeat clustering with different
   parameters.

Created on Fri Apr 24 2026

@author: Andrea
"""

import logging
import os
import shutil
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional

import autoemx.calibrations as calibs
import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
from autoemx.config.ledger_schemas import ClusteringConfig  # type: ignore
from autoemx.config.runtime_configs import (
    BulkMeasurementConfig,
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    QuantificationOptionsConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze
from autoemx.utils import print_double_separator

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

__all__ = ["quantify_external_spectra"]

_SUPPORTED_EXTENSIONS = set(cnst.EMSA_SPECTRUM_EXTENSIONS)


def _discover_spectra_files(spectra_dir: str) -> List[Path]:
    """Return sorted list of EMSA spectrum files found directly inside *spectra_dir*."""
    source = Path(spectra_dir)
    if not source.is_dir():
        raise FileNotFoundError(
            f"spectra_dir does not exist or is not a directory: '{spectra_dir}'"
        )
    return sorted(
        p for p in source.iterdir()
        if p.is_file() and p.suffix.lower() in _SUPPORTED_EXTENSIONS
    )


def _clear_destination_spectra_files(spectra_dest_dir: str) -> int:
    """Remove canonical imported spectra files from destination and return count."""
    dest = Path(spectra_dest_dir)
    if not dest.is_dir():
        return 0

    removed = 0
    for p in dest.iterdir():
        if not p.is_file():
            continue
        if p.suffix.lower() not in _SUPPORTED_EXTENSIONS:
            continue
        if not p.name.startswith(cnst.SPECTRUM_FILENAME_PREFIX):
            continue
        p.unlink()
        removed += 1
    return removed


def quantify_external_spectra(
    samples: List[Dict[str, Any]],
    microscope_ID: str = dflt.microscope_ID,
    microscope_type: str = dflt.microscope_type,
    measurement_type: str = dflt.measurement_type,
    measurement_mode: str = dflt.measurement_mode,
    beam_energy: float = 15.0,
    els_substrate: Optional[List[str]] = None,
    sample_substrate_type: str = "Ctape",
    sample_substrate_shape: str = "circle",
    sample_substrate_width_mm: float = 12,
    sample_type: str = "bulk",
    sample_halfwidth: float = 2.9,
    quantification_method: str = dflt.quantification_method,
    use_project_specific_std_dict: bool = False,
    use_instrument_background: bool = dflt.use_instrument_background,
    max_analytical_error_percent: float = 5.0,
    min_bckgrnd_cnts: float = 5.0,
    quant_flags_accepted: Optional[List[int]] = None,
    max_n_clusters: int = 6,
    show_unused_comps_clust: bool = True,
    interrupt_fits_bad_spectra: bool = False,
    num_CPU_cores: Optional[int] = None,
    run_analysis: bool = True,
    powder_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    bulk_meas_cfg_kwargs: Optional[Dict[str, Any]] = None,
    results_path: Optional[str] = None,
    overwrite_existing: bool = False,
    standards_dict: Optional[dict] = None,
    verbose: bool = True,
) -> List[EMXSp_Composition_Analyzer]:
    """
    Copy externally-acquired spectra into the AutoEMX sample folder layout,
    build a complete ``ledger.json``, then run quantification and (optionally)
    clustering analysis via :func:`batch_quantify_and_analyze`.

    Parameters
    ----------
    samples : list of dict
        Each entry describes one sample.  Recognised keys:

        - ``'ID'`` (str, required): Sample identifier; the output folder is named after this.
        - ``'els'`` (list of str, required): Element symbols expected in the sample.
        - ``'spectra_dir'`` (str, optional): Path to a folder containing ``.msa`` / ``.msg``
          / ``.emsa`` spectrum files.  Files are copied and renamed ``spectrum_0<ext>``,
          ``spectrum_1<ext>``, … in sort order.  When omitted (or set to ``None``), the
          spectra are expected to already be present at ``<results_path>/<ID>/spectra/``
          in the canonical naming.
        - ``'cnd'`` (list of str, optional): Candidate phase formulae for clustering
          (e.g. ``['PbMoO4']``).
        - ``'type'`` (str, optional): Override ``sample_type`` for this entry.
        - ``'half_width_mm'`` (float, optional): Override ``sample_halfwidth``.

    microscope_ID : str
        Identifier of the microscope calibration folder.  Default: ``'PhenomXL'``.
    microscope_type : str
        Microscope type (``'SEM'``).
    measurement_type : str
        Measurement type (``'EDS'``).
    measurement_mode : str
        Measurement mode (``'point'``); selects the detector-channel calibration.
    beam_energy : float
        Electron beam energy in keV.  Must match the energy used during acquisition.
    els_substrate : list of str, optional
        Elements in the sample substrate.  Default: ``['C', 'O', 'Al']``.
    sample_substrate_type : str
        Substrate material type (``'Ctape'`` or ``'None'``).
    sample_substrate_shape : str
        Stub shape (``'circle'`` or ``'square'``).
    sample_substrate_width_mm : float
        Stub lateral dimension in mm.
    sample_type : str
        Sample geometry for quantification corrections (``'powder'``, ``'bulk'``,
        ``'bulk_rough'``, ``'powder_continuous'``).
    sample_halfwidth : float
        Half-width of the sample area in mm.
    quantification_method : str
        Quantification algorithm.  Only ``'PB'`` is currently implemented.
    use_project_specific_std_dict : bool
        Load standards from the project folder rather than global calibrations.
    use_instrument_background : bool
        Use instrument background during fitting when available.
    max_analytical_error_percent : float
        Compositions above this analytical-error threshold (w%) are excluded
        from clustering.
    min_bckgrnd_cnts : int
        Minimum background counts under the reference peak for a spectrum to be
        accepted for quantification.
    quant_flags_accepted : list of int, optional
        Quantification flags included during clustering.  Default: ``[0, -1]``.
        See the user documentation page "Quantification flags (quant_flag)".
    max_n_clusters : int
        Maximum number of clusters.
    show_unused_comps_clust : bool
        Show discarded compositions in the clustering plot.
    interrupt_fits_bad_spectra : bool
        Abort fits early when quality indicators are poor (speeds up batch runs).
    num_CPU_cores : int, optional
        Parallel fitting cores.  ``None`` uses half the available cores.
    run_analysis : bool
        Run clustering / statistical analysis after quantification.
    powder_meas_cfg_kwargs : dict, optional
        Extra kwargs forwarded to :class:`PowderMeasurementConfig`.
    bulk_meas_cfg_kwargs : dict, optional
        Extra kwargs forwarded to :class:`BulkMeasurementConfig`.
    results_path : str, optional
        Root directory for sample sub-folders.  Defaults to ``<cwd>/Results``.
    overwrite_existing : bool
        Re-copy external spectra files into ``sample_id/spectra`` even when
        ``ledger.json`` already exists. The existing ledger is never deleted;
        quantification history is preserved.
    standards_dict : dict, optional
        Custom dictionary of reference PB values; ``None`` loads the defaults.
    verbose : bool
        Print progress information.

    Returns
    -------
    list of EMXSp_Composition_Analyzer
        One analyzer per successfully processed sample (from
        :func:`batch_quantify_and_analyze`).
    """
    if els_substrate is None:
        els_substrate = ["C", "O", "Al"]
    if quant_flags_accepted is None:
        quant_flags_accepted = [0, -1]
    if results_path is None:
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)

    # ------------------------------------------------------------------
    # Shared config objects (identical for every sample in the batch)
    # ------------------------------------------------------------------
    microscope_cfg = MicroscopeConfig(ID=microscope_ID, type=microscope_type)

    # Populate detector channel parameters so they are stored in the ledger.
    calibs.load_microscope_calibrations(
        microscope_ID, measurement_mode, load_detector_channel_params=True
    )
    meas_modes_calibs = calibs.detector_channel_params
    microscope_cfg.energy_zero = meas_modes_calibs[measurement_mode][cnst.OFFSET_KEY]
    microscope_cfg.bin_width   = meas_modes_calibs[measurement_mode][cnst.SCALE_KEY]

    measurement_cfg = MeasurementConfig(
        type=measurement_type,
        mode=measurement_mode,
        beam_energy_keV=beam_energy,
    )
    sample_substrate_cfg = SampleSubstrateConfig(
        elements=els_substrate,
        type=sample_substrate_type,
        shape=sample_substrate_shape,
        auto_detection=False,
        stub_w_mm=sample_substrate_width_mm,
    )
    quant_cfg = QuantificationOptionsConfig(
        method=quantification_method,
        use_project_specific_std_dict=use_project_specific_std_dict,
        use_instrument_background=use_instrument_background,
    )
    plot_cfg = PlotConfig(show_unused_comps_clust=show_unused_comps_clust)
    powder_meas_cfg = PowderMeasurementConfig(**(powder_meas_cfg_kwargs or {}))
    bulk_meas_cfg   = BulkMeasurementConfig(**(bulk_meas_cfg_kwargs or {}))

    # ------------------------------------------------------------------
    # Phase 1 — Per-sample: copy spectra + build ledger
    # ------------------------------------------------------------------
    ingested_ids: List[str] = []

    for sample in samples:
        sample_id: str        = sample["ID"]
        elements: List[str]   = sample["els"]
        spectra_dir: Optional[str] = sample.get("spectra_dir")
        ref_formulae: List[str]    = sample.get("cnd", [])
        smpl_type: str             = sample.get("type", sample_type)
        smpl_halfwidth: float      = sample.get("half_width_mm", sample_halfwidth)

        print_double_separator()
        logging.info("Preparing sample '%s' for quantification.", sample_id)

        sample_dir      = os.path.join(results_path, sample_id)
        ledger_path     = os.path.join(sample_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)
        spectra_dest_dir = os.path.join(sample_dir, cnst.SPECTRA_DIR)

        if os.path.exists(ledger_path) and not overwrite_existing:
            logging.info(
                "Ledger already exists at '%s'. Skipping ingestion "
                "(set overwrite_existing=True to replace it).",
                ledger_path,
            )
            ingested_ids.append(sample_id)
            continue

        # ---- Copy spectra (only when a source directory is given and differs) ----
        if spectra_dir is not None:
            src_resolved  = Path(spectra_dir).resolve()
            dest_resolved = Path(spectra_dest_dir).resolve()
            if src_resolved == dest_resolved:
                logging.info(
                    "spectra_dir is identical to the destination; skipping copy."
                )
            else:
                try:
                    source_files = _discover_spectra_files(spectra_dir)
                except Exception as exc:
                    logging.warning(
                        "Could not discover spectra for '%s': %s. Skipping.",
                        sample_id, exc,
                    )
                    continue

                if not source_files:
                    logging.warning(
                        "No .msa/.msg/.emsa files found in '%s'. Skipping sample '%s'.",
                        spectra_dir, sample_id,
                    )
                    continue

                os.makedirs(spectra_dest_dir, exist_ok=True)
                if overwrite_existing:
                    removed = _clear_destination_spectra_files(spectra_dest_dir)
                    if removed > 0:
                        logging.info(
                            "Removed %d existing spectrum file(s) from destination before copy.",
                            removed,
                        )
                for idx, src in enumerate(source_files):
                    dest_name = f"{cnst.SPECTRUM_FILENAME_PREFIX}{idx}{src.suffix.lower()}"
                    shutil.copy2(src, dest_resolved / dest_name)
                logging.info(
                    "Copied %d spectrum file(s) from '%s'.",
                    len(source_files), spectra_dir,
                )

        os.makedirs(spectra_dest_dir, exist_ok=True)

        # ---- Sample-specific configs ----
        sample_cfg = SampleConfig(
            elements=elements,
            type=smpl_type,
            half_width_mm=smpl_halfwidth,
        )

        n_els = len(elements)
        clustering_cfg = ClusteringConfig(
            method="kmeans",
            features=cnst.W_FR_CL_FEAT if n_els <= 2 else cnst.AT_FR_CL_FEAT,
            k_finding_method="calinski_harabasz" if n_els <= 2 else "silhouette",
            max_k=max_n_clusters,
            ref_formulae=ref_formulae,
            max_analytical_error_percent=max_analytical_error_percent,
            min_bckgrnd_cnts=min_bckgrnd_cnts,
            quant_flags_accepted=quant_flags_accepted,
        )

        # ---- Instantiate analyzer (non-acquisition mode) ----
        try:
            comp_analyzer = EMXSp_Composition_Analyzer(
                microscope_cfg=microscope_cfg,
                sample_id=sample_id,
                sample_cfg=sample_cfg,
                measurement_cfg=measurement_cfg,
                sample_substrate_cfg=sample_substrate_cfg,
                quant_cfg=quant_cfg,
                initial_clustering_cfg=clustering_cfg,
                powder_meas_cfg=powder_meas_cfg,
                bulk_meas_cfg=bulk_meas_cfg,
                plot_cfg=plot_cfg,
                is_acquisition=False,
                standards_dict=standards_dict,
                verbose=verbose,
                results_dir=sample_dir,
            )
        except Exception as exc:
            logging.warning(
                "Could not initialise analyzer for '%s': %s. Skipping.",
                sample_id, exc,
            )
            traceback.print_exc()
            continue

        # ---- Build and persist ledger ----
        # _load_or_create_ledger discovers every spectrum_*.msa in spectra/,
        # builds SpectrumEntry objects, seeds a QuantificationConfig, and writes
        # ledger.json to disk.  No manual QuantificationConfig creation needed here.
        try:
            comp_analyzer._load_or_create_ledger()
        except Exception as exc:
            logging.warning(
                "Could not build ledger for '%s': %s. Skipping.", sample_id, exc,
            )
            traceback.print_exc()
            continue

        logging.info("Ledger ready for sample '%s'.", sample_id)
        ingested_ids.append(sample_id)

    if not ingested_ids:
        logging.warning("No samples were successfully ingested. Quantification skipped.")
        return []

    # ------------------------------------------------------------------
    # Phase 2 — Quantification (and optional analysis) via batch runner
    # ------------------------------------------------------------------
    logging.info(
        "Starting quantification for %d sample(s): %s",
        len(ingested_ids), ingested_ids,
    )
    return batch_quantify_and_analyze(
        sample_IDs=ingested_ids,
        results_path=results_path,
        min_bckgrnd_cnts=min_bckgrnd_cnts,
        use_instrument_background=use_instrument_background,
        max_analytical_error=max_analytical_error_percent,
        interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
        num_CPU_cores=num_CPU_cores,
        run_analysis=run_analysis,
        use_project_specific_std_dict=use_project_specific_std_dict,
        standards_dict=standards_dict,
    )
