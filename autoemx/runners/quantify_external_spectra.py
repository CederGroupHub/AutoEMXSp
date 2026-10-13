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
from typing import Any, Dict, List, Optional, Tuple

import autoemx.calibrations as calibs
import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import create_sample_ledger, ingest_spectra, load_sample_ledger
from autoemx.config.ledger_schemas import ClusteringConfig, LedgerConfigs  # type: ignore
from autoemx.config.runtime_configs import (
    BulkMeasurementConfig,
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
from autoemx.runners.analyze_sample import analyze_sample
from autoemx.runners.batch_quantify_and_analyze import batch_quantify_and_analyze
from autoemx.utils import load_msa, print_double_separator

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


def read_energy_calibration(files: List[Path]) -> Tuple[float, float]:
    """
    Energy calibration ``(offset, width)`` in keV (energy of channel 0 and channel width),
    from the ``#OFFSET`` and ``#XPERCHAN`` lines of the EMSA headers of *files*,
    read as in the web app (``load_msa`` + ``parse_emsa_geometry``).

    Raises ValueError if the files have no calibration or different ones.
    """
    from autoemx.web.pipeline import parse_emsa_geometry

    if not files:
        raise ValueError("no spectrum files")
    calibrations = set()
    for path in files:
        _, _, metadata = load_msa(str(path))
        if "OFFSET" not in metadata or "XPERCHAN" not in metadata:
            raise ValueError(f"no #OFFSET or #XPERCHAN line in the header of '{path.name}'")
        geometry = parse_emsa_geometry(metadata)
        calibrations.add((round(geometry["det_ch_offset"], 9), round(geometry["det_ch_width"], 9)))
    if len(calibrations) > 1:
        raise ValueError(f"the spectra have {len(calibrations)} different energy calibrations")
    return calibrations.pop()


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
    beam_energy: float = dflt.beam_energy_keV,
    els_substrate: Optional[List[str]] = None,
    sample_substrate_type: str = "Ctape",
    sample_substrate_shape: str = "circle",
    sample_substrate_width_mm: float = dflt.substrate_width_mm,
    sample_type: str = "bulk",
    sample_halfwidth: float = dflt.sample_halfwidth_mm,
    quantification_method: str = dflt.quantification_method,
    use_project_specific_std_dict: bool = False,
    use_instrument_background: bool = dflt.use_instrument_background,
    max_analytical_error_percent: float = dflt.max_analytical_error_percent,
    min_bckgrnd_cnts: Optional[float] = dflt.min_bckgrnd_cnts,
    quant_flags_accepted: Optional[List[int]] = None,
    max_n_clusters: int = dflt.max_n_clusters,
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
    quantify: bool = True,
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
        Elements in the sample substrate.  Default: ``['C']`` (defaults.substrate_els).
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
    quantify : bool
        If False, only copies the spectra and builds the ledgers (no quantification nor analysis);
        the samples can then be quantified with :func:`batch_quantify_and_analyze`. Default: True.

    Returns
    -------
    list of EMXSp_Composition_Analyzer
        One analyzer per successfully quantified sample: from :func:`analyze_sample` if analysed, otherwise from
        :func:`batch_quantify_and_analyze`.
    """
    if els_substrate is None:
        els_substrate = list(dflt.substrate_els)
    if quant_flags_accepted is None:
        quant_flags_accepted = [0, -1]
    if results_path is None:
        results_path = os.path.join(os.getcwd(), cnst.RESULTS_DIR)

    # ------------------------------------------------------------------
    # Shared config objects (identical for every sample in the batch)
    # ------------------------------------------------------------------
    # The energy calibration of each sample is read from the headers of its spectra (below)
    microscope_cfg = MicroscopeConfig(ID=microscope_ID, type=microscope_type)

    measurement_cfg = MeasurementConfig(
        type=measurement_type,
        mode=measurement_mode,
        beam_energy_keV=beam_energy,
        emergence_angle=calibs.load_detector_model(microscope_ID).emergence_angle,  # Fixed by instrument geometry
    )
    sample_substrate_cfg = SampleSubstrateConfig(
        elements=els_substrate,
        type=sample_substrate_type,
        shape=sample_substrate_shape,
        auto_detection=False,
        stub_w_mm=sample_substrate_width_mm,
    )
    plot_cfg = PlotConfig(show_unused_comps_clust=show_unused_comps_clust)
    powder_meas_cfg = PowderMeasurementConfig(**(powder_meas_cfg_kwargs or {}))
    bulk_meas_cfg   = BulkMeasurementConfig(**(bulk_meas_cfg_kwargs or {}))

    # ------------------------------------------------------------------
    # Phase 1 — Per-sample: copy spectra + build ledger
    # ------------------------------------------------------------------
    ingested_ids: List[str] = []
    clustering_cfgs: Dict[str, ClusteringConfig] = {}  # Analysis settings of each sample

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

        # ---- Energy calibration, always read from the EMSA headers of the spectra ----
        try:
            energy_zero, bin_width = read_energy_calibration(_discover_spectra_files(spectra_dest_dir))
        except Exception as exc:
            logging.warning("Could not read the energy calibration of '%s': %s. Skipping.", sample_id, exc)
            continue
        sample_microscope_cfg = microscope_cfg.model_copy(
            update={"energy_zero": energy_zero, "bin_width": bin_width}
        )
        logging.info("Energy calibration from the spectra headers: offset %.4f keV, %.6f keV/channel.",
                     energy_zero, bin_width)

        # ---- Sample-specific configs ----
        sample_cfg = SampleConfig(
            elements=elements,
            type=smpl_type,
            half_width_mm=smpl_halfwidth,
        )

        n_els = len(elements)
        clustering_cfgs[sample_id] = ClusteringConfig(
            method="kmeans",
            features=cnst.W_FR_CL_FEAT if n_els <= 2 else cnst.AT_FR_CL_FEAT,
            k_finding_method="calinski_harabasz" if n_els <= 2 else "silhouette",
            max_k=max_n_clusters,
            ref_formulae=ref_formulae,
            max_analytical_error_percent=max_analytical_error_percent,
            min_bckgrnd_cnts=min_bckgrnd_cnts,
            quant_flags_accepted=quant_flags_accepted,
        )

        # ---- Build and persist ledger ----
        # As at the start of an acquisition: a new ledger with the sample configs, then every
        # spectrum_* file of spectra/ is ingested. An existing ledger keeps its quantification history.
        try:
            if os.path.exists(ledger_path):
                ledger = load_sample_ledger(ledger_path)
            else:
                ledger = create_sample_ledger(
                    sample_dir,
                    sample_id,
                    LedgerConfigs(
                        microscope_cfg=sample_microscope_cfg,
                        sample_cfg=sample_cfg,
                        measurement_cfg=measurement_cfg.model_copy(
                            update={"powder_meas_cfg": powder_meas_cfg, "bulk_meas_cfg": bulk_meas_cfg}
                        ),
                        sample_substrate_cfg=sample_substrate_cfg,
                        plot_cfg=plot_cfg,
                    ),
                )
            n_ingested, _ = ingest_spectra(ledger)
            if not ledger.spectra:
                raise ValueError(f"no spectrum files in '{spectra_dest_dir}'")
            ledger.to_json_file(ledger_path)
            logging.info("Ledger with %d spectra (%d new).", len(ledger.spectra), n_ingested)
        except Exception as exc:
            logging.warning(
                "Could not build ledger for '%s': %s. Skipping.", sample_id, exc,
            )
            traceback.print_exc()
            continue

        logging.info("Ledger ready for sample '%s'.", sample_id)
        ingested_ids.append(sample_id)

    if not quantify:
        logging.info("Ledgers ready for %d sample(s): %s. Quantification not requested.",
                     len(ingested_ids), ", ".join(ingested_ids))
        return []

    if not ingested_ids:
        logging.warning("No samples were successfully ingested. Quantification skipped.")
        return []

    # ------------------------------------------------------------------
    # Phase 2 — Quantification via batch runner
    # ------------------------------------------------------------------
    logging.info(
        "Starting quantification for %d sample(s): %s",
        len(ingested_ids), ingested_ids,
    )
    analyzers = batch_quantify_and_analyze(
        sample_IDs=ingested_ids,
        results_path=results_path,
        quantification_method=quantification_method,
        min_bckgrnd_cnts=min_bckgrnd_cnts,
        use_instrument_background=use_instrument_background,
        interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
        num_CPU_cores=num_CPU_cores,
        run_analysis=False,
        use_project_specific_std_dict=use_project_specific_std_dict,
        standards_dict=standards_dict,
    )
    if not run_analysis:
        return analyzers

    # ------------------------------------------------------------------
    # Phase 3 — Analysis with the clustering settings of each sample (candidate phases, clusters, flags)
    # ------------------------------------------------------------------
    analyzed = []
    for comp_analyzer in analyzers:
        clustering_cfg = clustering_cfgs[comp_analyzer.sample_id]
        try:
            analysis_analyzer = analyze_sample(
                comp_analyzer.sample_id,
                results_path=results_path,
                ref_formulae=clustering_cfg.ref_formulae,
                clustering_features=clustering_cfg.features,
                clustering_method=clustering_cfg.method,
                k_finding_method=clustering_cfg.k_finding_method,
                max_k=clustering_cfg.max_k,
                max_analytical_error_percent=max_analytical_error_percent,
                quant_flags_accepted=clustering_cfg.quant_flags_accepted,
                min_bckgrnd_cnts=min_bckgrnd_cnts,
                show_plots=False,
                show_unused_compositions_cluster_plot=show_unused_comps_clust,
            )
        except Exception as exc:
            logging.exception(f"Error during clustering analysis for '{comp_analyzer.sample_id}'. "
                              f"Rerun separately if needed: {exc}")
            analysis_analyzer = None
        analyzed.append(analysis_analyzer if analysis_analyzer is not None else comp_analyzer)
    return analyzed
