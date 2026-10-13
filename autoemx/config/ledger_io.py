#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convenience I/O helpers for sample ledgers."""

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, cast

import autoemx.utils.constants as cnst
from autoemx._logging import get_logger

from autoemx.config.ledger_schemas import (  # type: ignore
    AcquisitionDetails,
    ClusteringConfig,
    LedgerConfigs,
    QuantificationConfig,
    QuantificationResult,
    SampleLedger,
    SpectrumEntry,
)
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

logger = get_logger(__name__)

# Missing quantified spectra already reported, to warn once per file (the ledger is synchronized repeatedly)
_reported_missing_quantified_spectra: set[tuple[str, str]] = set()


def _extract_spectrum_id(pointer_file: Path) -> str:
    """Extract spectrum id from a pointer filename."""
    stem = pointer_file.stem
    if stem.startswith(cnst.SPECTRUM_FILENAME_PREFIX):
        return stem[len(cnst.SPECTRUM_FILENAME_PREFIX):]
    return stem


def _list_spectrum_pointer_files(spectra_dir: Path) -> List[Path]:
    """
    List the spectrum files of a spectra folder (``spectrum_<id>.<ext>``), one per spectrum id, in deterministic
    order (numeric ids first). Files of other names, and manufacturer background files, are ignored. When a spectrum
    has files of several formats, the one of highest priority (SPECTRUM_POINTER_EXT_PRIORITY) is used.
    """
    if not spectra_dir.exists():
        return []

    allowed_ext = set(cnst.SPECTRUM_POINTER_EXTENSIONS)
    ext_priority = dict(cnst.SPECTRUM_POINTER_EXT_PRIORITY)
    selected_by_id = {}

    for path in spectra_dir.iterdir():
        if not path.is_file() or path.suffix.lower() not in allowed_ext:
            continue
        if not path.stem.startswith(cnst.SPECTRUM_FILENAME_PREFIX) or path.stem.endswith(cnst.SPECTRUM_MAN_BACKGROUND_SUFFIX):
            continue

        spectrum_id = _extract_spectrum_id(path)
        existing = selected_by_id.get(spectrum_id)
        if existing is None:
            selected_by_id[spectrum_id] = path
            continue

        existing_priority = ext_priority.get(existing.suffix.lower(), 99)
        current_priority = ext_priority.get(path.suffix.lower(), 99)
        if current_priority < existing_priority:
            selected_by_id[spectrum_id] = path

    def _sort_key(path: Path):
        spectrum_id = _extract_spectrum_id(path)
        if spectrum_id.isdigit():
            return (0, int(spectrum_id), path.name)
        return (1, spectrum_id.lower(), path.name)

    return sorted(selected_by_id.values(), key=_sort_key)


def _background_relpath(spectrum_id: str) -> str:
    """Relative path of the manufacturer background vector file of a spectrum."""
    return str(Path(cnst.SPECTRA_DIR, f"{cnst.SPECTRUM_FILENAME_PREFIX}{spectrum_id}"
                                     f"{cnst.SPECTRUM_MAN_BACKGROUND_SUFFIX}{cnst.VECTOR_FILEEXT}"))


def _load_realtime_from_pointer_file(pointer_path: Path) -> Optional[float]:
    """Read REALTIME from EMSA-like headers when available."""
    if pointer_path.suffix.lower() not in cnst.EMSA_SPECTRUM_EXTENSIONS:
        return None

    try:
        with pointer_path.open("r", encoding="utf-8") as file_obj:
            for raw_line in file_obj:
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


def _build_spectrum_entry(
    sample_root: Path,
    pointer_file: Path,
    acquisition_details: Optional[AcquisitionDetails] = None,
    quantification_results: Optional[List[QuantificationResult]] = None,
) -> SpectrumEntry:
    """
    Build the ledger entry of a spectrum file of ``<sample_root>/spectra``.

    The live time is read from the file header (1.0 s if missing or not positive), and the manufacturer background
    vector file of the spectrum is referenced if it exists.
    """
    spectrum_id = _extract_spectrum_id(pointer_file)
    sample_root = sample_root.resolve()
    pointer_relpath = str(pointer_file.resolve().relative_to(sample_root).as_posix())

    try:
        counts = SampleLedger._load_counts_from_pointer_file(pointer_file.resolve())
        total_counts = int(round(sum(float(value) for value in counts)))
    except Exception:
        total_counts = 0

    live_time = _load_realtime_from_pointer_file(pointer_file)
    if live_time is not None and not (math.isfinite(live_time) and live_time > 0):
        live_time = None

    background_path = sample_root / _background_relpath(spectrum_id)
    background_relpath = background_path.relative_to(sample_root).as_posix() if background_path.exists() else None

    return SpectrumEntry(
        spectrum_id=spectrum_id,
        total_counts=total_counts,
        live_acquisition_time=live_time if live_time is not None else 1.0,
        acquisition_details=acquisition_details or AcquisitionDetails(frame_id=None, particle_id=None, spot_coordinates=None),
        spectrum_relpath=pointer_relpath,
        instrument_background_relpath=background_relpath,
        quantification_results=list(quantification_results or []),
    )


def create_sample_ledger(sample_dir: str | Path, sample_id: str, configs: LedgerConfigs) -> SampleLedger:
    """New ledger of a sample, with its configurations and no spectra, particles or quantification runs (not saved)."""
    return SampleLedger(
        sample_id=sample_id,
        sample_path=os.path.abspath(sample_dir),
        configs=configs,
        spectra=[],
        particles=[],
        quantifications=[],
        active_quant=None,
    )


def ingest_spectra(ledger: SampleLedger) -> tuple[int, int]:
    """Synchronize ledger spectra entries with files in ``<ledger.sample_path>/spectra``.

    Appends pointer files whose spectrum id is not already present in
    ``ledger.spectra`` and removes ledger entries whose spectrum id no longer
    exists on disk, so that wrong spectra can be deleted before quantification.
    Entries that already have quantification results are kept even when their
    file is missing, since deleting a quantified spectrum is assumed to be a mistake.

    Returns
    -------
    tuple[int, int]
        Number of newly ingested spectra and number of removed spectra.
    """
    sample_root = Path(ledger.sample_path)
    spectra_dir = sample_root / cnst.SPECTRA_DIR
    pointer_files = _list_spectrum_pointer_files(spectra_dir)
    available_ids = {_extract_spectrum_id(pointer_file) for pointer_file in pointer_files}
    available_relpaths = {
        pointer_file.resolve().relative_to(sample_root.resolve()).as_posix()
        for pointer_file in pointer_files
    }

    def _is_on_disk(entry: SpectrumEntry) -> bool:
        return entry.spectrum_relpath in available_relpaths or (
            entry.spectrum_relpath in (None, "")
            and entry.spectrum_id not in (None, "")
            and str(entry.spectrum_id) in available_ids
        )

    original_count = len(ledger.spectra)
    kept_entries = []
    for entry in ledger.spectra:
        if _is_on_disk(entry):
            kept_entries.append(entry)
        elif entry.quantification_results:
            kept_entries.append(entry)
            report_key = (str(sample_root), str(entry.spectrum_relpath or entry.spectrum_id))
            if report_key not in _reported_missing_quantified_spectra:
                _reported_missing_quantified_spectra.add(report_key)
                logger.warning(
                    f"⚠️ Spectrum file {entry.spectrum_relpath or entry.spectrum_id} not found, but spectrum "
                    f"#{entry.spectrum_id} has already been quantified: keeping it and its results in the ledger."
                )
    ledger.spectra[:] = kept_entries
    n_removed = original_count - len(ledger.spectra)

    existing_ids = {
        str(entry.spectrum_id)
        for entry in ledger.spectra
        if entry.spectrum_id not in (None, "")
    }

    n_ingested = 0
    for pointer_file in pointer_files:
        spectrum_id = _extract_spectrum_id(pointer_file)
        if spectrum_id in existing_ids:
            continue
        ledger.spectra.append(_build_spectrum_entry(sample_root, pointer_file))
        existing_ids.add(spectrum_id)
        n_ingested += 1

    return n_ingested, n_removed


def load_sample_ledger(file_path: str | Path) -> SampleLedger:
    """Load a sample ledger, with full legacy bootstrap when ledger.json is missing.

    If the ledger file does not exist, this uses legacy Data.csv + config JSON files
    to bootstrap a populated ledger (configs, spectra, and quantification payloads).
    """
    ledger_path = Path(file_path)

    if ledger_path.exists():
        ledger = SampleLedger.from_json_file(ledger_path)
        sample_root = str(ledger_path.parent.resolve())
        sample_path_changed = ledger.sample_path != sample_root
        if sample_path_changed:
            ledger.sample_path = sample_root

        n_ingested, n_removed = ingest_spectra(ledger)
        if sample_path_changed or n_ingested > 0 or n_removed > 0:
            ledger.to_json_file(ledger_path)

        if n_ingested > 0 or n_removed > 0:
            change_parts = []
            if n_ingested > 0:
                change_parts.append(f"ingested {n_ingested} new spectrum file{'s' if n_ingested != 1 else ''}")
            if n_removed > 0:
                change_parts.append(f"removed {n_removed} missing spectrum entr{'ies' if n_removed != 1 else 'y'}")
            print(
                f"Loaded sample ledger for sample '{ledger.sample_id}' from {ledger_path} "
                f"and {'; '.join(change_parts)}."
            )
        else:
            print(f"Loaded sample ledger for sample '{ledger.sample_id}' from {ledger_path}")
        return ledger

    sample_result_dir = ledger_path.parent
    from autoemx.utils.legacy.legacy_backfill import load_ledger_configs_from_legacy_json
    from autoemx.utils.legacy.legacy_config_loader import has_legacy_data_csv
    from autoemx.utils.legacy.ledger_bootstrap import ( # type: ignore
        build_legacy_background_pointer_writer, # type: ignore
        build_legacy_import_quantification_config, # type: ignore
        build_legacy_msa_pointer_writer, # type: ignore
        load_or_create_ledger_with_legacy_data_csv, # type: ignore
    )

    legacy_ledger_configs_raw = load_ledger_configs_from_legacy_json(str(sample_result_dir))
    if legacy_ledger_configs_raw is None:
        legacy_config_candidates = [
            sample_result_dir / f"{cnst.CONFIG_FILENAME}.json",
            sample_result_dir / f"{cnst.ACQUISITION_INFO_FILENAME}.json",
        ]
        if any(path.exists() for path in legacy_config_candidates):
            raise FileNotFoundError(
                f"No ledger file found at '{ledger_path}' and legacy config JSON in "
                f"'{sample_result_dir}' could not be parsed into LedgerConfigs."
            )
        raise FileNotFoundError(
            f"No ledger file found at '{ledger_path}' and no legacy config JSON "
            f"was found in '{sample_result_dir}'."
        )
    legacy_ledger_configs = cast(LedgerConfigs, legacy_ledger_configs_raw)

    if has_legacy_data_csv(str(sample_result_dir)):
        legacy_quant_config = build_legacy_import_quantification_config(
            sample_result_dir=str(sample_result_dir),
            ledger_configs=legacy_ledger_configs,
        )
        xperchan = legacy_ledger_configs.microscope_cfg.bin_width *1000
        offset = legacy_ledger_configs.microscope_cfg.energy_zero  *1000
        resolve_or_create_spectrum_pointer = build_legacy_msa_pointer_writer(str(sample_result_dir), xperchan=xperchan, offset=offset)
        write_background_pointer = build_legacy_background_pointer_writer(str(sample_result_dir))

        ledger = load_or_create_ledger_with_legacy_data_csv(
            sample_result_dir=str(sample_result_dir),
            sample_id=sample_result_dir.name,
            microscope_id=getattr(legacy_ledger_configs.microscope_cfg, "ID", None),
            use_instrument_background=bool(
                (legacy_quant_config.options or {}).get("use_instrument_background", False)
            ),
            default_ledger_configs=legacy_ledger_configs,
            resolve_or_create_spectrum_pointer=resolve_or_create_spectrum_pointer,
            write_background_pointer=write_background_pointer,
        )
        ledger.sample_path = str(sample_result_dir.resolve())
        print(
            f"Loaded sample ledger for sample '{ledger.sample_id}' from legacy Data.csv/configs "
            f"in '{sample_result_dir}' (ledger.json not found)."
        )
        return ledger

    legacy_quant_config = build_legacy_import_quantification_config(
        sample_result_dir=str(sample_result_dir),
        ledger_configs=legacy_ledger_configs,
    )
    ledger = SampleLedger(
        sample_id=sample_result_dir.name,
        sample_path=str(sample_result_dir.resolve()),
        configs=legacy_ledger_configs,
        spectra=[],
        quantifications=[legacy_quant_config],
        active_quant=int(legacy_quant_config.quantification_id),
    )
    print(
        f"Loaded legacy sample configuration for sample '{ledger.sample_id}' "
        f"from '{sample_result_dir}' (ledger.json and Data.csv not found)."
    )
    return ledger


def get_active_quantification_config(ledger: Optional[SampleLedger]) -> Optional[QuantificationConfig]:
    """Return the active quantification run of the ledger (the last one if none is active), or None if none."""
    if ledger is None or not ledger.quantifications:
        return None
    return next(
        (config for config in ledger.quantifications if config.quantification_id == ledger.active_quant),
        ledger.quantifications[-1],
    )


@dataclass
class LedgerRuntimeConfigs:
    """Configurations of a sample ledger, as used by the runners.

    The acquisition configurations are the ledger's own objects (not copies). The quantification options and the
    clustering configuration are those of the active quantification run, or defaults if the sample was never
    quantified or analysed.
    """

    microscope_cfg: MicroscopeConfig
    sample_cfg: SampleConfig
    measurement_cfg: MeasurementConfig
    sample_substrate_cfg: SampleSubstrateConfig
    plot_cfg: PlotConfig
    powder_meas_cfg: Optional[PowderMeasurementConfig]
    bulk_meas_cfg: Optional[BulkMeasurementConfig]
    quant_cfg: QuantificationOptionsConfig
    clustering_cfg: ClusteringConfig
    active_quant_config: Optional[QuantificationConfig]


def ledger_runtime_configs(ledger: SampleLedger) -> LedgerRuntimeConfigs:
    """Collect the configurations of a sample ledger needed to quantify or analyse it again."""
    active_quant_config = get_active_quantification_config(ledger)
    quant_cfg = (
        QuantificationOptionsConfig(**active_quant_config.options)
        if active_quant_config is not None
        else QuantificationOptionsConfig()
    )
    active_clustering_analysis = (
        active_quant_config.get_active_clustering_analysis() if active_quant_config is not None else None
    )
    clustering_cfg = (
        active_clustering_analysis.config if active_clustering_analysis is not None else ClusteringConfig()
    )
    configs = ledger.configs
    return LedgerRuntimeConfigs(
        microscope_cfg=configs.microscope_cfg,
        sample_cfg=configs.sample_cfg,
        measurement_cfg=configs.measurement_cfg,
        sample_substrate_cfg=configs.sample_substrate_cfg,
        plot_cfg=configs.plot_cfg,
        powder_meas_cfg=configs.measurement_cfg.powder_meas_cfg,
        bulk_meas_cfg=configs.measurement_cfg.bulk_meas_cfg,
        quant_cfg=quant_cfg,
        clustering_cfg=clustering_cfg,
        active_quant_config=active_quant_config,
    )


def save_sample_ledger(
    ledger: SampleLedger,
    file_path: str | Path,
    *,
    indent: int = 2,
) -> None:
    """Save a sample ledger JSON in one call."""
    ledger.to_json_file(file_path, indent=indent)
