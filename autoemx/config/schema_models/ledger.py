#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import time
from pathlib import Path, PureWindowsPath, PurePosixPath
from typing import Any, Dict, List, Literal, Optional

from packaging.version import InvalidVersion, Version
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from autoemx._version import __version__ as AUTOEMX_VERSION

import autoemx.utils.constants as cnst
from autoemx.config.runtime_configs import (
    BulkMeasurementConfig,
    ExpStandardsConfig,
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from .acquisition import ParticleInfo, SpectrumEntry
from .particle_stats import ParticleStatsProgress, ParticleStatsResult
from .quantification import QuantificationConfig, QuantificationResult


class LedgerConfigs(BaseModel):
    """Typed config bundle persisted in the ledger."""

    microscope_cfg: MicroscopeConfig
    sample_cfg: SampleConfig
    measurement_cfg: MeasurementConfig
    sample_substrate_cfg: SampleSubstrateConfig
    plot_cfg: PlotConfig

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def drop_legacy_quant_cfg(cls, data: Any) -> Any:
        if isinstance(data, dict):
            cleaned = dict(data)
            cleaned.pop("quant_cfg", None)
            return cleaned
        return data

    @model_validator(mode="before")
    @classmethod
    def move_legacy_exp_std_w_frs(cls, data: Any) -> Any:
        """Legacy compatibility: move the standard's mass fractions from exp_stds_cfg to sample_cfg.

        Old ledgers stored them in measurement_cfg.exp_stds_cfg.w_frs.
        """
        if not isinstance(data, dict):
            return data
        measurement_cfg = data.get("measurement_cfg")
        sample_cfg = data.get("sample_cfg")
        if not isinstance(measurement_cfg, dict) or not isinstance(sample_cfg, dict):
            return data
        exp_stds_cfg = measurement_cfg.get("exp_stds_cfg")
        if not isinstance(exp_stds_cfg, dict):
            return data

        exp_stds_cfg = dict(exp_stds_cfg)
        legacy_w_frs = exp_stds_cfg.pop("w_frs", None)
        if legacy_w_frs and exp_stds_cfg.get("is_exp_std_measurement") and not sample_cfg.get("w_frs"):
            sample_cfg = {**sample_cfg, "w_frs": legacy_w_frs}
        measurement_cfg = {**measurement_cfg, "exp_stds_cfg": exp_stds_cfg}
        return {**data, "measurement_cfg": measurement_cfg, "sample_cfg": sample_cfg}

    @model_validator(mode="after")
    def set_exp_std_w_frs(self) -> "LedgerConfigs":
        """Set the composition of an experimental standard from its formula, when not stored yet."""
        exp_stds_cfg = self.measurement_cfg.exp_stds_cfg
        if exp_stds_cfg is not None and exp_stds_cfg.is_exp_std_measurement and not self.sample_cfg.w_frs:
            self.sample_cfg.w_frs = exp_stds_cfg.formula_w_frs()
        return self

    @property
    def powder_meas_cfg(self) -> Optional[PowderMeasurementConfig]:
        return self.measurement_cfg.powder_meas_cfg

    @property
    def bulk_meas_cfg(self) -> Optional[BulkMeasurementConfig]:
        return self.measurement_cfg.bulk_meas_cfg

    @property
    def exp_stds_cfg(self) -> Optional[ExpStandardsConfig]:
        return self.measurement_cfg.exp_stds_cfg


class LedgerVersionError(ValueError):
    """Raised when a ledger written by a newer AutoEMX version cannot be read by the installed one."""


def _newer_ledger_version_message(payload: Any) -> Optional[str]:
    """Explain which AutoEMX version to use, if the ledger was written by a newer one than installed."""
    ledger_version = payload.get("autoemx_version") if isinstance(payload, dict) else None
    if not ledger_version:
        return None
    try:
        is_newer = Version(str(ledger_version)) > Version(AUTOEMX_VERSION)
    except InvalidVersion:
        return None
    if not is_newer:
        return None
    return (
        f"This ledger was written by AutoEMX {ledger_version}, but AutoEMX {AUTOEMX_VERSION} is installed. "
        f"Update AutoEMX, or read it in a separate environment: pip install autoemx=={ledger_version}"
    )


class SampleLedger(BaseModel):
    """Single-sample ledger with per-spectrum pointers and configuration bundle."""

    # Version of AutoEMX that last saved the ledger. Ledgers saved before this field existed get the current
    # version, since the current version reads them (through the legacy compatibility code).
    autoemx_version: str = AUTOEMX_VERSION
    sample_id: str
    sample_path: str
    configs: LedgerConfigs
    spectra: List[SpectrumEntry]
    particles: List[ParticleInfo] = Field(default_factory=list)
    particle_stats: Optional[ParticleStatsResult] = None
    particle_stats_progress: Optional[ParticleStatsProgress] = None
    active_quant: Optional[int] = None
    quantifications: List[QuantificationConfig] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def move_legacy_known_phase_mixture(cls, data: Any) -> Any:
        """Legacy compatibility: is_known_powder_mixture_meas (powder settings) is now the quantification setting
        is_known_phase_mixture, stored in the options of each quantification run."""
        # Only saved ledgers (plain JSON data) can contain the legacy setting
        if not isinstance(data, dict) or not isinstance(data.get("configs"), dict):
            return data
        measurement_cfg = data["configs"].get("measurement_cfg")
        powder_meas_cfg = measurement_cfg.get("powder_meas_cfg") if isinstance(measurement_cfg, dict) else None
        if not isinstance(powder_meas_cfg, dict) or "is_known_powder_mixture_meas" not in powder_meas_cfg:
            return data
        is_known_phase_mixture = bool(powder_meas_cfg["is_known_powder_mixture_meas"])
        quantifications = []
        for quant_config in data.get("quantifications") or []:
            if isinstance(quant_config, dict) and isinstance(quant_config.get("options"), dict):
                options = quant_config["options"]
                if "is_known_phase_mixture" not in options:
                    quant_config = {**quant_config, "options": {**options, "is_known_phase_mixture": is_known_phase_mixture}}
            quantifications.append(quant_config)
        return {**data, "quantifications": quantifications}

    @field_validator("sample_path")
    @classmethod
    def validate_sample_path(cls, v: str) -> str:
        # Check if it is a valid absolute path on EITHER Windows or POSIX
        is_windows_abs = PureWindowsPath(v).is_absolute()
        is_posix_abs = PurePosixPath(v).is_absolute()
        
        if not (is_windows_abs or is_posix_abs):
            raise ValueError("sample_path must be an absolute path")
            
        return v

    @field_validator("active_quant", mode="before")
    @classmethod
    def validate_active_quant_input(cls, v: Optional[Any]) -> Optional[Any]:
        if isinstance(v, bool):
            raise ValueError("active_quant must be a non-negative integer or None")
        return v

    @field_validator("active_quant")
    @classmethod
    def validate_active_quant(cls, v: Optional[int]) -> Optional[int]:
        if v is not None and v < 0:
            raise ValueError("active_quant must be a non-negative integer or None")
        return v

    @model_validator(mode="after")
    def validate_ledger_integrity(self) -> "SampleLedger":
        particle_ids = [particle.id for particle in self.particles]
        if len(particle_ids) != len(set(particle_ids)):
            raise ValueError("ParticleInfo.id values must be unique within a ledger")
        for name, stats in (("particle_stats", self.particle_stats), ("particle_stats_progress", self.particle_stats_progress)):
            if stats is None:
                continue
            band_ids = [pid for band in stats.bands for pid in band.particle_ids]
            unknown_ids = sorted(set(getattr(stats, "particle_ids", []) + band_ids) - set(particle_ids))
            if unknown_ids:
                raise ValueError(f"{name} references particle ids not in the ledger: {unknown_ids}")

        config_ids = [cfg.quantification_id for cfg in self.quantifications]
        if len(config_ids) != len(set(config_ids)):
            raise ValueError("quantification_id values must be unique within a ledger")

        if config_ids:
            if self.active_quant is None:
                self.active_quant = config_ids[-1]
            if self.active_quant not in set(config_ids):
                raise ValueError("active_quant must reference an existing quantification_id")
        elif self.active_quant is not None:
            raise ValueError("active_quant must be None when no quantification configs are available")

        config_id_set = set(config_ids)
        for spectrum_index, spectrum in enumerate(self.spectra):
            if spectrum.spectrum_relpath is not None:
                resolved = (Path(self.sample_path) / spectrum.spectrum_relpath).resolve()
                sample_root = Path(self.sample_path).resolve()
                try:
                    resolved.relative_to(sample_root)
                except ValueError as exc:
                    raise ValueError(
                        f"Spectrum at index {spectrum_index} has spectrum_relpath outside sample_path"
                    ) from exc

            if spectrum.instrument_background_relpath is not None:
                resolved = (Path(self.sample_path) / spectrum.instrument_background_relpath).resolve()
                sample_root = Path(self.sample_path).resolve()
                try:
                    resolved.relative_to(sample_root)
                except ValueError as exc:
                    raise ValueError(
                        f"Spectrum at index {spectrum_index} has instrument_background_relpath outside sample_path"
                    ) from exc

            result_ids = set()
            for result in spectrum.quantification_results:
                if result.quantification_id not in config_id_set:
                    raise ValueError(
                        f"Spectrum at index {spectrum_index} contains a quantification result "
                        f"with unknown quantification_id '{result.quantification_id}'"
                    )

                if result.quantification_id in result_ids:
                    raise ValueError(
                        f"Spectrum at index {spectrum_index} contains duplicate "
                        f"quantification_id '{result.quantification_id}'"
                    )
                result_ids.add(result.quantification_id)

        return self

    @staticmethod
    def _load_counts_from_pointer_file(file_path: Path) -> List[float]:
        """Load counts from an external spectrum pointer file."""
        def _read_once() -> List[float]:
            suffix = file_path.suffix.lower()

            if suffix in cnst.EMSA_SPECTRUM_EXTENSIONS:
                counts: List[float] = []
                in_data_section = False
                with file_path.open("r", encoding="utf-8") as f:
                    for raw_line in f:
                        line = raw_line.strip()
                        if not line:
                            continue
                        # Accept both '#SPECTRUM' and '##SPECTRUM' as start of data
                        if line.startswith("#SPECTRUM") or line.startswith("##SPECTRUM"):
                            in_data_section = True
                            continue
                        if line.startswith("#"):
                            continue
                        if not in_data_section:
                            continue

                        # Accept only counts (no index). Reject lines with commas.
                        token = line.rstrip(",")
                        if "," in token:
                            raise ValueError(f"Pointer file '{file_path}' contains invalid spectrum data line with comma: '{line}'")
                        try:
                            counts.append(float(token))
                        except ValueError:
                            continue

                if not counts:
                    raise ValueError(f"No spectrum counts found in pointer file '{file_path}'")
                return counts

            if suffix == ".json":
                with file_path.open("r", encoding="utf-8") as f:
                    payload = json.load(f)
                values = payload.get("spectrum_vals")
                if not isinstance(values, list) or not values:
                    raise ValueError(f"JSON pointer file '{file_path}' does not contain non-empty 'spectrum_vals'")
                return [float(v) for v in values]

            raise ValueError(f"Unsupported spectrum pointer extension '{suffix}' for file '{file_path}'")

        for attempt in range(2):
            try:
                return _read_once()
            except KeyboardInterrupt:
                if attempt == 0:
                    time.sleep(0.05)
                    continue
                raise
        raise RuntimeError(f"Failed to load spectrum counts from '{file_path}' after retries")

    def append_quantification_config(self, quant_config: QuantificationConfig) -> None:
        """Append a quantification config, enforcing unique config identifiers."""
        if any(existing.quantification_id == quant_config.quantification_id for existing in self.quantifications):
            raise ValueError(
                f"quantification_id '{quant_config.quantification_id}' already exists in this ledger"
            )
        self.quantifications.append(quant_config)
        self.active_quant = quant_config.quantification_id

    def append_quantification_result(self, spectrum_index: int, result: QuantificationResult) -> None:
        """Append a quantification result to the requested spectrum entry."""
        if spectrum_index < 0 or spectrum_index >= len(self.spectra):
            raise IndexError(f"spectrum_index {spectrum_index} is out of range")
        if result.quantification_id not in {cfg.quantification_id for cfg in self.quantifications}:
            raise ValueError(
                f"Quantification result references unknown quantification_id '{result.quantification_id}'"
            )
        if any(
            existing.quantification_id == result.quantification_id
            for existing in self.spectra[spectrum_index].quantification_results
        ):
            raise ValueError(
                f"Spectrum at index {spectrum_index} already contains quantification_id "
                f"'{result.quantification_id}'"
            )
        self.spectra[spectrum_index].quantification_results.append(result)

    @classmethod
    def from_json_file(cls, file_path: str | Path) -> "SampleLedger":
        """Load and validate a ledger from JSON in one call.

        Raises
        ------
        LedgerVersionError
            If the ledger cannot be read because it was written by a newer AutoEMX version.
        """
        path = Path(file_path)
        with path.open("r", encoding="utf-8") as f:
            payload = json.load(f)
        try:
            return cls.model_validate(payload)
        except ValidationError as e:
            message = _newer_ledger_version_message(payload)
            if message is not None:
                raise LedgerVersionError(f"{message}\n(ledger: {path})") from e
            raise

    def to_dict(
        self,
        config_key_style: Literal["cfg", "config"] = "cfg",
        *,
        include_configs: bool = True,
    ) -> Dict[str, Any]:
        """Serialize ledger to a dict and choose config key suffix style."""
        payload = self.model_dump(mode="json")
        if not include_configs:
            payload.pop("configs", None)
            return payload

        if config_key_style == "cfg":
            return payload
        if config_key_style != "config":
            raise ValueError("config_key_style must be either 'cfg' or 'config'")

        cfg_block = payload.get("configs", {})
        renamed_cfg_block: Dict[str, Any] = {}
        for key, value in cfg_block.items():
            if key.endswith("_cfg"):
                renamed_cfg_block[key[: -len("_cfg")] + "_config"] = value
            else:
                renamed_cfg_block[key] = value
        payload["configs"] = renamed_cfg_block
        return payload

    def to_json_file(
        self,
        file_path: str | Path,
        *,
        config_key_style: Literal["cfg", "config"] = "cfg",
        include_configs: bool = True,
        indent: int = 2,
    ) -> None:
        """Write ledger to JSON in one call, recording the AutoEMX version that saved it."""
        self.autoemx_version = AUTOEMX_VERSION
        path = Path(file_path)
        payload = self.to_dict(config_key_style=config_key_style, include_configs=include_configs)
        path.parent.mkdir(parents=True, exist_ok=True)
        # Per-process temporary file, so concurrent writers never share it. Atomic replace avoids
        # partial reads if another process opens the ledger mid-write.
        tmp_path = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        try:
            with tmp_path.open("w", encoding="utf-8") as f:
                json.dump(payload, f, indent=indent, allow_nan=False)
                f.write("\n")
            os.replace(tmp_path, path)
        except Exception:
            if tmp_path.exists():
                tmp_path.unlink()
            raise
