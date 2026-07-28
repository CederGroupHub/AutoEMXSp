#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import time
from pathlib import Path, PureWindowsPath, PurePosixPath
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

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

    @property
    def powder_meas_cfg(self) -> Optional[PowderMeasurementConfig]:
        return self.measurement_cfg.powder_meas_cfg

    @property
    def bulk_meas_cfg(self) -> Optional[BulkMeasurementConfig]:
        return self.measurement_cfg.bulk_meas_cfg

    @property
    def exp_stds_cfg(self) -> Optional[ExpStandardsConfig]:
        return self.measurement_cfg.exp_stds_cfg


class SampleLedger(BaseModel):
    """Single-sample ledger with per-spectrum pointers and configuration bundle."""

    sample_id: str
    sample_path: str
    configs: LedgerConfigs
    spectra: List[SpectrumEntry]
    particles: List[ParticleInfo] = Field(default_factory=list)
    active_quant: Optional[int] = None
    quantifications: List[QuantificationConfig] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")

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

            if suffix in {".msa", ".msg"}:
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
        """Load and validate a ledger from JSON in one call."""
        path = Path(file_path)
        with path.open("r", encoding="utf-8") as f:
            payload = json.load(f)
        return cls.model_validate(payload)

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
        """Write ledger to JSON in one call."""
        path = Path(file_path)
        payload = self.to_dict(config_key_style=config_key_style, include_configs=include_configs)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = path.with_suffix(path.suffix + ".tmp")
        try:
            with tmp_path.open("w", encoding="utf-8") as f:
                json.dump(payload, f, indent=indent, allow_nan=False)
                f.write("\n")
            os.replace(tmp_path, path)
        except Exception:
            if tmp_path.exists():
                tmp_path.unlink()
            raise
