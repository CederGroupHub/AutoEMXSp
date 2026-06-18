#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator

import autoemx.utils.constants as cnst
from autoemx.utils.legacy.legacy_standards import normalize_standards_file_payload



class StandardMeanZ(BaseModel):
    """Mean atomic number metrics stored for one standard measurement."""
    statham2016: float = Field(alias=cnst.Z_MEAN_STATHAM_KEY)
    markowicz1984: float = Field(alias=cnst.Z_MEAN_MARKOWICZ_KEY)
    mass_averaged: float = Field(alias=cnst.Z_MEAN_W_KEY)
    atomic_averaged: float = Field(alias=cnst.Z_MEAN_AT_KEY)

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    @field_validator("mass_averaged", "atomic_averaged", "statham2016", "markowicz1984", mode="before")
    @classmethod
    def round_to_2_decimals(cls, value):
        if value is None:
            return value
        try:
            return round(float(value), 2)
        except Exception:
            return value



class StandardPbSummary(BaseModel):
    """Shared PB summary fields used by fit results and persisted entries."""

    corrected_pb: float = Field(alias=cnst.COR_PB_DF_KEY)
    rel_stdev_pb_percent: float = Field(alias=cnst.REL_ER_PERCENT_PB_DF_KEY)
    measured_pb: Optional[float] = Field(default=None, alias=cnst.MEAS_PB_DF_KEY)
    stdev_pb: float = Field(alias=cnst.STDEV_PB_DF_KEY)

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    @field_validator("corrected_pb", "stdev_pb", "measured_pb", mode="before")
    @classmethod
    def round_to_1_decimal(cls, value):
        if value is None:
            return value
        try:
            return round(float(value), 1)
        except Exception:
            return value

    @field_validator("rel_stdev_pb_percent", mode="before")
    @classmethod
    def round_to_2_decimals(cls, value):
        if value is None:
            return value
        try:
            return round(float(value), 2)
        except Exception:
            return value


class StandardEntry(StandardPbSummary):
    """Single standards-library entry for one element reference line."""

    standard_id: str = Field(alias=cnst.STD_ID_KEY)
    datetime: str = Field(alias=cnst.DATETIME_KEY)
    formula: Optional[str] = Field(default=None, alias=cnst.STD_FORMULA_KEY)
    std_type: Optional[str] = Field(default=None, alias=cnst.STD_TYPE_KEY)
    use_for_mean_calc: Optional[bool] = Field(default=False, alias=cnst.STD_USE_FOR_MEAN_KEY)
    mean_z: Optional[StandardMeanZ] = Field(default=None, alias=cnst.STD_Z_KEY)

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    @field_validator("standard_id", "datetime")
    @classmethod
    def validate_non_empty_str(cls, value: str) -> str:
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("Standard entry identifiers and timestamps cannot be empty")
        return normalized


class ReferenceMean(StandardPbSummary):
    """Dedicated reference-mean values for one element reference line."""

    datetime: Optional[str] = Field(default=None, alias=cnst.DATETIME_KEY)

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    @field_validator("datetime")
    @classmethod
    def validate_optional_datetime(cls, value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("Reference mean timestamp cannot be empty when provided")
        return normalized


class StandardLine(BaseModel):
    """Structured standards payload for one element reference line."""

    entries: List[StandardEntry] = Field(default_factory=list)
    reference_mean: Optional[ReferenceMean] = None

    model_config = ConfigDict(extra="forbid", populate_by_name=True)


class StandardFitLineResult(StandardPbSummary):
    """Computed standards-fit metrics for one element reference line."""

    pb_ratios: List[float]
    peak_theoretical_energy_keV: float
    n_spectra_used: int

    model_config = ConfigDict(extra="forbid")
    def to_standard_entry(
        self,
        *,
        standard_id: str,
        datetime: str,
        formula: Optional[str] = None,
        std_type: Optional[str] = None,
        use_for_mean_calc: Optional[bool] = False,
        mean_z: Optional[StandardMeanZ] = None,
    ) -> StandardEntry:
        return StandardEntry.model_validate({
            cnst.STD_ID_KEY: standard_id,
            cnst.DATETIME_KEY: datetime,
            cnst.COR_PB_DF_KEY: self.corrected_pb,
            cnst.STDEV_PB_DF_KEY: self.stdev_pb,
            cnst.REL_ER_PERCENT_PB_DF_KEY: self.rel_stdev_pb_percent,
            cnst.MEAS_PB_DF_KEY: self.measured_pb,
            cnst.STD_FORMULA_KEY: formula,
            cnst.STD_TYPE_KEY: std_type,
            cnst.STD_USE_FOR_MEAN_KEY: use_for_mean_calc,
            cnst.STD_Z_KEY: mean_z.model_dump(mode="json", by_alias=True) if mean_z is not None else None,
        })


class StandardsFitResults(BaseModel):
    """Aggregate standards-fit output used for persistence and library updates."""

    lines: Dict[str, StandardFitLineResult] = Field(default_factory=dict)
    mean_z: Optional[StandardMeanZ] = None

    model_config = ConfigDict(extra="forbid")


class EDSStandardsFile(BaseModel):
    """Schema-backed representation of an EDS standards file for one beam energy."""

    schema_version: int = 1
    measurement_type: str
    beam_energy_keV: int
    standards_by_mode: Dict[str, Dict[str, StandardLine]] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")

    @field_validator("schema_version")
    @classmethod
    def validate_schema_version(cls, value: int) -> int:
        if value < 1:
            raise ValueError("schema_version must be >= 1")
        return value

    @field_validator("measurement_type")
    @classmethod
    def validate_measurement_type(cls, value: str) -> str:
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("measurement_type cannot be empty")
        return normalized

    @field_validator("beam_energy_keV", mode="before")
    @classmethod
    def validate_beam_energy(cls, value: Any) -> int:
        if isinstance(value, bool):
            raise ValueError("beam_energy_keV must be an integer")
        normalized = int(value)
        if normalized <= 0:
            raise ValueError("beam_energy_keV must be positive")
        return normalized

    @field_validator("standards_by_mode", mode="before")
    @classmethod
    def normalize_standards_by_mode_payload(cls, value: Any) -> Any:
        if not isinstance(value, dict):
            return value

        normalized_by_mode: Dict[str, Dict[str, Dict[str, Any]]] = {}
        for mode, lines_payload in value.items():
            if not isinstance(lines_payload, dict):
                raise ValueError("standards_by_mode entries must be dictionaries")

            normalized_lines: Dict[str, Dict[str, Any]] = {}
            for el_line, line_payload in lines_payload.items():
                if isinstance(line_payload, dict) and (
                    "entries" in line_payload or "reference_mean" in line_payload
                ):
                    normalized_lines[el_line] = line_payload
                    continue

                if not isinstance(line_payload, list):
                    raise ValueError(
                        f"Line payload for '{mode}:{el_line}' must be a list or structured dict"
                    )

                entries: List[Dict[str, Any]] = []
                reference_mean: Optional[Dict[str, Any]] = None
                for entry in line_payload:
                    if not isinstance(entry, dict):
                        raise ValueError(
                            f"Invalid standards entry for '{mode}:{el_line}': expected object"
                        )

                    if entry.get(cnst.STD_ID_KEY) == cnst.STD_MEAN_ID_KEY:
                        reference_mean = {
                            cnst.DATETIME_KEY: entry.get(cnst.DATETIME_KEY),
                            cnst.COR_PB_DF_KEY: entry.get(cnst.COR_PB_DF_KEY),
                            cnst.STDEV_PB_DF_KEY: entry.get(cnst.STDEV_PB_DF_KEY),
                            cnst.REL_ER_PERCENT_PB_DF_KEY: entry.get(cnst.REL_ER_PERCENT_PB_DF_KEY),
                        }
                        continue

                    entries.append(entry)

                line_model_payload: Dict[str, Any] = {"entries": entries}
                if reference_mean is not None and reference_mean.get(cnst.COR_PB_DF_KEY) is not None:
                    line_model_payload["reference_mean"] = reference_mean
                normalized_lines[el_line] = line_model_payload

            normalized_by_mode[mode] = normalized_lines

        return normalized_by_mode

    @classmethod
    def from_payload(
        cls,
        payload: Dict[str, Any],
        meas_type: str,
        beam_energy_keV: int,
    ) -> "EDSStandardsFile":
        normalized_payload = normalize_standards_file_payload(
            payload=payload,
            measurement_type=meas_type,
            beam_energy_keV=int(beam_energy_keV),
        )
        return cls.model_validate(normalized_payload)

    @classmethod
    def from_standards_dict(
        cls,
        standards_by_mode: Dict[str, Dict[str, List[Dict[str, Any]]]],
        meas_type: str,
        beam_energy_keV: int,
    ) -> "EDSStandardsFile":
        payload = {
            "schema_version": 1,
            "measurement_type": str(meas_type),
            "beam_energy_keV": int(beam_energy_keV),
            "standards_by_mode": standards_by_mode,
        }
        return cls.from_payload(payload, meas_type=meas_type, beam_energy_keV=beam_energy_keV)

    @classmethod
    def from_json_file(
        cls,
        file_path: str | Path,
        meas_type: str,
        beam_energy_keV: int,
    ) -> "EDSStandardsFile":
        path = Path(file_path)
        with path.open("r", encoding="utf-8") as file:
            payload = json.load(file)
        return cls.from_payload(payload, meas_type=meas_type, beam_energy_keV=beam_energy_keV)

    def to_standards_dict(self) -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
        """Return legacy-compatible per-mode dictionaries used by quantification code."""
        serialized: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
        for mode, lines in self.standards_by_mode.items():
            serialized[mode] = {}
            for el_line, line_payload in lines.items():
                line_entries = [
                    entry.model_dump(mode="json", by_alias=True, exclude_none=True)
                    for entry in line_payload.entries
                ]
                if line_payload.reference_mean is not None:
                    mean_entry = line_payload.reference_mean.model_dump(
                        mode="json",
                        by_alias=True,
                        exclude_none=True,
                    )
                    mean_entry[cnst.STD_ID_KEY] = cnst.STD_MEAN_ID_KEY
                    line_entries.append(mean_entry)

                serialized[mode][el_line] = line_entries
        return serialized

    def to_payload(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "measurement_type": self.measurement_type,
            "beam_energy_keV": self.beam_energy_keV,
            "standards_by_mode": {
                mode: {
                    el_line: line_payload.model_dump(
                        mode="json",
                        by_alias=True,
                        exclude_none=True,
                    )
                    for el_line, line_payload in lines.items()
                }
                for mode, lines in self.standards_by_mode.items()
            },
        }

    def to_json_file(self, file_path: str | Path, indent: int = 2) -> None:
        """
        Write the standards file to JSON, sorting el_line keys in each mode using el_line_sort_key logic.
        """
        from pymatgen.core.periodic_table import Element as PymatgenElement
        def el_line_sort_key(el_line):
            el, *rest = el_line.split("_")
            try:
                z = PymatgenElement(el).Z
            except Exception:
                z = 999  # Unknown elements go last
            line = rest[0] if rest else ""
            line_priority = {"Ka": 0, "La": 1, "Ma": 2, "Mz": 3}
            line_rank = line_priority.get(line, 99)
            return (z, line_rank, line)

        path = Path(file_path)
        payload = self.to_payload()
        # Sort el_line keys in standards_by_mode for each mode using el_line_sort_key
        if "standards_by_mode" in payload:
            for mode, lines in payload["standards_by_mode"].items():
                sorted_keys = sorted(lines, key=el_line_sort_key)
                sorted_lines = {k: lines[k] for k in sorted_keys}
                payload["standards_by_mode"][mode] = sorted_lines
        with path.open("w", encoding="utf-8") as file:
            json.dump(payload, file, indent=indent, allow_nan=False)
            file.write("\n")
