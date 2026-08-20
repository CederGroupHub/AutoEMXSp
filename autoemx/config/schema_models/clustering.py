#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import json
from typing import Any, ClassVar, Dict, List, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


def _canonicalize_json_value(value: Any) -> Any:
    """Recursively normalize values to deterministic JSON-compatible primitives."""
    if isinstance(value, dict):
        return {
            str(k): _canonicalize_json_value(v)
            for k, v in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_canonicalize_json_value(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _collect_payload_differences(
    *,
    left: Any,
    right: Any,
    path: str,
    out: Dict[str, Dict[str, Any]],
) -> None:
    """Recursively collect payload differences keyed by a dotted path."""
    if isinstance(left, dict) and isinstance(right, dict):
        all_keys = sorted(set(left) | set(right))
        for key in all_keys:
            next_path = f"{path}.{key}" if path else str(key)
            if key not in left:
                out[next_path] = {"old": None, "new": right[key]}
                continue
            if key not in right:
                out[next_path] = {"old": left[key], "new": None}
                continue
            _collect_payload_differences(
                left=left[key],
                right=right[key],
                path=next_path,
                out=out,
            )
        return

    if isinstance(left, list) and isinstance(right, list):
        if left != right:
            out[path or "root"] = {"old": left, "new": right}
        return

    if left != right:
        out[path or "root"] = {"old": left, "new": right}


class DBSCANParams(BaseModel):
    """Parameters for DBSCAN clustering.

    Only relevant when :attr:`ClusteringConfig.method` is ``"dbscan"``. Grouped in
    a dedicated sub-model so method-specific knobs do not pollute the top-level
    clustering config and remain individually typed and validated.
    """

    eps: float = 0.05
    min_samples: int = 3
    metric: str = "euclidean"

    model_config = ConfigDict(extra="forbid")

    @field_validator("eps")
    @classmethod
    def validate_eps(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0:
            raise ValueError("DBSCAN eps must be a positive, finite number")
        return float(value)

    @field_validator("min_samples")
    @classmethod
    def validate_min_samples(cls, value: int) -> int:
        if value < 1:
            raise ValueError("DBSCAN min_samples must be >= 1")
        return int(value)

    @field_validator("metric")
    @classmethod
    def validate_metric(cls, value: str) -> str:
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("DBSCAN metric cannot be empty")
        return normalized


class ClusteringConfig(BaseModel):
    """Configuration for clustering of compositions and their filtering."""

    ALLOWED_METHODS: ClassVar[Tuple[str, ...]] = ("kmeans", "dbscan")

    clustering_id: int = 0
    method: str = "kmeans"
    features: str = "at_fr"
    k_forced: Optional[int] = None
    k_resolved: Optional[int] = None
    k_finding_method: str = "silhouette"
    max_k: int = 6
    ref_formulae: List[str] = Field(default_factory=list)
    do_matrix_decomposition: bool = True
    max_analytical_error_percent: Optional[float] = 5.0
    min_bckgrnd_cnts: Optional[float] = 5.0
    quant_flags_accepted: List[int] = Field(
        default_factory=lambda: [0, -1],
        description=(
            "Quantification flags kept for clustering. "
            "See the user documentation page 'Quantification flags (quant_flag)'."
        ),
    )
    dbscan: DBSCANParams = Field(default_factory=DBSCANParams)

    model_config = ConfigDict(extra="forbid")

    @field_validator("clustering_id", mode="before")
    @classmethod
    def validate_clustering_id_input(cls, v: Any) -> Any:
        if isinstance(v, bool):
            raise ValueError("clustering_id must be a non-negative integer")
        return v

    @field_validator("clustering_id")
    @classmethod
    def validate_clustering_id(cls, value: int) -> int:
        if value < 0:
            raise ValueError("clustering_id must be non-negative")
        return value

    @field_validator("method", "features", "k_finding_method")
    @classmethod
    def validate_non_empty_strings(cls, value: str) -> str:
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("Clustering string fields cannot be empty")
        return normalized

    @field_validator("method")
    @classmethod
    def validate_method(cls, value: str) -> str:
        normalized = str(value).strip().lower()
        if normalized not in cls.ALLOWED_METHODS:
            raise ValueError(
                f"Clustering method must be one of {cls.ALLOWED_METHODS}, got '{value}'."
            )
        return normalized

    @field_validator("max_k")
    @classmethod
    def validate_max_k(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_k must be positive")
        return value

    @field_validator("k_forced", "k_resolved")
    @classmethod
    def validate_k(cls, value: Optional[int]) -> Optional[int]:
        if value is not None and value <= 0:
            raise ValueError("k values must be positive when provided")
        return value

    @model_validator(mode="after")
    def validate_k_semantics(self) -> "ClusteringConfig":
        forced_key = "forced"
        if self.k_forced is not None and self.k_finding_method != forced_key:
            self.k_finding_method = forced_key
        if self.k_forced is None and self.k_finding_method == forced_key:
            raise ValueError("k_finding_method='forced' requires k_forced to be set")
        return self

    @field_validator("max_analytical_error_percent")
    @classmethod
    def validate_max_analytical_error_percent(cls, value: Optional[float]) -> Optional[float]:
        if value is not None and not np.isfinite(value):
            raise ValueError("max_analytical_error_percent must be finite when provided")
        return value

    @field_validator("min_bckgrnd_cnts")
    @classmethod
    def validate_min_bckgrnd_cnts(cls, value: Optional[float]) -> Optional[float]:
        if value is not None and value < 0:
            raise ValueError("min_bckgrnd_cnts must be non-negative or None")
        return value

    @field_validator("ref_formulae")
    @classmethod
    def validate_ref_formulae(cls, values: List[str]) -> List[str]:
        normalized: List[str] = []
        seen = set()
        for value in values:
            formula = str(value).strip()
            if not formula:
                continue
            if formula not in seen:
                normalized.append(formula)
                seen.add(formula)
        return normalized

    @field_validator("quant_flags_accepted")
    @classmethod
    def validate_quant_flags_accepted(cls, values: List[int]) -> List[int]:
        normalized: List[int] = []
        seen = set()
        for value in values:
            numeric = int(value)
            if numeric not in seen:
                normalized.append(numeric)
                seen.add(numeric)
        return normalized

    def fingerprint_payload(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "method": self.method,
            "features": self.features,
            "k_forced": self.k_forced,
            "k_finding_method": self.k_finding_method,
            "max_k": self.max_k,
            "ref_formulae": sorted(self.ref_formulae),
            "do_matrix_decomposition": self.do_matrix_decomposition,
            "max_analytical_error_percent": self.max_analytical_error_percent,
            "min_bckgrnd_cnts": self.min_bckgrnd_cnts,
            "quant_flags_accepted": sorted(self.quant_flags_accepted),
        }
        # Only fold method-specific params into the fingerprint when they are
        # actually in use, so existing k-means configs keep their current hash.
        if self.method == "dbscan":
            payload["dbscan"] = self.dbscan.model_dump()
        return payload

    def fingerprint(self) -> str:
        canonical_payload = _canonicalize_json_value(self.fingerprint_payload())
        canonical_json = json.dumps(canonical_payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()

    def fingerprint_differences(self, other: "ClusteringConfig") -> Dict[str, Dict[str, Any]]:
        left_payload = _canonicalize_json_value(self.fingerprint_payload())
        right_payload = _canonicalize_json_value(other.fingerprint_payload())

        differences: Dict[str, Dict[str, Any]] = {}
        _collect_payload_differences(
            left=left_payload,
            right=right_payload,
            path="",
            out=differences,
        )
        return dict(sorted(differences.items()))


class ClusteringResult(BaseModel):
    """Cluster analysis artifacts produced by one analysis run."""

    centroids: List[List[float]] = Field(default_factory=list)
    els_std_dev_per_cluster: List[List[float]] = Field(default_factory=list)
    centroids_other_fr: List[List[float]] = Field(default_factory=list)
    els_std_dev_per_cluster_other_fr: List[List[float]] = Field(default_factory=list)
    n_points_per_cluster: List[int] = Field(default_factory=list)
    wcss_per_cluster: List[float] = Field(default_factory=list)
    rms_dist_cluster: List[float] = Field(default_factory=list)
    rms_dist_cluster_other_fr: List[float] = Field(default_factory=list)
    refs_assigned_rows: List[Dict[str, Any]] = Field(default_factory=list)
    wcss: float
    sil_score: float
    tot_n_points: int
    clusters_assigned_mixtures: List[Any] = Field(default_factory=list)

    model_config = ConfigDict(extra="forbid")


class ClusteringAnalysis(BaseModel):
    """Bundle of the clustering config and the resulting analysis artifacts."""

    config: ClusteringConfig
    result: Optional[ClusteringResult] = None

    model_config = ConfigDict(extra="forbid")
