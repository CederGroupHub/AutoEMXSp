#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hashlib
import json
from typing import Any, ClassVar, Dict, List, Optional, Tuple

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

import autoemx.config.defaults as dflt


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

    ``eps`` is a distance in the clustering geometry, whose scale differs by roughly an
    order of magnitude between fractions (Euclidean) and CLR coordinates (Aitchison).
    When left as ``None``, the geometry-specific default in ``DEFAULT_EPS`` is used.
    """

    DEFAULT_EPS: ClassVar[Dict[str, float]] = {"euclidean": 0.05, "aitchison": 0.3}

    eps: Optional[float] = None
    min_samples: int = 3
    metric: str = "euclidean"

    model_config = ConfigDict(extra="forbid")

    @field_validator("eps")
    @classmethod
    def validate_eps(cls, value: Optional[float]) -> Optional[float]:
        if value is None:
            return None
        if not np.isfinite(value) or value <= 0:
            raise ValueError("DBSCAN eps must be a positive, finite number")
        return float(value)

    def resolved_eps(self, geometry: str) -> float:
        """Return ``eps`` if set, otherwise the default for ``geometry`` ("euclidean" or "aitchison")."""
        return self.eps if self.eps is not None else self.DEFAULT_EPS[geometry]

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


class AitchisonParams(BaseModel):
    """Parameters for clustering in Aitchison geometry.

    Only relevant when :attr:`ClusteringConfig.geometry` is ``"aitchison"`` or ``"auto"``.
    Before the CLR transform, compositions are closed and low values are replaced
    multiplicatively (Martín-Fernández et al., 2003) with 0.65 x ``detection_limit_percent``:

    - exact zeros of every element (a measured zero means "below detection");
    - all values below the detection limit of *trace* elements, i.e. elements whose
      median fraction across the clustered compositions is below 2 x the detection limit.

    Measured sub-detection-limit values of major elements are kept. This prevents noisy
    trace elements from dominating log-ratio distances without collapsing real low
    concentrations of major elements.

    With ``geometry="auto"``, Euclidean geometry is used instead of Aitchison when either
    condition below holds on the clustered compositions, since log-ratios then tend to
    amplify noise or over-split:

    - at most two elements are present (log-ratio space is one-dimensional);
    - in at least ``auto_max_near_zero_fraction`` of the spectra, some major element
      (median >= 2 x the detection limit) is below ``auto_near_zero_percent``, e.g. phases
      with largely disjoint elements, or mixture lines reaching an end-member.
    """

    detection_limit_percent: float = 0.5
    auto_near_zero_percent: float = 1.0
    auto_max_near_zero_fraction: float = 0.10

    model_config = ConfigDict(extra="forbid")

    @field_validator("detection_limit_percent")
    @classmethod
    def validate_detection_limit_percent(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0 or value >= 100:
            raise ValueError("Aitchison detection_limit_percent must be a finite number in (0, 100)")
        return float(value)

    @field_validator("auto_near_zero_percent")
    @classmethod
    def validate_auto_near_zero_percent(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0 or value >= 100:
            raise ValueError("Aitchison auto_near_zero_percent must be a finite number in (0, 100)")
        return float(value)

    @field_validator("auto_max_near_zero_fraction")
    @classmethod
    def validate_auto_max_near_zero_fraction(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0 or value > 1:
            raise ValueError("Aitchison auto_max_near_zero_fraction must be a finite number in (0, 1]")
        return float(value)


class ClusterMergeParams(BaseModel):
    """Parameters for merging k-means clusters that are pieces of one continuous population.

    Only used when :attr:`ClusteringConfig.auto_merge_clusters` is ``True``, with k-means and an
    automatically chosen number of clusters. k-means can cut a continuous spread of compositions,
    e.g. a mixture line between two phases, into several clusters. After clustering, two clusters
    are merged when all of the following hold, in the clustering geometry:

    - their compositions projected onto the axis joining the two cluster centers form a single
      mode (Hartigan's dip test, p-value above ``dip_alpha``);
    - the empty gap between them along that axis is below ``max_gap_ratio`` times the spread
      (standard deviation) of the larger cluster;
    - DBSCAN connects them through dense neighbourhoods, with ``min_samples = dbscan_min_samples``
      and ``eps = dbscan_eps_factor`` times the median distance of the compositions to their
      ``dbscan_min_samples``-th nearest neighbour.

    Pairs are merged one at a time, most compatible first, and the criteria are re-evaluated after
    each merge. Clusters are never split.
    """

    dip_alpha: float = 0.05
    max_gap_ratio: float = 1.0
    dbscan_min_samples: int = 5
    dbscan_eps_factor: float = 3.0

    model_config = ConfigDict(extra="forbid")

    @field_validator("dip_alpha")
    @classmethod
    def validate_dip_alpha(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0 or value >= 1:
            raise ValueError("Cluster merge dip_alpha must be in (0, 1)")
        return float(value)

    @field_validator("max_gap_ratio", "dbscan_eps_factor")
    @classmethod
    def validate_positive(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Cluster merge parameters must be positive, finite numbers")
        return float(value)

    @field_validator("dbscan_min_samples")
    @classmethod
    def validate_min_samples(cls, value: int) -> int:
        if value < 2:
            raise ValueError("Cluster merge dbscan_min_samples must be >= 2")
        return int(value)


class MixtureParams(BaseModel):
    """Parameters for decomposing clusters into mixtures of candidate phases (``ref_formulae``).

    Clusters are fitted as non-negative combinations of candidate phases (NNLS), scored by the
    exponential reconstruction error ``mean(exp(recon_error_alpha * |X - WH|) - 1)``, which is
    mapped to a confidence score ``exp(-error^2 / (2 * conf_sigma^2))``.

    Adding phases always lowers the reconstruction error, so only the fewest phases that the
    shape of the cluster requires are used: a binary mixture spreads along a line between two
    phases (with a width set by measurement noise), a ternary one over a plane.
    Pairs of candidate phases are always tested, and listed when their error is below
    ``max_recon_error_binary`` (1, confidence ~0.14; only pairs below ``max_recon_error``
    explain the cluster, weaker pairs are listed for inspection). Combinations of 3 up to ``max_n_phases`` phases are only tested
    when no combination with fewer phases reaches an error below ``max_recon_error``, stopping
    at the first number of phases that does. ``max_n_phases = 2`` restricts the analysis to
    binary mixtures. The default ``max_recon_error`` (0.4, confidence ~0.73) separates
    single-phase standards, known binary mixtures and ternary mixtures measured on a Phenom XL.

    Combinations whose reconstruction errors differ by less than ``equivalent_recon_error_tol``
    explain the cluster equally well, e.g. collinear candidate phases such as CaO and Ca4Ta2O9
    in a CaO-Ta2O5 mixture. Among these, the combination with fewer phases and then with the
    phases closest to each other in composition space (smallest sum of pairwise distances) is
    ranked first.

    All mixtures are saved in the ledger. ``Clusters.csv`` shows a mixture only if its confidence
    is at least ``min_reported_conf_ratio`` of the best one in the cluster, and it is among the
    first ``max_reported_mixtures`` or its confidence is within ``report_within_conf_ratio`` of the
    best one. A note in ``Clusters.csv`` gives the number of mixtures only saved in the ledger. With
    ``collapse_equivalent_mixtures``, mixtures whose phases span the same mixing line or plane as a
    better-ranked mixture (within ``equivalent_span_tol``, in fractions) are equivalent decompositions
    and are not shown either, e.g. pairs of Sr-Ta oxides that all lie on the SrO-Ta2O5 line.

    Clusters are not decomposed if they are considered single-phase: their RMS distance from the
    centroid is below ``single_phase_max_rms_dist`` and they match a candidate phase with
    confidence above ``single_phase_min_ref_conf``. If no mixture of candidate phases reaches a
    confidence of ``nmf_min_mixture_conf``, the cluster is also decomposed into two phases of
    unknown composition with free NMF.
    """

    max_n_phases: int = 4
    max_recon_error: float = 0.4
    max_recon_error_binary: float = 1.0
    recon_error_alpha: float = 15.0
    conf_sigma: float = 0.5
    single_phase_max_rms_dist: float = 0.03
    single_phase_min_ref_conf: float = 0.5
    nmf_min_mixture_conf: float = 0.5
    equivalent_recon_error_tol: float = 0.01
    max_reported_mixtures: int = 5
    report_within_conf_ratio: float = 0.9
    min_reported_conf_ratio: float = 0.5
    collapse_equivalent_mixtures: bool = True
    equivalent_span_tol: float = 0.005

    model_config = ConfigDict(extra="forbid")

    @field_validator("max_n_phases")
    @classmethod
    def validate_max_n_phases(cls, value: int) -> int:
        if value < 2:
            raise ValueError("Mixture max_n_phases must be >= 2")
        return int(value)

    @field_validator("max_recon_error", "max_recon_error_binary", "recon_error_alpha", "conf_sigma",
                     "single_phase_max_rms_dist")
    @classmethod
    def validate_positive(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Mixture parameters must be positive, finite numbers")
        return float(value)

    @field_validator("max_reported_mixtures")
    @classmethod
    def validate_max_reported_mixtures(cls, value: int) -> int:
        if value < 1:
            raise ValueError("Mixture max_reported_mixtures must be >= 1")
        return int(value)

    @field_validator("report_within_conf_ratio", "min_reported_conf_ratio")
    @classmethod
    def validate_conf_ratio(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0 or value > 1:
            raise ValueError("Mixture confidence ratios must be in (0, 1]")
        return float(value)

    @field_validator("equivalent_span_tol")
    @classmethod
    def validate_equivalent_span_tol(cls, value: float) -> float:
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Mixture equivalent_span_tol must be a positive, finite number")
        return float(value)

    @field_validator("equivalent_recon_error_tol")
    @classmethod
    def validate_non_negative(cls, value: float) -> float:
        if not np.isfinite(value) or value < 0:
            raise ValueError("Mixture equivalent_recon_error_tol must be a non-negative, finite number")
        return float(value)

    @field_validator("single_phase_min_ref_conf", "nmf_min_mixture_conf")
    @classmethod
    def validate_confidence(cls, value: float) -> float:
        if not np.isfinite(value) or value < 0 or value > 1:
            raise ValueError("Mixture confidence thresholds must be in [0, 1]")
        return float(value)


class ClusteringConfig(BaseModel):
    """Configuration for clustering of compositions and their filtering."""

    ALLOWED_METHODS: ClassVar[Tuple[str, ...]] = ("kmeans", "dbscan")
    ALLOWED_GEOMETRIES: ClassVar[Tuple[str, ...]] = ("euclidean", "aitchison", "auto")
    LEGACY_GEOMETRY: ClassVar[str] = "euclidean"
    LEGACY_AUTO_MERGE_CLUSTERS: ClassVar[bool] = False

    clustering_id: int = 0
    method: str = "kmeans"
    # Geometry in which clustering distances are computed. "aitchison" clusters
    # CLR-transformed compositions; reported centroids/statistics stay in fraction space.
    # "auto" picks euclidean or aitchison from the measured compositions (see AitchisonParams).
    # New configs default to "auto"; configs saved before this field existed and configs
    # created from legacy data are "euclidean" (LEGACY_GEOMETRY), preserving their results.
    geometry: str = "auto"
    features: str = "at_fr"
    k_forced: Optional[int] = None
    k_resolved: Optional[int] = None
    k_finding_method: str = "silhouette"
    max_k: int = dflt.max_n_clusters
    ref_formulae: List[str] = Field(default_factory=list)
    do_matrix_decomposition: bool = True
    max_analytical_error_percent: Optional[float] = dflt.max_analytical_error_percent
    quant_flags_accepted: List[int] = Field(
        default_factory=lambda: [0, -1],
        description=(
            "Quantification flags kept for clustering. "
            "See the user documentation page 'Quantification flags (quant_flag)'."
        ),
    )
    dbscan: DBSCANParams = Field(default_factory=DBSCANParams)
    aitchison: AitchisonParams = Field(default_factory=AitchisonParams)
    # Merge k-means clusters that are pieces of one continuous population (see ClusterMergeParams).
    # New configs default to True; configs saved before this field existed and configs created from
    # legacy data are False (LEGACY_AUTO_MERGE_CLUSTERS), preserving their results.
    auto_merge_clusters: bool = True
    cluster_merge: ClusterMergeParams = Field(default_factory=ClusterMergeParams)
    mixture: MixtureParams = Field(default_factory=MixtureParams)

    model_config = ConfigDict(extra="forbid")

    @model_validator(mode="before")
    @classmethod
    def drop_legacy_min_bckgrnd_cnts(cls, data: Any) -> Any:
        """Legacy compatibility: min_bckgrnd_cnts is a quantification setting (QuantificationOptionsConfig).

        Its value in old clustering configs is moved to their quantification run by QuantificationConfig.
        """
        if isinstance(data, dict) and "min_bckgrnd_cnts" in data:
            data = dict(data)
            data.pop("min_bckgrnd_cnts")
        return data

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

    @field_validator("geometry")
    @classmethod
    def validate_geometry(cls, value: str) -> str:
        normalized = str(value).strip().lower()
        if normalized not in cls.ALLOWED_GEOMETRIES:
            raise ValueError(
                f"Clustering geometry must be one of {cls.ALLOWED_GEOMETRIES}, got '{value}'."
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
            "quant_flags_accepted": sorted(self.quant_flags_accepted),
        }
        # Only fold method-specific params into the fingerprint when they are
        # actually in use, so existing k-means configs keep their current hash.
        if self.method == "dbscan":
            payload["dbscan"] = self.dbscan.model_dump()
        if self.auto_merge_clusters:
            payload["auto_merge_clusters"] = True
            payload["cluster_merge"] = self.cluster_merge.model_dump()
        if self.geometry != "euclidean":
            payload["geometry"] = self.geometry
            payload["aitchison"] = self.aitchison.model_dump()
        if self.mixture != MixtureParams():
            payload["mixture"] = self.mixture.model_dump()
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

    @model_validator(mode="before")
    @classmethod
    def backfill_legacy_geometry(cls, data: Any) -> Any:
        """
        Configs saved before ``geometry`` existed were clustered in Euclidean geometry, and configs
        saved before ``auto_merge_clusters`` existed without merging clusters.
        """
        if isinstance(data, dict) and isinstance(data.get("config"), dict):
            config = dict(data["config"])
            config.setdefault("geometry", ClusteringConfig.LEGACY_GEOMETRY)
            config.setdefault("auto_merge_clusters", ClusteringConfig.LEGACY_AUTO_MERGE_CLUSTERS)
            data = {**data, "config": config}
        return data
