#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reference matching module for phase identification and mixture analysis."""

import itertools
from typing import Any, Dict, List, Optional, Tuple, cast

import cvxpy as cp
import numpy as np
import pandas as pd
from cvxpy.constraints.constraint import Constraint
from pymatgen.core.composition import Composition

import autoemx.utils.constants as cnst
from autoemx.core.composition_analysis.plotting import PlottingModule


# Transient key of mixture dictionaries: per-point molar fractions, used for plots and not saved
_POINT_MOLAR_FRS_KEY = '_point_molar_frs'


def _is_known_phase_mixture(analyzer: Any) -> bool:
	"""Whether the sample was quantified as a mixture of known phases (setting of the active quantification run)."""
	quant_config = getattr(analyzer, "current_quant_config", None)
	options = getattr(quant_config, "options", None) or {}
	if "is_known_phase_mixture" in options:
		return bool(options["is_known_phase_mixture"])
	return bool(getattr(getattr(analyzer, "quant_cfg", None), "is_known_phase_mixture", False))


def _draw_intermixing_violin_plots(analyzer: Any) -> bool:
	"""Whether to plot the intermixing of mixtures as violin plots (by default, for known phase mixtures)."""
	setting = getattr(getattr(analyzer, "plot_cfg", None), "intermixing_violin_plots", None)
	return _is_known_phase_mixture(analyzer) if setting is None else bool(setting)


class ReferenceMatchingModule:
	"""Container for reference matching algorithms extracted from the analyser."""

	# Attributes are populated by the orchestrating analyzer instance.
	ref_formulae: Optional[List[str]]
	ref_phases_df: pd.DataFrame
	ref_weights_in_mixture: np.ndarray
	detectable_els_sample: List[str]
	clustering_cfg: Any
	powder_meas_cfg: Any

	def _correlate_centroids_to_refs(
		self: Any,
		centroids: 'np.ndarray',
		cluster_radii: 'np.ndarray',
		ref_phases_df: 'pd.DataFrame'
	) -> Tuple[List[float], 'pd.DataFrame']:
		all_ref_phases = ref_phases_df.to_numpy()
		refs_dict = []
		max_raw_confs = []

		for centroid, radius in zip(centroids, cluster_radii):
			distances = np.linalg.norm(all_ref_phases - centroid, axis=1)
			indices = np.where(distances < max(0.1, 5 * radius))[0]
			ref_names = [self.ref_formulae[i] for i in indices]
			ref_phases = np.array([all_ref_phases[i] for i in indices])
			max_raw_conf, refs_dict_row = ReferenceMatchingModule._get_ref_confidences(
				centroid, ref_phases, ref_names
			)
			max_raw_confs.append(max_raw_conf)
			refs_dict.append(refs_dict_row)

		refs_assigned_df = pd.DataFrame(refs_dict)
		return max_raw_confs, refs_assigned_df

	def _assign_reference_phases(self, centroids, rms_dist_cluster):
		min_conf = None
		max_raw_confs = None
		refs_assigned_df = None
		if self.ref_formulae is not None:
			max_raw_confs, refs_assigned_df = ReferenceMatchingModule._correlate_centroids_to_refs(
				self,
				centroids,
				rms_dist_cluster,
				self.ref_phases_df,
			)
			if len(max_raw_confs) > 0:
				max_confs_num = [c for c in max_raw_confs if isinstance(c, (int, float))]
				if len(max_confs_num) > 0:
					min_conf = min(max_confs_num)
		return min_conf, max_raw_confs, refs_assigned_df

	@staticmethod
	def _get_ref_confidences(
		centroid: 'np.ndarray',
		ref_phases: 'np.ndarray',
		ref_names: List[str]
	) -> Tuple[Optional[float], Dict]:
		"""Compute raw and weighted confidence scores for candidate references."""
		if len(ref_phases) == 0:
			refs_dict = {
				f'{cnst.CND_DF_KEY}1': np.nan,
				f'{cnst.CS_RAW_CND_DF_KEY}1': np.nan,
				f'{cnst.CS_CND_DF_KEY}1': np.nan,
			}
			max_raw_conf = None
		else:
			distances = np.linalg.norm(ref_phases - centroid, axis=1)
			raw_confidences = np.exp(-distances**2 / (2 * 0.03**2))

			weights_conf = np.exp(-(1 - raw_confidences)**2 / (2 * 0.3**2))
			weights_conf /= np.sum(weights_conf)
			confidences = raw_confidences * weights_conf

			max_raw_conf = float(np.max(raw_confidences))

			sorted_indices = np.argsort(-confidences)
			sorted_ref_names = np.array(ref_names)[sorted_indices]
			sorted_confidences = confidences[sorted_indices]
			sorted_raw_confs = raw_confidences[sorted_indices]

			refs_dict = {}
			for i, (ref_name, conf, conf_raw) in enumerate(
				zip(sorted_ref_names, sorted_confidences, sorted_raw_confs)
			):
				if conf_raw > 0.05:
					refs_dict[f'{cnst.CND_DF_KEY}{i+1}'] = ref_name
					refs_dict[f'{cnst.CS_CND_DF_KEY}{i+1}'] = np.round(conf, 2)
					refs_dict[f'{cnst.CS_RAW_CND_DF_KEY}{i+1}'] = np.round(conf_raw, 2)

		return max_raw_conf, refs_dict

	# =========================================================================
	# Mixture analysis / NMF
	# =========================================================================

	def _assign_mixtures(self, k, labels, compositions_df, rms_dist_cluster, max_raw_confs, n_points_per_cluster):
		"""
		Determine if clusters are mixtures or single phases, using candidate phases and NMF if needed.

		Returns
		-------
		clusters_assigned_mixtures : list
			List of mixture assignments for each cluster.
			
		Potential improvements
		----------------------
		Instead of using the cluster standard deviation, use covariance of elemental fractions
		to discern clusters that may originate from binary phase mixtures or solid solutions.
		"""
		clusters_assigned_mixtures = []
		ref_formulae = self.ref_formulae or []
		mixture_cfg = self.clustering_cfg.mixture
		for i in range(k):
			# Get compositions of data points included in cluster as np.array (only detectable elements)
			cluster_data = compositions_df[self.detectable_els_sample].iloc[labels == i].values
			max_mix_conf = 0
			mixtures_dicts = []

			if rms_dist_cluster[i] < mixture_cfg.single_phase_max_rms_dist:
				if max_raw_confs is None or len(max_raw_confs) < 1:
					is_cluster_single_phase = n_points_per_cluster[i] > 3
				elif max_raw_confs[i] is not None and max_raw_confs[i] > mixture_cfg.single_phase_min_ref_conf:
					is_cluster_single_phase = True
				else:
					is_cluster_single_phase = False
			else:
				is_cluster_single_phase = False

			if is_cluster_single_phase and not _is_known_phase_mixture(self):
				# Cluster determined to stem from a single phase
				pass
			elif len(ref_formulae) > 1:
				max_mix_raw_conf, mixtures_dicts = ReferenceMatchingModule._identify_mixture_from_refs(self, cluster_data, cluster_ID=i)
				max_mix_conf = max(max_mix_conf, max_mix_raw_conf)
			if not is_cluster_single_phase and max_mix_conf < mixture_cfg.nmf_min_mixture_conf:
				mix_nmf_conf, mixture_dict = ReferenceMatchingModule._identify_mixture_nmf(self, cluster_data, cluster_ID=i)
				if mixture_dict is not None:
					mixtures_dicts.append(mixture_dict)
				max_mix_conf = max(max_mix_conf, mix_nmf_conf)
			mixtures_dicts = ReferenceMatchingModule._rank_mixtures(mixtures_dicts, mixture_cfg.equivalent_recon_error_tol)
			if _draw_intermixing_violin_plots(self):
				reported, _ = ReferenceMatchingModule._reported_mixtures(self, mixtures_dicts)
				for mixture in reported:
					if _POINT_MOLAR_FRS_KEY not in mixture:
						continue
					plot_violin = (
						PlottingModule._save_violin_plot_powder_mixture if len(mixture[cnst.REF_NAME_KEY]) == 2
						else PlottingModule._save_violin_plot_powder_mixture_multi
					)
					plot_violin(cast(Any, self), mixture[_POINT_MOLAR_FRS_KEY], mixture[cnst.REF_NAME_KEY], i)
			for mixture in mixtures_dicts:
				mixture.pop(_POINT_MOLAR_FRS_KEY, None)
			clusters_assigned_mixtures.append(mixtures_dicts)
		return clusters_assigned_mixtures


	@staticmethod
	def _rank_mixtures(mixtures: List[Dict], equivalent_recon_error_tol: float) -> List[Dict]:
		"""
		Rank the mixtures of a cluster, storing the rank in each mixture dictionary.

		Mixtures are ranked by confidence score (i.e., reconstruction error). Mixtures whose
		reconstruction error is within ``equivalent_recon_error_tol`` of the best remaining one
		explain the cluster equally well (e.g., collinear candidate phases), and are ranked by
		number of phases and then by the spread of their phases in composition space, so that
		the phases closest to the cluster are preferred.
		"""
		remaining = sorted(mixtures, key=lambda m: -m[cnst.CONF_SCORE_KEY])
		ranked: List[Dict] = []
		while remaining:
			best_recon_er = remaining[0].get(cnst.MIX_RECON_ERROR_KEY)
			if best_recon_er is None:
				group = remaining[:1]
			else:
				group = [
					m for m in remaining
					if m.get(cnst.MIX_RECON_ERROR_KEY) is not None
					and m[cnst.MIX_RECON_ERROR_KEY] <= best_recon_er + equivalent_recon_error_tol
				]
				group.sort(key=lambda m: (
					len(m[cnst.REF_NAME_KEY]),
					m.get(cnst.MIX_PHASES_SPREAD_KEY, np.inf),
					m[cnst.MIX_RECON_ERROR_KEY],
				))
			ranked.extend(group)
			remaining = [m for m in remaining if not any(m is g for g in group)]

		for rank, mixture in enumerate(ranked, start=1):
			mixture[cnst.MIX_RANK_KEY] = rank
		return ranked


	@staticmethod
	def _in_affine_span(points: 'np.ndarray', basis: 'np.ndarray', tol: float) -> bool:
		"""Whether all ``points`` lie within ``tol`` of the affine span (line, plane, ...) of the rows of ``basis``."""
		directions = (basis[1:] - basis[0]).T
		for p in points:
			offset = p - basis[0]
			if directions.size:
				coeffs = np.linalg.lstsq(directions, offset, rcond=None)[0]
				offset = offset - directions @ coeffs
			if np.linalg.norm(offset) > tol:
				return False
		return True


	@staticmethod
	def _collapse_equivalent_mixtures(
		sorted_mixtures: List[Dict],
		phase_compositions: Dict[str, 'np.ndarray'],
		tol: float,
	) -> Tuple[List[Dict], int]:
		"""
		Drop mixtures that are equivalent decompositions of a better-ranked one, from mixtures sorted by rank.

		Two mixtures are equivalent when their phases span the same mixing line or plane, and that span is
		narrower than the whole composition space (e.g., phases on the SrO-Ta2O5 line: any pair bracketing
		the cluster describes the same mixture). Mixtures spanning the whole composition space (e.g., three
		phases with three elements) are distinct decompositions and are kept. Mixtures with phases of unknown
		composition (free NMF) are kept.

		Returns the kept mixtures and the number of dropped equivalents.
		"""
		kept: List[Dict] = []
		kept_spans: List['np.ndarray'] = []
		n_equivalent = 0
		for mixture in sorted_mixtures:
			names = mixture.get(cnst.REF_NAME_KEY) or []
			if not names or not all(f in phase_compositions for f in names):
				kept.append(mixture)
				continue
			H = np.array([phase_compositions[f] for f in names], dtype=float)
			n_dims = np.linalg.matrix_rank(H[1:] - H[0], tol=tol) if len(H) > 1 else 0
			is_proper_span = n_dims < H.shape[1] - 1
			is_equivalent = is_proper_span and any(
				np.linalg.matrix_rank(K[1:] - K[0], tol=tol) == n_dims and ReferenceMatchingModule._in_affine_span(H, K, tol)
				for K in kept_spans
			)
			if is_equivalent:
				n_equivalent += 1
				continue
			kept.append(mixture)
			kept_spans.append(H)
		return kept, n_equivalent


	def _reported_mixtures(self, mixtures: List[Dict]) -> Tuple[List[Dict], Dict[str, int]]:
		"""
		Mixtures of a cluster reported in Clusters.csv (and plotted), in rank order.

		Equivalent decompositions of a better-ranked mixture are dropped (see ``_collapse_equivalent_mixtures``),
		then ``_select_reported_mixtures`` is applied. Also returns the number of mixtures not reported, by reason:
		'equivalent', 'below_min_conf' (confidence below min_reported_conf_ratio of the best) and 'beyond_max'
		(ranked after the first max_reported_mixtures and not within report_within_conf_ratio of the best).
		"""
		mixture_cfg = self.clustering_cfg.mixture
		sorted_mixtures = ReferenceMatchingModule._sort_mixtures(mixtures)
		candidates, n_equivalent = sorted_mixtures, 0
		phase_compositions = ReferenceMatchingModule._phase_compositions(self)
		if mixture_cfg.collapse_equivalent_mixtures and phase_compositions:
			candidates, n_equivalent = ReferenceMatchingModule._collapse_equivalent_mixtures(
				sorted_mixtures, phase_compositions, mixture_cfg.equivalent_span_tol
			)
		reported = ReferenceMatchingModule._select_reported_mixtures(candidates, mixture_cfg)
		best_conf = max((float(m.get(cnst.CONF_SCORE_KEY, 0.0)) for m in sorted_mixtures), default=0.0)
		n_below = sum(
			1 for m in candidates
			if float(m.get(cnst.CONF_SCORE_KEY, 0.0)) < mixture_cfg.min_reported_conf_ratio * best_conf
		)
		counts = {
			'equivalent': n_equivalent,
			'below_min_conf': n_below,
			'beyond_max': len(candidates) - len(reported) - n_below,
			'best_conf': best_conf,
		}
		return reported, counts


	def _phase_compositions(self) -> Dict[str, 'np.ndarray']:
		"""Compositions of the candidate phases (detectable elements), by formula."""
		els = list(getattr(self, 'detectable_els_sample', None) or [])
		ref_formulae = list(getattr(self, 'ref_formulae', None) or [])
		if not els or not ref_formulae or getattr(self, 'ref_phases_df', None) is None:
			return {}
		return {f: self.ref_phases_df[els].iloc[i].to_numpy(dtype=float) for i, f in enumerate(ref_formulae)}


	@staticmethod
	def _select_reported_mixtures(sorted_mixtures: List[Dict], mixture_cfg: Any) -> List[Dict]:
		"""
		Select the mixtures reported in Clusters.csv, from mixtures sorted by rank.

		A mixture is reported if its confidence is at least ``min_reported_conf_ratio`` times the best
		confidence of the cluster, and it is among the first ``max_reported_mixtures`` or its confidence
		is at least ``report_within_conf_ratio`` times the best one.
		"""
		if not sorted_mixtures:
			return []
		best_conf = max(float(m.get(cnst.CONF_SCORE_KEY, 0.0)) for m in sorted_mixtures)
		reported = []
		for i, mixture in enumerate(sorted_mixtures):
			conf = float(mixture.get(cnst.CONF_SCORE_KEY, 0.0))
			if conf < mixture_cfg.min_reported_conf_ratio * best_conf:
				continue
			if i < mixture_cfg.max_reported_mixtures or conf >= mixture_cfg.report_within_conf_ratio * best_conf:
				reported.append(mixture)
		return reported


	@staticmethod
	def _mixtures_not_reported_note(counts: Dict[str, float], mixture_cfg: Any) -> Optional[str]:
		"""Clusters.csv note on the mixtures saved in the ledger but not reported, with the reason and threshold."""
		reasons = []
		if counts['equivalent']:
			reasons.append(f"{counts['equivalent']} equivalent to a listed mixture (same mixing line/plane)")
		if counts['below_min_conf']:
			floor = mixture_cfg.min_reported_conf_ratio
			reasons.append(
				f"{counts['below_min_conf']} with confidence below {floor * 100:.0f}% of the best "
				f"(CS_mix < {floor * counts['best_conf']:.2f})"
			)
		if counts['beyond_max']:
			within = mixture_cfg.report_within_conf_ratio
			reasons.append(
				f"{counts['beyond_max']} ranked after the first {mixture_cfg.max_reported_mixtures} with confidence "
				f"more than {(1 - within) * 100:.0f}% below the best (CS_mix < {within * counts['best_conf']:.2f})"
			)
		if not reasons:
			return None
		n = int(counts['equivalent'] + counts['below_min_conf'] + counts['beyond_max'])
		return (f"{n} more mixture(s) saved in ledger.json (clusters_assigned_mixtures of this clustering analysis): "
				+ '; '.join(reasons))


	@staticmethod
	def _sort_mixtures(mixtures: List[Dict]) -> List[Dict]:
		"""Sort the mixtures of a cluster by their stored rank, or by confidence score for results saved without ranks."""
		if mixtures and all(cnst.MIX_RANK_KEY in m for m in mixtures):
			return sorted(mixtures, key=lambda m: m[cnst.MIX_RANK_KEY])
		return sorted(mixtures, key=lambda m: -float(m.get(cnst.CONF_SCORE_KEY, 0.0)))


	def _identify_mixture_from_refs(self, X: 'np.ndarray', cluster_ID: Optional[int] = None) -> Tuple[float, List[Dict]]:
		"""
		Identify mixtures within a cluster by testing combinations of candidate phases using constrained optimization.

		For each combination of candidate phases, tests if the cluster compositions (X)
		can be well described by a linear combination of those phases, using
		non-negative matrix factorization (NMF) with fixed bases.

		All pairs of candidate phases are tested. Since adding phases always lowers the
		reconstruction error, combinations of more phases (up to ``clustering_cfg.mixture.max_n_phases``)
		are only tested when no combination with fewer phases reconstructs the cluster with an error
		below ``clustering_cfg.mixture.max_recon_error``, stopping at the first number of phases that does.
		This way the number of phases follows the shape of the cluster: a line for binary mixtures,
		a plane for ternary ones.

		Parameters
		----------
		X : np.ndarray
			Cluster data (compositions), shape (n_samples, n_features).
		cluster_ID : int
			Current cluster ID. Used for violin plot name

		Returns
		-------
		max_confidence : float
			The highest confidence score among all tested mixtures.
		mixtures_dicts : list of Dict
			List of mixture descriptions for all successful combinations of candidate phases.
		"""
		ref_formulae = self.ref_formulae or []
		n_refs = len(self.ref_phases_df)
		mixture_cfg = self.clustering_cfg.mixture

		mixtures_dicts = []
		max_confidence = 0

		for n_phases in range(2, min(mixture_cfg.max_n_phases, n_refs) + 1):
			# Binary mixtures are listed with a looser threshold; mixtures of more phases
			# only if they reconstruct the cluster within noise
			max_listed_recon_error = mixture_cfg.max_recon_error_binary if n_phases == 2 else mixture_cfg.max_recon_error
			min_recon_er = np.inf

			for ref_comb in itertools.combinations(range(n_refs), n_phases):
				# Get the names of the candidate phases in this combination
				ref_names = [ref_formulae[ref_i] for ref_i in ref_comb]

				# Weights of references, for molar concentrations of parent phases
				ref_weights = [self.ref_weights_in_mixture[ref_i] for ref_i in ref_comb]

				# Get matrix of basis vectors (H) for the candidate phases
				H = np.array([
					self.ref_phases_df[self.detectable_els_sample].iloc[ref_i].values
					for ref_i in ref_comb
				])

				# Perform NMF with fixed H to fit the cluster data as a mixture of the candidate phases
				W, _ = ReferenceMatchingModule._nmf_with_constraints(self, X, n_components=n_phases, fixed_H=H)

				# Compute reconstruction error for the fit
				recon_er = ReferenceMatchingModule._calc_reconstruction_error(self, X, W, H)
				min_recon_er = min(min_recon_er, recon_er)

				# If the combination yields an acceptable reconstruction error, store the result
				mix_dict, conf = ReferenceMatchingModule._get_mixture_dict_with_conf(
					self, W, ref_weights, recon_er, ref_names, cluster_ID, max_listed_recon_error, H
				)
				if mix_dict is not None:
					mixtures_dicts.append(mix_dict)
					max_confidence = max(max_confidence, conf)

			if min_recon_er < mixture_cfg.max_recon_error:
				# Fewest phases that explain the cluster within noise
				break

		return max_confidence, mixtures_dicts


	def _calc_reconstruction_error(
		self,
		X: 'np.ndarray',
		W: 'np.ndarray',
		H: 'np.ndarray'
	) -> float:
		"""
		Calculate the reconstruction error for a matrix factorization X ≈ W @ H.

		The error metric is an exponential penalty (with parameter ``clustering_cfg.mixture.recon_error_alpha``)
		applied to the absolute difference between X and its reconstruction W @ H, normalized by the
		number of elements in X. This penalizes large deviations more strongly.
		"""
		WH = np.dot(W, H)
		alpha = self.clustering_cfg.mixture.recon_error_alpha
		norm = np.sum(np.exp(alpha * np.abs(X - WH)) - 1)
		m, n = X.shape
		normalized_norm = norm / (m * n)
		return normalized_norm


	def _get_mixture_dict_with_conf(
		self,
		W: 'np.ndarray',
		ref_weights: List[float],
		reconstruction_error: float,
		ref_names: List[str],
		cluster_ID: Optional[int] = None,
		max_recon_error: Optional[float] = None,
		H: Optional['np.ndarray'] = None
	) -> Tuple[Optional[Dict], float]:
		"""
		Evaluate if a cluster is a mixture of candidate phases, and compute a confidence score.

		If the reconstruction error is below ``max_recon_error`` (``clustering_cfg.mixture.max_recon_error_binary``
		if None), computes a confidence score and transforms the NMF coefficients into molar fractions.
		Returns a dictionary describing the mixture and the confidence score. If the phase compositions
		``H`` are given, their spread in composition space (sum of pairwise distances) is also stored,
		used to rank mixtures that reconstruct the cluster equally well.
		"""
		mixture_cfg = self.clustering_cfg.mixture
		if max_recon_error is None:
			max_recon_error = mixture_cfg.max_recon_error_binary

		# For known powder mixtures, all binary mixtures are kept (violin plots are drawn in
		# _assign_mixtures, for the mixtures reported in Clusters.csv)
		keep_binary = len(ref_names) == 2 and _is_known_phase_mixture(self)

		if reconstruction_error < max_recon_error or keep_binary:
			conf = np.exp(-reconstruction_error**2 / (2 * mixture_cfg.conf_sigma**2))

			# NMF coefficients are fractions of atoms (or mass) contributed by each phase. Dividing by
			# the atoms (or mass) per formula unit of each phase gives its molar fraction.
			W_mol_frs = np.asarray(W, dtype=float) / np.asarray(ref_weights, dtype=float)
			W_mol_frs /= W_mol_frs.sum(axis=1, keepdims=True)

			mol_frs_norm_means = np.mean(W_mol_frs, axis=0)
			mol_frs_norm_stddevs = np.std(W_mol_frs, axis=0)

			mixture_dict = {
				cnst.REF_NAME_KEY: ref_names,
				cnst.CONF_SCORE_KEY: conf,
				cnst.MOLAR_FR_MEAN_KEY: mol_frs_norm_means[0],
				cnst.MOLAR_FR_STDEV_KEY: mol_frs_norm_stddevs[0],
				cnst.MOLAR_FRS_MEAN_KEY: [float(v) for v in mol_frs_norm_means],
				cnst.MOLAR_FRS_STDEV_KEY: [float(v) for v in mol_frs_norm_stddevs],
				cnst.MIX_RECON_ERROR_KEY: float(reconstruction_error),
				# Per-point molar fractions, for violin plots; removed in _assign_mixtures (not saved in the ledger)
				_POINT_MOLAR_FRS_KEY: W_mol_frs,
			}
			if H is not None:
				mixture_dict[cnst.MIX_PHASES_SPREAD_KEY] = float(sum(
					np.linalg.norm(H[a] - H[b]) for a, b in itertools.combinations(range(len(H)), 2)
				))
		else:
			mixture_dict = None
			conf = 0

		return mixture_dict, conf


	def _nmf_with_constraints(
		self,
		X: 'np.ndarray',
		n_components: int,
		fixed_H: Optional['np.ndarray'] = None
	) -> Tuple['np.ndarray', 'np.ndarray']:
		"""
		Perform Non-negative Matrix Factorization (NMF) with optional constraints on the factor matrices.

		This function alternates between optimizing two non-negative matrices W and H, such that X ≈ W @ H:
		  - If H is fixed (provided via fixed_H), only W is updated.
		  - If H is not fixed, both W and H are updated via alternating minimization.

		Constraints:
		  - Both W and H are non-negative.
		  - The rows of both W (sum of coefficients) and H (sum of elemental fractions) sum to 1.
		  - Sparsity regularization (L1) is applied to H when it is updated.
		"""
		max_iter = 1000
		convergence_tol = 1e-3
		lambda_H = 0

		W = np.random.rand(X.shape[0], n_components)
		if fixed_H is None:
			H = np.random.rand(n_components, X.shape[1])
		else:
			H = fixed_H

		prev_W: Optional[np.ndarray] = None
		prev_H: Optional[np.ndarray] = None
		convergence = np.inf
		i = 0

		while convergence > convergence_tol and i < max_iter:
			W_var = cp.Variable((X.shape[0], n_components), nonneg=True)
			objective_W = cp.Minimize(cp.sum_squares(X - W_var @ H))
			constraints_W: List[Constraint] = [cast(Constraint, cp.sum(W_var, axis=1) == 1)]
			problem_W = cp.Problem(objective_W, constraints_W)
			problem_W.solve(solver=cp.ECOS)
			if W_var.value is None:
				raise RuntimeError('Constrained NMF failed to produce W values.')
			W = W_var.value

			if fixed_H is None:
				H_var = cp.Variable((n_components, X.shape[1]), nonneg=True)
				objective_H = cp.Minimize(
					cp.sum_squares(X - W @ H_var) + lambda_H * cp.norm1(H_var)
				)
				constraints_H: List[Constraint] = [cast(Constraint, cp.sum(H_var, axis=1) == 1)]
				problem_H = cp.Problem(objective_H, constraints_H)
				problem_H.solve(solver=cp.ECOS)
				if H_var.value is None:
					raise RuntimeError('Constrained NMF failed to produce H values.')
				H = H_var.value

			if prev_W is None:
				convergence_W = np.inf
			else:
				convergence_W = np.linalg.norm(W - prev_W, 'fro')

			if fixed_H is None and prev_H is not None:
				convergence_H = np.linalg.norm(H - prev_H, 'fro')
			else:
				convergence_H = 0.0
			convergence = max(convergence_W, convergence_H)

			prev_W, prev_H = W, H
			i += 1

		return W, H


	def _identify_mixture_nmf(
		self,
		X: 'np.ndarray',
		n_components: int = 2,
		cluster_ID: Optional[int] = None
	) -> Tuple[float, Optional[Dict]]:
		"""
		Identify a mixture within a cluster using unconstrained NMF (Non-negative Matrix Factorization).

		This method fits the cluster data X to n_components using NMF with constraints (rows of W and H sum to 1),
		evaluates the reconstruction error, and if acceptable, returns a dictionary describing the mixture and a confidence score.
		"""
		mixture_dict = None
		conf = 0

		W, H = ReferenceMatchingModule._nmf_with_constraints(self, X, n_components)
		recon_er = ReferenceMatchingModule._calc_reconstruction_error(self, X, W, H)
		ref_names, ref_weights = ReferenceMatchingModule._get_pretty_formulas_nmf(self, H, n_components)
		mixture_dict, conf = ReferenceMatchingModule._get_mixture_dict_with_conf(self, W, ref_weights, recon_er, ref_names, cluster_ID, H=H)

		return conf, mixture_dict


	def _get_pretty_formulas_nmf(
		self,
		phases: 'np.ndarray',
		n_components: int
	) -> Tuple[List[str], List[float]]:
		"""
		Generate human-readable (pretty) formulas from NMF bases, accounting for data noise.

		For each component, filters out small fractions, constructs a composition dictionary,
		and returns a formula string and a weight or atom count, depending on the clustering feature.
		"""
		ref_names = []
		ref_weights = []

		for i in range(n_components):
			frs = phases[i, :].copy()
			frs[frs < 0.005] = 0

			fr_dict = {el: float(fr) for el, fr in zip(self.detectable_els_sample, frs)}

			if self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
				comp = Composition().from_weight_dict(cast(Dict[Any, float], fr_dict))
			elif self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
				comp = Composition(fr_dict)
			else:
				raise ValueError(f'Unsupported clustering feature type: {self.clustering_cfg.features}')

			formula = comp.get_integer_formula_and_factor()[0]
			ref_integer_comp = Composition(formula)
			min_at_n = min(ref_integer_comp.get_el_amt_dict().values())
			pretty_at_frs = {el: round(n / min_at_n, 1) for el, n in ref_integer_comp.get_el_amt_dict().items()}
			pretty_comp = Composition(pretty_at_frs)
			pretty_formula = pretty_comp.formula
			ref_names.append(pretty_formula)

			if self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
				ref_weights.append(pretty_comp.weight)
			elif self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
				n_atoms_in_formula = sum(pretty_comp.get_el_amt_dict().values())
				ref_weights.append(n_atoms_in_formula)

		return ref_names, ref_weights


	def _build_mixtures_df(
		self,
		clusters_assigned_mixtures: List[List[Dict]]
	) -> 'pd.DataFrame':
		"""
		Build a DataFrame summarizing mixture assignments for each cluster.

		For each cluster, sorts mixture dictionaries by rank (see ``_rank_mixtures``) and extracts:
		  - candidate phase names (as a comma-separated string)
		  - Confidence score
		  - Molar ratio (mean / (1 - mean)), for binary mixtures only
		  - Mean and standard deviation of the main component's molar fraction
		  - Mean molar fractions of all phases, for mixtures of more than two phases

		Equivalent decompositions of a better-ranked mixture are dropped (see ``_collapse_equivalent_mixtures``),
		then only the mixtures selected by ``_select_reported_mixtures`` are included; if others were found,
		a note gives their number and where they are saved.
		"""
		mixture_cfg = self.clustering_cfg.mixture
		mixtures_strings_dict = []
		for mixtures_dict in clusters_assigned_mixtures:
			if mixtures_dict:
				reported_mixtures, counts = ReferenceMatchingModule._reported_mixtures(self, mixtures_dict)
				cluster_mix_dict = {}
				for i, mixture_dict in enumerate(reported_mixtures, start=1):
					is_binary = len(mixture_dict[cnst.REF_NAME_KEY]) == 2
					cluster_mix_dict[f'{cnst.MIX_DF_KEY}{i}'] = ', '.join(mixture_dict[cnst.REF_NAME_KEY])
					cluster_mix_dict[f'{cnst.CS_MIX_DF_KEY}{i}'] = float(f"{mixture_dict[cnst.CONF_SCORE_KEY]:.2f}")
					cluster_mix_dict[f'{cnst.MIX_MOLAR_RATIO_DF_KEY}{i}'] = (
						np.round(mixture_dict[cnst.MOLAR_FR_MEAN_KEY] / (1 - mixture_dict[cnst.MOLAR_FR_MEAN_KEY]), 2)
						if is_binary else np.nan
					)
					cluster_mix_dict[f'{cnst.MIX_FIRST_COMP_MEAN_DF_KEY}{i}'] = np.round(mixture_dict[cnst.MOLAR_FR_MEAN_KEY], 2)
					cluster_mix_dict[f'{cnst.MIX_FIRST_COMP_STDEV_DF_KEY}{i}'] = np.round(mixture_dict[cnst.MOLAR_FR_STDEV_KEY], 2)
					if not is_binary and cnst.MOLAR_FRS_MEAN_KEY in mixture_dict:
						cluster_mix_dict[f'{cnst.MIX_ALL_COMPS_MEAN_DF_KEY}{i}'] = '/'.join(
							f'{x:.2f}' for x in mixture_dict[cnst.MOLAR_FRS_MEAN_KEY]
						)
				note = ReferenceMatchingModule._mixtures_not_reported_note(counts, mixture_cfg)
				if note:
					cluster_mix_dict[cnst.MIX_MORE_DF_KEY] = note
				mixtures_strings_dict.append(cluster_mix_dict)
			else:
				mixtures_strings_dict.append({})

		mixtures_df = pd.DataFrame(mixtures_strings_dict)
		return mixtures_df
