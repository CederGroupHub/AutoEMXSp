#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Clustering mixin for EMXSp composition analysis."""

import os
import warnings
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from yellowbrick.cluster import KElbowVisualizer

import autoemx.utils.constants as cnst
from autoemx.utils.helper import print_single_separator

from autoemx._logging import get_logger
logger = get_logger(__name__)


class ClusteringModule:
    def _find_optimal_k(self, compositions_df, k, compute_k_only_once = False):
        """
        Determine the optimal number of clusters for k-means.
    
        Returns
        -------
        k : int
            Optimal number of clusters.
        """
        if not k:
            # Check if there is only one single cluster, or no clusters
            is_single_cluster = ClusteringModule._is_single_cluster(compositions_df, verbose=self.verbose)
            if is_single_cluster or self.clustering_cfg.max_k <= 1:
                k = 1
            elif compute_k_only_once:
                # Get number of clusters (k) and optionally save the plot
                results_dir = self.analysis_dir if self.plot_cfg.save_plots else None
                k = ClusteringModule._get_k(
                    compositions_df, self.clustering_cfg.max_k, self.clustering_cfg.k_finding_method,
                    show_plot=self.plot_cfg.show_plots, results_dir=results_dir
                )
            else:
                # Calculate most frequent number of clusters (k) with elbow method. Does not save the plot
                k = ClusteringModule._get_most_freq_k(
                    compositions_df, self.clustering_cfg.max_k, self.clustering_cfg.k_finding_method,
                    verbose=self.verbose
                )
        elif self.verbose:
            print_single_separator()
            logger.info(f"ℹ️ Number of clusters was forced to be {k}")
        return k
    
    
    @staticmethod
    def _get_most_freq_k(
        compositions_df: 'pd.DataFrame',
        max_k: int,
        k_finding_method: str,
        verbose: bool = False,
        show_plot: bool = False,
        results_dir: str = None
    ) -> int:
        """
        Determine the most frequent optimal number of clusters (k) for the given compositions.
    
        This method repeatedly runs the k-finding algorithm and selects the most robust k value.
        It loops until it finds a value of k that is at least twice as frequent as the second most frequent value,
        or until a maximum number of iterations is reached.
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame containing the compositions to cluster.
        max_k : int
            Maximum number of clusters to test.
        k_finding_method : str
            Method used to determine the optimal k (passed to _get_k).
        verbose : bool, optional
            If True, print progress and summary information.
        show_plot : bool, optional
            If True, display plots for each k-finding run.
        results_dir : str, optional
            Directory to save plots/results (if applicable).
    
        Returns
        -------
        k : int
            The most robustly determined number of clusters.
    
        Raises
        ------
        ValueError
            If there are not enough data points to determine clusters.
    
        Notes
        -----
        - The function tries up to 5 times to find a dominant k value.
        - If a tie or ambiguity remains, it picks the smallest k with frequency ≥ half of the most frequent.
    
        """
        if len(compositions_df) < 2:
            raise ValueError("Not enough data points to determine clusters (need at least 2).")
    
        if verbose:
            print_single_separator()
            logger.info("📊 Computing number of clusters k...")
    
        k_found = []
        k = None
        max_iter = 5
        i = 0
    
        max_allowed_k = len(compositions_df) - 1  # KElbowVisualizer throws error if max_k > n_samples - 1
        if max_k > max_allowed_k:
            max_k = max_allowed_k
            warnings.warn(
                f"Maximum number of clusters reduced to {max_allowed_k} because number of clustered points is {max_allowed_k + 1}.",
                UserWarning
            )
    
        while k is None:
            i += 1
            for n in range(20):
                k_val = ClusteringModule._get_k(
                    compositions_df, max_k, k_finding_method,
                    show_plot=show_plot, results_dir=results_dir
                )
                if not isinstance(k_val, int) or k_val < 1:
                    continue  # skip invalid k values
                k_found.append(k_val)
    
            if not k_found:
                raise ValueError("No valid cluster counts were found.")
    
            counts = np.bincount(k_found)
            total = counts.sum()
            sorted_k = np.argsort(-counts)  # descending
            first_k = sorted_k[0]
            first_count = counts[first_k]
    
            if len(sorted_k) > 1:
                second_k = sorted_k[1]
                second_count = counts[second_k]
            else:
                second_k = None
                second_count = 0
    
            # Check if first is at least twice as common as second
            if second_count == 0 or first_count >= 2 * second_count:
                k = first_k
            elif i >= max_iter:  # max_iter reached
                # Pick smallest of all k values whose frequency is ≥ half of the most frequent
                threshold = first_count / 2
                k = min([k_val for k_val, count in enumerate(counts) if count >= threshold])
            else:
                k = None
    
        if verbose:
            logger.info(f"📊 Most frequent k: {first_k} (count = {first_count}, frequency = {first_count / total:.2%})")
            if second_k is not None and second_k != 0:
                logger.info(f"📊 Second most frequent k: {second_k} (count = {second_count}, frequency = {second_count / total:.2%})")
            if len(np.where(counts == first_count)[0]) > 1:
                logger.info(f"ℹ️ Tie detected among: {np.where(counts == first_count)[0].tolist()} (choosing {k})")
    
        return int(k)
    
    
    @staticmethod
    def _get_k(
        compositions_df: 'pd.DataFrame',
        max_k: int = 6,
        method: str = 'silhouette',
        model: 'KMeans' = None,
        results_dir: str = None,
        show_plot: bool = False
    ) -> int:
        """
        Determine the optimal number of clusters for the data using visualizer methods.
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame containing the compositions to cluster.
        max_k : int, optional
            Maximum number of clusters to test (default: 6).
        method : str, optional
            Method for evaluating the number of clusters. One of 'elbow', 'silhouette', or 'calinski_harabasz' (default: 'silhouette').
        model : KMeans or compatible, optional
            Clustering model to use (default: KMeans(n_init='auto')).
        results_dir : str, optional
            Directory to save the plot (if provided).
        show_plot : bool, optional
            If True, show the plot interactively (default: False).
    
        Returns
        -------
        optimal_k : int
            The optimal number of clusters found.
    
        Raises
        ------
        ValueError
            If an unsupported method is provided.
    
        Notes
        -----
        - Uses yellowbrick's KElbowVisualizer.
        - For 'elbow', finds the inflection point; for 'silhouette' or 'calinski_harabasz', finds the k with the highest score.
        - If cluster finding fails, returns 1 and prints a warning.
        - If `show_plot` is True, the plot is shown interactively.
        - If `results_dir` is provided, the plot is saved as 'Elbow_plot.png'.
        """
        if model is None:
            from sklearn.cluster import KMeans
            model = KMeans(n_init='auto')
    
        # Map 'elbow' to 'distortion' for yellowbrick, but keep original for logic
        user_method = method
        if method == 'elbow':
            yb_method = 'distortion'
        elif method in ['silhouette', 'calinski_harabasz']:
            yb_method = method
        else:
            raise ValueError(f"Unsupported method '{method}' for evaluating number of clusters.")
    
        plt.figure(figsize=(10, 8))
        visualizer = KElbowVisualizer(model, k=max_k, metric=yb_method, timings=True, show=False)
    
        try:
            visualizer.fit(compositions_df)
        except ValueError as er:
            warnings.warn(f"Number of clusters could not be identified due to the following error:\n{er}\nForcing k = 1.", UserWarning)
            return 1
    
        # Get optimal number of clusters
        if user_method == 'elbow':
            optimal_k = visualizer.elbow_value_
        elif user_method in ['silhouette', 'calinski_harabasz']:
            # For silhouette and calinski_harabasz, k_scores_ is indexed from k=2
            optimal_k = np.argmax(visualizer.k_scores_) + 2
            visualizer.elbow_value_ = optimal_k  # For correct plotting
    
        # Add labels
        ax1, ax2 = visualizer.axes
        ax1.set_ylabel(f'{user_method} score')
        ax1.set_xlabel('k: number of clusters')
        ax2.set_ylabel('Fit time (sec)')
        ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    
        if show_plot:
            plt.ion()
            visualizer.show()
            plt.pause(0.001)
    
        if results_dir:
            fig = visualizer.fig
            fig.savefig(os.path.join(results_dir, 'Elbow_plot.png'))
    
        if not show_plot:
            plt.close(visualizer.fig)
    
        # Set k to 1 if elbow method was unsuccessful
        if optimal_k is None:
            optimal_k = 1
    
        return int(optimal_k)
    
    
    @staticmethod
    def _is_single_cluster(
        compositions_df: 'pd.DataFrame',
        verbose: bool = False
    ) -> bool:
        """
        Determine if the data effectively forms a single cluster using k-means and silhouette analysis.
    
        This method:
          - Fits k-means with k=1 and calculates the RMS distance from the centroid.
          - Fits k-means with k=2 multiple times, keeping the best silhouette score and inertia.
          - Uses empirically determined thresholds on silhouette score, centroid distance, and inertia ratio
            to decide if the data forms a single cluster or multiple clusters.
    
        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to analyze.
        verbose : bool, optional
            If True, print detailed output of the clustering metrics.
    
        Returns
        -------
        is_single_cluster : bool
            True if the data is best described as a single cluster, False otherwise.
    
        Notes
        -----
        - Uses silhouette score and inertia ratio as main criteria.
        - Empirical thresholds: mean centroid distance < 0.025, silhouette < 0.5, or inertia ratio < 1.5.
        """
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        import numpy as np
    
        if verbose:
            print_single_separator()
            logger.info('📊 Checking if more than 1 cluster is present...')
    
        # Fit k-means for k=1
        kmeans_1 = KMeans(n_clusters=1, random_state=0, n_init='auto')
        kmeans_1.fit(compositions_df)
        inertia_1 = kmeans_1.inertia_
        rms_distance_1 = np.sqrt(inertia_1 / len(compositions_df))
    
        # Fit k-means for k=2, keep best silhouette score and inertia
        best_silhouette_score_2 = -1
        best_inertia_2 = None
        for _ in range(10):
            kmeans_2 = KMeans(n_clusters=2, random_state=None, n_init='auto')
            labels_2 = kmeans_2.fit_predict(compositions_df)
            inertia_2 = kmeans_2.inertia_
            sil_score = silhouette_score(compositions_df, labels_2)
            if sil_score > best_silhouette_score_2:
                best_silhouette_score_2 = sil_score
                best_inertia_2 = inertia_2
    
        ratio_inertias = inertia_1 / best_inertia_2 if best_inertia_2 else float('inf')
    
        if verbose:
            logger.info(f"📊 RMS distance for k=1: {rms_distance_1*100:.1f}%")
            logger.info(f"📊 Inertia for k=1: {inertia_1:.3f}")
            logger.info(f"📊 Inertia for k=2: {best_inertia_2:.3f}")
            logger.info(f"📊 Ratio of inertia for k=1 over k=2: {ratio_inertias:.2f}")
            logger.info(f"📊 Silhouette Score for k=2: {best_silhouette_score_2:.2f}")
    
        # Empirical decision logic
        if rms_distance_1 < 0.03:
            is_single_cluster = True
            reason_str = 'd_rms < 3%'
        elif best_silhouette_score_2 < 0.5:
            is_single_cluster = True
            reason_str = 's < 0.5'
        elif best_silhouette_score_2 > 0.6:
            is_single_cluster = False
            reason_str = 's > 0.6'
        elif ratio_inertias < 1.5:
            is_single_cluster = True
            reason_str = 'ratio of inertias < 1.5'
        else:
            is_single_cluster = False
            reason_str = 'ratio of inertias > 1.5 and 0.5 < s < 0.6'
    
        if verbose:
            if is_single_cluster:
                logger.info("ℹ️ " + reason_str + ": The data effectively forms a single cluster.")
            else:
                logger.info("ℹ️ " + reason_str + ": The data forms multiple clusters.") 
    
        return is_single_cluster
    #%% Clustering operations
    # ============================================================================= 
    def _run_kmeans_clustering(self, k, compositions_df):
        """
        Run k-means clustering multiple times and select the best solution by silhouette score.
    
        Returns
        -------
        kmeans : KMeans
            The best fitted KMeans instance.
        labels : np.ndarray
            Cluster labels for each composition.
        sil_score : float
            Best silhouette score obtained.
        """
        if k > 1:
            n_clustering_eval = 20  # Number of clustering evaluations to run
            best_sil_score = -np.inf  # Initialise best silhouette score
            for _ in range(n_clustering_eval):
                # K-means is not ideal for clusters with varying sizes/densities. Consider alternatives (e.g., GMM).
                is_clustering_ok = False
                max_loops_nonneg_silh = 50  # Max loops to find clustering solutions with no negative silhouette values
                n_loop = 0
                while not is_clustering_ok and n_loop < max_loops_nonneg_silh:
                    n_loop += 1
                    kmeans, labels = ClusteringModule._get_clustering_kmeans(self, k, compositions_df)
                    silhouette_vals = silhouette_samples(compositions_df, labels)
                    if np.all(silhouette_vals > 0):
                        # Clustering is accepted only if all silhouette values are positive (no wrong clustering)
                        is_clustering_ok = True
                    sil_score = silhouette_score(compositions_df, labels)
                if sil_score > best_sil_score:
                    best_kmeans = kmeans
                    best_labels = labels
                    best_sil_score = sil_score
            return best_kmeans, best_labels, best_sil_score
        else:
            # Clustering with k = 1 is trivial, and has no silhouette score
            kmeans, labels = ClusteringModule._get_clustering_kmeans(self, k, compositions_df)
            return kmeans, labels, np.nan
    

    def _prepare_composition_dataframes(self, compositions_list_at, compositions_list_w):
        """
        Convert lists of compositions to DataFrames for clustering.
    
        Returns
        -------
        compositions_df : pd.DataFrame
            DataFrame of compositions for clustering (feature set selected).
        compositions_df_other_fr : pd.DataFrame
            DataFrame of compositions in the alternate fraction representation.
        """
        # Substitute nan with 0
        if self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
            compositions_df = pd.DataFrame(compositions_list_at).fillna(0)
            compositions_df_other_fr = (pd.DataFrame(compositions_list_w)).fillna(0) 
        elif self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
            compositions_df = pd.DataFrame(compositions_list_w).fillna(0)
            compositions_df_other_fr = (pd.DataFrame(compositions_list_at)).fillna(0)
            
        return compositions_df, compositions_df_other_fr
    
    
    def _get_clustering_kmeans(
        self,
        k: int,
        compositions_df: 'pd.DataFrame'
    ) -> Tuple['KMeans', 'np.ndarray']:
        """
        Perform k-means clustering on the given compositions.
    
        Parameters
        ----------
        k : int
            The number of clusters to find.
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to cluster.
    
        Returns
        -------
        kmeans : KMeans
            The fitted KMeans object.
        labels : np.ndarray
            Array of cluster (phase) labels for each composition point.
    
        Raises
        ------
        ValueError
            If clustering is unsuccessful due to invalid data or parameters.
    
        Notes
        -----
        - Uses k-means++ initialization and scikit-learn's default settings.
        - n_init='auto' requires scikit-learn >= 1.2.0.
        """
        try:
            kmeans = KMeans(n_clusters=k, init='k-means++', n_init='auto')
            # Perform clustering. Returns labels (array of cluster (= phase) ID each composition point belongs to)
            labels = kmeans.fit_predict(compositions_df)
        except Exception as e:
            raise ValueError(f"Clustering unsuccessful due to the following error:\n{e}")
    
        return kmeans, labels

    
    def _get_clustering_dbscan(
        self,
        compositions_df: 'pd.DataFrame'
    ) -> Tuple['np.ndarray', int]:
        """
        Perform DBSCAN clustering on the given compositions.

        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to cluster.

        Returns
        -------
        labels : np.ndarray
            Array of cluster labels for each composition point. Noise points are labeled as -1.
        num_labels : int
            Number of unique clusters found (excluding noise points).

        Raises
        ------
        ValueError
            If clustering is unsuccessful due to invalid data or parameters.

        Notes
        -----
        - eps, min_samples and metric are read from ``self.clustering_cfg.dbscan``.
        - The number of clusters excludes noise points (label -1).
        """
        dbscan_cfg = self.clustering_cfg.dbscan
        try:
            dbscan = DBSCAN(
                eps=dbscan_cfg.eps,
                min_samples=dbscan_cfg.min_samples,
                metric=dbscan_cfg.metric,
            )
            labels = dbscan.fit_predict(compositions_df)
        except Exception as e:
            raise ValueError(f"Clustering unsuccessful due to the following error:\n{e}")

        # Get the number of unique labels, excluding noise (-1)
        num_labels = len(set(labels)) - (1 if -1 in labels else 0)

        return labels, num_labels


    @staticmethod
    def _normalize_dbscan_labels(labels: 'np.ndarray') -> 'np.ndarray':
        """
        Remap DBSCAN cluster labels to a contiguous ``0..k-1`` range.

        DBSCAN may emit arbitrary, non-contiguous integer labels (e.g. ``{-1, 0, 3, 7}``)
        depending on density-reachability ordering. Downstream statistics, reference
        matching and mixture assignment assume contiguous labels indexed like k-means.
        Noise points (``-1``) are preserved as ``-1``.

        Parameters
        ----------
        labels : np.ndarray
            Raw DBSCAN labels.

        Returns
        -------
        np.ndarray
            Labels with non-noise clusters renumbered to ``0..k-1``.
        """
        labels = np.asarray(labels)
        unique_clusters = sorted(int(lbl) for lbl in set(labels.tolist()) if lbl != -1)
        remap = {old: new for new, old in enumerate(unique_clusters)}
        return np.array([remap[int(lbl)] if lbl != -1 else -1 for lbl in labels])


    def _run_dbscan_clustering(
        self,
        compositions_df: 'pd.DataFrame'
    ) -> Tuple['np.ndarray', 'np.ndarray', int, float, float]:
        """
        Run DBSCAN clustering and derive k-means-compatible artifacts.

        DBSCAN does not expose centroids or an inertia value, so they are computed
        here from the resulting clusters to keep downstream statistics, plotting and
        reference matching agnostic to the clustering method.

        Parameters
        ----------
        compositions_df : pd.DataFrame
            DataFrame of samples (rows) and features (columns) to cluster.

        Returns
        -------
        labels : np.ndarray
            Cluster labels normalized to ``0..k-1``; noise points are labeled ``-1``.
        centroids : np.ndarray
            Mean composition of each cluster, shape ``(k, n_features)``.
        k : int
            Number of clusters found (excluding noise).
        sil_score : float
            Silhouette score over non-noise points (``nan`` if fewer than 2 clusters).
        wcss : float
            Within-cluster sum of squares (analogous to k-means inertia).
        """
        labels_raw, _ = ClusteringModule._get_clustering_dbscan(self, compositions_df)
        labels = ClusteringModule._normalize_dbscan_labels(labels_raw)

        X = compositions_df.to_numpy()
        k = len({int(lbl) for lbl in labels if lbl != -1})

        if k > 0:
            centroids = np.array([X[labels == i].mean(axis=0) for i in range(k)])
        else:
            centroids = np.empty((0, X.shape[1]))

        wcss = 0.0
        for i in range(k):
            cluster_points = X[labels == i]
            wcss += float(np.sum((cluster_points - centroids[i]) ** 2))

        # Silhouette is only defined for 2+ clusters; restrict to non-noise points.
        non_noise_mask = labels != -1
        if k >= 2 and int(np.sum(non_noise_mask)) > k:
            sil_score = float(silhouette_score(compositions_df[non_noise_mask], labels[non_noise_mask]))
        else:
            sil_score = np.nan

        return labels, centroids, k, sil_score, wcss


    def _compute_cluster_statistics(self, compositions_df, compositions_df_other_fr, centroids, labels):
        """
        Compute statistics for each cluster, including WCSS, standard deviations, and centroids
        in terms of both atomic and mass fractions.
    
        Returns
        -------
        wcss_per_cluster : list
            Within-Cluster Sum of Squares for each cluster.
        rms_dist_cluster : list
            Standard deviation of distances to centroid for each cluster.
        rms_dist_cluster_other_fr : list
            Standard deviation of distances in alternate fraction representation.
        n_points_per_cluster : list
            Number of points in each cluster.
        els_std_dev_per_cluster : list
            Elemental standard deviations within each cluster.
        els_std_dev_per_cluster_other_fr : list
            Elemental standard deviations in alternate fraction representation.
        centroids_other_fr : list
            Centroids in alternate fraction representation.
        max_cl_rmsdist : float
            Maximum standard deviation across all clusters.
        """
        wcss_per_cluster = []
        rms_dist_cluster = []
        rms_dist_cluster_other_fr = []
        n_points_per_cluster = []
        els_std_dev_per_cluster = []
        els_std_dev_per_cluster_other_fr = []
        centroids_other_fr = []
        for i, centroid in enumerate(centroids):
            # Save data using the elemental fraction employed as feature
            cluster_points = compositions_df[labels == i].to_numpy()
            n_points_per_cluster.append(len(cluster_points))
            if len(cluster_points) > 1:
                els_std_dev_per_cluster.append(np.std(cluster_points, axis=0, ddof=1))
            else:
                # Append NaN or zero or skip
                els_std_dev_per_cluster.append(np.full(cluster_points.shape[1], np.nan))
            distances = np.linalg.norm(cluster_points - centroid, axis=1)
            wcss_per_cluster.append(np.sum(distances ** 2))
            rms_dist_cluster.append(np.sqrt(np.mean(distances ** 2)))
    
            # Save data also for the other elemental fraction type
            cluster_points_other = compositions_df_other_fr[labels == i].to_numpy()
            centroid_other_fr = np.mean(cluster_points_other, axis=0)
            centroids_other_fr.append(centroid_other_fr)
            if len(cluster_points) > 1:
                els_std_dev_per_cluster_other_fr.append(np.std(cluster_points_other, axis=0, ddof=1))
            else:
                # Append NaN or zero or skip
                els_std_dev_per_cluster_other_fr.append(np.full(cluster_points_other.shape[1], np.nan))
            distances_other_fr = np.linalg.norm(cluster_points_other - centroid_other_fr, axis=1)
            rms_dist_cluster_other_fr.append(np.sqrt(np.mean(distances_other_fr ** 2)))
            
        max_cl_rmsdist = max(rms_dist_cluster)
        
        return (wcss_per_cluster, rms_dist_cluster, rms_dist_cluster_other_fr, n_points_per_cluster,
                els_std_dev_per_cluster, els_std_dev_per_cluster_other_fr, centroids_other_fr, max_cl_rmsdist)


    def _select_good_compositions(self, max_analytical_error):
        """
        Select compositions for clustering, filtering out those with high analytical error or bad quantification flags.
    
        Returns
        -------
        compositions_list_at : list
            List of atomic fractions for good spectra.
        compositions_list_w : list
            List of mass fractions for good spectra.
        unused_compositions_list : list
            List of compositions not used for clustering (for plotting).
        df_indices : list
            Indices of rows used for phase identification.
        n_datapts : int
            Total number of spectra considered.
        """
        # Initialise counters for spectra filtered out
        self.n_sp_too_low_counts = 0
        self.n_sp_too_high_an_err = 0
        self.n_sp_bad_quant = 0
    
        compositions_list_at = []
        compositions_list_w = []
        unused_compositions_list = []
        df_indices = []
        n_datapts = len(self.spectra_quant_records)
    
        for i in range(n_datapts):
            record = self.spectra_quant_records[i]
            if record is not None and record.composition_atomic_fractions is not None:
                is_comp_ok = True
                spectrum_quant_result_at = dict(record.composition_atomic_fractions)
                spectrum_quant_result_w = dict(record.composition_weight_fractions)
                analytical_error = float(record.analytical_error)
                quant_flag = record.quant_flag

                # Check if composition was flagged as bad during quantification
                if quant_flag not in self.clustering_cfg.quant_flags_accepted:
                    is_comp_ok = False
                    self.n_sp_bad_quant += 1
    
                elif max_analytical_error is None:
                    # Analytical error check is disabled
                    is_comp_ok = True
                    pass
    
                # Check if analytical error is too high
                elif analytical_error < - (max_analytical_error + self.undetectable_an_er) or analytical_error > max_analytical_error:
                    is_comp_ok = False
                    self.n_sp_too_high_an_err += 1
    
                # Append composition to list of used or unused datapoints
                if is_comp_ok:
                    df_indices.append(i)
                    # Construct dictionary that includes all elements that are supposed to be in the sample.
                    comp_at = {el: spectrum_quant_result_at.get(el, 0) for el in self.all_els_sample}
                    compositions_list_at.append(comp_at)
                    comp_w = {el: spectrum_quant_result_w.get(el, 0) for el in self.all_els_sample}
                    compositions_list_w.append(comp_w)
                else:
                    # Collect unused data points to show them in the clustering plot
                    if self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
                        comp = [spectrum_quant_result_at.get(el, 0) for el in self.all_els_sample]
                    elif self.clustering_cfg.features == cnst.W_FR_CL_FEAT:
                        comp = [spectrum_quant_result_w.get(el, 0) for el in self.all_els_sample]
                    unused_compositions_list.append(comp)
            else:
                self.n_sp_too_low_counts += 1
    
        return compositions_list_at, compositions_list_w, unused_compositions_list, df_indices, n_datapts
