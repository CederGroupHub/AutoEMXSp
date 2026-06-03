#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plotting mixin for composition analysis outputs."""

import importlib.util
import os
import warnings
from typing import Any, List

import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans

import autoemx.calibrations as calibs
from autoemx.core.composition_analysis import custom_plotting_builtin as builtin_custom_plotting
import autoemx.utils.constants as cnst
from autoemx.utils import print_single_separator, to_latex_formula

from autoemx._logging import get_logger
logger = get_logger(__name__)


class PlottingModule:
    # Attributes are injected by the analyzer class during composition analysis.
    plot_cfg: Any
    clustering_cfg: Any
    ref_phases_df: Any
    ref_formulae: Any
    sample_cfg: Any
    sample_id: str
    analysis_dir: str
    detectable_els_sample: List[str]
    all_els_sample: List[str]
    verbose: bool

    @staticmethod
    def _find_ideal_3d_azimuth(points_xyz: 'np.ndarray') -> float:
        """Return azimuth that maximizes XY footprint after Z-axis rotation."""
        points_xyz = np.asarray(points_xyz, dtype=float)
        if points_xyz.size == 0:
            return 35.0
        centered = points_xyz - np.mean(points_xyz, axis=0, keepdims=True)
        best_azimuth = 35.0
        best_area = -1.0
        for azimuth in np.arange(0.0, 360.0, 15.0):
            theta = np.deg2rad(azimuth)
            cos_t, sin_t = np.cos(theta), np.sin(theta)
            x_rot = centered[:, 0] * cos_t - centered[:, 1] * sin_t
            y_rot = centered[:, 0] * sin_t + centered[:, 1] * cos_t
            area = (np.max(x_rot) - np.min(x_rot)) * (np.max(y_rot) - np.min(y_rot))
            if area > best_area:
                best_area = area
                best_azimuth = float(azimuth)
        return best_azimuth

    def _load_custom_plot_function(self):
        """Load a user-defined custom plotting callable from plot config."""
        custom_plot_file = getattr(self.plot_cfg, "custom_plot_file", None)
        if not custom_plot_file:
            return None

        custom_plot_file = os.path.abspath(custom_plot_file)
        if not os.path.exists(custom_plot_file):
            warnings.warn(
                f"Custom plot file not found: {custom_plot_file}. Falling back to default plot.",
                UserWarning,
            )
            return None

        module_name = f"autoemx_user_custom_plot_{abs(hash(custom_plot_file))}"
        try:
            spec = importlib.util.spec_from_file_location(module_name, custom_plot_file)
            if spec is None or spec.loader is None:
                warnings.warn(
                    f"Could not load custom plot module from {custom_plot_file}. Falling back to default plot.",
                    UserWarning,
                )
                return None

            user_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(user_module)
            return getattr(user_module, "_save_clustering_plot_custom_3D", None)
        except Exception as exc:
            warnings.warn(
                f"Failed to import custom plotting module '{custom_plot_file}': {exc}. Falling back.",
                UserWarning,
            )
            return None

    def _run_custom_clustering_plot(
        self,
        elements: List[str],
        els_comps_list: 'np.ndarray',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list,
    ) -> bool:
        """Run custom clustering plotting code and return True on success."""
        custom_plot_func = PlottingModule._load_custom_plot_function(self)

        if custom_plot_func is None:
            custom_plot_func = getattr(builtin_custom_plotting, "_save_clustering_plot_custom_3D", None)
            if custom_plot_func is None:
                return False

        ideal_elev = None
        ideal_azim = None
        if len(elements) == 3:
            base_points = np.asarray(els_comps_list, dtype=float).T
            base_azimuth = PlottingModule._find_ideal_3d_azimuth(
                base_points if base_points.size > 0 else np.empty((0, 3))
            )
            ideal_elev = 24.0
            ideal_azim = (base_azimuth + 180.0) % 360.0

        try:
            try:
                custom_plot_func(
                    elements,
                    els_comps_list,
                    centroids,
                    labels,
                    els_std_dev_per_cluster,
                    unused_compositions_list,
                    self.clustering_cfg.features,
                    self.ref_phases_df,
                    self.ref_formulae,
                    self.plot_cfg.show_plots,
                    self.sample_id,
                    analysis_dir=self.analysis_dir,
                    output_filename=cnst.CUSTOM_CLUSTERING_PLOT_FILENAME + cnst.CLUSTERING_PLOT_FILEEXT,
                    ideal_elev=ideal_elev,
                    ideal_azim=ideal_azim,
                )
            except TypeError:
                # Backward compatibility for legacy custom plotting signatures.
                custom_plot_func(
                    elements,
                    els_comps_list,
                    centroids,
                    labels,
                    els_std_dev_per_cluster,
                    unused_compositions_list,
                    self.clustering_cfg.features,
                    self.ref_phases_df,
                    self.ref_formulae,
                    self.plot_cfg.show_plots,
                    self.sample_id,
                )
            return True
        except Exception as exc:
            warnings.warn(
                f"Custom plotting failed with '{exc}'. Falling back to default plot.",
                UserWarning,
            )
            return False

    def _save_plots(
        self,
        kmeans: 'KMeans',
        compositions_df: 'pd.DataFrame',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list
    ) -> None:
        # Silhouette plot (only if more than one cluster)
        if len(centroids) > 1:
            PlottingModule._save_silhouette_plot(
                kmeans, compositions_df, self.analysis_dir, show_plot=self.plot_cfg.show_plots
            )

        can_plot_clustering = True
        els_for_plot = list(set(self.detectable_els_sample) - set(self.plot_cfg.els_excluded_clust_plot))
        els_excluded_clust_plot = list(set(self.all_els_sample) - set(els_for_plot))
        n_els = len(els_for_plot)

        if n_els == 1:
            can_plot_clustering = False
            print_single_separator()
            warnings.warn("Cannot generate clustering plot with a single element.", UserWarning)
            if len(self.detectable_els_sample) > 1:
                logger.warning('⚠️ Too many elements were excluded from the clustering plot via the use of "els_excluded_clust_plot".')
                logger.info(f'ℹ️ Consider removing one or more among the list: {self.plot_cfg.els_excluded_clust_plot}')
        elif n_els > 3:
            els_excluded_clust_plot += els_for_plot[3:]
            els_for_plot = els_for_plot[:3]

        indices_to_remove = [self.all_els_sample.index(el) for el in els_excluded_clust_plot]
        els_for_plot = [el for i, el in enumerate(self.all_els_sample) if i not in indices_to_remove]
        centroids = np.array([[coord for i, coord in enumerate(row) if i not in indices_to_remove] for row in centroids])
        els_std_dev_per_cluster = [[stddev for i, stddev in enumerate(row) if i not in indices_to_remove] for row in els_std_dev_per_cluster]
        unused_compositions_list = [[fr for i, fr in enumerate(row) if i not in indices_to_remove] for row in unused_compositions_list]

        if can_plot_clustering:
            els_comps_list = compositions_df[els_for_plot].to_numpy().T
            if self.plot_cfg.use_custom_plots:
                custom_successful = PlottingModule._run_custom_clustering_plot(self,
                    els_for_plot,
                    els_comps_list,
                    centroids,
                    labels,
                    els_std_dev_per_cluster,
                    unused_compositions_list,
                )
                if not custom_successful:
                    PlottingModule._save_clustering_plot(self,
                        els_for_plot, els_comps_list, centroids, labels,
                        els_std_dev_per_cluster, unused_compositions_list
                    )
            else:
                PlottingModule._save_clustering_plot(self,
                    els_for_plot, els_comps_list, centroids, labels,
                    els_std_dev_per_cluster, unused_compositions_list
                )
        elif self.verbose:
            logger.warning('⚠️ Clusters were not plotted because only one detectable element was present in the sample.')
            undetectable = getattr(calibs, 'undetectable_els', [])
            logger.warning(f"⚠️ Elements {undetectable} cannot be detected at the employed instrument.")

    def _save_clustering_plot(
        self,
        elements: List[str],
        els_comps_list: 'np.ndarray',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list
    ) -> None:
        plt.rcParams['font.family'] = 'Arial'
        fontsize = 14
        labelpad = 12
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.titlesize'] = fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize

        axis_label_add = ' (w%)' if self.clustering_cfg.features == cnst.W_FR_CL_FEAT else ' (at%)'
        ticks = np.arange(0, 1, 0.1)
        ticks_labels = [f"{x*100:.0f}" for x in ticks]

        def _compute_zoom_limits(
            values: 'np.ndarray',
            margin_ratio: float = 0.08,
            min_span: float = 0.10,
            central_fraction: float = 0.90,
        ) -> tuple[float, float]:
            values = np.asarray(values, dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                return 0.0, 1.0

            lower_q = max(0.0, (1.0 - central_fraction) / 2.0)
            upper_q = min(1.0, 1.0 - lower_q)
            v_min = float(np.quantile(values, lower_q))
            v_max = float(np.quantile(values, upper_q))
            span = max(v_max - v_min, min_span)
            margin = span * margin_ratio
            low = max(0.0, v_min - margin)
            high = min(1.0, v_max + margin)

            if high <= low:
                center = 0.5 * (v_min + v_max)
                half = max(min_span * 0.5, 0.02)
                low = max(0.0, center - half)
                high = min(1.0, center + half)

            return low, high

        def _plot_clustering_scene(ax: Any, title_suffix: str = "", show_legend: bool = True) -> None:
            ax.scatter(*els_comps_list, c=labels, cmap='viridis', marker='o')
            ax.scatter(*centroids.T, c='red', marker='x', s=100, label='Centroids')

            first_ellipse = True
            for centroid, stdevs in zip(centroids, els_std_dev_per_cluster):
                if not np.any(np.isnan(stdevs)):
                    if len(elements) == 3:
                        x_c, y_c, z_c = centroid
                        rx, ry, rz = stdevs
                        u = np.linspace(0, 2 * np.pi, 100)
                        v = np.linspace(0, np.pi, 100)
                        x = x_c + rx * np.outer(np.cos(u), np.sin(v))
                        y = y_c + ry * np.outer(np.sin(u), np.sin(v))
                        z = z_c + rz * np.outer(np.ones_like(u), np.cos(v))
                        ax.plot_surface(x, y, z, color='red', alpha=0.1, edgecolor='none')
                        if first_ellipse:
                            first_ellipse = False
                            ax.plot([], [], [], color='red', alpha=0.1, label='Stddev')
                    else:
                        x_c, y_c = centroid
                        rx, ry = stdevs
                        ellipse = patches.Ellipse((x_c, y_c), rx, ry, edgecolor='red', facecolor='red', linestyle='--', alpha=0.2)
                        if first_ellipse:
                            ellipse.set_label('Stddev')
                            first_ellipse = False
                        ax.add_patch(ellipse)

            if unused_compositions_list and self.plot_cfg.show_unused_comps_clust:
                ax.scatter(*np.array(unused_compositions_list).T, c='grey', marker='^', label='Discarded comps.')

            if self.ref_formulae is not None:
                first_ref = True
                ref_phases_df = self.ref_phases_df[elements]
                for index, row in ref_phases_df.iterrows():
                    label = 'Candidate phases' if first_ref else None
                    ax.scatter(*row.values, c='blue', marker='*', s=100, label=label)
                    ref_label = to_latex_formula(self.ref_formulae[index])
                    ax.text(*row.values, ref_label, color='black', fontsize=fontsize, ha='left', va='bottom')
                    first_ref = False

            for i, centroid in enumerate(centroids):
                ax.text(*centroid, str(i), color='black', fontsize=fontsize, ha='right', va='bottom')

            ax.set_xlabel(elements[0] + axis_label_add, labelpad=labelpad)
            ax.set_ylabel(elements[1] + axis_label_add, labelpad=labelpad)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticks(ticks)
            ax.set_xticklabels(ticks_labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(ticks_labels)
            if len(elements) == 3:
                ax.set_zlabel(elements[2] + axis_label_add, labelpad=labelpad * 0.95)
                # Keep (0,0,0) at the back while preserving the chosen camera angle.
                ax.set_xlim(1, 0)
                ax.set_ylim(1, 0)
                ax.set_zlim(0, 1)
                ax.set_zticks(ticks)
                ax.set_zticklabels(ticks_labels)
            ax.set_title(f'{self.clustering_cfg.method} clustering {self.sample_id}{title_suffix}')

            if show_legend and getattr(self.plot_cfg, 'show_legend_clustering', None):
                ax.legend(fontsize=fontsize, loc='best')

        fig = plt.figure(figsize=(6, 6))
        full_view_elev = 5.0
        full_view_azim = None
        if len(elements) == 3:
            ax: Any = fig.add_subplot(111, projection='3d')
            base_points = np.asarray(els_comps_list, dtype=float).T
            base_azimuth = PlottingModule._find_ideal_3d_azimuth(base_points if base_points.size > 0 else np.empty((0, 3)))
            # Keep 3D axes on the far side for a clearer foreground view of clusters.
            full_view_azim = (base_azimuth + 180.0) % 360.0
            ax.view_init(elev=full_view_elev, azim=full_view_azim)
        else:
            ax: Any = fig.add_subplot(111)
        _plot_clustering_scene(ax, show_legend=True)
        if self.plot_cfg.show_plots:
            plt.show()
        fig.savefig(
            os.path.join(self.analysis_dir, cnst.CLUSTERING_PLOT_FILENAME + cnst.CLUSTERING_PLOT_FILEEXT),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.1,
        )

        fig_zoomed = plt.figure(figsize=(6, 6))
        if len(elements) == 3:
            ax_zoomed: Any = fig_zoomed.add_subplot(111, projection='3d')
        else:
            ax_zoomed: Any = fig_zoomed.add_subplot(111)
        _plot_clustering_scene(ax_zoomed, title_suffix=' (zoomed)', show_legend=False)

        # Zoom includes most of total sample points, including discarded compositions.
        zoom_points = [np.asarray(els_comps_list, dtype=float).T]
        if unused_compositions_list:
            zoom_points.append(np.asarray(unused_compositions_list, dtype=float))
        zoom_points.append(np.asarray(centroids, dtype=float))
        # Add reference phases within 10% of any outlier to zoom extent
        if self.ref_formulae is not None and unused_compositions_list:
            ref_phases_df_zoom = self.ref_phases_df[elements]
            unused_arr = np.asarray(unused_compositions_list, dtype=float)
            threshold = 20  # 20% distance
            for _, row in ref_phases_df_zoom.iterrows():
                ref_point = np.array(row.values, dtype=float)
                dists = np.linalg.norm(unused_arr - ref_point, axis=1)
                if np.any(dists < threshold):
                    zoom_points.append(ref_point.reshape(1, -1))
        all_points = np.vstack([pts for pts in zoom_points if pts.size > 0]) if zoom_points else np.empty((0, len(elements)))

        # Compute zoom limits and expand by 20% of the span for better visibility
        def expand_limits(low, high, margin_ratio=0.20):
            span = high - low
            margin = span * margin_ratio
            return max(0.0, low - margin), min(1.0, high + margin)

        x_low, x_high = _compute_zoom_limits(all_points[:, 0] if all_points.size > 0 else np.array([]))
        y_low, y_high = _compute_zoom_limits(all_points[:, 1] if all_points.size > 0 else np.array([]))
        x_low, x_high = expand_limits(x_low, x_high)
        y_low, y_high = expand_limits(y_low, y_high)
        if len(elements) == 3:
            ax_zoomed.set_xlim(x_high, x_low)
            ax_zoomed.set_ylim(y_high, y_low)
        else:
            ax_zoomed.set_xlim(x_low, x_high)
            ax_zoomed.set_ylim(y_low, y_high)

        int_percent_formatter = FuncFormatter(lambda value, _: f"{int(round(value * 100.0))}")
        ax_zoomed.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax_zoomed.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ax_zoomed.xaxis.set_major_formatter(int_percent_formatter)
        ax_zoomed.yaxis.set_major_formatter(int_percent_formatter)

        if len(elements) == 3:
            z_low, z_high = _compute_zoom_limits(all_points[:, 2] if all_points.size > 0 else np.array([]))
            z_low, z_high = expand_limits(z_low, z_high)
            ax_zoomed.set_zlim(z_low, z_high)
            ax_zoomed.zaxis.set_major_locator(MaxNLocator(nbins=6))
            ax_zoomed.zaxis.set_major_formatter(int_percent_formatter)
            # Keep identical orientation to the full plot for direct visual comparison.
            if full_view_azim is not None:
                ax_zoomed.view_init(elev=full_view_elev, azim=full_view_azim)

        fig_zoomed.savefig(
            os.path.join(
                self.analysis_dir,
                cnst.CLUSTERING_PLOT_FILENAME + '_zoomed' + cnst.CLUSTERING_PLOT_FILEEXT,
            ),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.1,
        )

    def _save_violin_plot_powder_mixture(
        self,
        W_mol_frs: 'np.ndarray',
        ref_names: List[str],
        cluster_ID: int
    ) -> None:
        plt.rcParams['font.family'] = 'Arial'
        fontsize = 17
        labelpad = 0
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.titlesize'] = fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize
        purple_cmap = cm.get_cmap('Purples')
        yellow_cmap = cm.get_cmap('autumn')

        y_vals = np.asarray(W_mol_frs, dtype=float)[:, 0]
        fig, ax_left = plt.subplots(figsize=(4, 4))
        mean = np.mean(y_vals)
        std = np.std(y_vals)

        ax_left = sns.violinplot(data=y_vals, inner=None, color=purple_cmap(0.3), linewidth=1.5, density_norm='area', width=1, zorder=1)
        sns.swarmplot(data=y_vals, color=purple_cmap(0.8), edgecolor=purple_cmap(1.0), linewidth=2, size=5, label='data', zorder=2)
        ax_left.errorbar(0, mean, yerr=std / 2, fmt='none', color=yellow_cmap(0.9), label='Mean ±1 Std Dev', capsize=5, elinewidth=1, zorder=4, markerfacecolor=yellow_cmap(0.9), markeredgecolor='black', markeredgewidth=1, marker='o', linestyle='none')
        ax_left.errorbar(0, mean, yerr=std / 2, fmt='none', color='none', label='_nolegend_', capsize=6, elinewidth=2, zorder=3, markerfacecolor='none', markeredgecolor='black', markeredgewidth=2, marker='o', linestyle='none', ecolor='black')
        ax_left.scatter(0, mean, color=yellow_cmap(0.9), marker='o', s=50, edgecolors='k', linewidths=1, label='Mean', zorder=10)

        ax_left.set_xticks([])
        ax_left.set_yticks([0, 1])
        ax_left.set_frame_on(True)
        for spine in ax_left.spines.values():
            spine.set_color('black')
            spine.set_linewidth(0.5)
        plt.grid(False)

        plt.xlim(-0.5, 0.5)
        ylim_bottom = 0
        ylim_top = 1
        ax_left.set_ylim(ylim_bottom, ylim_top)

        left_formula = to_latex_formula(ref_names[0], include_dollar_signs=False)
        ax_left.set_ylabel(rf"$x_{{\mathrm{{{left_formula}}}}}$", labelpad=labelpad)
        ax_right = ax_left.twinx()
        ax_right.set_ylim(ylim_top, ylim_bottom)
        ax_right.set_yticks([1, 0])
        right_formula = to_latex_formula(ref_names[1], include_dollar_signs=False)
        ax_right.set_ylabel(rf"$x_{{\mathrm{{{right_formula}}}}}$", labelpad=labelpad)
        ax_left.text(0.03, 0.03, rf"$\sigma_x = {std*100:.1f}$%", fontsize=fontsize, color='black', ha='left', va='bottom', transform=ax_left.transAxes)
        ax_left.set_title(f'Violin plot {self.sample_id}')

        fig.savefig(
            os.path.join(self.analysis_dir, cnst.POWDER_MIXTURE_PLOT_FILENAME + f"_cl{cluster_ID}_{ref_names[0]}_{ref_names[1]}" + cnst.CLUSTERING_PLOT_FILEEXT),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0,
        )

    @staticmethod
    def _save_silhouette_plot(
        model: 'KMeans',
        compositions_df: 'pd.DataFrame',
        results_dir: str,
        show_plot: bool
    ) -> None:
        try:
            yellowbrick_cluster = importlib.import_module('yellowbrick.cluster')
            silhouette_visualizer_cls = getattr(yellowbrick_cluster, 'SilhouetteVisualizer', None)
        except Exception:
            silhouette_visualizer_cls = None

        if silhouette_visualizer_cls is None:
            warnings.warn(
                "yellowbrick is not available; skipping silhouette plot generation.",
                UserWarning,
            )
            return

        plt.figure(figsize=(10, 8))
        sil_visualizer = silhouette_visualizer_cls(model, colors='yellowbrick')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            sil_visualizer.fit(compositions_df)

        plt.ylabel('Cluster label')
        plt.xlabel('Silhouette coefficient values')
        plt.legend(loc='upper right', frameon=True)

        if show_plot:
            plt.ion()
            sil_visualizer.show()
            plt.pause(0.001)
            plt.ioff()

        fig = sil_visualizer.fig
        fig.savefig(os.path.join(results_dir, 'Silhouette_plot.png'))
        if not show_plot:
            plt.close(fig)
