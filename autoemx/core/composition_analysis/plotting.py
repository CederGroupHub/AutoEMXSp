#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plotting mixin for composition analysis outputs."""

import importlib.util
import os
import warnings
from typing import Any, List, Optional, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans

import autoemx.calibrations as calibs
from autoemx.core.composition_analysis import custom_plotting_builtin as builtin_custom_plotting
from autoemx.core.composition_analysis.clustering_plot_axes import (
    CLUSTERING_3D_VIEW_AZIM,
    CLUSTERING_3D_VIEW_ELEV,
    apply_data_driven_axis_limits,
    apply_fixed_full_range_ticks,
    compute_data_driven_axis_limits,
    configure_interactive_clustering_axes,
    gather_clustering_zoom_points,
    point_within_composition_limits,
)
import autoemx.utils.constants as cnst
from autoemx.utils.helper import print_single_separator, to_latex_formula

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
        kmeans: 'Optional[KMeans]',
        compositions_df: 'pd.DataFrame',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list,
        silhouette_df: 'Optional[pd.DataFrame]' = None
    ) -> None:
        # Silhouette plot (only if more than one cluster). The Yellowbrick
        # visualizer requires a fitted KMeans model, so it is skipped for methods
        # (e.g. DBSCAN) that do not provide one. It refits the model, so it must use
        # the features clustering operated on (silhouette_df, e.g. CLR coordinates).
        if kmeans is not None and len(centroids) > 1:
            PlottingModule._save_silhouette_plot(
                kmeans, compositions_df if silhouette_df is None else silhouette_df,
                self.analysis_dir, show_plot=self.plot_cfg.show_plots
            )

        can_plot_clustering = True
        els_to_plot = list(dict.fromkeys(getattr(self.plot_cfg, "els_to_plot", None) or []))
        excluded = list(dict.fromkeys(self.plot_cfg.els_excluded_clust_plot or []))

        if els_to_plot:
            unknown = [el for el in els_to_plot if el not in self.detectable_els_sample]
            if unknown:
                raise ValueError(
                    f'els_to_plot contains elements not among the sample\'s detectable elements '
                    f'{list(self.detectable_els_sample)}: {unknown}'
                )
            if len(els_to_plot) > 3:
                raise ValueError(
                    f'els_to_plot must contain at most 3 elements for the clustering plot, '
                    f'got {len(els_to_plot)}: {els_to_plot}'
                )

            conflict = [el for el in els_to_plot if el in excluded]
            if conflict:
                logger.warning(
                    '⚠️ Removing %s from els_excluded_clust_plot because els_to_plot forces them onto the axes.',
                    conflict,
                )
                excluded = [el for el in excluded if el not in set(els_to_plot)]

            if len(els_to_plot) in (2, 3):
                # Authoritative selection: plot exactly these elements, in the given order.
                els_for_plot = list(els_to_plot)
            else:
                # A single forced element: keep it and fill from remaining detectable elements.
                els_for_plot = list(els_to_plot)
                for el in self.detectable_els_sample:
                    if el not in els_for_plot and el not in excluded:
                        els_for_plot.append(el)
                if len(els_for_plot) > 3:
                    els_for_plot = els_for_plot[:3]

            if len(els_for_plot) not in (2, 3):
                can_plot_clustering = False
                print_single_separator()
                warnings.warn(
                    f"Cannot generate clustering plot: els_to_plot={els_to_plot} leaves "
                    f"{len(els_for_plot)} axis element(s) ({els_for_plot}); need 2 or 3.",
                    UserWarning,
                )
                logger.warning(
                    '⚠️ els_to_plot must result in 2 or 3 plot axes. '
                    f'Got {len(els_for_plot)} element(s): {els_for_plot}. '
                    'Provide 2 or 3 elements in els_to_plot, or remove exclusions that leave too few.'
                )
        else:
            els_for_plot = [el for el in self.all_els_sample if el in self.detectable_els_sample and el not in excluded]
            n_els = len(els_for_plot)

            if n_els == 1:
                can_plot_clustering = False
                print_single_separator()
                warnings.warn("Cannot generate clustering plot with a single element.", UserWarning)
                if len(self.detectable_els_sample) > 1:
                    logger.warning('⚠️ Too many elements were excluded from the clustering plot via the use of "els_excluded_clust_plot".')
                    logger.info(f'ℹ️ Consider removing one or more among the list: {self.plot_cfg.els_excluded_clust_plot}')
            elif n_els > 3:
                els_for_plot = els_for_plot[:3]

        # 2D clustering plot (ternary diagram or two most representative elements), from the full compositions
        PlottingModule._save_clustering_plot_2d(
            self, compositions_df, centroids, labels, els_std_dev_per_cluster, unused_compositions_list
        )

        indices_to_keep = [self.all_els_sample.index(el) for el in els_for_plot]
        centroids = np.array([[row[i] for i in indices_to_keep] for row in centroids])
        els_std_dev_per_cluster = [[row[i] for i in indices_to_keep] for row in els_std_dev_per_cluster]
        unused_compositions_list = [[row[i] for i in indices_to_keep] for row in unused_compositions_list]

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
        is_3d = len(elements) == 3

        def _plot_clustering_scene(
            ax: Any,
            title_suffix: str = "",
            show_legend: bool = True,
            use_fixed_full_range_ticks: bool = True,
            ref_phase_limits: Optional[tuple[tuple[float, float], tuple[float, float], Optional[tuple[float, float]]]] = None,
        ) -> None:
            labels_arr = np.asarray(labels)
            noise_mask = labels_arr == -1
            comps_arr = np.asarray(els_comps_list)
            if np.any(noise_mask):
                # DBSCAN noise points: render in grey, separate from clustered points.
                ax.scatter(*comps_arr[:, ~noise_mask], c=labels_arr[~noise_mask], cmap='viridis', marker='o')
                ax.scatter(*comps_arr[:, noise_mask], c='lightgrey', marker='o', label='Noise (unclustered)')
            else:
                ax.scatter(*comps_arr, c=labels_arr, cmap='viridis', marker='o')
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
                ref_xlim = ref_ylim = ref_zlim = None
                if ref_phase_limits is not None:
                    ref_xlim, ref_ylim, ref_zlim = ref_phase_limits
                for index, row in ref_phases_df.iterrows():
                    if ref_xlim is not None and ref_ylim is not None and not point_within_composition_limits(
                        row.values, ref_xlim, ref_ylim, ref_zlim
                    ):
                        continue
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
            if is_3d:
                ax.set_zlabel(elements[2] + axis_label_add, labelpad=labelpad * 0.95)
                # Keep (0,0,0) at the back while preserving the chosen camera angle.
                ax.set_xlim(1, 0)
                ax.set_ylim(1, 0)
                ax.set_zlim(0, 1)
            if use_fixed_full_range_ticks:
                apply_fixed_full_range_ticks(ax, is_3d=is_3d)
            ax.set_title(f'{self.clustering_cfg.method} clustering {self.sample_id}{title_suffix}')

            if show_legend and getattr(self.plot_cfg, 'show_legend_clustering', None):
                ax.legend(fontsize=fontsize, loc='best')

        fig = plt.figure(figsize=(6, 6))
        if len(elements) == 3:
            ax: Any = fig.add_subplot(111, projection='3d')
            ax.view_init(elev=CLUSTERING_3D_VIEW_ELEV, azim=CLUSTERING_3D_VIEW_AZIM)
        else:
            ax: Any = fig.add_subplot(111)
        _plot_clustering_scene(ax, show_legend=True)
        if is_3d:
            # Reserve extra space on the right so the z-axis label is visible in exports.
            fig.subplots_adjust(right=0.88)
        fig.savefig(
            os.path.join(self.analysis_dir, cnst.CLUSTERING_PLOT_FILENAME + cnst.CLUSTERING_PLOT_FILEEXT),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.2 if is_3d else 0.1,
        )
        if self.plot_cfg.show_plots:
            configure_interactive_clustering_axes(
                ax,
                is_3d=is_3d,
                reversed_xy=is_3d,
            )
        else:
            plt.close(fig)

        all_points = gather_clustering_zoom_points(
            els_comps_list,
            centroids,
            unused_compositions_list,
            elements,
            ref_phases_df=self.ref_phases_df if self.ref_formulae is not None else None,
        )
        zoom_limits = compute_data_driven_axis_limits(all_points, is_3d=is_3d)

        fig_zoomed = plt.figure(figsize=(6, 6))
        if len(elements) == 3:
            ax_zoomed: Any = fig_zoomed.add_subplot(111, projection='3d')
        else:
            ax_zoomed: Any = fig_zoomed.add_subplot(111)
        _plot_clustering_scene(
            ax_zoomed,
            title_suffix=' (zoomed)',
            show_legend=False,
            use_fixed_full_range_ticks=False,
            ref_phase_limits=zoom_limits,
        )
        apply_data_driven_axis_limits(
            ax_zoomed,
            all_points,
            is_3d=is_3d,
            reversed_xy=is_3d,
        )
        if is_3d:
            ax_zoomed.view_init(elev=CLUSTERING_3D_VIEW_ELEV, azim=CLUSTERING_3D_VIEW_AZIM)
            # Match base-plot spacing to avoid clipping the z-axis label.
            fig_zoomed.subplots_adjust(right=0.88)

        fig_zoomed.savefig(
            os.path.join(
                self.analysis_dir,
                cnst.CLUSTERING_PLOT_FILENAME + '_zoomed' + cnst.CLUSTERING_PLOT_FILEEXT,
            ),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.2 if is_3d else 0.1,
        )
        if not self.plot_cfg.show_plots:
            plt.close(fig_zoomed)

    def _clustering_plot_2d_elements(self, compositions_df: 'pd.DataFrame', centroids: 'np.ndarray') -> Tuple[List[str], bool]:
        """
        Elements of the 2D clustering plot, and whether they are shown in a ternary diagram.

        With exactly 3 detectable elements (after exclusions) and atomic fractions, all 3 are shown in a
        ternary diagram. Otherwise the 2 most representative elements are shown: those listed in
        ``els_to_plot`` if it has exactly 2 entries, else those whose fractions differ most between the
        cluster centroids (or vary most across compositions, with a single cluster).
        """
        excluded = set(self.plot_cfg.els_excluded_clust_plot or [])
        els = [el for el in self.all_els_sample if el in self.detectable_els_sample and el not in excluded]
        forced = [el for el in (getattr(self.plot_cfg, 'els_to_plot', None) or []) if el in self.detectable_els_sample]
        if len(forced) == 3 and self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
            return forced, True
        if len(forced) == 2:
            return forced, False
        if len(els) == 3 and self.clustering_cfg.features == cnst.AT_FR_CL_FEAT:
            return els, True
        C = np.atleast_2d(np.asarray(centroids, dtype=float))
        idx = [self.all_els_sample.index(el) for el in els]
        spread = C[:, idx].std(axis=0) if len(C) > 1 else compositions_df[els].to_numpy(dtype=float).std(axis=0)
        top = sorted(np.argsort(-spread, kind='stable')[:2])
        return [els[i] for i in top], False


    def _save_clustering_plot_2d(
        self,
        compositions_df: 'pd.DataFrame',
        centroids: 'np.ndarray',
        labels: 'np.ndarray',
        els_std_dev_per_cluster: list,
        unused_compositions_list: list,
    ) -> None:
        """
        2D clustering plot, for samples with at least 3 detectable elements (full and zoomed).

        Ternary diagram for 3 elements in atomic fractions; otherwise the 2 most representative elements,
        in fractions (see ``_clustering_plot_2d_elements``). Saved as Clustering_plot_2D(_zoomed).png.
        """
        n_detectable = len([el for el in self.all_els_sample if el in self.detectable_els_sample])
        if n_detectable < 3:
            return
        els, is_ternary = PlottingModule._clustering_plot_2d_elements(self, compositions_df, centroids)
        if len(els) < 2:
            return
        fontsize = PlottingModule._mixture_style()
        unit = 'w%' if self.clustering_cfg.features == cnst.W_FR_CL_FEAT else 'at%'
        idx = [self.all_els_sample.index(el) for el in els]
        X = compositions_df[els].to_numpy(dtype=float)
        C = np.atleast_2d(np.asarray(centroids, dtype=float))[:, idx]
        S = np.atleast_2d(np.asarray(els_std_dev_per_cluster, dtype=float))[:, idx] if len(els_std_dev_per_cluster) else None
        U = np.atleast_2d(np.asarray(unused_compositions_list, dtype=float))[:, idx] if unused_compositions_list else None
        R = None
        if self.ref_formulae:
            R = self.ref_phases_df[els].to_numpy(dtype=float)
        if is_ternary:
            proj = PlottingModule._ternary_xy
        else:
            proj = lambda A: np.atleast_2d(A) * 100
        labels_arr = np.asarray(labels)
        P, Pc = proj(X), proj(C)
        Pu = proj(U) if U is not None and len(U) and self.plot_cfg.show_unused_comps_clust else None
        Pr = proj(R) if R is not None and len(R) else None

        for zoomed in (False, True):
            fig, ax = plt.subplots(figsize=(6, 6))
            if zoomed:
                view = PlottingModule._square_view(np.vstack([P, Pc]))
            else:
                view = ((-0.12, 1.12), (-0.12, 0.98)) if is_ternary else ((0, 100), (0, 100))
            corner_labels = []
            if is_ternary:
                PlottingModule._draw_ternary_grid(ax, 0.05 if zoomed else 0.1)
                if not zoomed:
                    corners = PlottingModule._ternary_xy(np.eye(3))
                    for (x, y), el, (dx, dy) in zip(corners, els, [(-0.03, -0.06), (0.03, -0.06), (0, 0.05)]):
                        corner_labels.append(ax.text(x + dx, y + dy, f'{el} ({unit})', ha='center', va='center', fontsize=fontsize))
                ax.axis('off')
            else:
                ax.set_xlabel(f'{els[0]} ({unit})')
                ax.set_ylabel(f'{els[1]} ({unit})')
            noise = labels_arr == -1
            ax.scatter(P[~noise, 0], P[~noise, 1], c=labels_arr[~noise], cmap='viridis', marker='o', s=25, zorder=2)
            if np.any(noise):
                ax.scatter(P[noise, 0], P[noise, 1], c='lightgrey', marker='o', s=25, label='Noise (unclustered)', zorder=2)
            if Pu is not None:
                ax.scatter(Pu[:, 0], Pu[:, 1], c='grey', marker='^', s=25, label='Discarded comps.', zorder=1)
            if not is_ternary and S is not None:
                first = True
                for c_xy, sd in zip(Pc, S * 100):
                    if not np.any(np.isnan(sd)):
                        ax.add_patch(patches.Ellipse(tuple(c_xy), sd[0], sd[1], edgecolor='red', facecolor='red',
                                                     linestyle='--', alpha=0.2, label='Stddev' if first else None))
                        first = False
            ax.scatter(Pc[:, 0], Pc[:, 1], c='red', marker='x', s=100, label='Centroids', zorder=4)
            texts = [ax.text(*c_xy, str(i), color='black', fontsize=fontsize - 2, ha='right', va='bottom',
                             bbox=PlottingModule.LABEL_BOX, zorder=3) for i, c_xy in enumerate(Pc)]
            if Pr is not None:
                for f, r in zip(self.ref_formulae, Pr):
                    if not (view[0][0] <= r[0] <= view[0][1] and view[1][0] <= r[1] <= view[1][1]):
                        continue
                    ax.scatter(*r, c='blue', marker='*', s=100, label='Candidate phases' if len(texts) == len(Pc) else None, zorder=4)
                    start = PlottingModule._label_start(r, PlottingModule._ternary_xy(np.eye(3)) if is_ternary else None)
                    texts.append(ax.text(*start, to_latex_formula(f), color='black', fontsize=fontsize - 3, ha='left',
                                         va='bottom', bbox=PlottingModule.LABEL_BOX, zorder=3))
            ax.set_xlim(*view[0])
            ax.set_ylim(*view[1])
            ax.set_aspect('equal')
            avoid = [P, Pc] + ([Pr] if Pr is not None else [])
            PlottingModule._adjust_labels(ax, texts, np.vstack(avoid), corner_labels)
            ax.set_title(f'{self.clustering_cfg.method} clustering {self.sample_id}' + (' (zoomed)' if zoomed else ''))
            if not zoomed and getattr(self.plot_cfg, 'show_legend_clustering', None):
                ax.legend(fontsize=fontsize - 3, loc='best')
            fig.savefig(
                os.path.join(self.analysis_dir, cnst.CLUSTERING_PLOT_FILENAME + '_2D' + ('_zoomed' if zoomed else '')
                             + cnst.CLUSTERING_PLOT_FILEEXT),
                dpi=300, bbox_inches='tight', pad_inches=0.1,
            )
            if not self.plot_cfg.show_plots:
                plt.close(fig)


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
        purple_cmap = plt.get_cmap('Purples')
        yellow_cmap = plt.get_cmap('autumn')

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

    def _save_violin_plot_powder_mixture_multi(
        self,
        W_mol_frs: 'np.ndarray',
        ref_names: List[str],
        cluster_ID: int
    ) -> None:
        """Violin plots of the molar fractions of each phase of a mixture of 3 or more phases, side by side."""
        plt.rcParams['font.family'] = 'Arial'
        fontsize = 17
        plt.rcParams['font.size'] = fontsize
        plt.rcParams['axes.titlesize'] = fontsize
        plt.rcParams['axes.labelsize'] = fontsize
        plt.rcParams['xtick.labelsize'] = fontsize
        plt.rcParams['ytick.labelsize'] = fontsize
        purple_cmap = plt.get_cmap('Purples')
        yellow_cmap = plt.get_cmap('autumn')

        x = np.asarray(W_mol_frs, dtype=float)
        n_phases = x.shape[1]
        long_df = pd.DataFrame({
            'phase': np.repeat(np.arange(n_phases), len(x)),
            'x': x.T.ravel(),
        })
        fig, ax = plt.subplots(figsize=(1.8 * n_phases + 1, 4))
        sns.violinplot(data=long_df, x='phase', y='x', inner=None, color=purple_cmap(0.3), linewidth=1.5,
                       density_norm='area', width=0.9, cut=0, zorder=1, ax=ax)
        sns.swarmplot(data=long_df, x='phase', y='x', color=purple_cmap(0.8), edgecolor=purple_cmap(1.0),
                      linewidth=1, size=3, zorder=2, ax=ax)
        means, stds = x.mean(axis=0), x.std(axis=0)
        ax.errorbar(np.arange(n_phases), means, yerr=stds / 2, fmt='none', ecolor='black', capsize=6, elinewidth=2, zorder=3)
        ax.errorbar(np.arange(n_phases), means, yerr=stds / 2, fmt='none', ecolor=yellow_cmap(0.9), capsize=5, elinewidth=1, zorder=4)
        ax.scatter(np.arange(n_phases), means, color=yellow_cmap(0.9), marker='o', s=50, edgecolors='k', linewidths=1, zorder=10)

        ax.set_xticks(np.arange(n_phases))
        ax.set_xticklabels([f"{to_latex_formula(f)}\n$\\sigma_x = {sd * 100:.1f}$%" for f, sd in zip(ref_names, stds)])
        ax.set_xlabel('')
        ax.set_ylabel('$x$')
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1])
        ax.set_frame_on(True)
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(0.5)
        ax.grid(False)
        ax.set_title(f'Violin plot {self.sample_id}')

        fig.savefig(
            os.path.join(self.analysis_dir, cnst.POWDER_MIXTURE_PLOT_FILENAME + f"_cl{cluster_ID}_" + '_'.join(ref_names) + cnst.CLUSTERING_PLOT_FILEEXT),
            dpi=300,
            bbox_inches='tight',
            pad_inches=0.05,
        )
        plt.close(fig)

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
            # Finalize without displaying; all figures are shown together at the end of the analysis
            sil_visualizer.finalize()

        fig = sil_visualizer.fig
        fig.savefig(os.path.join(results_dir, 'Silhouette_plot.png'))
        if not show_plot:
            plt.close(fig)


    # =========================================================================
    # Best mixture of each cluster
    # =========================================================================

    def _save_best_mixture_plots(
        self,
        compositions_df: 'pd.DataFrame',
        labels: 'np.ndarray',
        clusters_assigned_mixtures: List[List[dict]],
    ) -> None:
        """
        Plot the top-ranked mixture of each cluster, if it combines 2 or 3 candidate phases.

        - 3 detectable elements: ternary diagram of the elements (full and zoomed).
        - 2 phases, other numbers of elements: 2D plot of the two elements that differ most
          between the phases, or of ``plot_cfg.els_to_plot`` if it lists two elements (full and zoomed).
        - 3 phases, 4+ elements: ternary diagram of the molar fractions of the three phases.

        Mixtures from free NMF (phases not among ``ref_formulae``) are not plotted.
        """
        # Imported here: reference_matching imports this module
        from autoemx.core.composition_analysis.reference_matching import ReferenceMatchingModule

        els = list(self.detectable_els_sample)
        ref_formulae = list(self.ref_formulae or [])
        unit = 'w%' if self.clustering_cfg.features == cnst.W_FR_CL_FEAT else 'at%'
        labels_arr = np.asarray(labels)

        for cluster_id, mixtures in enumerate(clusters_assigned_mixtures or []):
            if not mixtures:
                continue
            top = ReferenceMatchingModule._sort_mixtures(mixtures)[0]
            phases = list(top[cnst.REF_NAME_KEY])
            if len(phases) not in (2, 3) or not all(f in ref_formulae for f in phases):
                continue
            X = compositions_df[els].to_numpy()[labels_arr == cluster_id]
            if len(X) == 0:
                continue
            H = np.array([self.ref_phases_df[els].iloc[ref_formulae.index(f)].to_numpy() for f in phases], dtype=float)
            W, _ = ReferenceMatchingModule._nmf_with_constraints(self, X, len(phases), fixed_H=H)
            other_refs = {
                f: self.ref_phases_df[els].iloc[i].to_numpy(dtype=float)
                for i, f in enumerate(ref_formulae) if f not in phases
            }
            info = PlottingModule._mixture_info_text(top, phases)
            title = f'Mixture {self.sample_id} cluster {cluster_id}'
            filename = cnst.BEST_MIXTURE_PLOT_FILENAME + f'_cl{cluster_id}'

            try:
                if len(els) == 3:
                    for zoomed in (False, True):
                        fig = PlottingModule._plot_mixture_element_ternary(
                            X, W, H, phases, other_refs, els, unit, title, info, zoomed
                        )
                        PlottingModule._finish_mixture_figure(self, fig, filename + ('_zoomed' if zoomed else ''))
                elif len(phases) == 2:
                    els_xy = PlottingModule._mixture_plot_elements(self, H, els)
                    idx = [els.index(el) for el in els_xy]
                    for zoomed in (False, True):
                        fig = PlottingModule._plot_mixture_element_2d(
                            X[:, idx], (W @ H)[:, idx], H[:, idx],
                            {f: v[idx] for f, v in other_refs.items()}, phases, els_xy, unit, title, info, zoomed
                        )
                        PlottingModule._finish_mixture_figure(self, fig, filename + ('_zoomed' if zoomed else ''))
                elif len(els) > 3:
                    mol_frs = W / np.asarray([self.ref_weights_in_mixture[ref_formulae.index(f)] for f in phases], dtype=float)
                    mol_frs /= mol_frs.sum(axis=1, keepdims=True)
                    dist = np.linalg.norm(X - W @ H, axis=1) * 100
                    fig = PlottingModule._plot_mixture_phase_ternary(mol_frs, dist, phases, unit, title, info)
                    PlottingModule._finish_mixture_figure(self, fig, filename)
            except Exception as exc:
                logger.warning(f"⚠️ Could not plot the mixture of cluster {cluster_id}: {exc}")

    @staticmethod
    def _mixture_style() -> int:
        fontsize = 14
        plt.rcParams['font.family'] = 'Arial'
        for key in ('font.size', 'axes.titlesize', 'axes.labelsize', 'xtick.labelsize', 'ytick.labelsize'):
            plt.rcParams[key] = fontsize
        return fontsize

    @staticmethod
    def _mixture_info_text(mixture: dict, phases: List[str]) -> str:
        fractions = mixture.get(cnst.MOLAR_FRS_MEAN_KEY) or []
        frs = ', '.join(f'{to_latex_formula(f)} {x * 100:.0f}%' for f, x in zip(phases, fractions))
        return f"$CS_{{mix}}$ = {mixture[cnst.CONF_SCORE_KEY]:.2f}   |   {frs}"

    def _mixture_plot_elements(self, H: 'np.ndarray', els: List[str]) -> List[str]:
        """Axes of the 2D mixture plot: plot_cfg.els_to_plot if it lists 2 elements, else the 2 elements differing most between the phases."""
        forced = [el for el in (getattr(self.plot_cfg, 'els_to_plot', None) or []) if el in els]
        if len(forced) == 2:
            return forced
        order = np.argsort(-np.abs(H[0] - H[1]), kind='stable')
        return [els[i] for i in sorted(order[:2])]

    def _finish_mixture_figure(self, fig: Any, filename: str) -> None:
        fig.savefig(
            os.path.join(self.analysis_dir, filename + cnst.CLUSTERING_PLOT_FILEEXT),
            dpi=300, bbox_inches='tight', pad_inches=0.1,
        )
        if not self.plot_cfg.show_plots:
            plt.close(fig)

    @staticmethod
    def _ternary_xy(A: 'np.ndarray') -> 'np.ndarray':
        """Fractions of (a, b, c) -> 2D coordinates: a bottom-left, b bottom-right, c top."""
        A = np.atleast_2d(np.asarray(A, dtype=float))
        A = A / A.sum(axis=1, keepdims=True)
        return np.c_[A[:, 1] + 0.5 * A[:, 2], np.sqrt(3) / 2 * A[:, 2]]

    @staticmethod
    def _draw_ternary_grid(ax: Any, step: float, edge_color: str = 'k') -> None:
        C = PlottingModule._ternary_xy(np.eye(3))
        ax.plot(*np.vstack([C, C[:1]]).T, '-', color=edge_color, lw=1.5 if edge_color != 'k' else 1, zorder=1)
        for f in np.arange(step, 1.0 - 1e-9, step):
            for i in range(3):
                j, k = (i + 1) % 3, (i + 2) % 3
                a = np.zeros(3); a[i] = f; a[j] = 1 - f
                b = np.zeros(3); b[i] = f; b[k] = 1 - f
                ax.plot(*PlottingModule._ternary_xy(np.array([a, b])).T, '-', color='#DDDDDD', lw=0.7, zorder=0)

    # Label styles of the 2D plots: small, on a light background, placed to avoid overlaps
    LABEL_BOX = dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.6)

    @staticmethod
    def _adjust_labels(ax: Any, texts: List[Any], avoid_points: Optional['np.ndarray'] = None,
                       avoid_artists: Optional[List[Any]] = None) -> None:
        """
        Move labels so they overlap neither each other, the plotted points, nor fixed artists such as the
        corner labels of a ternary diagram (adjustText), linking moved labels to their marker with a thin line. Axis limits must be set before calling. Labels keep their
        initial position if adjustText is not installed.
        """
        if not texts:
            return
        try:
            from adjustText import adjust_text
        except ImportError:
            return
        kwargs = dict(ax=ax, expand=(1.15, 1.3), ensure_inside_axes=True,
                      arrowprops=dict(arrowstyle='-', color='#777777', lw=0.6))
        pts = [np.atleast_2d(avoid_points)] if avoid_points is not None and len(avoid_points) else []
        if avoid_artists:
            # Fixed artists (e.g. corner labels): obstacle points covering their extent on the canvas
            renderer = ax.figure.canvas.get_renderer()
            to_data = ax.transData.inverted()
            for artist in avoid_artists:
                bb = artist.get_window_extent(renderer)
                gx, gy = np.meshgrid(np.linspace(bb.x0, bb.x1, 7), np.linspace(bb.y0, bb.y1, 3))
                pts.append(to_data.transform(np.c_[gx.ravel(), gy.ravel()]))
        if pts:
            pts = np.vstack(pts)
            kwargs.update(x=pts[:, 0], y=pts[:, 1])
        try:
            adjust_text(texts, **kwargs)
        except Exception as exc:  # label placement must never break plotting
            logger.debug(f"Label placement failed: {exc}")

    @staticmethod
    @staticmethod
    def _label_start(q: 'np.ndarray', corners: Optional['np.ndarray']) -> 'np.ndarray':
        """Initial label position: at the marker, or moved inwards for a marker on a ternary corner."""
        if corners is not None and np.min(np.linalg.norm(corners - q, axis=1)) < 0.05:
            return q + 0.06 * (corners.mean(axis=0) - q) / np.linalg.norm(corners.mean(axis=0) - q)
        return q

    @staticmethod
    def _label_mixture_phases(ax: Any, Q: 'np.ndarray', phases: List[str], fontsize: int,
                              corners: Optional['np.ndarray'] = None) -> List[Any]:
        """Mixture phases as blue stars, with labels (to be placed by ``_adjust_labels``)."""
        texts = []
        for k, (f, q) in enumerate(zip(phases, Q)):
            ax.scatter(*q, color='blue', marker='*', s=160, zorder=4, label='Mixture phases' if k == 0 else None)
            texts.append(ax.text(*PlottingModule._label_start(q, corners), to_latex_formula(f), fontsize=fontsize - 2,
                                 ha='left', va='bottom', bbox=PlottingModule.LABEL_BOX, zorder=3))
        return texts

    @staticmethod
    def _plot_other_refs(ax: Any, points: dict, view: tuple, fontsize: int) -> List[Any]:
        """Other candidate phases in view, as grey stars, with labels (to be placed by ``_adjust_labels``)."""
        (x0, x1), (y0, y1) = view
        texts = []
        for f, a in points.items():
            if not (x0 <= a[0] <= x1 and y0 <= a[1] <= y1):
                continue
            ax.scatter(*a, color='#9A9A9A', marker='*', s=80, zorder=4, label='Other candidate phases' if not texts else None)
            texts.append(ax.text(*a, to_latex_formula(f), color='#777777', fontsize=fontsize - 4, ha='left', va='top',
                                 bbox=PlottingModule.LABEL_BOX, zorder=3))
        return texts

    @staticmethod
    def _plot_reconstruction(ax: Any, P: 'np.ndarray', R: 'np.ndarray') -> None:
        for p, r in zip(P, R):
            ax.plot([p[0], r[0]], [p[1], r[1]], '-', color='grey', lw=0.5, alpha=0.6, zorder=1)
        ax.scatter(R[:, 0], R[:, 1], color='red', marker='x', s=15, linewidths=1, label='Position on mixture line', zorder=3)

    @staticmethod
    def _square_view(points: 'np.ndarray', margin: float = 0.15) -> tuple:
        lo, hi = points.min(axis=0), points.max(axis=0)
        half = (hi - lo).max() / 2 * (1 + margin) + 1e-3
        mid = (lo + hi) / 2
        return (mid[0] - half, mid[0] + half), (mid[1] - half, mid[1] + half)

    @staticmethod
    def _plot_mixture_element_ternary(X, W, H, phases, other_refs, els, unit, title, info, zoomed):
        """Cases with 3 elements: points, mixture phases (line or triangle), other candidates, in the element ternary."""
        fontsize = PlottingModule._mixture_style()
        fig, ax = plt.subplots(figsize=(6, 6.6))
        txy = PlottingModule._ternary_xy
        P, Q = txy(X), txy(H)
        corners = txy(np.eye(3))
        view = PlottingModule._square_view(np.vstack([P, Q])) if zoomed else ((-0.12, 1.12), (-0.12, 0.98))
        PlottingModule._draw_ternary_grid(ax, 0.05 if zoomed else 0.1)
        corner_labels = []
        if not zoomed:
            offsets = [(-0.03, -0.06), (0.03, -0.06), (0, 0.05)]
            for (x, y), el, (dx, dy) in zip(corners, els, offsets):
                corner_labels.append(ax.text(x + dx, y + dy, f'{el} ({unit})', ha='center', va='center', fontsize=fontsize))
        if len(phases) == 2:
            PlottingModule._plot_reconstruction(ax, P, txy(W @ H))
        ax.scatter(P[:, 0], P[:, 1], color=plt.get_cmap('viridis')(0.0), marker='o', s=25, alpha=0.8, label='Compositions', zorder=2)
        ax.plot(*(np.vstack([Q, Q[:1]]) if len(Q) > 2 else Q).T, '-', color='blue', lw=1.5, zorder=2)
        texts = PlottingModule._plot_other_refs(ax, {f: txy(v)[0] for f, v in other_refs.items()}, view, fontsize)
        texts += PlottingModule._label_mixture_phases(ax, Q, phases, fontsize, corners)
        ax.set_xlim(*view[0]); ax.set_ylim(*view[1])
        ax.set_aspect('equal'); ax.axis('off')
        PlottingModule._adjust_labels(ax, texts, np.vstack([P, Q]), corner_labels)
        ax.set_title(title + (' (zoomed)' if zoomed else ''))
        ax.legend(fontsize=fontsize - 3, loc='upper center', bbox_to_anchor=(0.5, -0.01), ncol=2, frameon=False)
        fig.text(0.5, 0.02, info, ha='center', va='bottom', fontsize=fontsize - 2)
        fig.subplots_adjust(bottom=0.16)
        return fig

    @staticmethod
    def _plot_mixture_element_2d(X, R, H, other_refs, phases, els_xy, unit, title, info, zoomed):
        """2 phases, not 3 elements: points, mixture line and other candidates in the plane of the two most representative elements."""
        fontsize = PlottingModule._mixture_style()
        fig, ax = plt.subplots(figsize=(6, 6.6))
        P, Rp, Q = X * 100, R * 100, H * 100
        if zoomed:
            # Fractions are non-negative: keep the zoomed view within the positive quadrant
            (x0, x1), (y0, y1) = PlottingModule._square_view(np.vstack([P, Q]))
            view = ((max(x0, 0.0) - 1, x1 - min(x0, 0.0)), (max(y0, 0.0) - 1, y1 - min(y0, 0.0)))
        else:
            view = ((0, 100), (0, 100))
        PlottingModule._plot_reconstruction(ax, P, Rp)
        ax.scatter(P[:, 0], P[:, 1], color=plt.get_cmap('viridis')(0.0), marker='o', s=25, alpha=0.8, label='Compositions', zorder=2)
        ax.plot(*Q.T, '-', color='blue', lw=1.5, zorder=2)
        texts = PlottingModule._plot_other_refs(ax, {f: v * 100 for f, v in other_refs.items()}, view, fontsize)
        texts += PlottingModule._label_mixture_phases(ax, Q, phases, fontsize)
        ax.set_xlim(*view[0]); ax.set_ylim(*view[1]); ax.set_aspect('equal')
        PlottingModule._adjust_labels(ax, texts, np.vstack([P, Q]))
        ax.set_xlabel(f'{els_xy[0]} ({unit})'); ax.set_ylabel(f'{els_xy[1]} ({unit})')
        ax.set_title(title + (' (zoomed)' if zoomed else ''))
        ax.legend(fontsize=fontsize - 3, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
        fig.text(0.5, 0.01, info, ha='center', va='bottom', fontsize=fontsize - 2)
        fig.subplots_adjust(bottom=0.22)
        return fig

    @staticmethod
    def _plot_mixture_phase_ternary(mol_frs, dist, phases, unit, title, info):
        """3 phases, 4+ elements: molar fractions of the three phases, points coloured by their distance from the mixture plane."""
        fontsize = PlottingModule._mixture_style()
        fig, ax = plt.subplots(figsize=(6.6, 6))
        corners = PlottingModule._ternary_xy(np.eye(3))
        PlottingModule._draw_ternary_grid(ax, 0.1, edge_color='blue')
        P = PlottingModule._ternary_xy(mol_frs)
        sc = ax.scatter(P[:, 0], P[:, 1], c=dist, cmap='viridis', marker='o', s=30, alpha=0.9, zorder=2)
        offsets = [(-0.02, -0.07), (0.02, -0.07), (0, 0.06)]
        for (x, y), f, (dx, dy) in zip(corners, phases, offsets):
            ax.scatter(x, y, color='blue', marker='*', s=160, zorder=4)
            ax.text(x + dx, y + dy, to_latex_formula(f), ha='center', va='center', fontsize=fontsize)
        cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label(f'Distance from mixture plane ({unit})')
        ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.15, 1.0)
        ax.set_aspect('equal'); ax.axis('off')
        ax.set_title(title + '\n(molar fractions of the mixture phases)')
        fig.text(0.45, 0.01, info, ha='center', va='bottom', fontsize=fontsize - 4)
        return fig
