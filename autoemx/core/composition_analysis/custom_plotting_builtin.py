import os

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

import autoemx.utils.constants as cnst
from autoemx.core.composition_analysis.clustering_plot_axes import (
    CLUSTERING_3D_VIEW_AZIM,
    CLUSTERING_3D_VIEW_ELEV,
    apply_data_driven_axis_limits,
    gather_clustering_zoom_points,
)
from autoemx.utils.helper import to_latex_formula

custom_dir = ''

#%% Customize 3D clustering plot
def _save_clustering_plot_custom_3D(elements, els_comps_list, centroids, labels, els_std_dev_per_cluster,
                                    unused_compositions_list,
                                    clustering_features,
                                    ref_phases_df,
                                    ref_formulae,
                                    show_plots,
                                    sample_ID,
                                    analysis_dir=None,
                                    output_filename=None,
                                    ):

    plot_file_title = output_filename or cnst.CUSTOM_CLUSTERING_PLOT_FILENAME + cnst.CLUSTERING_PLOT_FILEEXT

    plt.rcParams['font.family'] = 'Arial'
    fontsize = 14
    labelpad = 12
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize

    axis_label_add = ' (w%)' if str(clustering_features).startswith('w') else ' (at%)'
    is_3d = len(elements) == 3

    fig = plt.figure(figsize=(6, 6))
    if is_3d:
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)

    ax.scatter(*els_comps_list, c=labels, cmap='viridis', marker='o', label='Measured comps.')
    ax.scatter(*centroids.T, c='red', marker='x', s=100, label='Centroids')

    if unused_compositions_list:
        ax.scatter(*np.array(unused_compositions_list).T, c='grey', marker='^', label='Discarded comps.')

    first_ellipse = True
    for centroid, stdevs in zip(centroids, els_std_dev_per_cluster):
        if np.any(np.isnan(stdevs)):
            continue

        if len(elements) == 3: # 3D plot
            x_c, y_c, z_c = centroid
            rx, ry, rz = stdevs

            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = x_c + rx * np.outer(np.cos(u), np.sin(v))
            y = y_c + ry * np.outer(np.sin(u), np.sin(v))
            z = z_c + rz * np.outer(np.ones_like(u), np.cos(v))

            ax.plot_surface(x, y, z, color='red', alpha=0.1, edgecolor='none')
            if first_ellipse:
                first_ellipse= False
                ax.plot([], [], [], color='red', alpha=0.1, label='Stddev')
        else: # 2D plot
            x_c, y_c = centroid
            rx, ry = stdevs

            ellipse = patches.Ellipse((x_c, y_c), rx, ry, edgecolor='red', facecolor='red', linestyle='--', alpha=0.2)
            if first_ellipse:
                ellipse.set_label('Stddev')
                first_ellipse = False
            ax.add_patch(ellipse)

    if ref_phases_df is not None:
        first_ref = True
        ref_phases_df = ref_phases_df[elements] if ref_phases_df is not None else None
        for index, row in ref_phases_df.iterrows():
            label = 'Candidate phases' if first_ref else None
            ref_formula = ref_formulae[index]
            ax.scatter(*row.values, c='blue', marker='*', s=100, label=label)

            if len(elements) == 3:
                dx = 0.05
                x_label, y_label, z_label = row.values
                ax.text(x_label + dx, y_label + dx, z_label + dx,
                        to_latex_formula(ref_formula), color='black', fontsize=fontsize, ha='left', va='bottom')
            else:
                dx = 0.002
                x_label, y_label = row.values
                ax.text(x_label + dx, y_label + dx,
                        to_latex_formula(ref_formula), color='black', fontsize=fontsize, ha='left', va='bottom')
            first_ref = False

    for i, centroid in enumerate(centroids):
        ax.text(*centroid, str(i), color='black', fontsize=fontsize, ha='right', va='bottom')

    ax.set_xlabel(elements[0] + axis_label_add, labelpad=labelpad)
    ax.set_ylabel(elements[1] + axis_label_add, labelpad=labelpad)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if is_3d:
        ax.set_zlabel(elements[2] + axis_label_add, labelpad=labelpad)
        ax.set_zlim(0, 1)
        ax.view_init(elev=CLUSTERING_3D_VIEW_ELEV, azim=CLUSTERING_3D_VIEW_AZIM)

    all_points = gather_clustering_zoom_points(
        els_comps_list,
        centroids,
        unused_compositions_list,
        elements,
        ref_phases_df=ref_phases_df,
    )

    apply_data_driven_axis_limits(ax, all_points, is_3d=is_3d)

    ax.set_title(f'Custom clustering {sample_ID}')
    ax.legend(fontsize=fontsize)

    output_dir = analysis_dir if analysis_dir else custom_dir
    fig.savefig(os.path.join(output_dir, plot_file_title), dpi=300, bbox_inches='tight', pad_inches=0.1)

    if show_plots:
        plt.ion()
        plt.show()
        plt.pause(0.001)

    if not show_plots:
        plt.close(fig)
