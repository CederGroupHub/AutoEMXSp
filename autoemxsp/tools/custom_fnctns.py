import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import numpy as np

from autoemxsp.tools.utils import to_latex_formula

custom_dir = ''

#%% Edit plot colors and fonts
c_SEMEDS_paper = "#519fa7"
plt.rcParams['font.family'] = 'Arial'
paper_fontsize = 15

# Function to truncate colormap
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    return mcolors.LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )

# Create truncated viridis (e.g. darker half)
custom_cmap = truncate_colormap(plt.cm.viridis, 0.4, 0.9)

# Okabe-Ito color palette (8-color version, colorblind-friendly)
okabe_ito_colors = [
    "#009E73",  # green
    "#F0E442",  # yellow
    "#56B4E9",  # sky blue
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#000000",  # black
]
# Create a ListedColormap using the first few (e.g., 5) colors
custom_cmap = mcolors.ListedColormap(okabe_ito_colors[:5])

#%% Customize 3D clustering plot
def _save_clustering_plot_custom_3D(elements, els_comps_list, centroids, labels, els_std_dev_per_cluster,
                                    unused_compositions_list,
                                    clustering_features,
                                    ref_phases_df,
                                    ref_formulae,
                                    show_plots,
                                    sample_ID
                                    ):
    
    plot_file_title = f'{sample_ID}_3Dclustering.pdf'
    
    # Plot options
    plt.rcParams['font.family'] = 'Arial'
    fontsize = paper_fontsize
    labelpad = 15
    plt.rcParams['font.size'] = fontsize         # General font size (default 10)
    plt.rcParams['axes.titlesize'] = fontsize    # Title font size
    plt.rcParams['axes.labelsize'] = fontsize    # Axis label font size
    plt.rcParams['xtick.labelsize'] = fontsize   # X-axis tick label font size
    plt.rcParams['ytick.labelsize'] = fontsize   # Y-axis tick label font size
    
    # x_NaP_vals = np.array(els_comps_list[0])/(2* 0.308) + np.array(els_comps_list[1])/(2* 0.154)
    # x_Sn_vals = np.array(els_comps_list[2])/0.333
    # #Select points that fit within a certain analytical error!
    # valid_pts = [i for i in range(len(x_NaP_vals)) if np.abs(x_NaP_vals[i] + x_Sn_vals[i] - 1) <0.5]
    
    # Define labels
    if clustering_features[0] == 'w':
        axis_label_add = ' (w%)'
    else:
        axis_label_add = ' (at%)'
    
    ### Figure
    # fig = plt.figure(figsize=(8, 8))
    fig = plt.figure(figsize=(5, 5))
    # Add 3rd dimension if 3D plot
    if len(elements) == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.set_zlabel(elements[2] + axis_label_add, labelpad = labelpad)
        ax.set_zlim(0, 1)
        
        # Clean background and axes
        ax.set_facecolor('white')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
    else:
        ax = fig.add_subplot(111)
    
    
    
    
    # x_vals = [els_comps_list[0]]#[i] for i in valid_pts]
    # y_vals = [els_comps_list[1]]#[i] for i in valid_pts]
    # z_vals = [els_comps_list[2]]#[i] for i in valid_pts]
    # labels = [labels[i] for i in valid_pts]
    # #Select points that fit within a certain analytical error!
    # ax.scatter(*els_comps_list, c=labels, cmap=custom_cmap, marker='o')

    
    # Plot datapoints used for clustering, coloring them based on the cluster
    # ax.scatter(x_vals,y_vals,z_vals, c=labels, s=50, cmap=custom_cmap, marker='o')
    ax.scatter(*els_comps_list, c=labels, cmap=custom_cmap, s=30, marker='o', label = 'Measured comps.')
    # ax.scatter(*els_comps_list, c= c_SEMEDS_paper, s=40, marker='o')
    
    
    # Add datapoints that were not used for clustering
    if unused_compositions_list != []:
        ax.scatter(*np.array(unused_compositions_list).T, c='black', marker='^', label = 'Discarded comps.')

    
    # Plot centroids
    ax.scatter(*centroids.T, c = 'red', marker='x', s=100, label= 'Centroid')
    
    # Plot standard deviation of centroids
    first_ellipse = True
    for centroid, stdevs in zip(centroids, els_std_dev_per_cluster):
        if len(elements) == 3: # 3D plot
            x_c, y_c, z_c = centroid
            rx, ry, rz = stdevs
            
            # Create the ellipsoid
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)
            x = x_c + rx * np.outer(np.cos(u), np.sin(v))
            y = y_c + ry * np.outer(np.sin(u), np.sin(v))
            z = z_c + rz * np.outer(np.ones_like(u), np.cos(v))

            # Plot the surface with transparency
            ax.plot_surface(x, y, z, color='red', alpha=0.06, edgecolor='none')
            
            # Add legend entry
            if first_ellipse:
                first_ellipse= False
                # Plot a dummy point with the label
                ax.plot([], [], [], color='red', alpha=0.06, label='Std. dev.')
        
            
        else: # 2D plot
            x_c, y_c = centroid
            rx, ry = stdevs
            
            # Plot the ellipse with transparency
            ellipse = patches.Ellipse((x_c, y_c), rx, ry, edgecolor='red', facecolor='red', linestyle='--', alpha=0.2)
            
            # Add label
            if first_ellipse:
                ellipse.set_label('Stddev')
                first_ellipse = False
            
            # Add ellipse to plot
            ax.add_patch(ellipse)
    
    
    # Add reference phases
    if ref_phases_df is not None:
        first_ref = True
        ref_phases_df = ref_phases_df[elements]
        for index, row in ref_phases_df.iterrows():
            if first_ref:
                first_ref = False
                label = 'Reference phases'
            else:
                label = None
            ref_formula = ref_formulae[index]
            ax.scatter(*row.values, c = 'blue', marker='*', s=100, label= label)
            dx = 0.05  # fixed horizontal offset in data units
            
            if len(elements) == 3:
                x_label, y_label, z_label = row.values
                ax.text(x_label + dx, y_label + dx, z_label + dx,
                        to_latex_formula(ref_formula), color='black', fontsize=fontsize, ha='left', va='bottom')
            else:
                dx = 0.002
                x_label, y_label = row.values
                ax.text(x_label + dx, y_label + dx,
                        to_latex_formula(ref_formula), color='black', fontsize=fontsize, ha='left', va='bottom')

    
    ax.set_xlabel(elements[0] + axis_label_add, labelpad=labelpad)
    ax.set_ylabel(elements[1] + axis_label_add, labelpad=labelpad)
    
    ### Customisation
    lims_x = (0.05, 0.2501)
    lims_y = (0, 0.2001)
    lims_z = (0.05, 0.2501)
    lower_lim_x, upper_lim_x = lims_x
    lower_lim_y, upper_lim_y = lims_y
    lower_lim_z, upper_lim_z = lims_z

    delta_grid = 0.1

    
        
    
    ax.set_xlim(lower_lim_x, upper_lim_x)
    ax.set_ylim(lower_lim_y, upper_lim_y)
    
    if len(elements) == 3:
        ax.set_zlim(lower_lim_z, upper_lim_z)
    
    x_ticks = np.arange(lower_lim_x, upper_lim_x, delta_grid)
    y_ticks = np.arange(lower_lim_y, upper_lim_y, delta_grid)
    if len(elements) == 3: z_ticks = np.arange(lower_lim_z, upper_lim_z, delta_grid)
    
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f"{x*100:.0f}" for x in x_ticks])
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f"{x*100:.0f}" for x in y_ticks])
    if len(elements) == 3:
        ax.set_zticks(z_ticks)
        ax.set_zticklabels([f"{x*100:.0f}" for x in z_ticks])
    
    # ref_phases_df = ref_phases_df[elements]
    # for index, row in ref_phases_df.iterrows():
    #     ref_formula = ref_formulae[index]
    #     ax.scatter(*row.values, c = 'blue', marker='*', s=100, label= 'Reference phases')
        # if ref_formula == 'PbMO2':
        #     ax.scatter(*row.values, c = 'blue', marker='o', s=70, label= 'Precursors')
        #     ax.text(*row.values + np.array([0.01,0,0]), 'SnO$_2$', color='black', fontsize=fontsize, ha='left', va='top')
        # if ref_formula == 'PbMoO4':
        #     ax.scatter(*row.values, c = 'green', marker='*', s=70, label= 'Phase')
            # ax.text(*row.values + np.array([-0.04,0.01,0]), 'PbMoO$_4$', color='black', fontsize=fontsize, ha='right', va='top')
            
    # Plot projection
    if len(elements) == 3:
        ax.view_init(elev=20, azim=40)
    
    # legend = ax.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9), fontsize=fontsize, 
    #   frameon=True, facecolor='white', edgecolor='black')    
    # frame = legend.get_frame()
    # frame.set_alpha(0.9)  # Set the alpha value to 0.9 (you can adjust this value as needed)
        
    # ax.set_title(f'K-Means Clustering {self.sample_ID}')
    # ax.legend(fontsize = fontsize)
    
    
    if show_plots:
        plt.ion()
        plt.show()
        plt.pause(0.001)
    
    # plt.tight_layout()
    
    plt.savefig(os.path.join(custom_dir, plot_file_title), dpi=300, bbox_inches='tight', pad_inches=0.5)
    
        
    # # Extract legend handles and labels
    # handles, labels = ax.get_legend_handles_labels()
    
    # # Create a new figure just for the legend
    # fig_legend = plt.figure(figsize=(3, 2))
    # fig_legend.legend(handles, labels, loc='center')
    # fig_legend.tight_layout()
    
    # # Save just the legend
    # fig_legend.savefig(os.path.join(custom_dir, "Clustering_plot_legend.pdf"), dpi=300, transparent=True)
    
    # Close the figure if plot is False
    if not show_plots:
        plt.close(fig)
