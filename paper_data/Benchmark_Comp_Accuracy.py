#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for plotting histograms of relative and absolute deviations from expected values for composition accuracy evaluation.
Plots and prints statistics.

Select samples to analyse in Sample selection, and play with the Options to visualize different histograms.

Created on Fri Nov 22 14:17:43 2024

@author: Andrea
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from pymatgen.core.composition import Composition
from pymatviz.ptable import ptable_heatmap
from collections import Counter
from scipy.optimize import curve_fit

import autoemxsp
from autoemxsp.tools.utils import print_single_separator, print_double_separator, get_sample_dir
import autoemxsp.tools.constants as cnst

okabeito_green = '#009E73'
okabeito_orange = '#E69F00'
okabeito_blue = '#0072B2'
okabeito_lightblue = '#5bb4e5'
frame_grey = '#cccccc'
text_grey = '#4f4e50'

color_SEMEDS_paper_teal = "#519fa7"
color_SEMEDS_paper_lighter_teal = '#c4dfe1'  # lighter shade of #519fa7
color_SEMEDS_paper_gray = '#4E4E50'
    
added_text = ''
min_fr = 0
max_fr = 100
min_abs_er = 0
min_rdev = 0
n_samples = 0
els_not_to_consider = None
els_to_stack_hist = []
add_fit = False

undetectable_els = ['H', 'He', 'Li', 'Be']

#%% Sample lists
synthetic_samples = [
        # NASICONS
        {'sample_ID' : 'NaGe2(PO4)3', 'formula' : 'NaGe2(PO4)3'},
        {'sample_ID' : 'NaSn2(PO4)3', 'formula' : 'NaSn2(PO4)3'},
        {'sample_ID' : 'Na0.4Zr1.4Ta0.6(PO4)3', 'formula' : 'Na0.4Zr1.4Ta0.6(PO4)3', 'clusters_to_ignore' : [1]}, # half-reacted ZrO2 leftover
        {'sample_ID' : 'NaZrTi(PO4)3', 'formula' : 'NaZrTi(PO4)3', 'clusters_to_ignore' : [1]}, # TiO2 leftover
        {'sample_ID' : 'NaTiSn(PO4)3', 'formula' : 'NaTiSn(PO4)3'},
        
        # Other
        {'sample_ID' : 'MnAgO2', 'formula' : 'MnAgO2', 'clusters_to_ignore' : [1]}, # Mn-O impurity, forced k=2
        {'sample_ID' : 'CaCo(PO3)4', 'formula' : 'CaCo(PO3)4', 'clusters_to_ignore' : [1]}, # Ca(PO3)2 impurity
        {'sample_ID' : 'MgTi2NiO6', 'formula' : 'MgTi2NiO6', 'clusters_to_ignore' : [1]}, # Some half-reacted mix of precursors left
        {'sample_ID' : 'K4MgFe3(PO4)5', 'formula' : 'K4MgFe3(PO4)5', 'clusters_to_ignore' : [1]}, # Some unreacted mix of precursors left
        {'sample_ID' : 'KNaTi2(PO5)2', 'formula' : 'KNaTi2(PO5)2', 'clusters_to_ignore' : [1]}, # TiO2 impurity
        {'sample_ID' : 'Hf2Sb2Pb4O13', 'formula' : 'Hf2Sb2Pb4O13', 'clusters_to_ignore' : [1]}, # Hf-rich impurity leftover as expected from XRD
        {'sample_ID' : 'K2TiCr(PO4)3', 'formula' : 'K2TiCr(PO4)3', 'clusters_to_ignore' : [1]}, # Cr2O3 leftover
        {'sample_ID' : 'MgCuP2O7', 'formula' : 'MgCuP2O7'},
        {'sample_ID' : 'MgTi4(PO4)6', 'formula' : 'MgTi4(PO4)6', 'clusters_to_ignore' : [1]}, # A few Mg-rich points
        {'sample_ID' : 'NaCaMgFe(SiO3)4', 'formula' : 'NaCaMgFe(SiO3)4', 'clusters_to_ignore' : [1]}, # Some CaSiO3 present
        {'sample_ID' : 'Bi2Fe4O9', 'formula' : 'Bi2Fe4O9'},
        {'sample_ID' : 'Bi25FeO39', 'formula' : 'Bi25FeO39'}, # Clusters forced to be 1 instead of 2. Not necessary if using w_fr as features
        {'sample_ID' : 'LaNbO4', 'formula' : 'LaNbO4', 'clusters_to_ignore' : [1]}
    ]

eds_standards = [
        {'sample_ID' : 'ScPO4_mineral', 'formula' : 'ScPO4'},
        {'sample_ID' : 'YIG_mineral', 'formula' : 'Y3Fe5O12'},
        {'sample_ID' : 'YPO4_mineral', 'formula' : 'YPO4'},
        {'sample_ID' : 'K-412_NISTstd_mineral', 'formula' : 'Fe0.156Mg0.538Ca0.305Al0.204Si0.847O3', 'clusters_to_ignore' : [1]}, # SiO2 point

        {'sample_ID' : 'Alamosite_mineral', 'formula' : 'PbSiO3', 'clusters_to_ignore' : [1]}, # There is clearly a particle of different composition
        {'sample_ID' : 'Albite_mineral', 'formula' : 'NaAlSi3O8'},
        {'sample_ID' : 'Anhydrite_mineral', 'formula' : 'CaSO4'},
        {'sample_ID' : 'Anorthite_mineral', 'formula' : 'CaAl2Si2O8'},
        {'sample_ID' : 'Benitoite_mineral', 'formula' : 'BaTiSi3O9'},
        {'sample_ID' : 'Bornite_mineral', 'formula' : 'Cu5FeS4'},
        {'sample_ID' : 'Chalcopyrite_mineral', 'formula' : 'CuFeS2'},
        {'sample_ID' : 'CoOlivine_mineral', 'formula' : 'Co2SiO4'},
        {'sample_ID' : 'FeOlivine_mineral', 'formula' : 'Fe2SiO4'},
        {'sample_ID' : 'Fluorphlogopite_mineral', 'formula' : 'KMg3(AlSi3O10)F2'},
        {'sample_ID' : 'Jadeite_mineral', 'formula' : 'NaAlSi2O6'},
        {'sample_ID' : 'Labradorite_mineral', 'formula' : '(Ca0.65Na0.34)(Al1.66Si2.33)O8'},
        {'sample_ID' : 'MnOlivine_mineral', 'formula' : 'Mn2SiO4', 'clusters_to_ignore' : [1]}, # Single particle with a different composition, clearly seen from all its spectra
        {'sample_ID' : 'Nepheline_mineral', 'formula' : '(Na0.79K0.17)(Al0.95Si1.04)O4'},
        {'sample_ID' : 'Orthoclase_mineral', 'formula' : '(K0.915Na0.082)(Al0.912Fe0.073)Si3O8'},
        {'sample_ID' : 'Rhodonite_mineral', 'formula' : '(Mn0.589Fe0.196Ca0.178)SiO3'},
        {'sample_ID' : 'Wulfenite_mineral', 'formula' : 'PbMoO4'}
    ]


commercial_precursors = [
            {'sample_ID' : 'Al2O3_precursor', 'formula' : 'Al2O3'},
            {'sample_ID' : 'AlPO4_precursor', 'formula' : 'AlPO4'},
            {'sample_ID' : 'Ba3(PO4)2_precursor', 'formula' : 'Ba3(PO4)2'},
            {'sample_ID' : 'Co3O4_precursor', 'formula' : 'Co3O4'},
            {'sample_ID' : 'CuO_precursor', 'formula' : 'CuO'},
            {'sample_ID' : 'Fe2O3_precursor', 'formula' : 'Fe2O3'},
            {'sample_ID' : 'Ga2O3_precursor', 'formula' : 'Ga2O3'},
            {'sample_ID' : 'GeO2_precursor', 'formula' : 'GeO2', 'clusters_to_ignore' : [1]}, # Extra point, automatically excluded 
            {'sample_ID' : 'HfO2_precursor', 'formula' : 'HfO2'},
            {'sample_ID' : 'In2O3_precursor', 'formula' : 'In2O3'},
            {'sample_ID' : 'KCl_precursor', 'formula' : 'KCl', 'clusters_to_ignore' : [1]}, # Extra point, automatically excluded
            {'sample_ID' : 'Li2WO4_precursor', 'formula' : 'Li2WO4'}, # Forced k =1
            {'sample_ID' : 'LiCoPO4_precursor', 'formula' : 'LiCoPO4'},
            {'sample_ID' : 'LiNiCoMnO2_precursor', 'formula' : 'Li3NiCoMnO6'},
            {'sample_ID' : 'MgF2_precursor', 'formula' : 'MgF2'},
            {'sample_ID' : 'MgO_precursor', 'formula' : 'MgO'},
            {'sample_ID' : 'MnO_precursor', 'formula' : 'MnO'},
            {'sample_ID' : 'MnO2_precursor', 'formula' : 'MnO2'},
            {'sample_ID' : 'Mn2O3_precursor', 'formula' : 'Mn2O3'},
            {'sample_ID' : 'MoO2_precursor', 'formula' : 'MoO2'},
            {'sample_ID' : 'Na2MoO4_precursor', 'formula' : 'Na2MoO4'},
            {'sample_ID' : 'Na4P2O7_precursor', 'formula' : 'Na4P2O7'}, # Hygroscopic?
            {'sample_ID' : 'NaNO3_precursor', 'formula' : 'NaNO3'},
            {'sample_ID' : 'Ni(OH)2_precursor', 'formula' : 'Ni(OH)2', 'clusters_to_ignore' : [1]}, # Some NiO present
            {'sample_ID' : 'NiO_precursor', 'formula' : 'NiO'},
            {'sample_ID' : 'PbO_precursor', 'formula' : 'PbO'},
            {'sample_ID' : 'Sb2O3_precursor', 'formula' : 'Sb2O3'},
            {'sample_ID' : 'SiO2_precursor', 'formula' : 'SiO2'},
            {'sample_ID' : 'SnO2_precursor', 'formula' : 'SnO2'},
            {'sample_ID' : 'Ta2O5_precursor', 'formula' : 'Ta2O5'},
            {'sample_ID' : 'TiN_precursor', 'formula' : 'TiN'},
            {'sample_ID' : 'TiO2_precursor', 'formula' : 'TiO2'},
            {'sample_ID' : 'WO3_precursor', 'formula' : 'WO3'},
            {'sample_ID' : 'ZnO_precursor', 'formula' : 'ZnO'},
            {'sample_ID' : 'ZrO2_precursor', 'formula' : 'Zr0.985O2'}, #Removed measured mass fraction of Hf from reference formula
]

def get_n_els_in_formula(formula):
    # Automatically remove undetectable elements from count
    return len(set(Composition(formula).as_dict().keys()) - set(undetectable_els))

all_samples = commercial_precursors + synthetic_samples + eds_standards
binary_samples = [s for s in all_samples if get_n_els_in_formula(s['formula']) == 2]
ternary_samples = [s for s in all_samples if get_n_els_in_formula(s['formula']) == 3]
quaternary_plus_samples = [s for s in all_samples if get_n_els_in_formula(s['formula']) >= 4]


#%% Sample selection

samples = all_samples
# samples = synthetic_samples
# samples = eds_standards
# samples = commercial_precursors
# samples = commercial_precursors + eds_standards
# samples = binary_samples
# samples = ternary_samples
# samples = quaternary_plus_samples

# Get list of elements
all_elements = set()
for sample in samples:
    comp = Composition(sample['formula'])
    for el in comp.elements:
        all_elements.add(el.symbol)
#print(element_symbols)

#%% Options
results_dir = os.path.dirname(os.path.abspath(__file__))
save_dir = 'Quant Analysis'

comp_type = 'at_fr'
max_analytical_errors = [5] # Plots histograms of all spectra with different max analytical error values
# Does not affect the plotting of clustering values, which need to be re-analysed separately, and saved with a different name, e.g. "For Histogram 50an" for a max 50% analytical error

hist_bin = 1 # Binning of histogram counts

n_pts_thresh = 5 # Checks if all clusters have at least n_pts_thresh points
f_string_hist = "For Histogram 5an" # String defining Analysis folder to use in the sample
added_text += '_5an'

### Perform analysis on a limited number of elements
# -- Select a few elements
# selected_elements_to_consider = ['P']
# els_not_to_consider = all_elements - set(selected_elements_to_consider)
# added_text += f'_Only{selected_elements_to_consider}'

# -- Exclude a few elements
# els_not_to_consider = ['N','O','F']
# added_text += f'_No{els_not_to_consider}'


### Elements highlighted in lighter color in histogram
els_to_stack_hist = ['N','O','F']


### Plot only elemental fraction larger or smaller than a specific value
# min_fr = 10
# added_text += f'_LargerThan{min_fr}'

# max_fr = 10
# added_text += f'_SmallerThan{max_fr}'

### Plot only ADEV larger than a specific value
# min_abs_er = 0.01
# added_text += f'_minadev{min_abs_er}'

### Plot only RDEV larger than a specific value
# min_rdev = 10
# added_text += f'_minrdev{min_rdev}'


### Add string to file in function of set of selected samples
# added_text += '_AllSamples'
# added_text += '_OnlyMinerals&Prec'
# added_text += '_OnlySynthetic'
# added_text += '_OnlyEDSstds'
# added_text += '_OnlyCommercialPrecursors'

# added_text += '_OnlyBinaryComps'
# added_text += '_OnlyTernaryComps'
# added_text += '_MoreThan3ElsComps'


# Define quant_flags to consider in the analysis
acceptable_quant_flags = [-1, 0] # Only valid compositions

# acceptable_quant_flags = [-1, 0 ,4, 5, 6,7,8] # All quantified compositions
# added_text += '_WithAllBadFlags'

# Set discarded_spectra_no_aner_considered=True to plot spectra discarded due to filters other than analytical error
discarded_spectra_no_aner_considered = False
# added_text += '_OnlyBadFlags'

# Set all_discarded_spectra=True to plot all discarded spectra
all_discarded_spectra = False
# added_text += '_AllDiscardedSpectra'

### Plot spectra discarded due to specific quant_flags 
# acceptable_quant_flags= [8]
# added_text += f'_OnlyFlag{acceptable_quant_flags[0]}'


### Add Lorentian fit to the histogram
# add_fit = True
# added_text += '_WithFit'


#%% Functions
# Define the Lorentzian function
def lorentzian(x, A, x0, gamma):
    return A / (1 + ((x - x0) / gamma) ** 2)

# # Define the Gaussian function
# def gaussian(x, A, x0, sigma):
#     return A * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))

def fit_histogram(n, bins):
   if isinstance(n, np.ndarray) and n.ndim == 2:
       n = np.sum(n, axis=0)  # Collapse into a single histogram
    
    # Fit the histogram with both Lorentzian and Gaussian
   bin_centers = (bins[:-1] + bins[1:]) / 2  # Compute bin centers

   # Initial guesses for fitting
   initial_guess_lorentzian = [1, 0, 1]  # A, x0, gamma
   # initial_guess_gaussian = [1, 0, 1]    # A, x0, sigma
   
   # Fit the histogram with the Lorentzian
   bounds = (0, [np.inf, np.inf, np.inf])
   popt_lorentzian, _ = curve_fit(lorentzian, bin_centers, n, p0=initial_guess_lorentzian, bounds=bounds)
   A_lor, x0_lor, gamma_lor = popt_lorentzian
   fwhm_lor = 2 * gamma_lor  # FWHM for Lorentzian
   
   # # Fit the histogram with the Gaussian
   # popt_gaussian, _ = curve_fit(gaussian, bin_centers, n, p0=initial_guess_gaussian)
   # A_gauss, x0_gauss, sigma_gauss = popt_gaussian
   # fwhm_gauss = 2.355 * sigma_gauss  # FWHM for Gaussian

   # Plot the fits
   x_fit = np.linspace(bins[0], bins[-1], 500)
   plt.plot(x_fit, lorentzian(x_fit, *popt_lorentzian), 'c--', label=f"Lorentzian Fit\nFWHM = {fwhm_lor:.2f}")
   # plt.plot(x_fit, gaussian(x_fit, *popt_gaussian), 'g--', label=f"Gaussian Fit\nFWHM = {fwhm_gauss:.2f}")
   
   print(f"Lorentzian Parameters: A = {A_lor:.2f}, x0 = {x0_lor:.2f}, gamma = {gamma_lor:.2f}")
   print(f"Lorentzian FWHM = {fwhm_lor:.2f}")
   # print(f"Gaussian Parameters: A = {A_gauss:.2f}, x0 = {x0_gauss:.2f}, sigma = {sigma_gauss:.2f}")
   # print(f"Gaussian FWHM = {fwhm_gauss:.2f}")
   
   
def make_histogram(data, is_cluster_plot, is_rdev_plot, is_abs_er_plot, is_stdev_plot = False, is_drms_plot = False,
                   rdevs_to_stack = [],
                   adevs_to_stack = [], cluster_rstdevs_to_stack = [], cluster_astdevs_to_stack = [],
                   n_undiscarded_spectra = 0, plot_title = 'Histogram', file_name = 'hist', add_fit = False):
    all_data = np.array(data)

    # Determine the range of the data
    min_val = np.min(all_data)
    max_val = np.max(all_data)
    
            
    if comp_type == 'at_fr':
        at_percent_label = '($at$%)'
    elif comp_type == 'w_fr':
        at_percent_label = '($w$%)'
    

    # Extend the range to be symmetric around 0
    range_limit = max(abs(min_val), abs(max_val))
    
    if (is_stdev_plot and is_abs_er_plot) or is_drms_plot:
        plot_hist_bin = hist_bin / 4
    elif is_abs_er_plot:
        plot_hist_bin = hist_bin / 2
    else:
        plot_hist_bin = hist_bin
    
    # Create bins symmetrically around 0 with a step size of 0.01
    bins = np.arange(-range_limit, range_limit + plot_hist_bin, plot_hist_bin)
    # bins = bins - (bins[0] % plot_hist_bin)  # Adjust to align bins with zero
    
    # Round range_limit up to nearest multiple of bin width
    range_limit = np.ceil(range_limit / plot_hist_bin) * plot_hist_bin
    
    # Create bins centered on 0, aligned to integer multiples of plot_hist_bin
    bins = np.arange(-range_limit, range_limit + plot_hist_bin, plot_hist_bin)
    
    # Ensure the rightmost bin includes max value
    bins[-1] += 1e-8  # Tiny offset so max(data) is included
    
    # Plot the histogram
    plt.figure(figsize=(5, 4.5))
        
    print_single_separator()

    # Remove duplicate entries in data, and stack values of low-Z elements
    tolerance = 1e-5
    bins_linewidth = 0.3
    bins_edgecolor = text_grey
    
    plt.grid(axis='y', color=frame_grey, linewidth = bins_linewidth)
    
    
    ### Plot histogram
    def stack_and_plot(data, stack, label=None):
        data = np.array(data)
        stack = np.array(stack)
        mask = np.zeros(len(data), dtype=bool)
        for val in stack:
            idx = np.where((np.abs(data - val) < tolerance) & ~mask)[0]
            if len(idx): mask[idx[0]] = True
        data_filtered = data[~mask]
        if len(els_to_stack_hist) > 0 and label:
            print(f"{np.sum(mask)} fractions from low-Z elements")
        return plt.hist(
            [stack, data_filtered], bins=bins, alpha=1,
            color=[color_SEMEDS_paper_lighter_teal, color_SEMEDS_paper_teal],
            stacked=True, linewidth=bins_linewidth,
            edgecolor=bins_edgecolor, zorder=5,
            label=[", ".join(els_to_stack_hist), "Other els"]  # legend labels
        )
        
    
    data = np.array(data)
    
    if is_rdev_plot and len(rdevs_to_stack)>0:
        n, plotted_bins, _ = stack_and_plot(data, rdevs_to_stack, label="rdev")
    elif is_abs_er_plot and len(adevs_to_stack)>0:
        n, plotted_bins, _ = stack_and_plot(data, adevs_to_stack, label="adev")
    elif is_stdev_plot and len(cluster_astdevs_to_stack)>0:
        n, plotted_bins, _ = stack_and_plot(data, cluster_astdevs_to_stack, label="astdev")
    elif is_stdev_plot and len(cluster_rstdevs_to_stack)>0:
        n, plotted_bins, _ = stack_and_plot(data, cluster_rstdevs_to_stack, label="rstdev")
    else:
        data = data
        if is_drms_plot:
            color = okabeito_lightblue
        else:
            color = color_SEMEDS_paper_teal
        n, plotted_bins, _ = plt.hist(
            data, bins=bins, color=color,
            edgecolor=bins_edgecolor, linewidth=bins_linewidth, zorder=5
        )


    ### Determine band positions
    if is_stdev_plot:
        if is_rdev_plot:
            green_line = 5
            red_line = 10
            # plt.xlabel('$\sigma_{rel}=\sigma_{i}/n_i$ (%)', fontsize = fontsize, color = text_grey)
            plt.xlabel('$\sigma_{i,rel}$ (%)', fontsize = fontsize, color = text_grey)
        elif is_abs_er_plot:
            green_line = 1
            red_line = 2
            plt.xlabel('$\sigma_{i}$ '+at_percent_label, fontsize = fontsize, color = text_grey)
    elif is_drms_plot:
        green_line = 2
        red_line = 3
        percent_label = 'at%' if 'at%' in plot_title else 'w%'
        plt.xlabel('$d_{RMS}$ '+percent_label, fontsize = fontsize, color = text_grey)
    elif is_rdev_plot:
        green_line = 5
        red_line = 10
        plt.xlabel('RDEV (%)', fontsize = fontsize, color = text_grey, fontweight='bold')
    elif is_abs_er_plot:
        green_line = 1
        red_line = 2
        plt.xlabel('ADEV ' + at_percent_label, fontsize = fontsize, color = text_grey, fontweight='bold')
    else:
        plt.xlabel('Error (%)', fontsize = fontsize, color = text_grey)
        
    plt.ylabel('Counts', fontsize = fontsize, color = text_grey)#', binning: {plot_hist_bin}%')
    # Tick marks for both axes (all ticks in frame_grey)
    plt.tick_params(axis='x', which='both', color=frame_grey)
    plt.tick_params(axis='y', which='both', color=frame_grey)
    
    # Major tick labels only in text_grey
    plt.tick_params(axis='x', which='major', labelsize=fontsize, labelcolor=text_grey)
    plt.tick_params(axis='y', which='major', labelsize=fontsize, labelcolor=text_grey)
    ax = plt.gca()
    # ax.set_aspect('equal')  # force square plot area

    
    # Change axes color
    # Alternatively, you can use a loop to change all spines at once
    for spine in ax.spines.values():
        spine.set_color(frame_grey)
        
    
    # Set x range limits
    if not is_cluster_plot:
        if is_rdev_plot:
            # plt.xlim(-100, 100)
            plt.xlim(-30, 30)
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.xaxis.set_major_locator(MultipleLocator(10))
        elif is_abs_er_plot:
            plt.xlim(-20, 20)
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.xaxis.set_major_locator(MultipleLocator(5))

    elif is_stdev_plot:
        if is_rdev_plot:
            plt.xlim(0, 50)
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.xaxis.set_major_locator(MultipleLocator(10))

        elif is_abs_er_plot:
            plt.xlim(0, 5)
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.xaxis.set_major_locator(MultipleLocator(1))
    elif is_drms_plot:
        plt.xlim(0, 8)
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.xaxis.set_major_locator(MultipleLocator(1))
    else:
        if is_rdev_plot:
            plt.xlim(-30, 30)
            # plt.ylim(0, 30)
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.xaxis.set_major_locator(MultipleLocator(10))

        elif is_abs_er_plot:
            plt.xlim(-5, 5)
            # plt.ylim(0, 50)
            ax.xaxis.set_minor_locator(MultipleLocator(1))
            ax.xaxis.set_major_locator(MultipleLocator(2))

    
    # Determine the maximum frequency
    max_height = n.max()
    
    # Conditional tick spacing
    if max_height > 35:
        # Round up to the next multiple of 10
        y_max = int(np.ceil(max_height / 10.0)) * 10 * 1.075
        ax.set_ylim(0, y_max)

        # Choose a tick spacing that gives at most 4 ticks
        for spacing in [10, 20, 25, 30, 40,  50, 75, 100, 200, 300, 400, 500]:
            n_ticks = y_max / spacing + 1  # +1 because ticks include 0

            if 3 <= n_ticks <= 4:
                tick_spacing = spacing
                break
        else:
            tick_spacing = 100  # fallback if all spacings would give too many ticks
        
        ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))
    elif 9 < max_height <=35:
        # Round up to the next multiple of 5
        y_max = int(np.ceil(max_height / 5.0)) * 5
        if y_max == max_height:
            y_max+=1
        ax.set_ylim(0, y_max)
        
        # Choose a tick spacing that gives at most 4 ticks
        for spacing in [3, 5, 2, 10]:
            n_ticks = y_max / spacing + 1  # +1 because ticks include 0

            if 2 <= n_ticks <= 5:
                tick_spacing = spacing
                break
        ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))

    else:
        # Round up to the next multiple of 2
        y_max = 10#int(np.ceil(max_height / 2)) * 2
        ax.set_ylim(0, y_max)
    
        tick_spacing = 2
    
        ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # Use integer ticks
    # Green bands
    # Force autoscale update
    plt.gcf().canvas.draw()
    ax.figure.canvas.draw()
    ymin, ymax = ax.get_ylim()
    collapsed_n = n[1] if n.ndim == 2 else n
    
    for x_pos in [-red_line, red_line, -green_line, green_line]:
        # Find nearest bin index instead of exact match
        idx = np.searchsorted(plotted_bins, x_pos) - 1
        idx = np.clip(idx, 0, len(collapsed_n) - 1)
    
        idx_next = min(idx + 1, len(collapsed_n) - 1)
        height = max(collapsed_n[idx], collapsed_n[idx_next])
    
        plt.vlines(x_pos, 0, height, color='black', linestyle='--', linewidth=0.75, zorder=10)
    
    if is_drms_plot:
        plt.axvline(0, color=okabeito_blue, linestyle='-', linewidth=0.75, zorder = 3)

        plt.fill_between(np.array([-red_line, red_line]), ymin, y_max,
                     color='#e2f1ff', step='mid', label=f'±{red_line}%', zorder=1)
        plt.fill_between(np.array([-green_line, green_line]), ymin, ymax,
                     color='#c5e2ff', step='mid', label=f'±{green_line}%', zorder=2)
    else:
        plt.axvline(0, color=okabeito_green, linestyle='-', linewidth=0.75, zorder = 3)

        plt.fill_between(np.array([-red_line, red_line]), ymin, y_max,
                     color='#E6F5F1', step='mid', label=f'±{red_line}%', zorder=1)
        plt.fill_between(np.array([-green_line, green_line]), ymin, ymax,
                     color='#B3E2D5', step='mid', label=f'±{green_line}%', zorder=2)
    
    print(plot_title)
    n_pts = len(all_data)
    print('# data points: ', n_pts)
    if not is_cluster_plot and not is_stdev_plot:
        print('# spectra: ', n_undiscarded_spectra)
    if is_rdev_plot:
        less_than_5 = len([er for er in all_data if np.abs(er) < green_line])
        less_than_10 = len([er for er in all_data if np.abs(er) < red_line])
        less_than_15 = len([er for er in all_data if np.abs(er) < 15])
        print(f'<{green_line}%: {less_than_5/n_pts*100:.1f}%')
        print(f'<{red_line}%: {less_than_10/n_pts*100:.1f}%')
        print(f'<15%: {less_than_15/n_pts*100:.1f}%')
    elif is_abs_er_plot:
        less_than_1 = len([er for er in all_data if np.abs(er) < green_line])
        less_than_red = len([er for er in all_data if np.abs(er) < red_line])
        less_than_3 = len([er for er in all_data if np.abs(er) < 3])
        print(f'<{green_line}%: {less_than_1/n_pts*100:.1f}%')
        print(f'<{red_line}%: {less_than_red/n_pts*100:.1f}%')
        print(f'<3%: {less_than_3/n_pts*100:.1f}%')
    elif is_drms_plot:
        less_than_1 = len([er for er in all_data if np.abs(er) < 1])
        less_than_2 = len([er for er in all_data if np.abs(er) < 2])
        less_than_3 = len([er for er in all_data if np.abs(er) < 3])
        print(f'<{1}%: {less_than_1/n_pts*100:.1f}%')
        print(f'<{2}%: {less_than_2/n_pts*100:.1f}%')
        print(f'<3%: {less_than_3/n_pts*100:.1f}%')
    if add_fit:
        fit_histogram(n, bins)
    
    # Add legend
    plt.legend(
        loc='upper right',
        fontsize=fontsize,
        frameon=True,                # Enables the frame around the legend
        framealpha=1.0,              # Makes the frame fully opaque
        fancybox=False,              # Ensures a rectangular box
        edgecolor='black',           # Border color
        bbox_transform=ax.transAxes,
        facecolor='white'            # Background color (use bbox_to_anchor if needed)
    )
        
    # plt.show()
    plt.savefig(os.path.join(save_dir, file_name), dpi = 300)
    plt.close()

def extract_errors(df, ref_comp_dict, elements, pd_comp_str, min_fr = 0, max_fr = 100, min_abs_er = 0, min_rdev = 0, is_cluster = False):
    rdevs = []
    adevs = []
    abs_stdevs = []
    rel_stdevs = []
    adevs_to_stack = []
    rdevs_to_stack = []
    abs_stdevs_to_stack = []
    rel_stdevs_to_stack = []
    stdev_str = cnst.STDEV_DF_KEY
    stdevs_extracted = True
    is_first_eval = True
    n_spectra = 0
    for el in elements:
        if els_not_to_consider is None or el not in els_not_to_consider:
            exp_value = ref_comp_dict[el]
            # Cut data with fractions lower than min_fr
            if exp_value > min_fr and exp_value < max_fr:
                at_frs = df[el + pd_comp_str].to_numpy()
                if is_cluster: # Clusters are saved in %
                    stdev = float(clusters[el + stdev_str + pd_comp_str].iloc[0])
                    abs_stdevs.append(stdev)
                    
                    rel_stdev = stdev / at_frs *100
                    rel_stdevs.append(float(rel_stdev[0]))
                    
                    if el in els_to_stack_hist:
                        abs_stdevs_to_stack.append(float(stdev))
                        rel_stdevs_to_stack.append(float(rel_stdev[0]))
                            
                                
                el_adevs = (at_frs - exp_value)

                # Filter by minimal absolute error
                el_adevs = np.array([er for er in el_adevs if np.abs(er)>=min_abs_er])
                
                # Filter by minimal relative error
                el_adevs = np.array([er for er in el_adevs if np.abs(er/exp_value)>=min_rdev/100])
                
                el_rdevs =  el_adevs / exp_value * 100
                if el in els_to_stack_hist:
                    rdevs_to_stack += list(el_rdevs)
                    adevs_to_stack += list(el_adevs)

                rdevs += list(el_rdevs)
                adevs += list(el_adevs)
                
                if is_first_eval and not is_cluster:
                    is_first_eval = False
                    n_spectra += len(el_rdevs)
                
                if is_cluster and any(np.abs(el_er) >= 9.95 for el_er in el_rdevs):
                    print(f"{el}: {el_rdevs[0]:.1f}% RDEV from {at_frs[0]:.1f}% instead of {exp_value:.1f}%")
                    print(f"ADEV: {el_adevs[0]:.1f}%")
                    
                if is_cluster and any(np.abs(el_er) >= 2.95 for el_er in el_adevs):
                    print(f"{el}: {el_adevs[0]:.1f}% ADEV from {at_frs[0]:.1f}% instead of {exp_value:.1f}%")
                    print(f"RDEV: {el_rdevs[0]:.1f}%")
                
    return rdevs, adevs, abs_stdevs, rel_stdevs, stdevs_extracted, adevs_to_stack, rdevs_to_stack, n_spectra, abs_stdevs_to_stack, rel_stdevs_to_stack


def get_cluster_info(analysis_dir, clusters_to_ignore):
    clusters = None
    
    # Load obtained clusters
    try:
        clusters = pd.read_csv(os.path.join(analysis_dir, 'Clusters.csv'), index_col=0)
    except: #old files
        clusters = pd.read_csv(os.path.join(analysis_dir, 'Phases.csv'), index_col=0)
    
    # Filter out clusters to ignore
    try:
        for cl in clusters_to_ignore:
            clusters = clusters.drop(cl)
    except:
        pass
    
    try:
        n_pts_cluster = clusters['n_points'][0]
    except:
        ValueError(f"No Cluster file present in {analysis_dir}")
        
    if len(clusters) > 1:
        raise ValueError("Missing value for clusters to ignore")
    
    return clusters, n_pts_cluster


def get_analysis_dir(sample_dir, folders):
    analysis_dir = os.path.join(sample_dir, sorted(folders)[-1])
    return analysis_dir

#%% Code
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


# Plot config
fontsize = 15
labelpad = 12
global_lw = 0.5

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = fontsize         # General font size (default 10)
plt.rcParams['axes.titlesize'] = fontsize    # Title font size
plt.rcParams['axes.labelsize'] = fontsize    # Axis label font size
plt.rcParams['xtick.labelsize'] = fontsize   # X-axis tick label font size
plt.rcParams['ytick.labelsize'] = fontsize   # Y-axis tick label font size
plt.rcParams.update({
    'lines.linewidth': global_lw,           # plot lines
    'axes.linewidth': global_lw,            # axes border
    'xtick.major.width': global_lw,         # major ticks on x-axis
    'ytick.major.width': global_lw,         # major ticks on y-axis
    'xtick.minor.width': global_lw,         # minor ticks on x-axis
    'ytick.minor.width': global_lw,         # minor ticks on y-axis
    'grid.linewidth': global_lw,            # grid lines
    'patch.linewidth': global_lw            # patch objects (e.g., legend frame)
})

image_format = ".png"

if min_rdev > 0 or min_abs_er > 0:
    # Code not adapted to plot stdevs with any of these options
    plot_stdevs = False
else:
    plot_stdevs = True

frequency_els_list = []

n_unfiltered_spectra = 0
n_undiscarded_spectra = [0 for i in range(len(max_analytical_errors))]

all_rdevs_to_stack = [[] for i in range(len(max_analytical_errors))]
all_adevs_to_stack = [[] for i in range(len(max_analytical_errors))]
cluster_rdevs_to_stack = []
cluster_adevs_to_stack = []
cluster_rstdevs_to_stack = []
cluster_astdevs_to_stack = []
all_rdevs = [[] for i in range(len(max_analytical_errors))]
all_adevs = [[] for i in range(len(max_analytical_errors))]
all_cluster_rdevs = []
all_cluster_adevs = []
all_abs_stdevs = []
all_rel_stdevs = []
all_cluster_drms_at = []
all_cluster_drms_w = []
all_adevsstdevs_diffs = []
all_adevsstdevs_rel_diffs = []
folders_not_found = []
n_adevs_larger_stdev = 0
samples_with_insufficient_points = []
samples_with_no_stddev = []
for sample in samples:
    # Get sample info
    sample_ID = sample['sample_ID']
    formula = sample['formula']
    
    clusters_to_ignore = []
    if sample.get('clusters_to_ignore') is not None:
        clusters_to_ignore = sample.get('clusters_to_ignore')
        
    ref_comp = Composition(formula)
    filtered_ref_comp = Composition({el: amt for el, amt in ref_comp.items() if el.symbol not in undetectable_els})
    elements = [el.symbol for el in filtered_ref_comp.elements]
    
    # Get dictionary of reference composition
    if comp_type == 'at_fr':
        pd_comp_str = cnst.AT_FR_DF_KEY
        pd_phases_str = cnst.AT_FR_DF_KEY
        comp = Composition(filtered_ref_comp.fractional_composition).get_el_amt_dict()
        ref_comp_dict = {el : float(at_fr) *100 for el, at_fr in comp.items()}
    elif comp_type == 'w_fr':
        pd_comp_str = cnst.W_FR_DF_KEY
        pd_phases_str = cnst.W_FR_DF_KEY
        comp = filtered_ref_comp.to_weight_dict
        ref_comp_dict = {el : float(w_fr) *100 for el, w_fr in comp.items()}

    
    print_single_separator()
    print(f"Extracting data from sample {sample_ID}, with formula {formula}") 
    
    sample_dir = get_sample_dir(results_dir, sample_ID)
    all_folders = [f for f in os.listdir(sample_dir) if os.path.isdir(os.path.join(sample_dir, f)) and 'Analysis' in f]
    
    n_pts_cluster = 0 # initialise
    # Check if any folder has the piece of string to identify which one should be used for this analysis
    if any(f_string_hist in folder for folder in all_folders):
        folders = [f for f in all_folders if f_string_hist in f and 'OLD' not in f]
        analysis_dir = get_analysis_dir(sample_dir, folders)
        clusters, n_pts_cluster = get_cluster_info(analysis_dir, clusters_to_ignore)
        
        
    if n_pts_cluster == 0:
        folders_not_found.append(sample_ID)
        print(f"No Analysis folder with string {f_string_hist} were found")
        continue
    
    print(f"Using compositions from folder {analysis_dir}")
    n_samples += 1
    
    # Load measured compositions
    comps = pd.read_csv(os.path.join(analysis_dir, 'Compositions.csv'), index_col = 0)
    
    n_unfiltered_spectra += len(comps)
    
    # Filter out compositions to ignore
    for cl_n in clusters_to_ignore:
        comps = comps[comps["Cluster_ID"] != cl_n]
    
    # Histogram plots of all measurements
    for i, an_er in enumerate(max_analytical_errors):
        mask_ok_quant_flags = []
        if all_discarded_spectra:
            mask_ok_comps = comps[cnst.QUANT_FLAG_DF_KEY].isin([4,5,6,7,8])  | (
            comps[cnst.QUANT_FLAG_DF_KEY].isin([0,-1]) &
            (np.abs(comps[cnst.AN_ER_DF_KEY]) >= an_er)
            )
        elif discarded_spectra_no_aner_considered:
            mask_ok_comps = comps[cnst.QUANT_FLAG_DF_KEY].isin([4,5,6,7,8])
        else:
            mask_ok_comps = (
                    comps[cnst.QUANT_FLAG_DF_KEY].isin(acceptable_quant_flags) &
                    (np.abs(comps[cnst.AN_ER_DF_KEY]) < an_er)
                    )   
        
        
        comps_filtered = comps[mask_ok_comps]
        
        rdevs, adevs, _, _, _, adevs_to_stack, rdevs_to_stack, n_spectra, _, _ \
            = extract_errors(comps_filtered, ref_comp_dict, elements, pd_comp_str, 
                             min_fr, max_fr, min_abs_er, min_rdev)
    
        all_rdevs[i] += rdevs
        all_adevs[i] += adevs
        
        n_undiscarded_spectra[i] += n_spectra
        
        all_rdevs_to_stack[i] += rdevs_to_stack
        all_adevs_to_stack[i] += adevs_to_stack
    
    # Histogram plots of cluster averages
    if n_pts_cluster >= n_pts_thresh:
        frequency_els_list += list(ref_comp.get_el_amt_dict().keys())   # Update list of elements
        
        # extract rdev and adev values
        cl_rdevs, cl_adevs, abs_stdevs, rel_stdevs, stdevs_extracted, adevs_to_stack, rdevs_to_stack, _, \
            abs_stdevs_to_stack, rel_stdevs_to_stack \
            = extract_errors(clusters, ref_comp_dict, elements, pd_phases_str, min_fr, 
                             max_fr, min_abs_er, min_rdev, is_cluster = True)
        all_cluster_rdevs += cl_rdevs
        all_cluster_adevs += cl_adevs
        
        cluster_rdevs_to_stack += rdevs_to_stack
        cluster_adevs_to_stack += adevs_to_stack 
        cluster_rstdevs_to_stack += rel_stdevs_to_stack 
        cluster_astdevs_to_stack += abs_stdevs_to_stack
        
        all_cluster_drms_w.append(clusters[cnst.RMS_DIST_DF_KEY+ cnst.W_FR_DF_KEY])
        drms_at = clusters[cnst.RMS_DIST_DF_KEY+ cnst.AT_FR_DF_KEY][0]
        all_cluster_drms_at.append(drms_at)
        if drms_at > 3:
            print(f"!!!High d_RMS of {drms_at:.1f}")

        # extract measurement uncertainties
        if plot_stdevs:
            
            if stdevs_extracted:
                all_abs_stdevs += abs_stdevs
                
                # Count how many fractions have ADEV value larger than the measured stdev
                all_adevsstdevs_diffs += list(np.abs(cl_adevs) - np.array(abs_stdevs))
                all_adevsstdevs_rel_diffs += list(np.abs(cl_adevs)/np.array(abs_stdevs))
                n_adevs_larger_stdev_sample = np.sum(np.abs(cl_adevs) > 1* np.array(abs_stdevs))
                n_adevs_larger_stdev += n_adevs_larger_stdev_sample
                
                all_rel_stdevs += rel_stdevs
            else:
                samples_with_no_stddev.append(sample_ID)
    else:
        samples_with_insufficient_points.append(sample_ID)


print_double_separator()
if f'_minrdev{min_abs_er}' not in added_text:
    for i, an_er in enumerate(max_analytical_errors):
        make_histogram(all_rdevs[i], is_cluster_plot = False, is_rdev_plot = True, is_abs_er_plot = False,
                       rdevs_to_stack = all_rdevs_to_stack[i],
                       adevs_to_stack = all_adevs_to_stack[i],
                       n_undiscarded_spectra = n_undiscarded_spectra[i],
                       plot_title = f'RDEV, max an er {an_er} w%', file_name = f"bin{hist_bin}_rdev_aner{int(an_er)}"+added_text+image_format, add_fit = add_fit)
        make_histogram(all_adevs[i], is_cluster_plot = False, is_rdev_plot = False, is_abs_er_plot = True,
                       rdevs_to_stack = all_rdevs_to_stack[i],
                       adevs_to_stack = all_adevs_to_stack[i],
                       n_undiscarded_spectra = n_undiscarded_spectra[i],
                       plot_title = f'ADEV, max an er {an_er} w%', file_name = f"bin{hist_bin}_adev_aner{int(an_er)}"+added_text+image_format, add_fit =add_fit)

make_histogram(all_cluster_rdevs, is_cluster_plot = True, is_rdev_plot = True, is_abs_er_plot = False,
               rdevs_to_stack = cluster_rdevs_to_stack,
               plot_title = 'RDEV, clusters', file_name = f"bin{hist_bin}_rdev_clusters"+added_text+image_format,add_fit = add_fit)

make_histogram(all_cluster_adevs, is_cluster_plot = True, is_rdev_plot = False, is_abs_er_plot = True,
               adevs_to_stack = cluster_adevs_to_stack,
               plot_title = 'ADEV, clusters', file_name = f"bin{hist_bin}_adev_clusters"+added_text+image_format, add_fit =add_fit)

if plot_stdevs:
    make_histogram(all_rel_stdevs, is_cluster_plot = True, is_rdev_plot = True, is_abs_er_plot = False, is_stdev_plot=True, 
                   cluster_rstdevs_to_stack = cluster_rstdevs_to_stack, plot_title = 'Relative stdevs, clusters', file_name = f"bin{hist_bin}_rstdev_clusters"+added_text+image_format, add_fit = add_fit)
    make_histogram(all_abs_stdevs, is_cluster_plot = True, is_rdev_plot = False, is_abs_er_plot = True, is_stdev_plot=True,
                   cluster_astdevs_to_stack = cluster_astdevs_to_stack, plot_title = 'Absolute stdevs, clusters', file_name = f"bin{hist_bin/4}_astdev_clusters"+added_text+image_format, add_fit =add_fit)
    make_histogram(all_cluster_drms_at, is_cluster_plot = True, is_rdev_plot = False, is_abs_er_plot = False, is_stdev_plot=False, is_drms_plot = True, 
               cluster_rstdevs_to_stack = [], plot_title = 'd_RMS at%, clusters', file_name = f"bin{hist_bin}_drmsat_clusters"+added_text+image_format, add_fit = add_fit)
    make_histogram(all_cluster_drms_w, is_cluster_plot = True, is_rdev_plot = False, is_abs_er_plot = False, is_stdev_plot=False, is_drms_plot = True,
                   cluster_rstdevs_to_stack = [], plot_title = 'd_RMS w%, clusters', file_name = f"bin{hist_bin}_drmsw_clusters"+added_text+image_format, add_fit = add_fit)
    make_histogram(all_adevsstdevs_diffs, is_cluster_plot = True, is_rdev_plot = False, is_abs_er_plot = True, is_stdev_plot=False, is_drms_plot = False,
                   cluster_rstdevs_to_stack = [], plot_title = 'ADEV - stdev, clusters', file_name = f"bin{hist_bin}_diffADEVstdev_clusters"+added_text+image_format, add_fit = add_fit)
    make_histogram(all_adevsstdevs_rel_diffs, is_cluster_plot = True, is_rdev_plot = False, is_abs_er_plot = True, is_stdev_plot=False, is_drms_plot = True,
                   cluster_rstdevs_to_stack = [], plot_title = 'ADEV / stdev, clusters', file_name = f"bin{hist_bin}_relADEVstdev_clusters"+added_text+image_format, add_fit = add_fit)

#%% Plot periodic table with counts for each element analysed
fontsize_periodictable = 14.5
# Count occurrences of each element
element_counts = Counter(frequency_els_list)

# Convert counts to a format suitable for pymatviz
element_data = {el: count for el, count in element_counts.items()}
added_text
if f'_minrdev{min_abs_er}' not in added_text:
    # Generate the periodic table heatmap
    fig = ptable_heatmap(data = element_data, return_type = 'figure', log = True)
    # Adjust figure size to make squares bigger
    fig.set_size_inches(10, 6)  # Increase width and height
    
    # Get the colorbar from the figure and modify its ticks
    cbar = fig.axes[-1]  # Colorbar is usually the last axis in Matplotlib figures
    max_tick = int(10 * np.floor(max(element_data.values())/ 10))
    cbar.set_xticks([1, 10, max_tick],[1, 10, max_tick])  # Ensure 60 is displayed
    # cbar.tick_params(axis='both', colors=color_SEMEDS_paper_gray) 
    cbar.tick_params(axis='both', colors=color_SEMEDS_paper_gray, labelsize=fontsize_periodictable, which='both')
    for spine in cbar.spines.values():
        spine.set_color(color_SEMEDS_paper_gray)
        # spine.set_linewidth(0.8)


    # Identify and modify the colorbar title
    for text in cbar.get_children():
        if isinstance(text, plt.Text) and text.get_text() == "Element Count":
            text.set_fontweight('normal')  # Remove bold formatting
            text.set_color(color_SEMEDS_paper_gray)  # Change text color
            text.set_fontsize(fontsize_periodictable)  # Change font size (adjust the value as needed)
            break  # Stop once the correct text is found


    # Save the plot to a file
    file_path = os.path.join(save_dir, f"Els_counts{added_text}" + image_format)
    fig.savefig(file_path, dpi=300, bbox_inches="tight")
    plt.close()

    
# Print samples whose correct folder was not found
print_double_separator()
if len(folders_not_found) > 0:
    print(f"The following samples did not have an Analysis folder containing '{f_string_hist}'")
    for sample in folders_not_found:
        print(sample)
else:
    print("Analysis folders were correctly found for all samples")

print_double_separator()
print(f"{n_samples} samples were analysed")
print(f"{n_unfiltered_spectra} spectra were analysed")

print_single_separator()
print(f"{n_adevs_larger_stdev} elemental fractions had ADEV larger than the measurement stdev")

print_double_separator()
if plot_stdevs:
    if len(samples_with_no_stddev) > 0:
        print(f"Could not extract atomic standard deviations from the following {len(samples_with_no_stddev)} samples:")
        for sample in samples_with_no_stddev:
            print(sample)
    else:
        print("Could extract atomic standard deviations from all analysed samples")
    print_single_separator()

# Print samples that did not have enough points making up the cluster
if len(samples_with_insufficient_points) > 0:
    print(f"The following samples did not have at least {n_pts_thresh} measurement points")
    for sample in samples_with_insufficient_points:
        print(sample)
else:
    print(f"All analysed samples had at least {n_pts_thresh} measurement points")


    
    
    
    
    
    
    
    
    
    
    