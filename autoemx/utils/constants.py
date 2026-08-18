#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 11:44:19 2025

This file contains the strings used across the various AutoEMX modules

@author: Andrea
"""

### DICTIONARY KEYS

# Fitting results
PEAK_AREA_KEY = 'area'
PEAK_SIGMA_KEY = 'sigma'
PEAK_CENTER_KEY = 'center'
PEAK_FWHM_KEY = 'fwhm'
PEAK_INTENSITY_KEY = 'peak_intensity'
BACKGROUND_INT_KEY = 'background_intensity'
PEAK_TH_ENERGY_KEY = 'th_energy'
PEAK_HEIGHT_KEY = 'height'
PB_RATIO_KEY = 'PB_ratio'

# Quantification results
R_SQ_KEY = 'r_squared'
REDCHI_SQ_KEY = 'redchi_sq'
AN_ER_KEY = 'Analytical error'
COMP_AT_FR_KEY = 'Comp_at_fr'
COMP_W_FR_KEY = 'Comp_w_fr'
Z_MEAN_W_KEY = 'mass-averaged'
Z_MEAN_AT_KEY = 'atomic-averaged'
Z_MEAN_STATHAM_KEY = 'Statham2016'
Z_MEAN_MARKOWICZ_KEY = 'Markowicz1984'

# Spectrum collection info dictionary keys
DATETIME_KEY = 'datetime'
MICROSCOPE_CFG_KEY = 'microscope_cfg'
SAMPLE_CFG_KEY = 'sample_cfg'
MEASUREMENT_CFG_KEY = 'measurement_cfg'
SAMPLESUBSTRATE_CFG_KEY = 'sample_substrate_cfg'
QUANTIFICATION_CFG_KEY = 'quant_cfg'
ACQUISITION_CFG_KEY = 'acquisition_cfg'
POWDER_MEASUREMENT_CFG_KEY = 'powder_meas_cfg'
BULK_MEASUREMENT_CFG_KEY = 'bulk_meas_cfg'
EXP_STD_MEASUREMENT_CFG_KEY = 'exp_stds_cfg'
CLUSTERING_CFG_KEY = 'clustering_cfg'
PLOT_CFG_KEY = 'plot_cfg'

# Image annotation dictionary keys
ANNOTATION_CIRCLE_KEY = 'circle'
ANNOTATION_TEXT_KEY = 'text'

# Clustering info dictionary keys
N_SP_ACQUIRED_KEY = 'n_spectra_collected'
N_SP_USED_KEY = 'n_spectra_used'
N_CLUST_KEY = 'n_clusters'
WCSS_KEY = 'wcss'
SIL_SCORE_KEY = 'sil_score'
REF_NAME_KEY = 'refs'
CONF_SCORE_KEY = 'conf_score'
MOLAR_FR_MEAN_KEY = 'mean'
MOLAR_FR_STDEV_KEY = 'stddev'

# Experimental standard measurements
MEAN_PB_KEY = 'PB_mean'
MEAS_PB_DF_KEY = 'Measured_PB'
STDEV_PB_DF_KEY = 'Stdev_PB'
COR_PB_DF_KEY =  'Corrected_PB'
REL_ER_PERCENT_PB_DF_KEY = 'Rel_stdev_PB (%)'
STD_ID_KEY = 'ID'
STD_FORMULA_KEY = 'Formula'
STD_TYPE_KEY = 'Std_type'
STD_Z_KEY = 'Mean_Z'
STD_MEAN_ID_KEY = 'Mean'
STD_USE_FOR_MEAN_KEY = 'Use_for_mean_calc'



### DATAFRAME HEADERS

# Headers of Data.csv and Compositions.csv files, containing the collected spectra and their compositions
SP_ID_DF_KEY = 'Spectrum #'
PAR_ID_DF_KEY = 'Particle #'
FRAME_ID_DF_KEY = 'Frame ID'
SP_X_COORD_DF_KEY = 'x'
SP_Y_COORD_DF_KEY = 'y'
SP_X_PIXEL_COORD_DF_KEY = 'x_px'
SP_Y_PIXEL_COORD_DF_KEY = 'y_px'
SPECTRUM_DF_KEY = 'Spectrum'
BACKGROUND_DF_KEY = 'Background'
REAL_TIME_DF_KEY = 'Real_time'
LIVE_TIME_DF_KEY = 'Live_time'
QUANT_FLAG_DF_KEY = 'Quant_flag'
COMMENTS_DF_KEY = 'Comments'
AN_ER_DF_KEY = 'An er w%'#'Analytical error' #'An er w%'
CND_DF_KEY = 'Cnd'
CS_RAW_CND_DF_KEY = 'CS_raw'
CS_CND_DF_KEY = 'CS_cnd'
MIX_DF_KEY = 'Mix'
CS_MIX_DF_KEY = 'CS_mix'
MIX_MOLAR_RATIO_DF_KEY = 'Mol_Ratio'
MIX_FIRST_COMP_MEAN_DF_KEY = 'X1_mean'
MIX_FIRST_COMP_STDEV_DF_KEY = 'X1_stdev'


# Headers of Clusters.csv files
CL_ID_DF_KEY = 'Cluster_ID'
N_PTS_DF_KEY = 'n_points'
RMS_DIST_DF_KEY = 'RMS_dist'
WCSS_DF_KEY = 'wcss'

# Headers common to both types of files
AT_FR_DF_KEY = '_at%' # Added to elements
W_FR_DF_KEY = '_w%' # Added to elements
STDEV_DF_KEY = '_std' # Added to elements

# Particle statistics keys
PAR_AREA_UM_KEY = 'Area (μm²)'
PAR_EQ_D_KEY = 'Equivalent Diameter (μm)'



### FILES AND DIRECTORIES

# Directories
RESULTS_DIR = 'Results'
STDS_DIR = 'Std_measurements'
ANALYSIS_DIR = 'Analysis'
ANALYSIS_SUBDIR = 'analysis'  # new lowercase subdirectory for timestamped analysis folders
SPECTRA_DIR = 'spectra'  # directory containing per-spectrum pointer files
IMAGES_DIR = 'SEM images'
CALIBS_DIR = 'calibrations'
XRAY_SPECTRA_CALIBS_DIR = 'calibrations'
DETECTOR_CHANNEL_PARAMS_CALIBR_DIR = 'Detector_channel_params_calibs'
SDD_CALIBS_MEAS_DIR = 'SDD calibrations'
PAR_SEGMENTATION_MODELS_DIR = 'particle_segmentation_models'


# File names
STD_FILENAME = 'Stds'
ACQUISITION_INFO_FILENAME = 'Comp_analysis_configs'  # legacy config filename (kept during transition)
CONFIG_FILENAME = 'config'  # new config filename
LEDGER_FILENAME = 'ledger'  # new per-sample spectral ledger filename
LEDGER_FILEEXT = '.json'
ANALYSIS_TIMESTAMP_FORMAT = '%Y-%m-%d_%H-%M'  # timestamp format for analysis subdirectory names
DATA_FILENAME = 'Data'
MSA_SP_FILENAME = "EM_metadata.msa"
SPECTRUM_FILENAME_PREFIX = "spectrum_"
SPECTRUM_MAN_BACKGROUND_SUFFIX = "_man_bckgrnd"
# EMSA/MAS spectrum files accepted for ingest (same text format; .emsa is a common alias).
EMSA_SPECTRUM_EXTENSIONS = {".msa", ".msg", ".emsa"}
SPECTRUM_POINTER_EXTENSIONS = {".msa", ".msg", ".emsa", ".json"}
SPECTRUM_POINTER_EXT_PRIORITY = {".msa": 0, ".emsa": 0, ".msg": 1, ".json": 2}
VECTOR_FILEEXT = ".npy"
PARTICLE_STATS_FILENAME = 'Par_size_stats'
PARTICLE_SIZES_FILENAME = 'Par_sizes'
PARTICLE_STAT_HIST_FILENAME = 'Par_size_distribution_hist'
STDS_MEAS_FILENAME = 'Std_measurements'
STDS_RESULT_FILENAME = 'Std_results'
DATA_FILEEXT = '.csv'
COMPOSITIONS_FILENAME = 'Compositions'
CLUSTERS_FILENAME = 'Clusters'
CLUSTERING_PLOT_FILENAME = 'Clustering_plot'
ANALYSIS_CONFIG_SUMMARY_FILENAME = 'Analysis_config_summary'
CUSTOM_CLUSTERING_PLOT_FILENAME = 'Clustering_plot_custom'
CUSTOM_PLOT_FILENAME = 'custom_plot.py'
POWDER_MIXTURE_PLOT_FILENAME = 'Mixture_decomposition_plot'
CLUSTERING_PLOT_FILEEXT = '.png'
NAVCAM_IM_FILENAME = 'Analysed_region'
INITIAL_SEM_IM_FILENAME = 'Initial_Position'
DETECTOR_CONV_MATRICES_FILENAME = 'detector_response_convolution_matrices.json'
DETECTOR_EFFICIENCY_FILENAME ='detector_efficiency.msa'
DETECTOR_CHANNEL_PARAMS_CALIBR_FILENAME = 'detector_channel_params_calibs'


# json files keys
LIST_SPECTRAL_DATA_QUANT_KEYS = [
        SPECTRUM_DF_KEY,
        BACKGROUND_DF_KEY,
        REAL_TIME_DF_KEY,
        LIVE_TIME_DF_KEY,
    ]

LIST_SPECTRAL_DATA_KEYS = LIST_SPECTRAL_DATA_QUANT_KEYS + [
        COMMENTS_DF_KEY,
        QUANT_FLAG_DF_KEY,
    ]

LIST_SPECTRUM_COORDINATES_KEYS = [
    SP_ID_DF_KEY,
    FRAME_ID_DF_KEY,
    PAR_ID_DF_KEY,
    SP_X_COORD_DF_KEY,
    SP_Y_COORD_DF_KEY,
    SP_X_PIXEL_COORD_DF_KEY,
    SP_Y_PIXEL_COORD_DF_KEY,
    ] # Update this list if adding keys
# This allows correct loading when quantifying or analysing spectra after acquisition



### CONFIGURATION DATACLASSES
S_POWDER_SAMPLE_TYPE = "powder"
S_POWDER_CONTINUOUS_SAMPLE_TYPE = "powder_continuous"
S_BULK_SAMPLE_TYPE = "bulk"
S_BULK_ROUGH_SAMPLE_TYPE = "bulk_rough"
S_FILM_SAMPLE_TYPE = "film"

SQUARE_SUBSTRATE_SHAPE = 'square'
CIRCLE_SUBSTRATE_SHAPE = 'circle'
CTAPE_SUBSTRATE_TYPE = 'Ctape'
NONE_SUBSTRATE_TYPE = 'None'


### OTHER

# Clustering features
AT_FR_CL_FEAT = 'at_fr'
W_FR_CL_FEAT = 'w_fr'

# X-ray detector channel calibration keys
BEAM_CURRENT_KEY = 'spot_size'
SCALE_KEY = 'scale'
OFFSET_KEY = 'offset'
