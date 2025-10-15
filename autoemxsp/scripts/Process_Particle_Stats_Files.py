#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process Particle Size Distribution File and Recalculate Particle Statistics

Choose particles to filter out from {Sample_ID_Par_Sizes.csv file} either by
particle ID ('Particle #' header) or filtering out all particles in a specific frame ('Frame ID' header)

Created on Wed Jul 31 10:40:26 2024

@author: Andrea
"""

import os
import pandas as pd

from autoemxsp.core.EM_particle_finder import EM_Particle_Finder
import autoemxsp.tools.constants as cnst
from autoemxsp.tools.utils import get_sample_dir, print_single_separator

# --------------------------
# Define sample
# --------------------------
sample_ID = 'example_particle_stats'  # Name of folder where data has been saved
input_dir = 'input'  # Base directory where data is saved. It can find sample_ID also within its subdirectories

# --------------------------
# Define particles and frames to filter out
# --------------------------
particles_IDs_to_filter = [1, 3]        # Particle IDs to remove
frame_IDs_to_filter = []            # Example frame IDs to remove (strings match the Frame ID column)

#%%
# --------------------------
# Code
# --------------------------
# Check if sample exists
sample_dir = get_sample_dir(input_dir, sample_ID)
if not os.path.exists(sample_dir):
    raise FileNotFoundError(f"Could not find sample at {sample_dir}. Please check 'sample_ID' and 'input_dir'.")

# Load particle size data
par_data_path = os.path.join(sample_dir, f"{sample_ID}_{cnst.PARTICLE_SIZES_FILENAME}.csv")
if not os.path.exists(par_data_path):
    raise FileNotFoundError(f"Particle size file not found at {par_data_path}")

par_data = pd.read_csv(par_data_path)

# Store original IDs for reference before filtering
original_particle_ids = set(par_data[cnst.PAR_ID_DF_KEY])
original_frame_ids = set(par_data[cnst.FRAME_ID_DF_KEY])

# Determine which IDs actually exist in the dataset
missing_particles = [pid for pid in particles_IDs_to_filter if pid not in original_particle_ids]
missing_frames = [fid for fid in frame_IDs_to_filter if fid not in original_frame_ids]

# Perform the filtering
if particles_IDs_to_filter is not None:
    par_data = par_data[~par_data[cnst.PAR_ID_DF_KEY].isin(particles_IDs_to_filter)]

if frame_IDs_to_filter is not None:
    par_data = par_data[~par_data[cnst.FRAME_ID_DF_KEY].isin(frame_IDs_to_filter)]

# --- Report results ---
if not missing_particles and not missing_frames:
    print("✅ All requested particle IDs and frame IDs were present in the dataset and have been removed.")
else:
    print("⚠️ The following IDs were not found in the dataset (and therefore were not removed):")
    if missing_particles:
        print(f" - Missing particle IDs: {missing_particles}")
    if missing_frames:
        print(f" - Missing frame IDs: {missing_frames}")

# Re-calculate statistics and save filtered results
calculator = EM_Particle_Finder(None, None, results_dir=sample_dir, verbose=True)
calculator._sample_ID = sample_ID
calculator.analyzed_pars = list(
    zip(par_data[cnst.PAR_ID_DF_KEY], par_data[cnst.FRAME_ID_DF_KEY], par_data[cnst.PAR_AREA_UM_KEY])
)

if len(par_data) > 0:
    calculator.save_particle_statistics(output_file_suffix='_processed')
else:
    print("No particles left to build statistics. Ensure you do not filter out all particles through: "
          "particles_IDs_to_filter and frame_IDs_to_filter")