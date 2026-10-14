#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process Particle Size Distribution File and Recalculate Particle Statistics

Choose particles to filter out from {Sample_ID_Par_Sizes.csv file} either by
particle ID ('Particle #' header) or filtering out all particles in a specific frame ('Frame ID' header).
The particles are read from the sample ledger (ledger.json), where the recalculated statistics replace those of the
measurement (recording the particles and frames left out; rerun without filters to restore them), or from the
particle sizes file for samples measured with older versions of AutoEMX.
The statistics are also saved with the suffix '_processed'.

Created on Wed Jul 31 10:40:26 2024

@author: Andrea
"""

# --------------------------
# Define sample
# --------------------------
sample_ID = 'example_particle_stats'  # Name of folder where data has been saved
results_path = 'input'  # Base directory where data is saved. It can find sample_ID also within its subdirectories

# --------------------------
# Define particles and frames to filter out
# --------------------------
particles_IDs_to_filter = [1, 3]        # Particle IDs to remove
frame_IDs_to_filter = []            # Example frame IDs to remove (strings match the Frame ID column)

# --------------------------
# Run
# --------------------------
from autoemx.runners.process_particle_statistics import process_particle_statistics

process_particle_statistics(
    sample_ID=sample_ID,
    results_path=results_path,
    particle_ids_to_filter=particles_IDs_to_filter,
    frame_ids_to_filter=frame_IDs_to_filter,
)
