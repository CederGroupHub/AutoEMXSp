#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process Particle Size Distribution File and Recalculate Particle Statistics

Choose the statistics result to start from (each segmentation of the frames, e.g. by another segmentation model,
gives its own result: run once with show_results_only = True to list them), then the particles to filter out from
{Sample_ID_Par_Sizes.csv file} either by frame and particle number in the frame ('Frame ID' and 'Particle #'
headers) or filtering out all particles in a specific frame ('Frame ID' header). The particles are read from the
frames of the sample ledger (ledger.json), as segmented by the chosen result: the recalculated statistics are added
as a new result, made the active one (recording the particles and frames left out; earlier results are kept, and
rerunning without filters gives back the statistics of the measurement). Samples measured with older versions of
AutoEMX, without particle sizes in their ledger, are read from their particle sizes file, by particle ID.
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
# Choose the statistics result
# --------------------------
show_results_only = False  # True: only list the results of the sample (id, segmentation, particles...)
result_id = None           # Result to recalculate (its measurement and segmentation). None: the active one

# --------------------------
# Define particles and frames to filter out (all of them: those of the chosen result are not kept)
# --------------------------
particles_to_filter = [("A0", 1), ("A2", 3)]  # Particles to remove: (Frame ID, Particle #)
particles_IDs_to_filter = []        # Particle IDs to remove (older samples read from their particle sizes file)
frame_IDs_to_filter = []            # Example frame IDs to remove (strings match the Frame ID column)

# --------------------------
# Run
# --------------------------
from autoemx.runners.process_particle_statistics import list_size_results, process_particle_statistics

list_size_results(sample_ID, results_path)
if not show_results_only:
    process_particle_statistics(
        sample_ID=sample_ID,
        results_path=results_path,
        particles_to_filter=particles_to_filter,
        particle_ids_to_filter=particles_IDs_to_filter,
        frame_ids_to_filter=frame_IDs_to_filter,
        result_id=result_id,
    )
