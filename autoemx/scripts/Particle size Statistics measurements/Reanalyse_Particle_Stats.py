#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reanalyse a particle size statistics measurement from its saved frame images

Segments again the frame images saved in 'SEM images' with other settings (e.g. brightness threshold), and rebuilds
the statistics of the sample (in ledger.json) with the same size bands and frames. The frames searched by EDS on the
sample are reanalysed too. The statistics are also saved with the suffix '_reanalysed'.

Created on Fri Oct 16 2026

@author: Andrea
"""

# --------------------------
# Define sample
# --------------------------
sample_ID = 'example_particle_stats'  # Name of folder where data has been saved
results_path = 'input'  # Base directory where data is saved. It can find sample_ID also within its subdirectories

# --------------------------
# Segmentation settings (None: those of the measurement)
# --------------------------
par_brightness_thresh = None   # 8-bit brightness above which pixels belong to particles
par_segmentation_model = None  # e.g. 'threshold_bright'
min_par_pixels = None          # Minimum pixels of a particle for its area to be measured (only raised: bands stay)
include_eds_frames = True      # Also reanalyse the frames searched by EDS on the sample

# --------------------------
# Run
# --------------------------
from autoemx.runners.reanalyze_particle_statistics import reanalyze_particle_statistics

reanalyze_particle_statistics(
    sample_ID=sample_ID,
    results_path=results_path,
    par_brightness_thresh=par_brightness_thresh,
    par_segmentation_model=par_segmentation_model,
    min_par_pixels=min_par_pixels,
    include_eds_frames=include_eds_frames,
)
