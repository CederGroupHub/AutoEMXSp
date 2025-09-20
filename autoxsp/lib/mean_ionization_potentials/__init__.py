#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 14:14:20 2025

@author: Andrea
"""
import os
import pandas as pd

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Directory of file
mean_ionization_potentials_dir = os.path.join(script_dir, 'MeanIonizationPotentials.csv')

# Load into pandas dataframe
J_df = pd.read_csv(mean_ionization_potentials_dir, header=0)

J_df.index = range(1, len(J_df) + 1)  # Set indices to atomic number (Z), starting from 1 for H

