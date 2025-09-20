#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:54:00 2025

Useful in case you need to modify the detector calibrations retroactively

@author: Andrea
"""
import os,json

def list_all_directories(root_dir):
    directories = []
    for dirpath, dirnames, _ in os.walk(root_dir):
        for dirname in dirnames:
            directories.append(os.path.join(dirpath, dirname))
    return directories


#%% For standard files
path_to_measurements = '/Users/Andrea_1/Desktop/Work/Codes/EDX/Std measurements'
all_directories = list_all_directories(path_to_measurements)

# Modify calibration of all measurements done with name High Cnts
dirs = [d for d in all_directories if any(s in d for s in ['50k', '250k', '1000k'])]

for d in dirs:
    with open(os.path.join(d,"EDX_info_point_mode.json"), 'r') as f:
        meas_file = json.load(f)
        print(meas_file['energy_zero'], meas_file['bin_width'])
        meas_file['energy_zero'] = -0.034245
        meas_file['bin_width'] = 0.009972
        print(meas_file['energy_zero'], meas_file['bin_width'])
    with open(os.path.join(d,"EDX_info_point_mode.json"), 'w') as f: 
        json.dump(meas_file, f, ensure_ascii=False, indent=4)
        
        
#%% For result files
path_to_measurements = '/Users/Andrea_1/Desktop/Work/Codes/EDX/Results'
all_directories = [d for d in list_all_directories(path_to_measurements) if not any(s in d for s in ['SEM images', 'Analysis'])]

# Modify calibration of all measurements done with name High Cnts
dirs = [d for d in all_directories if any(s in d for s in ['mineral'])]

spectra_info_file_name = 'Spectra_collection_info.json'

for d in dirs:
    with open(os.path.join(d,spectra_info_file_name), 'r') as f:
        meas_file = json.load(f)
        print(meas_file['energy_zero'], meas_file['bin_width'])
        meas_file['energy_zero'] = -0.034245
        meas_file['bin_width'] = 0.009972
        print(meas_file['energy_zero'], meas_file['bin_width'])
    with open(os.path.join(d, spectra_info_file_name), 'w') as f: 
        json.dump(meas_file, f, ensure_ascii=False, indent=4)