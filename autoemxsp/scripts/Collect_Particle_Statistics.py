#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 10:40:26 2024

@author: Andrea
"""

import os
import numpy as np
import tkinter as tk
from tkinter import messagebox, filedialog
import json

from autoemxsp.runners import collect_particle_statistics
from autoemxsp.tools.config_classes import PowderMeasurementConfig
AVAILABLE_SEGMENTATION_MODELS = PowderMeasurementConfig.AVAILABLE_PAR_SEGMENTATION_MODELS

#%% Functions

# File to store default values
defaults_file = 'defaults_particle_stats.json'

def load_defaults():
    if os.path.exists(defaults_file):
        with open(defaults_file, 'r') as f:
            return json.load(f)
    return {
        'sample_ID': '',
        'sample_center_x': '',
        'sample_center_y': '',
        'num_particles_to_analyse': '200',
        'par_segmentation_model' : PowderMeasurementConfig.DEFAULT_PAR_SEGMENTATION_MODEL,
        'min_particle_diameter_um': '0.3',
        'max_particle_diameter_um': '20',
        'working_distance_mm' : '5',
        'carbon_tape_diameter_mm': '5',
        'is_manual_navigation' : 'False',
        'auto_detect_Ctape' : 'no',
        'auto_contrast_brightness': 'no',
        'brightness': '0.4107',
        'contrast': '0.4221',
        'results_dir': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Results')
    }

def save_defaults(user_input):
    with open(defaults_file, 'w') as f:
        json.dump(user_input, f)
    

def get_user_input():
    defaults = load_defaults()
    
    def submit():
        try:
            sample_ID = sample_ID_entry.get()
            if not sample_ID:
                raise ValueError("Sample ID is required.")

            sample_center_x = sample_center_x_entry.get()
            sample_center_y = sample_center_y_entry.get()
            if not sample_center_x or not sample_center_y:
                raise ValueError("Sample Center X and Y are required.")

            num_particles_to_analyse = num_particles_entry.get()
            if not num_particles_to_analyse:
                raise ValueError("Number of Particles to Analyse is required.")

            min_diameter = min_diameter_entry.get()
            max_diameter = max_diameter_entry.get()
            if not min_diameter or not max_diameter:
                raise ValueError("Min and Max Particle Diameter are required.")

            carbon_tape_diameter = carbon_tape_diameter_entry.get()
            if not carbon_tape_diameter:
                raise ValueError("Carbon Tape Diameter is required.")

            working_distance_val = working_distance_entry.get()
            try:
                working_distance_val = float(working_distance_val)
            except ValueError:
                raise ValueError("Working Distance must be a valid number.")

            auto_contrast_brightness = auto_contrast_brightness_var.get()
            user_input['auto_contrast_brightness'] = auto_contrast_brightness
            
            if auto_contrast_brightness == 'no':
                brightness = brightness_entry.get()
                contrast = contrast_entry.get()
                if not brightness:
                    raise ValueError("Brightness is required.")
                if not contrast:
                    raise ValueError("Contrast is required.")
                user_input['brightness'] = float(brightness)
                user_input['contrast'] = float(contrast)
            else:
                user_input['brightness'] = defaults['brightness']
                user_input['contrast'] = defaults['contrast']

            # Collect values into user_input
            user_input['sample_ID'] = sample_ID
            user_input['sample_center_x'] = float(sample_center_x)
            user_input['sample_center_y'] = float(sample_center_y)
            user_input['num_particles_to_analyse'] = int(num_particles_to_analyse)
            user_input['min_particle_diameter_um'] = float(min_diameter)
            user_input['max_particle_diameter_um'] = float(max_diameter)
            user_input['carbon_tape_diameter_mm'] = float(carbon_tape_diameter)
            user_input['working_distance_mm'] = working_distance_val
            
            user_input['par_segmentation_model'] = par_segmentation_var.get()
            user_input['is_manual_navigation'] = bool(is_manual_navigation_var.get())
            user_input['auto_detect_Ctape'] = bool(auto_detect_Ctape_var.get())

            user_input['results_dir'] = defaults['results_dir']

            save_defaults(user_input)
            root.destroy()

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))

    user_input = {}
    root = tk.Tk()
    root.title("User Options")

    def show_info(message):
        messagebox.showinfo(message=message)

    # Sample ID
    tk.Label(root, text="Sample ID").grid(row=0, column=0)
    sample_ID_entry = tk.Entry(root, width=20)
    sample_ID_entry.grid(row=0, column=1, columnspan=2)
    sample_ID_entry.insert(0, defaults['sample_ID'])

    # Sample Center X/Y
    tk.Label(root, text="Sample Center X, Y").grid(row=1, column=0)
    sample_center_x_entry = tk.Entry(root, width=10)
    sample_center_x_entry.grid(row=1, column=1)
    sample_center_y_entry = tk.Entry(root, width=10)
    sample_center_y_entry.grid(row=1, column=2)
    sample_center_x_entry.insert(0, defaults['sample_center_x'])
    sample_center_y_entry.insert(0, defaults['sample_center_y'])

    # Info button
    info_message = 'Ensure one or more particles are present at these coordinates for auto-focus.'
    info_button = tk.Button(root, text="?", command=lambda: show_info(info_message), padx=2)
    info_button.grid(row=1, column=3, sticky='w')

    # Number of Particles
    tk.Label(root, text="Number of Particles to Analyse").grid(row=2, column=0)
    num_particles_entry = tk.Entry(root, width=10)
    num_particles_entry.insert(0, defaults['num_particles_to_analyse'])
    num_particles_entry.grid(row=2, column=1, columnspan=2)

    # Min/Max Particle Diameter
    tk.Label(root, text="Min & Max Particle Diameter (µm)").grid(row=3, column=0)
    min_diameter_entry = tk.Entry(root, width=10)
    min_diameter_entry.grid(row=3, column=1)
    min_diameter_entry.insert(0, defaults['min_particle_diameter_um'])
    max_diameter_entry = tk.Entry(root, width=10)
    max_diameter_entry.grid(row=3, column=2)
    max_diameter_entry.insert(0, defaults['max_particle_diameter_um'])

    # Carbon Tape Diameter
    tk.Label(root, text="Carbon Tape Diameter (mm)").grid(row=4, column=0)
    carbon_tape_diameter_entry = tk.Entry(root, width=10)
    carbon_tape_diameter_entry.insert(0, defaults['carbon_tape_diameter_mm'])
    carbon_tape_diameter_entry.grid(row=4, column=1, columnspan=3)

    # Working Distance (mm)
    tk.Label(root, text="Working Distance (mm)").grid(row=5, column=0)
    working_distance_entry = tk.Entry(root, width=10)
    working_distance_entry.insert(0, defaults['working_distance_mm'])
    working_distance_entry.grid(row=5, column=1, columnspan=3)

    # Segmentation Model Dropdown
    tk.Label(root, text="Segmentation Model").grid(row=6, column=0)
    par_segmentation_var = tk.StringVar(value=defaults['par_segmentation_model'])
    par_model_menu = tk.OptionMenu(root, par_segmentation_var, *AVAILABLE_SEGMENTATION_MODELS)
    par_model_menu.grid(row=6, column=1, columnspan=3)
    if len(AVAILABLE_SEGMENTATION_MODELS) == 1:
        par_model_menu.config(state="disabled")

    # Manual Navigation Checkbox
    is_manual_navigation_var = tk.IntVar(value=int(defaults['is_manual_navigation'] == 'True'))
    tk.Checkbutton(root, text="Manual Navigation", variable=is_manual_navigation_var).grid(row=7, column=0, columnspan=2)

    # Auto-detect Carbon Tape Checkbox
    auto_detect_Ctape_var = tk.IntVar(value=int(defaults['auto_detect_Ctape'] == 'yes'))
    tk.Checkbutton(root, text="Auto Detect Carbon Tape", variable=auto_detect_Ctape_var).grid(row=7, column=2, columnspan=2)

    # Auto Contrast / Brightness
    tk.Label(root, text="Auto Adjust Brightness and Contrast").grid(row=8, column=0)
    auto_contrast_brightness_var = tk.StringVar(value=defaults['auto_contrast_brightness'])
    tk.Radiobutton(root, text="Yes", variable=auto_contrast_brightness_var, value='yes').grid(row=8, column=1)
    tk.Radiobutton(root, text="No", variable=auto_contrast_brightness_var, value='no').grid(row=8, column=2)

    # Contrast
    tk.Label(root, text="Contrast").grid(row=9, column=0)
    contrast_entry = tk.Entry(root, width=10)
    contrast_entry.insert(0, defaults['contrast'])
    contrast_entry.grid(row=9, column=1, columnspan=3)

    # Brightness
    tk.Label(root, text="Brightness").grid(row=10, column=0)
    brightness_entry = tk.Entry(root, width=10)
    brightness_entry.insert(0, defaults['brightness'])
    brightness_entry.grid(row=10, column=1, columnspan=3)

    def update_fields():
        if auto_contrast_brightness_var.get() == 'yes':
            contrast_entry.grid_forget()
            brightness_entry.grid_forget()
        else:
            contrast_entry.grid(row=9, column=1, columnspan=3)
            brightness_entry.grid(row=10, column=1, columnspan=3)

    update_fields()
    auto_contrast_brightness_var.trace("w", lambda *args: update_fields())

    # Submit Button
    submit_button = tk.Button(root, text="Submit", command=submit)
    submit_button.grid(row=11, columnspan=4)

    root.mainloop()
    return user_input

def select_save_directory(default_directory):
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    directory = filedialog.askdirectory(initialdir=default_directory, title="Select Directory to Save Data")
    root.destroy()
    return directory


#%% Unpack tkinter
user_inputs = get_user_input()

# Unpack values from the dictionary
sample_ID = user_inputs['sample_ID']
sample_center_position = (user_inputs['sample_center_x'], user_inputs['sample_center_y'])
samples = [{'ID': sample_ID, 'pos' : sample_center_position}]

num_particles_to_analyse = user_inputs['num_particles_to_analyse']
min_particle_diameter_um = user_inputs['min_particle_diameter_um']
max_particle_diameter_um = user_inputs['max_particle_diameter_um']
carbon_tape_diameter_mm = user_inputs['carbon_tape_diameter_mm']
auto_contrast_brightness = user_inputs['auto_contrast_brightness']
if auto_contrast_brightness == 'no':
    auto_contrast_brightness = False
else:
    auto_contrast_brightness = True
brightness = user_inputs['brightness']
contrast = user_inputs['contrast']
par_segmentation_model = user_inputs['par_segmentation_model']
working_distance = user_inputs['working_distance_mm']
is_manual_navigation = user_inputs['is_manual_navigation']
is_auto_substrate_detection = user_inputs['auto_detect_Ctape']

# Directory for saving data
defaults = load_defaults()
default_directory = defaults['results_dir']
results_dir = select_save_directory(default_directory)
if not results_dir:
    print("No directory selected. Exiting.")
    exit()

# Update defaults file with the new directory if a valid one was selected
user_inputs['results_dir'] = results_dir
save_defaults(user_inputs)

# Define options based on user
max_area_par = (max_particle_diameter_um/2)**2 * np.pi
min_area_par = (min_particle_diameter_um/2)**2 * np.pi
C_tape_r = (carbon_tape_diameter_mm - 0.2 *carbon_tape_diameter_mm)/ 2

#%% General  Configuration
# =============================================================================
# Microscope Configuration
# =============================================================================
microscope_ID = 'PhenomXL'
microscope_type = 'SEM'
sample_substrate_type = 'Ctape'
sample_substrate_shape = 'circle'
sample_substrate_width_mm = 12 # Al stub

# =============================================================================
# Powder sample options
# =============================================================================
powder_meas_cfg_kwargs = dict(
    is_manual_particle_selection = False,
    is_known_powder_mixture_meas = False,
    par_search_frame_width_um = None,
    max_n_par_per_frame=np.inf,
    max_area_par=max_area_par,
    min_area_par=min_area_par,
    par_segmentation_model = par_segmentation_model,
    par_brightness_thresh=100,
)


collect_particle_statistics(
    samples = samples,
    n_par_target = num_particles_to_analyse,
    microscope_ID = microscope_ID,
    microscope_type = microscope_type,
    sample_halfwidth = C_tape_r,
    sample_substrate_type = sample_substrate_type,
    sample_substrate_shape = sample_substrate_shape,
    sample_substrate_width_mm = sample_substrate_width_mm,
    working_distance = working_distance,
    is_manual_navigation = is_manual_navigation,
    is_auto_substrate_detection = is_auto_substrate_detection,
    auto_adjust_brightness_contrast = auto_contrast_brightness,
    contrast = contrast,
    brightness = brightness,
    powder_meas_cfg_kwargs = powder_meas_cfg_kwargs,
    development_mode = False,
    verbose = True,
    results_dir = results_dir
)

