#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 18:10:58 2024

@author: Andrea
"""
import json
import numpy as np

current_offset = -0.032415
current_scale = 0.009972

print(f"Current scale: {current_scale:.6f}")
print(f"Current offset: {current_offset:.6f}")

# ### Using Mn_Ka1 and Al_Ka1
# # We employ Mn_Ka1 and Cu_La1 peaks so that they are at sufficient distance:
# x_measured_en = np.mean([5.86996, 5.87007, 5.870139]) # Mn_Ka1
# y_measured_en = np.mean([1.479459, 1.4794615, 1.47952]) # Al_Ka1
# # y_measured_en = np.mean([1.4793158212574276, 1.4794232914235788, 1.4793835847103245]) # Al_Ka1

# x_th_en = 5.898 # Mn_Ka1
# y_th_en = 1.4869 # Al_Ka1

# ### Using Si_Ka1 and Fe_Ka1 from NIST standards
# # We employ Mn_Ka1 and Cu_La1 peaks so that they are at sufficient distance:
# x_measured_en = np.mean([6.37062, 6.36638, 6.37112]) # Fe_Ka1
# y_measured_en = np.mean([1.73316, 1.73215, 1.73136]) # Si_Ka1

# x_th_en = 6.403 # Fe_Ka1
# y_th_en = 1.740 # Si_Ka1

### Using Cu_Ka1 and Al_Ka1
# We employ Cu_Ka1 and Cu_La1 peaks so that they are at sufficient distance:
x_measured_en = np.mean([8.03903634, 8.04116283, 8.03926358, 8.03848308]) # Cu_Ka1
y_measured_en = np.mean([1.48622366, 1.48614360, 1.48607739]) # Al_Ka1

x_th_en = 8.046 # Cu_Ka1
y_th_en = 1.4869 # Al_Ka1


#%%
i_x = (x_measured_en - current_offset) / current_scale
i_y = (y_measured_en - current_offset) / current_scale

Dx = x_th_en - x_measured_en
Dy = y_th_en - y_measured_en

new_scale = (Dx - Dy + x_measured_en - y_measured_en) / (i_x - i_y)

new_offset = Dy + y_measured_en - i_y * new_scale

print(f"New scale: {new_scale:.6f}")
print(f"New offset: {new_offset:.6f}")

# # EDX modes employed to collect standards
# # For each mode, we assign (offset, scale) to calculate energy_vals
# EDX_modes = {'point' : {'spot_size' : 4.3, 'offset': -0.033, 'scale': 0.009957},
#                  'map' : {'spot_size' : 5.1, 'offset': -0.033, 'scale': 0.009921} }
#                 # 'image' : {'spot_size' : 3.3, 'offset': , 'scale': }
#                 # 'low' : {'spot_size' : 2.1, 'offset': , 'scale': }


# # Save EDX modes in json file, accessible when measuring spectra for quantification or when making standards        
# with open('phenom_EDX_modes.json', 'w') as file:
#     json.dump(EDX_modes, file)