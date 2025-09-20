#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:31:09 2025

@author: Andrea
"""

import pandas as pd
import numpy as np
import os
from pymatgen.core import Element
from scipy import constants


# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

### Import libraries
# Scattering factors from Henke database: https://henke.lbl.gov/optical_constants/asf.html
# Accessed: January 15th, 2025

# Directory of scattering factor files
xray_sf_dir = os.path.join(script_dir, 'Xray_scattering_factors')

def _xray_wavelength_from_energy(energy_eV):
    # Convert energy from eV to Joules
    energy_J = energy_eV * constants.e  # 1 eV = 1.60218e-19 J
    
    # Calculate wavelength using the formula λ = hc / E
    wavelength_m = (constants.h * constants.c) / energy_J
    
    return wavelength_m


# Returns dictionary of xray lines with their energy (in keV) and weights for element el
def xray_mass_absorption_coeff(element, energies):
    '''
    Returns mass absorption coefficient of element at given energies.
    Calculated from atomic form factors from Henke database (https://henke.lbl.gov/optical_constants/asf.html)

        Parameters
    ----------
    el : str or int
        The element symbol (e.g., 'O' for oxygen) or atomic number (e.g., 8 for oxygen).

    energy : float, list, or numpy.array
        The X-ray energy values in keV at which to compute the mass absorption coefficient.
        The energy can be provided as a single value (float), a list of energies, 
        or a numpy array. All input energies will be converted into a numpy array 
        for consistent processing.

    Returns
    -------
    mass_abs : numpy.array
        The mass absorption coefficient (in cm²/g) at the provided energies.
        The result is computed using the interpolation of atomic scattering factors 
        (`f2`) for the given energies and based on the atomic form factor data.
    '''
    
    # Ensure energy is a NumPy array, even if it's provided as a list or single value
    if isinstance(energies, list):
        energies = np.array(energies)
    elif isinstance(energies, (int, float)):
        energies = np.array([energies])  # Convert single value to array
    elif isinstance(energies, (np.ndarray)):    
        energies = energies.copy()
    else:
        raise ValueError("energies must be float, int or np.array")
    
    # Convert to eV
    energies *= 1000
    
    # Ensure element is Element object from pymatgen
    if isinstance(element, str):
        el_obj = Element(element)
    elif isinstance(element, int):
        el_obj = Element.from_Z(element)
    elif isinstance(element, Element):
        el_obj = element
    else:
        raise ValueError(f"You entered el={element}. The value of 'el' must be the element abbreviation or its atomic number")
    
    # Load scattering factors for element
    el_sf_dir = os.path.join(xray_sf_dir, el_obj.symbol.lower() + '.nff')
    data = np.loadtxt(el_sf_dir, delimiter='\t', skiprows=1, usecols=[0,1,2])
    df = pd.DataFrame(data, columns = ('E(eV)', 'f1', 'f2'))
    
    # Replace missing data with np.nan
    df.replace(-9999., np.nan, inplace=True)
    
    # Interpolate f2 values for the provided energy (single value, list, or np.array)
    energy_values = df['E(eV)'].values  # The energy values from the DataFrame
    f2_values = df['f2'].values        # The corresponding f2 values from the DataFrame
    f2_interp = np.interp(energies, energy_values, f2_values)
    
    # Calculate mass absorption coefficient from f2
    # https://gisaxs.com/index.php/Absorption_length#:~:text=The%20absorption%20length%20arises%20from,displaystyle%20%5Cmu%20/%5Crho%20%7D%20.
    atomic_mass = el_obj.atomic_mass  # in amu
    e_r = constants.physical_constants['classical electron radius'][0]
    Na = constants.Avogadro
    xray_lambda = _xray_wavelength_from_energy(energies)
    mass_abs = Na * 2 * f2_interp * e_r * xray_lambda / atomic_mass * 10**4 # convert from m^2 to cm^2
    
    # If the input was a single value, return a single value
    if energies.size == 1:
        mass_abs = mass_abs[0]  # Return the single value
    
    return mass_abs
    
# print(xray_mass_absorption_coeff('Ca', [100,200]))
