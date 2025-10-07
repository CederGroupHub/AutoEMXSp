#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility Functions for EDS Spectrum Analysis and Data Handling

Created on Fri Jun 28 11:50:53 2024

@author: Andrea Giunto

This module provides a collection of utility functions and classes to support 
X-ray energy-dispersive spectroscopy (EDS) analysis workflows, especially in 
combination with the main spectrum modeling and fitting modules.

Main Features
-------------
- **Compositional Conversions:**  
  Functions to convert between atomic and weight fractions, parse chemical formulas, 
  and generate LaTeX representations for chemical formulas, using `pymatgen`.

- **String and Table Utilities:**  
  Functions for formatted printing of tables and separators for clear console output.
  Includes an `AlphabetMapper` class for converting integer indices to Excel-style 
  column labels (e.g., A, B, ..., Z, AA, AB, ...).
  
- **Image Utilities:**  
  Function to draw a scale bar on an image

- **File and Directory Handling:**  
  Functions to locate sample directories within a results folder and to load 
  .msa/.msg spectral data files (especially from Thermo Fisher Phenom systems), 
  including metadata parsing and energy scale calculation.

- **Miscellaneous Tools:**  
  - A simple Tkinter-based GUI prompt for user interaction during workflows.
  - Custom exception classes for clearer error handling in EDS/SEM applications.

Typical Usage
-------------
This module is intended to be imported by higher-level scripts or modules that 
perform EDS spectrum modeling and data analysis. Example usages include:

- Converting between atomic and weight fractions for quantification:
    >>> atomic_to_weight_fr([0.5, 0.5], ['K', 'Cl'])
    >>> weight_to_atomic_fr([0.5, 0.5], ['K', 'Cl'])

- Parsing and standardizing chemical formulas:
    >>> get_std_comp_from_formula('LiCoO2')

- Formatting output in the console:
    >>> print_nice_1d_row('W_fr', [0.123, 0.456, 0.789])
    >>> print_single_separator()

- Loading EDS spectra from MSA files:
    >>> energy, counts, meta = load_msa('my_spectrum.msa')

- Prompting the user for action during a workflow:
    >>> Prompt_User("Pause", "Please reposition the sample and press OK.").run()

- Handling custom errors in EDS/SEM processing:
    >>> raise EDSError("Detector not found", code=404)

Dependencies
------------
- numpy
- pymatgen
- tkinter (for GUI prompts)
- Standard Python libraries: os, sys, re

Notes
-----
This module is designed for flexibility and ease of use in scientific data analysis 
and should be adapted as needed for specific workflows or instrument environments.
"""
#TODO update docs with new functions

# Standard library imports
import os
import sys
import re
import json
import warnings
import tkinter as tk
import pandas as pd
import ast

# Third-party imports
import cv2
import numpy as np
from pymatgen.core import Element, Composition

# Typing imports
from typing import Any, Sequence, List, Optional, Tuple, Dict, Union

import autoemxsp.tools.constants as cnst


#%% Compositions and formulas
def atomic_to_weight_fr(
    atomic_fractions: Union[Sequence[float], np.ndarray],
    elements: List[str],
    verbose: bool = True
) -> np.ndarray:
    """
    Convert atomic fractions to weight fractions for a given set of elements.

    Parameters
    ----------
    atomic_fractions : Sequence[float] or np.ndarray
        The atomic fractions of the elements. Should sum to approximately 1.
    elements : List[str]
        List of element symbols (e.g., ["Si", "O", "Al"]).
    verbose : bool, optional
        If True, prints warnings about normalization. Default is True.

    Returns
    -------
    weight_fractions : np.ndarray
        The corresponding weight fractions, summing to 1.

    Raises
    ------
    ValueError
        If the lengths of atomic_fractions and elements do not match.

    Examples
    --------
    >>> atomic_to_weight_fr([0.333, 0.667], ["Si", "O"])
    array([0.467..., 0.532...])
    """
    atomic_fractions = np.asarray(atomic_fractions, dtype=np.float64)
    if len(atomic_fractions) != len(elements):
        raise ValueError("Length of atomic_fractions and elements must match.")

    total_atomic = np.sum(atomic_fractions)
    if not np.isclose(total_atomic, 1.0, atol=1e-3):
        if verbose:
            print(f"Warning: atomic_fractions sum to {total_atomic:.4f}. Normalizing to 1.")
        atomic_fractions = atomic_fractions / total_atomic

    molar_masses = np.array([Element(el).atomic_mass for el in elements], dtype=np.float64)
    if np.any(molar_masses == 0):
        raise ValueError("One or more elements have zero molar mass. Check element symbols.")

    # Calculate raw weights
    raw_weights = atomic_fractions * molar_masses
    total_weight = np.sum(raw_weights)
    weight_fractions = raw_weights / total_weight

    return weight_fractions

# Test
if __name__ == "__main__":
    print(atomic_to_weight_fr([0.5,0.5], ['K', 'Cl']))


def weight_to_atomic_fr(
    weight_fractions: Union[Sequence[float], np.ndarray],
    elements: List[str],
    verbose: bool = True
) -> np.ndarray:
    """
    Convert weight fractions to atomic fractions for a given set of elements.

    Parameters
    ----------
    weight_fractions : Sequence[float] or np.ndarray
        The weight fractions of the elements. Should sum to approximately 1.
    elements : List[str]
        List of element symbols (e.g., ["Si", "O", "Al"]).
    verbose : bool, optional
        If True, prints warnings about normalization. Default is True.

    Returns
    -------
    atomic_fractions : np.ndarray
        The corresponding atomic fractions, summing to 1.

    Raises
    ------
    ValueError
        If the lengths of weight_fractions and elements do not match.

    Examples
    --------
    >>> weight_to_atomic_fr([0.467, 0.533], ["Si", "O"])
    array([0.333..., 0.666...])
    """
    weight_fractions = np.asarray(weight_fractions, dtype=np.float64)
    if len(weight_fractions) != len(elements):
        raise ValueError("Length of weight_fractions and elements must match.")

    total_weight = np.sum(weight_fractions)
    if not np.isclose(total_weight, 1.0, atol=1e-3):
        if verbose:
            print(f"Warning: weight_fractions sum to {total_weight:.4f}. Normalizing to 1.")
        weight_fractions = weight_fractions / total_weight

    molar_masses = np.array([Element(el).atomic_mass for el in elements], dtype=np.float64)
    # Avoid division by zero in case of invalid element symbols
    if np.any(molar_masses == 0):
        raise ValueError("One or more elements have zero molar mass. Check element symbols.")

    # Calculate atomic amounts (number of moles per element)
    raw_atomics = weight_fractions / molar_masses
    total_atoms = np.sum(raw_atomics)
    atomic_fractions = raw_atomics / total_atoms

    return atomic_fractions


# Test
if __name__ == "__main__":    
    print(weight_to_atomic_fr([0.5,0.5], ['K', 'Cl']))


def to_latex_formula(formula: str, include_dollar_signs: bool = True) -> str:
    """
    Convert a chemical formula string into its LaTeX representation.

    Supports nested parentheses, element subscripts, and group multipliers.
    For example:
        'Al2(SO4)3' -> 'Al$_{2}$(SO$_{4}$)$_{3}$'

    Args:
        formula (str): The chemical formula as a string.
        include_dollar_signs (bool): Whether to wrap LaTeX in $...$.

    Returns:
        str: The LaTeX-formatted formula.
    """
    def convert(s: str) -> str:
        result = ''
        i = 0
        while i < len(s):
            if s[i] == '(':
                # Find the matching closing parenthesis for the group
                depth = 1
                j = i + 1
                while j < len(s) and depth > 0:
                    if s[j] == '(':
                        depth += 1
                    elif s[j] == ')':
                        depth -= 1
                    j += 1
                group = s[i+1:j-1]  # Content inside the parentheses

                # Check for a multiplier after the group (e.g., (SO4)3)
                multiplier = ''
                while j < len(s) and (s[j].isdigit() or s[j] == '.'):
                    multiplier += s[j]
                    j += 1

                group_latex = convert(group)
                if multiplier:
                    result += f"({group_latex})$_{{{multiplier}}}$"
                else:
                    result += f"({group_latex})"
                i = j
            else:
                # Match element symbol and optional subscript
                m = re.match(r'([A-Z][a-z]*)([0-9.]+)?', s[i:])
                if m:
                    element, subscript = m.groups()
                    if subscript and subscript != '1':
                        result += f"{element}$_{{{subscript}}}$"
                    else:
                        result += element
                    i += len(m.group(0))
                else:
                    # Skip invalid characters (should not occur in well-formed formulas)
                    i += 1
        return result
    
    latex_formula = convert(formula)
    if include_dollar_signs:
        return latex_formula
    else:
        # Remove all $ signs (in case someone wants a pure LaTeX string)
        return latex_formula.replace('$', '')

# Test
if __name__ == "__main__":
    print(to_latex_formula('Bi2(Fe4O9)0.33'))     

    

#%% Strings
def print_nice_1d_row(
    first_col: Any,
    row: Sequence[Any],
    floatfmt: str = ".2f",
    first_col_width: int = 10,
    col_width: int = 10
) -> None:
    """
    Print a single row for a table: first entry (label) left-aligned, rest as formatted floats.
    Used to print quantification results during ZAF corrections

    Args:
        first_col: The value for the first (label) column.
        row: Sequence of values for the remaining columns (numbers or strings).
        floatfmt: Format for floating point numbers (default: '.3f').
        first_col_width: Width for the first (label) column.
        col_width: Width for each numeric column.
    """
    first = str(first_col).ljust(first_col_width)
    rest = []
    for val in row:
        try:
            rest.append(f"{float(val):{col_width}{floatfmt}}")
        except (ValueError, TypeError):
            rest.append(str(val).rjust(col_width))
    print(first + "".join(rest))

# Test
if __name__ == "__main__":
    print_nice_1d_row('', ['Fe', 'Mn', 'O'])  # header row with element names
    print_nice_1d_row('W_fr', [1.23456, 23.4, 0.005678])
    print_nice_1d_row('Z_vals', [0.123, 2.3456, 0.0005678])


def print_element_fractions_table(formula):
    """
    Given a chemical formula, print a table with elements, atomic %, and weight % using pymatgen.
    """
    comp = Composition(formula)
    elements = [el.symbol for el in comp.elements]
    at_fr = [comp.get_atomic_fraction(el) for el in comp.elements]
    w_fr = [comp.get_wt_fraction(el) for el in comp.elements]
    
    print_nice_1d_row(formula, elements)
    print_nice_1d_row('at%', [x*100 for x in at_fr])
    print_nice_1d_row('w%',  [x*100 for x in w_fr])
    
# Test
if __name__ == "__main__":
    print_element_fractions_table('Mn2SiO4')  # header row with element names
    
    

def print_single_separator():
    """Print a single-line separator (50 dashes) for visual clarity."""
    print('-' * 50)
    sys.stdout.flush()
# print_single_separator()




def print_double_separator():
    """Print a double-line separator (50 equals signs) for visual clarity."""
    print('=' * 50)
    sys.stdout.flush()
# print_double_separator()


class AlphabetMapper:
    """
    Maps a zero-based integer index to Excel-style column letters.
    0 -> 'A', 1 -> 'B', ..., 25 -> 'Z', 26 -> 'AA', etc.
    
    Used for labeling frames analyzed during particle search.
    """
    def __init__(self):
        self.alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self.n_letters = len(self.alphabet)

    def get_letter(self, index: int) -> str:
        """
        Convert a zero-based index to Excel-style column letters.

        Parameters
        ----------
        index : int
            Zero-based index (e.g., 0 for 'A', 25 for 'Z', 26 for 'AA', ...)

        Returns
        -------
        str
            Excel-style column letters.

        Raises
        ------
        ValueError
            If index is negative.
        """
        if index < 0:
            raise ValueError("Index must be non-negative.")

        letters = ''
        while True:
            index, remainder = divmod(index, self.n_letters)
            letters = self.alphabet[remainder] + letters
            if index == 0:
                break
            index -= 1  # Excel-style: after 'Z' comes 'AA', not 'BA'
        return letters

# Test the AlphabetMapper
if __name__ == "__main__":
    mapper = AlphabetMapper()
    test_indices = [0, 1, 25, 26, 27, 51, 52, 701, 702, 703, 1378]
    for i in test_indices:
        print(f"{i}: {mapper.get_letter(i)}")

#%% Files
def load_configurations_from_json(json_path, config_classes_dict):
    """
    Load configuration dataclasses and metadata from a spectrum collection info JSON file.

    Parameters
    ----------
    json_path : str
        Path to the JSON file saved by EMXSp_Composition_Analyzer._save_spectrum_collection_info.
    config_classes_dict : dict
        Mapping from JSON keys to dataclass types, e.g.:
            {'sample_cfg': SampleConfig, ...}
        See configurations in tools/config_classes.py:
            - MicroscopeConfig: Settings for microscope hardware, calibration, and imaging parameters.
            - SampleConfig: Defines the sample’s identity, elements, and spatial properties.
            - SampleSubstrateConfig: Specifies the substrate composition and geometry supporting the sample.
            - MeasurementConfig: Controls measurement type, beam parameters, and acquisition settings.
            - FittingConfig: Parameters for spectral fitting and background handling.
            - QuantConfig: Options for quantification and filtering of X-ray spectra.
            - PowderMeasurementConfig: Settings for analyzing powder samples and particle selection.
            - BulkMeasurementConfig: Settings for analyzing bulk samples.
            - ExpStandardsConfig: Settings for experimental standard measurements
            - ClusteringConfig: Configures clustering algorithms and feature selection for data grouping.
            - PlotConfig: Options for saving, displaying, and customizing plots.

    Returns
    -------
    configs : dict
        Dictionary of configuration objects reconstructed from JSON, keyed by their JSON key.
        If a key from config_classes_dict is missing in the JSON, it will not be present in configs.
    metadata : dict
        Dictionary of any additional metadata (e.g., timestamp) found in the JSON.

    Raises
    ------
    FileNotFoundError
        If the JSON file does not exist.

    Example
    -------
    >>> config_classes = {'sample_cfg': SampleConfig, ...}
    >>> configs, metadata = load_configurations_from_json('acquisition_info.json', config_classes)
    """
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    configs = {}
    metadata = {}
    for key, cls in config_classes_dict.items():
        if key in data:
            configs[key] = cls(**data[key])
        else:
            configs[key] = None
            # warnings.warn(f"Configuration key '{key}' not found in JSON file '{json_path}'.")

    # Any other keys in JSON are treated as metadata
    for key, value in data.items():
        if key not in config_classes_dict:
            metadata[key] = value

    return configs, metadata


def extract_spectral_data(data_csv_path):
    """
    Extract spectra quantification, spectral data, and coordinates from Data.csv file.

    Parameters
    ----------
    data_csv_path : str
        Path to the Data.csv file.
    cnst : module or object
        Should provide all the necessary attribute keys as in your code.

    Returns
    -------
    spectra_quant : list of dict or None
        List of quantification results per spectrum (None if not quantified).
        If all entries are None, returns None.
    spectral_data : dict
        Dictionary of lists for each spectral data column (e.g., spectrum, background, real_time, live_time, etc.).
        If a column is missing, the value is an empty list.
    sp_coords : list of dict
        List of dicts for each spectrum's coordinates, keys as in cnst.LIST_SPECTRUM_COORDINATES_KEYS.
    df : pandas.DataFrame
        The loaded DataFrame from the CSV file.
    """
    df = pd.read_csv(data_csv_path)

    # --- Extract spectral data ---
    spectral_keys = cnst.LIST_SPECTRAL_DATA_KEYS
    spectral_data = {}

    for key in spectral_keys:
        if key in df.columns:
            # For spectrum and background, convert string "[1.1,2.2,...]" to list of floats
            if key in [cnst.SPECTRUM_DF_KEY, cnst.BACKGROUND_DF_KEY]:
                spectral_data[key] = [
                    ast.literal_eval(val) if pd.notnull(val) else None
                    for val in df[key]
                ]
            else:
                spectral_data[key] = df[key].to_numpy()
        else:
            spectral_data[key] = [None] * len(df)

    # --- Extract spectral coordinates ---
    available_cols = [col for col in cnst.LIST_SPECTRUM_COORDINATES_KEYS if col in df.columns]
    sp_coords = df[available_cols].to_dict(orient='records')

    # --- Extract quantification results ---
    spectra_quant = []
    el_atfr_cols = [c for c in df.columns if c.endswith(cnst.AT_FR_DF_KEY)]
    el_wfr_cols = [c for c in df.columns if c.endswith(cnst.W_FR_DF_KEY)]
    elements = [c.replace(cnst.AT_FR_DF_KEY, '') for c in el_atfr_cols]
    for _, row in df.iterrows():
        # Only consider quantified if at least one atfr or wfr is not null
        if any(pd.notnull(row[c]) for c in el_atfr_cols + el_wfr_cols):
            comp_atfr = {el: row[el + cnst.AT_FR_DF_KEY]/100 for el in elements if pd.notnull(row[el + cnst.AT_FR_DF_KEY])}
            comp_wfr  = {el: row[el + cnst.W_FR_DF_KEY]/100 for el in elements if pd.notnull(row[el + cnst.W_FR_DF_KEY])}
            quant_dict = {
                cnst.COMP_AT_FR_KEY: comp_atfr,
                cnst.COMP_W_FR_KEY: comp_wfr,
                cnst.AN_ER_KEY: row[cnst.AN_ER_DF_KEY]/100 if cnst.AN_ER_DF_KEY in df.columns else [],
                cnst.R_SQ_KEY: row[cnst.R_SQ_KEY] if cnst.R_SQ_KEY in df.columns else [],
                cnst.REDCHI_SQ_KEY: row[cnst.REDCHI_SQ_KEY] if cnst.REDCHI_SQ_KEY in df.columns else [],
            }
            spectra_quant.append(quant_dict)
        else:
            spectra_quant.append(None)

    if all(x is None for x in spectra_quant):
        spectra_quant = None

    return spectra_quant, spectral_data, sp_coords, df


def make_unique_path(parent_dir: str, base_name: str, extension: str = None) -> str:
    """
    Generate a unique file or directory path inside parent_dir based on base_name.
    If a path with the base name exists, appends a counter (e.g., 'Sample1_2').
    Optionally, add an extension for files.

    Parameters
    ----------
    parent_dir : str
        The parent directory in which to generate the new path.
    base_name : str
        The base name for the new file or directory.
    extension : str, optional
        The file extension (e.g., 'txt' or '.txt'). If None, treat as directory.

    Returns
    -------
    unique_path : str
        The full, unique path (not created).

    Raises
    ------
    ValueError
        If `parent_dir` or `base_name` is invalid.

    Example
    -------
    >>> make_unique_path('./results', 'Sample1')
    './results/Sample1'
    >>> make_unique_path('./results', 'Sample1', extension='txt')
    './results/Sample1.txt'
    >>> make_unique_path('./results', 'Sample1', extension='txt')  # If exists
    './results/Sample1_2.txt'
    """
    if not isinstance(parent_dir, str) or not isinstance(base_name, str):
        raise ValueError("Both `parent_dir` and `base_name` must be strings.")

    # Normalize extension
    ext = ''
    if extension:
        ext = extension if extension.startswith('.') else f'.{extension}'

    counter = 1
    if ext:
        unique_path = os.path.join(parent_dir, f"{base_name}{ext}")
    else:
        unique_path = os.path.join(parent_dir, base_name)

    while os.path.exists(unique_path):
        counter += 1
        if ext:
            unique_path = os.path.join(parent_dir, f"{base_name}_{counter}{ext}")
        else:
            unique_path = os.path.join(parent_dir, f"{base_name}_{counter}")

    return unique_path


def get_sample_dir(results_path: str, sample_ID: str) -> str:
    """
    Find the directory for a given sample_ID inside results_path or any of its subdirectories.

    If the sample_ID folder exists both directly under results_path and in a subdirectory,
    raises a RuntimeError due to ambiguity.

    Parameters
    ----------
    results_path : str
        The root directory to search in.
    sample_ID : str
        The folder name to search for.

    Returns
    -------
    sample_dir : str
        The full path to the sample_ID directory.

    Raises
    ------
    RuntimeError
        If the sample_ID folder exists both in root and a subdirectory.
    FileNotFoundError
        If the sample_ID folder is not found.
    """
    sample_dir = os.path.join(results_path, sample_ID)
    found_in_root = os.path.isdir(sample_dir)
    found_in_subdir: Optional[str] = None

    for root, dirs, _ in os.walk(results_path):
        if root == results_path:
            continue  # Skip top-level
        if sample_ID in dirs:
            found_in_subdir = os.path.join(root, sample_ID)
            break  # Stop at first occurrence

    if found_in_root and found_in_subdir:
        print(f"WARNING: '{sample_ID}' folder exists both in '{results_path}' and in subdirectory '{found_in_subdir}'", file=sys.stderr)
        raise RuntimeError(f"Ambiguous '{sample_ID}' folder location.")

    if found_in_root:
        return sample_dir
    elif found_in_subdir:
        return found_in_subdir
    else:
        raise FileNotFoundError(f"'{sample_ID}' folder not found in '{os.path.abspath(results_path)}' or its subdirectories.")



def load_msa(filepath: str) -> Tuple[np.ndarray, np.ndarray, Dict[str, str]]:
    """
    Load a .msa or .msg file containing Y-only spectral data and compute the energy scale.
    Designed for raw spectra exported from Thermo Fisher Phenom systems.
    May work with other EMSA/MAS format files, but minor variations can occur.

    Parameters
    ----------
    filepath : str
        Path to the .msa or .msg file.

    Returns
    -------
    energy : np.ndarray
        Energy values computed from OFFSET and XPERCHAN.
    counts : np.ndarray
        Measured counts per energy channel.
    metadata : dict
        Parsed metadata from the header (all values as strings).
    """
    metadata = {}
    counts = []
    in_data_section = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('#SPECTRUM'):
                in_data_section = True
                continue

            if line.startswith('#'):
                # Parse header metadata
                parts = line[1:].split(':', maxsplit=1)
                if len(parts) == 2:
                    key, value = parts
                    metadata[key.strip()] = value.strip()
                continue

            if in_data_section:
                # Parse counts (allow comma or space separated)
                token = line.replace(',', '').strip()
                try:
                    counts.append(float(token))
                except ValueError:
                    continue  # Skip lines that aren't valid numbers

    # Retrieve required metadata, with sensible defaults
    try:
        npoints = int(float(metadata.get("NPOINTS", len(counts))))
        offset = float(metadata.get("OFFSET", 0.0))
        xperchan = float(metadata.get("XPERCHAN", 1.0))
    except (TypeError, ValueError) as e:
        raise ValueError(f"Error reading required metadata fields: {e}")

    # Adjust counts to match NPOINTS
    if len(counts) > npoints:
        counts = counts[:npoints]
    elif len(counts) < npoints:
        counts += [0.0] * (npoints - len(counts))

    energy = offset + np.arange(npoints) * xperchan
    counts = np.array(counts)

    return energy, counts, metadata

# Test
if __name__ == "__main__":
    energy, counts, meta = load_msa(os.path.join('lib','PhenomXL_detector efficiency.msa'))
    print(energy)
    print(counts)
    print(meta)
    

#%% Images
def draw_scalebar(image, pixel_size_um, bar_width = 0.25):
    """
    Draw a scale bar on the given image.

    The scale bar is drawn as a filled white rectangle in the bottom-left corner of the image,
    with a label indicating its length in micrometers (um). The actual scale bar length is chosen
    from a set of standard values (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200 um) to best match
    the 'bar_width', expressed as fraction of image width.

    Parameters
    ----------
    image : np.ndarray
        Input image (grayscale or color, as a NumPy array). The scale bar will be drawn directly on this image.
    pixel_size_um : float
        Size of a pixel in micrometers (um/pixel).
    bar_width : float, optional
        Target width of scale bar, in fraction of image width. Default: 0.25

    Returns
    -------
    image : np.ndarray
        The same image with the scale bar and label drawn on it.

    Notes
    -----
    - The function detects if the image is grayscale or color and draws the scale bar in white accordingly.
    - The scale bar is placed in the bottom-left corner with a margin from the edges.
    - The label is centered above the scale bar rectangle.
    """
    
    # Check if image is color or grey scale
    if len(image.shape) == 3:
        im_height, im_width, _ = image.shape
        white_color = (255, 255, 255)
    else:
        im_height, im_width = image.shape
        white_color = 255

    # Calculate the desired scale bar length in pixels
    desired_length_pixels = im_width * bar_width

    # Define the possible scale bar lengths in um
    possible_lengths_um = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]

    # Find the closest possible scale bar length to the desired length
    best_length_um = min(possible_lengths_um, key=lambda x: abs(x / pixel_size_um - desired_length_pixels))

    # Set the scale bar length and label
    bar_length_pixels = int(best_length_um / pixel_size_um)
    bar_label = f"{best_length_um} um"

    # Set the position for the scalebar (bottom-left corner)
    bar_thickness = int(im_height / 100)  # Height of the scalebar in pixels
    distance_from_border = int(im_width / 40)
    bottom_left_corner = (distance_from_border, im_height - distance_from_border)

    # Calculate the top-right corner of the rectangle
    top_right_corner = (bottom_left_corner[0] + bar_length_pixels,
                        bottom_left_corner[1] - bar_thickness)

    # Draw the rectangular scalebar
    cv2.rectangle(image, bottom_left_corner, top_right_corner, white_color, -1)  # Thickness of -1 fills the rectangle

    # Set position of label
    label_x = bottom_left_corner[0] + int(bar_length_pixels / 2) - 60
    label_y = top_right_corner[1] - bar_thickness * 2

    # Add the label for the scalebar
    cv2.putText(image, bar_label, (label_x, label_y), cv2.FONT_HERSHEY_SIMPLEX, 1.5, white_color, 2, cv2.LINE_AA)

    # cv2.imshow('Added scalebar', image)
    return image 
#%% Prompt    
class Prompt_User:
    """
    A simple GUI prompt using Tkinter to display a message and wait for user confirmation.

    This is useful for pausing execution and prompting the user to take action,
    such as selecting a position at the micriscope for manual EDS spectrum collection.

    Attributes
    ----------
    title : str
        The window title.
    message : str
        The message to display in the prompt.
    execution_stopped : bool
        True if the user closed the window or pressed Esc.
    ok_pressed : bool
        True if the user pressed OK or Return.
    root : tk.Tk or None
        The Tkinter root window.
    """
    def __init__(self, title: str, message: str):
        self.title = title
        self.message = message
        self.execution_stopped = False
        self.ok_pressed = False
        self.root = None

    def press_ok(self):
        """Handle OK button or Return key press."""
        self.ok_pressed = True
        self.root.quit()
        self.root.destroy()

    def stop_execution(self):
        """Handle window close (X) or Esc key press."""
        self.execution_stopped = True
        self.root.quit()
        self.root.destroy()

    def run(self):
        """
        Start the prompt window and wait for user interaction.
        Sets ok_pressed or execution_stopped depending on user action.
        """
        self.root = tk.Tk()
        self.root.title(self.title)

        # Create a label
        label = tk.Label(self.root, text=self.message)
        label.pack(pady=20)

        # Create an OK button
        ok_button = tk.Button(self.root, text="OK", command=self.press_ok)
        ok_button.pack(pady=20)

        # Bind the Return key to OK
        self.root.bind('<Return>', lambda event: self.press_ok())

        # Handle window close (X)
        self.root.protocol("WM_DELETE_WINDOW", self.stop_execution)

        # Bind the Esc key to stop execution
        self.root.bind('<Escape>', lambda event: self.stop_execution())

        # Start the main loop
        self.root.mainloop()
        
        
        
#%% Error handling

class RefLineError(Exception):
    """Exception raised for errors related to reference lines."""
    pass

class MissingHintError(Exception):
    """Exception raised when a required hint is missing."""
    pass

class EMError(Exception):
    """
    Custom exception class for electron microscope (EM)-related errors.

    Parameters
    ----------
    message : str
        Description of the error.
    code : int, optional
        Optional error code.
    """
    def __init__(self, message: str, code: Optional[int] = None) -> None:
        super().__init__(message)
        self.code = code


class EDSError(EMError):
    """
    Custom exception class for EDS-related errors.

    Parameters
    ----------
    message : str
        Description of the error.
    code : int, optional
        Optional error code.
    """
    def __init__(self, message: str, code: Optional[int] = None) -> None:
        super().__init__(message, code)