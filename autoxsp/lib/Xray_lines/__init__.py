#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
X-ray Line Data Utilities
=========================

This module provides functionality to retrieve X-ray emission line energies 
and relative weights for chemical elements, with optional conversion to 
Siegbahn notation.

Data Sources:
    - Line energies: LineEnergies.csv
    - Line weights: LineWeights.csv
    - Line nomenclature conversion: lines_nomenclature_conversion.json

Dependencies:
    - pandas
    - pyatomgen
    - json
    - os

Author:
    Andrea Giunto
Date:
    2024-09-18
"""

import os
import json
from typing import Dict, Union

import pandas as pd
from pymatgen.core import Element


# -------------------------------------------------------------------------
# Locate data files relative to this script location
# -------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

LINE_ENERGIES_FILE = os.path.join(SCRIPT_DIR, "LineEnergies.csv")
LINE_WEIGHTS_FILE = os.path.join(SCRIPT_DIR, "LineWeights.csv")
LINE_NOMENCL_FILE = os.path.join(SCRIPT_DIR, "lines_nomenclature_conversion.json")


# -------------------------------------------------------------------------
# Load datasets
# -------------------------------------------------------------------------
def _load_csv(filepath: str) -> pd.DataFrame:
    """Load a CSV file into a pandas DataFrame, ensuring index starts at 1."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Required file not found: {filepath}")
    df = pd.read_csv(filepath)
    df.index = df.index + 1
    return df


def _load_json(filepath: str) -> dict:
    """Load a JSON file into a Python dictionary."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Required file not found: {filepath}")
    with open(filepath, "r", encoding="utf-8") as file:
        return json.load(file)


# Load the data files
LINE_ENERGIES_DF = _load_csv(LINE_ENERGIES_FILE)
LINE_WEIGHTS_DF = _load_csv(LINE_WEIGHTS_FILE)
LINE_NOMENCL_DICT = _load_json(LINE_NOMENCL_FILE)


# -------------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------------
def get_el_xray_lines(el: Union[str, int]) -> Dict[str, Dict[str, float]]:
    """
    Retrieve X-ray emission lines for a given element.

    Parameters
    ----------
    el : str | int
        The element symbol (e.g., "Fe") or its atomic number (e.g., 26).

    Returns
    -------
    Dict[str, Dict[str, float]]
        Dictionary mapping line names to sub-dictionaries containing:
            - "energy (keV)": float
                Line energy in keV.
            - "weight": float
                Relative intensity weight.

    Raises
    ------
    ValueError
        If `el` is neither a valid element symbol nor an integer atomic number.
    KeyError
        If the element is not found in the datasets.
    """
    # Determine atomic number using pyatomgen
    if isinstance(el, str):
        try:
            atomic_n = Element(el).Z
        except Exception:
            raise ValueError(f"Invalid element symbol: '{el}'")
    elif isinstance(el, int):
        if el < 1 or el > 118:
            raise ValueError(f"Invalid atomic number: {el}")
        atomic_n = el
    else:
        raise ValueError(
            f"Invalid type for 'el': {type(el)}. Must be str or int."
        )

    # Retrieve element data
    try:
        line_en_row = LINE_ENERGIES_DF.loc[atomic_n]
        line_w_vals = LINE_WEIGHTS_DF.loc[atomic_n].tolist()
    except KeyError:
        raise KeyError(f"No data found for atomic number {atomic_n}")

    # Combine energies and weights into a DataFrame
    line_w_row = pd.DataFrame([line_w_vals], columns=LINE_ENERGIES_DF.columns).loc[0]
    el_lines_df = pd.concat([line_en_row, line_w_row], axis=1)

    el_lines_dict = {}
    for index, row in el_lines_df.iterrows():
        line_en = row.iloc[0]  # Energy in eV
        line_w = row.iloc[1]   # Weight

        if line_en > 1 and line_w >= 0.0001:
            # Convert line name to Siegbahn notation if available
            line_name = LINE_NOMENCL_DICT.get(index, "") or index
            el_lines_dict[line_name] = {
                "energy (keV)": line_en / 1000,
                "weight": line_w
            }

    return el_lines_dict


# -------------------------------------------------------------------------
# Script execution example
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Example usage
    element = "Fe"
    lines = get_el_xray_lines(element)
    print(f"X-ray lines for {element}:\n{json.dumps(lines, indent=4)}")