#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility Functions for EDS Spectrum Analysis and Data Handling
==================================================================

A collection of general-purpose utility functions and lightweight classes for data handling,
visualization, and file management within the EDS analysis and modeling framework.

This module is designed to be self-contained and importable across different parts
of the project, covering tasks from compositional conversions to structured file
loading and formatted console output.

-------------------------------------------------------------------------------

Main Features
--------------

**Compositional Conversions**
- `atomic_to_weight_fr()`: Convert atomic fractions to weight fractions.
- `weight_to_atomic_fr()`: Convert weight fractions to atomic fractions.
  Both functions use `pymatgen.Element` for accurate atomic mass values.

**Formula Handling**
- `get_std_comp_from_formula()`: Parse chemical formulas into standardized element dictionaries.
- `to_latex_formula()`: Convert a chemical formula into LaTeX format (e.g., ``Fe2O3`` → ``Fe$_2$O$_3$``).

**String and Table Utilities**
- `print_nice_1d_row()`: Print a formatted 1D table row (with adjustable width and alignment).
- `print_element_fractions_table()`: Display element names with their atomic and weight percentages.
- `print_single_separator()`: Print a single horizontal line separator.
- `print_double_separator()`: Print a double horizontal line separator for clear sectioning.
- `AlphabetMapper`: Map integers to alphabetic letters (0 → 'A', 25 → 'Z', etc.).

**Image Utilities**
- `draw_scalebar()`: Draw a scale bar on an image using OpenCV, with optional text label and styling.

**File and Directory Handling**
- `get_sample_dir()`: Locate a directory named after a given sample ID, with optional recursive search.
- `make_unique_path()`: Generate a unique file or directory path if one already exists (e.g., ``file_1.png``).
- `load_msa()`: Load and parse a `.msa` spectral file, returning energy and intensity arrays.

**User Interaction**
- `Prompt_User`: Display a Tkinter prompt window for manual user confirmation (e.g., proceed/stop execution).

**Error Handling**
- `EDSError`: Custom exception class for handling EDS-related errors gracefully.

Created on Fri Jun 28 11:50:53 2024

@author: Andrea Giunto
"""

# Standard library imports
import os
import sys
import re
import logging
import warnings
import unicodedata
import difflib

# Third-party imports
import cv2
import numpy as np
from pathlib import Path
from pymatgen.core import Element, Composition

# Typing imports
from typing import Any, Sequence, List, Optional, Tuple, Dict, Union


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
        include_dollar_signs (bool): Whether to wrap the LaTeX in $...$.

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
                group = s[i + 1:j - 1]  # Content inside the parentheses

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
    
    

def _safe_autoemx_log_info(msg: str) -> None:
    """Log to autoemx logger; if stream handlers break, fall back to stdout."""
    try:
        logging.getLogger("autoemx").info(msg)
    except (KeyboardInterrupt, OSError, ValueError):
        # Some Windows environments can surface stream/handler issues as KeyboardInterrupt
        # while writing logs (without an actual user Ctrl+C). Keep the workflow running.
        try:
            print(msg, flush=True)
        except Exception:
            pass


def print_single_separator(message: Optional[str] = None):
    """Log a single-line separator and optionally print a message right below it."""
    _safe_autoemx_log_info('-' * 50)
    if message:
        print(message)




def print_double_separator(message: Optional[str] = None):
    """Log a double-line separator (50 equals signs) for visual clarity."""
    _safe_autoemx_log_info('=' * 50)
    if message:
        print(message)


class AlphabetMapper:
    """
    Maps a zero-based integer index to Excel-style column letters.
    0 -> 'A', 1 -> 'B', ..., 25 -> 'Z', 26 -> 'AA', etc.
    
    Used for labeling frames analyzed during particle search.
    Pass a lowercase alphabet for compact child-frame labels in multi-scale scans.
    """
    def __init__(self, alphabet: str = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'):
        if not alphabet:
            raise ValueError("alphabet must be non-empty")
        self.alphabet = alphabet
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


def make_unique_path(
    parent_dir: Union[str, Path],
    base_name: Union[str, Path],
    extension: Optional[str] = None
) -> Path:
    """
    Generate a unique file or directory path inside parent_dir based on base_name.
    If a path with the base name exists, appends a counter (e.g., 'Sample1_2').
    Optionally, add an extension for files.

    Parameters
    ----------
    parent_dir : str or Path
        The parent directory in which to generate the new path.
    base_name : str or Path
        The base name for the new file or directory.
    extension : str, optional
        The file extension (e.g., 'txt' or '.txt'). If None, treat as directory.

    Returns
    -------
    unique_path : Path
        The full, unique path (not created).

    Raises
    ------
    ValueError
        If `parent_dir` or `base_name` is invalid.
    """

    if not isinstance(parent_dir, (str, Path)):
        raise ValueError("`parent_dir` must be a string or Path.")
    if not isinstance(base_name, (str, Path)):
        raise ValueError("`base_name` must be a string or Path.")

    parent_dir = Path(parent_dir)
    base_name = Path(base_name).name  # Ensure only the name is used

    # Normalize extension
    ext = ''
    if extension:
        ext = extension if extension.startswith('.') else f'.{extension}'

    counter = 1
    if ext:
        unique_path = parent_dir / f"{base_name}{ext}"
    else:
        unique_path = parent_dir / base_name

    while unique_path.exists():
        counter += 1
        if ext:
            unique_path = parent_dir / f"{base_name}_{counter}{ext}"
        else:
            unique_path = parent_dir / f"{base_name}_{counter}"

    return unique_path


def get_sample_dir(
    results_path: str,
    sample_ID: str,
    case_insensitive: bool = True,
    verbose: bool = False,
) -> str:
    """
    Find a directory named `sample_ID` under `results_path` or its subdirectories.

    Strategy:
      1. Walk the entire directory tree under results_path and collect exact matches.
      2. If multiple matches -> raise RuntimeError (ambiguous).
      3. If none -> show close matches and raise FileNotFoundError.

    Parameters
    ----------
    results_path : str
        Root directory to search from.
    sample_ID : str
        Directory name to search for (exact match).
    case_insensitive : bool, optional
        Whether to ignore case when matching. Default is True.
    verbose : bool, optional
        Print additional debug information. Default is False.

    Returns
    -------
    sample_dir : str
        Full path to the matched directory.

    Raises
    ------
    RuntimeError
        If multiple matches are found.
    FileNotFoundError
        If no match is found.
    """
    results_path = os.path.abspath(os.path.expanduser(results_path))
    if verbose:
        print(f"[get_sample_dir] Searching for '{sample_ID}' under '{results_path}'", file=sys.stderr)

    if not os.path.isdir(results_path):
        raise FileNotFoundError(f"Search root does not exist or is not a directory: '{results_path}'")

    def _norm_name(s: str) -> str:
        """Normalize unicode and collapse case for robust name matching."""
        return unicodedata.normalize("NFKC", s).strip()

    sample_norm = _norm_name(sample_ID)

    def match_name(name: str) -> bool:
        name_norm = _norm_name(name)
        return name_norm.lower() == sample_norm.lower() if case_insensitive else name_norm == sample_norm

    # Fast-path: direct child under results_path.
    direct_candidate = os.path.join(results_path, sample_ID)
    if os.path.isdir(direct_candidate):
        return os.path.abspath(direct_candidate)

    # One-level check (case-insensitive if requested) before deeper traversal.
    try:
        with os.scandir(results_path) as it:
            for entry in it:
                if entry.is_dir(follow_symlinks=False) and match_name(entry.name):
                    return os.path.abspath(entry.path)
    except PermissionError:
        pass

    # Skip known expensive/irrelevant directory names while walking.
    skip_dir_names = {
        ".git", "__pycache__", ".ipynb_checkpoints", ".DS_Store",
        "node_modules", "Library", "Applications", "System",
    }

    matches: List[str] = []
    seen_names = set()
    verbose_listed_paths: List[str] = []

    for root, dirs, _ in os.walk(results_path, topdown=True, followlinks=False):
        # In-place prune: hidden dirs and known large/system trees.
        dirs[:] = [
            d for d in dirs
            if d not in skip_dir_names and not d.startswith('.')
        ]

        for d in dirs:
            seen_names.add(d)
            if verbose and len(verbose_listed_paths) < 50:
                verbose_listed_paths.append(os.path.join(root, d))
            if match_name(d):
                matches.append(os.path.abspath(os.path.join(root, d)))
                # We only need to know if there is ambiguity, not all possible matches.
                if len(matches) > 1:
                    break
        if len(matches) > 1:
            break

    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError(f"Ambiguous '{sample_ID}' — found in multiple locations:\n" + "\n".join(matches))

    close = difflib.get_close_matches(sample_ID, sorted(seen_names), n=5, cutoff=0.6)
    debug_msg_lines = [
        f"'{sample_ID}' folder not found in '{results_path}' or its subdirectories."
    ]
    if close:
        debug_msg_lines.append(f"Did you mean one of: {close}?")
    else:
        debug_msg_lines.append("No similar directory names found.")

    if verbose and verbose_listed_paths:
        debug_msg_lines.append("\nSome directories under the search root (first 50):")
        debug_msg_lines.extend(verbose_listed_paths)

    raise FileNotFoundError("\n".join(debug_msg_lines))
        

def load_msa(filepath: str) -> Tuple[np.ndarray, np.ndarray, Dict[str, str]]:
    """
    Load an EMSA/MAS spectrum file (``.msa``, ``.msg``, or ``.emsa``) containing
    Y-only spectral data and compute the energy scale.
    Designed for raw spectra exported from Thermo Fisher Phenom systems.
    May work with other EMSA/MAS format files, but minor variations can occur.

    Parameters
    ----------
    filepath : str
        Path to the ``.msa``, ``.msg``, or ``.emsa`` file.

    Returns
    -------
    energy : np.ndarray
        Energy values computed from OFFSET and XPERCHAN (in eV).
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

            if line.startswith('#SPECTRUM') or line.startswith('##SPECTRUM'):
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
                # Handle both single and double column formats
                # Removing trailing commas before splitting
                line = line.rstrip(',')
                if ',' in line:
                    parts = line.split(',', maxsplit=1)
                    token = parts[1].strip() if len(parts) > 1 else parts[0].strip()
                else:
                    token = line.strip()
                    
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
def draw_scalebar(image, pixel_size_um, bar_width = 0.25, color = None):
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
    if color is None:
        color = white_color
    cv2.rectangle(image, bottom_left_corner, top_right_corner, color, -1)  # Thickness of -1 fills the rectangle

    # Set position of label
    label_x = bottom_left_corner[0] + int(bar_length_pixels / 2) - 60
    label_y = top_right_corner[1] - bar_thickness * 2

    # Add the label for the scalebar
    cv2.putText(image, bar_label, (label_x, label_y), cv2.FONT_HERSHEY_SIMPLEX, 1.5, color, 2, cv2.LINE_AA)

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
        if self.root:
            self.root.quit()
            self.root.destroy()

    def stop_execution(self):
        """Handle window close (X) or Esc key press."""
        self.execution_stopped = True
        if self.root:
            self.root.quit()
            self.root.destroy()

    def run(self):
        """
        Start the prompt window and wait for user interaction.
        Sets ok_pressed or execution_stopped depending on user action.
        """
        import tkinter as tk
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
