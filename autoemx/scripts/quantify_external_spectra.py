#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quantify External Spectra with AutoEMX

Use this script when you have .msa, .msg, or .emsa spectra already collected outside
AutoEMX (e.g. exported from a third-party EDS software) and you want to
quantify them with the standard AutoEMX batch workflow.

This script copies spectra into the canonical folder layout, builds a
complete ledger.json, and immediately runs quantification and (optionally)
clustering analysis — all in a single call.

Typical usage
-------------
1. Edit the ``samples`` list below to specify each sample's ID, elements, and
   (optionally) the folder containing its spectra.  If ``spectra_dir`` is
   omitted, the spectra must already be present at
   ``<results_dir>/<sample_ID>/spectra/``.
2. Set the instrument and analysis parameters to match the acquisition
   conditions under which the spectra were collected.
3. Run this script.

Created on Fri Apr 24 2026

@author: Andrea
"""

# =============================================================================
# Sample Definitions
# =============================================================================

samples = [
    {
        # ---- Required ----
        'ID':  'Wulfenite_imported',           # Sample name; output folder is named after this
        'els': ['Pb', 'Mo', 'O'],              # Elements expected in the sample
        # 'spectra_dir': '/path/to/spectra',   # Omit if spectra are already in results_dir

        # ---- Optional ----
        'cnd': ['PbMoO4'],                     # Candidate phases for clustering
    },
    # Add more samples as needed, e.g.:
    # {
    #     'ID': 'K-412',
    #     'els': ['Mg', 'Al', 'Si', 'Ca', 'Fe', 'O'],
    #     'spectra_dir': '/path/to/K412_spectra',
    #     'cnd': [],
    # },
]

results_dir = None  # Set to an absolute path to save elsewhere; None → <cwd>/Results

# =============================================================================
# Sample Options
# =============================================================================

sample_type = 'bulk'     # 'powder' | 'bulk' | 'bulk_rough' | 'powder_continuous'
                         # Affects geometry correction factors during quantification.
                         # Use 'powder' for particles on a substrate,
                         # 'bulk' for a flat, polished surface.

sample_halfwidth = 2.9   # Half-width of the sample area in mm (informational only).

# =============================================================================
# Instrument Options
# =============================================================================
# Must match the acquisition conditions under which the spectra were collected.

# Microscope calibration folder must exist under autoemx/calibrations/<microscope_ID>/
microscope_ID    = 'PhenomXL'
measurement_mode = 'point'   # Selects detector-channel calibration and reference PB values

beam_energy = 15  # keV.  Critical: must match the beam energy used for acquisition.

# =============================================================================
# Substrate Options
# =============================================================================

els_substrate            = ['C', 'O', 'Al']  # Elements in the substrate; excluded from
                                              # quantification unless also in sample elements
sample_substrate_type    = 'Ctape'            # 'Ctape' or 'None'
sample_substrate_shape   = 'circle'           # 'circle' or 'square'
sample_substrate_width_mm = 12               # Stub lateral dimension in mm

# =============================================================================
# Quantification Options
# =============================================================================

use_project_specific_std_dict = False  # If True, loads standards from results_dir

interrupt_fits_bad_spectra = True  # Abort fits early when quality indicators are poor.

num_CPU_cores = None  # None → use half the available cores

use_instrument_background = False  # Use instrument background during fitting when available.

max_analytical_error = 5    # w%. Compositions above this threshold are excluded during clustering.

min_bckgrnd_cnts = 5        # Minimum counts under the reference peak for a spectrum to be
                             # accepted for quantification and clustering.

quant_flags_accepted = [0, -1]  # Quantification flags accepted during clustering.
                                   # See docs: Quantification flags (quant_flag).

max_n_clusters = 6          # Maximum number of clusters during compositional clustering.

show_unused_comps_clust = True  # Whether to plot discarded compositions in the clustering plot.

# =============================================================================
# Run
# =============================================================================

from autoemx.runners.quantify_external_spectra import quantify_external_spectra

quant_results = quantify_external_spectra(
    samples=samples,
    microscope_ID=microscope_ID,
    measurement_mode=measurement_mode,
    beam_energy=beam_energy,
    els_substrate=els_substrate,
    sample_substrate_type=sample_substrate_type,
    sample_substrate_shape=sample_substrate_shape,
    sample_substrate_width_mm=sample_substrate_width_mm,
    sample_type=sample_type,
    sample_halfwidth=sample_halfwidth,
    use_project_specific_std_dict=use_project_specific_std_dict,
    use_instrument_background=use_instrument_background,
    max_analytical_error_percent=max_analytical_error,
    min_bckgrnd_cnts=min_bckgrnd_cnts,
    quant_flags_accepted=quant_flags_accepted,
    max_n_clusters=max_n_clusters,
    show_unused_comps_clust=show_unused_comps_clust,
    interrupt_fits_bad_spectra=interrupt_fits_bad_spectra,
    num_CPU_cores=num_CPU_cores,
    results_path=results_dir,
)
