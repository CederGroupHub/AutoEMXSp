#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publish experimental standards to the package.

Copies the standards file of your instrument (~/.autoemx/calibrations/<microscope_ID>/EDS_Stds_<beam_energy>keV.json)
into the package (autoemx/calibrations/<microscope_ID>/), and regenerates its details report there. Commit both files
to distribute the new values as the initial standards of the microscope.

Only runs in a git checkout of AutoEMX.

Created on Mon Oct 12 2026

@author: Andrea
"""

# =============================================================================
# Input
# =============================================================================

microscope_ID = 'PhenomXL'
beam_energy = 15  # keV

# =============================================================================
# Run
# =============================================================================

from autoemx.runners.publish_standards import publish_standards

publish_standards(microscope_ID=microscope_ID, beam_energy=beam_energy)
