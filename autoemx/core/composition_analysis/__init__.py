#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Composition analysis package.

This package encapsulates all composition analysis functionality including:
- Clustering of compositional data
- Visualization and plotting
- Reference material standards handling
- Phase identification and reference matching
- Main analysis orchestration

Modules
-------
- clustering: K-means and DBSCAN clustering algorithms
- plotting: Visualization of analysis results
- standards: Reference material and correction factor management
- ledger: Sample ledger (ledger.json) of the analyzer: spectra, particles, quantification runs, clustering analyses
- reference_matching: Phase identification and correlation
- analyser: Main EMXSp_Composition_Analyzer orchestrator

Public API
----------
EMXSp_Composition_Analyzer: Main analysis class

Examples
--------
    >>> from autoemx.core.composition_analysis import EMXSp_Composition_Analyzer
    >>> analyzer = EMXSp_Composition_Analyzer(...)
    >>> analyzer.run_collection_and_quantification()

    >>> # Sample already acquired
    >>> analyzer = EMXSp_Composition_Analyzer.from_ledger('Results/sample_001')

@author: Andrea
Created during refactoring of EMXSp_Composition_Analyzer
"""

from .clustering import ClusteringModule
from .plotting import PlottingModule
from .reference_matching import ReferenceMatchingModule
from .standards import StandardsModule
from .analyser import EMXSp_Composition_Analyzer

__all__ = [
    'EMXSp_Composition_Analyzer',
    'ClusteringModule',
    'PlottingModule',
    'ReferenceMatchingModule',
    'StandardsModule',
]
