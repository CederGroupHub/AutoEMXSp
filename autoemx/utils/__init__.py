#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 18:06:04 2025

@author: Andrea
"""

# NOTE:
# Avoid importing helper.py eagerly here because helper imports config schemas,
# and config initialization imports utils.constants, which would create a
# circular import chain at package import time.
from importlib import import_module


def __getattr__(name):
	_helper = import_module(f"{__name__}.helper")
	if name == "helper":
		return _helper
	return getattr(_helper, name)


def __dir__():
	_helper = import_module(f"{__name__}.helper")
	return sorted(set(globals().keys()) | set(dir(_helper)))