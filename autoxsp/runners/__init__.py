#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:30:06 2025

@author: Andrea
"""

import pkgutil
import importlib
from pathlib import Path

# Automatically import all modules in the current package
package_dir = Path(__file__).parent
for _, module_name, _ in pkgutil.iter_modules([str(package_dir)]):
    module = importlib.import_module(f".{module_name}", package=__name__)
    # Import all functions from the module into this namespace
    for attr in dir(module):
        if not attr.startswith("_"):
            globals()[attr] = getattr(module, attr)
