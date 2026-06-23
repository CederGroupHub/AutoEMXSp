"""Backward-compatibility shim for the renamed package.

This package was renamed from ``autoemxsp`` to ``autoemx``.
Legacy imports are kept working temporarily and emit a warning.
"""

from __future__ import annotations

import importlib
import sys
import warnings

from autoemx import *  # noqa: F401,F403
import autoemx as _autoemx
from autoemx._compat import warn_if_stale_autoemxsp_install

warnings.warn(
    "Package 'autoemxsp' has been renamed to 'autoemx'. "
    "Please update imports and use: pip install autoemx",
    FutureWarning,
    stacklevel=2,
)

warn_if_stale_autoemxsp_install()

# Reuse the new package path so imports like 'autoemxsp.runners' resolve.
__path__ = _autoemx.__path__

# Register top-level subpackage aliases to support imports like:
# from autoemxsp.runners.analyze_sample import analyze_sample
for _subpkg in (
    "calibrations",
    "config",
    "core",
    "data",
    "microscope_drivers",
    "runners",
    "utils",
):
    sys.modules[f"{__name__}.{_subpkg}"] = importlib.import_module(f"autoemx.{_subpkg}")
