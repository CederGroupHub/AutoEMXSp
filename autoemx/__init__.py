#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 13:30:06 2025

@author: Andrea
"""

from ._version import __version__

import logging as _logging
# Library best-practice: ship a NullHandler so users never see
# "No handlers could be found for logger 'autoemx'" warnings.
_logging.getLogger("autoemx").addHandler(_logging.NullHandler())
# Route warnings.warn() calls into the logging system.
_logging.captureWarnings(True)

# After captureWarnings, so that the warning also reaches the log
from ._compat import warn_if_stale_autoemxsp_install
warn_if_stale_autoemxsp_install()

