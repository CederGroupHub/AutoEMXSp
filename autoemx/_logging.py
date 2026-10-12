#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Centralized logging utilities for the autoemx package.

Library best-practice: the 'autoemx' logger carries only a NullHandler by
default (installed in __init__.py).  Users opt in by adding their own handlers,
e.g. ``logging.basicConfig(level=logging.INFO)``.
"""

import logging


def get_logger(name: str) -> logging.Logger:
    """Return a module-level logger.  Typically called as ``get_logger(__name__)``."""
    return logging.getLogger(name)
