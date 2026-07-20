#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared fixtures for CI-safe AutoEMX tests."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from ci_paths import (
    INPUTS_DIR,
    K412_CLUSTER_MINI_ID,
    OUTPUTS_DIR,
    WULFENITE_MINI_ID,
)

# Non-interactive plotting for headless CI runners.
os.environ.setdefault("MPLBACKEND", "Agg")


@pytest.fixture
def inputs_dir() -> Path:
    return INPUTS_DIR


@pytest.fixture
def outputs_dir() -> Path:
    OUTPUTS_DIR.mkdir(exist_ok=True)
    return OUTPUTS_DIR


@pytest.fixture
def results_path(tmp_path: Path) -> str:
    """Copy mini sample ledgers into a temp results root for isolation."""
    for sample_id in (WULFENITE_MINI_ID, K412_CLUSTER_MINI_ID):
        src = INPUTS_DIR / sample_id
        dst = tmp_path / sample_id
        shutil.copytree(src, dst)
    return str(tmp_path)


@pytest.fixture(autouse=True)
def _patch_opencv_gui(monkeypatch: pytest.MonkeyPatch):
    """Prevent OpenCV GUI calls from hanging headless CI."""
    try:
        import cv2
    except ImportError:
        return

    monkeypatch.setattr(cv2, "imshow", lambda *args, **kwargs: None, raising=False)
    monkeypatch.setattr(cv2, "waitKey", lambda *args, **kwargs: 1, raising=False)
    monkeypatch.setattr(cv2, "destroyAllWindows", lambda *args, **kwargs: None, raising=False)
