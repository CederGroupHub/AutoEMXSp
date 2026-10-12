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


@pytest.fixture(scope="session")
def _shared_autoemx_cache(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("autoemx_cache")


@pytest.fixture(autouse=True)
def autoemx_home(monkeypatch: pytest.MonkeyPatch, tmp_path_factory: pytest.TempPathFactory,
                 _shared_autoemx_cache: Path) -> Path:
    """Local AutoEMX folder (~/.autoemx) in a temporary folder, so tests never read or change the user's one.

    The cache (~/.autoemx/cache) is shared by all tests of the session, so detector convolution matrices are computed
    once: cached files are identified by the detector model they are computed from, so tests changing it get their own.
    """
    import autoemx._paths as paths

    home = tmp_path_factory.mktemp("autoemx_home")
    monkeypatch.setattr(paths, "AUTOEMX_HOME", home)

    def shared_cache_dir(microscope_ID: str) -> Path:
        path = _shared_autoemx_cache / microscope_ID
        path.mkdir(parents=True, exist_ok=True)
        return path

    monkeypatch.setattr(paths, "user_cache_dir", shared_cache_dir)
    return home
