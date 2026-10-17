#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle segmentation with SAM2 (Segment Anything Model 2.1), each particle segmented on its own, so that touching
particles are measured separately:

    N. Ravi et al., "SAM 2: Segment Anything in Images and Videos", arXiv:2408.00714 (2024).
    https://github.com/facebookresearch/sam2

SAM2 segments the whole frame (automatic mask generation, from a grid of points), through Hugging Face
``transformers`` with Meta's weights (``facebook/sam2.1-hiera-*``, Apache-2.0, downloaded from the Hugging Face Hub
on first use). The masks are then filtered (substrate, agglomerates, surface features, occluded particles; see
``_instance_filters``). The brightness threshold is not used.

Install (optional): ``pip install "autoemx[segmentation]"`` (torch, torchvision, transformers), or, for an editable
installation, ``pip install -e ".[segmentation]"`` from the repository folder.

Options (``PowderMeasurementConfig.par_segmentation_params``), with their defaults:

- ``checkpoint``: ``"facebook/sam2.1-hiera-large"`` (or ``-base-plus``, ``-small``, ``-tiny``: faster, less accurate);
- ``device``: ``"auto"`` (CUDA if available, else CPU), or a torch device such as ``"cpu"``, ``"cuda"``, ``"mps"``;
- ``points_per_batch`` (64), ``points_per_crop`` (32, grid of 32 × 32 points), ``crops_n_layers`` (0),
  ``pred_iou_thresh`` (0.8), ``stability_score_thresh`` (0.9), ``crops_nms_thresh`` (0.7): options of the
  transformers mask-generation pipeline;
- the filter options of ``_instance_filters.FILTER_DEFAULTS`` (e.g. ``min_visible_fraction``).

Created on Sat Oct 17 2026
@author: Andrea
"""

import importlib.util
from pathlib import Path
from typing import Any, Dict

import numpy as np

from autoemx._logging import get_logger
from autoemx.core.em_runtime.particle_segmentation_models._instance_filters import filter_options, particle_labels

logger = get_logger(__name__)



def _install_hint() -> str:
    """The pip command installing SAM2: from the repository folder for an editable installation of AutoEMX (a PyPI
    install would replace it), else from PyPI."""
    repo = Path(__file__).resolve().parents[4]
    if (repo / "pyproject.toml").exists():
        return f'pip install -e "{repo}[segmentation]"'
    return 'pip install "autoemx[segmentation]"'


INSTALL_HINT = _install_hint()

#: Options of the mask generation and their defaults
GENERATOR_DEFAULTS: Dict[str, Any] = {
    "points_per_batch": 64,
    "points_per_crop": 32,
    "crops_n_layers": 0,
    "pred_iou_thresh": 0.8,
    "stability_score_thresh": 0.9,
    "crops_nms_thresh": 0.7,
}
DEFAULT_CHECKPOINT = "facebook/sam2.1-hiera-large"

_generators: Dict[tuple, Any] = {}  # mask generators loaded in this process, by (checkpoint, device)


def is_available() -> bool:
    """Whether torch, torchvision and transformers are installed."""
    return all(importlib.util.find_spec(name) is not None for name in ("torch", "torchvision", "transformers"))


def _device(name: str) -> str:
    if name != "auto":
        return name
    import torch

    return "cuda" if torch.cuda.is_available() else "cpu"


def get_mask_generator(checkpoint: str = DEFAULT_CHECKPOINT, device: str = "auto"):
    """The transformers mask-generation pipeline of a SAM2 checkpoint (loaded once per process)."""
    key = (checkpoint, _device(device))
    if key not in _generators:
        from transformers import pipeline

        logger.info(f"Loading SAM2 checkpoint {checkpoint} on {key[1]} (its weights, up to ~900 MB, are downloaded "
                    f"on first use)…")
        _generators[key] = pipeline("mask-generation", model=checkpoint, device=key[1])
    return _generators[key]


def generate_masks(frame_image: np.ndarray, params: Dict[str, Any]) -> list:
    """Masks of every object SAM2 finds in a greyscale frame (bool arrays of the frame's shape)."""
    import torch
    from PIL import Image

    generator = get_mask_generator(params.get("checkpoint", DEFAULT_CHECKPOINT), params.get("device", "auto"))
    options = {key: params.get(key, default) for key, default in GENERATOR_DEFAULTS.items()}
    gray = frame_image if frame_image.ndim == 2 else frame_image[..., 0]
    with torch.inference_mode():
        output = generator(Image.fromarray(np.ascontiguousarray(gray)).convert("RGB"), **options)
    return [np.asarray(m.cpu() if hasattr(m, "cpu") else m, dtype=bool) for m in output["masks"]]


def segment_particles(frame_image: np.ndarray,
                      powder_meas_config: 'PowderMeasurementConfig' = None,
                      save_image: bool = False,
                      EM: 'EM_controller' = None) -> np.ndarray:
    """
    Label image of the particles of a frame segmented by SAM2 (int32, background 0, one label per particle).

    Parameters
    ----------
    frame_image : ndarray
        Greyscale frame image.
    powder_meas_config : PowderMeasurementConfig, optional
        Its ``par_segmentation_params`` are the options of the model (see the module documentation).
    save_image : bool, optional
        Save the label image as a mask through ``EM`` (default: False).
    EM : EM_controller, optional
        Used to save the mask image.
    """
    if not is_available():
        raise ImportError(f"SAM2 needs torch, torchvision and transformers: {INSTALL_HINT}")
    params = dict(powder_meas_config.par_segmentation_params) if powder_meas_config is not None else {}
    masks = generate_masks(frame_image, params)
    labels = particle_labels(masks, frame_image, **filter_options(params))
    logger.debug(f"SAM2: {len(masks)} masks, {int(labels.max())} particles kept.")
    if save_image and EM:
        EM.save_frame_image(f"{EM.sample_id}_fr{EM.current_frame_label}_ml_mask",
                            frame_image=np.where(labels > 0, 255, 0).astype(np.uint8))
    return labels
