#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Particle segmentation with micro-sam (Segment Anything for Microscopy), each particle segmented on its own:

    A. Archit et al., "Segment Anything for Microscopy", Nature Methods 22, 579–591 (2025).
    https://doi.org/10.1038/s41592-024-02580-4, https://github.com/computational-cell-analytics/micro-sam

Automatic instance segmentation (AIS) of the frame with a micro-sam model (default ``vit_b_lm``, fine-tuned for
light microscopy, which separates touching objects well), then the filters of ``_instance_filters`` (substrate,
occluded particles). The brightness threshold is not used.

micro-sam is MIT-licensed; its model weights (CC-BY-4.0, by the micro-sam authors) are downloaded from Zenodo on
first use. Install (optional): ``python -m autoemx.utils.install_segmentation micro-sam`` (installs micro-sam
without PyQt6, the GPL Qt binding of its napari plugin, not needed here).

Options (``PowderMeasurementConfig.par_segmentation_params``), with their defaults:

- ``model_type``: ``"vit_b_lm"`` (other micro-sam models, e.g. ``"vit_l_lm"``, ``"vit_b_em_organelles"``);
- ``device``: ``"auto"`` (micro-sam's choice: CUDA, MPS or CPU), or ``"cpu"``, ``"cuda"``, ``"mps"``;
- ``tile_shape``, ``halo``: tiled segmentation of large frames (e.g. ``[512, 512]`` and ``[64, 64]``; default:
  whole frame);
- the filter options of ``_instance_filters.FILTER_DEFAULTS`` (e.g. ``min_visible_fraction``).

Created on Sat Oct 17 2026
@author: Andrea
"""

import importlib.util
from typing import Any, Dict

import numpy as np

from autoemx._logging import get_logger
from autoemx.core.em_runtime.particle_segmentation_models._instance_filters import filter_options, particle_labels

logger = get_logger(__name__)

INSTALL_HINT = "python -m autoemx.utils.install_segmentation micro-sam"
DEFAULT_MODEL_TYPE = "vit_b_lm"

_segmenters: Dict[tuple, Any] = {}  # (predictor, segmenter) loaded in this process, by (model, device, tiled)


def is_available() -> bool:
    """Whether micro-sam (and torch) is installed."""
    return all(importlib.util.find_spec(name) is not None for name in ("torch", "micro_sam"))


def get_predictor_and_segmenter(model_type: str = DEFAULT_MODEL_TYPE, device: str = "auto", is_tiled: bool = False):
    """The micro-sam predictor and AIS segmenter of a model (loaded once per process)."""
    key = (model_type, device, is_tiled)
    if key not in _segmenters:
        from micro_sam.automatic_segmentation import get_predictor_and_segmenter as load

        logger.info(f"Loading micro-sam model {model_type}…")
        _segmenters[key] = load(model_type=model_type, device=None if device == "auto" else device,
                                is_tiled=is_tiled)
    return _segmenters[key]


def instance_labels(frame_image: np.ndarray, params: Dict[str, Any]) -> np.ndarray:
    """Label image of the objects micro-sam finds in a greyscale frame."""
    from micro_sam.automatic_segmentation import automatic_instance_segmentation

    tile_shape, halo = params.get("tile_shape"), params.get("halo")
    predictor, segmenter = get_predictor_and_segmenter(params.get("model_type", DEFAULT_MODEL_TYPE),
                                                       params.get("device", "auto"), tile_shape is not None)
    gray = frame_image if frame_image.ndim == 2 else frame_image[..., 0]
    kwargs = {"tile_shape": tuple(tile_shape), "halo": tuple(halo or (64, 64))} if tile_shape is not None else {}
    return np.asarray(automatic_instance_segmentation(predictor=predictor, segmenter=segmenter,
                                                      input_path=np.ascontiguousarray(gray), ndim=2,
                                                      verbose=False, **kwargs))


def segment_particles(frame_image: np.ndarray,
                      powder_meas_config: 'PowderMeasurementConfig' = None,
                      save_image: bool = False,
                      EM: 'EM_controller' = None) -> np.ndarray:
    """
    Label image of the particles of a frame segmented by micro-sam (int32, background 0, one label per particle).

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
        raise ImportError(f"micro-sam is not installed: {INSTALL_HINT}")
    params = dict(powder_meas_config.par_segmentation_params) if powder_meas_config is not None else {}
    instances = instance_labels(frame_image, params)
    labels = particle_labels(instances, frame_image, **filter_options(params))
    logger.debug(f"micro-sam: {int(instances.max())} objects, {int(labels.max())} particles kept.")
    if save_image and EM:
        EM.save_frame_image(f"{EM.sample_id}_fr{EM.current_frame_label}_ml_mask",
                            frame_image=np.where(labels > 0, 255, 0).astype(np.uint8))
    return labels
