#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optional particle segmentation models giving one mask per particle (SAM2, micro-sam): filters of their masks, label
images in the particle finder, registry of installed models, the plugins (with fake models, no download), each
model's result chosen and its frames and particles left out, and the install helper. The real models are run only
where they are installed with their weights.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import cv2
import numpy as np
import pytest

import autoemx.core.em_runtime.particle_segmentation_models as par_seg_models
from autoemx.config.ledger_io import load_sample_ledger
from autoemx.config.ledger_schemas import SegmentationConfig
from autoemx.config.runtime_configs import PowderMeasurementConfig
from autoemx.core.em_runtime.particle_finder import is_binary_mask, label_image_stats
from autoemx.core.em_runtime.particle_segmentation_models import SAM2, microSAM
from autoemx.core.em_runtime.particle_segmentation_models._instance_filters import particle_labels
from autoemx.runners.process_particle_statistics import list_size_results, process_particle_statistics
from autoemx.runners.reanalyze_particle_statistics import reanalyze_particle_statistics
from autoemx.utils import install_segmentation
from _psd_helpers import make_measurement_sample

H, W = 200, 300


def _disk(cx: int, cy: int, r: int) -> np.ndarray:
    mask = np.zeros((H, W), np.uint8)
    cv2.circle(mask, (cx, cy), r, 1, -1)
    return mask.astype(bool)


def _labels_of(labels: np.ndarray, mask: np.ndarray) -> set:
    return set(np.unique(labels[mask]).tolist()) - {0}


# =============================================================================
# Filters of the masks
# =============================================================================
@pytest.fixture
def frame():
    """Bright particles on a dark substrate, and the masks a promptable model would give."""
    image = np.full((H, W), 40, np.uint8)
    masks = {
        "grain": _disk(60, 60, 30), "pit": _disk(60, 60, 5),  # a pit on a grain
        "a1": _disk(175, 60, 15), "a2": _disk(204, 60, 15),  # two touching particles...
        "front": _disk(100, 150, 15), "behind": _disk(112, 150, 16),  # a particle partly hidden behind another
        "dark": _disk(250, 160, 15),  # a dark region of the substrate
    }
    masks["agglomerate"] = masks["a1"] | masks["a2"]  # ...and their agglomerate
    masks["substrate"] = np.ones((H, W), bool)
    for name in ("grain", "a1", "a2", "front", "behind"):
        image[masks[name]] = 200
    return image, masks


def test_masks_filtered(frame):
    image, masks = frame
    labels = particle_labels(list(masks.values()), image)
    assert labels.dtype == np.int32
    assert sorted(np.unique(labels).tolist()) == list(range(labels.max() + 1))  # consecutive labels
    assert len(_labels_of(labels, masks["grain"])) == 1  # the pit is a surface feature: dropped
    assert _labels_of(labels, masks["pit"]) == _labels_of(labels, masks["grain"])
    # The agglomerate is dropped, its particles kept
    a1, a2 = _labels_of(labels, masks["a1"] & ~masks["a2"]), _labels_of(labels, masks["a2"] & ~masks["a1"])
    assert len(a1) == len(a2) == 1 and a1 != a2
    # The particle hidden behind another is dropped; the substrate and the dark region too
    assert _labels_of(labels, masks["front"]) and not _labels_of(labels, masks["behind"] & ~masks["front"])
    assert not _labels_of(labels, masks["dark"])
    assert labels.max() == 4  # grain, a1, a2, front


def test_filter_options(frame):
    image, masks = frame
    everything = list(masks.values())
    # Occluded particles kept with a lower minimum visible fraction; dark masks kept if asked
    assert _labels_of(particle_labels(everything, image, min_visible_fraction=0.3), masks["behind"] & ~masks["front"])
    assert _labels_of(particle_labels(everything, image, drop_dark_masks=False), masks["dark"])
    # Irregular particles dropped with a minimum solidity
    crescent = _disk(250, 60, 25) & ~_disk(265, 60, 22)
    image[crescent] = 200
    assert _labels_of(particle_labels([crescent], image), crescent)
    assert not _labels_of(particle_labels([crescent], image, min_solidity=0.8), crescent)


def test_label_image_input(frame):
    image, masks = frame
    instances = np.zeros((H, W), np.int32)
    for i, name in enumerate(("grain", "a1", "a2", "dark"), start=1):
        instances[masks[name]] = i
    instances[_disk(60, 60, 5)] = 9  # a pit, as a separate instance (label 9: labels need not be consecutive)
    labels = particle_labels(instances, image)
    assert labels.max() == 3 and not _labels_of(labels, masks["dark"])
    assert _labels_of(labels, masks["pit"]) == _labels_of(labels, masks["grain"])


# =============================================================================
# Label images in the particle finder
# =============================================================================
def test_label_image_stats_as_connected_components():
    rng = np.random.default_rng(0)
    labels = np.zeros((120, 160), np.int32)
    for i in range(1, 40):
        y, x = rng.integers(0, 110), rng.integers(0, 150)
        labels[y:y + rng.integers(1, 10), x:x + rng.integers(1, 10)] = i
    labels[labels == 7] = 0  # a label absent from the image
    num_labels, stats, centroids = label_image_stats(labels)
    assert num_labels == 40
    for label in range(num_labels):
        ys, xs = np.nonzero(labels == label)
        if not ys.size:
            assert not stats[label].any()
            continue
        assert stats[label].tolist() == [xs.min(), ys.min(), xs.max() - xs.min() + 1, ys.max() - ys.min() + 1, ys.size]
        assert centroids[label] == pytest.approx([xs.mean(), ys.mean()])


def test_binary_masks_and_label_images_told_apart():
    assert is_binary_mask(np.array([[0, 255]], np.uint8)) and is_binary_mask(np.array([[True]]))
    assert is_binary_mask(np.zeros((2, 2), np.uint8))
    assert not is_binary_mask(np.array([[0, 255]], np.int32))  # two particles labelled 0 and 255, not binary
    assert not is_binary_mask(np.array([[0, 3]], np.uint8))


# =============================================================================
# Registry of the models
# =============================================================================
def test_models_registered_even_if_not_installed(monkeypatch):
    monkeypatch.setattr(SAM2, "is_available", lambda: False)
    assert "SAM2" in par_seg_models.AVAILABLE_SEGMENTATION_MODELS
    assert "SAM2" not in par_seg_models.available_models()
    # A configuration (e.g. of a ledger) recorded with the model stays valid
    assert PowderMeasurementConfig(par_segmentation_model="SAM2").par_segmentation_model == "SAM2"
    with pytest.raises(ImportError, match="\\[segmentation\\]"):
        par_seg_models.require_model("SAM2")
    monkeypatch.setattr(SAM2, "is_available", lambda: True)
    assert "SAM2" in par_seg_models.available_models()
    assert par_seg_models.require_model("SAM2") is SAM2


def test_model_options_part_of_the_segmentation():
    cfg = PowderMeasurementConfig(par_segmentation_model="SAM2", par_segmentation_params={"pred_iou_thresh": 0.7})
    segmentation = SegmentationConfig.from_powder_meas_cfg(cfg, 0)
    assert segmentation.model_params == {"pred_iou_thresh": 0.7}
    assert segmentation.label == "SAM2 · 80 px · pred_iou_thresh=0.7"
    default = SegmentationConfig.from_powder_meas_cfg(PowderMeasurementConfig(par_segmentation_model="SAM2"), 1)
    assert default.fingerprint() != segmentation.fingerprint()


# =============================================================================
# The plugins, with fake models
# =============================================================================
def _component_masks(image: np.ndarray) -> list:
    """What a promptable model would find: each bright object, the substrate, and an agglomerate of two objects."""
    n, labels = cv2.connectedComponents((image > 100).astype(np.uint8))
    masks = [labels == i for i in range(1, n)]
    return masks + [np.ones(image.shape, bool)] + ([masks[0] | masks[1]] if len(masks) > 1 else [])


@pytest.fixture
def fake_models(monkeypatch):
    """SAM2 and micro-sam 'installed', with fake models; the options they were given."""
    calls = {"SAM2": [], "microSAM": []}

    def sam2_masks(image, params):
        calls["SAM2"].append(dict(params))
        return _component_masks(image)

    def msam_instances(image, params):
        calls["microSAM"].append(dict(params))
        return cv2.connectedComponents((image > 100).astype(np.uint8))[1]

    for module, name, fake in ((SAM2, "generate_masks", sam2_masks), (microSAM, "instance_labels", msam_instances)):
        monkeypatch.setattr(module, "is_available", lambda: True)
        monkeypatch.setattr(module, name, fake)
    return calls


@pytest.mark.parametrize("model", ["SAM2", "microSAM"])
def test_plugin_label_image(fake_models, model):
    image = np.full((H, W), 40, np.uint8)
    for cx in (50, 120, 190):
        image[_disk(cx, 100, 20)] = 200
    cfg = PowderMeasurementConfig(par_segmentation_model=model, par_segmentation_params={"checkpoint": "x"})
    labels = par_seg_models.require_model(model).segment_particles(image, cfg)
    assert labels.dtype == np.int32 and labels.max() == 3
    assert fake_models[model] == [{"checkpoint": "x"}]


def test_plugin_not_installed(monkeypatch):
    monkeypatch.setattr(microSAM, "is_available", lambda: False)
    with pytest.raises(ImportError, match="install_segmentation micro-sam"):
        microSAM.segment_particles(np.zeros((10, 10), np.uint8), PowderMeasurementConfig())


# =============================================================================
# One result per model; frames and particles of the chosen result left out
# =============================================================================
@pytest.fixture
def sample(tmp_path: Path) -> Path:
    make_measurement_sample(tmp_path / "results" / "S", images=True)
    return tmp_path / "results"


def test_each_model_gives_its_result_and_its_particles_are_left_out(sample: Path, fake_models, capsys):
    reanalyze_particle_statistics("S", str(sample), par_segmentation_model="SAM2", verbose=False)
    reanalyze_particle_statistics("S", str(sample), par_segmentation_model="microSAM",
                                  par_segmentation_params={"model_type": "vit_l_lm"}, verbose=False)
    ledger = load_sample_ledger(str(sample / "S" / "ledger.json"))
    models = {s.segmentation_id: (s.par_segmentation_model, s.model_params) for s in ledger.segmentations}
    assert sorted(models.values(), key=str) == sorted(
        [("threshold_bright", {}), ("SAM2", {}), ("microSAM", {"model_type": "vit_l_lm"})], key=str)
    assert fake_models["microSAM"][0] == {"model_type": "vit_l_lm"}
    # A result per model (none combined), each from its own segmentation of every frame
    sam2_result, msam_result = ledger.size_results[1], ledger.size_results[2]
    assert models[sam2_result.segmentation_id][0] == "SAM2" and models[msam_result.segmentation_id][0] == "microSAM"
    frame = ledger.frames_of("size_stats", 0)[0]
    assert {s.segmentation_id for s in frame.segmentations} == set(models)
    assert sam2_result.size_stats.n_par_analysed > 0

    # The user lists the results, chooses the SAM2 one, and leaves out a frame and a particle of it
    labels = list_size_results("S", str(sample))
    assert len(labels) == 3 and "SAM2" in labels[1] and "SAM2" in capsys.readouterr().out
    other = ledger.frames_of("size_stats", 0)[1]
    process_particle_statistics("S", str(sample), result_id=sam2_result.result_id, frame_ids_to_filter=[frame.frame_id],
                                particles_to_filter=[(other.frame_id, 0)], verbose=False)
    ledger = load_sample_ledger(str(sample / "S" / "ledger.json"))
    chosen = ledger.size_result
    assert chosen.segmentation_id == sam2_result.segmentation_id  # the SAM2 segmentation, not the active one before
    assert chosen.excluded_frame_ids == [frame.frame_id] and chosen.excluded_particles == [(other.frame_id, 0)]
    n_left_out = frame.segmentation(sam2_result.segmentation_id).n_particles + 1
    assert chosen.size_stats.n_par_analysed == sam2_result.size_stats.n_par_analysed - n_left_out

    # Segmented again by another model: the particles left out (numbered per segmentation) are not carried over,
    # the frames left out are
    reanalyze_particle_statistics("S", str(sample), par_segmentation_model="microSAM", verbose=False)
    ledger = load_sample_ledger(str(sample / "S" / "ledger.json"))
    assert ledger.size_result.excluded_particles == [] and ledger.size_result.excluded_frame_ids == [frame.frame_id]


# =============================================================================
# Install helper
# =============================================================================
def test_install_commands(capsys):
    assert install_segmentation.main(["micro-sam", "--dry-run"]) == 0
    first, second = capsys.readouterr().out.strip().splitlines()
    assert "--no-deps micro_sam==1.8.14 git+https://github.com/facebookresearch/segment-anything.git" in first
    assert "torch-em" in second and "napari" in second
    for package in ("PyQt6", "PyQt5", "segment-anything-py", "magicgui", "superqt"):
        assert package not in first + second
    [sam2] = install_segmentation.install_commands("sam2")
    assert sam2[-3:] == ["torch>=2.1", "torchvision>=0.16", "transformers>=4.56,<6"]


# =============================================================================
# Real models (where installed, with their weights)
# =============================================================================
def _sam2_weights() -> bool:
    try:
        from huggingface_hub import try_to_load_from_cache
    except ImportError:
        return False
    return isinstance(try_to_load_from_cache(SAM2.DEFAULT_CHECKPOINT, "config.json"), str)


def _micro_sam_weights() -> bool:
    if not microSAM.is_available():
        return False
    from micro_sam.util import get_cache_directory

    return any((Path(get_cache_directory()) / "models").glob(f"{microSAM.DEFAULT_MODEL_TYPE}*"))


@pytest.mark.skipif(not (SAM2.is_available() and _sam2_weights()), reason="SAM2 or its weights not installed")
def test_real_sam2():
    image = np.full((H, W), 40, np.uint8)
    for cx in (60, 150, 240):
        image[_disk(cx, 100, 25)] = 200
    labels = SAM2.segment_particles(image, PowderMeasurementConfig(
        par_segmentation_model="SAM2", par_segmentation_params={"points_per_crop": 8}))
    assert labels.dtype == np.int32 and labels.max() >= 1


@pytest.mark.skipif(not _micro_sam_weights(),
                    reason="micro-sam or its weights not installed")
def test_real_micro_sam():
    image = np.full((H, W), 40, np.uint8)
    for cx in (60, 150, 240):
        image[_disk(cx, 100, 25)] = 200
    labels = microSAM.segment_particles(image, PowderMeasurementConfig(par_segmentation_model="microSAM"))
    assert labels.dtype == np.int32 and labels.max() >= 1


# =============================================================================
# Particle outlines (shown in the GUI frame browser)
# =============================================================================
def test_frame_outlines():
    from autoemx.core.particle_outlines import frame_outlines, load_outlines, save_outlines

    labels = np.zeros((H, W), np.int32)
    labels[_disk(60, 60, 20)] = 3
    labels[_disk(200, 100, 15)] = 1
    outlines = frame_outlines(labels, counted=[3], others=[1])
    [counted], [other] = outlines["counted"], outlines["others"]
    xy = np.asarray(counted).reshape(-1, 2)
    assert np.allclose(xy.mean(axis=0), (60, 60), atol=2) and np.ptp(xy[:, 0]) == pytest.approx(40, abs=2)
    assert np.asarray(other).reshape(-1, 2)[:, 0].min() == pytest.approx(185, abs=1)


def test_outlines_saved_with_each_segmentation(sample: Path, fake_models):
    from autoemx.core.particle_outlines import load_outlines

    sample_dir = sample / "S"
    reanalyze_particle_statistics("S", str(sample), par_segmentation_model="SAM2", verbose=False)
    ledger = load_sample_ledger(str(sample_dir / "ledger.json"))
    seg_id = ledger.size_result.segmentation_id
    for frame in ledger.frames_of("size_stats", 0):
        outlines = load_outlines(str(sample_dir), seg_id, frame.frame_id)
        # One outline per particle counted, in the order of the frame's particles
        assert len(outlines["counted"]) == frame.segmentation(seg_id).n_particles
        for (x, y, _), poly in zip(frame.segmentation(seg_id).particle_marks_px, outlines["counted"]):
            xy = np.asarray(poly).reshape(-1, 2)
            assert xy[:, 0].min() <= x <= xy[:, 0].max() and xy[:, 1].min() <= y <= xy[:, 1].max()


def test_missing_outlines_added_only_if_the_same_particles_are_found(sample: Path):
    import json
    import shutil

    from autoemx.core.particle_outlines import OUTLINES_DIR, load_outlines

    sample_dir = sample / "S"
    reanalyze_particle_statistics("S", str(sample), par_brightness_thresh=120, verbose=False)
    ledger_path = sample_dir / "ledger.json"
    seg_id = load_sample_ledger(str(ledger_path)).size_result.segmentation_id
    # As if segmented before outlines were saved; one frame recorded with other particles
    shutil.rmtree(sample_dir / OUTLINES_DIR)
    raw = json.loads(ledger_path.read_text())
    changed = next(f for f in raw["frames"] if f["segmentations"][-1]["particle_areas_um2"])
    changed["segmentations"][-1]["particle_areas_um2"][0] *= 1.5
    ledger_path.write_text(json.dumps(raw))

    reanalyze_particle_statistics("S", str(sample), par_brightness_thresh=120, verbose=False)
    ledger = load_sample_ledger(str(ledger_path))
    assert len(ledger.segmentations) == 2  # same segmentation, none added
    for frame in ledger.frames_of("size_stats", 0):
        outlines = load_outlines(str(sample_dir), seg_id, frame.frame_id)
        if frame.frame_id == changed["frame_id"]:
            assert outlines is None  # other particles than recorded: circles kept
        else:
            assert len(outlines["counted"]) == frame.segmentation(seg_id).n_particles
