#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sample folders with particle size data only (no spectra), built without acquisition, for the GUI tests."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence, Tuple

import numpy as np

import autoemx.utils.constants as cnst

from autoemx.config.ledger_io import create_sample_ledger
from autoemx.config.ledger_schemas import (
    FrameRecord,
    FrameSegmentation,
    LedgerConfigs,
    SampleLedger,
    SegmentationConfig,
    SizeBand,
    SizeMeasurement,
)
from autoemx.config.runtime_configs import (
    MeasurementConfig,
    MicroscopeConfig,
    PlotConfig,
    PowderMeasurementConfig,
    SampleConfig,
    SampleSubstrateConfig,
)
from autoemx.core.particle_statistics import recalculate_particle_stats


def area_um2(d_um: float) -> float:
    return float(np.pi * d_um**2 / 4)


FRAME_PX = (640, 480)


def _frame_image(sample_dir: Path, frame_id: str, width_um: float, diameters_um, rng) -> Tuple[str, list]:
    """Synthetic frame image (bright discs on a dark background) in the 'SEM images' folder, and its particle marks
    (x, y, radius in pixels)."""
    from PIL import Image, ImageDraw

    images_dir = sample_dir / cnst.IMAGES_DIR
    images_dir.mkdir(exist_ok=True)
    px_um = width_um / FRAME_PX[0]
    im = Image.fromarray(rng.normal(40, 8, FRAME_PX[::-1]).clip(0, 255).astype(np.uint8))
    draw = ImageDraw.Draw(im)
    marks = []
    for d in diameters_um:
        r = max(2, int(round(d / px_um / 2)))
        x = int(rng.integers(r + 1, max(r + 2, FRAME_PX[0] - r - 1)))
        y = int(rng.integers(r + 1, max(r + 2, FRAME_PX[1] - r - 1)))
        draw.ellipse([x - r, y - r, x + r, y + r], fill=int(rng.integers(150, 230)))
        marks.append((x, y, r))
    name = f"{frame_id}_particles.png"
    im.save(images_dir / name)
    return name, marks


def _ledger(sample_dir: Path, elements: Sequence[str] = ()) -> SampleLedger:
    sample_dir.mkdir(parents=True)
    configs = LedgerConfigs(
        microscope_cfg=MicroscopeConfig(ID="PhenomXL", type="SEM"),
        sample_cfg=SampleConfig(elements=list(elements), type="powder"),
        sample_substrate_cfg=SampleSubstrateConfig(), measurement_cfg=MeasurementConfig(), plot_cfg=PlotConfig(),
    )
    return create_sample_ledger(str(sample_dir), sample_dir.name, configs)


def _segmentation(ledger: SampleLedger, cfg: PowderMeasurementConfig) -> int:
    return ledger.append_segmentation(SegmentationConfig.from_powder_meas_cfg(cfg, 0))


def _frame(ledger, frame_id, purpose, width_um, window, diameters_um, segmentation_id, rng, images, sample_dir,
           measurement_id=None, center=None):
    """A frame with its particles of the given diameters (and a synthetic image if ``images``)."""
    image_file, marks = _frame_image(sample_dir, frame_id, width_um, diameters_um, rng) if images else (None, [])
    return ledger.add_frame(FrameRecord(
        frame_id=frame_id, label=frame_id.removeprefix("EDS-"), purpose=purpose, measurement_id=measurement_id,
        center_mm=center, frame_width_um=width_um, pixel_size_um=width_um / FRAME_PX[0],
        area_um2=width_um * width_um * FRAME_PX[1] / FRAME_PX[0], search_window_um=window, image_file=image_file,
        segmentations=[FrameSegmentation(segmentation_id=segmentation_id, d_range_um=window,
                                         particle_areas_um2=[area_um2(d) for d in diameters_um],
                                         particle_edge_corrections=[1.0] * len(diameters_um), particle_marks_px=marks)],
    ))


def make_measurement_sample(sample_dir: Path, seed: int = 0, median_um: float = 4.0, sigma: float = 0.5,
                            images: bool = False) -> Path:
    """Sample measured only by a particle size statistics measurement: two size bands of 6 frames each (with
    synthetic frame images and particle marks if ``images``), and its statistics."""
    rng = np.random.default_rng(seed)
    ledger = _ledger(sample_dir)
    cfg = PowderMeasurementConfig()
    measurement = SizeMeasurement(measurement_id=0, status="completed", n_par_target=40, powder_meas_cfg=cfg,
                                  segmentation_id=_segmentation(ledger, cfg))
    ledger.size_measurements.append(measurement)
    for b, (d_min, d_max, width) in enumerate(((median_um / 2, median_um * 10, 120.0),
                                               (median_um / 10, median_um / 2, 30.0))):
        band = SizeBand(d_min_um=d_min, d_max_um=d_max, frame_width_um=width)
        measurement.bands.append(band)
        for f in range(6):
            ds = rng.lognormal(np.log(median_um), sigma, 40)
            ds = ds[(ds >= d_min) & (ds < d_max)][:15]
            frame = _frame(ledger, f"B{b}F{f}", "size_stats", width, (d_min, d_max), ds, measurement.segmentation_id,
                           rng, images, sample_dir, measurement_id=0)
            band.frame_ids.append(frame.frame_id)
    ledger.append_size_result(recalculate_particle_stats(ledger))
    ledger.to_json_file(str(sample_dir / "ledger.json"))
    return sample_dir


def make_eds_frames_sample(sample_dir: Path, seed: int = 1, median_um: float = 8.0, images: bool = False) -> Path:
    """Sample with only frames searched by EDS (no spectra, no statistics saved in the ledger; with synthetic frame
    images and particle marks if ``images``)."""
    rng = np.random.default_rng(seed)
    ledger = _ledger(sample_dir, ["Si"])
    segmentation_id = _segmentation(ledger, ledger.configs.measurement_cfg.powder_meas_cfg or PowderMeasurementConfig())
    for f in range(6):
        ds = rng.lognormal(np.log(median_um), 0.3, 10)
        _frame(ledger, f"EDS-A{f}", "eds_search", 50.0, (3.0, 30.0), ds, segmentation_id, rng, images, sample_dir)
    ledger.to_json_file(str(sample_dir / "ledger.json"))
    return sample_dir


def make_empty_sample(sample_dir: Path) -> Path:
    """Sample with neither spectra nor particle size data."""
    _ledger(sample_dir, ["Si"]).to_json_file(str(sample_dir / "ledger.json"))
    return sample_dir
