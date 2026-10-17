#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Legacy compatibility: particle size statistics of ledgers saved before frames and segmentations were recorded
(``particle_stats``, ``particle_stats_progress`` and ``eds_size_scans``), converted to ``frames``,
``segmentations``, ``size_measurements`` and ``size_results`` when the ledger is loaded.

- The settings of the measurement (``particle_stats.powder_meas_cfg``) and of the EDS acquisition give the
  segmentations (one if they are the same).
- The frames searched by EDS, and the frames of the size bands of the measurement (and of the measurement in
  progress), become frames with the result of their segmentation; the particles of the size bands, recorded as
  ``ParticleInfo`` (selection 'size_stats'), become the particles of their frame (those measured by EDS stay).
- ``particle_stats`` becomes measurement 0 (completed) and result 0 (the active one, statistics kept), the particles
  left out becoming (frame_id, index) pairs.

Frames of earlier measurements, replaced in the ledger by a later one, were not kept and cannot be recovered.
"""

from typing import Any, Dict, List, Optional, Set, Tuple

from autoemx.config.runtime_configs import PowderMeasurementConfig

from .segmentation import SegmentationConfig

LEGACY_KEYS = ("particle_stats", "particle_stats_progress", "eds_size_scans")


def _defaults_if_missing(record: Dict[str, Any]) -> Dict[str, Any]:
    """Without version or date recorded, the defaults of the model."""
    return {k: v for k, v in record.items() if not (k in ("autoemx_version", "datetime") and v is None)}


def migrate_legacy_size_statistics(data: Any) -> Any:
    """Convert the legacy particle size statistics of a saved ledger (plain JSON data); other data is returned as is."""
    if not isinstance(data, dict) or not any(key in data for key in LEGACY_KEYS):
        return data
    data = dict(data)
    stats, progress = data.pop("particle_stats", None), data.pop("particle_stats_progress", None)
    eds_scans = data.pop("eds_size_scans", None) or []
    if not (stats or progress or eds_scans):
        return data

    segmentations: List[Dict[str, Any]] = list(data.get("segmentations") or [])
    frames: List[Dict[str, Any]] = list(data.get("frames") or [])
    measurements: List[Dict[str, Any]] = list(data.get("size_measurements") or [])
    results: List[Dict[str, Any]] = list(data.get("size_results") or [])
    used_frame_ids: Set[str] = {f["frame_id"] for f in frames}

    def segmentation_id(powder_meas_cfg: Optional[dict]) -> int:
        cfg = PowderMeasurementConfig.model_validate(powder_meas_cfg or {})
        candidate = SegmentationConfig.from_powder_meas_cfg(cfg, len(segmentations))
        for existing in segmentations:
            if SegmentationConfig.model_validate(existing).fingerprint() == candidate.fingerprint():
                return existing["segmentation_id"]
        segmentations.append(candidate.model_dump())
        return candidate.segmentation_id

    def unique_frame_id(frame_id: str) -> str:
        new_id, n = frame_id, 2
        while new_id in used_frame_ids:
            new_id, n = f"{frame_id}#{n}", n + 1
        used_frame_ids.add(new_id)
        return new_id

    def frame_record(frame: dict, frame_id: str, purpose: str, measurement_id: Optional[int], window: Tuple[float, float],
                     result: dict) -> dict:
        return {
            "frame_id": frame_id, "label": frame["frame_id"].removeprefix("EDS-").split("#")[0], "purpose": purpose,
            "measurement_id": measurement_id, "datetime": None, "center_mm": frame.get("center_mm"),
            "frame_width_um": frame["frame_width_um"], "pixel_size_um": frame["pixel_size_um"],
            "area_um2": frame["area_um2"], "search_window_um": list(window), "image_file": frame.get("image_file"),
            "segmentations": [result],
        }

    # Frames searched by EDS
    measurement_cfg = (data.get("configs") or {}).get("measurement_cfg") or {}
    eds_segmentation = segmentation_id(measurement_cfg.get("powder_meas_cfg")) if eds_scans else None
    for scan in eds_scans:
        window = (scan["d_min_um"], scan["d_max_um"])
        for frame in scan.get("frames") or []:
            frames.append(frame_record(frame, unique_frame_id(frame["frame_id"]), "eds_search", None, window, {
                "segmentation_id": eds_segmentation, "d_range_um": list(window),
                "particle_areas_um2": frame.get("particle_areas_um2") or [],
                "particle_edge_corrections": frame.get("particle_edge_corrections") or [],
                "particle_marks_px": [list(m[:3]) for m in frame.get("particle_marks_px") or []],
            }))

    # Frames of the size bands, with their particles (ParticleInfo 'size_stats')
    particles: List[Dict[str, Any]] = list(data.get("particles") or [])
    particles_by_id = {p["id"]: p for p in particles}
    moved_particles: Set[int] = set()

    def migrate_bands(record: dict, measurement_id: int, segmentation: int):
        bands, key_of, frame_map = [], {}, {}
        for band in record.get("bands") or []:
            window = (band["d_min_um"], band["d_max_um"])
            band_ids = [pid for pid in band.get("particle_ids") or [] if pid in particles_by_id]
            frame_ids = []
            for frame in band.get("frames") or []:
                frame_id = unique_frame_id(frame["frame_id"])
                frame_map[frame["frame_id"]] = frame_id
                frame_particles = [particles_by_id[pid] for pid in band_ids
                                   if particles_by_id[pid].get("frame_id") == frame["frame_id"]]
                marks = {m[3]: list(m[:3]) for m in frame.get("particle_marks_px") or [] if m[3] is not None}
                result = {"segmentation_id": segmentation, "d_range_um": list(window),
                          "particle_areas_um2": [p["area_um"] for p in frame_particles],
                          "particle_edge_corrections": [p.get("edge_correction") or 1.0 for p in frame_particles],
                          "particle_marks_px": [marks[p["id"]] for p in frame_particles if p["id"] in marks]}
                if len(result["particle_marks_px"]) != len(frame_particles):
                    result["particle_marks_px"] = []
                for index, particle in enumerate(frame_particles):
                    key_of[particle["id"]] = (frame_id, index)
                    moved_particles.add(particle["id"])
                frames.append(frame_record(frame, frame_id, "size_stats", measurement_id, window, result))
                frame_ids.append(frame_id)
            bands.append({"d_min_um": band["d_min_um"], "d_max_um": band["d_max_um"],
                          "frame_width_um": band.get("frame_width_um"), "frame_ids": frame_ids})
        return bands, key_of, frame_map

    def measurement(record: dict, measurement_id: int, status: str):
        segmentation = segmentation_id(record.get("powder_meas_cfg"))
        bands, key_of, frame_map = migrate_bands(record, measurement_id, segmentation)
        measurements.append(_defaults_if_missing({
            "measurement_id": measurement_id, "datetime": record.get("datetime"), "status": status,
            "autoemx_version": record.get("autoemx_version"), "n_par_target": record.get("n_par_target") or 0,
            "powder_meas_cfg": record.get("powder_meas_cfg") or {}, "segmentation_id": segmentation,
            "bands": bands, "precision_reached": record.get("precision_reached"),
        }))
        return segmentation, key_of, frame_map

    next_measurement = max([m["measurement_id"] for m in measurements], default=-1) + 1
    stats_measurement, segmentation, key_of, frame_map = None, None, {}, {}
    if stats and stats.get("bands"):
        stats_measurement = next_measurement
        segmentation, key_of, frame_map = measurement(stats, stats_measurement, "completed")
        next_measurement += 1
    if progress:
        measurement(progress, next_measurement, "in_progress")

    if stats:
        def summary(band: dict) -> dict:
            band_frames = band.get("frames") or []
            return {
                "d_min_um": band["d_min_um"], "d_max_um": band["d_max_um"], "frame_width_um": band.get("frame_width_um"),
                "n_frames": band.get("n_frames", len(band_frames)), "scanned_area_um2": band.get("scanned_area_um2", 0.0),
                "n_particles": len(band.get("particle_ids") or []) or sum(f.get("n_particles", 0) for f in band_frames),
                "number_density_per_mm2": band.get("number_density_per_mm2"),
                "n_eds_particles": band.get("n_eds_particles", 0),
                "frame_ids": [frame_map.get(f["frame_id"], f["frame_id"]) for f in band_frames],
            }

        result_id = max([r["result_id"] for r in results], default=-1) + 1
        results.append(_defaults_if_missing({
            "result_id": result_id, "autoemx_version": stats.get("autoemx_version"), "datetime": stats.get("datetime"),
            "measurement_id": stats_measurement, "segmentation_id": segmentation,
            "eds_segmentation_id": eds_segmentation if stats.get("eds_scans") else None,
            "n_par_target": stats.get("n_par_target"), "powder_meas_cfg": stats.get("powder_meas_cfg") or {},
            "n_par_by_selection": stats.get("n_par_by_selection") or {}, "size_stats": stats["size_stats"],
            "unweighted_size_stats": stats.get("unweighted_size_stats"),
            "n_par_without_weight": stats.get("n_par_without_weight", 0),
            "confidence_intervals_95": stats.get("confidence_intervals_95") or {},
            "precision_reached": stats.get("precision_reached"),
            "bands": [summary(b) for b in stats.get("bands") or []],
            "eds_scans": [summary(s) for s in stats.get("eds_scans") or []],
            "excluded_frame_ids": [frame_map.get(f, f) for f in stats.get("excluded_frame_ids") or []],
            "excluded_particles": [list(key_of[pid]) for pid in stats.get("excluded_particle_ids") or []
                                   if pid in key_of],
            "overlapping_frame_ids": [frame_map.get(f, f) for f in stats.get("overlapping_frame_ids") or []],
        }))
        if data.get("active_size_result") is None:
            data["active_size_result"] = result_id

    # Particles of the size bands now in their frames (those measured by EDS, with spectra, stay)
    measured_by_eds = {(s.get("acquisition_details") or {}).get("particle_id") for s in data.get("spectra") or []}
    data["particles"] = [p for p in particles if p["id"] not in moved_particles or p["id"] in measured_by_eds]
    data.update(frames=frames, segmentations=segmentations, size_measurements=measurements, size_results=results)
    return data
