#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample ledger module of the composition analyzer.

Reads and writes the sample ledger (ledger.json) of the analyzer: ledger path and files, spectrum entries, particles,
quantification runs and clustering analyses. Like the other modules of the analyzer, its methods take the analyzer as
first argument and are called as ``LedgerModule._method(analyzer, ...)``.
"""

import functools
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
from autoemx.config.ledger_io import _background_relpath, ingest_spectra
from autoemx.config.schema_models import (
    ClusteringConfig as LedgerClusteringConfig,
    ClusteringResult,
    Coordinate2D,
    LedgerConfigs,
    LedgerVersionError,
    ParticleInfo,
    QuantificationConfig,
    SampleLedger,
    SpectrumEntry,
    SpotCoordinates,
)
from autoemx.utils.legacy.ledger_bootstrap import create_analyzer_ledger_from_legacy_data_csv
from autoemx.utils.legacy.legacy_config_loader import has_legacy_data_csv
from .reference_matching import ReferenceMatchingModule

from autoemx._logging import get_logger
logger = get_logger(__name__)


class LedgerModule:
    """Sample ledger methods of the composition analyzer."""

    def _get_ledger_path(self: Any) -> str:
        """Return the ledger path for the current sample result directory."""
        return os.path.join(self.sample_result_dir, cnst.LEDGER_FILENAME + cnst.LEDGER_FILEEXT)

    def _build_spectrum_relpath(self: Any, spectrum_id: str) -> str:
        """Build the relative path for one raw spectrum pointer file."""
        # Extension of the format used by the EDS software to export spectra, defined by the microscope driver
        spectrum_ext = dflt.RAW_SPECTRUM_EXT
        em_controller = getattr(self, "EM_controller", None)
        if em_controller is not None:
            try:
                spectrum_ext = getattr(em_controller.EM_driver, "RAW_SPECTRUM_EXT", spectrum_ext)
            except Exception:  # No driver loaded
                pass
        filename = f"{cnst.SPECTRUM_FILENAME_PREFIX}{spectrum_id}{spectrum_ext}"
        return os.path.join(cnst.SPECTRA_DIR, filename)

    @staticmethod
    def _load_background_vector_from_file(background_path: Path) -> Optional[np.ndarray]:
        """Load an instrument background vector from a sidecar file when available."""
        if not background_path.exists() or background_path.suffix.lower() != cnst.VECTOR_FILEEXT:
            return None
        try:
            background_vals = np.load(background_path, allow_pickle=False)
        except Exception:
            return None
        if background_vals is None:
            return None
        return np.asarray(background_vals, dtype=float)

    @staticmethod
    def _parse_optional_float(value: Any) -> Optional[float]:
        if value is None or pd.isna(value) or value == "":
            return None
        try:
            return float(value)
        except Exception:
            return None

    @staticmethod
    def _eq_diameter_um_from_area(area_um: float) -> float:
        """Equivalent circular diameter in µm for a given area in µm²."""
        return float(2.0 * np.sqrt(area_um / np.pi))

    @staticmethod
    def _format_particle_composition_string(
        cluster_ids: List[int],
        clustering_result: ClusteringResult,
    ) -> Optional[str]:
        """Build the Candidate phases / mixtures string for a particle's clusters."""
        phase_parts: List[str] = []
        mixture_parts: List[str] = []
        refs_rows = clustering_result.refs_assigned_rows or []
        mixtures_by_cluster = clustering_result.clusters_assigned_mixtures or []

        for cluster_id in cluster_ids:
            if cluster_id < 0:
                continue
            if cluster_id < len(refs_rows):
                row = refs_rows[cluster_id] or {}
                phase_name = row.get(f"{cnst.CND_DF_KEY}1")
                phase_conf = row.get(f"{cnst.CS_CND_DF_KEY}1")
                if phase_name:
                    if phase_conf is None:
                        phase_parts.append(str(phase_name))
                    else:
                        phase_parts.append(f"{phase_name} (CScnd: {float(phase_conf):.2f})")

            if cluster_id < len(mixtures_by_cluster):
                mixtures = mixtures_by_cluster[cluster_id] or []
                if mixtures:
                    top_mix = ReferenceMatchingModule._sort_mixtures(mixtures)[0]
                    refs = top_mix.get(cnst.REF_NAME_KEY) or []
                    conf = top_mix.get(cnst.CONF_SCORE_KEY)
                    if refs:
                        mix_name = "+".join(str(r) for r in refs)
                        if conf is None:
                            mixture_parts.append(mix_name)
                        else:
                            mixture_parts.append(f"{mix_name} (CSmix: {float(conf):.2f})")

        sections: List[str] = []
        if phase_parts:
            sections.append("Candidate phases: " + ", ".join(phase_parts))
        if mixture_parts:
            sections.append("candidate mixtures: " + ", ".join(mixture_parts))
        if not sections:
            return None
        return "; ".join(sections)

    @staticmethod
    def _update_non_fingerprint_options(quant_config: QuantificationConfig, candidate: QuantificationConfig) -> None:
        """Store the current values of the options that do not affect the results (not in the fingerprint) in a
        reused quantification run, e.g. is_known_phase_mixture, which the analysis reads."""
        for key in QuantificationConfig.NON_FINGERPRINT_OPTIONS:
            if key in candidate.options:
                quant_config.options[key] = candidate.options[key]

    @staticmethod
    def _quantification_configs_match(left: QuantificationConfig, right: QuantificationConfig) -> bool:
        """Compare quantification configs by deterministic scientific-input fingerprint."""
        return left.fingerprint() == right.fingerprint()

    @staticmethod
    def _clustering_configs_match(left: LedgerClusteringConfig, right: LedgerClusteringConfig) -> bool:
        """Compare clustering configs by deterministic scientific-input fingerprint."""
        return left.fingerprint() == right.fingerprint()

    @staticmethod
    def _format_quantification_config_changes(
        changes: Dict[str, Dict[str, Any]],
        max_items: int = 8,
        max_value_length: int = 80,
    ) -> str:
        """Return a concise deterministic summary of config changes, one field per line."""
        if not changes:
            return "  • none"

        sorted_items = sorted(changes.items(), key=lambda item: item[0])
        formatted_entries: List[str] = []
        for path, change in sorted_items[:max_items]:
            old_val = LedgerModule._short_repr(change.get("old"), max_value_length)
            new_val = LedgerModule._short_repr(change.get("new"), max_value_length)
            formatted_entries.append(f"  • {path}: {old_val} -> {new_val}")

        remaining = len(sorted_items) - len(formatted_entries)
        if remaining > 0:
            formatted_entries.append(f"  • ... and {remaining} more")

        return "\n".join(formatted_entries)

    @staticmethod
    def _short_repr(value: Any, max_length: int) -> str:
        """Return a stable shortened representation suitable for warning messages."""
        text = repr(value)
        if len(text) <= max_length:
            return text
        return text[: max_length - 3] + "..."

    @staticmethod
    def _runtime_clustering_cfg_payload(clustering_cfg: LedgerClusteringConfig) -> Dict[str, Any]:
        """Convert a ledger clustering config to the runtime config payload shape."""
        payload = clustering_cfg.model_dump(mode="json")
        payload.pop("clustering_id", None)
        return payload

    @staticmethod
    def _next_clustering_config_id(quant_config: QuantificationConfig) -> int:
        """Return the next clustering config id local to the provided quantification config."""
        if not quant_config.clustering_analyses:
            return 0
        return max(analysis.config.clustering_id for analysis in quant_config.clustering_analyses) + 1

    @staticmethod
    def _next_quantification_id(existing_ledger: Optional[SampleLedger]) -> int:
        """Return the next ledger-local integer quantification identifier."""
        if existing_ledger is None or not existing_ledger.quantifications:
            return 0
        return max(config.quantification_id for config in existing_ledger.quantifications) + 1

    @staticmethod
    def _ledger_has_quantification_results(
        existing_ledger: Optional[SampleLedger],
        quantification_id: int,
    ) -> bool:
        """Return True when at least one spectrum contains results for the given quantification id."""
        if existing_ledger is None:
            return False
        for spectrum in existing_ledger.spectra:
            for result in spectrum.quantification_results:
                if result.quantification_id == quantification_id:
                    return True
        return False

    def _write_manufacturer_background_vector(self: Any, spectrum_id: str, background_vals: Optional[List[float]]) -> Optional[str]:
        """Persist manufacturer background counts as a companion vector file when enabled."""
        if background_vals is None or not self.quant_cfg.use_instrument_background:
            return None

        background_relpath = _background_relpath(spectrum_id)
        background_abs_path = os.path.join(self.sample_result_dir, background_relpath)
        os.makedirs(os.path.dirname(background_abs_path), exist_ok=True)
        np.save(background_abs_path, np.asarray(list(map(float, background_vals)), dtype=float))
        return background_relpath

    def _load_existing_ledger(self: Any) -> Optional[SampleLedger]:
        """Load an existing ledger if present and valid."""
        ledger_path = LedgerModule._get_ledger_path(self)
        if not os.path.exists(ledger_path):
            return None
        try:
            return SampleLedger.from_json_file(ledger_path)
        except LedgerVersionError:
            # Never replace a ledger written by a newer AutoEMX version with a new one
            raise
        except Exception:
            return None

    def _get_ingested_spectrum_ids(self: Any) -> set[str]:
        """Return the spectrum ids already ingested for the current sample directory."""
        cached_sample_dir = getattr(self, "_ingested_spectrum_ids_sample_result_dir", None)
        ingested_spectrum_ids = getattr(self, "_ingested_spectrum_ids", None)

        if ingested_spectrum_ids is None or cached_sample_dir != self.sample_result_dir:
            ingested_spectrum_ids = set()
            ledger = LedgerModule._load_existing_ledger(self)
            if ledger is not None:
                ingested_spectrum_ids.update(
                    str(entry.spectrum_id)
                    for entry in ledger.spectra
                    if entry.spectrum_id not in (None, "")
                )
            self._ingested_spectrum_ids = ingested_spectrum_ids
            self._ingested_spectrum_ids_sample_result_dir = self.sample_result_dir

        return ingested_spectrum_ids

    def _get_spectra_dir(self: Any) -> str:
        """Return the absolute path to the spectra pointer directory."""
        return os.path.join(self.sample_result_dir, cnst.SPECTRA_DIR)

    def _build_spot_coordinates(
        self: Any,
        machine_x: Any = None,
        machine_y: Any = None,
        pixel_x: Any = None,
        pixel_y: Any = None,
    ) -> Optional[SpotCoordinates]:
        machine_coords = None
        pixel_coords = None

        x_pixel = LedgerModule._parse_optional_float(pixel_x)
        y_pixel = LedgerModule._parse_optional_float(pixel_y)
        if x_pixel is not None and y_pixel is not None:
            pixel_coords = (int(round(x_pixel)), int(round(y_pixel)))

        if pixel_coords is not None and hasattr(self, "EM_controller"):
            try:
                em_driver = self.EM_controller.EM_driver
                rel_coords = em_driver.frame_pixel_to_rel_coords(
                    np.asarray([pixel_coords], dtype=float),
                    int(self.EM_controller.im_width),
                    int(self.EM_controller.im_height),
                )
                machine_coords = Coordinate2D(
                    x=float(rel_coords[0][0]),
                    y=float(rel_coords[0][1]),
                )
            except Exception:
                machine_coords = None

        if machine_coords is None:
            x_machine = LedgerModule._parse_optional_float(machine_x)
            y_machine = LedgerModule._parse_optional_float(machine_y)
            if x_machine is not None and y_machine is not None:
                machine_coords = Coordinate2D(x=x_machine, y=y_machine)

        if machine_coords is None and pixel_coords is None:
            return None

        return SpotCoordinates(
            machine_coordinates=machine_coords,
            pixel_coordinates=pixel_coords,
        )

    def _load_or_create_ledger(self: Any) -> SampleLedger:
        """Load a ledger, using legacy Data.csv bootstrap only for non-acquisition workflows."""
        if self.is_acquisition:
            self._last_acq_ledger_created = False
            self._last_acq_ledger_ingested_spectra_count = 0

            spectra_dir = LedgerModule._get_spectra_dir(self)
            os.makedirs(spectra_dir, exist_ok=True)

            ledger = LedgerModule._load_existing_ledger(self)
            ledger_changed = False
            if ledger is None:
                ledger = SampleLedger(
                    sample_id=self.sample_id,
                    sample_path=os.path.abspath(self.sample_result_dir),
                    configs=LedgerModule._build_ledger_configs(self),
                    spectra=[],
                    quantifications=[],
                    active_quant=None,
                )
                self._last_acq_ledger_created = True
                ledger_changed = True

            n_ingested, n_removed = ingest_spectra(ledger)
            self._last_acq_ledger_ingested_spectra_count = n_ingested
            if n_ingested > 0 or n_removed > 0:
                ledger_changed = True

            if ledger.sample_path != os.path.abspath(self.sample_result_dir):
                ledger.sample_path = os.path.abspath(self.sample_result_dir)
                ledger_changed = True

            if ledger_changed:
                ledger.to_json_file(LedgerModule._get_ledger_path(self))

            return ledger

        ledger = LedgerModule._load_existing_ledger(self)
        if ledger is not None:
            return ledger

        if has_legacy_data_csv(self.sample_result_dir):
            return create_analyzer_ledger_from_legacy_data_csv(
                self,
                default_ledger_configs=LedgerModule._build_ledger_configs(self),
                write_background_pointer=functools.partial(LedgerModule._write_manufacturer_background_vector, self),
            )

        raise FileNotFoundError(
            f"No ledger.json or Data.csv found for sample '{self.sample_id}' "
            f"in '{self.sample_result_dir}'. Please ensure the sample directory contains "
            "either ledger.json or Data.csv."
        )

    def _build_ledger_configs(self: Any) -> LedgerConfigs:
        """Build inline ledger configs from current analyzer configuration objects."""
        self.measurement_cfg = self.measurement_cfg.model_copy(
            update={
                "powder_meas_cfg": self.powder_meas_cfg,
                "bulk_meas_cfg": self.bulk_meas_cfg,
                "exp_stds_cfg": self.exp_stds_cfg if self.exp_stds_cfg.is_exp_std_measurement else None,
            }
        )
        return LedgerConfigs(
            microscope_cfg=self.microscope_cfg,
            sample_cfg=self.sample_cfg,
            measurement_cfg=self.measurement_cfg,
            sample_substrate_cfg=self.sample_substrate_cfg,
            plot_cfg=self.plot_cfg,
        )

    def _extract_spectrum_info(
        self: Any,
        spectrum: 'SpectrumEntry',
        index: int,
        particles_by_id: Optional[Dict[int, ParticleInfo]] = None,
    ) -> dict:
        """Extract spectrum metadata from a ledger SpectrumEntry for use in CSV export."""
        par_id = frame_id = ""
        particle_id_int: Optional[int] = None
        acq = spectrum.acquisition_details
        if acq is not None:
            particle_id_int = acq.particle_id
            par_id = str(acq.particle_id) if acq.particle_id is not None else ""
            frame_id = str(acq.frame_id or "")

        row = {
            cnst.SP_ID_DF_KEY: str(spectrum.spectrum_id) if spectrum.spectrum_id is not None else str(index),
            cnst.PAR_ID_DF_KEY: par_id,
        }
        if particles_by_id:
            particle = particles_by_id.get(particle_id_int) if particle_id_int is not None else None
            if particle is not None and particle.area_um is not None:
                row[cnst.PAR_AREA_UM_KEY] = round(float(particle.area_um), 3)
            else:
                row[cnst.PAR_AREA_UM_KEY] = ""
            if particle is not None and particle.eq_diameter_um is not None:
                row[cnst.PAR_EQ_D_KEY] = round(float(particle.eq_diameter_um), 3)
            else:
                row[cnst.PAR_EQ_D_KEY] = ""
        row[cnst.FRAME_ID_DF_KEY] = frame_id
        return row

    def _upsert_particle_info(
        self: Any,
        ledger: SampleLedger,
        particle_id: int,
        area_um: Optional[float],
        frame_id: Optional[str] = None,
        coordinates: Optional[Coordinate2D] = None,
        image_pixel_size_um: Optional[float] = None,
    ) -> None:
        """Create or update the ParticleInfo entry for a particle in the ledger."""
        eq_diameter_um = (
            LedgerModule._eq_diameter_um_from_area(area_um) if area_um is not None else None
        )
        for particle in ledger.particles:
            if particle.id == particle_id:
                if area_um is not None:
                    particle.area_um = area_um
                    particle.eq_diameter_um = eq_diameter_um
                if frame_id is not None:
                    particle.frame_id = frame_id
                if coordinates is not None:
                    particle.coordinates = coordinates
                if image_pixel_size_um is not None:
                    particle.image_pixel_size_um = image_pixel_size_um
                return
        ledger.particles.append(
            ParticleInfo(
                id=particle_id,
                area_um=area_um,
                eq_diameter_um=eq_diameter_um,
                frame_id=frame_id,
                coordinates=coordinates,
                image_pixel_size_um=image_pixel_size_um,
            )
        )

    def _update_particles_from_clustering(
        self: Any,
        ledger: SampleLedger,
        labels: Optional[List],
        df_indices: Optional[List],
        clustering_result: ClusteringResult,
    ) -> bool:
        """
        Fill ParticleInfo.clusters and .composition from the active clustering run.

        Returns True if the ledger was modified.
        """
        if not ledger.particles or labels is None or df_indices is None:
            return False

        clusters_by_particle: Dict[int, List[int]] = {
            particle.id: [] for particle in ledger.particles
        }
        for label_pos, spectrum_index in enumerate(df_indices):
            if label_pos >= len(labels):
                continue
            cluster_id = labels[label_pos]
            if cluster_id is None or pd.isna(cluster_id):
                continue
            cluster_id_int = int(cluster_id)
            if cluster_id_int < 0:
                continue
            if spectrum_index < 0 or spectrum_index >= len(ledger.spectra):
                continue
            acq = ledger.spectra[spectrum_index].acquisition_details
            if acq is None or acq.particle_id is None:
                continue
            particle_id = int(acq.particle_id)
            if particle_id not in clusters_by_particle:
                continue
            if cluster_id_int not in clusters_by_particle[particle_id]:
                clusters_by_particle[particle_id].append(cluster_id_int)

        changed = False
        for particle in ledger.particles:
            ordered_clusters = sorted(clusters_by_particle.get(particle.id, []))
            composition = LedgerModule._format_particle_composition_string(
                ordered_clusters, clustering_result
            )
            if particle.clusters != ordered_clusters or particle.composition != composition:
                particle.clusters = ordered_clusters
                particle.composition = composition
                changed = True
        return changed

    def _attach_particle_stats_ledger(self: Any) -> None:
        """Create or load a SampleLedger and attach it to the particle finder for PSD."""
        ledger = LedgerModule._load_existing_ledger(self)
        if ledger is None:
            ledger = SampleLedger(
                sample_id=self.sample_id,
                sample_path=os.path.abspath(self.sample_result_dir),
                configs=LedgerModule._build_ledger_configs(self),
                spectra=[],
                particles=[],
                quantifications=[],
                active_quant=None,
            )
        else:
            # Fresh PSD run replaces prior particle geometry records.
            ledger.particles = []
            ledger.configs = LedgerModule._build_ledger_configs(self)
        ledger_path = LedgerModule._get_ledger_path(self)
        ledger.to_json_file(ledger_path)

        particle_finder = None
        if getattr(self, "EM_controller", None) is not None:
            particle_finder = getattr(self.EM_controller, "particle_finder", None)
        if particle_finder is not None:
            particle_finder.attach_stats_ledger(ledger, ledger_path=ledger_path)

    def _configure_list_particle_selection(self: Any) -> None:
        """Validate and attach list-mode particle targets to the particle finder."""
        if self._particle_ids is not None and self._particle_coords_mm is not None:
            raise ValueError(
                "Provide exactly one of particle_ids or particle_coords_mm for "
                "par_selection_mode='list'"
            )
        if self._particle_ids is None and self._particle_coords_mm is None:
            raise ValueError(
                "par_selection_mode='list' requires particle_ids or particle_coords_mm"
            )

        particle_finder = None
        if getattr(self, "EM_controller", None) is not None:
            particle_finder = getattr(self.EM_controller, "particle_finder", None)
        if particle_finder is None:
            raise RuntimeError(
                "Particle finder is not available for list particle selection"
            )

        ledger = LedgerModule._load_or_create_ledger(self)
        particle_finder.set_particle_targets(
            particle_ids=self._particle_ids,
            particle_coords_mm=self._particle_coords_mm,
            ledger=ledger,
        )

    def _record_acquired_spectrum(
        self: Any,
        spectrum_entry: SpectrumEntry,
        particle_id: Optional[int] = None,
        particle_area_um: Optional[float] = None,
        frame_id: Optional[str] = None,
        particle_coordinates: Optional[Coordinate2D] = None,
        image_pixel_size_um: Optional[float] = None,
    ) -> None:
        """
        Add a newly acquired spectrum to the ledger (created if missing), with its particle if particle_id is
        given, and save the ledger.
        """
        try:
            ledger = LedgerModule._load_existing_ledger(self)
        except Exception:
            ledger = None
        if ledger is None:
            ledger = SampleLedger(
                sample_id=self.sample_id,
                sample_path=os.path.abspath(self.sample_result_dir),
                configs=LedgerModule._build_ledger_configs(self),
                spectra=[],
                particles=[],
                quantifications=[],
                active_quant=None,
            )
        ledger.spectra.append(spectrum_entry)
        if particle_id is not None:
            LedgerModule._upsert_particle_info(
                self,
                ledger,
                particle_id=particle_id,
                area_um=particle_area_um,
                frame_id=frame_id,
                coordinates=particle_coordinates,
                image_pixel_size_um=image_pixel_size_um,
            )
        ledger.to_json_file(LedgerModule._get_ledger_path(self))

    def _backfill_particles_from_clustering(
        self: Any,
        labels: Optional[List],
        df_indices: Optional[List],
        clustering_result: ClusteringResult,
    ) -> None:
        """Save the clusters and composition of each particle of the ledger, from the clustering result."""
        ledger = LedgerModule._load_or_create_ledger(self)
        if LedgerModule._update_particles_from_clustering(self, ledger, labels, df_indices, clustering_result):
            ledger.to_json_file(LedgerModule._get_ledger_path(self))
