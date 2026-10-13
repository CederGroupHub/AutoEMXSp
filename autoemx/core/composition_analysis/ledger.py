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
import re
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import autoemx.calibrations as calibs
import autoemx.config.defaults as dflt
import autoemx.utils.constants as cnst
from autoemx._version import __version__ as AUTOEMX_VERSION
from autoemx.config.ledger_io import (
    _background_relpath,
    create_sample_ledger,
    get_active_quantification_config,
    ingest_spectra,
)
from autoemx.config.schema_models import (
    ClusteringAnalysis,
    ClusteringConfig as LedgerClusteringConfig,
    ClusteringResult,
    Coordinate2D,
    LedgerConfigs,
    LedgerVersionError,
    EDSStandardsFile,
    ParticleInfo,
    QuantificationConfig,
    QuantificationResult,
    SampleLedger,
    SpectrumEntry,
    SpotCoordinates,
)
from autoemx.utils.legacy.ledger_bootstrap import create_analyzer_ledger_from_legacy_data_csv
from autoemx.utils.legacy.legacy_config_loader import has_legacy_data_csv
from autoemx.core.quantifier import XSp_Quantifier
from .reference_matching import ReferenceMatchingModule
from .standards import StandardsModule

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
                ledger = create_sample_ledger(self.sample_result_dir, self.sample_id, LedgerModule._build_ledger_configs(self))
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
            ledger = create_sample_ledger(self.sample_result_dir, self.sample_id, LedgerModule._build_ledger_configs(self))
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
            ledger = create_sample_ledger(self.sample_result_dir, self.sample_id, LedgerModule._build_ledger_configs(self))
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

    def _ensure_quant_tracking_length(self: Any, total_spectra: int) -> None:
        """Ensure in-memory quantification tracking lists are indexable for all spectra."""
        if len(self.spectra_quant_records) < total_spectra:
            self.spectra_quant_records.extend([None] * (total_spectra - len(self.spectra_quant_records)))

    def _ensure_current_quantification_run(self: Any, force_new: bool = False) -> None:
        """Resolve the active quantification run for this launch.

        Parameters
        ----------
        force_new : bool
            If True, always create a new quantification config id.
        """
        existing_ledger = LedgerModule._load_or_create_ledger(self)
        active_quant_config = LedgerModule._get_active_quantification_config(existing_ledger)
        active_quant_has_results = (
            active_quant_config is not None
            and LedgerModule._ledger_has_quantification_results(
                existing_ledger,
                int(active_quant_config.quantification_id),
            )
        )

        current_sample_elements = list(self.all_els_sample)
        current_substrate_elements = list(self.all_els_substrate)
        current_options = LedgerModule._build_quantification_options(self)
        current_reference_values = LedgerModule._get_reference_values_by_el_line(self, active_quant_config=active_quant_config)

        candidate_id = (
            active_quant_config.quantification_id
            if active_quant_config is not None
            else LedgerModule._next_quantification_id(existing_ledger)
        )
        candidate_config = LedgerModule._build_quantification_config(self, 
            quantification_id=candidate_id,
            sample_elements=current_sample_elements,
            substrate_elements=current_substrate_elements,
            options=current_options,
            reference_values_by_el_line=current_reference_values,
        )

        if (
            not force_new
            and active_quant_config is not None
            and LedgerModule._quantification_configs_match(active_quant_config, candidate_config)
        ):
            self.current_quant_config = active_quant_config
            self.current_quantification_id = active_quant_config.quantification_id
            LedgerModule._update_non_fingerprint_options(self.current_quant_config, candidate_config)
            # Even when reusing quantification settings, clustering inputs (e.g., ref_formulae
            # passed via runners as sample['cnd']) may have changed for this run.
            LedgerModule._ensure_current_clustering_run(self, 
                self.current_quant_config,
            )
            LedgerModule._apply_active_clustering_config(self, self.current_quant_config)
            return

        # Same science as a prior (non-active) run: reuse the latest matching fingerprint.
        if not force_new and existing_ledger is not None:
            matching_prior = [
                config
                for config in existing_ledger.quantifications
                if LedgerModule._quantification_configs_match(config, candidate_config)
            ]
            if matching_prior:
                latest_match = max(matching_prior, key=lambda config: config.quantification_id)
                self.current_quant_config = latest_match
                self.current_quantification_id = latest_match.quantification_id
                LedgerModule._update_non_fingerprint_options(self.current_quant_config, candidate_config)
                LedgerModule._ensure_current_clustering_run(self, self.current_quant_config)
                LedgerModule._apply_active_clustering_config(self, self.current_quant_config)
                return

        if not force_new and active_quant_config is not None:
            changes = active_quant_config.fingerprint_differences(candidate_config)
            if changes:
                changed_summary = LedgerModule._format_quantification_config_changes(changes)
                if not active_quant_has_results:
                    logger.info(
                        "🆕 Quantification scientific fingerprint changed; "
                        "replacing unused active quantification config.\n%s",
                        changed_summary,
                    )
                    self.current_quant_config = candidate_config
                    self.current_quantification_id = candidate_config.quantification_id
                    LedgerModule._ensure_current_clustering_run(self, 
                        self.current_quant_config,
                    )
                    LedgerModule._apply_active_clustering_config(self, self.current_quant_config)
                    return
                logger.info(
                    "🆕 Quantification scientific fingerprint changed; "
                    "creating a new quantification config.\n%s",
                    changed_summary,
                )

        quantification_id = LedgerModule._next_quantification_id(existing_ledger)
        self.current_quant_config = LedgerModule._build_quantification_config(self, 
            quantification_id=quantification_id,
            sample_elements=current_sample_elements,
            substrate_elements=current_substrate_elements,
            options=current_options,
            reference_values_by_el_line=current_reference_values,
        )
        self.current_quantification_id = quantification_id
        LedgerModule._ensure_current_clustering_run(self, self.current_quant_config)
        LedgerModule._apply_active_clustering_config(self, self.current_quant_config)

    def _build_quantification_options(self: Any) -> Dict[str, Any]:
        """Build the subset of quantification options that defines result reuse."""
        options: Dict[str, Any] = {
            "method": self.quant_cfg.method,
            "spectrum_lims": [
                float(self.quant_cfg.spectrum_lims[0]),
                float(self.quant_cfg.spectrum_lims[1]),
            ],
            "fit_tolerance": float(self.quant_cfg.fit_tolerance),
            "use_instrument_background": bool(self.quant_cfg.use_instrument_background),
            # Effective quantifier state values that materially affect quantification output.
            "is_particle": bool(self._apply_geom_factors),
            "beam_energy_keV": float(self.measurement_cfg.beam_energy_keV),
            "emergence_angle": float(self.measurement_cfg.emergence_angle),
            "det_ch_offset": float(self.det_ch_offset),
            "det_ch_width": float(self.det_ch_width),
            "min_total_counts_fraction": self._min_total_counts_fraction(),
            "is_known_phase_mixture": bool(getattr(self.quant_cfg, "is_known_phase_mixture", False)),
            "identify_missing_elements": bool(getattr(self.quant_cfg, "identify_missing_elements", dflt.identify_missing_elements)),
            "identification_seed_spectra": int(
                getattr(self.quant_cfg, "identification_seed_spectra", dflt.identification_seed_spectra)
            ),
            "identification_z_range": [
                int(z) for z in getattr(self.quant_cfg, "identification_z_range", dflt.identification_z_range)
            ],
            "identification_excluded_elements": sorted(
                getattr(self.quant_cfg, "identification_excluded_elements", dflt.identification_excluded_elements)
            ),
        }
        return options

    def _get_reference_values_by_el_line(
        self: Any,
        active_quant_config: Optional[QuantificationConfig] = None,
    ) -> Dict[str, Any]:
        """Return method-dependent reference values for the scientific fingerprint.

        Only sample-element reference lines are stored and compared. Substrate
        elements are excluded from ``reference_values_by_el_line``.

        Always loads standards when needed so live values for sample elements are
        compared against any cached active-config values. Cached values are reused
        only when they match the live standards within tolerance.
        """
        # load_if_missing=True is required: the default quant path leaves
        # standards_dict=None until XSp_Quantifier loads standards, so a False
        # flag would skip the comparison and always reuse stale cached refs.
        sample_elements = set(self.all_els_sample)
        current_reference_values = LedgerModule._extract_reference_values_from_standards(self, load_if_missing=True)

        if active_quant_config is not None and active_quant_config.reference_values_by_el_line:
            cached_reference_values = LedgerModule._filter_reference_values_for_elements(
                active_quant_config.reference_values_by_el_line,
                relevant_elements=sample_elements,
            )
            if not current_reference_values:
                return cached_reference_values

            if LedgerModule._reference_values_changed(
                current_reference_values,
                cached_reference_values,
                decimals=2,
                relevant_elements=sample_elements,
            ):
                logger.info(
                    "🆕 Reference values in standards differ from active quantification config "
                    "(beyond 2-decimal tolerance) for sample elements; "
                    "a new quantification config will be created."
                )
                return current_reference_values
            return cached_reference_values

        return current_reference_values

    @staticmethod
    def _filter_reference_values_for_elements(
        reference_values_by_el_line: Dict[str, Any],
        relevant_elements: Optional[set[str]] = None,
    ) -> Dict[str, Any]:
        """Keep only reference-line entries whose element is in ``relevant_elements``."""
        if relevant_elements is None:
            return dict(sorted((str(k), v) for k, v in reference_values_by_el_line.items()))

        filtered: Dict[str, Any] = {}
        for el_line, value in reference_values_by_el_line.items():
            element = str(el_line).split("_", maxsplit=1)[0]
            if element in relevant_elements:
                filtered[str(el_line)] = value
        return dict(sorted(filtered.items()))

    @staticmethod
    def _rounded_reference_values(
        reference_values_by_el_line: Dict[str, Any],
        decimals: int,
        relevant_elements: Optional[set[str]] = None,
    ) -> Dict[str, float]:
        """Return deterministic rounded reference values keyed by element-line."""
        filtered = LedgerModule._filter_reference_values_for_elements(
            reference_values_by_el_line,
            relevant_elements=relevant_elements,
        )
        rounded: Dict[str, float] = {}
        for el_line, value in filtered.items():
            rounded[str(el_line)] = round(float(value), decimals)
        return dict(sorted(rounded.items()))

    @staticmethod
    def _reference_values_changed(
        current_reference_values: Dict[str, Any],
        cached_reference_values: Dict[str, Any],
        decimals: int,
        relevant_elements: Optional[set[str]] = None,
    ) -> bool:
        """Return True when reference values differ beyond configured decimal tolerance."""
        rounded_current = LedgerModule._rounded_reference_values(
            current_reference_values,
            decimals,
            relevant_elements=relevant_elements,
        )
        rounded_cached = LedgerModule._rounded_reference_values(
            cached_reference_values,
            decimals,
            relevant_elements=relevant_elements,
        )
        return rounded_current != rounded_cached

    def _extract_reference_values_from_standards(self: Any, load_if_missing: bool) -> Dict[str, Any]:
        """Extract per-reference-line values from standards for the configured quantification method."""
        standards_by_line = None

        if self.quant_cfg.use_project_specific_std_dict:
            if self.standards_dict is not None and not load_if_missing:
                standards_by_line = self.standards_dict
            else:
                standards, _ = StandardsModule._load_standards(self)
                standards_by_line = standards.standards_by_mode.get(self.measurement_cfg.mode, {})
        else:
            standards_by_line = self.standards_dict
            if standards_by_line is None and self.standards is not None:
                if isinstance(self.standards, EDSStandardsFile):
                    standards_by_line = self.standards.standards_by_mode.get(
                        self.measurement_cfg.mode,
                        {},
                    )
                elif self.measurement_cfg.mode in self.standards:
                    standards_by_line = self.standards[self.measurement_cfg.mode]
                else:
                    standards_by_line = self.standards
            if standards_by_line is None and load_if_missing:
                standards, _ = StandardsModule._load_standards(self)
                standards_by_line = standards.standards_by_mode.get(self.measurement_cfg.mode, {})

        if standards_by_line is None:
            return {}

        # Persist/compare only sample-element reference lines (never substrate).
        relevant_elements = set(self.detectable_els_sample)
        reference_line_suffixes = tuple(XSp_Quantifier.xray_quant_ref_lines)
        reference_values_by_el_line: Dict[str, Any] = {}
        if self.quant_cfg.method == "PB":
            for el_line, std_values in standards_by_line.items():
                element, _, line_suffix = str(el_line).partition("_")
                if element not in relevant_elements or line_suffix not in reference_line_suffixes:
                    continue

                if hasattr(std_values, "reference_mean"):
                    reference_mean = getattr(std_values, "reference_mean", None)
                    if reference_mean is None:
                        continue
                    reference_values_by_el_line[str(el_line)] = float(reference_mean.corrected_pb)
                    continue

                if isinstance(std_values, dict):
                    reference_mean = std_values.get("reference_mean")
                    if isinstance(reference_mean, dict) and cnst.COR_PB_DF_KEY in reference_mean:
                        reference_values_by_el_line[str(el_line)] = float(reference_mean[cnst.COR_PB_DF_KEY])
                        continue
                    std_values = std_values.get("entries", [])

                mean_std = next(
                    (std for std in std_values if std.get(cnst.STD_ID_KEY) == cnst.STD_MEAN_ID_KEY),
                    None,
                )
                if mean_std is None or cnst.COR_PB_DF_KEY not in mean_std:
                    continue
                reference_values_by_el_line[str(el_line)] = float(mean_std[cnst.COR_PB_DF_KEY])

        return dict(sorted(reference_values_by_el_line.items()))

    @staticmethod
    def _get_active_quantification_config(existing_ledger: Optional[SampleLedger]) -> Optional[QuantificationConfig]:
        """Return the active quantification config from the ledger when available."""
        return get_active_quantification_config(existing_ledger)

    def _build_quantification_config(
        self: Any,
        quantification_id: int,
        sample_elements: List[str],
        substrate_elements: List[str],
        options: Dict[str, Any],
        reference_values_by_el_line: Dict[str, Any],
    ) -> QuantificationConfig:
        """Build the persisted descriptor for the current quantification run."""
        label = f"Quant {quantification_id} {self.output_filename_suffix}"
        if isinstance(self.output_filename_suffix, str) and self.output_filename_suffix:
            label += self.output_filename_suffix
        # Never persist substrate reference lines in the quantification config.
        sample_reference_values = LedgerModule._filter_reference_values_for_elements(
            reference_values_by_el_line,
            relevant_elements=set(sample_elements),
        )
        reference_lines_by_element = QuantificationConfig.derive_reference_lines_by_element(
            reference_values_by_el_line=sample_reference_values,
            preferred_lines=XSp_Quantifier.xray_quant_ref_lines,
        )
        return QuantificationConfig(
            quantification_id=quantification_id,
            label=label,
            autoemx_version=AUTOEMX_VERSION,
            sample_elements=sample_elements,
            substrate_elements=substrate_elements,
            els_w_fr=LedgerModule._get_forced_mass_fractions(self),
            options=options,
            reference_values_by_el_line=sample_reference_values,
            reference_lines_by_element=reference_lines_by_element,
            clustering_analyses=[],
            active_clustering_analysis_index=None,
        )

    def _get_forced_mass_fractions(self: Any) -> Optional[Dict[str, float]]:
        """Return forced elemental mass fractions for the active run, if any."""
        forced = getattr(self.sample_cfg, "w_frs", None)
        if isinstance(forced, dict) and forced:
            return {str(el): float(fr) for el, fr in forced.items()}

        return None

    def _upsert_current_quantification_config_on_ledger(self: Any, ledger: SampleLedger) -> None:
        """Insert or replace the current quantification config inside the provided ledger.

        Same scientific fingerprint with the same quantification_id updates that slot in place.
        Same fingerprint with a different id is treated as a distinct run (no id remapping),
        so force_requantification can keep a newly allocated id.
        """
        if self.current_quant_config is None:
            return

        # Same science and same id: refresh that slot (clustering history, etc.).
        for idx, config in enumerate(ledger.quantifications):
            if (
                config.quantification_id == self.current_quant_config.quantification_id
                and LedgerModule._quantification_configs_match(config, self.current_quant_config)
            ):
                ledger.quantifications[idx] = self.current_quant_config
                self.current_quantification_id = self.current_quant_config.quantification_id
                ledger.active_quant = self.current_quant_config.quantification_id
                return

        # Upsert by quantification_id only (do not collapse distinct ids onto fingerprint matches).
        existing_index = next(
            (
                i
                for i, config in enumerate(ledger.quantifications)
                if config.quantification_id == self.current_quant_config.quantification_id
            ),
            None,
        )

        if existing_index is None:
            ledger.append_quantification_config(self.current_quant_config)
        else:
            ledger.quantifications[existing_index] = self.current_quant_config
        ledger.active_quant = self.current_quant_config.quantification_id

    def _persist_current_quantification_config(self: Any) -> None:
        """Persist the currently selected quantification config as the ledger active run."""
        if self.current_quant_config is None or self.current_quantification_id is None:
            return

        ledger = LedgerModule._load_or_create_ledger(self)
        LedgerModule._upsert_current_quantification_config_on_ledger(self, ledger)

        # Ensure the full current config payload (including clustering history/index)
        # is written for the active quantification id.
        for idx, config in enumerate(ledger.quantifications):
            if config.quantification_id == self.current_quantification_id:
                ledger.quantifications[idx] = self.current_quant_config.model_copy(deep=True)
                break

        ledger.active_quant = self.current_quantification_id
        ledger.to_json_file(LedgerModule._get_ledger_path(self))

    def _sync_existing_quantification_from_ledger(self: Any) -> None:
        """Populate in-memory quantification slots from ledger results with the current id."""
        if self.current_quantification_id is None:
            return
        existing_ledger = LedgerModule._load_or_create_ledger(self)
        if existing_ledger is None:
            return

        total_spectra = len(existing_ledger.spectra)
        LedgerModule._ensure_quant_tracking_length(self, total_spectra)

        for index, spectrum in enumerate(existing_ledger.spectra):
            matching = next(
                (
                    result for result in spectrum.quantification_results
                    if result.quantification_id == self.current_quantification_id
                ),
                None,
            )
            if matching is None:
                continue

            self.spectra_quant_records[index] = matching

    def _persist_quantification_record(self: Any, spectrum_index: int, quant_record: QuantificationResult, overwrite: bool = False) -> None:
        """Write one quantification result to the ledger immediately for interruption-safe progress.

        Parameters
        ----------
        overwrite : bool
            If True and a result with the same quantification_id already exists for this spectrum,
            replace it (used for force_requantification, requantify_only_unquantified_spectra,
            and interrupted re-runs).
        """
        existing_ledger = LedgerModule._load_or_create_ledger(self)
        ledger = existing_ledger

        if self.current_quant_config is None:
            raise ValueError("Current quantification config is not initialized")

        LedgerModule._upsert_current_quantification_config_on_ledger(self, ledger)

        ledger.active_quant = self.current_quantification_id

        existing_ids = {
            existing.quantification_id
            for existing in ledger.spectra[spectrum_index].quantification_results
        }
        if quant_record.quantification_id not in existing_ids:
            ledger.append_quantification_result(spectrum_index, quant_record)
        elif overwrite:
            ledger.spectra[spectrum_index].quantification_results = [
                quant_record if r.quantification_id == quant_record.quantification_id else r
                for r in ledger.spectra[spectrum_index].quantification_results
            ]

        ledger.to_json_file(LedgerModule._get_ledger_path(self))

    def _prepare_acquisition_resume_state(self: Any, log_resume_info: bool = True) -> dict:
        """Load the ledger, reconcile spectra pointers, and resolve acquisition resume counters.

        Parameters
        ----------
        log_resume_info : bool, optional
            If True, emit resume/progress messages when prior spectra are detected.

        Returns
        -------
        dict
            Keys: ``tot_n_spectra``, ``next_spectrum_id``, ``particle_id_offset``,
            ``had_ledger_before``, and ``acquisition_needed``.
        """
        _had_ledger_before = LedgerModule._load_existing_ledger(self) is not None
        _existing_ledger = LedgerModule._load_or_create_ledger(self)
        if _existing_ledger is not None and _existing_ledger.spectra:
            tot_n_spectra = len(_existing_ledger.spectra)
            _numeric_ids = [
                int(e.spectrum_id)
                for e in _existing_ledger.spectra
                if e.spectrum_id is not None and str(e.spectrum_id).isdigit()
            ]
            next_spectrum_id = (max(_numeric_ids) + 1) if _numeric_ids else 0
            _particle_ids = [
                e.acquisition_details.particle_id
                for e in _existing_ledger.spectra
                if e.acquisition_details is not None
                and e.acquisition_details.particle_id is not None
            ]
            _particle_id_offset = 0
            self.particle_cntr = -1
            try:
                if _particle_ids:
                    _particle_id_offset = max(int(p) for p in _particle_ids)
                    self.particle_cntr = _particle_id_offset
                else:
                    inferred_particle_id = LedgerModule._infer_max_particle_id_from_saved_images(self)
                    if inferred_particle_id is not None:
                        _particle_id_offset = inferred_particle_id
                        self.particle_cntr = inferred_particle_id
            except (ValueError, TypeError):
                _particle_id_offset = 0
                self.particle_cntr = -1
            next_particle_id = _particle_id_offset + 1 if self.sample_cfg.is_particle_acquisition else "n/a"
            if log_resume_info and self.verbose:
                n_ingested = getattr(self, "_last_acq_ledger_ingested_spectra_count", 0)
                if (
                    not _had_ledger_before
                    and getattr(self, "_last_acq_ledger_created", False)
                    and n_ingested > 0
                ):
                    logger.info(
                        "ℹ️ No ledger found; %d pre-existing spectrum file/s ingested into a new ledger "
                        "(assuming the same acquisition configurations). "
                        "%d spectrum/spectra available; restarting from spectrum #%d, particle #%s. "
                        "If a fresh acquisition is desired, delete the sample folder at: %s",
                        n_ingested,
                        tot_n_spectra,
                        next_spectrum_id,
                        next_particle_id,
                        os.path.abspath(self.sample_result_dir),
                    )
                else:
                    ledger_note = (
                        "Existing ledger detected; acquisition will resume"
                        if _had_ledger_before
                        else "Prior spectra detected"
                    )
                    logger.info(
                        "ℹ️ %s. %d spectrum/spectra available; restarting from spectrum #%d, particle #%s. "
                        "If a fresh acquisition is desired, delete the sample folder at: %s",
                        ledger_note,
                        tot_n_spectra,
                        next_spectrum_id,
                        next_particle_id,
                        os.path.abspath(self.sample_result_dir),
                    )
        else:
            tot_n_spectra = 0
            next_spectrum_id = 0
            _particle_id_offset = 0

        return {
            "tot_n_spectra": tot_n_spectra,
            "next_spectrum_id": next_spectrum_id,
            "particle_id_offset": _particle_id_offset,
            "had_ledger_before": _had_ledger_before,
            "acquisition_needed": tot_n_spectra < self.max_n_spectra,
        }

    def _initialise_acquisition_ledger(self: Any) -> None:
        """Initialize acquisition ledger and reconcile any existing spectra pointers.

        This guarantees that if a ledger is missing but spectra files already exist
        under ``spectra/``, those spectra are ingested into the newly created ledger
        before acquisition restarts. Also determines whether new spectra must be
        acquired so microscope initialization can be skipped when unnecessary.
        """
        resume_state = LedgerModule._prepare_acquisition_resume_state(self)
        self._acquisition_resume_state = resume_state
        self._spectrum_acquisition_needed = resume_state["acquisition_needed"]
        self.EM_controller = None

        if self.verbose and not self._spectrum_acquisition_needed:
            logger.info(
                "ℹ️ %d spectrum/spectra already acquired (maximum %d). Skipping acquisition.",
                "If more spectra are desired, increase max_n_spectra and restart acquisition.",
                resume_state["tot_n_spectra"],
                self.max_n_spectra,
            )

    def _resolve_active_analysis_config_ids(self: Any) -> Tuple[int, int]:
        """Resolve active quantification and clustering config ids for analysis folder naming."""
        quant_config = self.current_quant_config

        if quant_config is None:
            ledger = LedgerModule._load_or_create_ledger(self)
            quant_config = LedgerModule._get_active_quantification_config(ledger)

        if quant_config is None:
            raise RuntimeError("No active quantification config is available to name analysis directory")

        active_clustering_analysis = quant_config.get_active_clustering_analysis()
        if active_clustering_analysis is None:
            raise RuntimeError("No active clustering config is available to name analysis directory")

        return int(quant_config.quantification_id), int(active_clustering_analysis.config.clustering_id)

    def _infer_max_particle_id_from_saved_images(self: Any) -> Optional[int]:
        """Infer the maximum particle id from saved acquisition image filenames."""
        images_dir = Path(self.sample_result_dir, cnst.IMAGES_DIR)
        if not images_dir.exists():
            return None

        max_particle_id = None
        particle_pattern = re.compile(r"_par(\d+)_")
        for image_path in images_dir.iterdir():
            if not image_path.is_file():
                continue
            match = particle_pattern.search(image_path.name)
            if match is None:
                continue
            try:
                particle_id = int(match.group(1))
            except (TypeError, ValueError):
                continue
            if max_particle_id is None or particle_id > max_particle_id:
                max_particle_id = particle_id

        return max_particle_id

    def _apply_active_clustering_config(self: Any, quant_config: QuantificationConfig) -> None:
        """Load runtime clustering options from the active clustering config of a quantification config."""
        active_clustering_analysis = quant_config.get_active_clustering_analysis()
        if active_clustering_analysis is None:
            raise ValueError(
                "Active quantification config does not contain an active clustering config"
            )

        self.clustering_cfg = active_clustering_analysis.config.model_copy(deep=True)
        self.ref_formulae = list(dict.fromkeys(self.clustering_cfg.ref_formulae))
        self._calc_reference_phases_df()

    def _build_clustering_config_descriptor(self: Any, clustering_id: int) -> LedgerClusteringConfig:
        """Build a persisted clustering config descriptor for the current analysis settings."""
        return LedgerClusteringConfig(
            clustering_id=clustering_id,
            method=self.clustering_cfg.method,
            features=self.clustering_cfg.features,
            k_forced=self.clustering_cfg.k_forced,
            k_resolved=self.clustering_cfg.k_resolved,
            k_finding_method=self.clustering_cfg.k_finding_method,
            max_k=self.clustering_cfg.max_k,
            ref_formulae=list(self.clustering_cfg.ref_formulae),
            do_matrix_decomposition=self.clustering_cfg.do_matrix_decomposition,
            max_analytical_error_percent=self.clustering_cfg.max_analytical_error_percent,
            min_bckgrnd_cnts=self.clustering_cfg.min_bckgrnd_cnts,
            quant_flags_accepted=list(self.clustering_cfg.quant_flags_accepted),
            dbscan=self.clustering_cfg.dbscan.model_copy(deep=True),
            geometry=self.clustering_cfg.geometry,
            aitchison=self.clustering_cfg.aitchison.model_copy(deep=True),
            mixture=self.clustering_cfg.mixture.model_copy(deep=True),
            auto_merge_clusters=self.clustering_cfg.auto_merge_clusters,
            cluster_merge=self.clustering_cfg.cluster_merge.model_copy(deep=True),
        )

    def _ensure_current_clustering_run(
        self: Any,
        quant_config: QuantificationConfig,
    ) -> None:
        """Ensure active quantification config tracks clustering run history and active config."""
        candidate_clustering_config = LedgerModule._build_clustering_config_descriptor(self, 
            clustering_id=LedgerModule._next_clustering_config_id(quant_config)
        )
        matching_index: Optional[int] = None
        runtime_clustering_id = getattr(self.clustering_cfg, "clustering_id", None)
        for idx, analysis in enumerate(quant_config.clustering_analyses):
            if not LedgerModule._clustering_configs_match(analysis.config, candidate_clustering_config):
                continue
            if runtime_clustering_id is not None and analysis.config.clustering_id == runtime_clustering_id:
                matching_index = idx
                break
            if matching_index is None:
                matching_index = idx

        if matching_index is not None:
            quant_config.active_clustering_analysis_index = matching_index
            return

        active_clustering_analysis = quant_config.get_active_clustering_analysis()
        active_clustering_config = (
            active_clustering_analysis.config if active_clustering_analysis is not None else None
        )
        if active_clustering_config is not None:
            changes = active_clustering_config.fingerprint_differences(candidate_clustering_config)
            if changes:
                changed_summary = LedgerModule._format_quantification_config_changes(changes)
                logger.info(
                    "🆕 Clustering scientific fingerprint changed; "
                    "appending a new clustering config to the active quantification config.\n%s",
                    changed_summary,
                )

        quant_config.clustering_analyses.append(
            ClusteringAnalysis(config=candidate_clustering_config, result=None)
        )
        quant_config.active_clustering_analysis_index = len(quant_config.clustering_analyses) - 1

    def _persist_resolved_k_on_active_clustering_config(self: Any, resolved_k: Optional[int]) -> None:
        """Persist resolved k for the active clustering config of the active quantification run."""
        if resolved_k is None:
            return

        resolved_k_int = int(resolved_k)
        self.clustering_cfg.k_resolved = resolved_k_int

        if self.current_quant_config is not None:
            active_clustering_analysis = self.current_quant_config.get_active_clustering_analysis()
            active_clustering_config = (
                active_clustering_analysis.config if active_clustering_analysis is not None else None
            )
            if active_clustering_config is not None and active_clustering_config.k_resolved != resolved_k_int:
                active_clustering_config.k_resolved = resolved_k_int
                LedgerModule._persist_current_quantification_config(self)
            return

        ledger = LedgerModule._load_or_create_ledger(self)
        active_quant_config = LedgerModule._get_active_quantification_config(ledger)
        if active_quant_config is None:
            return

        active_clustering_analysis = active_quant_config.get_active_clustering_analysis()
        active_clustering_config = (
            active_clustering_analysis.config if active_clustering_analysis is not None else None
        )
        if active_clustering_config is None or active_clustering_config.k_resolved == resolved_k_int:
            return

        active_clustering_config.k_resolved = resolved_k_int
        ledger.to_json_file(LedgerModule._get_ledger_path(self))

    def _persist_clustering_result_on_active_analysis(self: Any, clustering_result: ClusteringResult) -> None:
        """Persist the latest clustering result on the active clustering analysis entry."""
        ledger = LedgerModule._load_or_create_ledger(self)
        active_quant_config = LedgerModule._get_active_quantification_config(ledger)
        if active_quant_config is None:
            return

        # Ensure the current runtime clustering settings map to a concrete analysis entry
        # on the same ledger object being persisted.
        LedgerModule._ensure_current_clustering_run(self, active_quant_config)

        active_quant_config.set_active_clustering_result(clustering_result.model_copy(deep=True))

        if (
            self.current_quant_config is not None
            and self.current_quant_config.quantification_id == active_quant_config.quantification_id
        ):
            self.current_quant_config.set_active_clustering_result(clustering_result.model_copy(deep=True))

        ledger.to_json_file(LedgerModule._get_ledger_path(self))

    def _save_experimental_config(self: Any, is_XSp_measurement) -> None:
        """
        Save all relevant configuration dataclasses and metadata related to the
        current spectrum collection/acquisition to a JSON file.
    
        The saved file includes:
            - Timestamp of saving
            - All configuration dataclasses
    
        This function is intended to be called after acquisition to ensure
        reproducibility and traceability of the experimental configuration.
    
        Raises
        ------
        OSError
            If the output directory cannot be created or file cannot be written.
        """
    
        now = datetime.now()

        def _serialize_config(config_obj: Any) -> Dict[str, Any]:
            if hasattr(config_obj, "model_dump"):
                return config_obj.model_dump(mode="json")
            return asdict(config_obj)
    
        # Gather configuration dataclasses as dictionaries
        self.measurement_cfg = self.measurement_cfg.model_copy(
            update={
                "powder_meas_cfg": self.powder_meas_cfg,
                "bulk_meas_cfg": self.bulk_meas_cfg,
                "exp_stds_cfg": self.exp_stds_cfg if self.exp_stds_cfg.is_exp_std_measurement else None,
            }
        )
        cfg_dataclasses = {
            cnst.SAMPLE_CFG_KEY: _serialize_config(self.sample_cfg),
            cnst.MICROSCOPE_CFG_KEY: _serialize_config(self.microscope_cfg),
            cnst.MEASUREMENT_CFG_KEY: _serialize_config(self.measurement_cfg),
            cnst.SAMPLESUBSTRATE_CFG_KEY: _serialize_config(self.sample_substrate_cfg),
        }
        
        if is_XSp_measurement:
            cfg_dataclasses[cnst.QUANTIFICATION_CFG_KEY] = _serialize_config(self.quant_cfg)
    
        # Include dataclasses corresponding to measurement type
        if is_XSp_measurement and not self.exp_stds_cfg.is_exp_std_measurement:
            cfg_dataclasses[cnst.CLUSTERING_CFG_KEY] = LedgerModule._runtime_clustering_cfg_payload(self.clustering_cfg)
            cfg_dataclasses[cnst.PLOT_CFG_KEY] = _serialize_config(self.plot_cfg)


    
        # Compose the metadata dictionary for saving
        spectrum_collection_info = {
            cnst.DATETIME_KEY: now.strftime("%Y-%m-%d %H:%M:%S"),
            **cfg_dataclasses
        }
    
        # Ensure the output directory exists before saving
        try:
            os.makedirs(self.sample_result_dir, exist_ok=True)
        except Exception as e:
            raise OSError(f"Could not create output directory '{self.sample_result_dir}': {e}")
    
        try:
            ledger = LedgerModule._load_or_create_ledger(self)
            ledger.configs = LedgerConfigs.model_validate({
                cnst.MICROSCOPE_CFG_KEY: spectrum_collection_info[cnst.MICROSCOPE_CFG_KEY],
                cnst.SAMPLE_CFG_KEY: spectrum_collection_info[cnst.SAMPLE_CFG_KEY],
                cnst.MEASUREMENT_CFG_KEY: spectrum_collection_info[cnst.MEASUREMENT_CFG_KEY],
                cnst.SAMPLESUBSTRATE_CFG_KEY: spectrum_collection_info[cnst.SAMPLESUBSTRATE_CFG_KEY],
                cnst.PLOT_CFG_KEY: spectrum_collection_info.get(cnst.PLOT_CFG_KEY, self.plot_cfg.model_dump(mode="json")),
            })

            # Persist an initial quantification+clustering descriptor at acquisition start
            # so user-provided clustering inputs (e.g., sample['cnd'] -> ref_formulae)
            # are saved to ledger even before any quantification run is executed.
            if is_XSp_measurement and not self.exp_stds_cfg.is_exp_std_measurement and not ledger.quantifications:
                quantification_id = LedgerModule._next_quantification_id(ledger)
                self.current_quant_config = LedgerModule._build_quantification_config(self, 
                    quantification_id=quantification_id,
                    sample_elements=list(self.all_els_sample),
                    substrate_elements=list(self.all_els_substrate),
                    options=LedgerModule._build_quantification_options(self),
                    reference_values_by_el_line=LedgerModule._get_reference_values_by_el_line(self, active_quant_config=None),
                )
                self.current_quantification_id = quantification_id
                LedgerModule._ensure_current_clustering_run(self, self.current_quant_config)
                ledger.append_quantification_config(self.current_quant_config)
                ledger.active_quant = quantification_id

            ledger.to_json_file(LedgerModule._get_ledger_path(self))
        except Exception as e:
            raise OSError(f"Could not persist configurations to ledger '{LedgerModule._get_ledger_path(self)}': {e}")
            

    def _sync_active_quantification_from_ledger(self: Any) -> None:
        """
        Sync the in-memory state of the analyzer with the active quantification run of the ledger, before analysis.

        This is essential when analyse_data is called without a preceding run_quantification (i.e. analysis-only
        path): the typed QuantificationResult records must be populated in memory, and the elements and spectrum
        limits are those recorded for the active quantification run.
        """
        if self.sample_result_dir is None:
            return
        existing_ledger = LedgerModule._load_or_create_ledger(self)
        # Resolve current_quantification_id from the ledger when not already set
        # (analysis-only callers skip run_quantification so it is never assigned).
        if self.current_quantification_id is None and existing_ledger is not None:
            if existing_ledger.active_quant is not None:
                self.current_quantification_id = existing_ledger.active_quant

        # Analysis-only runs must also register the current clustering inputs on the
        # active quantification config so clustering history/index stays in sync.
        active_quant_config = LedgerModule._get_active_quantification_config(existing_ledger)
        if active_quant_config is not None:
            self.current_quant_config = active_quant_config
            self.current_quantification_id = active_quant_config.quantification_id
            # Use the elements and spectrum limits recorded for the active quantification
            # run, falling back to the base sample config when they are not stored.
            self.all_els_sample = list(dict.fromkeys(
                active_quant_config.sample_elements or self.sample_cfg.elements
            ))
            self.detectable_els_sample = [el for el in self.all_els_sample if el not in calibs.undetectable_els]
            self.all_els_substrate = list(dict.fromkeys(
                active_quant_config.substrate_elements or self.sample_substrate_cfg.elements
            ))
            detectable_els_substrate = [el for el in self.all_els_substrate if el not in calibs.undetectable_els]
            self.detectable_els_substrate = [el for el in detectable_els_substrate if el not in self.detectable_els_sample]
            active_spectrum_lims = (active_quant_config.options or {}).get("spectrum_lims")
            if active_spectrum_lims is not None:
                self.sp_start = int(round(active_spectrum_lims[0]))
                self.sp_end = int(round(active_spectrum_lims[1]))
                if self.det_ch_offset is not None and self.det_ch_width:
                    self.energy_vals = np.array([
                        self.det_ch_offset + self.det_ch_width * i
                        for i in range(self.sp_start, self.sp_end)
                    ])
            LedgerModule._ensure_current_clustering_run(self, self.current_quant_config)
            LedgerModule._persist_current_quantification_config(self)
            LedgerModule._apply_active_clustering_config(self, self.current_quant_config)

        tot = len(existing_ledger.spectra) if existing_ledger is not None else 0
        if tot > 0:
            LedgerModule._ensure_quant_tracking_length(self, tot)
            LedgerModule._sync_existing_quantification_from_ledger(self)
