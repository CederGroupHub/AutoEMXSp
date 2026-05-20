#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Experimental standards mixin for PB ratio workflows."""

import os
import shutil
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from pymatgen.core.composition import Composition

import autoemx.calibrations as calibs
import autoemx.utils.constants as cnst
from autoemx.config.schema_models import (
    EDSStandardsFile,  # type: ignore
        ReferenceMean,  # type: ignore
    StandardEntry,  # type: ignore
    StandardFitLineResult,  # type: ignore
    StandardsFitResults,  # type: ignore
    StandardLine,  # type: ignore
    StandardMeanZ,  # type: ignore
)
from autoemx.core.quantifier import Quant_Corrections, XSp_Quantifier
from autoemx.utils import print_double_separator, print_single_separator, weight_to_atomic_fr

from autoemx._logging import get_logger  # type: ignore
logger = get_logger(__name__)


class StandardsModule:
    def _get_std_exports_dir(self: Any) -> str:
        """Return the folder used for user-facing standards CSV exports."""
        output_dir = getattr(self, "analysis_dir", None) or self.sample_result_dir
        os.makedirs(output_dir, exist_ok=True)
        return output_dir

    def _save_std_measurements_from_records(
        self: Any,
    ) -> Optional[pd.DataFrame]:
        """Write per-spectrum standards measurements CSV from quant records."""
        n_records = len(getattr(self, "spectra_quant_records", []))
        n_coords = len(getattr(self, "sp_coords", []))
        n_spectra = max(n_records, n_coords)
        if n_spectra < 1:
            return None

        ref_lines = set(XSp_Quantifier.xray_quant_ref_lines)

        rows: List[Dict[str, Any]] = []
        for i in range(n_spectra):
            record = self.spectra_quant_records[i] if i < n_records else None
            coords = self.sp_coords[i] if i < n_coords else {cnst.SP_ID_DF_KEY: str(i)}

            data_row: Dict[str, Any] = dict(coords)
            fit_pb_data: Dict[str, float] = {}

            if record is not None and record.fit_result is not None:
                for peak_key, peak in record.fit_result.fitted_peaks.items():
                    if peak.pb_ratio is None:
                        continue
                    if peak.line in ref_lines:
                        fit_pb_data[peak_key] = round(peak.pb_ratio,1)

            data_row.update(fit_pb_data)

            if record is not None and record.fit_result is not None:
                if record.quant_flag is not None:
                    data_row[cnst.QUANT_FLAG_DF_KEY] = record.quant_flag
                if record.comment is not None:
                    data_row[cnst.COMMENTS_DF_KEY] = record.comment
                if record.fit_result.r_squared is not None:
                    data_row[cnst.R_SQ_KEY] = float(f"{record.fit_result.r_squared:.5f}")
                if record.fit_result.reduced_chi_squared is not None:
                    data_row[cnst.REDCHI_SQ_KEY] = float(f"{record.fit_result.reduced_chi_squared:.1f}")

            rows.append(data_row)

        data_df = pd.DataFrame(rows)

        output_dir = StandardsModule._get_std_exports_dir(self)
        suffix = self.output_filename_suffix if isinstance(self.output_filename_suffix, str) else ""
        filename = f"{cnst.STDS_MEAS_FILENAME}" + suffix
        out_path = os.path.join(output_dir, filename + cnst.DATA_FILEEXT)

        data_df.to_csv(out_path, index=False, header=True)
        return data_df

    @staticmethod
    def _serialize_standard_mean_z(z_sample: Any) -> Optional[Dict[str, float]]:
        """Convert runtime Z summary payload to schema-compatible mean_z mapping."""
        if z_sample is None:
            return None

        if isinstance(z_sample, dict):
            required = [
                cnst.Z_MEAN_W_KEY,
                cnst.Z_MEAN_AT_KEY,
                cnst.Z_MEAN_STATHAM_KEY,
                cnst.Z_MEAN_MARKOWICZ_KEY,
            ]
            if all(key in z_sample for key in required):
                return {
                    cnst.Z_MEAN_W_KEY: float(z_sample[cnst.Z_MEAN_W_KEY]),
                    cnst.Z_MEAN_AT_KEY: float(z_sample[cnst.Z_MEAN_AT_KEY]),
                    cnst.Z_MEAN_STATHAM_KEY: float(z_sample[cnst.Z_MEAN_STATHAM_KEY]),
                    cnst.Z_MEAN_MARKOWICZ_KEY: float(z_sample[cnst.Z_MEAN_MARKOWICZ_KEY]),
                }
            return None

        if isinstance(z_sample, (list, tuple, np.ndarray)) and len(z_sample) >= 4:
            return {
                cnst.Z_MEAN_W_KEY: float(z_sample[0]),
                cnst.Z_MEAN_AT_KEY: float(z_sample[1]),
                cnst.Z_MEAN_STATHAM_KEY: float(z_sample[2]),
                cnst.Z_MEAN_MARKOWICZ_KEY: float(z_sample[3]),
            }

        return None

    def _compile_standards_from_references(self: Any) -> dict:
        standards, _ = StandardsModule._load_standards(self)
        std_dict_all_lines = standards.standards_by_mode[self.measurement_cfg.mode]
        ref_lines = XSp_Quantifier.xray_quant_ref_lines
        ref_formulae = self.clustering_cfg.ref_formulae

        filtered_std_dict = {}

        for el in self.detectable_els_sample:
            for line in ref_lines:
                el_line = f"{el}_{line}"
                if el_line not in std_dict_all_lines:
                    continue

                line_payload = std_dict_all_lines[el_line]
                ref_entries = []
                for i, std_entry in enumerate(line_payload.entries):
                    try:
                        if std_entry.formula is None:
                            continue
                        std_comp = Composition(std_entry.formula)
                        for ref_formula in ref_formulae:
                            if std_comp.reduced_formula == Composition(ref_formula).reduced_formula:
                                ref_entries += [i]
                    except Exception:
                        pass

                mean_payload = getattr(line_payload, "reference_mean", None)
                std_mean_value = mean_payload.corrected_pb if (mean_payload is not None and hasattr(mean_payload, "corrected_pb")) else None

                # Only use mean if it exists
                if len(ref_entries) < 1 and not self.exp_stds_cfg.is_exp_std_measurement:
                    if std_mean_value is None:
                        continue
                    ref_value = std_mean_value
                else:
                    new_std_ref_list = [std_d for i, std_d in enumerate(line_payload.entries) if i in ref_entries]
                    if len(new_std_ref_list) < 1:
                        if std_mean_value is None:
                            continue
                        ref_value = std_mean_value
                        filtered_std_dict[el_line] = [{
                            cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY,
                            cnst.COR_PB_DF_KEY: ref_value,
                        }]
                        continue
                    list_pb = [ref_line.corrected_pb for ref_line in new_std_ref_list]
                    ref_value = float(np.mean(list_pb))

                filtered_std_dict[el_line] = [{
                    cnst.STD_ID_KEY: cnst.STD_MEAN_ID_KEY,
                    cnst.COR_PB_DF_KEY: ref_value,
                }]

        return filtered_std_dict

    def _fit_stds_and_save_results(
        self: Any,
        run_fitting: bool = True,
    ) -> Optional[StandardsFitResults]:
        fit_results = None

        if run_fitting:
            self._fit_and_quantify_spectra(quantify=False)
        StandardsModule._save_std_measurements_from_records(self)

        std_ref_lines = StandardsModule._assemble_std_PB_data_from_records(self)
        if std_ref_lines:
            pb_corrected, z_sample = StandardsModule._calc_corrected_PB(self, std_ref_lines)
            fit_results = StandardsModule._save_std_results(self, std_ref_lines, pb_corrected, z_sample)
        return fit_results

    def _assemble_std_PB_data_from_records(self: Any) -> Dict[str, StandardFitLineResult]:
        """Assemble standards PB statistics directly from quantification records.

        This reads fit results from in-memory quantification records (which are
        ledger-backed) and no longer relies on a temporary DataFrame export.
        """
        accepted_flags = set(self.exp_stds_cfg.quant_flags_accepted)
        ref_lines = set(XSp_Quantifier.xray_quant_ref_lines)
        detectable = set(self.detectable_els_sample)

        pb_by_line: Dict[str, List[float]] = {}
        th_energy_by_line: Dict[str, float] = {}

        for record in getattr(self, "spectra_quant_records", []):
            if record is None:
                continue
            if record.quant_flag not in accepted_flags:
                continue
            if record.fit_result is None or not record.fit_result.fitted_peaks:
                continue

            for peak_key, peak in record.fit_result.fitted_peaks.items():
                if peak.line not in ref_lines:
                    continue
                if peak.element not in detectable:
                    continue
                if peak.pb_ratio is None:
                    continue

                if peak_key not in pb_by_line:
                    pb_by_line[peak_key] = []
                pb_by_line[peak_key].append(float(peak.pb_ratio))

                if peak.theoretical_energy is not None and peak_key not in th_energy_by_line:
                    th_energy_by_line[peak_key] = float(peak.theoretical_energy)

        std_ref_lines: Dict[str, StandardFitLineResult] = {}
        for el_line, meas_pb_ratios in pb_by_line.items():
            if len(meas_pb_ratios) < 1:
                continue

            if el_line in th_energy_by_line:
                th_energy = th_energy_by_line[el_line]
            elif el_line in self._th_peak_energies:
                th_energy = float(self._th_peak_energies[el_line])
            else:
                # Keep compatibility with previous behavior if energy cache is absent.
                th_energy = float("nan")

            std_ref_lines[el_line] = StandardFitLineResult.model_validate({
                "pb_ratios": meas_pb_ratios,
                "measured_pb": float(np.nanmean(meas_pb_ratios)),
                "stdev_pb": float(np.nanstd(meas_pb_ratios)),
                "peak_theoretical_energy_keV": float(th_energy),
                "corrected_pb": float("nan"),
                "rel_stdev_pb_percent": float("nan"),
                "n_spectra_used": 0,
            })

        return std_ref_lines

    def _evaluate_exp_std_fit(self: Any, tot_n_spectra: int) -> Tuple[int, bool, bool]:
        is_fit_successful = False
        is_converged = False
        num_valid_spectra = 0
        
        if self.verbose:
            print_double_separator()
            logger.info(f"🔬 Fitting after collection of {tot_n_spectra} spectra...")

        fit_results = StandardsModule._fit_stds_and_save_results(self)
        if fit_results is not None and fit_results.lines:
            is_fit_successful = True
            num_valid_spectra = int(np.min([line.n_spectra_used for line in fit_results.lines.values()]))
            print_single_separator()
            logger.info(f"📈 Acquired {num_valid_spectra} valid spectra out of {self.min_n_spectra} requested")
            is_converged = num_valid_spectra >= self.min_n_spectra
            if is_converged:
                logger.info(f"✅ Convergence reached, acquisition completed.")
            else:
                n_remaining = self.min_n_spectra - num_valid_spectra
                logger.info(f"⏳ {n_remaining} more valid spectra required.")

        return num_valid_spectra, is_fit_successful, is_converged

    def _assemble_std_PB_data(
        self: Any,
        data_df: "pd.DataFrame"
    ) -> Dict[str, Dict[str, Union[float, List[float]]]]:
        data_df = data_df.dropna(axis=1, how="all")
        if cnst.QUANT_FLAG_DF_KEY not in data_df.columns:
            raise RuntimeError(f"Missing required column '{cnst.QUANT_FLAG_DF_KEY}' in input DataFrame.")
        data_filtered_df = data_df[data_df[cnst.QUANT_FLAG_DF_KEY].isin(self.exp_stds_cfg.quant_flags_accepted)]

        all_fitted_el_lines = [el_line for el_line in self._th_peak_energies.keys() if el_line in data_filtered_df.columns]
        fitted_std_el_lines = [el_line for el_line in all_fitted_el_lines if el_line.split("_")[0] in self.detectable_els_sample]

        std_ref_lines = {}
        for el_line in fitted_std_el_lines:
            meas_pb_ratios = data_filtered_df[el_line].tolist()
            if len(meas_pb_ratios) > 0:
                std_ref_lines[el_line] = {
                    cnst.PB_RATIO_KEY: meas_pb_ratios,
                    cnst.MEAN_PB_KEY: float(np.nanmean(meas_pb_ratios)),
                    cnst.STDEV_PB_DF_KEY: float(np.nanstd(meas_pb_ratios)),
                    cnst.PEAK_TH_ENERGY_KEY: self._th_peak_energies[el_line],
                }

        return std_ref_lines

    def _calc_corrected_PB(
        self: Any,
        std_ref_lines: Dict[str, StandardFitLineResult]
    ) -> Tuple[np.ndarray, Any]:
        peak_energies_dict: Dict[str, float] = {}
        means_pb: List[float] = []
        w_frs: List[float] = []

        for el_line, line_result in std_ref_lines.items():
            peak_energies_dict[el_line] = float(line_result.peak_theoretical_energy_keV)
            means_pb.append(float(line_result.measured_pb))
            el = el_line.split('_')[0]
            if el not in self.exp_stds_cfg.w_frs:
                raise RuntimeError(f"Mass fraction for element '{el}' not found in exp_stds_cfg.w_frs.")
            w_frs.append(float(self.exp_stds_cfg.w_frs[el]))

        zaf_calculator = Quant_Corrections(
            elements=self.detectable_els_sample,
            beam_energy=self.measurement_cfg.beam_energy_keV,
            emergence_angle=self.measurement_cfg.emergence_angle,
            meas_mode=self.measurement_cfg.mode,
            verbose=False,
        )

        missing_elements = [el for el in self.detectable_els_sample if el not in self.exp_stds_cfg.w_frs]
        if missing_elements:
            missing_txt = ", ".join(missing_elements)
            raise RuntimeError(f"Missing mass fraction(s) for detectable element(s): {missing_txt}.")
        nominal_w_frs = np.asarray([self.exp_stds_cfg.w_frs[el] for el in self.detectable_els_sample], dtype=float)
        zaf_pb, z_sample = zaf_calculator.get_ZAF_mult_f_pb(nominal_w_frs, peak_energies_dict)
        pb_corrected = zaf_pb * np.array(means_pb) / np.array(w_frs)
        return pb_corrected, z_sample

    def _save_std_results(
        self: Any,
        std_ref_lines: Dict[str, StandardFitLineResult],
        pb_corrected: Union[List[float], np.ndarray],
        z_sample: Any,
    ) -> Optional[StandardsFitResults]:
        if not std_ref_lines:
            return None

        means_pb = []
        stdevs_pb = []
        at_percent_values = []
        w_percent_values = []
        n_spectra_per_line = []
        corrected_pb_values = []
        rel_errors_percent = []
        fit_lines: Dict[str, StandardFitLineResult] = {}
        line_keys = list(std_ref_lines.keys())
        if len(pb_corrected) != len(line_keys):
            raise ValueError("Length of pb_corrected does not match number of reference lines.")

        std_w_frs = {el: float(w) for el, w in dict(self.exp_stds_cfg.w_frs).items()}
        std_elements = list(std_w_frs.keys())
        std_w_fr_list = [std_w_frs[el] for el in std_elements]
        std_at_fr_list = weight_to_atomic_fr(std_w_fr_list, std_elements, verbose=False)
        at_percent_by_el = {el: float(at_fr) * 100 for el, at_fr in zip(std_elements, std_at_fr_list)}
        w_percent_by_el = {el: float(w_fr) * 100 for el, w_fr in std_w_frs.items()}

        for index, el_line in enumerate(line_keys):
            line_result = std_ref_lines[el_line]
            means_pb.append(float(line_result.measured_pb))
            stdevs_pb.append(float(line_result.stdev_pb))
            el = el_line.split('_')[0]
            at_percent = at_percent_by_el.get(el, float("nan"))
            w_percent = w_percent_by_el.get(el, float("nan"))
            if not np.isnan(at_percent):
                at_percent = round(at_percent, 3)
            if not np.isnan(w_percent):
                w_percent = round(w_percent, 3)
            at_percent_values.append(at_percent)
            w_percent_values.append(w_percent)
            pb_ratios = line_result.pb_ratios
            n_spectra_used = sum((x is not None) and (not (isinstance(x, float) and np.isnan(x))) for x in pb_ratios)
            n_spectra_per_line.append(n_spectra_used)
            corrected_value = float(pb_corrected[index])
            corrected_value = round(corrected_value, 1) if not np.isnan(corrected_value) else corrected_value
            corrected_pb_values.append(corrected_value)

            mean_pb = float(line_result.measured_pb)
            mean_pb = round(mean_pb, 1) if not np.isnan(mean_pb) else mean_pb
            stdev_pb = float(line_result.stdev_pb)
            stdev_pb = round(stdev_pb, 1) if not np.isnan(stdev_pb) else stdev_pb
            rel_error = (stdev_pb / mean_pb * 100) if mean_pb else float("nan")
            rel_error = round(rel_error, 2) if not np.isnan(rel_error) else rel_error
            rel_errors_percent.append(rel_error)

            fit_lines[el_line] = StandardFitLineResult.model_validate({
                "pb_ratios": list(pb_ratios),
                "measured_pb": mean_pb,
                "stdev_pb": stdev_pb,
                "peak_theoretical_energy_keV": float(line_result.peak_theoretical_energy_keV),
                "corrected_pb": corrected_value,
                "rel_stdev_pb_percent": rel_error,
                "n_spectra_used": int(n_spectra_used),
            })

        results_df = pd.DataFrame({
            "at%": at_percent_values,
            "w%": w_percent_values,
            cnst.MEAS_PB_DF_KEY: means_pb,
            cnst.STDEV_PB_DF_KEY: stdevs_pb,
            cnst.COR_PB_DF_KEY: corrected_pb_values,
            cnst.REL_ER_PERCENT_PB_DF_KEY: rel_errors_percent,
            cnst.N_SP_USED_KEY: n_spectra_per_line,
        }, index=line_keys)

        suffix = self.output_filename_suffix if isinstance(self.output_filename_suffix, str) else ""
        filename = f"{cnst.STDS_RESULT_FILENAME}" + suffix
        results_path = os.path.join(StandardsModule._get_std_exports_dir(self), filename + cnst.DATA_FILEEXT)
        results_df.to_csv(results_path, index=True, index_label="peak", header=True)

        mean_z_payload = StandardsModule._serialize_standard_mean_z(z_sample)
        mean_z = StandardMeanZ.model_validate(mean_z_payload) if mean_z_payload is not None else None
        return StandardsFitResults.model_validate({"lines": fit_lines, "mean_z": mean_z})

    def _load_xsp_standards(self: Any) -> Tuple[EDSStandardsFile, str]:
        """Backward-compatible alias returning typed standards model."""
        return StandardsModule._load_standards(self)

    def _load_standards(self: Any) -> Tuple[EDSStandardsFile, str]:
        """Load standards as structured schema model for internal processing."""
        meas_mode = self.measurement_cfg.mode
        update_separate_std_dict = self.exp_stds_cfg.is_exp_std_measurement and self.exp_stds_cfg.generate_separate_std_dict
        project_dir: Optional[str] = None

        if self.standards is None:
            std_f_dir = None
            if update_separate_std_dict or self.quant_cfg.use_project_specific_std_dict:
                project_dir = os.path.dirname(self.sample_result_dir)
                std_f_dir = project_dir

            try:
                standards, stds_filepath = calibs.load_standards_model(
                    self.measurement_cfg.type,
                    self.measurement_cfg.beam_energy_keV,
                    std_f_dir=std_f_dir,
                )
            except FileNotFoundError:
                stds_filepath = calibs.standards_dir
                standards = EDSStandardsFile(
                    schema_version=1,
                    measurement_type=self.measurement_cfg.type,
                    beam_energy_keV=int(self.measurement_cfg.beam_energy_keV),
                    standards_by_mode={meas_mode: {}},
                )
            else:
                if update_separate_std_dict and project_dir and os.path.dirname(stds_filepath) != project_dir:
                    stds_filepath = shutil.copy(stds_filepath, project_dir)
        else:
            standards = calibs.ensure_standards_model(
                payload=self.standards,
                meas_type=self.measurement_cfg.type,
                beam_energy=self.measurement_cfg.beam_energy_keV,
            )
            stds_filepath = ''

        if meas_mode not in standards.standards_by_mode:
            standards.standards_by_mode[meas_mode] = {}

        return standards, stds_filepath

    def _update_standard_library(
        self: Any,
        fit_results: StandardsFitResults,
    ) -> None:
        from pymatgen.core.periodic_table import Element as PymatgenElement
        meas_mode = self.measurement_cfg.mode
        if self.standards is not None:
            self.standards = None
        standards, stds_filepath = StandardsModule._load_standards(self)
        std_lib = standards.standards_by_mode[meas_mode]

        for el_line, line_payload in list(std_lib.items()):
            line_payload.entries = [
                entry
                for entry in line_payload.entries
                if entry.standard_id != self.sample_id
            ]

        now = datetime.now()
        std_mean_z = fit_results.mean_z

        # Determine which elements/lines should be used for mean PB calculation
        els_to_use = self.exp_stds_cfg.els_to_use_for_mean_PB_calc
        def _should_use_for_mean(el_line):
            if els_to_use is None:
                return False
            if els_to_use == ["all"]:
                return True
            el = el_line.split("_")[0]
            return el in els_to_use


        for el_line in fit_results.lines.keys():
            line_result = fit_results.lines[el_line]
            if el_line not in std_lib:
                std_lib[el_line] = StandardLine()

            use_for_mean_calc = _should_use_for_mean(el_line)
            line_payload = std_lib[el_line]
            line_payload.entries.append(line_result.to_standard_entry(
                standard_id=self.sample_id,
                datetime=now.strftime("%Y-%m-%d %H:%M:%S"),
                formula=self.exp_stds_cfg.formula,
                std_type=self.sample_cfg.type,
                use_for_mean_calc=use_for_mean_calc,
                mean_z=std_mean_z,
            ))

            list_pb_for_mean = [
                entry.corrected_pb
                for entry in line_payload.entries
                if getattr(entry, "use_for_mean_calc", False)
            ]
            if len(list_pb_for_mean) > 0:
                mean_pb = float(np.mean(list_pb_for_mean))
                stddev_mean_pb = float(np.std(list_pb_for_mean))
                error_mean_pb = (stddev_mean_pb / mean_pb * 100) if mean_pb else float("nan")
                line_payload.reference_mean = ReferenceMean.model_validate({
                    cnst.DATETIME_KEY: now.strftime("%Y-%m-%d %H:%M:%S"),
                    cnst.COR_PB_DF_KEY: mean_pb,
                    cnst.STDEV_PB_DF_KEY: stddev_mean_pb,
                    cnst.REL_ER_PERCENT_PB_DF_KEY: error_mean_pb,
                    cnst.MEAS_PB_DF_KEY: None,
                })
            else:
                line_payload.reference_mean = None

        standards.to_json_file(stds_filepath, indent=2)
        if stds_filepath:
            logger.info(
                "✅ Standards library updated successfully for sample '%s' (%d peaks): %s",
                self.sample_id,
                len(fit_results.lines),
                stds_filepath,
            )
        else:
            logger.info(
                "✅ Standards library updated successfully for sample '%s' (%d peaks).",
                self.sample_id,
                len(fit_results.lines),
            )
