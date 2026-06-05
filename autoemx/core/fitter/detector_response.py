#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Detector Response Function module.

Provides handling of EDS detector response, including convolution matrices for energy resolution
and incomplete charge collection (ICC).

This module is used by both peaks and background models to accurately simulate detector effects.

Main Class:
    - DetectorResponseFunction: Manages detector response and convolution operations
"""

import os
import json
import time
import warnings
import numpy as np
from pathlib import Path
from scipy.integrate import quad, trapezoid
from scipy.optimize import root_scalar
from pydantic import BaseModel

from autoemx.utils.helper import print_single_separator, load_msa
import autoemx.utils.constants as cnst
import autoemx.calibrations as calibs
from autoemx.data.Xray_absorption_coeffs import xray_mass_absorption_coeff

from autoemx._logging import get_logger
logger = get_logger(__name__)


def _safe_log(level: str, message: str, *, verbose: bool = True) -> None:
    """Best-effort logging that never interrupts numerical routines."""
    if not verbose:
        return
    try:
        getattr(logger, level)(message)
    except BaseException:
        # On some microscope PCs, stream writes can fail and bubble up as
        # KeyboardInterrupt/lost stderr from worker contexts.
        try:
            print(message, flush=True)
        except Exception:
            pass

parent_dir = str(Path(__file__).resolve().parent.parent.parent)


class SparseConvolutionRow(BaseModel):
    """Sparse representation of a matrix row."""

    start_index: int
    values: list[float]


class DetectorConvolutionMatricesCache(BaseModel):
    """Typed cache payload for detector convolution matrices."""

    det_ch_offset: float
    det_ch_width: float
    detector_ch_n: int
    energy_vals_padding: int
    matrix_size: int
    zero_tol: float
    det_res_rows: list[SparseConvolutionRow]
    icc_rows: list[SparseConvolutionRow]


class DetectorResponseFunction:
    """
    Provides static and class methods for handling the detector's instrumental response.
    
    Includes convolution matrices for energy resolution and incomplete charge collection (ICC).
    This class is initialized with calibration data and is used by both background and peak 
    models to accurately simulate detector effects.
    
    Class Attributes
    ----------------
    det_res_conv_matrix : np.ndarray or None
        Detector resolution convolution matrix (cropped to spectrum limits).
    icc_conv_matrix : np.ndarray or None
        ICC convolution matrix (cropped to spectrum limits).
    det_eff_energy_vals : np.ndarray
        Energy values for detector efficiency curve.
    det_eff_vals : np.ndarray
        Detector efficiency values.
    energy_vals_padding : int
        Padding added to energy_vals for convolution operations.
    """
    
    det_res_conv_matrix = None
    icc_conv_matrix = None
    energy_vals_padding = 30  # Padding added to energy_vals to ensure correct functioning of convolution operation
    conv_zero_tol = 1e-14

    @staticmethod
    def _conv_matrices_key(det_ch_offset, det_ch_width):
        return f"O{det_ch_offset},W{det_ch_width}"

    @staticmethod
    def _sparsify_row(row_vals, zero_tol):
        """Trim leading and trailing zero values from a dense row."""
        row_arr = np.asarray(row_vals, dtype=float)
        nz_idx = np.where(np.abs(row_arr) > zero_tol)[0]
        if nz_idx.size == 0:
            return SparseConvolutionRow(start_index=0, values=[])
        start_index = int(nz_idx[0])
        end_index = int(nz_idx[-1]) + 1
        return SparseConvolutionRow(start_index=start_index, values=row_arr[start_index:end_index].tolist())

    @staticmethod
    def _dense_from_sparse_rows(sparse_rows, matrix_size):
        """Rebuild a dense matrix from sparse rows, padding missing values with zeros."""
        dense = np.zeros((matrix_size, matrix_size), dtype=float)
        for row_i, row in enumerate(sparse_rows[:matrix_size]):
            if not row.values:
                continue
            col_start = max(0, int(row.start_index))
            if col_start >= matrix_size:
                continue
            row_vals = np.asarray(row.values, dtype=float)
            col_end = min(matrix_size, col_start + row_vals.size)
            dense[row_i, col_start:col_end] = row_vals[:col_end - col_start]
        return dense

    @classmethod
    def _load_conv_matrices_cache(cls, file_path, det_ch_offset, det_ch_width, matrix_size):
        """Load typed cache matrices padded to matrix_size."""
        if not os.path.exists(file_path):
            return None
        detector_ch_n = getattr(calibs, 'detector_ch_n', None)
        if detector_ch_n is None:
            return None
        try:
            with open(file_path, 'r') as file:
                payload = json.load(file)
        except Exception:
            return None

        cache_key = cls._conv_matrices_key(det_ch_offset, det_ch_width)
        if not isinstance(payload, dict) or cache_key not in payload:
            return None

        try:
            cache_payload = DetectorConvolutionMatricesCache.model_validate(payload[cache_key])
        except Exception:
            return None

        same_settings = (
            cache_payload.detector_ch_n == detector_ch_n and
            cache_payload.energy_vals_padding == cls.energy_vals_padding and
            np.isclose(cache_payload.det_ch_offset, det_ch_offset) and
            np.isclose(cache_payload.det_ch_width, det_ch_width)
        )
        if not same_settings:
            return None

        det_res_dense = cls._dense_from_sparse_rows(cache_payload.det_res_rows, matrix_size)
        icc_dense = cls._dense_from_sparse_rows(cache_payload.icc_rows, matrix_size)
        if det_res_dense.shape != (matrix_size, matrix_size) or icc_dense.shape != (matrix_size, matrix_size):
            return None
        return det_res_dense, icc_dense

    @classmethod
    def setup_detector_response_vars(cls, det_ch_offset, det_ch_width, spectrum_lims, microscope_ID, verbose=True):
        """
        Initialize detector response variables for the EDS system.
    
        Loads detector efficiency and convolution matrices for the specified microscope.
        If convolution matrices for the given channel settings do not exist, they are calculated and saved.
    
        Parameters
        ----------
        det_ch_offset : float
            Energy offset for detector channels (in keV).
        det_ch_width : float
            Channel width (in keV).
        spectrum_lims : tuple of int
            (low, high) indices for the usable spectrum region.
        microscope_ID : str
            Identifier for the microscope/calibration directory.
        verbose : bool, optional
            If True, print status messages.
    
        Sets Class Attributes
        ---------------------
        det_eff_energy_vals : np.ndarray
            Energy values for detector efficiency.
        det_eff_vals : np.ndarray
            Detector efficiency values.
        det_res_conv_matrix : np.ndarray
            Detector response convolution matrix (cropped to spectrum_lims).
        icc_conv_matrix : np.ndarray
            ICC convolution matrix (cropped to spectrum_lims).
    
        Notes
        -----
        - Convolution matrices are detector-dependent and cached in JSON for efficiency.
        - If multiple EDS detectors are used, this code should be adapted.
        """
    
        # --- Load EDS detector efficiency spectrum ---
        detector_efficiency_path = os.path.join(
            parent_dir, cnst.XRAY_SPECTRA_CALIBS_DIR, microscope_ID, cnst.DETECTOR_EFFICIENCY_FILENAME
        )
        det_eff_energy_vals, det_eff_vals, metadata = load_msa(detector_efficiency_path)
        if metadata['XUNITS'] == 'eV':
            det_eff_energy_vals /= 1000  # Convert to keV
    
        cls.det_eff_energy_vals = det_eff_energy_vals
        cls.det_eff_vals = det_eff_vals
        detector_ch_n = getattr(calibs, 'detector_ch_n', None)
        if detector_ch_n is None:
            raise AttributeError("calibrations.detector_ch_n is not initialized")
    
        # --- Load or calculate convolution matrices ---
        conv_matrices_file_path = os.path.join(
            parent_dir, cnst.XRAY_SPECTRA_CALIBS_DIR, microscope_ID, cnst.DETECTOR_CONV_MATRICES_FILENAME
        )
        lock_file_path = conv_matrices_file_path + ".lock"
        matrix_size = detector_ch_n + cls.energy_vals_padding + cls.energy_vals_padding // 2 - 1
        
        conv_matrices = None

        # 1. FAST PATH: Load cache without lock if settings match.
        conv_matrices = cls._load_conv_matrices_cache(
            conv_matrices_file_path,
            det_ch_offset,
            det_ch_width,
            matrix_size,
        )

        # 2. SLOW PATH: We need to compute it (or wait for another core to compute it)
        if conv_matrices is None:
            start_wait_time = time.time()
            max_wait_time = 600  # 10 minute timeout for stale locks
            stale_lock_max_age_s = 120.0

            def _read_lock_pid(path):
                try:
                    with open(path, 'r') as f:
                        txt = f.read().strip()
                    if txt.startswith('PID '):
                        return int(txt.split(' ', 1)[1])
                except Exception:
                    return None
                return None

            def _is_process_alive(pid):
                if pid is None or pid <= 0:
                    return False
                try:
                    os.kill(pid, 0)
                except OSError:
                    return False
                return True

            def _compute_conv_matrices_without_lock():
                compute_start_time = time.time()
                if verbose:
                    print('-' * 50, flush=True)
                    print("ℹ️ No cached convolution matrices found for this detector setup.", flush=True)
                    print("🔬 Calculating detector convolution matrices...", flush=True)

                full_en_vector = [det_ch_offset + j * det_ch_width for j in range(detector_ch_n)]
                det_res_conv_matrix = cls._calc_det_res_conv_matrix(full_en_vector, verbose)
                icc_conv_matrix = cls._calc_icc_conv_matrix(full_en_vector, verbose)

                cache_payload = DetectorConvolutionMatricesCache(
                    det_ch_offset=det_ch_offset,
                    det_ch_width=det_ch_width,
                    detector_ch_n=detector_ch_n,
                    energy_vals_padding=cls.energy_vals_padding,
                    matrix_size=matrix_size,
                    zero_tol=cls.conv_zero_tol,
                    det_res_rows=[cls._sparsify_row(row, cls.conv_zero_tol) for row in det_res_conv_matrix],
                    icc_rows=[cls._sparsify_row(row, cls.conv_zero_tol) for row in icc_conv_matrix],
                )

                cache_key = cls._conv_matrices_key(det_ch_offset, det_ch_width)
                cache_file_payload = {}
                if os.path.exists(conv_matrices_file_path):
                    try:
                        with open(conv_matrices_file_path, 'r') as file:
                            existing_payload = json.load(file)
                        if isinstance(existing_payload, dict):
                            cache_file_payload = existing_payload
                    except Exception:
                        cache_file_payload = {}

                cache_file_payload[cache_key] = cache_payload.model_dump()
                with open(conv_matrices_file_path, 'w') as file:
                    json.dump(cache_file_payload, file)

                if verbose:
                    compute_time = time.time() - compute_start_time
                    print(f"✅ Finished computing convolution matrices in {compute_time:.1f} s", flush=True)
                return det_res_conv_matrix, icc_conv_matrix

            # Loop until conv_matrices is successfully populated
            while conv_matrices is None:
                try:
                    try:
                        # Attempt to exclusively create the lock file.
                        with open(lock_file_path, 'x') as f:
                            f.write(f"PID {os.getpid()}")

                        # === WE HAVE THE LOCK ===
                        try:
                            # Reload cache: another process might have finished while we waited.
                            conv_matrices = cls._load_conv_matrices_cache(
                                conv_matrices_file_path,
                                det_ch_offset,
                                det_ch_width,
                                matrix_size,
                            )

                            # If it's STILL missing, we actually do the heavy lifting.
                            if conv_matrices is None:
                                conv_matrices = _compute_conv_matrices_without_lock()

                        finally:
                            # === RELEASE THE LOCK ===
                            if os.path.exists(lock_file_path):
                                os.remove(lock_file_path)

                    except FileExistsError:
                        # === ANOTHER CORE HAS THE LOCK ===
                        lock_age = None
                        if os.path.exists(lock_file_path):
                            try:
                                lock_age = time.time() - os.path.getmtime(lock_file_path)
                            except OSError:
                                lock_age = None
                        lock_pid = _read_lock_pid(lock_file_path) if os.path.exists(lock_file_path) else None
                        lock_owner_dead = lock_pid is not None and not _is_process_alive(lock_pid)
                        lock_too_old = lock_age is not None and lock_age > stale_lock_max_age_s

                        if lock_owner_dead or lock_too_old:
                            try:
                                os.remove(lock_file_path)
                                start_wait_time = time.time()
                                if verbose:
                                    reason = []
                                    if lock_owner_dead:
                                        reason.append(f"owner PID {lock_pid} is not alive")
                                    if lock_too_old:
                                        reason.append(f"lock age {lock_age:.1f} s > {stale_lock_max_age_s:.0f} s")
                                    print_single_separator()
                                    print(f"⚠️ Removed stale convolution lock ({'; '.join(reason)}).", flush=True)
                                continue
                            except OSError:
                                pass

                        if time.time() - start_wait_time > max_wait_time:
                            # Lock is stale (the other core crashed). Force delete it and compute locally.
                            try:
                                os.remove(lock_file_path)
                                start_wait_time = time.time()
                                if verbose:
                                    print_single_separator()
                                    print("⚠️ Removed stale convolution lock file; falling back to local computation.", flush=True)
                            except OSError:
                                pass
                            conv_matrices = _compute_conv_matrices_without_lock()
                        else:
                            # Wait patiently.
                            time.sleep(3)

                            # Peek at the file to see if the other core finished our key.
                            conv_matrices = cls._load_conv_matrices_cache(
                                conv_matrices_file_path,
                                det_ch_offset,
                                det_ch_width,
                                matrix_size,
                            )

                except KeyboardInterrupt:
                    # Some environments surface interrupted child work as KeyboardInterrupt.
                    # Fall back to a local, non-locking computation so the run can continue.
                    if verbose:
                        print_single_separator()
                        print("⚠️ Convolution cache setup was interrupted; falling back to local computation.", flush=True)
                    try:
                        if os.path.exists(lock_file_path):
                            os.remove(lock_file_path)
                    except OSError:
                        pass
                    conv_matrices = _compute_conv_matrices_without_lock()
                
        det_res_conv_matrix, icc_conv_matrix = conv_matrices
    
        # --- Crop matrices to match spectrum limits ---
        low_l = spectrum_lims[0] + 1
        high_l = spectrum_lims[1] + cls.energy_vals_padding // 2 + cls.energy_vals_padding - 1
        cls.det_res_conv_matrix = np.array(det_res_conv_matrix)[low_l:high_l, low_l:high_l]
        cls.icc_conv_matrix = np.array(icc_conv_matrix)[low_l:high_l, low_l:high_l] 

    # =============================================================================
    # Convolution of signal with detector response function
    # =============================================================================
    @classmethod
    def _apply_padding_with_fit(cls, signal):
        """
        Pad a signal at both ends by linear extrapolation, using a fit to the first and last few points.
    
        This method is used to extend the signal array, avoiding edge effects in convolution.
        Padding values are clipped to be non-negative.
    
        Parameters
        ----------
        signal : np.ndarray or list
            The input 1D signal to pad.
    
        Returns
        -------
        padded_signal : np.ndarray
            The input signal with linear-extrapolated padding added at both ends.
    
        Notes
        -----
        - The number of padding points is determined by cls.energy_vals_padding.
        - Linear fits are performed on the first and last 4 points of the signal.
        - Extrapolated values are clipped at zero to avoid negative padding.
        """
        n_pts_fitted = 4  # Number of points used for linear fit at each end
    
        # --- Padding at the beginning (head) ---
        x_head = signal[:n_pts_fitted]
        x_indices_head = np.arange(len(x_head))
        # Linear fit to the first few points
        slope_head, intercept_head = np.polyfit(x_indices_head, x_head, 1)
        # Indices for extrapolation (negative indices for padding before signal)
        extrapolation_indices_head = np.arange(-cls.energy_vals_padding // 2 + 1, 0)
        # Linear extrapolation for padding
        extrapolated_values_head = slope_head * extrapolation_indices_head + intercept_head
        # Ensure no negative values in padding
        extrapolated_values_head = np.clip(extrapolated_values_head, 0, None)
    
        # --- Padding at the end (tail) ---
        x_tail = signal[-n_pts_fitted:]
        x_indices_tail = np.arange(len(x_tail))
        # Linear fit to the last few points
        slope_tail, intercept_tail = np.polyfit(x_indices_tail, x_tail, 1)
        # Indices for extrapolation (beyond signal end)
        extrapolation_indices_tail = np.arange(len(x_tail) + 1, len(x_tail) + cls.energy_vals_padding)
        # Linear extrapolation for padding
        extrapolated_values_tail = slope_tail * extrapolation_indices_tail + intercept_tail
        # Ensure no negative values in padding
        extrapolated_values_tail = np.clip(extrapolated_values_tail, 0, None)
    
        # --- Combine padding and original signal ---
        padded_signal = np.concatenate([extrapolated_values_head, signal, extrapolated_values_tail])
    
        return padded_signal
    
    
    @classmethod
    def _apply_det_response_fncts(cls, signal):
        """
        Apply both detector resolution and ICC convolution functions to a signal,
        with edge padding by linear extrapolation.
    
        The signal is padded at both ends using a linear fit, then convolved sequentially
        with the detector resolution and ICC convolution matrices. The padding is removed
        from the final result.
    
        Parameters
        ----------
        signal : np.ndarray or list
            The input 1D signal to be convolved.
    
        Returns
        -------
        processed_model : np.ndarray
            The signal after both convolutions, trimmed to exclude the padding.
    
        Notes
        -----
        - Padding is performed by fitting and extrapolating the first and last few points.
        - The convolution matrices must be initialized in the class before use.
        - This approach is more accurate than simple replicative padding, but slightly slower.
        """
        # Pad the signal at both ends using linear fit/extrapolation
        padded_signal = cls._apply_padding_with_fit(signal)
    
        # First, convolve with the detector resolution matrix
        model = np.sum(cls.det_res_conv_matrix * padded_signal, axis=1)
    
        # Then, convolve with the ICC convolution matrix
        model = np.sum(cls.icc_conv_matrix * model, axis=1)
    
        # Remove padding from the result to return only the original signal region
        processed_model = model[cls.energy_vals_padding // 2 - 1 : -cls.energy_vals_padding + 1]
    
        return processed_model
    
    # =============================================================================
    # Detector resolution convolution
    # =============================================================================
    @staticmethod
    def _det_sigma(E):
        """
        Calculate the detector Gaussian sigma for a given X-ray energy.
        
        Requires calibration file to be correctly loaded via setup_detector_response_vars.
    
        Based on:
        N.W.M. Ritchie, "Spectrum Simulation in DTSA-II", Microsc. Microanal. 15 (2009) 454–468.
        https://doi.org/10.1017/S1431927609990407
    
        Parameters
        ----------
        E : float or array-like
            X-ray energy in keV.
    
        Returns
        -------
        sigma : float or np.ndarray
            Standard deviation (sigma) of the detector response at energy E.
    
        Notes
        -----
        - Parameters (conv_eff, elec_noise, F) are calibrated on several elements in bulk standards.
        - See Calculation_pars_peak_fwhm.py for details.
        """
        # Calculate detector sigma using calibrated parameters
        sigma = calibs.conv_eff * np.sqrt(calibs.elec_noise**2 + E * calibs.F / calibs.conv_eff)
        return sigma
        
    
    @classmethod
    def _calc_det_res_conv_matrix(cls, energy_vals, verbose = True):
        """
        Calculate the detector resolution convolution matrix.
    
        Each row of the matrix represents the probability distribution (Gaussian)
        for an input energy, accounting for the detector's energy resolution.
        Padding is added to minimize edge effects during convolution.
    
        Parameters
        ----------
        energy_vals : array-like
            Array of energy values (in keV) for which to compute the convolution matrix.
        verbose : bool, optional
            If True, print status messages.
    
        Returns
        -------
        det_res_conv_matrix : np.ndarray
            The detector resolution convolution matrix (energy_vals x energy_vals).
    
        Notes
        -----
        - Padding is applied to account for convolution spillover at the spectrum edges.
        - The Gaussian sigma is energy-dependent and calculated with _det_sigma().
        - Integration is performed for each matrix element to ensure normalization.
        """
    
        start_time = time.time()
        if verbose:
            _safe_log("info", "🔬 Calculating convolution matrix for detector resolution")
    
        deltaE = energy_vals[5] - energy_vals[4]
        n_intervals = cls.energy_vals_padding
    
        # Extend the energy axis with padding on both sides to avoid edge effects
        left_pad = [energy_vals[0] - deltaE * i for i in range(n_intervals // 2, 0, -1)]
        right_pad = [energy_vals[-1] + deltaE * i for i in range(1, n_intervals)]
        energy_vals_extended = left_pad + list(energy_vals) + right_pad
    
        def gaussian(E, E0, sigma):
            """Normalized Gaussian function."""
            return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 / sigma**2 * (E - E0)**2)

        len_axis = len(energy_vals_extended)
        sparse_rows = []
        prev_start_index = 0
        for i, en in enumerate(energy_vals_extended):
            sigma = cls._det_sigma(en)
            row_start_index = None
            row_values = []
            seen_nonzero = False
            zero_streak = 0

            for idx in range(prev_start_index, len_axis):
                cen_E = energy_vals_extended[idx]
                try:
                    int_E, _ = quad(
                        lambda E: gaussian(E, en, sigma),
                        cen_E - deltaE / 2,
                        cen_E + deltaE / 2,
                    )
                except Exception:
                    int_E = 0.0

                if abs(int_E) > cls.conv_zero_tol:
                    if row_start_index is None:
                        row_start_index = idx
                    row_values.append(int_E)
                    seen_nonzero = True
                    zero_streak = 0
                elif seen_nonzero:
                    zero_streak += 1
                    if zero_streak >= 2:
                        break

            if row_start_index is None:
                sparse_rows.append(SparseConvolutionRow(start_index=0, values=[]))
            else:
                sparse_rows.append(
                    SparseConvolutionRow(start_index=row_start_index, values=row_values)
                )
                prev_start_index = row_start_index

        det_res_conv_matrix = cls._dense_from_sparse_rows(sparse_rows, len_axis).T
        
        if verbose:
            process_time = time.time() - start_time
            _safe_log("info", f"✅ Calculation executed in {process_time:.1f} s")
    
        return det_res_conv_matrix
    
    
    # =============================================================================
    # Incomplete Charge Collection Spectrum Calculations
    # =============================================================================
    @staticmethod
    def get_icc_spectrum(energy_vals, line_en, R_e=50e-7, F_loss=0.27):
        """
        Generate the ICC (Incomplete Charge Collection) smearing function for a given line energy.
    
        Parameters
        ----------
        energy_vals : array-like
            Energy axis (in keV) of detector channels.
        line_en : float
            X-ray line energy (in keV).
        R_e : float, optional
            Effective recombination parameter (cm).
        F_loss : float, optional
            Fractional charge loss parameter.
    
        Returns
        -------
        icc_n_vals_distr : list
            ICC distribution, mapped to detector channel energies.
        """
        icc_e_vals, icc_n_vals = DetectorResponseFunction._icc_fnct(line_en, R_e, F_loss)
        icc_n_vals_distr = DetectorResponseFunction._distribute_icc_over_EDS_channels(
            energy_vals, icc_e_vals, icc_n_vals
        )
        return icc_n_vals_distr
    
    
    @staticmethod
    def _icc_fnct(line_en, R_e=50e-7, F_loss=0.27):
        """
        Calculate the ICC function n(E) for a given X-ray line energy.
        
        ICC model as described in:
        Redus, R. H., & Huber, A. C. (2015). Response Function of Silicon Drift Detectors for Low Energy X-rays.
        In Advances in X-ray Analysis (AXA) (pp. 274–282). International Centre for Diffraction Data (ICDD).
    
        Parameters
        ----------
        line_en : float
            X-ray line energy (in keV).
        R_e : float, optional
            Effective recombination parameter (cm).
        F_loss : float, optional
            Fractional charge loss parameter.
    
        Returns
        -------
        e_vals : list of float
            Energy values (in keV) for the ICC function.
        n_vals : list of float
            ICC function values at those energies.
        """
        # Absorption coefficient of Si at energy line_en
        alpha = xray_mass_absorption_coeff(element='Si', energies=line_en) * calibs.Si_density  # cm^-1
    
        V_tot = 4 * np.pi / 3 * (R_e) ** 3
    
        def V_1(z):
            return np.pi / 3 * (R_e - z) ** 2 * (2 * R_e + z)
    
        def Q(z):
            return (V_tot - F_loss * V_1(z)) / V_tot
    
        def dQ_dz(z):
            return -np.pi * F_loss / V_tot * (z ** 2 - R_e ** 2)
    
        def N(z):
            return 1 - np.exp(-alpha * z)
    
        def dN_dz(z):
            return alpha * np.exp(-alpha * z)
    
        def get_z(Q_val, z_min=0.0):
            Q_val_rnd = np.clip(Q_val, Q_min, 1)
            solution = root_scalar(lambda z: Q(z) - Q_val_rnd, method='brentq', bracket=[z_min, R_e])
            return solution.root
    
        def n(x, z_min=0.0):
            Q_val = x / line_en
            z_val = get_z(Q_val, z_min=z_min)
            n_val = dN_dz(z_val) * dQ_dz(z_val) ** -1 / line_en
            return n_val, z_val
    
        def get_n_at_line_en(integral_rest, last_E, last_n):
            def n_fnct(n_):
                n_ = np.float64(n_)
                res = np.trapz([last_n, n_], [last_E, line_en])
                return res
            guess = (1 - integral_rest) / (line_en - last_E)
            guess_2 = guess * 2
            solution = root_scalar(lambda n_: n_fnct(n_) - (1 - integral_rest), x0=guess, x1=guess_2, method='secant')
            return solution.root
    
        # Calculate left boundary of ICC smearing function
        Q_min = Q(0)
        E_min = line_en * Q_min
    
        e_vals = list(np.linspace(E_min, line_en, 1000))
        e_vals.pop()  # Remove last energy value corresponding to line_en
        n_vals = []
        prev_z = 0.0
        for en in e_vals:
            n_val, prev_z = n(en, z_min=prev_z)
            n_vals.append(n_val)
        signal_integral = trapezoid(n_vals, e_vals)
        n_val_at_E = get_n_at_line_en(signal_integral, e_vals[-1], n_vals[-1])
    
        # Update lists with values at line_en
        e_vals.append(line_en)
        n_vals.append(n_val_at_E)
    
        return e_vals, n_vals


    @staticmethod
    def _distribute_icc_over_EDS_channels(eds_en_vals, icc_en_vals, icc_n_vals):
        """
        Distribute the ICC function over EDS detector channels.
    
        Parameters
        ----------
        eds_en_vals : array-like
            Detector channel energy values (in keV).
        icc_en_vals : list
            ICC function energy values (in keV).
        icc_n_vals : list
            ICC function values.
    
        Returns
        -------
        eds_icc_n_vals : list
            ICC values distributed over the detector channels.
        """
        ch_width = eds_en_vals[1] - eds_en_vals[0]
        icc_en_spacing = icc_en_vals[1] - icc_en_vals[0]
        icc_en_vals = np.asarray(icc_en_vals, dtype=float)
        icc_n_vals = np.asarray(icc_n_vals, dtype=float)
    
        # Determine which channels are affected by ICC
        left_bound = icc_en_vals[0] - ch_width / 2
        right_bound = icc_en_vals[-1] + ch_width / 2
        first_affected = int(np.searchsorted(eds_en_vals, left_bound, side='right'))
        last_affected = int(np.searchsorted(eds_en_vals, right_bound, side='right'))
        indices_affected = list(range(first_affected, last_affected))
        # Calculate number of points to add on the right side of the list to make it symmetrical. Needed to avoid shifts during convolution
        n_pts_to_center_data = len(indices_affected) - 1
        n_pts_added = 20  # Pad array of energy values for full overlap during convolution
        
        # Check if reference energy is outside the range of energy values. This can happen with characteristic X-rays outside the energy range
        if len(indices_affected) == 0:
            return [1] # ICC convolution is not applied if peak is outside energy range
    
        indices_affected = (
            list(range(indices_affected[0] - n_pts_added, indices_affected[0])) +
            indices_affected +
            list(range(indices_affected[-1] + 1, indices_affected[-1] + n_pts_added + n_pts_to_center_data + 1))
        )
    
        # Remove negative indices, maintaining symmetry
        if indices_affected[0] < 0:
            index_zero = indices_affected.index(0)
            indices_affected = indices_affected[index_zero: len(indices_affected) - index_zero]
    
        # Remove indices beyond dimension of eds_en_vals if needed
        len_eds_en_vals = len(eds_en_vals)
        if indices_affected[-1] >= len_eds_en_vals:
            index_last = indices_affected.index(len_eds_en_vals)
            n_pts_to_remove = len(indices_affected) - index_last + 1
            indices_affected = indices_affected[n_pts_to_remove: len(indices_affected) - n_pts_to_remove]
    
        # Distribute ICC function over the detector channels
        eds_icc_en_vals = [eds_en_vals[i] for i in indices_affected]
        eds_icc_n_vals = []
        for index, en in enumerate(eds_icc_en_vals):
            interval_boundary_left = en - ch_width / 2
            interval_boundary_right = en + ch_width / 2
            left_idx = int(np.searchsorted(icc_en_vals, interval_boundary_left, side='right'))
            right_idx = int(np.searchsorted(icc_en_vals, interval_boundary_right, side='right'))
            if left_idx < right_idx:
                e_vals_to_int = icc_en_vals[left_idx:right_idx]
                n_vals_to_int = icc_n_vals[left_idx:right_idx]
                if len(e_vals_to_int) > 1:
                    # Integrate ICC function over interval corresponding to energy value en
                    eds_icc_n_val = trapezoid(n_vals_to_int, e_vals_to_int)
                else: # Case of only 1 point within the detector channel
                    # There is no full interval of en_spacing width within the current detector channel
                    # The portion of interval within this channel is added on the next steps
                    eds_icc_n_val = 0
                # Add portion of interval shared with left of en, unless at boundary
                if interval_boundary_left < e_vals_to_int[0] and e_vals_to_int[0] > 0 and left_idx != 0:
                    extra_i_left = left_idx - 1
                    left_int = trapezoid([icc_n_vals[extra_i_left], n_vals_to_int[0]],
                                         [icc_en_vals[extra_i_left], e_vals_to_int[0]])
                    eds_icc_n_val += left_int * (e_vals_to_int[0] - interval_boundary_left) / icc_en_spacing
                # Add portion of interval shared with right of en, unless at boundary
                if interval_boundary_right > e_vals_to_int[-1] and right_idx < len(icc_n_vals):
                    extra_i_right = right_idx
                    right_int = trapezoid([n_vals_to_int[-1], icc_n_vals[extra_i_right]],
                                          [e_vals_to_int[-1], icc_en_vals[extra_i_right]])
                    eds_icc_n_val += right_int * (interval_boundary_right - e_vals_to_int[-1]) / icc_en_spacing
            else:
                eds_icc_n_val = 0
            eds_icc_n_vals.append(eds_icc_n_val)
    
        return eds_icc_n_vals


    @classmethod
    def _calc_icc_conv_matrix(cls, energy_vals, verbose = True):
        """
        Calculate the ICC convolution matrix for all detector channels.
    
        Parameters
        ----------
        energy_vals : array-like
            Array of energy values (in keV) for which to compute the convolution matrix.
        verbose : bool, optional
            If True, print status messages.
            
        Returns
        -------
        icc_conv_matrix : np.ndarray
            The ICC convolution matrix (energy_vals x energy_vals).
        """
        start_time = time.time()
        if verbose:
            _safe_log("info", "🔬 Calculating convolution matrix for incomplete charge collection")
    
        deltaE = energy_vals[5] - energy_vals[4]
        n_intervals = cls.energy_vals_padding
    
        # Extend the energy axis with padding on both sides to avoid edge effects
        left_pad = [energy_vals[0] - deltaE * i for i in range(n_intervals // 2, 0, -1)]
        right_pad = [energy_vals[-1] + deltaE * i for i in range(1, n_intervals)]
        energy_vals_extended = left_pad + list(energy_vals) + right_pad
    
        len_row = len(energy_vals_extended)
        sparse_rows = []
        prev_start_index = 0
        seen_nonzero = False
        trailing_zero_rows = 0

        for i, en in enumerate(energy_vals_extended):
            if verbose:
                _safe_log("debug", f'  {i}\tEnergy: {en * 1000:.1f} eV')

            if en > 0:
                icc_spec = DetectorResponseFunction.get_icc_spectrum(
                    energy_vals_extended, en, calibs.R_e_background, calibs.F_loss_background
                )
                if len(icc_spec) == 0:
                    icc_spec = [0]

                icc_n_vals = np.zeros([len_row])
                icc_n_vals[i] = 1
                icc_n_vals = np.convolve(icc_n_vals, icc_spec, mode='same')

                row_start_index = None
                row_values = []
                seen_row_nonzero = False
                zero_streak = 0
                for idx in range(prev_start_index, len_row):
                    val = float(icc_n_vals[idx])
                    if abs(val) > cls.conv_zero_tol:
                        if row_start_index is None:
                            row_start_index = idx
                        row_values.append(val)
                        seen_row_nonzero = True
                        zero_streak = 0
                    elif seen_row_nonzero:
                        zero_streak += 1
                        if zero_streak >= 2:
                            break

                if row_start_index is None:
                    sparse_rows.append(SparseConvolutionRow(start_index=0, values=[]))
                    if seen_nonzero:
                        trailing_zero_rows += 1
                        if trailing_zero_rows >= 2:
                            sparse_rows.extend(
                                [SparseConvolutionRow(start_index=0, values=[]) for _ in range(len_row - i - 1)]
                            )
                            break
                else:
                    sparse_rows.append(
                        SparseConvolutionRow(start_index=row_start_index, values=row_values)
                    )
                    prev_start_index = row_start_index
                    seen_nonzero = True
                    trailing_zero_rows = 0
            else:
                sparse_rows.append(SparseConvolutionRow(start_index=0, values=[]))
                if seen_nonzero:
                    trailing_zero_rows += 1
                    if trailing_zero_rows >= 2:
                        sparse_rows.extend(
                            [SparseConvolutionRow(start_index=0, values=[]) for _ in range(len_row - i - 1)]
                        )
                        break

        icc_conv_matrix = cls._dense_from_sparse_rows(sparse_rows, len_row).T
        
        if verbose:
            process_time = time.time() - start_time
            _safe_log("info", f"✅ Calculation executed in {process_time:.1f} s")
    
        return icc_conv_matrix
