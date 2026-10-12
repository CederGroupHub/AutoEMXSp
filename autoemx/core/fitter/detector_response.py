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
import hashlib
import time
import warnings
import numpy as np
from pathlib import Path
from scipy import sparse
from scipy.integrate import quad, trapezoid
from pydantic import BaseModel

from autoemx.utils.helper import load_msa
import autoemx.utils.constants as cnst
import autoemx.calibrations as calibs
import autoemx._paths as paths
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

class SparseConvolutionRow(BaseModel):
    """Sparse representation of a matrix row."""

    start_index: int
    values: list[float]


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
    det_response_op : scipy.sparse.csr_matrix or None
        Sparse product icc_conv_matrix @ det_res_conv_matrix, applied in _apply_det_response_fncts.
    det_eff_energy_vals : np.ndarray
        Energy values for detector efficiency curve.
    det_eff_vals : np.ndarray
        Detector efficiency values.
    energy_vals_padding : int
        Padding added to energy_vals for convolution operations.
    """
    
    det_res_conv_matrix = None
    icc_conv_matrix = None
    det_response_op = None
    energy_vals_padding = 30  # Padding added to energy_vals to ensure correct functioning of convolution operation
    conv_zero_tol = 1e-14

    @classmethod
    def _conv_matrices_path(cls, det_ch_offset, det_ch_width, microscope_ID, matrix_size) -> Path:
        """
        Cache file of the convolution matrices of a detector channel setup, in ~/.autoemx/cache/<microscope_ID>/.

        The name includes a fingerprint of the detector model and of the matrix settings the matrices are computed
        from, so matrices computed with other values (e.g. after recalibrating the peak width) are never used.
        """
        if calibs.detector_model is None:
            raise AttributeError("calibrations.detector_model is not initialized")
        settings = f"{calibs.detector_model.fingerprint()},{cls.energy_vals_padding},{matrix_size},{cls.conv_zero_tol}"
        fingerprint = hashlib.sha256(settings.encode("utf-8")).hexdigest()[:12]
        cache_dir = paths.user_cache_dir(microscope_ID) / "detector_response"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir / f"conv_matrices_O{det_ch_offset}_W{det_ch_width}_{fingerprint}.npz"

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
    def _load_conv_matrices_cache(cls, file_path, matrix_size):
        """Load cached matrices (None if not cached, or unreadable)."""
        if not os.path.exists(file_path):
            return None
        try:
            with np.load(file_path) as cache:
                matrices = tuple(
                    sparse.csr_matrix(
                        (cache[f"{name}_data"], cache[f"{name}_indices"], cache[f"{name}_indptr"]),
                        shape=(matrix_size, matrix_size),
                    ).toarray()
                    for name in ("det_res", "icc")
                )
        except Exception:
            return None
        return matrices

    @classmethod
    def _compute_and_cache_conv_matrices(
        cls,
        *,
        conv_matrices_file_path,
        det_ch_offset: float,
        det_ch_width: float,
        detector_ch_n: int,
        matrix_size: int,
        verbose: bool = True,
    ):
        """
        Compute detector convolution matrices and save them to the cache.

        Leading and trailing values of each row below conv_zero_tol are set to zero, so the matrices are the same
        whether computed or loaded from the cache.
        """
        compute_start_time = time.time()
        if verbose:
            print('-' * 50, flush=True)
            print("ℹ️ No cached convolution matrices found for this detector setup.", flush=True)
            print("🔬 Calculating detector convolution matrices...", flush=True)

        full_en_vector = [det_ch_offset + j * det_ch_width for j in range(detector_ch_n)]
        matrices = tuple(
            cls._dense_from_sparse_rows([cls._sparsify_row(row, cls.conv_zero_tol) for row in matrix], matrix_size)
            for matrix in (cls._calc_det_res_conv_matrix(full_en_vector, verbose),
                           cls._calc_icc_conv_matrix(full_en_vector, verbose))
        )

        arrays = {}
        for name, matrix in zip(("det_res", "icc"), matrices):
            csr = sparse.csr_matrix(matrix)
            arrays.update({f"{name}_data": csr.data, f"{name}_indices": csr.indices, f"{name}_indptr": csr.indptr})
        # Atomic replace avoids partial reads if another process opens the cache mid-write.
        tmp_path = f"{conv_matrices_file_path}.{os.getpid()}.tmp"
        try:
            with open(tmp_path, 'wb') as file:
                np.savez_compressed(file, **arrays)
            os.replace(tmp_path, conv_matrices_file_path)
        finally:
            if os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

        if verbose:
            compute_time = time.time() - compute_start_time
            print(f"✅ Finished computing convolution matrices in {compute_time:.1f} s", flush=True)
        return matrices

    @classmethod
    def ensure_conv_matrices_cached(
        cls,
        det_ch_offset,
        det_ch_width,
        microscope_ID,
        verbose=True,
    ) -> bool:
        """
        Ensure convolution matrices for this detector channel setup exist on disk.

        Returns
        -------
        bool
            True if matrices were newly computed, False if already cached.
        """
        detector_ch_n = getattr(calibs, 'detector_ch_n', None)
        if detector_ch_n is None:
            raise AttributeError("calibrations.detector_ch_n is not initialized")

        matrix_size = detector_ch_n + cls.energy_vals_padding + cls.energy_vals_padding // 2 - 1
        conv_matrices_file_path = cls._conv_matrices_path(det_ch_offset, det_ch_width, microscope_ID, matrix_size)

        if cls._load_conv_matrices_cache(conv_matrices_file_path, matrix_size) is not None:
            return False

        cls._compute_and_cache_conv_matrices(
            conv_matrices_file_path=conv_matrices_file_path,
            det_ch_offset=det_ch_offset,
            det_ch_width=det_ch_width,
            detector_ch_n=detector_ch_n,
            matrix_size=matrix_size,
            verbose=verbose,
        )
        return True

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
        - Convolution matrices are detector-dependent and cached in ~/.autoemx/cache/<microscope_ID>/ for efficiency.
        - Prefer calling ``ensure_conv_matrices_cached`` once on the main process before
          parallel quantification so workers only hit the cached fast path.
        - If multiple EDS detectors are used, this code should be adapted.
        """
    
        # --- Load EDS detector efficiency spectrum ---
        detector_efficiency_path = str(paths.package_calibration_dir(microscope_ID) / cnst.DETECTOR_EFFICIENCY_FILENAME)
        det_eff_energy_vals, det_eff_vals, metadata = load_msa(detector_efficiency_path)
        if metadata['XUNITS'] == 'eV':
            det_eff_energy_vals /= 1000  # Convert to keV
    
        cls.det_eff_energy_vals = det_eff_energy_vals
        cls.det_eff_vals = det_eff_vals
        detector_ch_n = getattr(calibs, 'detector_ch_n', None)
        if detector_ch_n is None:
            raise AttributeError("calibrations.detector_ch_n is not initialized")
    
        # --- Load or calculate convolution matrices ---
        matrix_size = detector_ch_n + cls.energy_vals_padding + cls.energy_vals_padding // 2 - 1
        conv_matrices_file_path = cls._conv_matrices_path(det_ch_offset, det_ch_width, microscope_ID, matrix_size)

        conv_matrices = cls._load_conv_matrices_cache(conv_matrices_file_path, matrix_size)
        if conv_matrices is None:
            conv_matrices = cls._compute_and_cache_conv_matrices(
                conv_matrices_file_path=conv_matrices_file_path,
                det_ch_offset=det_ch_offset,
                det_ch_width=det_ch_width,
                detector_ch_n=detector_ch_n,
                matrix_size=matrix_size,
                verbose=verbose,
            )

        det_res_conv_matrix, icc_conv_matrix = conv_matrices
    
        # --- Crop matrices to match spectrum limits ---
        low_l = spectrum_lims[0] + 1
        high_l = spectrum_lims[1] + cls.energy_vals_padding // 2 + cls.energy_vals_padding - 1
        cls.det_res_conv_matrix = np.array(det_res_conv_matrix)[low_l:high_l, low_l:high_l]
        cls.icc_conv_matrix = np.array(icc_conv_matrix)[low_l:high_l, low_l:high_l] 
        # Combined (ICC ∘ resolution) operator, stored sparse since both matrices are banded
        cls.det_response_op = sparse.csr_matrix(cls.icc_conv_matrix @ cls.det_res_conv_matrix)

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
    
        # Convolve with the detector resolution matrix, then with the ICC matrix,
        # using the precomputed combined operator
        model = cls.det_response_op @ padded_signal
    
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
        return DetectorResponseFunction.det_sigma(E, calibs.conv_eff, calibs.elec_noise, calibs.F)

    @staticmethod
    def det_sigma(E, conv_eff, elec_noise, F):
        """
        Detector Gaussian sigma (keV) at X-ray energy E (keV), for given detector parameters.

        sigma = conv_eff * sqrt(elec_noise^2 + E * F / conv_eff), with conv_eff the energy per
        electron-hole pair (keV), elec_noise the electronic noise (in electrons) and F the Fano factor.
        See _det_sigma for the version using the loaded microscope calibrations.
        """
        return conv_eff * np.sqrt(elec_noise**2 + E * F / conv_eff)
        
    
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
    
        def dQ_dz(z):
            return -np.pi * F_loss / V_tot * (z ** 2 - R_e ** 2)
    
        def dN_dz(z):
            return alpha * np.exp(-alpha * z)
    
        # Calculate left boundary of ICC smearing function. Q(0) = 1 - F_loss / 2, since V_1(0) = V_tot / 2
        Q_min = 1 - F_loss / 2
        E_min = line_en * Q_min
    
        e_vals = np.linspace(E_min, line_en, 1000)[:-1]  # Remove last energy value corresponding to line_en
    
        # Invert Q(z) = E / line_en analytically, with Q(z) = (V_tot - F_loss * V_1(z)) / V_tot and
        # V_1(z) = pi / 3 * (R_e - z)^2 * (2 * R_e + z). With u = z / R_e, Q(z) = q reduces to the cubic
        # u^3 - 3u + (2 - c) = 0, with c = 4 * (1 - q) / F_loss, whose root in [0, 1] is given below
        Q_vals = np.clip(e_vals / line_en, Q_min, 1)
        c = 4 * (1 - Q_vals) / F_loss
        u = 2 * np.cos(np.arccos(np.clip((c - 2) / 2, -1, 1)) / 3 - 2 * np.pi / 3)
        z_vals = u * R_e
        n_vals = dN_dz(z_vals) * dQ_dz(z_vals) ** -1 / line_en
    
        # Value at line_en such that the trapezoid over the last interval brings the total integral to 1
        signal_integral = trapezoid(n_vals, e_vals)
        n_val_at_E = 2 * (1 - signal_integral) / (line_en - e_vals[-1]) - n_vals[-1]
    
        # Update lists with values at line_en
        e_vals = list(e_vals) + [line_en]
        n_vals = list(n_vals) + [n_val_at_E]
    
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
