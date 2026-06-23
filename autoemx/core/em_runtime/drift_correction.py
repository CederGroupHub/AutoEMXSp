import numpy as np
import cv2
from skimage.registration import phase_cross_correlation
from skimage.filters import window
from typing import Tuple


def preprocess_image(img: np.ndarray) -> np.ndarray:
    """
    Normalize and apply a window function to reduce FFT edge effects.
    Critical for non-periodic real-world images.
    """
    img_float = img.astype(np.float32)

    denom = img_float.max() - img_float.min()
    if denom > 1e-6:
        img_float = (img_float - img_float.min()) / denom

    w = window('hann', img_float.shape)
    return img_float * w


def shift_image(img: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Apply an affine shift (translation) to an image."""
    M = np.float32([[1, 0, dx], [0, 1, dy]])
    return cv2.warpAffine(
        img,
        M,
        (img.shape[1], img.shape[0]),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_REFLECT,
    )


def estimate_drift_of(img1: np.ndarray, img2: np.ndarray) -> Tuple[float, float]:
    """
    Robust estimate of the translation (dx, dy) of img2 relative to img1.

    Uses a 2-stage hierarchical approach:
    1. Coarse estimation on downscaled images to handle large displacements.
    2. Fine estimation using an ensemble of crops on the full-resolution images,
       corrected by the coarse shift.
    3. Robust averaging (IQR) to reject outliers.
    """

    def get_phase_corr_shift(im1, im2):
        if np.std(im1) < 1e-3 or np.std(im2) < 1e-3:
            return None
        try:
            detected_shift, _error, _diffphase = phase_cross_correlation(
                preprocess_image(im1),
                preprocess_image(im2),
                upsample_factor=20,
                normalization=None,
            )
            return -detected_shift[1], -detected_shift[0]
        except Exception:
            return None

    scale = 0.25
    h, w = img1.shape
    small_h, small_w = int(h * scale), int(w * scale)

    img1_s = cv2.resize(img1, (small_w, small_h), interpolation=cv2.INTER_AREA)
    img2_s = cv2.resize(img2, (small_w, small_h), interpolation=cv2.INTER_AREA)

    coarse_shift = get_phase_corr_shift(img1_s, img2_s)

    if coarse_shift is None:
        coarse_dx, coarse_dy = 0.0, 0.0
    else:
        coarse_dx = coarse_shift[0] / scale
        coarse_dy = coarse_shift[1] / scale

    img2_aligned = shift_image(img2, -coarse_dx, -coarse_dy)

    shifts = []

    crop_h = int(h * 0.75)
    crop_w = int(w * 0.75)

    crops = [
        (h // 2 - crop_h // 2, h // 2 + crop_h // 2, w // 2 - crop_w // 2, w // 2 + crop_w // 2),
        (0, crop_h, 0, crop_w),
        (0, crop_h, w - crop_w, w),
        (h - crop_h, h, 0, crop_w),
        (h - crop_h, h, w - crop_w, w),
    ]

    for (r1, r2, c1, c2) in crops:
        r1, r2 = max(0, r1), min(h, r2)
        c1, c2 = max(0, c1), min(w, c2)

        try:
            detected_shift, _error, _diffphase = phase_cross_correlation(
                preprocess_image(img1[r1:r2, c1:c2]),
                preprocess_image(img2_aligned[r1:r2, c1:c2]),
                upsample_factor=100,
                normalization=None,
            )
            res_dx = -detected_shift[1]
            res_dy = -detected_shift[0]
            shifts.append((coarse_dx + res_dx, coarse_dy + res_dy))
        except Exception:
            pass

    if not shifts:
        return float(coarse_dx), float(coarse_dy)

    shifts_arr = np.array(shifts)

    def robust_mean(values):
        if len(values) < 3:
            return np.mean(values)
        q1 = np.percentile(values, 25)
        q3 = np.percentile(values, 75)
        iqr = q3 - q1
        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr
        filtered = [v for v in values if lower <= v <= upper]
        if not filtered:
            return np.median(values)
        return np.mean(filtered)

    final_dx = robust_mean(shifts_arr[:, 0])
    final_dy = robust_mean(shifts_arr[:, 1])

    return float(final_dx), float(final_dy)
