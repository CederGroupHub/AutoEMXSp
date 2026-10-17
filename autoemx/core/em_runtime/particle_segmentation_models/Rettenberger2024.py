"""
Created on Thu Oct  9 09:34:39 2025

Particle segmentation model from:
    
    Rettenberger, L., Szymanski, N.J., Zeng, Y. et al. Uncertainty-aware particle segmentation for electron microscopy at varied length scales. npj Comput Mater 10, 124 (2024). https://doi.org/10.1038/s41524-024-01302-w

To use this model for particle segmentation, clone the repo from https://github.com/lrettenberger/Uncertainty-Aware-Particle-Segmentation-for-SEM to the same directory as this module, and rename it Rettenberg2024_model_data

@author: Andrea
"""

import os
from pathlib import Path

import cv2
import numpy as np

from autoemx._logging import get_logger
logger = get_logger(__name__)

INSTALL_HINT = ("pip install onnxruntime, then clone https://github.com/lrettenberger/"
                "Uncertainty-Aware-Particle-Segmentation-for-SEM into this module's directory as "
                "Rettenberg2024_model_data")


def is_available() -> bool:
    """Whether onnxruntime and the model files are installed."""
    import importlib.util
    model_path = Path(__file__).resolve().parent / "Rettenberg2024_model_data" / "highmag" / "model_high_mag_maskedrcnn.onnx"
    return importlib.util.find_spec("onnxruntime") is not None and model_path.exists()

def segment_particles(frame_image: np.ndarray,
                      powder_meas_config : 'PowderMeasurementConfig' = None,
                      save_image: bool = False,
                      EM: 'EM_controller' = None) -> np.ndarray:
    """
    Segments particles in the given frame image using a Mask R-CNN ONNX model, then
    returns an 8-bit *index map* image where each detected particle is assigned a
    unique brightness value based on its index.

    Parameters
    ----------
    frame_image : ndarray
        A grayscale (or RGB) input image of the current frame containing particles.

    save_image : bool
        Optionally save the 8-bit index map through EM_controller.

    EM : EM_controller object
        Used to optionally save the 8-bit index map image.

    Returns
    -------
    par_mask : ndarray (uint8)
        An 8-bit index image where:

        - Background pixels are 0.
        - Each particle (up to 255 particles) receives an index 1..255 and all of
          its pixels are set to that index value. If >255 particles are found,
          the extras are set to 255.

    Note
    ----
    - All configuration is hardcoded inside this function.
    - The ONNX model path is resolved *relative* to this file.
    """
    import onnxruntime as ort
    # --- Replace the failing import with this local import ---
    from pathlib import Path
    import sys

    HERE = Path(__file__).resolve().parent
    # Ensure this file's directory is importable
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))

    # Import the model helper relative to THIS file
    from Rettenberg2024_model_data.highmag.model import get_segmentation_mask


    # -----------------------
    # Hardcoded configuration
    # -----------------------
    USE_CLAHE: bool      = True
    CLAHE_CLIP: float    = 2.0
    CLAHE_TILE: tuple    = (8, 8)

    UPSCALE_FACTOR: float = 1.0  # >1.0 upscales before inference

    MIN_AREA_PX: int      = 5    # drop instances smaller than this
    MERGE_UNSURE_PX: int  = 40   # merge tiny unsure instances into confident set if non-overlapping

    # -----------------------
    # Resolve model path (relative to this file)
    # Expected layout:
    #   <this_file_dir>/
    #       Rettenberg2024_model_data/highmag/model_high_mag_maskedrcnn.onnx
    # -----------------------
    here = Path(__file__).resolve().parent
    model_path = (here / "Rettenberg2024_model_data" / "highmag" / "model_high_mag_maskedrcnn.onnx").resolve()

    if not model_path.exists():
        raise FileNotFoundError(
            f"segment_particles: Could not find ONNX model at:\n{model_path}\n"
            "Ensure the model is placed relative to this script as shown."
        )

    # -----------------------
    # Helpers (minimal, local)
    # -----------------------
    def _ensure_gray_uint8(img: np.ndarray) -> np.ndarray:
        if img.ndim == 3:
            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        if img.dtype != np.uint8:
            img = cv2.normalize(img, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        return img

    def _apply_clahe_if_needed(img: np.ndarray) -> np.ndarray:
        if not USE_CLAHE:
            return img
        clahe = cv2.createCLAHE(clipLimit=float(CLAHE_CLIP), tileGridSize=tuple(CLAHE_TILE))
        return clahe.apply(img)

    def _upscale(img: np.ndarray):
        f = float(UPSCALE_FACTOR)
        if abs(f - 1.0) < 1e-6:
            return img, (lambda m: m)
        h, w = img.shape[:2]
        up_w = int(round(w * f))
        up_h = int(round(h * f))
        img_up = cv2.resize(img, (up_w, up_h), interpolation=cv2.INTER_CUBIC)

        def down_fn(m: np.ndarray):
            return cv2.resize(m.astype(np.int32), (w, h), interpolation=cv2.INTER_NEAREST).astype(m.dtype)

        return img_up, down_fn

    def _drop_small_instances(label_img: np.ndarray, min_area_px: int) -> np.ndarray:
        if label_img is None or label_img.max() == 0 or min_area_px <= 1:
            return label_img
        out = np.zeros_like(label_img)
        next_id = 1
        for lbl in np.unique(label_img):
            if lbl == 0:
                continue
            mask = (label_img == lbl)
            if int(mask.sum()) >= min_area_px:
                out[mask] = next_id
                next_id += 1
        return out

    def _merge_tiny_unsure(confident_lbl: np.ndarray, unsure_lbl: np.ndarray, max_area_px: int) -> np.ndarray:
        if max_area_px <= 0:
            return confident_lbl
        out = confident_lbl.copy()
        next_id = int(out.max()) + 1
        for lbl in np.unique(unsure_lbl):
            if lbl == 0:
                continue
            mask = (unsure_lbl == lbl)
            a = int(mask.sum())
            if a <= max_area_px and not np.any(out[mask] > 0):
                out[mask] = next_id
                next_id += 1
        return out

    # -----------------------
    # Prepare image
    # -----------------------
    gray = _ensure_gray_uint8(frame_image)
    gray = _apply_clahe_if_needed(gray)
    img_for_net, down_fn = _upscale(gray)

    # -----------------------
    # ONNX Runtime session (cached per model path)
    # -----------------------
    need_new_session = (
        not hasattr(segment_particles, "_ort_sess") or
        not hasattr(segment_particles, "_ort_model_path") or
        segment_particles._ort_model_path != str(model_path)
    )
    if need_new_session:
        available = ort.get_available_providers()
        preferred = ['CUDAExecutionProvider', 'DmlExecutionProvider', 'CPUExecutionProvider']
        providers = [p for p in preferred if p in available] or ['CPUExecutionProvider']
        segment_particles._ort_sess = ort.InferenceSession(str(model_path), providers=providers)
        segment_particles._ort_model_path = str(model_path)

    ort_sess = segment_particles._ort_sess

    # -----------------------
    # Inference (repo helper)
    # -----------------------
    m_not_conf, m_conf, m_comb = get_segmentation_mask(img_for_net, ort_sess)

    # Downscale masks back if we upscaled
    if abs(float(UPSCALE_FACTOR) - 1.0) > 1e-6:
        m_not_conf = down_fn(m_not_conf)
        m_conf     = down_fn(m_conf)
        m_comb     = down_fn(m_comb)

    # Clean up (tiny specks)
    m_not_conf = _drop_small_instances(m_not_conf, MIN_AREA_PX)
    m_conf     = _drop_small_instances(m_conf,     MIN_AREA_PX)
    m_comb     = _drop_small_instances(m_comb,     MIN_AREA_PX)

    # Optionally “rescue” tiny unsure instances into confident set
    m_conf = _merge_tiny_unsure(m_conf, m_not_conf, MERGE_UNSURE_PX)

    # -----------------------
    # Build 8-bit index map from combined mask
    # -----------------------
    labels_comb = np.unique(m_comb)
    labels_comb = labels_comb[labels_comb != 0]

    index_map = np.zeros(m_comb.shape, dtype=np.uint8)
    if labels_comb.size > 0:
        if labels_comb.size > 255:
            logger.warning(f"⚠️ segment_particles: {labels_comb.size} particles detected; "
                  f"capping index map at 255 (extra particles set to 255).")

        for idx, lbl in enumerate(labels_comb[:255], start=1):
            index_map[m_comb == lbl] = idx
        if labels_comb.size > 255:
            for lbl in labels_comb[255:]:
                index_map[m_comb == lbl] = 255

       # --- Multicolor EDGE overlay: draw only particle boundaries in unique colors ---
    base = cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR)  # keep grayscale background for context
    num_labels = int(index_map.max())

    # Color lookup table (0..255) in BGR; label 0 stays black (unused)
    color_lut = np.zeros((256, 3), dtype=np.uint8)
    if num_labels > 0:
        # Evenly spaced hues in OpenCV HSV space (H:0..179)
        hues = (np.linspace(0, 179, num=num_labels, endpoint=False).astype(np.uint8) + 37) % 180
        sats = np.full(num_labels, 200, dtype=np.uint8)
        vals = np.full(num_labels, 255, dtype=np.uint8)
        hsv = np.stack([hues, sats, vals], axis=1)
        bgr = cv2.cvtColor(hsv.reshape(-1, 1, 3), cv2.COLOR_HSV2BGR).reshape(-1, 3)
        color_lut[1:num_labels + 1] = bgr

    # Start from the grayscale image and draw only edges (no interior alpha-blend)
    overlay = base.copy()
    EDGE_THICKNESS = 2  # change to 1 for thinner lines

    # For each label, find outer contours and draw them in that label's color
    for lbl in range(1, num_labels + 1):
        m = (index_map == lbl).astype(np.uint8) * 255
        if not m.any():
            continue
        cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        c = tuple(int(v) for v in color_lut[lbl])  # BGR tuple
        cv2.drawContours(overlay, cnts, -1, c, EDGE_THICKNESS)



    if save_image and EM:
        # Ensure m_comb is OpenCV-friendly (uint8 or uint16)
        if m_comb.dtype != np.uint8:
            # If labels are <=255, cast to uint8; else cast to uint16
            max_val = m_comb.max()
            if max_val <= 255:
                m_comb = m_comb.astype(np.uint8)
            elif max_val <= 65535:
                m_comb = m_comb.astype(np.uint16)
            else:
                raise ValueError(f"m_comb contains label {max_val} > 65535, unsupported for direct saving.")

        filename = f"{EM.sample_id}_fr{EM.current_frame_label}_Rettenberger2024_mask"
        EM.save_frame_image(filename, frame_image = overlay)

    # Return the 8-bit index map
    return index_map


if __name__ == "__main__":
    # ===== TESTING CONTINUATION =====
    # Resolve a bundled example image relative to this file:
    # Expected at:
    #   <this_file_dir>/example_figure.tiff
    test_here = Path(__file__).resolve().parent
    TEST_IMAGE =  r"C:\Users\ThermoFisher\AppData\Local\Programs\Python\Python311\Lib\site-packages\autoemx\core\particle_segmentation_models/example_figure.tiff"

    import matplotlib.pyplot as plt

    # Load test image (TIFF). If it’s RGB, convert to grayscale.
    img = cv2.imread(str(TEST_IMAGE), cv2.IMREAD_UNCHANGED)
    if img is None:
        raise FileNotFoundError(f"Could not read image: {TEST_IMAGE}")
    if img.ndim == 3:
        img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    else:
        img_gray = img

    # Run segmentation → returns 8-bit index map (0 = background, 1..255 = particles)
    index_map = segment_particles(img_gray, powder_meas_config=None, save_image=False, EM=None)


    # Display the result
    plt.figure(figsize=(8, 6))
    plt.imshow(index_map, cmap="gray", vmin=0, vmax=255)
    plt.title("Particle Index Map (8-bit)")
    plt.axis("off")
    plt.show()
