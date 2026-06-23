#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np

from autoemx.core.em_runtime.drift_correction import estimate_drift_of, shift_image


def load_grayscale_image(image_path: Path) -> np.ndarray:
    image = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise FileNotFoundError(f"Could not load image: {image_path}")
    return image


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    inputs_dir = script_dir / "inputs"
    outputs_dir = script_dir / "outputs"
    outputs_dir.mkdir(exist_ok=True)

    original = load_grayscale_image(inputs_dir / "pre_drift.png")
    shifted = load_grayscale_image(inputs_dir / "post_drift.png")

    estimated_dx, estimated_dy = estimate_drift_of(original, shifted)
    reshifted = shift_image(shifted, -estimated_dx, -estimated_dy)
    difference = np.abs(original.astype(np.float32) - reshifted.astype(np.float32))

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    fig.suptitle(
        "Image drift correction algorithm\n"
        f"estimated shift ({estimated_dx:.2f}, {estimated_dy:.2f})",
        fontsize=14,
    )

    panels = [
        (axes[0, 0], original, "Original"),
        (axes[0, 1], shifted, "Shifted"),
        (axes[1, 0], reshifted, "Reshifted"),
        (axes[1, 1], difference, "Pixelwise Error"),
    ]

    for axis, image, title in panels:
        cmap = "inferno" if title == "Pixelwise Error" else "gray"
        axis.imshow(image, cmap=cmap)
        axis.set_title(title)
        axis.axis("off")

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    figure_path = outputs_dir / "drift_correction_visualization.png"
    fig.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
