#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Install the optional particle segmentation models in the current Python environment.

    python -m autoemx.utils.install_segmentation sam2 [--dry-run]
    python -m autoemx.utils.install_segmentation micro-sam [--dry-run]

``sam2``: SAM2 through Hugging Face transformers, the ``segmentation`` extra of AutoEMX (torch, torchvision,
transformers; Apache-2.0 / BSD).

``micro-sam``: micro-sam (MIT; its model weights are CC-BY-4.0, downloaded from BioImage.IO on first use) is installed
without its declared dependencies, which include PyQt6 (GPL) for its napari plugin, not needed to segment images,
and ``segment-anything-py``, a third-party repackaging of Meta's ``segment_anything`` (installed instead from
Meta's repository, Apache-2.0). Its other dependencies are then installed (``MICRO_SAM_REQUIREMENTS``).

Created on Sat Oct 17 2026
@author: Andrea
"""

import argparse
import subprocess
import sys
from typing import List

MICRO_SAM_VERSION = "1.8.14"
SEGMENT_ANYTHING = "git+https://github.com/facebookresearch/segment-anything.git"

#: Dependencies of micro-sam used to segment images: its declared ones without PyQt6, segment-anything-py, and
#: trackastra (tracking only); magicgui and superqt (MIT/BSD) come with napari, without any Qt binding
MICRO_SAM_REQUIREMENTS = [
    "torch>=2.5",
    "torchvision",
    "torch-em>=0.9",
    "python-elf>=0.9",
    "bioimage-cpp>=0.3",
    "bioimageio.core>=0.11",
    "napari>=0.7.1,<0.10",
    "h5py",
    "imagecodecs",
    "imageio",
    "joblib",
    "kornia",
    "natsort",
    "networkx",
    "pooch",
    "timm",
    "tqdm",
    "xarray",
    "xxhash",
    "zarr",
]

#: Packages that must not be installed (GPL)
EXCLUDED = ("PyQt6", "PyQt5")


def install_commands(model: str) -> List[List[str]]:
    """The pip commands installing a model."""
    pip = [sys.executable, "-m", "pip", "install"]
    if model == "sam2":
        return [pip + ["torch>=2.1", "torchvision>=0.16", "transformers>=4.56,<6"]]
    if model == "micro-sam":
        return [
            pip + ["--no-deps", f"micro_sam=={MICRO_SAM_VERSION}", SEGMENT_ANYTHING],
            pip + MICRO_SAM_REQUIREMENTS,
        ]
    raise ValueError(f"Unknown model '{model}': 'sam2' or 'micro-sam'")


def gpl_packages_installed() -> List[str]:
    """Installed packages among those that must not be (EXCLUDED)."""
    from importlib.metadata import PackageNotFoundError, version

    found = []
    for name in EXCLUDED:
        try:
            version(name)
            found.append(name)
        except PackageNotFoundError:
            pass
    return found


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Install an optional particle segmentation model of AutoEMX.")
    parser.add_argument("model", choices=["sam2", "micro-sam"])
    parser.add_argument("--dry-run", action="store_true", help="print the pip commands without running them")
    args = parser.parse_args(argv)
    for command in install_commands(args.model):
        print(" ".join(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, check=True)
    if not args.dry_run:
        found = gpl_packages_installed()
        if found:
            print(f"Warning: {', '.join(found)} installed (GPL), by another package of this environment.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
