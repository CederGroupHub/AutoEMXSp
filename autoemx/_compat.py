"""Compatibility helpers for package-rename transition."""

from __future__ import annotations

import glob
import os
import shlex
import site
import sys
import warnings


def _candidate_site_packages() -> list[str]:
    paths: set[str] = set()
    for getter in (site.getsitepackages,):
        try:
            for p in getter():
                if p:
                    paths.add(os.path.abspath(p))
        except Exception:
            continue

    try:
        user_site = site.getusersitepackages()
        if user_site:
            paths.add(os.path.abspath(user_site))
    except Exception:
        pass

    return sorted(paths)


def _find_stale_autoemxsp_artifacts() -> list[str]:
    stale: list[str] = []
    for site_pkg in _candidate_site_packages():
        ns_dir = os.path.join(site_pkg, "autoemxsp")
        init_file = os.path.join(ns_dir, "__init__.py")
        if os.path.isdir(ns_dir) and not os.path.exists(init_file):
            stale.append(ns_dir)

        stale.extend(glob.glob(os.path.join(site_pkg, "__editable__.autoemxsp-*.pth")))
        stale.extend(glob.glob(os.path.join(site_pkg, "__editable___autoemxsp_*_finder.py")))

    # De-duplicate while preserving order.
    unique: list[str] = []
    seen: set[str] = set()
    for path in stale:
        norm = os.path.abspath(path)
        if norm not in seen:
            seen.add(norm)
            unique.append(norm)
    return unique


def _build_cleanup_message(paths: list[str]) -> str:
    exe = shlex.quote(sys.executable)
    listed_paths = "\n".join(f"  {p}" for p in paths)
    return (
        "Leftover files from an old installation of 'autoemxsp' (the former name of AutoEMX) were found:\n"
        f"{listed_paths}\n"
        "They only affect scripts that still use 'import autoemxsp' (e.g. 'autoemxsp.runners'), which may fail "
        "to import. AutoEMX itself is not affected.\n"
        "To remove them:\n"
        f"  1. {exe} -m pip uninstall -y autoemxsp   (skip if pip reports it is not installed)\n"
        "  2. Delete the files and folders listed above."
    )


def warn_if_stale_autoemxsp_install() -> None:
    """Warn users when stale legacy package artifacts may shadow the shim."""
    stale_paths = _find_stale_autoemxsp_artifacts()
    if not stale_paths:
        return

    warnings.warn(_build_cleanup_message(stale_paths), UserWarning, stacklevel=2)
