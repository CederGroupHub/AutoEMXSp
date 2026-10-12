#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Launch the AutoEMX GUI: ``python -m autoemx.gui [results_folder] [--samples samples.json]``."""

from __future__ import annotations

from typing import Optional

import argparse
import hashlib
import json
import os
import socket
import stat
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
from pathlib import Path
from urllib.parse import urlencode

from autoemx._paths import AUTOEMX_HOME

LAUNCHER_NAME = "AutoEMX"
OLD_LAUNCHER_NAMES = ("AutoEMX GUI",)  # replaced by LAUNCHER_NAME when the launcher is created again
LAUNCHER_ICON = Path(__file__).with_name("assets") / "autoemx-icon-512.png"
LAUNCHER_ICON_WINDOWS = Path(__file__).with_name("assets") / "autoemx-icon.ico"
# On Windows, the .bat run by the Desktop shortcut is kept here (a .bat cannot have its own icon)
WINDOWS_LAUNCHERS_DIR = AUTOEMX_HOME / "launchers"
# Output of the GUIs started in the background by open_acquisition (macOS, Linux)
BACKGROUND_LOGS_DIR = AUTOEMX_HOME / "logs"


def acquisition_url(samples: Optional[list] = None, results_folder: Optional[str] = None, port: int = 8050,
                    run_id: Optional[str] = None) -> str:
    """
    URL of the GUI opening the Acquisition tab with these samples and results folder filled in.

    *samples* are as in ``Run_Acquisition.py``: ``{'ID', 'els', 'pos': (x, y), 'cnd'}``. With a *run_id*,
    the runs started from that page are reported (see ``acquisition_report``).
    """
    url = f"http://127.0.0.1:{port}/"
    acq: dict = {}
    if results_folder:
        acq["folder"] = os.path.abspath(os.path.expanduser(results_folder))
    if samples is not None:
        acq["samples"] = list(samples)
    if run_id is not None:
        from autoemx.gui.run_report import clean_run_id

        acq["run_id"] = clean_run_id(run_id)
    return f"{url}?{urlencode({'acq': json.dumps(acq)})}" if acq else url


def _read_samples_file(path: str) -> tuple:
    """Samples, results folder and run ID (None if not given) of a JSON file: a list of samples, or
    ``{"folder", "samples", "run_id"}``."""
    with open(os.path.expanduser(path), encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list):
        return data, None, None
    if isinstance(data, dict) and isinstance(data.get("samples"), list):
        return data["samples"], data.get("folder"), data.get("run_id")
    raise ValueError('expected a list of samples, or {"folder": ..., "samples": [...], "run_id": ...}')


def _port_in_use(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        return sock.connect_ex(("127.0.0.1", port)) == 0


def open_acquisition(samples: list, results_folder: Optional[str] = None, port: int = 8050,
                     run_id: Optional[str] = None, open_browser: bool = True, timeout: float = 120.0) -> str:
    """
    Open the Acquisition tab of the GUI with these samples and results folder filled in, and return the ID
    of the run (*run_id*, or a new one), to get its report with ``acquisition_report`` or
    ``wait_for_acquisition``. Meant for external scripts: it returns once the GUI is ready, without waiting
    for the acquisition.

    *samples* are as in ``Run_Acquisition.py``: ``{'ID', 'els', 'pos': (x, y), 'cnd'}``. *run_id*: letters,
    digits, '-' or '_'. If the user starts several acquisitions from the page, the report is of the last one.

    If a GUI already runs on *port*, it is used. Otherwise one is started in the background, with the
    current Python interpreter, and keeps running after the calling script ends:

    - Windows: in its own console window; close it to stop the GUI.
    - macOS, Linux: without a window, its output written to ``~/.autoemx/logs/gui_<port>.log``;
      stop it with ``kill`` (the PID is in the log).

    Raises ``RuntimeError`` if the GUI could not be started within *timeout* seconds.
    """
    run_id = run_id if run_id is not None else uuid.uuid4().hex[:12]
    url = acquisition_url(samples, results_folder, port, run_id)  # checks the run ID and that samples are JSON
    if not _port_in_use(port):
        cmd = [sys.executable, "-m", "autoemx.gui", "--port", str(port), "--no-browser"]
        if results_folder:
            cmd.append(os.path.abspath(os.path.expanduser(results_folder)))
        env = dict(os.environ)
        package_parent = str(Path(__file__).resolve().parents[2])  # also works if AutoEMX is not installed
        env["PYTHONPATH"] = os.pathsep.join(p for p in (package_parent, env.get("PYTHONPATH")) if p)
        if os.name == "nt":
            log = None
            proc = subprocess.Popen(cmd, env=env, cwd=str(Path.home()),
                                    creationflags=subprocess.CREATE_NEW_CONSOLE)
        else:
            BACKGROUND_LOGS_DIR.mkdir(parents=True, exist_ok=True)
            log = BACKGROUND_LOGS_DIR / f"gui_{port}.log"
            with open(log, "w", encoding="utf-8") as out:
                proc = subprocess.Popen(cmd, env=env, cwd=str(Path.home()), stdin=subprocess.DEVNULL,
                                        stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
                out.write(f"AutoEMX GUI, PID {proc.pid}\n")
        deadline = time.monotonic() + timeout
        while not _port_in_use(port):
            if proc.poll() is not None or time.monotonic() > deadline:
                if proc.poll() is None:
                    proc.kill()
                where = f" See {log}." if log else ""
                raise RuntimeError(f"The AutoEMX GUI could not be started on port {port}.{where}")
            time.sleep(0.3)
    if open_browser:
        webbrowser.open(url)
    return run_id


def _set_macos_icon(path: Path) -> bool:
    """Give a file the AutoEMX icon in Finder (best effort; uses the AppKit bridge of osascript)."""
    script = (
        "ObjC.import('AppKit');"
        f"var img = $.NSImage.alloc.initWithContentsOfFile({json.dumps(str(LAUNCHER_ICON))});"
        f"$.NSWorkspace.sharedWorkspace.setIconForFileOptions(img, {json.dumps(str(path))}, 0);"
    )
    try:
        res = subprocess.run(["osascript", "-l", "JavaScript", "-e", script],
                             capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return False
    return res.returncode == 0 and res.stdout.strip() == "true"


def _hide_macos_extension(path: Path) -> bool:
    """Hide the extension of a file in Finder, so that the launcher shows as an app name (best effort)."""
    script = (
        "ObjC.import('Foundation');"
        "$.NSFileManager.defaultManager.setAttributesOfItemAtPathError("
        f"$({{NSFileExtensionHidden: true}}), {json.dumps(str(path))}, null);"
    )
    try:
        res = subprocess.run(["osascript", "-l", "JavaScript", "-e", script],
                             capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return False
    return res.returncode == 0 and res.stdout.strip() == "true"


def _powershell(script: str) -> Optional[str]:
    """Run a PowerShell script (Windows); its output, or None if it failed."""
    try:
        res = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", script],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return res.stdout.strip() if res.returncode == 0 else None


def _ps_quote(text: str) -> str:
    """Single-quoted PowerShell string literal."""
    return "'" + str(text).replace("'", "''") + "'"


def _windows_desktop() -> Path:
    """The user's Desktop folder on Windows (it may be redirected, e.g. to OneDrive)."""
    out = _powershell("[Environment]::GetFolderPath('Desktop')")
    return Path(out) if out and Path(out).is_dir() else Path.home() / "Desktop"


def _windows_shortcut_script(shortcut: Path, target: Path, icon: Path, workdir: Path) -> str:
    """PowerShell script creating a shortcut (.lnk) to *target*, with the icon *icon*."""
    return (
        "$s = (New-Object -ComObject WScript.Shell).CreateShortcut(" + _ps_quote(shortcut) + "); "
        "$s.TargetPath = " + _ps_quote(target) + "; "
        "$s.WorkingDirectory = " + _ps_quote(workdir) + "; "
        "$s.IconLocation = " + _ps_quote(f"{icon},0") + "; "
        "$s.Description = 'AutoEMX GUI'; "
        "$s.Save()"
    )


def _remove_old_launchers(dest: Path, extension: str, names=OLD_LAUNCHER_NAMES) -> None:
    """Delete launchers written under a previous name in the same folder (only if they launch the GUI)."""
    for name in names:
        old = dest / f"{name}{extension}"
        try:
            if old.is_file() and "-m autoemx.gui" in old.read_text(encoding="utf-8", errors="replace"):
                old.unlink()
        except OSError:
            pass


def create_launcher(dest_dir: Optional[str] = None, results_folder: Optional[str] = None) -> Path:
    """
    Write a double-clickable launcher of the GUI that runs it with the current Python interpreter
    and AutoEMX installation, and return its path.

    - macOS: ``AutoEMX.command``, with the AutoEMX icon and its extension hidden.
    - Windows: an ``AutoEMX`` shortcut (``.lnk``) with the AutoEMX icon, running a ``.bat`` kept in
      ``~/.autoemx/launchers`` (a ``.bat`` cannot have its own icon). If the shortcut cannot be
      created, the ``.bat`` is written in the destination folder instead.
    - Linux: ``AutoEMX.sh``.

    Launchers written under a previous name in the same folder are replaced.
    """
    if dest_dir:
        dest = Path(dest_dir).expanduser()
    else:
        dest = _windows_desktop() if os.name == "nt" else Path.home() / "Desktop"
    if not dest.is_dir():
        dest = Path.home()
    python = sys.executable
    package_parent = str(Path(__file__).resolve().parents[2])
    args = f' "{os.path.abspath(results_folder)}"' if results_folder else ""
    extension = ".bat" if os.name == "nt" else ".command" if sys.platform == "darwin" else ".sh"
    _remove_old_launchers(dest, extension)
    path = dest / f"{LAUNCHER_NAME}{extension}"
    if os.name == "nt":
        bat_text = (
            "@echo off\r\n"
            "rem Double-click to open the AutoEMX GUI. Close this window to stop it.\r\n"
            f'set "PYTHONPATH={package_parent};%PYTHONPATH%"\r\n'
            f'"{python}" -m autoemx.gui{args} %*\r\n'
            "pause\r\n"
        )
        # One .bat per shortcut (destination and results folder), so that shortcuts do not overwrite each other
        key = hashlib.sha1(f"{dest.resolve()}|{args}".encode()).hexdigest()[:8]
        bat = WINDOWS_LAUNCHERS_DIR / f"{LAUNCHER_NAME}_{key}.bat"
        shortcut = dest / f"{LAUNCHER_NAME}.lnk"
        try:
            WINDOWS_LAUNCHERS_DIR.mkdir(parents=True, exist_ok=True)
            bat.write_text(bat_text, encoding="utf-8")
            created = _powershell(_windows_shortcut_script(
                shortcut, bat, LAUNCHER_ICON_WINDOWS, Path.home())) is not None and shortcut.exists()
        except OSError:
            created = False
        if created:
            _remove_old_launchers(dest, ".bat", OLD_LAUNCHER_NAMES + (LAUNCHER_NAME,))
            return shortcut
        path.write_text(bat_text, encoding="utf-8")  # no shortcut: plain .bat, without icon
    else:
        path.write_text(
            "#!/bin/bash\n"
            "# Double-click to open the AutoEMX GUI. Close this window to stop it.\n"
            f'export PYTHONPATH="{package_parent}${{PYTHONPATH:+:$PYTHONPATH}}"\n'
            f'"{python}" -m autoemx.gui{args} "$@"\n',
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        if sys.platform == "darwin":
            _set_macos_icon(path)
            _hide_macos_extension(path)
    return path


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m autoemx.gui",
        description="AutoEMX GUI: quantification, clustering analysis and single-spectrum fits of AutoEMX samples.",
    )
    parser.add_argument("results_folder", nargs="?", default=None,
                        help="Folder containing the sample folders (default: ./Results if present, else none).")
    parser.add_argument("--samples", default=None, metavar="JSON",
                        help="JSON file of samples to acquire, as in Run_Acquisition.py (a list, or "
                             '{"folder": ..., "samples": [...], "run_id": ...}): opens the Acquisition tab with '
                             "them filled in. With a run_id, each run is reported (see autoemx.gui.acquisition_report).")
    parser.add_argument("--port", type=int, default=8050, help="Local port (default: 8050).")
    parser.add_argument("--no-browser", action="store_true", help="Do not open a browser tab.")
    parser.add_argument("--debug", action="store_true", help="Dash debug mode (auto-reload, error pop-ups).")
    parser.add_argument("--create-launcher", nargs="?", const="", default=None, metavar="FOLDER",
                        help="Write a double-clickable launcher (default: on the Desktop) and exit. "
                             "If results_folder is given, the launcher opens it.")
    args = parser.parse_args(argv)

    if args.create_launcher is not None:
        path = create_launcher(args.create_launcher or None, args.results_folder)
        print(f"Launcher written to: {path}\nDouble-click it to open the AutoEMX GUI.")
        return

    import importlib.util

    if importlib.util.find_spec("dash") is None:
        print(
            "The AutoEMX GUI requires Dash.\n"
            "Install it with:  pip install dash\n"
            "Then run:         python -m autoemx.gui",
            file=sys.stderr,
        )
        sys.exit(1)

    url = f"http://127.0.0.1:{args.port}/"
    open_url = url
    if args.samples:
        try:
            samples, samples_folder, run_id = _read_samples_file(args.samples)
            open_url = acquisition_url(samples, args.results_folder or samples_folder, args.port, run_id)
        except (OSError, ValueError) as exc:
            print(f"Cannot read the samples file {args.samples}: {exc}", file=sys.stderr)
            sys.exit(1)
    if _port_in_use(args.port):
        # Most likely the GUI is already running (e.g. launcher double-clicked twice).
        # The samples, if any, are passed to it in the URL.
        print(f"Port {args.port} is already in use; opening {url}. Use --port to start another instance.")
        if not args.no_browser:
            webbrowser.open(open_url)
        elif open_url != url:
            print(f"Open the samples with: {open_url}")
        return

    os.environ.setdefault("MPLBACKEND", "Agg")
    from autoemx.gui.app import JOBS, create_app

    folder = args.results_folder
    if folder is None and os.path.isdir(os.path.join(os.getcwd(), "Results")):
        folder = os.path.join(os.getcwd(), "Results")
    app = create_app(os.path.abspath(folder) if folder else None)

    print(f"AutoEMX GUI running at {url}  (Ctrl+C to stop)")
    if args.no_browser and open_url != url:
        print(f"Open the samples with: {open_url}")
    if not args.no_browser:
        threading.Timer(1.2, lambda: webbrowser.open(open_url)).start()
    try:
        # Bound to localhost only: the GUI reads and writes your local files.
        app.run(host="127.0.0.1", port=args.port, debug=args.debug, use_reloader=False)
    finally:
        JOBS.shutdown()


if __name__ == "__main__":
    main()
