#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Launch the AutoEMX spectrum GUI: ``python -m autoemx.web``."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        print(
            "The AutoEMX GUI requires Streamlit.\n"
            "Install it with:  pip install streamlit\n"
            "Then run:         python -m autoemx.web",
            file=sys.stderr,
        )
        sys.exit(1)

    app_path = Path(__file__).with_name("app.py")
    sys.argv = [
        "streamlit",
        "run",
        str(app_path),
        "--browser.gatherUsageStats=false",
    ]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
