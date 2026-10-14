"""AutoEMX version.

Single source of the version: read by pyproject.toml (setuptools dynamic version), docs/conf.py, and stored in
the ledgers AutoEMX writes. Unlike the installed package metadata, it stays correct in editable installs.
"""

__version__ = "1.6.1"

# Version recorded for quantification runs computed before AutoEMX recorded its version
LEGACY_VERSION = "<=1.4.0"
