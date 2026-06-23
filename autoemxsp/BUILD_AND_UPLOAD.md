# Upload `autoemxsp` redirect release to PyPI

Build and upload from this directory (`autoemxsp/`), not from the repo root.

## Prerequisites

1. **`autoemx` 1.0.0** must already be on PyPI (the redirect depends on it).
2. PyPI API token for the **`autoemxsp`** project owner account.

## Build

```bash
cd autoemxsp
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
```

Expected artifacts:

- `dist/autoemxsp-0.1.8-py3-none-any.whl`
- `dist/autoemxsp-0.1.8.tar.gz`

## Upload

```bash
export TWINE_USERNAME='__token__'
export TWINE_PASSWORD='pypi-...'   # your API token
python -m twine upload dist/*
```

On Windows PowerShell:

```powershell
$env:TWINE_USERNAME = "__token__"
$env:TWINE_PASSWORD = "pypi-..."
python -m twine upload dist/*
```

## Verify after upload

```bash
pip install --upgrade autoemxsp
python -c "import autoemxsp; import autoemx; print('autoemx', autoemx.__file__)"
```

You should see a `FutureWarning` about the rename, and `autoemx` should be installed as a dependency.

## Notes

- Version **0.1.8** is the compatibility redirect release after **0.1.7**.
- Do **not** delete the old `autoemxsp` releases on PyPI; leave them for pinned installs.
- `autoemxsp/autoemxsp/__init__.py` is the single compatibility shim used by both the redirect wheel and the main `autoemx` source tree.
