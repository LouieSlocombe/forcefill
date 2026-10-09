#!/usr/bin/env bash
# Run in the checkout with the scientific dependencies already installed.
set -euo pipefail

python -m build
python -m pip install --force-reinstall --no-deps dist/*.whl
(
    cd "$(mktemp -d)"  # import from outside the source tree
    python -c "import forcefill; print(forcefill.__version__)"
    python -c "import importlib.resources as r; assert r.files('forcefill').joinpath('py.typed').is_file()"
)
