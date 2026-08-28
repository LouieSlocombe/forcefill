# Contributing

## Development environment

AmberTools is conda-only, so the conda-forge environment is the supported route
for working on forcefill — the integration tests run the real `antechamber` and
`parmchk2`.

```bash
conda env create -f environment.yml
conda activate forcefill
pip install -e . --no-deps
```

## Tests

```bash
pytest -m "not integration and not smirnoff"   # fast hermetic tests only
pytest                                         # everything, including real antechamber and OpenFF
```

Three markers select the slow legs: `integration` (the real AmberTools
executables), `smirnoff` (real OpenFF parameterization) and `espaloma` (the real
Espaloma model, which needs the optional `espaloma` package and downloads a model
on first use). CI runs all of them and enforces a coverage floor.

## Style

`ruff` owns both layout and lint, and gates CI:

```bash
pip install -e '.[dev]' && pre-commit install
```

The ruff version is pinned in three places that must agree — `.pre-commit-config.yaml`,
the `dev` extra in `pyproject.toml`, and `.github/workflows/ci.yml` — because the
formatter's output shifts between releases, so a skew lets a locally clean tree
fail the CI format gate.

Docstrings are Google-style (`Args:` / `Returns:` / `Raises:`) with reST roles in
the prose, and dataclass fields and module constants are documented with Sphinx
`#:` comments. That is what the API reference is built from, so a new public
function wants the same treatment.

## Documentation

The docs build is pip-only and does **not** use `environment.yml`. Two runtime
dependencies cannot be installed from PyPI — `openmmforcefields` stops at 0.15.1
there, below the 0.16 floor, and `parmed` ships no wheels — so Read the Docs
installs forcefill with `--no-deps` and `docs/conf.py` mocks them. `openmm` is
the one dependency installed for real, because `forcefill/checks.py` divides two
units at import time and a mock cannot do that.

```bash
pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs docs/_build/html
python -m http.server -d docs/_build/html 8000
```

`-W` matches `fail_on_warning: true` in `.readthedocs.yaml` and the `docs` job in
CI, which reproduces the Read the Docs environment exactly — mocks included — so
that a mock which stops working fails on a pull request rather than on the docs
site.
