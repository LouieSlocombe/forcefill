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

## Releasing

Versions are read from one place — `version` in `pyproject.toml`.
`forcefill/__init__.py` and `docs/conf.py` both read it back through the
installed metadata, so there is nothing else to edit.

1. Add the new section to `CHANGELOG.md`, moving anything under `[Unreleased]`
   into it, and update the two link definitions at the bottom of the file.
2. Bump `version` in `pyproject.toml`. From 1.0.0 on this follows semantic
   versioning, where the public API is everything exported from
   `forcefill/__init__.py`.
3. Commit, push, and let CI go green on `main`.
4. Tag and push:

   ```bash
   git tag -a v1.2.3 -m "forcefill 1.2.3"
   git push origin v1.2.3
   ```

The tag triggers `.github/workflows/release.yml`, which refuses to go any
further unless the tag, `pyproject.toml` and `CHANGELOG.md` agree on the version
— a version uploaded to PyPI can never be re-uploaded, so the wrong number
published is the wrong number forever. It then builds the sdist and wheel in the
conda environment, installs the wheel and imports it from outside the source
tree, uploads to PyPI, and creates a GitHub release whose notes are the
changelog section.

`workflow_dispatch` runs the same checks and the same build without publishing
anything, which is the way to try a change to the workflow itself.

### One-time PyPI setup

Uploads use [Trusted Publishing](https://docs.pypi.org/trusted-publishers/), so
no API token is stored in the repository. Before the first tag, on PyPI:

1. Create a [pending publisher](https://pypi.org/manage/account/publishing/) for
   the project name `forcefill` — repository owner `LouieSlocombe`, repository
   `forcefill`, workflow `release.yml`, environment `pypi`. "Pending" is the
   right choice while the project does not exist yet; PyPI creates it on the
   first successful upload.
2. On GitHub, create an environment named `pypi` under **Settings →
   Environments**. A required-reviewer rule there is worth adding: it turns
   every upload into something a human approves.

If this has not been done when the tag is pushed, the `pypi` job fails at the
upload step and no GitHub release is created. The built artifacts are still
attached to the workflow run, so the fix is to configure PyPI and re-run the
failed jobs — the tag does not need deleting and re-pushing.
