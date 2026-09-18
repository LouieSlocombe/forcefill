# Installation

AmberTools is conda-only, so conda-forge is the recommended route:

```bash
conda env create -f environment.yml
conda activate forcefill
pip install -e . --no-deps
```

Or into an existing environment:

```bash
conda install -c conda-forge openmm parmed ambertools
pip install forcefill
```

Requires Python ≥ 3.12 and, at import time, `openmm >= 7.6`, `parmed >= 3.4`,
`rdkit`, `openff-toolkit >= 0.16` and `openmmforcefields >= 0.16` — all ordinary
dependencies, with no extras to pick and nothing imported lazily.

The openmmforcefields floor is 0.16 and not lower: that release is where
`smirnoff_filenames`, a multi-file `forcefield=` selection, and constraints and
virtual sites in the generated template all arrive. On 0.15 and earlier the
smirnoff backend does not work at all.

Two things are not ordinary dependencies. AmberTools is not a Python package:
the `antechamber` and `parmchk2` executables must be on `PATH` at run time for
the `gaff` backend, and no other backend needs them. And `espaloma` is optional,
because it pulls in PyTorch — install it with
`conda install -c conda-forge espaloma` if you want `backend="espaloma"`;
forcefill says so by name if you ask for it without.

