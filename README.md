# forcefill

[![ci](https://github.com/LouieSlocombe/forcefill/actions/workflows/ci.yml/badge.svg)](https://github.com/LouieSlocombe/forcefill/actions/workflows/ci.yml)
[![docs](https://readthedocs.org/projects/forcefill/badge/?version=latest)](https://forcefill.readthedocs.io/en/latest/)

Turn ligands into a ready-to-use [OpenMM](https://openmm.org) force-field XML —
either the non-standard residues (ligands, cofactors, hetero molecules) found in
a PDB, or ligand files on their own. Parameters come from AmberTools
(`antechamber` + `parmchk2`) via [ParmEd](https://github.com/ParmEd/ParmEd),
from [OpenFF](https://openforcefield.org) Sage or
[Espaloma](https://github.com/choderalab/espaloma) via
[openmmforcefields](https://github.com/openmm/openmmforcefields), or — for
CHARMM — by converting the CGenFF stream file
[ParamChem](https://cgenff.paramchem.org) gave you.

The output is a plain ffxml file you load alongside the standard force fields:

```python
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
system = ff.createSystem(pdb.topology)
```

Two ways in, depending on whether you have a structure:

```python
from forcefill import build_forcefield_xml, build_ligand_xml

# A complex: find what amber14 cannot match and parameterize it
build_forcefield_xml("complex.pdb", "extras.xml")

# Just the ligand, no structure anywhere
build_ligand_xml("benzamidinium.sdf", "ben.xml")
```

## What it does

1. **Identify** — `ForceField.getUnmatchedResidues` finds every residue the
   base force field (default: `amber14-all.xml` + `amber14/tip3p.xml`) has no
   template for.
2. **Classify** — only chemistry a stand-alone GAFF treatment is actually
   *valid* for gets parameterized; see
   [What gets skipped, and why](https://forcefill.readthedocs.io/en/latest/guide/what-gets-parameterized.html).
3. **Check** — before anything expensive: the net charge is read from the ligand
   file, a supplied file is confirmed to be the same molecule as the residue, and
   the geometry is checked for the faults that produce NaN energies.
4. **Parameterize** — each unique residue through `antechamber` (GAFF2 atom
   types, AM1-BCC charges → `.mol2`) and `parmchk2` (missing parameters →
   `.frcmod`), through OpenFF with `backend="smirnoff"`, through the Espaloma
   graph network with `backend="espaloma"`, or from a CGenFF stream file with
   `backend="charmm"`.
5. **Assemble** — one XML per residue, plus one combined XML.
6. **Validate** — an `openmm.System` is built from `base force field + new XML`
   for every parameterized residue on its own (and for the whole input when
   nothing was skipped), so a template that does not match its residue fails
   loudly here instead of at simulation time. `minimize=True` additionally
   catches parameters that are unphysical rather than merely absent (see
   [Checks](https://forcefill.readthedocs.io/en/latest/guide/checks.html)).

## Documentation

Full documentation is at **[forcefill.readthedocs.io](https://forcefill.readthedocs.io/en/latest/)**:

| | |
|---|---|
| [Installation](https://forcefill.readthedocs.io/en/latest/installation.html) | conda-forge, and the two dependencies that are not ordinary ones |
| [Quickstart](https://forcefill.readthedocs.io/en/latest/quickstart.html) | a complex in, an ffxml out, and simulating with it |
| [What gets skipped, and why](https://forcefill.readthedocs.io/en/latest/guide/what-gets-parameterized.html) | the residues forcefill refuses, and what to do instead |
| [Cleaning the structure](https://forcefill.readthedocs.io/en/latest/guide/cleaning.html) | water, buffer ions and crystallization additives — and the metals it keeps |
| [Ligand input](https://forcefill.readthedocs.io/en/latest/guide/ligands.html) | SDF, MOL2 or SMILES; `LigandSpec`; ligands with no structure at all |
| [Four backends](https://forcefill.readthedocs.io/en/latest/guide/backends.html) | GAFF, OpenFF Sage, Espaloma, CGenFF — and which mix |
| [CHARMM and CGenFF](https://forcefill.readthedocs.io/en/latest/guide/charmm.html) | converting a ParamChem stream file, and the three ways doing it by hand goes wrong |
| [Bespoke torsions](https://forcefill.readthedocs.io/en/latest/guide/bespokefit.html) | a BespokeFit OFFXML, and the checks that protect the QC you paid for |
| [Checks](https://forcefill.readthedocs.io/en/latest/guide/checks.html) | what runs before the expensive step, and what runs after |
| [Things to get right](https://forcefill.readthedocs.io/en/latest/guide/gotchas.html) | explicit hydrogens, net charge, virtual sites, the periodic box |
| [Relation to `openmmforcefields`](https://forcefill.readthedocs.io/en/latest/guide/openmmforcefields.html) | which template generators forcefill uses, and why not `GAFFTemplateGenerator` |
| [API reference](https://forcefill.readthedocs.io/en/latest/api/index.html) | every public function, class and constant |

## Installation

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

Requires Python ≥ 3.10. Two things are not ordinary dependencies: AmberTools is
not a Python package (the `antechamber` and `parmchk2` executables must be on
`PATH` for the `gaff` backend), and `espaloma` is optional because it pulls in
PyTorch. The full list, and the reason the `openmmforcefields` floor is 0.16, is
in [the installation guide](https://forcefill.readthedocs.io/en/latest/installation.html).
## Quickstart

```python
from forcefill import build_forcefield_xml

result = build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    net_charges={"LIG": -1},  # essential for sensible AM1-BCC charges
    clean_structure=True,  # drop water, buffer ions and crystallization additives
)
print(result.parameterized)  # ['LIG']
print(result.skipped)  # {'ZN': 'monatomic species - ...'}
```

then simulate with:

```python
from openmm import app

pdb = app.PDBFile("complex.pdb")
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
system = ff.createSystem(pdb.topology)
```

`result` also reports the per-residue XML files (`result.residue_xmls`), the
skip reasons (`result.skipped`), and the directory holding every intermediate
file (`result.workdir`) for inspection — pass `cleanup=True` to remove it on
success.

More, including per-ligand settings and the other three backends, in
[the guide](https://forcefill.readthedocs.io/en/latest/guide/index.html).

## Development

```bash
conda env create -f environment.yml && conda activate forcefill
pip install -e . --no-deps
pytest -m "not integration and not smirnoff"   # fast hermetic tests only
pytest                                         # everything, including real antechamber and OpenFF
```

Style is enforced by ruff (`pip install -e '.[dev]' && pre-commit install`).

The docs build is pip-only and separate from the conda environment — see
[Contributing](https://forcefill.readthedocs.io/en/latest/contributing.html).

## Roadmap

- A `forcefill` command-line interface — `build_ligand_xml` is the shape one
  wants.
- Caching, so re-runs into the same workdir skip finished antechamber jobs —
  and finished SMIRNOFF/Espaloma templates, which openmmforcefields' own
  `cache=` cannot help with on the code path forcefill uses.
- Covalently bound ligands. Still skipped, and deliberately: a stand-alone
  treatment of a polymer-linked residue is wrong whichever backend produces it.

## License

MIT — see [LICENSE](LICENSE).
