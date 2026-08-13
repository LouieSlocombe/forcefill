# forcefill

[![ci](https://github.com/LouieSlocombe/forcefill/actions/workflows/ci.yml/badge.svg)](https://github.com/LouieSlocombe/forcefill/actions/workflows/ci.yml)

Identify the non-standard residues (ligands, cofactors, hetero molecules) in a
PDB file and turn them into a ready-to-use [OpenMM](https://openmm.org)
force-field XML, using AmberTools (`antechamber` GAFF/GAFF2 atom types +
AM1-BCC charges, `parmchk2` for missing parameters) and
[ParmEd](https://github.com/ParmEd/ParmEd) for the assembly.

The output is a plain ffxml file you load alongside the standard force fields:

```python
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
system = ff.createSystem(pdb.topology)
```

## What it does

1. **Identify** — `ForceField.getUnmatchedResidues` finds every residue the
   base force field (default: `amber14-all.xml` + `amber14/tip3p.xml`) has no
   template for.
2. **Classify** — unmatched residues are triaged; only chemistry that a
   stand-alone GAFF treatment is actually *valid* for gets parameterized
   (see the table below).
3. **Parameterize** — each unique residue is written to its own PDB and run
   through `antechamber` (GAFF2 atom types, AM1-BCC charges → `.mol2`) and
   `parmchk2` (missing parameters → `.frcmod`).
4. **Assemble** — ParmEd merges the GAFF database, the frcmods and the mol2
   templates into one XML per residue plus one combined XML.
5. **Validate** — an `openmm.System` is built from `base force field + new
   XML` for every parameterized residue on its own (and for the whole input
   when nothing was skipped), so a template that does not match its residue
   fails loudly here instead of at simulation time.

### What gets skipped, and why

Most of the value is in what forcefill *refuses* to parameterize. Running
antechamber on the wrong kind of residue produces plausible-looking but
physically wrong parameters, so these are reported and skipped:

| Unmatched residue | Action | Do this instead |
|---|---|---|
| Standard residue (e.g. `ALA`, `HOH`) that failed to match | skip | It is missing atoms or has non-standard atom names — repair the structure with [PDBFixer](https://github.com/openmm/pdbfixer) or `Modeller.addHydrogens` |
| Monatomic species (ions such as `ZN`, `NA`) | skip | Load an ion parameter file; GAFF/antechamber cannot treat bare ions |
| Residue covalently bonded to its neighbours (modified amino acids, glycans) | skip | Stand-alone GAFF is not valid for polymer-linked residues; cap the fragment and derive charges consistently with the backbone force field (pyRED- or ffparam-style workflows) |
| Free-standing hetero molecule (ligand, cofactor) | **parameterize** | — |

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

Requires Python ≥ 3.10, `openmm >= 7.6`, `parmed >= 3.4`, and the AmberTools
executables (`antechamber`, `parmchk2`) on `PATH` at run time.

## Quickstart

```python
from forcefill import build_forcefield_xml

result = build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    net_charges={"LIG": -1},  # essential for sensible AM1-BCC charges
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
file (`result.workdir`) for inspection.

### Things to get right

- **Explicit hydrogens.** Ligands must contain all hydrogens with reasonable
  geometry; AM1-BCC charges are meaningless otherwise. forcefill warns when a
  ligand has none.
- **Net charge.** Pass `net_charges={"RES": q}` for every charged ligand.
  A wrong net charge is the classic source of plausible-but-wrong charges;
  forcefill warns about `net_charges` keys that match no residue (typos,
  case mismatches).
- **Connectivity.** Element columns and (for hetero groups) CONECT records
  should be present in the PDB.
- **One XML at a time.** Load either the combined XML *or* the per-residue
  XMLs, never both — the duplicated GAFF atom-type definitions would collide.

## Relation to `openmmforcefields`

If you would rather not manage XML files at all,
[`openmmforcefields.generators.GAFFTemplateGenerator`](https://github.com/openmm/openmmforcefields)
does the same antechamber job on the fly at `createSystem` time. forcefill is
for when you want the opposite trade-off: explicit, inspectable, versionable
XML artifacts, produced once, with the skip-classification above telling you
which residues need a different treatment entirely.

## Development

```bash
conda env create -f environment.yml && conda activate forcefill
pip install -e . --no-deps
pytest -m "not integration"   # fast hermetic tests, no AmberTools needed
pytest                        # includes end-to-end antechamber runs
```

Style is enforced by ruff (`pip install -e '.[dev]' && pre-commit install`).

## Roadmap

- Direct ligand-file input (SDF/MOL2) so antechamber sees the drawn bond
  orders instead of re-perceiving them from PDB geometry.
- A worked example under `examples/`.
- A `forcefill` command-line interface.
- Caching, so re-runs into the same workdir skip finished antechamber jobs.

## License

MIT — see [LICENSE](LICENSE).
