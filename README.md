# forcefill

[![ci](https://github.com/LouieSlocombe/forcefill/actions/workflows/ci.yml/badge.svg)](https://github.com/LouieSlocombe/forcefill/actions/workflows/ci.yml)

Turn ligands into a ready-to-use [OpenMM](https://openmm.org) force-field XML —
either the non-standard residues (ligands, cofactors, hetero molecules) found in
a PDB, or ligand files on their own. Parameters come from AmberTools
(`antechamber` GAFF/GAFF2 atom types + AM1-BCC charges, `parmchk2` for missing
parameters) via [ParmEd](https://github.com/ParmEd/ParmEd), or from
[OpenFF](https://openforcefield.org) Sage via
[openmmforcefields](https://github.com/openmm/openmmforcefields).

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
2. **Classify** — unmatched residues are triaged; only chemistry that a
   stand-alone GAFF treatment is actually *valid* for gets parameterized
   (see the table below).
3. **Check** — before anything expensive: the net charge is read from the ligand
   file, a supplied file is confirmed to be the same molecule as the residue, and
   the geometry is checked for the faults that produce NaN energies.
4. **Parameterize** — each unique residue goes through `antechamber` (GAFF2 atom
   types, AM1-BCC charges → `.mol2`) and `parmchk2` (missing parameters →
   `.frcmod`), or through OpenFF with `backend="smirnoff"`.
5. **Assemble** — ParmEd merges the GAFF database, the frcmods and the mol2
   templates into one XML per residue plus one combined XML.
6. **Validate** — an `openmm.System` is built from `base force field + new
   XML` for every parameterized residue on its own (and for the whole input
   when nothing was skipped), so a template that does not match its residue
   fails loudly here instead of at simulation time. With `minimize=True` each
   of those is also energy-minimized, which catches unphysical parameters
   that build a perfectly valid System (see below).

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

Note the last row cuts both ways: a glycerol or a sulfate left over from the
crystallization drop *is* a free-standing hetero molecule, so it gets
parameterized too. Strip those first — see below.

### Cleaning the structure first

A structure straight from the PDB carries water, buffer ions and whatever was
in the crystallization drop. `clean_pdb` removes them:

```python
from forcefill import clean_pdb

result = clean_pdb("3ptb.pdb", "3ptb_clean.pdb")
print(result.removed)  # {'HOH': ('water', 62)}
print(result.retained)  # {'CA': 'structural metal retained by default ...'}
```

or in memory, as part of the pipeline:

```python
result = build_forcefield_xml("complex.pdb", "extras.xml", clean_structure=True)
print(result.cleaning.n_atoms_removed)
```

| Category | Examples | Default |
|---|---|---|
| Water | `HOH`, `WAT`, `SOL`, `DOD` | **removed** |
| Bulk counter-ions | `NA`, `CL`, `K`, `BR`, `IOD` | **removed** |
| Crystallization additives | `GOL`, `EDO`, `PEG`, `DMS`, `SO4`, `EPE`, `BME` | **removed** |
| Structural metals | `CA`, `ZN`, `MG`, `MN`, `FE`, `CU` | **kept**, and reported |

The split between the last two ion rows is the point. Bulk ions come from the
buffer or from neutralizing the box — they occupy no defined site and you re-add
them with `Modeller.addSolvent` anyway. Structural metals are buried,
directionally coordinated, and often required for the fold or the chemistry:
trypsin's Ca²⁺ (3PTB residue `CA` 480) rigidifies the calcium-binding loop.
Deleting a needed metal is silent and wrong; keeping an unwanted one shows up in
`retained` and in the log, and goes away with
`remove_structural_metals=True`. Inspect the tables directly
(`forcefill.ADDITIVE_RESIDUES`, `STRUCTURAL_METAL_RESIDUES`, …) before you run.

Two escape hatches for the long tail, since these lists always drift:
`keep=("IMD",)` protects an additive that is really your ligand, and
`extra_remove=("HEM",)` drops something the tables leave alone. `extra_remove`
refuses standard residue names, so a typo cannot shred a protein.

**What it will not do.** The cleaner is subtractive only. It never adds missing
atoms, models missing loops, protonates anything, selects chains or strips
hydrogens — repairing a structure is
[PDBFixer's](https://github.com/openmm/pdbfixer) job and forcefill does not
duplicate it. Clean *after* you repair, not instead. It also refuses to delete
any residue that is covalently bonded to a neighbour: OpenMM's `Modeller` drops
the bonds along with the atoms and never says so, leaving the partner with an
unsatisfied valence.

One honest caveat: cleaning does not by itself make the full-structure checks
run. On a *raw* crystal structure the protein is unmatched too — every residue
is missing its hydrogens — so `skipped` stays non-empty either way. Cleaning
removes one common reason the check is suppressed (un-hydrogenated
crystallographic water); the rest is PDBFixer's job.

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

Requires Python ≥ 3.10 and, at import time, `openmm >= 7.6`, `parmed >= 3.4`,
`rdkit`, `openff-toolkit >= 0.16` and `openmmforcefields >= 0.14`. All five are
ordinary dependencies — there are no extras to pick and nothing is imported
lazily.

AmberTools is the exception, because it is not a Python package: the
`antechamber` and `parmchk2` executables must be on `PATH` at run time for the
`gaff` backend. `backend="smirnoff"` does not need them.

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
file (`result.workdir`) for inspection (pass `cleanup=True` to remove it on
success).

### Checking the parameters, not just the templates

Building a `System` proves the templates match and no parameter is missing. It
says nothing about whether the numbers are physical: a NaN charge or a zero
force constant survives it and only shows up later as an exploding simulation.
`minimize=True` adds an energy evaluation and a short minimization — of each
parameterized residue in vacuum, and of the whole input when nothing was
skipped — and raises if the potential energy is not finite at either end:

```python
result = build_forcefield_xml("complex.pdb", "extras.xml", minimize=True)
lig = result.minimizations["LIG"]
print(f"{lig.initial_energy:.0f} -> {lig.final_energy:.0f} kJ/mol")
print(result.full_minimization.max_force)  # kJ/mol/nm
```

`max_force` and `energy_change` are reported for inspection, not enforced —
what counts as converged depends on the system. The same check is available
on its own as `minimize_with_forcefield_xml(topology, positions, xml)`, which
takes the OpenMM knobs (`nonbonded_method`, `max_iterations`, `platform_name`)
that the pipeline leaves at their defaults.

### Supplying the ligand as drawn (SDF, MOL2 or SMILES)

Extracting a ligand from a PDB forces antechamber to re-perceive bond orders
from geometry — a classic source of silently wrong atom types for aromatics
and charged groups. If you have the ligand as an SDF or MOL2 with explicit
bonds and protonation, pass it directly:

```python
result = build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    residue_files={"LIG": "lig.sdf"},  # used instead of PDB extraction
)
```

A SMILES works too. When the residue is also in the structure,
the coordinates stay as deposited and only the bond orders come from the
SMILES — the crystal geometry is better than anything embedding produces:

```python
from forcefill import LigandSpec

result = build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    ligands={"BEN": LigandSpec(smiles="NC(=[NH2+])c1ccccc1")},
)
```

Either way the molecule must be the same one as the residue in the PDB,
hydrogens included — the generated template is matched against the PDB's bond
graph. forcefill checks that before running anything expensive; see below.

### Per-ligand settings

`LigandSpec` carries everything about one ligand: where it comes from and how to
treat it. Anything it does not set falls back to the call-level default, so a
spec states only what it overrides.

```python
result = build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    atom_type="gaff2",  # the default for everything...
    ligands={
        "BEN": LigandSpec(file="ben.sdf", backend="smirnoff"),
        "ATP": LigandSpec(file="atp.mol2", atom_type="gaff", net_charge=-4),
        "GOL": LigandSpec(smiles="OCC(O)CO", antechamber_args=("-dr", "no")),
    },
)
```

The older `net_charges`, `multiplicities` and `residue_files` mappings still
work and are folded in. Setting the same thing both ways raises rather than
silently picking a winner.

### Two backends

| | `backend="gaff"` (default) | `backend="smirnoff"` |
|---|---|---|
| Parameters | GAFF/GAFF2 atom types, AM1-BCC charges | OpenFF Sage, SMARTS-matched |
| Needs | AmberTools on `PATH` | nothing beyond the install |
| Ligand source | PDB residue, SDF, MOL2 or SMILES | **SDF, MOL2 or SMILES only** |

SMIRNOFF assigns parameters by matching SMARTS against the chemical graph, so it
has no way to work from a PDB residue — a PDB records no bond orders. A ligand on
that backend needs a `file` or a `smiles`, and says so if it has neither.

Backends can be mixed in one call: forcefill writes one combined XML and OpenMM
loads it. That works because SMIRNOFF names its atom types by a hash of the
molecule, so nothing collides — and because the merge keeps the two
`<PeriodicTorsionForce>` sections apart, since GAFF and SMIRNOFF impropers use
different `ordering` conventions.

### Ligands without a structure

`build_ligand_xml` is the same pipeline with the ligand as the whole input:

```python
from forcefill import build_ligand_xml, LigandSpec

build_ligand_xml("benzamidinium.sdf", "ben.xml")  # name from the file: BEN
build_ligand_xml(["a.sdf", "b.sdf"], "ligs.xml")  # several at once
build_ligand_xml({"LIG": LigandSpec(smiles="CO")}, "l.xml")  # named explicitly
```

Residue names not given explicitly come from the file name
(`benzamidinium.sdf` → `BEN`); a bare string is always a path, never a SMILES.
Validation still runs — it proves the template covers every atom and that no
parameter is missing — but with no structure there is no bond graph to match
against, so the molecule supplies its own topology. If you *do* have the complex,
`build_forcefield_xml` checks the thing that actually matters.

`examples/parameterize_ligand_standalone.py` runs this through both backends.

### Checks that run before the expensive step

antechamber's AM1-BCC can take an hour on a drug-sized ligand. Three mistakes
that used to cost that hour — or worse, silently produce wrong numbers — are now
caught in the first second:

- **Net charge read from the file.** An SDF or MOL2 states its own formal
  charge, so `net_charge` no longer defaults to a silent 0 for those. A supplied
  value that contradicts the file raises rather than picking one.
- **The ligand file must be the residue in the PDB.** A mismatch used to appear
  only at the end, as an opaque "no template matched". Now it names the
  difference: `ben.sdf has C7H9N2 (18 atoms), residue BEN has C7H8N2 (17)`.
- **Geometry sanity.** Coincident atoms, non-finite coordinates or a molecule
  written with no conformer — the standard causes of the NaN energies that
  `minimize=True` otherwise only finds at the very end.

`strict=False` downgrades the last two to warnings; these are heuristics and the
long tail is real.

**One behaviour change to know about.** If you already call
`build_forcefield_xml(..., residue_files={"BEN": "ben.sdf"})` without
`net_charges`, forcefill used to send `-nc 0` to antechamber. It now sends the
charge the SDF states, and logs that it did. For a charged ligand that changes
the output — to the correct answer. It applies only to ligands supplied as a
file or SMILES, never to a residue extracted from a PDB, and an explicit
`net_charge` always wins.

### Things to get right

- **Explicit hydrogens.** Ligands must contain all hydrogens with reasonable
  geometry; AM1-BCC charges are meaningless otherwise. forcefill warns when a
  ligand has none.
- **Net charge.** Read automatically from a supplied SDF/MOL2/SMILES. For a
  ligand extracted straight from a PDB there is nothing to read it from, so pass
  `net_charges={"RES": q}` yourself — a wrong net charge is the classic source
  of plausible-but-wrong charges. forcefill warns about keys that match no
  residue (typos, case mismatches).
- **Connectivity.** Element columns and (for hetero groups) CONECT records
  should be present in the PDB.
- **One XML at a time.** Load either the combined XML *or* the per-residue
  XMLs, never both — the duplicated GAFF atom-type definitions would collide.
- **Cleaning changes what the checks describe.** With `clean_structure=True`
  the full-structure `validate`/`minimize` results, and so
  `full_minimization.n_atoms`, refer to the *cleaned* topology, not the file on
  disk. Reconcile against `cleaning.n_atoms_after`.
- **The periodic box survives the strip.** A de-solvated structure keeps the
  box vectors of the solvated one, so a later PME run would use a
  mostly-empty box. Reset them yourself if you intend to simulate the cleaned
  structure directly.

## Relation to `openmmforcefields`

If you would rather not manage XML files at all,
[`openmmforcefields`](https://github.com/openmm/openmmforcefields)'
`GAFFTemplateGenerator` and `SMIRNOFFTemplateGenerator` do the same
parameterization on the fly at `createSystem` time — and forcefill uses the
latter for its own `smirnoff` backend. forcefill is for when you want the
opposite trade-off: explicit, inspectable, versionable XML artifacts, produced
once, with the skip-classification and the preflight checks above telling you
which residues need a different treatment entirely.

## Development

```bash
conda env create -f environment.yml && conda activate forcefill
pip install -e . --no-deps
pytest -m "not integration and not smirnoff"   # fast hermetic tests only
pytest                                         # everything, including real antechamber and OpenFF
```

Style is enforced by ruff (`pip install -e '.[dev]' && pre-commit install`).

## Roadmap

- A `forcefill` command-line interface — `build_ligand_xml` is the shape one
  wants.
- Caching, so re-runs into the same workdir skip finished antechamber jobs.
- Covalently bound ligands. Still skipped, and deliberately: a stand-alone
  treatment of a polymer-linked residue is wrong whichever backend produces it.

## License

MIT — see [LICENSE](LICENSE).
