# forcefill

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
   {doc}`guide/what-gets-parameterized`.
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
{doc}`guide/checks`).

## Where to next

- {doc}`installation` — conda-forge is the supported route; AmberTools is not
  a Python package.
- {doc}`quickstart` — a complex in, an ffxml out, and simulating with it.
- {doc}`guide/index` — the four backends, CHARMM/CGenFF, BespokeFit, cleaning,
  and what forcefill refuses to parameterize.
- {doc}`api/index` — every public function, class and constant.


```{toctree}
:hidden:
:maxdepth: 2

installation
quickstart
guide/index
examples
api/index
contributing
changelog
```
