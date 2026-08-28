# Ligand input

## Supplying the ligand as drawn (SDF, MOL2 or SMILES)

Extracting a ligand from a PDB forces antechamber to re-perceive bond orders
from geometry — a classic source of silently wrong atom types for aromatics and
charged groups. If you have the ligand as an SDF or MOL2 with explicit bonds and
protonation, pass it directly:

```python
result = build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    residue_files={"LIG": "lig.sdf"},  # used instead of PDB extraction
)
```

A SMILES works too. When the residue is also in the structure, the coordinates
stay as deposited and only the bond orders come from the SMILES — the crystal
geometry is better than anything embedding produces:

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
graph. forcefill checks that before running anything expensive; see {doc}`checks`.

## Per-ligand settings

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

## Ligands without a structure

`build_ligand_xml` is the same pipeline with the ligand as the whole input:

```python
from forcefill import build_ligand_xml, LigandSpec

build_ligand_xml("benzamidinium.sdf", "ben.xml")  # name from the file: BEN
build_ligand_xml(["a.sdf", "b.sdf"], "ligs.xml")  # several at once
build_ligand_xml({"LIG": LigandSpec(smiles="CO")}, "l.xml")  # named explicitly
```

Residue names not given explicitly come from the file name
(`benzamidinium.sdf` → `BEN`); a bare string is always a path, never a SMILES.
Validation still runs, but with no structure there is no bond graph to match
against — the molecule supplies its own topology. If you *do* have the complex,
`build_forcefield_xml` checks the thing that actually matters.

[`examples/parameterize_ligand_standalone.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/parameterize_ligand_standalone.py) runs this
through both backends.

