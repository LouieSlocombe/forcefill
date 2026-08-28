# Cleaning the structure first

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

The split between the two ion rows is the point. Bulk ions come from the buffer
or from neutralizing the box — they occupy no defined site and you re-add them
with `Modeller.addSolvent` anyway. Structural metals are buried, directionally
coordinated and often required for the fold or the chemistry: trypsin's Ca²⁺
(3PTB residue `CA` 480) rigidifies the calcium-binding loop. Deleting a needed
metal is silent and wrong; keeping an unwanted one shows up in `retained` and
goes away with `remove_structural_metals=True`.

For the long tail, `keep=("IMD",)` protects an additive that is really your
ligand and `extra_remove=("HEM",)` drops something the tables leave alone;
`extra_remove` refuses standard residue names, so a typo cannot shred a protein.
The tables themselves are importable (`forcefill.ADDITIVE_RESIDUES`,
`STRUCTURAL_METAL_RESIDUES`, …).

The cleaner is **subtractive only**: it never adds missing atoms, models loops,
protonates, selects chains or strips hydrogens — that is
[PDBFixer's](https://github.com/openmm/pdbfixer) job, so clean *after* you
repair, not instead. It also refuses to delete a residue covalently bonded to a
neighbour: `Modeller` drops the bonds along with the atoms and never says so,
leaving the partner with an unsatisfied valence. And it does not by itself make
the full-structure checks run — on a *raw* crystal structure the protein is
unmatched too, every residue missing its hydrogens, so `skipped` stays
non-empty.

