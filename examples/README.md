# Example: benzamidine bound to trypsin (PDB 3PTB)

A complete, honest walkthrough of the part of ligand parameterization that
tools cannot do for you — structure preparation — followed by the part
forcefill does: turning the prepared complex into a working force field.

```bash
conda activate forcefill        # needs ambertools on PATH
python parameterize_ligand.py
```

`parameterize_ligand.py` runs `build_forcefield_xml` on the prepared
structure, writes `ben_ff.xml`, then proves the result by building an
`openmm.System` with `amber14 + ben_ff.xml` and minimizing the complex.

## Where the prepared structure came from

`data/trypsin_ben_prepared.pdb` and `data/benzamidinium.sdf` were produced by
`prepare_trypsin_ben.py` from the deposited 3PTB entry. The script is short,
and every step in it is a decision you will face with your own systems:

1. **The protein and waters** get missing atoms and hydrogens from
   [PDBFixer](https://github.com/openmm/pdbfixer) at pH 7. PDBFixer *cannot*
   protonate the ligand — it has hydrogen templates only for standard
   residues — so benzamidine comes out of this step still bare. This is the
   step people miss: an X-ray structure has no hydrogens anywhere, and
   `Modeller.addHydrogens` quietly fixes only the residues it knows.
2. **The ligand's ionization state is your call, not software's.**
   Benzamidine's amidine group (pKa ≈ 11.6) is protonated at physiological
   pH — that is precisely why it binds the S1 pocket of trypsin. The deposited
   formula (`C7 H8 N2`, neutral) is *not* what you should simulate. RDKit
   assigns bond orders from the SMILES template `NC(=[NH2+])c1ccccc1` and adds
   hydrogens with 3D coordinates, giving benzamidinium (+1, 18 atoms).
3. **The ligand file carries the drawn bonds.** The protonated molecule is
   written both into the PDB (with CONECT records) and as
   `data/benzamidinium.sdf`. The example passes the SDF to
   `build_forcefield_xml(residue_files={"BEN": ...})`, so antechamber reads
   explicit bond orders instead of re-perceiving them from geometry — the
   aromatic ring and the delocalized amidinium are exactly the kind of
   chemistry geometry-based perception gets wrong.
4. **The structural Ca²⁺ keeps its position, loses its CONECT records.** The
   deposited entry draws four coordination "bonds" from the calcium to protein
   oxygens. Force fields model ions nonbonded; an ion with bonds can never
   match an ion template. amber14 then handles Ca²⁺ itself — forcefill never
   sees it. (Had it stayed bonded, forcefill would have refused it as
   "covalently bonded to neighbouring residues", which is the correct answer
   for a mis-drawn ion.)

With that preparation, `amber14-all.xml + amber14/tip3p.xml` matches
everything except `BEN`, and the run reduces to:

```python
result = build_forcefield_xml(
    "data/trypsin_ben_prepared.pdb",
    "ben_ff.xml",
    net_charges={"BEN": 1},  # benzamidinium, not benzamidine
    residue_files={"BEN": "data/benzamidinium.sdf"},
)
```

Expected output: `parameterized: ['BEN']`, nothing skipped, validation OK
(forcefill builds a System for BEN alone and for the full complex), and a
large negative energy after minimization.

## Files

| File | Committed? | What it is |
|---|---|---|
| `prepare_trypsin_ben.py` | yes | One-time preparation script (downloads 3PTB, needs pdbfixer + rdkit) |
| `parameterize_ligand.py` | yes | The actual example: forcefill + OpenMM minimization |
| `data/trypsin_ben_prepared.pdb` | yes | Prepared complex: protonated protein/waters/benzamidinium + bond-less Ca²⁺ |
| `data/benzamidinium.sdf` | yes | The ligand as drawn: bond orders, +1 charge, 3D hydrogens |
| `data/3PTB.pdb` | no (downloaded) | The deposited entry, fetched by the prep script |
| `ben_ff.xml`, `wd/` | no (generated) | Output of the example run |
