# Example: benzamidine bound to trypsin (PDB 3PTB)

A walkthrough of the part of ligand parameterization that tools cannot do for
you — structure preparation — followed by the part forcefill does: turning the
prepared complex into a working force field.

```bash
conda activate forcefill        # needs ambertools on PATH
python parameterize_ligand.py
```

`parameterize_ligand.py` runs `build_forcefield_xml(minimize=True)` on the
prepared structure and writes `ben_ff.xml`. The `minimize=True` proves the result
from inside the library — a finite energy that a minimizer can lower, for
benzamidinium alone and for the whole complex — and the script then loads
`amber14 + ben_ff.xml` by hand and runs a few steps of dynamics, which is how you
would actually use the file.

Two more scripts cover the other entry points:

```bash
python parameterize_ligand_standalone.py   # no structure: the ligand file is the input
python parameterize_ligand_charmm.py       # no AmberTools, no toppar download
```

The standalone script runs `build_ligand_xml` on the same
`data/benzamidinium.sdf` through both backends (GAFF2 and OpenFF Sage), then
through both at once merged into a single XML. It never states the net charge:
the SDF says `+1` and forcefill reads it, which is what separates it from the
script above.

The CHARMM script converts `data/benzamidinium_cgenff.str` into an ffxml, loads
it on top of `charmm36.xml`, minimizes benzamidinium at the crystal geometry, and
then shows the two combinations OpenMM could never load — a CHARMM ligand with an
Amber base force field, and a CHARMM ligand mixed with a GAFF one — being refused
by name. That refusal is the point: the two conventions scale 1-4 interactions
differently, so mixing them is an error, not a trade-off.

## Where the prepared structure came from

`data/trypsin_ben_prepared.pdb` and `data/benzamidinium.sdf` were produced by
`prepare_trypsin_ben.py` from the deposited 3PTB entry. Every step in it is a
decision you will face with your own systems:

1. **The protein and waters** get missing atoms and hydrogens from
   [PDBFixer](https://github.com/openmm/pdbfixer) at pH 7. PDBFixer *cannot*
   protonate the ligand — it has hydrogen templates only for standard residues —
   so benzamidine comes out still bare. This is the step people miss: an X-ray
   structure has no hydrogens anywhere, and `Modeller.addHydrogens` quietly
   fixes only the residues it knows.
2. **The ligand's ionization state is your call, not software's.** Benzamidine's
   amidine group (pKa ≈ 11.6) is protonated at physiological pH — which is
   precisely why it binds the S1 pocket of trypsin. The deposited formula
   (`C7 H8 N2`, neutral) is *not* what you should simulate. RDKit assigns bond
   orders from the SMILES template `NC(=[NH2+])c1ccccc1` and adds hydrogens with
   3D coordinates, giving benzamidinium (+1, 18 atoms).
3. **The ligand file carries the drawn bonds.** The protonated molecule is
   written both into the PDB (with CONECT records) and as
   `data/benzamidinium.sdf`. The example passes the SDF to
   `build_forcefield_xml(residue_files={"BEN": ...})`, so antechamber reads
   explicit bond orders instead of re-perceiving them from geometry — the
   aromatic ring and the delocalized amidinium are exactly what geometry-based
   perception gets wrong.
4. **The structural Ca²⁺ keeps its position, loses its CONECT records.** The
   deposited entry draws four coordination "bonds" from the calcium to protein
   oxygens. Force fields model ions nonbonded, and an ion with bonds can never
   match an ion template; amber14 then handles Ca²⁺ itself, so forcefill never
   sees it. (Had it stayed bonded, forcefill would have refused it as
   "covalently bonded to neighbouring residues", which is the correct answer for
   a mis-drawn ion.)

   This calcium is also why `clean_pdb` keeps structural metals by default.
   Running it on the prepared structure removes the 62 waters (3425 → 3239
   atoms) and leaves the Ca²⁺ in place, saying so in `result.retained`:

   ```python
   from forcefill import clean_pdb

   result = clean_pdb("data/trypsin_ben_prepared.pdb", "data/trypsin_ben_dry.pdb")
   result.removed  # {'HOH': ('water', 62)}
   result.retained  # {'CA': 'structural metal retained by default ...'}
   ```

   Pass `remove_structural_metals=True` if you really do want it gone — but for
   trypsin you do not: the calcium-binding loop needs it.

With that preparation, `amber14-all.xml + amber14/tip3p.xml` matches everything
except `BEN`, and the run reduces to:

```python
result = build_forcefield_xml(
    "data/trypsin_ben_prepared.pdb",
    "ben_ff.xml",
    net_charges={"BEN": 1},  # benzamidinium, not benzamidine
    residue_files={"BEN": "data/benzamidinium.sdf"},
)
```

Expected output: `parameterized: ['BEN']`, nothing skipped, validation OK (a
System for BEN alone and for the full complex), and a large negative energy after
minimization for both.

## Files

| File | Committed? | What it is |
|---|---|---|
| `prepare_trypsin_ben.py` | yes | One-time preparation script (downloads 3PTB, needs pdbfixer + rdkit) |
| `parameterize_ligand.py` | yes | The actual example: forcefill + OpenMM minimization and dynamics |
| `parameterize_ligand_standalone.py` | yes | The ligand on its own, no PDB: both backends, then both merged |
| `parameterize_ligand_charmm.py` | yes | The same ligand under CHARMM: converting a CGenFF stream file |
| `data/trypsin_ben_prepared.pdb` | yes | Prepared complex: protonated protein/waters/benzamidinium + bond-less Ca²⁺ |
| `data/benzamidinium.sdf` | yes | The ligand as drawn: bond orders, +1 charge, 3D hydrogens |
| `data/benzamidinium_cgenff.str` | yes | The same ligand's CGenFF topology, in the format ParamChem returns |
| `data/3PTB.pdb` | no (downloaded) | The deposited entry, fetched by the prep script |
| `ben_ff.xml`, `wd/` | no (generated) | Output of the example run |
| `ben_standalone_*.xml`, `wd_standalone_*/` | no (generated) | Output of the standalone run |
| `ben_charmm.xml`, `wd_charmm/` | no (generated) | Output of the CHARMM run |

`data/benzamidinium_cgenff.str` is not a `cgenff` run — the program is licensed
and its output is not redistributable. Its atom types and charges are taken
verbatim from `RESI BAMI` in the published CGenFF topology file, renamed and
reordered to match `data/benzamidinium.sdf`. For your own ligand, use what
[ParamChem](https://cgenff.paramchem.org) gives you.
