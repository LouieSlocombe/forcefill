# Examples

Five runnable scripts live in
[`examples/`](https://github.com/LouieSlocombe/forcefill/tree/main/examples),
all built around one system: benzamidine bound to trypsin, PDB entry
[3PTB](https://www.rcsb.org/structure/3PTB). The prepared structure and the
ligand files are committed, so every script runs from a checkout without a
download.

```bash
conda activate forcefill        # needs ambertools on PATH
python examples/parameterize_ligand.py
```

| Script | What it shows |
|---|---|
| [`parameterize_ligand.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/parameterize_ligand.py) | The main path: `build_forcefield_xml(minimize=True)` on the prepared complex, then loading `amber14 + ben_ff.xml` by hand and running dynamics |
| [`parameterize_ligand_standalone.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/parameterize_ligand_standalone.py) | {doc}`No structure at all <guide/ligands>`: the same ligand through GAFF2 and OpenFF Sage, then both merged into one XML. Never states the net charge — the SDF says `+1` and forcefill reads it |
| [`parameterize_ligand_charmm.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/parameterize_ligand_charmm.py) | {doc}`CHARMM <guide/charmm>` without AmberTools or a toppar download: a CGenFF stream file converted, loaded on top of `charmm36.xml` and minimized — then the two combinations OpenMM could never load, being refused by name |
| [`parameterize_ligand_bespoke.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/parameterize_ligand_bespoke.py) | {doc}`BespokeFit <guide/bespokefit>` end to end, including all three refusals, against a stand-in OFFXML synthesized from the installed release — so it needs no BespokeFit and no quantum chemistry |
| [`prepare_trypsin_ben.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/prepare_trypsin_ben.py) | Where the prepared structure came from: the preparation forcefill does *not* do for you |

## The part that is not forcefill's job

`prepare_trypsin_ben.py` is worth reading even if you never run it, because every
step in it is a decision you will face with your own systems and none of them are
automatic:

- **[PDBFixer](https://github.com/openmm/pdbfixer) cannot protonate your ligand.**
  It has hydrogen templates for standard residues only, so an X-ray structure
  comes back with a fully protonated protein and a still-bare ligand. This is the
  step people miss.
- **The ligand's ionization state is your call.** Benzamidine's amidine group
  (pKa ≈ 11.6) is protonated at physiological pH, which is exactly why it binds
  the S1 pocket — but the deposited formula (`C7 H8 N2`) is neutral. Simulating
  what was deposited would be simulating the wrong molecule.
- **Draw the bonds once, in a file.** The protonated molecule is written out as
  `benzamidinium.sdf` and passed in with `residue_files=`, so antechamber reads
  explicit bond orders instead of re-perceiving them from geometry — which is
  what gets the aromatic ring and the delocalized amidinium wrong.
- **A structural metal keeps its position and loses its CONECT records.** Force
  fields model ions nonbonded, and an ion with bonds can never match an ion
  template. This calcium is also why {func}`~forcefill.clean_pdb` keeps
  structural metals by default; see {doc}`guide/cleaning`.

The full walkthrough, including the file inventory and the provenance of the
CGenFF stream file, is in
[`examples/README.md`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/README.md).
