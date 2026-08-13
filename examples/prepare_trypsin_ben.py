"""Prepare the benzamidine-trypsin example structure (PDB 3PTB) from scratch.

This is the one-time preparation that produced the committed files in
``examples/data/``; run it only if you want to regenerate them. It documents
every decision the tutorial depends on:

1. Download 3PTB (bovine trypsin + benzamidine, CC0 data from the RCSB).
2. PDBFixer adds the missing protein atoms and all protein/water hydrogens at
   pH 7. PDBFixer cannot protonate the ligand: BEN has no hydrogen template,
   so it comes out of this step still bare.
3. RDKit protonates benzamidine: bond orders are assigned from the SMILES
   template ``NC(=[NH2+])c1ccccc1`` (benzamidinium, the +1 species at pH 7 -
   its amidine pKa is ~11.6), then hydrogens are added with 3D coordinates.
   The result is written both into the prepared PDB and as an SDF, so the
   parameterization example can hand antechamber the drawn bond orders via
   ``residue_files`` instead of letting it re-perceive them from geometry.
4. The protonated ligand replaces the bare one and the merged structure is
   written with CONECT records.

Requires: pdbfixer, rdkit (both on conda-forge).
"""

import logging
import urllib.request
from pathlib import Path

from openmm import Vec3, app, unit
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("prepare_trypsin_ben")

HERE = Path(__file__).parent
DATA = HERE / "data"
RAW_PDB = DATA / "3PTB.pdb"  # not committed; downloaded on demand
PREPARED_PDB = DATA / "trypsin_ben_prepared.pdb"
BEN_SDF = DATA / "benzamidinium.sdf"

BENZAMIDINIUM_SMILES = "NC(=[NH2+])c1ccccc1"


def download_3ptb() -> None:
    if RAW_PDB.is_file():
        log.info("Using existing %s", RAW_PDB)
        return
    url = "https://files.rcsb.org/download/3PTB.pdb"
    log.info("Downloading %s", url)
    urllib.request.urlretrieve(url, RAW_PDB)


def fix_protein() -> PDBFixer:
    """Repair the protein and protonate protein + waters (not the ligand)."""
    fixer = PDBFixer(filename=str(RAW_PDB))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # 3PTB has no internal gaps; do not model termini
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def protonate_benzamidine() -> Chem.Mol:
    """RDKit: bare BEN (heavy atoms + CONECT from the deposited PDB) -> benzamidinium."""
    lines = RAW_PDB.read_text().splitlines()
    ben_serials = set()
    block = []
    for line in lines:
        if line.startswith("HETATM") and line[17:20] == "BEN":
            ben_serials.add(int(line[6:11]))
            block.append(line)
    # The deposited CONECT records carry the ligand's bond graph.
    block += [line for line in lines if line.startswith("CONECT") and int(line[6:11]) in ben_serials]
    block.append("END")

    mol = Chem.MolFromPDBBlock("\n".join(block), removeHs=False, proximityBonding=False)
    template = Chem.MolFromSmiles(BENZAMIDINIUM_SMILES)
    mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    mol_h = Chem.AddHs(mol, addCoords=True)
    Chem.MolToMolFile(mol_h, str(BEN_SDF))
    log.info("Wrote %s (%d atoms, net charge %+d)", BEN_SDF, mol_h.GetNumAtoms(), Chem.GetFormalCharge(mol_h))
    return mol_h


def rdkit_to_openmm(mol: Chem.Mol) -> tuple[app.Topology, unit.Quantity]:
    """Convert the protonated ligand into an OpenMM topology + positions."""
    top = app.Topology()
    chain = top.addChain("L")
    res = top.addResidue("BEN", chain)
    omm_atoms = []
    n_h = 0
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is not None and info.GetName().strip():
            name = info.GetName().strip()
        else:  # RDKit-added hydrogen
            n_h += 1
            name = f"H{n_h}"
        element = app.element.Element.getByAtomicNumber(atom.GetAtomicNum())
        omm_atoms.append(top.addAtom(name, element, res))
    for bond in mol.GetBonds():
        top.addBond(omm_atoms[bond.GetBeginAtomIdx()], omm_atoms[bond.GetEndAtomIdx()])
    conf = mol.GetConformer()
    # Positions must be Vec3, not plain tuples: openmm's unit math scales each
    # element, and a tuple silently breaks that (Vec3 supports scalar multiply).
    positions = unit.Quantity(
        [Vec3(p.x, p.y, p.z) for p in (conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms()))],
        unit.angstrom,
    ).in_units_of(unit.nanometer)
    return top, positions


def main() -> None:
    DATA.mkdir(exist_ok=True)
    download_3ptb()
    fixer = fix_protein()
    mol_h = protonate_benzamidine()

    modeller = app.Modeller(fixer.topology, fixer.positions)

    # The deposited CONECT records give the structural Ca2+ four coordination
    # "bonds" to protein oxygens. Force fields model the ion nonbonded (and a
    # bonded ion could never match an ion template), so re-add it bond-less.
    ca_residues = [r for r in modeller.topology.residues() if r.name == "CA"]
    ca_positions = [modeller.positions[next(r.atoms()).index] for r in ca_residues]
    modeller.delete(ca_residues)
    for ca_pos in ca_positions:
        ion_top = app.Topology()
        ion_res = ion_top.addResidue("CA", ion_top.addChain("I"))
        ion_top.addAtom("CA", app.element.calcium, ion_res)
        modeller.add(ion_top, unit.Quantity([ca_pos.value_in_unit(unit.nanometer)], unit.nanometer))

    # Swap the bare BEN for the protonated one.
    bare_ben = [r for r in modeller.topology.residues() if r.name == "BEN"]
    modeller.delete(bare_ben)
    ben_top, ben_pos = rdkit_to_openmm(mol_h)
    modeller.add(ben_top, ben_pos)

    with open(PREPARED_PDB, "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    log.info("Wrote %s", PREPARED_PDB)

    # Sanity check: against amber14, only BEN should now be unmatched.
    from forcefill import find_nonstandard_residues

    unmatched = find_nonstandard_residues(app.PDBFile(str(PREPARED_PDB)).topology)
    log.info("Unmatched residues in the prepared structure: %s", sorted({r.name for r in unmatched}))


if __name__ == "__main__":
    main()
