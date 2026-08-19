"""Turn a CGenFF stream file into an OpenMM force field, and prove it works.

The CHARMM counterpart of parameterize_ligand_standalone.py. The difference is
where the parameters come from: there forcefill derives them, here it converts
parameters you already have. There is no CHARMM equivalent of antechamber to
call - CGenFF parameters come from https://cgenff.paramchem.org or the licensed
cgenff program, both of which hand you a .str file.

Worth watching in the output:

  * the net charge is never stated. The RESI block's charges sum to +1 and
    forcefill reads that, as it reads an SDF's M CHG record.
  * the generated XML contains a residue template and nothing else. charmm36.xml
    already carries all 412 CGenFF atom types and their parameters, so forcefill
    names them rather than redefining them - which is not merely redundant to
    do: ParmEd writes epsilon=0 for types it only knows the mass of, and OpenMM
    lets a later file override an atom type without a word, so the real
    Lennard-Jones parameters would quietly disappear.
  * the two combinations OpenMM could never load are refused by name, up front.

Run from this directory:

    python parameterize_ligand_charmm.py

Needs nothing beyond an ordinary forcefill install - no AmberTools, no CHARMM
toppar download. data/benzamidinium_cgenff.str and OpenMM's own charmm36.xml
are the whole input.
"""

import logging
import xml.etree.ElementTree as ET
from pathlib import Path

import openmm
from openmm import app, unit
from rdkit import Chem

from forcefill import CHARMM_BASE_FORCEFIELD, LigandSpec, build_ligand_xml

HERE = Path(__file__).parent
BEN_STR = HERE / "data" / "benzamidinium_cgenff.str"
BEN_SDF = HERE / "data" / "benzamidinium.sdf"


def show_refusals():
    """The two things a CHARMM user hits first, and what forcefill says about them."""
    print("\n--- what is refused, and why ---")
    attempts = {
        "the Amber base force field": {"base_forcefield": ("amber14-all.xml", "amber14/tip3p.xml")},
        "a gaff ligand in the same call": {
            "base_forcefield": CHARMM_BASE_FORCEFIELD,
            "ligands": {
                "BEN": LigandSpec(backend="charmm", charmm_files=[BEN_STR]),
                "MOL": LigandSpec(backend="gaff", smiles="CO"),
            },
        },
    }
    for label, kwargs in attempts.items():
        ligands = kwargs.pop("ligands", {"BEN": LigandSpec(backend="charmm", charmm_files=[BEN_STR])})
        try:
            build_ligand_xml(ligands, HERE / "_never_written.xml", **kwargs)
        except ValueError as exc:
            print(f"{label}:\n  {str(exc).splitlines()[0]}")


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    print("--- converting the stream file ---")
    result = build_ligand_xml(
        {"BEN": LigandSpec(backend="charmm", charmm_files=[BEN_STR])},
        HERE / "ben_charmm.xml",
        base_forcefield=CHARMM_BASE_FORCEFIELD,
        workdir=HERE / "wd_charmm",
    )
    print(f"parameterized: {result.parameterized}")
    print(f"combined XML:  {result.forcefield_xml}")

    # What came out: a residue template, and only what the stream file itself
    # added on top of charmm36. For benzamidinium that is nothing - every term
    # it needs is in the released CGenFF set already.
    root = ET.parse(result.forcefield_xml).getroot()
    print(f"sections:      {[section.tag for section in root]}")
    print(f"atom types redefined: {len(root.findall('./AtomTypes/Type'))}")

    # An ordinary force-field file, loaded next to the standard CHARMM set - and
    # meaningless without it, since it refers to atom types it does not define.
    print("\n--- using it ---")
    forcefield = app.ForceField(*CHARMM_BASE_FORCEFIELD, result.forcefield_xml)

    # build_ligand_xml cannot minimize a CHARMM ligand: a stream file records
    # internal coordinates, so the template's positions are all zero. The
    # geometry has to come from elsewhere - here the SDF, whose atom order the
    # stream file was written to match.
    molecule = Chem.MolFromMolFile(str(BEN_SDF), removeHs=False)
    conformer = molecule.GetConformer()
    topology = app.Topology()
    residue = topology.addResidue("BEN", topology.addChain("A"))
    atoms = {}
    for name, atom in zip(_stream_atom_names(result), molecule.GetAtoms(), strict=True):
        atoms[atom.GetIdx()] = topology.addAtom(name, app.Element.getByAtomicNumber(atom.GetAtomicNum()), residue)
    for bond in molecule.GetBonds():
        topology.addBond(atoms[bond.GetBeginAtomIdx()], atoms[bond.GetEndAtomIdx()])
    positions = unit.Quantity(
        [openmm.Vec3(*conformer.GetAtomPosition(i)) for i in range(molecule.GetNumAtoms())], unit.angstrom
    )

    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(positions)
    before = context.getState(getEnergy=True).getPotentialEnergy()
    openmm.LocalEnergyMinimizer.minimize(context, 10.0, 200)
    after = context.getState(getEnergy=True).getPotentialEnergy()
    print(f"BEN at the crystal geometry: {before.value_in_unit(unit.kilojoule_per_mole):>7.0f} kJ/mol")
    print(f"after 200 minimizer steps:   {after.value_in_unit(unit.kilojoule_per_mole):>7.0f} kJ/mol")

    show_refusals()

    print(
        "\nNote: build_forcefield_xml would find nothing to do for this ligand.\n"
        "charmm36.xml ships 814 residue templates, including CGenFF's own\n"
        "benzamidinium (RESI BAMI), so the base force field matches it already.\n"
        "The backend earns its keep on a drug-sized ligand, which no released\n"
        "parameter set contains - and which is exactly what ParamChem is for."
    )


def _stream_atom_names(result):
    """Atom names from the generated template, in order - the order the SDF is in."""
    root = ET.parse(result.forcefield_xml).getroot()
    return [atom.get("name") for atom in root.findall("./Residues/Residue/Atom")]


if __name__ == "__main__":
    main()
