"""Parameterize benzamidine bound to trypsin (PDB 3PTB), then build and minimize the complex.

The prepared structure and the ligand SDF come from prepare_trypsin_ben.py;
see examples/README.md for the preparation story. Run from this directory:

    python parameterize_ligand.py

Requires AmberTools on PATH (conda install -c conda-forge ambertools).
"""

import logging
from pathlib import Path

from openmm import LangevinMiddleIntegrator, app, unit

from forcefill import build_forcefield_xml

HERE = Path(__file__).parent
PREPARED_PDB = HERE / "data" / "trypsin_ben_prepared.pdb"
BEN_SDF = HERE / "data" / "benzamidinium.sdf"


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    # Step 1: build the force-field XML for everything amber14 cannot match.
    # - net_charges: benzamidinium is protonated (+1) at pH 7; without this,
    #   AM1-BCC would happily produce plausible-looking charges for the wrong
    #   ionization state.
    # - residue_files: antechamber reads the drawn SDF (explicit bond orders,
    #   aromatic ring, amidinium) instead of re-perceiving bonds from the
    #   PDB geometry.
    result = build_forcefield_xml(
        PREPARED_PDB,
        HERE / "ben_ff.xml",
        net_charges={"BEN": 1},
        residue_files={"BEN": BEN_SDF},
        workdir=HERE / "wd",
    )
    print(f"parameterized: {result.parameterized}")
    print(f"skipped:       {result.skipped or 'nothing'}")
    print(f"combined XML:  {result.forcefield_xml}")

    # Step 2: the produced XML is a normal force-field file - load it next to
    # the standard ones and simulate.
    pdb = app.PDBFile(str(PREPARED_PDB))
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=200)
    e1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    kj = unit.kilojoule_per_mole
    print(f"potential energy: {e0.value_in_unit(kj):.0f} -> {e1.value_in_unit(kj):.0f} kJ/mol after minimization")


if __name__ == "__main__":
    main()
