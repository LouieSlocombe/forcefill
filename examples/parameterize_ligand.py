"""Parameterize benzamidine bound to trypsin (PDB 3PTB), minimize it, then simulate.

The prepared structure and the ligand SDF come from prepare_trypsin_ben.py;
see examples/README.md for the preparation story. Run from this directory:

    python parameterize_ligand.py

Requires AmberTools on PATH (conda install -c conda-forge ambertools).
"""

from __future__ import annotations

import logging
from pathlib import Path

from openmm import LangevinMiddleIntegrator, app, unit

from forcefill import build_forcefield_xml

HERE = Path(__file__).parent
PREPARED_PDB = HERE / "data" / "trypsin_ben_prepared.pdb"
BEN_SDF = HERE / "data" / "benzamidinium.sdf"


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    # Step 1: build the force-field XML for everything amber14 cannot match.
    # - net_charges: benzamidinium is protonated (+1) at pH 7; without this,
    #   AM1-BCC gives plausible charges for the wrong ionization state.
    # - residue_files: antechamber reads the drawn SDF (explicit bond orders,
    #   aromatic ring, amidinium) instead of re-perceiving bonds from geometry.
    # - minimize: checks the numbers are physical, not just that a System
    #   builds - BEN alone in vacuum, then the whole complex.
    result = build_forcefield_xml(
        PREPARED_PDB,
        HERE / "ben_ff.xml",
        net_charges={"BEN": 1},
        residue_files={"BEN": BEN_SDF},
        workdir=HERE / "wd",
        minimize=True,
    )
    print(f"parameterized: {result.parameterized}")
    print(f"skipped:       {result.skipped or 'nothing'}")
    print(f"combined XML:  {result.forcefield_xml}")

    for label, report in [("BEN alone", result.minimizations["BEN"]), ("complex", result.full_minimization)]:
        print(
            f"{label:>10}: {report.initial_energy:>10.0f} -> {report.final_energy:>10.0f} kJ/mol "
            f"({report.n_atoms} atoms, max force {report.max_force:.0f} kJ/mol/nm)"
        )

    # Step 2: the produced XML is a normal force-field file - load it next to
    # the standard ones and simulate.
    pdb = app.PDBFile(str(PREPARED_PDB))
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.step(10)
    print(f"ran 10 steps of Langevin dynamics on {system.getNumParticles()} particles")


if __name__ == "__main__":
    main()
