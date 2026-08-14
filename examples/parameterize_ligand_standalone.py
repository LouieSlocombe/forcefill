"""Parameterize benzamidinium on its own - no protein, no PDB - through both backends.

The companion to parameterize_ligand.py, which needs the prepared complex. Here
the ligand file is the entire input, which is the common case when you are
building a ligand library or preparing a molecule before you have a structure to
put it in.

Two things worth watching in the output:

  * the net charge is never stated. data/benzamidinium.sdf carries an M CHG
    record saying +1, and forcefill reads it - passing net_charges={"BEN": 1} by
    hand, as parameterize_ligand.py does, is now belt and braces rather than the
    difference between right and wrong answers.
  * the same molecule goes through GAFF2 and through OpenFF Sage, and the two
    disagree about its energy. That is expected: they are different force
    fields. What matters is that each is internally consistent, which is what
    minimize=True checks.

Run from this directory:

    python parameterize_ligand_standalone.py

Requires AmberTools on PATH for the GAFF half; the OpenFF stack the SMIRNOFF
half needs is an ordinary forcefill dependency, so it is already there.
"""

import logging
import shutil
from pathlib import Path

from openmm import app

from forcefill import LigandSpec, build_ligand_xml

HERE = Path(__file__).parent
BEN_SDF = HERE / "data" / "benzamidinium.sdf"


def run(backend, output_xml, workdir):
    """Parameterize the ligand with one backend and report what came out."""
    print(f"\n--- {backend} ---")
    result = build_ligand_xml(
        BEN_SDF,  # the file name supplies the residue name: benzamidinium -> BEN
        output_xml,
        backend=backend,
        workdir=workdir,
        minimize=True,
    )
    report = result.minimizations["BEN"]
    print(f"parameterized: {result.parameterized}")
    print(f"combined XML:  {result.forcefield_xml}")
    print(
        f"BEN: {report.initial_energy:>9.0f} -> {report.final_energy:>9.0f} kJ/mol "
        f"({report.n_atoms} atoms, max force {report.max_force:.0f} kJ/mol/nm)"
    )

    # The product is an ordinary force-field file: load it next to the standard set.
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    print(f"loads alongside amber14; BEN template present: {'BEN' in forcefield._templates}")
    return result


def main():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    run("smirnoff", HERE / "ben_standalone_smirnoff.xml", HERE / "wd_standalone_smirnoff")

    # AmberTools is the one piece that cannot be a Python dependency, so it is
    # the only thing worth checking for.
    if not shutil.which("antechamber"):
        print("\nskipping the gaff backend: antechamber is not on PATH")
        return
    run("gaff", HERE / "ben_standalone_gaff.xml", HERE / "wd_standalone_gaff")

    # Mixing them in one call: each ligand names its own backend, and the two
    # force fields are merged into a single XML that OpenMM will load.
    print("\n--- mixed, in one XML ---")
    result = build_ligand_xml(
        {
            "BEN": LigandSpec(file=BEN_SDF, backend="smirnoff"),
            "MOL": LigandSpec(smiles="CO", backend="gaff"),
        },
        HERE / "ben_standalone_mixed.xml",
        workdir=HERE / "wd_standalone_mixed",
        minimize=True,
    )
    for name, report in sorted(result.minimizations.items()):
        print(f"{name}: {report.initial_energy:>9.0f} -> {report.final_energy:>9.0f} kJ/mol")
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    print(f"one XML, both force fields: {result.forcefield_xml}")


if __name__ == "__main__":
    main()
