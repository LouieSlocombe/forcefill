"""Use a bespoke SMIRNOFF force field - BespokeFit's output - for one ligand.

`BespokeFit <https://github.com/openforcefield/openff-bespokefit>`_ fits torsion
parameters to quantum-chemical torsion drives for one specific molecule and
writes an OFFXML: a stock OpenFF release with the bespoke torsions layered on
top. That file is not something OpenMM can load. This is the last mile - OFFXML
in, ffxml out, validated on the way through:

    build_ligand_xml(
        {"BEN": LigandSpec(file="ben.sdf", forcefield="ben_bespoke.offxml")},
        "ben.xml",
        backend="smirnoff",
    )

forcefill neither runs nor requires BespokeFit. To produce the input for real:

    openff-bespoke executor run --file ben.sdf \
                               --workflow default \
                               --output-force-field ben_bespoke.offxml

That wants psi4 or xtb, torsiondrive and ForceBalance, and takes hours. So this
script *synthesizes* a file of the same shape from the installed release - one
extra, molecule-specific proper torsion appended to Sage - which is enough to
show the plumbing and all three refusals.

Run from this directory:

    python parameterize_ligand_bespoke.py

Needs no AmberTools: the smirnoff backend is pure OpenFF.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from openff.toolkit import ForceField, Molecule
from openff.units import unit
from openmm import app

from forcefill import LigandSpec, build_ligand_xml

HERE = Path(__file__).parent
BEN_SDF = HERE / "data" / "benzamidinium.sdf"
WORKDIR = HERE / "bespoke_work"

#: The unconstrained build, which is what BespokeFit itself starts from. The
#: constrained one is used below to show why.
RELEASE = "openff_unconstrained-2.2.1.offxml"


def fake_bespoke_offxml(path: Path, smiles: str, release: str = RELEASE) -> Path:
    """Write an OFFXML shaped like BespokeFit's: a release plus one bespoke torsion."""
    forcefield = ForceField(release)
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    bond = next(
        b
        for b in molecule.bonds
        if any(n is not b.atom2 for n in b.atom1.bonded_atoms) and any(n is not b.atom1 for n in b.atom2.bonded_atoms)
    )
    atoms = [
        next(n for n in bond.atom1.bonded_atoms if n is not bond.atom2),
        bond.atom1,
        bond.atom2,
        next(n for n in bond.atom2.bonded_atoms if n is not bond.atom1),
    ]
    # Every atom decorated with element and connectivity: stock Sage torsions
    # wildcard the outer two, so nothing in the release is written this way.
    smirks = "~".join(f"[#{a.atomic_number}X{len(list(a.bonded_atoms))}:{i + 1}]" for i, a in enumerate(atoms))
    forcefield.get_parameter_handler("ProperTorsions").add_parameter(
        {
            "smirks": smirks,
            "periodicity1": 1,
            "phase1": 0.0 * unit.degree,
            "k1": 1.234 * unit.kilocalorie_per_mole,
            "idivf1": 1.0,
        }
    )
    forcefield.to_file(str(path))
    print(f"  wrote {path.name}: {release} + 1 bespoke torsion {smirks}")
    return path


def expect_refusal(label: str, **kwargs: object) -> None:
    """Run a build that should be refused, and print the reason."""
    print(f"\n--- {label} ---")
    try:
        build_ligand_xml(output_xml=WORKDIR / "refused.xml", backend="smirnoff", workdir=WORKDIR / "wd", **kwargs)
    except ValueError as exc:
        print(f"  refused: {str(exc).splitlines()[0]}")
    else:
        print("  NOT refused - this example is out of date with the library")


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    WORKDIR.mkdir(exist_ok=True)

    print("\n--- building a stand-in for BespokeFit's output ---")
    bespoke = fake_bespoke_offxml(WORKDIR / "ben_bespoke.offxml", "NC(=[NH2+])c1ccccc1")

    print("\n--- parameterizing BEN against it ---")
    result = build_ligand_xml(
        {"BEN": LigandSpec(file=BEN_SDF, forcefield=bespoke)},
        WORKDIR / "ben_bespoke.xml",
        backend="smirnoff",
        workdir=WORKDIR / "wd",
        minimize=True,
    )
    report = result.minimizations["BEN"]
    print(f"  {result.forcefield_xml}")
    print(f"  {report.initial_energy:.1f} -> {report.final_energy:.1f} kJ/mol over {report.n_atoms} atoms")

    # The whole point: it is an ordinary ffxml, loadable next to the standard set.
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    print("  loads alongside amber14-all.xml")

    # The refusals are the feature. Each of these otherwise costs you the QC.
    expect_refusal(
        "a constrained OFFXML (the obvious-looking 'openff-2.2.1.offxml')",
        ligands={
            "BEN": LigandSpec(
                file=BEN_SDF,
                forcefield=fake_bespoke_offxml(
                    WORKDIR / "constrained.offxml", "NC(=[NH2+])c1ccccc1", "openff-2.2.1.offxml"
                ),
            )
        },
    )
    expect_refusal(
        "an OFFXML fitted for a different molecule",
        ligands={
            "BEN": LigandSpec(
                file=BEN_SDF,
                forcefield=fake_bespoke_offxml(WORKDIR / "hexanol.offxml", "CCCCCCO"),
            )
        },
    )
    expect_refusal(
        "a path that does not exist",
        ligands={"BEN": LigandSpec(file=BEN_SDF, forcefield=WORKDIR / "typo.offxml")},
    )

    shutil.rmtree(WORKDIR, ignore_errors=True)


if __name__ == "__main__":
    main()
