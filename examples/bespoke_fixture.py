"""Synthetic BespokeFit-shaped input shared by the runnable example and tests.

These parameters demonstrate conversion and validation; they are not a fitted
force field for scientific use.
"""

from __future__ import annotations

import logging
from pathlib import Path

from openff.toolkit import ForceField, Molecule
from openff.units import unit

log = logging.getLogger(__name__)

#: The unconstrained build, which is what BespokeFit itself starts from.
RELEASE = "openff_unconstrained-2.2.1.offxml"


def fake_bespoke_offxml(path: Path, smiles: str, release: str = RELEASE) -> Path:
    """Write an installed release plus one matching, molecule-specific torsion."""
    forcefield = ForceField(release)
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    smirks = _bespoke_torsion_smirks(molecule)
    # Appended, so it wins: SMIRNOFF assignment is last match wins, which is how
    # BespokeFit's parameters override the general ones they sit on top of.
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
    log.info("Wrote %s: %s + 1 bespoke torsion %s", path.name, release, smirks)
    return path


def _bespoke_torsion_smirks(molecule: Molecule) -> str:
    """Select a matching torsion decorated unlike the stock Sage patterns."""
    for bond in molecule.bonds:
        outer_first = [n for n in bond.atom1.bonded_atoms if n is not bond.atom2]
        outer_last = [n for n in bond.atom2.bonded_atoms if n is not bond.atom1]
        if not (outer_first and outer_last):
            continue
        atoms = [outer_first[0], bond.atom1, bond.atom2, outer_last[0]]
        # Stock Sage torsions wildcard the outer atoms. Decorate every atom with
        # element and connectivity; "~" also matches double and aromatic bonds.
        smirks = "~".join(
            f"[#{atom.atomic_number}X{len(list(atom.bonded_atoms))}:{i + 1}]" for i, atom in enumerate(atoms)
        )
        if molecule.chemical_environment_matches(smirks):
            return smirks
    raise AssertionError(f"no usable torsion in {molecule.to_smiles()}")
