"""Find, classify and extract the residues a base force field cannot match.

The first half of the pipeline, and the only part that reads the input
structure. Everything here works on ``openmm.app.Topology`` objects:

    * :func:`find_nonstandard_residues` asks OpenMM which residues have no
      template;
    * :func:`_classify_unmatched` triages those into "parameterize" and "skip
      and say why", because a stand-alone GAFF treatment is valid for a
      free-standing hetero molecule and wrong for everything else;
    * :func:`extract_residue_to_pdb` and the two slicing helpers cut a single
      residue out of the structure, which is what antechamber and the
      per-residue checks are given.

Residue objects returned from here stay bound to the topology they came from -
they index into its coordinate array. Never mix them with a topology that was
rebuilt (by cleaning, say) in the meantime.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence

from openmm import app, unit

from ._residue_names import STANDARD_RESIDUES
from ._spec import DEFAULT_BASE_FORCEFIELD, LigandSpec, PathLike

log = logging.getLogger(__name__)

__all__ = [
    "extract_residue_to_pdb",
    "find_nonstandard_residues",
]


# --------------------------------------------------------------------------
# Identification
# --------------------------------------------------------------------------


def find_nonstandard_residues(
    topology: app.Topology,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
) -> list[app.topology.Residue]:
    """Return every residue that *base_forcefield* has no template for.

    This is OpenMM's own definition of "non-standard": a residue whose
    element/bond graph matches no registered template. Note that standard
    residues with missing atoms (e.g. a protein without hydrogens) also fail
    to match; :func:`_classify_unmatched` filters those out separately.
    """
    forcefield = app.ForceField(*base_forcefield)
    return forcefield.getUnmatchedResidues(topology)


def _classify_unmatched(
    unmatched: Sequence[app.topology.Residue],
) -> tuple[dict[str, app.topology.Residue], dict[str, str]]:
    """Split unmatched residues into those to parameterize and those to skip.

    Returns ``({name: representative_residue}, {name: skip_reason})``.
    """
    groups: dict[str, list[app.topology.Residue]] = defaultdict(list)
    for res in unmatched:
        groups[res.name.strip()].append(res)

    to_param: dict[str, app.topology.Residue] = {}
    skipped: dict[str, str] = {}
    for name, residues in groups.items():
        counts = sorted({sum(1 for _ in r.atoms()) for r in residues})
        rep = max(residues, key=lambda r: sum(1 for _ in r.atoms()))
        n_atoms = sum(1 for _ in rep.atoms())
        if len(counts) > 1:
            log.warning(
                "Copies of residue %s differ in atom count (%s); "
                "using the most complete copy (%d atoms) as the template.",
                name,
                counts,
                n_atoms,
            )
        if name in STANDARD_RESIDUES:
            skipped[name] = (
                f"standard residue that failed to match ({len(residues)} "
                "copies) - it is probably missing atoms or has non-standard "
                "atom names; repair the structure (e.g. with PDBFixer or "
                "Modeller.addHydrogens) instead of reparameterizing it"
            )
        elif n_atoms == 1:
            skipped[name] = (
                "monatomic species - use an ion parameter file for it (GAFF/antechamber cannot treat bare ions)"
            )
        elif n_linked := sum(1 for r in residues if any(True for _ in r.external_bonds())):
            skipped[name] = (
                f"covalently bonded to neighbouring residues ({n_linked} of "
                f"{len(residues)} copies) - a stand-alone GAFF "
                "parameterization is not valid for polymer-linked residues; "
                "cap the fragment and derive charges consistently with the "
                "backbone force field instead"
            )
        else:
            to_param[name] = rep
    return to_param, skipped


def _warn_unused_overrides(
    to_param: Mapping[str, app.topology.Residue],
    skipped: Mapping[str, str],
    net_charges: Mapping[str, int],
    multiplicities: Mapping[str, int],
    residue_files: Mapping[str, PathLike] | None = None,
    removed: Iterable[str] = (),
    ligands: Mapping[str, LigandSpec] | None = None,
) -> None:
    """Warn about ligands/net_charges/multiplicities/residue_files keys with no effect.

    A typo'd or case-mismatched key silently leaves the defaults (net
    charge 0, multiplicity 1, PDB extraction), which yields plausible but
    wrong AM1-BCC charges - the worst failure mode.

    *removed* names the residues cleaning deleted, so an override aimed at one
    of them reports the real cause rather than "matches no residue".
    """
    removed = set(removed)
    for label, mapping in (
        ("ligands", ligands or {}),
        ("net_charges", net_charges),
        ("multiplicities", multiplicities),
        ("residue_files", residue_files or {}),
    ):
        for key in mapping:
            if key in to_param:
                continue
            if key in removed:
                log.warning(
                    "%s[%r] has no effect: residue %s was removed from the structure by clean_structure=True.",
                    label,
                    key,
                    key,
                )
            elif key in skipped:
                log.warning(
                    "%s[%r] has no effect: residue %s is being skipped, not parameterized.",
                    label,
                    key,
                    key,
                )
            else:
                log.warning(
                    "%s[%r] does not match any residue selected for "
                    "parameterization %s; check the spelling and case of "
                    "the residue name.",
                    label,
                    key,
                    sorted(to_param),
                )


# --------------------------------------------------------------------------
# Extraction
# --------------------------------------------------------------------------


def _residue_subtopology(residue: app.topology.Residue) -> app.Topology:
    """Copy *residue* (atoms and internal bonds) into a fresh Topology."""
    sub_top = app.Topology()
    chain = sub_top.addChain("A")
    new_res = sub_top.addResidue(residue.name, chain)
    atom_map = {}
    for atom in residue.atoms():
        atom_map[atom] = sub_top.addAtom(atom.name, atom.element, new_res)
    for bond in residue.internal_bonds():
        sub_top.addBond(atom_map[bond.atom1], atom_map[bond.atom2])
    return sub_top


def _residue_positions(positions: unit.Quantity, residue: app.topology.Residue) -> unit.Quantity:
    """Slice *positions* (indexed by global atom index) down to *residue*, in its own atom order.

    The order must match :func:`_residue_subtopology`, which iterates the same
    ``residue.atoms()``.
    """
    return unit.Quantity(
        [positions[a.index].value_in_unit(unit.nanometer) for a in residue.atoms()],
        unit.nanometer,
    )


def extract_residue_to_pdb(
    positions: unit.Quantity,
    residue: app.topology.Residue,
    out_pdb: PathLike,
) -> str:
    """Write a single residue's atoms, internal bonds and coordinates to *out_pdb* and return the path."""
    for atom in residue.atoms():
        if atom.element is None:
            log.warning(
                "Atom %s in residue %s has no element assigned; antechamber "
                "may misread it. Check the element columns of the PDB.",
                atom.name,
                residue.name,
            )
    sub_top = _residue_subtopology(residue)

    with open(out_pdb, "w") as fh:
        app.PDBFile.writeFile(sub_top, _residue_positions(positions, residue), fh)
    return str(out_pdb)


def _describe_topology(topology: app.Topology) -> str:
    """Name *topology* for an error message: a lone residue by name, anything else by size."""
    n_atoms = topology.getNumAtoms()
    n_residues = topology.getNumResidues()
    if n_residues == 1:
        return f"residue {next(topology.residues()).name} ({n_atoms} atoms)"
    return f"the topology ({n_atoms} atoms, {n_residues} residues)"
