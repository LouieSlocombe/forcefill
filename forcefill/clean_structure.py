"""Strip crystallographic water, bulk counter-ions and crystallization additives from a structure.

What it removes, by default:

    * **water** - ``HOH``/``WAT``/``SOL`` and the other model aliases.
    * **bulk counter-ions** - group-1 cations and group-17 anions (``NA``,
      ``CL``, ``K``, ...). These come from the buffer or from neutralizing the
      box; they occupy no defined site and you re-add them with
      ``Modeller.addSolvent(ionicStrength=...)`` anyway.
    * **crystallization additives** - cryoprotectants (``GOL``, ``EDO``,
      ``PEG``), solvents (``DMS``), buffers (``EPE``, ``TRS``), precipitants
      (``SO4``, ``PO4``, ``ACT``) and reductants (``BME``, ``DTT``). Left in
      place these are exactly what :func:`~forcefill.build_forcefield_xml`
      sends to antechamber - a free-standing hetero molecule looks like a
      ligand - so a raw crystal structure burns AM1-BCC cycles on glycerol,
      and since X-ray additives carry no hydrogens the resulting charges are
      meaningless anyway.

What it does **not** remove:

    * **Structural metals** (``CA``, ``ZN``, ``MG``, ``MN``, ``FE``, ...) are
      kept by default and reported in :attr:`CleaningResult.retained`. Trypsin's
      Ca2+ (PDB 3PTB, residue ``CA`` 480) rigidifies the calcium-binding loop;
      deleting it silently changes the science. The asymmetry decides the
      default: dropping a needed metal is silent and wrong, keeping an unwanted
      one is visible in the log and reversible with
      ``remove_structural_metals=True``.
    * **Anything covalently bonded to a neighbouring residue.** ``Modeller``
      drops the bonds along with the atoms and never says so, leaving the
      surviving partner with an unsatisfied valence. Those residues are kept
      and reported instead.

This module is **subtractive only**. It does not add missing atoms, model
missing loops, protonate anything, select chains or strip hydrogens. Repairing
a structure is [PDBFixer's](https://github.com/openmm/pdbfixer) job and
forcefill deliberately does not duplicate it - clean *after* you repair.

The residue-name tables live in :mod:`forcefill._residue_names`, which also
documents the het codes deliberately left out (cofactors, glycans, ``BEN``) and
the borderline ones (``IMD``, ``AZI``, ``SO4``). ``keep=`` and ``extra_remove=``
are the escape hatches for both directions.

Example:
    >>> from forcefill import clean_pdb
    >>> result = clean_pdb("3ptb.pdb", "3ptb_clean.pdb")
    >>> result.removed["HOH"]
    ('water', 62)
    >>> result.retained["CA"]
    'structural metal retained by default ...'
"""

from __future__ import annotations

import logging
import os
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from openmm import app, unit

from ._residue_names import (
    ADDITIVE_RESIDUES,
    BULK_ION_RESIDUES,
    ION_ELEMENTS,
    STANDARD_RESIDUES,
    STRUCTURAL_METAL_RESIDUES,
    WATER_RESIDUES,
)

log = logging.getLogger(__name__)

__all__ = [
    "ADDITIVE_RESIDUES",
    "BULK_ION_RESIDUES",
    "STRUCTURAL_METAL_RESIDUES",
    "WATER_RESIDUES",
    "CleaningResult",
    "clean_pdb",
    "clean_topology",
]

PathLike = str | os.PathLike

#: Largest residue still accepted as a water molecule. Three real atoms, plus
#: the virtual site of a 4- or 5-site model (TIP4P's ``M``, TIP5P's ``EP``),
#: whose element is None.
_MAX_WATER_ATOMS = 5


@dataclass
class CleaningResult:
    """What one cleaning pass removed, and what it deliberately did not."""

    #: Atom count of the input topology.
    n_atoms_before: int
    #: Atom count after cleaning; equal to *n_atoms_before* when nothing matched.
    n_atoms_after: int
    #: Deleted residue names, mapped to ``(category, n_copies)``. The category
    #: is one of ``"water"``, ``"bulk_ion"``, ``"structural_metal"``,
    #: ``"additive"`` or ``"requested"`` (from ``extra_remove``).
    removed: dict[str, tuple[str, int]] = field(default_factory=dict)
    #: Residue names that matched a removal rule but were kept, mapped to the
    #: reason - the audit trail for "why is the calcium still in my file".
    #: A name can appear here *and* in :attr:`removed` when only some copies
    #: were kept.
    retained: dict[str, str] = field(default_factory=dict)
    #: Where the cleaned structure was written; None for :func:`clean_topology`.
    output_pdb: str | None = None

    @property
    def n_atoms_removed(self) -> int:
        """Atoms deleted by this pass."""
        return self.n_atoms_before - self.n_atoms_after

    @property
    def n_residues_removed(self) -> int:
        """Residue copies deleted by this pass, summed over every name."""
        return sum(n for _, n in self.removed.values())


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------


def _linked_residues(topology: app.Topology) -> set[int]:
    """Indices of residues with at least one bond to a *different* residue.

    One pass over the bonds. ``Residue.external_bonds()`` re-scans every bond in
    the topology per residue and tests membership against a list, which is
    quadratic in the residue count - 250x slower than this on a 287-residue
    structure, and far worse on a solvated box.
    """
    linked: set[int] = set()
    for atom1, atom2 in topology.bonds():
        if atom1.residue is not atom2.residue:
            linked.add(atom1.residue.index)
            linked.add(atom2.residue.index)
    return linked


def _ion_mismatch(residue: app.topology.Residue, name: str, n_atoms: int) -> str | None:
    """Why *residue* cannot be treated as the ion its name claims, or None if it can.

    A name alone is never enough to delete something. ``I`` is iodide but also
    inosine; a residue named ``CA`` holding a carbon is an alpha carbon that
    lost its residue name somewhere upstream. Both are caught by requiring a
    single atom carrying the expected element.
    """
    if n_atoms != 1:
        return (
            f"named like the ion {name} but has {n_atoms} atoms; a monatomic "
            "ion has exactly one, so this is a different molecule with a "
            "colliding name"
        )
    atom = next(residue.atoms())
    if atom.element is None:
        return (
            f"named like the ion {name} but its atom has no element assigned, "
            "so the name cannot be confirmed; fill in the element columns "
            "(77-78) of the PDB"
        )
    expected = ION_ELEMENTS[name]
    if atom.element.symbol != expected:
        return (
            f"named like the ion {name} but its atom is {atom.element.symbol}, "
            f"not {expected}; the residue name and the element disagree"
        )
    return None


def _classify_residue(
    residue: app.topology.Residue,
    name: str,
    *,
    linked: set[int],
    remove_water: bool,
    remove_ions: bool,
    remove_additives: bool,
    remove_structural_metals: bool,
    keep: frozenset[str],
    extra_remove: frozenset[str],
) -> tuple[str | None, str | None]:
    """Decide the fate of one residue.

    Returns ``(category, reason)``: a category names the rule that removes it,
    a reason explains why a residue that matched a removal rule is being kept
    anyway. Both are None when no rule applied at all.

    Precedence is ``keep`` > covalent bonds > ``extra_remove`` > the categories,
    so an explicit request outranks the built-in tables but never the safety
    invariant.
    """
    if name in keep:
        return None, None

    is_requested = name in extra_remove
    is_water = remove_water and name in WATER_RESIDUES
    is_bulk_ion = remove_ions and name in BULK_ION_RESIDUES
    # Metals are considered whenever either ion flag is on: with the defaults
    # that is what produces the "kept your calcium" audit trail, and with both
    # off the caller has said hands-off, so there is nothing to report.
    is_metal = (remove_ions or remove_structural_metals) and name in STRUCTURAL_METAL_RESIDUES
    is_additive = remove_additives and name in ADDITIVE_RESIDUES
    if not (is_requested or is_water or is_bulk_ion or is_metal or is_additive):
        return None, None

    # Checked before every removal rule: Modeller drops the bonds along with
    # the atoms, silently leaving the partner with an unsatisfied valence.
    if residue.index in linked:
        n_linked = sum(1 for _ in residue.external_bonds())
        return None, (
            f"covalently bonded to {n_linked} neighbouring atom(s) - deleting "
            "it would leave the partner with an unsatisfied valence and no "
            "warning; cap or repair the structure instead"
        )

    if is_requested:
        return "requested", None

    n_atoms = len(residue)

    # Water first: the virtual site of a 4- or 5-site model has no element, so
    # it must not reach the element checks below.
    if is_water:
        if n_atoms > _MAX_WATER_ATOMS:
            return None, f"named like water but has {n_atoms} atoms; not removing it"
        return "water", None

    if is_bulk_ion or is_metal:
        mismatch = _ion_mismatch(residue, name, n_atoms)
        if mismatch is not None:
            return None, mismatch
        if is_bulk_ion:
            return "bulk_ion", None
        if remove_structural_metals:
            return "structural_metal", None
        return None, (
            "structural metal retained by default (it may be a cofactor or a "
            "fold-stabilising site, as trypsin's Ca2+ is); pass "
            "remove_structural_metals=True to strip it"
        )

    return "additive", None


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------


def clean_topology(
    topology: app.Topology,
    positions: unit.Quantity,
    *,
    remove_water: bool = True,
    remove_ions: bool = True,
    remove_additives: bool = True,
    remove_structural_metals: bool = False,
    keep: Iterable[str] = (),
    extra_remove: Iterable[str] = (),
) -> tuple[app.Topology, unit.Quantity, CleaningResult]:
    """Remove solvent, bulk ions and crystallization additives from a topology.

    Args:
        topology: Structure to clean. Not modified; a new Topology is returned.
        positions: Coordinates, indexed by atom index and the same length as
            *topology*.
        remove_water: Remove water (:data:`~forcefill.WATER_RESIDUES`).
        remove_ions: Remove bulk counter-ions
            (:data:`~forcefill.BULK_ION_RESIDUES`).
        remove_additives: Remove crystallization additives
            (:data:`~forcefill.ADDITIVE_RESIDUES`).
        remove_structural_metals: Also remove
            :data:`~forcefill.STRUCTURAL_METAL_RESIDUES`. Off by default: those
            metals are frequently catalytic or fold-stabilising, and dropping
            one silently is a much worse failure than keeping it.
        keep: Residue names never to remove, whatever category they fall in.
            The escape hatch for an additive that is really your ligand.
        extra_remove: Additional residue names to remove. Rejected with
            ValueError for standard residues, so a typo cannot shred a protein.

    Returns:
        ``(topology, positions, CleaningResult)``. Rebind both halves together -
        the returned positions are indexed by the *new* topology, and the old
        Residue objects belong to a discarded topology.
    """
    n_atoms_before = topology.getNumAtoms()
    if len(positions) != n_atoms_before:
        raise ValueError(
            f"positions has {len(positions)} entries but the topology has "
            f"{n_atoms_before} atoms. They must correspond one to one; pass "
            "the positions that came with this topology."
        )

    keep = frozenset(n.strip().upper() for n in keep)
    extra_remove = frozenset(n.strip().upper() for n in extra_remove)
    if bad := sorted(extra_remove & STANDARD_RESIDUES):
        raise ValueError(
            f"extra_remove names standard residue(s) {bad}, which would delete "
            "part of the protein or nucleic acid. Remove them from the list; "
            "if a standard residue really is junk, delete it yourself with "
            "openmm.app.Modeller."
        )

    if not (remove_water or remove_ions or remove_additives or remove_structural_metals or extra_remove):
        log.info("Nothing enabled for removal; returning the structure unchanged.")
        return topology, positions, CleaningResult(n_atoms_before, n_atoms_before)

    linked = _linked_residues(topology)
    to_delete: list[app.topology.Residue] = []
    counts: dict[str, int] = defaultdict(int)
    categories: dict[str, str] = {}
    retained: dict[str, str] = {}
    kept_counts: dict[str, int] = defaultdict(int)
    seen: set[str] = set()

    for residue in topology.residues():
        name = (residue.name or "").strip().upper()
        seen.add(name)
        category, reason = _classify_residue(
            residue,
            name,
            linked=linked,
            remove_water=remove_water,
            remove_ions=remove_ions,
            remove_additives=remove_additives,
            remove_structural_metals=remove_structural_metals,
            keep=keep,
            extra_remove=extra_remove,
        )
        if category is not None:
            to_delete.append(residue)
            counts[name] += 1
            categories[name] = category
        elif reason is not None:
            retained[name] = reason
            kept_counts[name] += 1

    for name, reason in sorted(retained.items()):
        log.warning("Keeping %s (%d copies): %s", name, kept_counts[name], reason)
    for name in sorted(keep | extra_remove):
        if name not in seen:
            log.warning(
                "Residue name %r matches nothing in the structure; check its spelling and case.",
                name,
            )

    if not to_delete:
        log.info("Nothing to remove: the structure is already clean.")
        return topology, positions, CleaningResult(n_atoms_before, n_atoms_before, retained=retained)

    modeller = app.Modeller(topology, positions)
    modeller.delete(to_delete)
    n_atoms_after = modeller.topology.getNumAtoms()

    removed = {name: (categories[name], counts[name]) for name in counts}
    for name in sorted(removed):
        # Same "N of M copies" shape as the parameterization skip messages: a
        # name lands in both maps when only some of its copies were removed.
        if name in kept_counts:
            log.info(
                "Removed %d of %d %s (%s); %d kept, see the warning above.",
                counts[name],
                counts[name] + kept_counts[name],
                name,
                categories[name],
                kept_counts[name],
            )
        else:
            log.info("Removed %d %s (%s).", counts[name], name, categories[name])
    log.info(
        "Cleaned structure: %d -> %d atoms (%d residues removed).",
        n_atoms_before,
        n_atoms_after,
        sum(counts.values()),
    )
    if n_atoms_after == 0:
        log.warning("Cleaning removed every atom; the structure was nothing but solvent.")

    return (
        modeller.topology,
        modeller.positions,
        CleaningResult(n_atoms_before, n_atoms_after, removed=removed, retained=retained),
    )


def clean_pdb(
    pdb_file: PathLike,
    output_pdb: PathLike,
    *,
    remove_water: bool = True,
    remove_ions: bool = True,
    remove_additives: bool = True,
    remove_structural_metals: bool = False,
    keep: Iterable[str] = (),
    extra_remove: Iterable[str] = (),
) -> CleaningResult:
    """Read *pdb_file*, clean it with :func:`clean_topology`, write *output_pdb*.

    The keyword arguments are :func:`clean_topology`'s; see it for what each
    category covers and why structural metals are kept.

    Only the first MODEL is read and only the primary altLoc is kept - that is
    OpenMM's PDB reader, and the same is true of
    :func:`~forcefill.build_forcefield_xml`.

    Returns:
        CleaningResult, with ``output_pdb`` set to the path written.
    """
    pdb_file = Path(pdb_file)
    output_pdb = Path(output_pdb)
    # A subtractive tool that overwrites its own source is unrecoverable.
    if output_pdb.resolve() == pdb_file.resolve():
        raise ValueError(
            f"output_pdb and pdb_file are the same file ({pdb_file}). Cleaning "
            "deletes atoms, so overwriting the input would lose them for good; "
            "write to a different path."
        )

    pdb = app.PDBFile(str(pdb_file))
    topology, positions, result = clean_topology(
        pdb.topology,
        pdb.positions,
        remove_water=remove_water,
        remove_ions=remove_ions,
        remove_additives=remove_additives,
        remove_structural_metals=remove_structural_metals,
        keep=keep,
        extra_remove=extra_remove,
    )

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as fh:
        app.PDBFile.writeFile(topology, positions, fh, keepIds=True)
    log.info("Wrote cleaned structure: %s", output_pdb)
    result.output_pdb = str(output_pdb)
    return result
