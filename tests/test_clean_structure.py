"""Tests for :mod:`forcefill.clean_structure`.

Hermetic: no AmberTools executables and no PDB fixtures are required. Every
structure is built in memory by :mod:`tests.helpers`.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Sequence
from pathlib import Path

import pytest

pytest.importorskip("openmm")

from openmm import Vec3, app, unit
from openmm.app import element

from forcefill import clean_structure
from forcefill._residue_names import (
    ADDITIVE_RESIDUES,
    BULK_ION_RESIDUES,
    ION_ELEMENTS,
    STANDARD_RESIDUES,
    STRUCTURAL_METAL_RESIDUES,
    WATER_RESIDUES,
)
from forcefill.clean_structure import clean_pdb, clean_topology
from tests.helpers import (
    add_glycerol_residue,
    add_ion_residue,
    add_methanol_residue,
    add_water_residue,
    bond_across_residues,
    write_dirty_pdb,
)


def _build(
    *builders: Callable[[app.Topology], Sequence[tuple[float, float, float]]],
) -> tuple[app.Topology, unit.Quantity]:
    """Assemble a topology from ``add_*`` helpers; returns ``(topology, positions)``."""
    top = app.Topology()
    xyz = []
    for builder in builders:
        xyz += builder(top)
    return top, unit.Quantity([Vec3(*p) for p in xyz], unit.angstrom)


def _names(topology: app.Topology) -> list[str]:
    return sorted({res.name for res in topology.residues()})


# -- residue-name tables ---------------------------------------------------


def test_category_sets_are_disjoint() -> None:
    sets = {
        "water": WATER_RESIDUES,
        "bulk_ion": BULK_ION_RESIDUES,
        "structural_metal": STRUCTURAL_METAL_RESIDUES,
        "additive": ADDITIVE_RESIDUES,
    }
    for a_name, a in sets.items():
        for b_name, b in sets.items():
            if a_name < b_name:
                assert not (a & b), f"{a_name} and {b_name} share {sorted(a & b)}"


def test_removal_sets_never_touch_protein_or_nucleic_acids() -> None:
    # Water is in STANDARD_RESIDUES on purpose; nothing else may be.
    protein_and_nucleic = STANDARD_RESIDUES - WATER_RESIDUES
    for names in (BULK_ION_RESIDUES, STRUCTURAL_METAL_RESIDUES, ADDITIVE_RESIDUES):
        assert not (names & protein_and_nucleic)


def test_standard_residues_contains_water() -> None:
    # Locks the refactor: _classify_unmatched still recognises water as standard.
    assert WATER_RESIDUES <= STANDARD_RESIDUES


def test_known_ligands_are_not_removable() -> None:
    # BEN is this package's own example ligand; the rest are common cofactors.
    for name in ("BEN", "HEM", "NAD", "FAD", "ATP", "ADP", "NAG"):
        assert name not in BULK_ION_RESIDUES | STRUCTURAL_METAL_RESIDUES | ADDITIVE_RESIDUES


def test_every_ion_code_maps_to_a_real_element() -> None:
    ions = BULK_ION_RESIDUES | STRUCTURAL_METAL_RESIDUES
    assert set(ION_ELEMENTS) == ions
    for name, symbol in ION_ELEMENTS.items():
        assert element.get_by_symbol(symbol) is not None, name


# -- clean_topology --------------------------------------------------------


def test_positions_stay_aligned_with_the_topology() -> None:
    # The invariant everything else depends on: a surviving atom must keep its
    # own coordinates, not inherit a deleted neighbour's.
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0)),
        lambda t: add_glycerol_residue(t, origin=(0.0, -8.0, 0.0)),
    )
    before = {(atom.residue.chain.id, atom.residue.name, atom.name): positions[atom.index] for atom in top.atoms()}

    new_top, new_positions, _ = clean_topology(top, positions)

    assert len(new_positions) == new_top.getNumAtoms()
    assert new_positions.unit.is_compatible(unit.nanometer)
    for atom in new_top.atoms():
        key = (atom.residue.chain.id, atom.residue.name, atom.name)
        expected = before[key].value_in_unit(unit.nanometer)
        actual = new_positions[atom.index].value_in_unit(unit.nanometer)
        assert actual == pytest.approx([expected.x, expected.y, expected.z])


def test_removes_water_and_keeps_the_ligand() -> None:
    top, positions = _build(add_methanol_residue, lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0)))
    new_top, _, result = clean_topology(top, positions)
    assert _names(new_top) == ["LIG"]
    assert result.removed == {"HOH": ("water", 1)}
    assert result.n_atoms_removed == 3
    assert result.n_residues_removed == 1


def test_removes_bulk_ion_but_keeps_structural_metal(caplog: pytest.LogCaptureFixture) -> None:
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_ion_residue(t, "NA", element.sodium, chain_id="N", origin=(-5.0, 0.0, 0.0)),
        lambda t: add_ion_residue(t, "CA", element.calcium, chain_id="C", origin=(-8.0, 0.0, 0.0)),
    )
    with caplog.at_level(logging.WARNING):
        new_top, _, result = clean_topology(top, positions)

    assert _names(new_top) == ["CA", "LIG"]
    assert result.removed == {"NA": ("bulk_ion", 1)}
    assert "structural metal retained by default" in result.retained["CA"]
    assert "remove_structural_metals=True" in caplog.text


def test_strips_structural_metal_when_asked() -> None:
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_ion_residue(t, "CA", element.calcium, chain_id="C", origin=(-8.0, 0.0, 0.0)),
    )
    new_top, _, result = clean_topology(top, positions, remove_structural_metals=True)
    assert _names(new_top) == ["LIG"]
    assert result.removed == {"CA": ("structural_metal", 1)}
    assert result.retained == {}


def test_metals_are_silent_when_ion_handling_is_off() -> None:
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_ion_residue(t, "CA", element.calcium, chain_id="C", origin=(-8.0, 0.0, 0.0)),
    )
    _, _, result = clean_topology(top, positions, remove_ions=False)
    assert result.retained == {}


def test_removes_crystallization_additives() -> None:
    top, positions = _build(add_methanol_residue, lambda t: add_glycerol_residue(t, origin=(0.0, -8.0, 0.0)))
    new_top, _, result = clean_topology(top, positions)
    assert _names(new_top) == ["LIG"]
    assert result.removed == {"GOL": ("additive", 1)}


def test_additives_survive_when_disabled() -> None:
    top, positions = _build(add_methanol_residue, lambda t: add_glycerol_residue(t, origin=(0.0, -8.0, 0.0)))
    new_top, _, result = clean_topology(top, positions, remove_additives=False)
    assert _names(new_top) == ["GOL", "LIG"]
    assert result.removed == {}


def test_never_deletes_a_covalently_bonded_residue(caplog: pytest.LogCaptureFixture) -> None:
    # Modeller drops the bonds along with the atoms and says nothing, leaving
    # the partner short a valence.
    top, positions = _build(add_methanol_residue, lambda t: add_glycerol_residue(t, origin=(0.0, -8.0, 0.0)))
    bond_across_residues(top, "GOL", "LIG")

    with caplog.at_level(logging.WARNING):
        new_top, _, result = clean_topology(top, positions)

    assert _names(new_top) == ["GOL", "LIG"]
    assert result.removed == {}
    assert "covalently bonded" in result.retained["GOL"]
    assert "covalently bonded" in caplog.text


def test_bonded_copy_is_kept_while_free_copies_go() -> None:
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_glycerol_residue(t, chain_id="G1", origin=(0.0, -8.0, 0.0)),
        lambda t: add_glycerol_residue(t, chain_id="G2", origin=(0.0, -16.0, 0.0)),
    )
    # bond_across_residues links the first GOL only.
    bond_across_residues(top, "GOL", "LIG")

    new_top, _, result = clean_topology(top, positions)

    assert result.removed == {"GOL": ("additive", 1)}
    assert "covalently bonded" in result.retained["GOL"]
    assert sorted(res.name for res in new_top.residues()) == ["GOL", "LIG"]


def test_ion_name_needs_a_single_atom() -> None:
    # "I" is iodide, but it is also inosine. Atom count settles it.
    def add_inosine(top: app.Topology) -> list[tuple[float, float, float]]:
        chain = top.addChain("R")
        res = top.addResidue("I", chain)
        for i in range(12):
            top.addAtom(f"C{i}", element.carbon, res)
        return [(float(i), 0.0, 0.0) for i in range(12)]

    top, positions = _build(add_methanol_residue, add_inosine)
    new_top, _, result = clean_topology(top, positions)

    assert "I" in _names(new_top)
    assert result.removed == {}
    assert "12 atoms" in result.retained["I"]


def test_ion_name_needs_a_matching_element() -> None:
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_ion_residue(t, "NA", element.carbon, chain_id="N", origin=(-5.0, 0.0, 0.0)),
    )
    new_top, _, result = clean_topology(top, positions)

    assert "NA" in _names(new_top)
    assert "its atom is C, not Na" in result.retained["NA"]


def test_ion_with_unknown_element_is_kept() -> None:
    # Deleting on the strength of a name alone is the one real silent-error risk.
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_ion_residue(t, "NA", None, chain_id="N", origin=(-5.0, 0.0, 0.0)),
    )
    new_top, _, result = clean_topology(top, positions)

    assert "NA" in _names(new_top)
    assert "no element assigned" in result.retained["NA"]


def test_four_site_water_is_still_removed() -> None:
    # The TIP4P M particle has element None, so water must be classified before
    # any element check runs.
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0), virtual_site=True),
    )
    new_top, _, result = clean_topology(top, positions)
    assert _names(new_top) == ["LIG"]
    assert result.removed == {"HOH": ("water", 1)}


def test_oversized_residue_named_like_water_is_kept() -> None:
    def add_fat_hoh(top: app.Topology) -> list[tuple[float, float, float]]:
        chain = top.addChain("X")
        res = top.addResidue("HOH", chain)
        for i in range(8):
            top.addAtom(f"O{i}", element.oxygen, res)
        return [(float(i), 5.0, 0.0) for i in range(8)]

    top, positions = _build(add_methanol_residue, add_fat_hoh)
    new_top, _, result = clean_topology(top, positions)

    assert "HOH" in _names(new_top)
    assert "8 atoms" in result.retained["HOH"]


def test_residue_names_are_normalized() -> None:
    # Topology.addResidue does not strip, unlike PDBFile: a converter-built
    # topology can carry " hoh" or "gol".
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_water_residue(t, name=" hoh ", origin=(5.0, 0.0, 0.0)),
        lambda t: add_glycerol_residue(t, name="gol", origin=(0.0, -8.0, 0.0)),
    )
    new_top, _, result = clean_topology(top, positions)

    assert _names(new_top) == ["LIG"]
    assert sorted(result.removed) == ["GOL", "HOH"]


def test_keep_overrides_the_category() -> None:
    top, positions = _build(
        add_methanol_residue,
        lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0)),
        lambda t: add_glycerol_residue(t, origin=(0.0, -8.0, 0.0)),
    )
    new_top, _, result = clean_topology(top, positions, keep=["gol"])

    assert _names(new_top) == ["GOL", "LIG"]
    assert result.removed == {"HOH": ("water", 1)}
    # keep is an explicit decision, not something to warn about.
    assert "GOL" not in result.retained


def test_extra_remove_deletes_a_name_the_tables_do_not_know() -> None:
    top, positions = _build(add_methanol_residue)
    new_top, _, result = clean_topology(top, positions, extra_remove=["lig"])
    assert new_top.getNumAtoms() == 0
    assert result.removed == {"LIG": ("requested", 1)}


def test_extra_remove_refuses_standard_residues() -> None:
    top, positions = _build(add_methanol_residue)
    with pytest.raises(ValueError, match="standard residue"):
        clean_topology(top, positions, extra_remove=["ALA"])


def test_warns_about_names_that_match_nothing(caplog: pytest.LogCaptureFixture) -> None:
    top, positions = _build(add_methanol_residue)
    with caplog.at_level(logging.WARNING):
        clean_topology(top, positions, keep=["ZZZ"])
    assert "'ZZZ'" in caplog.text
    assert "spelling and case" in caplog.text


def test_all_flags_off_returns_the_inputs_unchanged() -> None:
    top, positions = _build(add_methanol_residue, lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0)))
    new_top, new_positions, result = clean_topology(
        top,
        positions,
        remove_water=False,
        remove_ions=False,
        remove_additives=False,
    )
    # Identity, not equality: nothing was copied.
    assert new_top is top
    assert new_positions is positions
    assert result.removed == {}
    assert result.n_atoms_removed == 0


def test_clean_structure_is_a_no_op_on_an_already_clean_topology() -> None:
    top, positions = _build(add_methanol_residue)
    new_top, new_positions, result = clean_topology(top, positions)
    assert new_top is top
    assert new_positions is positions
    assert result.n_atoms_before == result.n_atoms_after


def test_emptied_chains_are_dropped() -> None:
    top, positions = _build(add_methanol_residue, lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0)))
    assert top.getNumChains() == 2
    new_top, _, _ = clean_topology(top, positions)
    assert new_top.getNumChains() == 1


def test_ids_bonds_and_box_survive() -> None:
    top, positions = _build(add_methanol_residue, lambda t: add_water_residue(t, origin=(5.0, 0.0, 0.0)))
    box = unit.Quantity([Vec3(4.0, 0, 0), Vec3(0, 4.0, 0), Vec3(0, 0, 4.0)], unit.nanometer)
    top.setPeriodicBoxVectors(box)

    new_top, _, _ = clean_topology(top, positions)

    lig = next(new_top.residues())
    assert lig.chain.id == "A"
    assert sorted(atom.name for atom in lig.atoms()) == ["C1", "H1", "H2", "H3", "H4", "O1"]
    assert new_top.getNumBonds() == 5
    assert new_top.getPeriodicBoxVectors() is not None


def test_removing_everything_warns_and_still_reports(caplog: pytest.LogCaptureFixture) -> None:
    top, positions = _build(add_water_residue)
    with caplog.at_level(logging.WARNING):
        new_top, _, result = clean_topology(top, positions)
    assert new_top.getNumAtoms() == 0
    assert result.removed == {"HOH": ("water", 1)}
    assert "nothing but solvent" in caplog.text


def test_mismatched_positions_are_rejected() -> None:
    top, positions = _build(add_methanol_residue)
    with pytest.raises(ValueError, match="positions has 3 entries"):
        clean_topology(top, positions[:3])


# -- clean_pdb -------------------------------------------------------------


def test_clean_pdb_roundtrip(tmp_path: Path) -> None:
    dirty = write_dirty_pdb(tmp_path / "dirty.pdb")
    out = tmp_path / "clean.pdb"

    result = clean_pdb(dirty, out)

    assert result.output_pdb == str(out)
    assert sorted(result.removed) == ["GOL", "HOH", "NA"]
    assert "CA" in result.retained

    reread = app.PDBFile(str(out))
    assert _names(reread.topology) == ["CA", "LIG"]
    assert reread.topology.getNumAtoms() == result.n_atoms_after
    # The ligand's CONECT records must survive, or antechamber loses the bonds.
    assert reread.topology.getNumBonds() == 5


def test_clean_pdb_creates_the_parent_directory(tmp_path: Path) -> None:
    dirty = write_dirty_pdb(tmp_path / "dirty.pdb")
    out = tmp_path / "nested" / "deeper" / "clean.pdb"
    clean_pdb(dirty, out)
    assert out.is_file()


def test_clean_pdb_refuses_to_overwrite_its_input(tmp_path: Path) -> None:
    dirty = write_dirty_pdb(tmp_path / "dirty.pdb")
    with pytest.raises(ValueError, match="same file"):
        clean_pdb(dirty, dirty)
    # And the input is still intact.
    assert app.PDBFile(str(dirty)).topology.getNumResidues() == 7


def test_clean_pdb_keeps_a_requested_additive(tmp_path: Path) -> None:
    dirty = write_dirty_pdb(tmp_path / "dirty.pdb")
    out = tmp_path / "clean.pdb"
    clean_pdb(dirty, out, keep=["GOL"])
    assert "GOL" in _names(app.PDBFile(str(out)).topology)


def test_module_exports_the_name_tables() -> None:
    # A user must be able to inspect exactly what will be deleted before running.
    for name in clean_structure.__all__:
        assert hasattr(clean_structure, name), name
    assert "GOL" in clean_structure.ADDITIVE_RESIDUES
    assert "CA" in clean_structure.STRUCTURAL_METAL_RESIDUES
