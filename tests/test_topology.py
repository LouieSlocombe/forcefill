"""Tests for :mod:`forcefill.topology`: finding, classifying and extracting residues.

Hermetic: no AmberTools executables and no PDB fixtures are required.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path

import pytest

pytest.importorskip("openmm")

from openmm import Vec3, app, unit

from forcefill import clean_topology, topology
from tests.helpers import METHANOL_ATOMS, METHANOL_XYZ, methanol_positions, methanol_residue


def _residue(name: str, n_atoms: int = 2, external_bonds: int = 0) -> app.topology.Residue:
    """Build a residue and, optionally, bonds to atoms in a neighbouring residue."""
    top = app.Topology()
    chain = top.addChain()
    residue = top.addResidue(name, chain)
    atoms = [top.addAtom(f"C{i}", app.element.carbon, residue) for i in range(n_atoms)]
    if external_bonds:
        neighbour = top.addResidue("PRT", chain)
        for i in range(external_bonds):
            top.addBond(atoms[0], top.addAtom(f"N{i}", app.element.nitrogen, neighbour))
    return residue


# -- classification --------------------------------------------------------


def test_classify_ligand_is_parameterized() -> None:
    to_param, skipped = topology._classify_unmatched([_residue("LIG", n_atoms=10)])
    assert list(to_param) == ["LIG"]
    assert skipped == {}


def test_classify_standard_residue_is_skipped() -> None:
    to_param, skipped = topology._classify_unmatched([_residue("ALA", n_atoms=5)])
    assert to_param == {}
    assert "repair the structure" in skipped["ALA"]


@pytest.mark.parametrize("name", ["HSD", "HSE", "HSP", "ADE", "CYT", "GUA", "THY", "URA", "TIP3"])
def test_classify_charmm_spellings_are_standard(name: str) -> None:
    # A structure prepared for CHARMM names histidine HSD/HSE/HSP and the bases
    # ADE/CYT/GUA/THY/URA. Without those, an incomplete protein residue would be
    # sent to a backend instead of reported as needing repair.
    to_param, skipped = topology._classify_unmatched([_residue(name, n_atoms=8)])
    assert to_param == {}
    assert "repair the structure" in skipped[name]


def test_classify_monatomic_is_skipped() -> None:
    to_param, skipped = topology._classify_unmatched([_residue("ZN", n_atoms=1)])
    assert to_param == {}
    assert "monatomic" in skipped["ZN"]


def test_classify_covalently_linked_is_skipped() -> None:
    to_param, skipped = topology._classify_unmatched([_residue("PTM", n_atoms=8, external_bonds=1)])
    assert to_param == {}
    assert "covalently bonded" in skipped["PTM"]


def test_classify_picks_most_complete_copy() -> None:
    small = _residue("LIG", n_atoms=4)
    big = _residue("LIG", n_atoms=9)
    to_param, skipped = topology._classify_unmatched([small, big])
    assert to_param["LIG"] is big
    assert skipped == {}


def test_classify_skips_if_any_copy_linked() -> None:
    # The representative (most atoms) is free-standing, but another copy is
    # covalently linked: the name must still be skipped.
    free_big = _residue("SUG", n_atoms=12)
    linked_small = _residue("SUG", n_atoms=11, external_bonds=1)
    to_param, skipped = topology._classify_unmatched([free_big, linked_small])
    assert to_param == {}
    assert "covalently bonded" in skipped["SUG"]
    assert "1 of 2 copies" in skipped["SUG"]


def test_external_bond_counts_ignore_internal_bonds_and_count_both_ends() -> None:
    top = app.Topology()
    chain = top.addChain()
    first, second, third, isolated = [top.addResidue(name, chain) for name in ("LIG", "LIG", "OTH", "FREE")]
    atoms = [[top.addAtom(f"C{i}", app.element.carbon, res) for i in range(2)] for res in top.residues()]
    for pair in atoms:
        top.addBond(*pair)
    top.addBond(atoms[0][0], atoms[1][0])
    top.addBond(atoms[0][1], atoms[1][1])
    top.addBond(atoms[1][0], atoms[2][0])

    counts = topology._external_bond_counts(top)

    assert counts == {first: 2, second: 3, third: 1}
    assert isolated not in counts
    assert topology._external_bond_counts(app.Topology()) == {}


def test_classify_does_not_confuse_indices_from_different_topologies() -> None:
    linked = _residue("LINK", external_bonds=2)
    free = _residue("FREE")
    assert linked.index == free.index

    to_param, skipped = topology._classify_unmatched([linked, free])

    assert to_param == {"FREE": free}
    assert list(skipped) == ["LINK"]
    assert "1 of 1 copies" in skipped["LINK"]


def test_classify_scans_each_distinct_topology_once(monkeypatch: pytest.MonkeyPatch) -> None:
    first = _residue("LIG", external_bonds=1)
    top = first.chain.topology
    neighbour = list(top.residues())[1]
    original_bonds = top.bonds
    traversals = 0

    def counted_bonds() -> Iterator[app.topology.Bond]:
        nonlocal traversals
        traversals += 1
        return original_bonds()

    monkeypatch.setattr(top, "bonds", counted_bonds)
    topology._classify_unmatched([first, neighbour])

    assert traversals == 1


def test_cleaning_keeps_linked_copies_while_parameterization_skips_the_name() -> None:
    linked = _residue("GOL", external_bonds=2)
    top = linked.chain.topology
    free = top.addResidue("GOL", linked.chain)
    for i in range(2):
        top.addAtom(f"C{i}", app.element.carbon, free)
    positions = unit.Quantity([Vec3(i, 0, 0) for i in range(top.getNumAtoms())], unit.angstrom)

    _, _, cleaning = clean_topology(top, positions)
    to_param, skipped = topology._classify_unmatched([linked, free])

    assert cleaning.removed == {"GOL": ("additive", 1)}
    assert "2 neighbouring atom(s)" in cleaning.retained["GOL"]
    assert to_param == {}
    assert "1 of 2 copies" in skipped["GOL"]


def test_warn_unused_overrides(caplog: pytest.LogCaptureFixture) -> None:
    to_param = {"LIG": _residue("LIG")}
    skipped = {"ZN": "monatomic species - ..."}
    with caplog.at_level(logging.WARNING):
        topology._warn_unused_overrides(
            to_param,
            skipped,
            net_charges={"lig": -1, "ZN": 2, "LIG": 0},
            multiplicities={"XYZ": 3},
            residue_files={"NOPE": "lig.sdf", "LIG": "ok.sdf"},
        )
    assert len(caplog.records) == 4  # 'lig', 'ZN', 'XYZ', 'NOPE'; 'LIG' is fine
    text = caplog.text
    assert "'lig'" in text and "spelling" in text
    assert "'ZN'" in text and "skipped" in text
    assert "'XYZ'" in text
    assert "residue_files['NOPE']" in text


# -- extraction ------------------------------------------------------------


def test_extract_residue_to_pdb_roundtrip(tmp_path: Path) -> None:
    residue = methanol_residue()
    out = tmp_path / "LIG.pdb"
    path = topology.extract_residue_to_pdb(methanol_positions(), residue, out)
    assert Path(path) == out

    reread = app.PDBFile(path)
    atoms = list(reread.topology.atoms())
    assert [a.name for a in atoms] == [name for name, _ in METHANOL_ATOMS]
    assert [a.element.symbol for a in atoms] == ["C", "O", "H", "H", "H", "H"]
    assert reread.topology.getNumBonds() == 5  # CONECT records survive
    new = reread.positions.value_in_unit(unit.angstrom)
    for (x, y, z), b in zip(METHANOL_XYZ, new, strict=True):
        assert max(abs(x - b.x), abs(y - b.y), abs(z - b.z)) < 1e-2


def test_extract_residue_warns_on_missing_element(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    top = app.Topology()
    chain = top.addChain("A")
    res = top.addResidue("UNK", chain)
    top.addAtom("X1", None, res)
    positions = unit.Quantity([Vec3(0.0, 0.0, 0.0)], unit.angstrom)
    with caplog.at_level(logging.WARNING):
        topology.extract_residue_to_pdb(positions, next(top.residues()), tmp_path / "unk.pdb")
    assert "no element" in caplog.text


def test_describe_topology_names_a_lone_residue() -> None:
    assert topology._describe_topology(methanol_residue().chain.topology) == "residue LIG (6 atoms)"
