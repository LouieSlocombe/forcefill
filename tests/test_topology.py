"""Tests for :mod:`forcefill.topology`: finding, classifying and extracting residues.

Hermetic: no AmberTools executables and no PDB fixtures are required.
"""

import logging
from pathlib import Path

import pytest

pytest.importorskip("openmm")

from openmm import Vec3, app, unit

from forcefill import topology
from tests.helpers import METHANOL_ATOMS, METHANOL_XYZ, methanol_positions, methanol_residue


class StubResidue:
    """Duck-types the Residue surface ``_classify_unmatched`` uses: ``.name``, ``.atoms()``, ``.external_bonds()``."""

    def __init__(self, name, n_atoms=2, external_bonds=0):
        self.name = name
        self._atoms = [object() for _ in range(n_atoms)]
        self._external = [object() for _ in range(external_bonds)]

    def atoms(self):
        return iter(self._atoms)

    def external_bonds(self):
        return iter(self._external)


# -- classification --------------------------------------------------------


def test_classify_ligand_is_parameterized():
    to_param, skipped = topology._classify_unmatched([StubResidue("LIG", n_atoms=10)])
    assert list(to_param) == ["LIG"]
    assert skipped == {}


def test_classify_standard_residue_is_skipped():
    to_param, skipped = topology._classify_unmatched([StubResidue("ALA", n_atoms=5)])
    assert to_param == {}
    assert "repair the structure" in skipped["ALA"]


def test_classify_monatomic_is_skipped():
    to_param, skipped = topology._classify_unmatched([StubResidue("ZN", n_atoms=1)])
    assert to_param == {}
    assert "monatomic" in skipped["ZN"]


def test_classify_covalently_linked_is_skipped():
    to_param, skipped = topology._classify_unmatched([StubResidue("PTM", n_atoms=8, external_bonds=1)])
    assert to_param == {}
    assert "covalently bonded" in skipped["PTM"]


def test_classify_picks_most_complete_copy():
    small = StubResidue("LIG", n_atoms=4)
    big = StubResidue("LIG", n_atoms=9)
    to_param, skipped = topology._classify_unmatched([small, big])
    assert to_param["LIG"] is big
    assert skipped == {}


def test_classify_skips_if_any_copy_linked():
    # The representative (most atoms) is free-standing, but another copy is
    # covalently linked: the name must still be skipped.
    free_big = StubResidue("SUG", n_atoms=12)
    linked_small = StubResidue("SUG", n_atoms=11, external_bonds=1)
    to_param, skipped = topology._classify_unmatched([free_big, linked_small])
    assert to_param == {}
    assert "covalently bonded" in skipped["SUG"]
    assert "1 of 2 copies" in skipped["SUG"]


def test_warn_unused_overrides(caplog):
    to_param = {"LIG": StubResidue("LIG")}
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


def test_extract_residue_to_pdb_roundtrip(tmp_path):
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


def test_extract_residue_warns_on_missing_element(tmp_path, caplog):
    top = app.Topology()
    chain = top.addChain("A")
    res = top.addResidue("UNK", chain)
    top.addAtom("X1", None, res)
    positions = unit.Quantity([Vec3(0.0, 0.0, 0.0)], unit.angstrom)
    with caplog.at_level(logging.WARNING):
        topology.extract_residue_to_pdb(positions, next(top.residues()), tmp_path / "unk.pdb")
    assert "no element" in caplog.text


def test_describe_topology_names_a_lone_residue():
    assert topology._describe_topology(methanol_residue().chain.topology) == "residue LIG (6 atoms)"
