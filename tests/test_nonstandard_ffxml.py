"""Tests for :mod:`forcefill.nonstandard_ffxml`.

Hermetic: no AmberTools executables and no PDB fixtures are required.
openmm and parmed must be importable (they are package dependencies);
everything here is skipped cleanly when they are not.
"""

import shutil

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

import forcefill
from forcefill import nonstandard_ffxml


class StubResidue:
    """Duck-types the Residue surface used by ``_classify_unmatched``:
    ``.name``, ``.atoms()`` and ``.external_bonds()``."""

    def __init__(self, name, n_atoms=2, external_bonds=0):
        self.name = name
        self._atoms = [object() for _ in range(n_atoms)]
        self._external = [object() for _ in range(external_bonds)]

    def atoms(self):
        return iter(self._atoms)

    def external_bonds(self):
        return iter(self._external)


def test_public_api_resolves():
    for name in forcefill.__all__:
        assert hasattr(forcefill, name), name
    assert forcefill.__all__ == nonstandard_ffxml.__all__
    assert forcefill.build_forcefield_xml is nonstandard_ffxml.build_forcefield_xml
    assert isinstance(forcefill.__version__, str)


def test_parameterization_result_defaults():
    result = forcefill.ParameterizationResult(forcefield_xml=None)
    assert result.residue_xmls == {}
    assert result.parameterized == []
    assert result.skipped == {}
    assert result.workdir is None


def test_classify_ligand_is_parameterized():
    to_param, skipped = nonstandard_ffxml._classify_unmatched(
        [StubResidue("LIG", n_atoms=10)]
    )
    assert list(to_param) == ["LIG"]
    assert skipped == {}


def test_classify_standard_residue_is_skipped():
    to_param, skipped = nonstandard_ffxml._classify_unmatched(
        [StubResidue("ALA", n_atoms=5)]
    )
    assert to_param == {}
    assert "repair the structure" in skipped["ALA"]


def test_classify_monatomic_is_skipped():
    to_param, skipped = nonstandard_ffxml._classify_unmatched(
        [StubResidue("ZN", n_atoms=1)]
    )
    assert to_param == {}
    assert "monatomic" in skipped["ZN"]


def test_classify_covalently_linked_is_skipped():
    to_param, skipped = nonstandard_ffxml._classify_unmatched(
        [StubResidue("PTM", n_atoms=8, external_bonds=1)]
    )
    assert to_param == {}
    assert "covalently bonded" in skipped["PTM"]


def test_classify_picks_most_complete_copy():
    small = StubResidue("LIG", n_atoms=4)
    big = StubResidue("LIG", n_atoms=9)
    to_param, skipped = nonstandard_ffxml._classify_unmatched([small, big])
    assert to_param["LIG"] is big
    assert skipped == {}


def test_require_executable_raises_when_missing(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="AmberTools"):
        nonstandard_ffxml._require_executable("antechamber")
