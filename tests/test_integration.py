"""End-to-end integration tests for :func:`forcefill.build_forcefield_xml`.

These run antechamber/parmchk2 for real and are skipped when AmberTools
(or its GAFF data files) is not available.
"""

import shutil
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from openmm import Vec3, app, unit
from openmm.app import element

from forcefill import build_forcefield_xml
from forcefill.nonstandard_ffxml import locate_gaff_dat


def _ambertools_available():
    if shutil.which("antechamber") is None or shutil.which("parmchk2") is None:
        return False
    try:
        locate_gaff_dat("gaff2")
    except FileNotFoundError:
        return False
    return True


pytestmark = pytest.mark.skipif(
    not _ambertools_available(), reason="AmberTools is not on PATH"
)


def _write_test_pdb(path, broken_gly=False):
    """Methanol as hetero residue LIG; optionally plus a hydrogen-stripped
    free glycine (standard name, missing atoms -> gets skipped)."""
    top = app.Topology()
    chain = top.addChain("A")
    res = top.addResidue("LIG", chain)
    atoms = {}
    for name, elem in [("C1", element.carbon), ("O1", element.oxygen),
                       ("H1", element.hydrogen), ("H2", element.hydrogen),
                       ("H3", element.hydrogen), ("H4", element.hydrogen)]:
        atoms[name] = top.addAtom(name, elem, res)
    for a, b in [("C1", "O1"), ("C1", "H1"), ("C1", "H2"), ("C1", "H3"),
                 ("O1", "H4")]:
        top.addBond(atoms[a], atoms[b])
    xyz = [
        (0.000, 0.000, 0.000),
        (1.410, 0.000, 0.000),
        (-0.360, 1.030, 0.000),
        (-0.360, -0.515, 0.892),
        (-0.360, -0.515, -0.892),
        (1.730, 0.890, 0.000),
    ]
    if broken_gly:
        chain2 = top.addChain("B")
        gly = top.addResidue("GLY", chain2)
        g = {}
        for name, elem in [("N", element.nitrogen), ("CA", element.carbon),
                           ("C", element.carbon), ("O", element.oxygen)]:
            g[name] = top.addAtom(name, elem, gly)
        top.addBond(g["N"], g["CA"])
        top.addBond(g["CA"], g["C"])
        top.addBond(g["C"], g["O"])
        xyz += [
            (8.000, 8.000, 8.000),
            (9.450, 8.000, 8.000),
            (10.050, 9.350, 8.000),
            (11.280, 9.400, 8.000),
        ]
    positions = unit.Quantity([Vec3(*p) for p in xyz], unit.angstrom)
    with open(path, "w") as fh:
        app.PDBFile.writeFile(top, positions, fh)
    return path


def test_methanol_end_to_end(tmp_path):
    pdb = _write_test_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml",
                                  workdir=tmp_path / "wd")
    assert result.parameterized == ["LIG"]
    assert result.skipped == {}
    assert result.residue_xmls.keys() == {"LIG"}

    xml_text = Path(result.forcefield_xml).read_text()
    assert '<Residue name="LIG">' in xml_text

    # The combined XML must work alongside the base force field.
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml",
                        result.forcefield_xml)
    ff.createSystem(app.PDBFile(str(pdb)).topology)


def test_relative_workdir_and_output(tmp_path, monkeypatch):
    # Regression: relative paths used to break antechamber, whose cwd is
    # changed to the per-residue directory.
    monkeypatch.chdir(tmp_path)
    _write_test_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml("in.pdb", "extras.xml", workdir="wd")
    assert result.parameterized == ["LIG"]
    assert Path(result.workdir).is_absolute()
    assert (tmp_path / "extras.xml").is_file()


def test_skipped_residue_still_returns_result(tmp_path):
    # Regression: a skipped residue used to make default validation raise
    # after all the work was done, losing the result.
    pdb = _write_test_pdb(tmp_path / "in.pdb", broken_gly=True)
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml",
                                  workdir=tmp_path / "wd")
    assert result.parameterized == ["LIG"]
    assert "GLY" in result.skipped
    assert Path(result.forcefield_xml).is_file()
