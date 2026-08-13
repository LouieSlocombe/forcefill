"""End-to-end integration tests for :func:`forcefill.build_forcefield_xml`.

These run antechamber/parmchk2 for real and are skipped when AmberTools
(or its GAFF data files) is not available.
"""

import math
import shutil
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from openmm import app

from forcefill import build_forcefield_xml
from forcefill.nonstandard_ffxml import locate_gaff_dat
from tests.helpers import write_methanol_pdb, write_methanol_sdf


def _ambertools_available():
    if shutil.which("antechamber") is None or shutil.which("parmchk2") is None:
        return False
    try:
        locate_gaff_dat("gaff2")
    except FileNotFoundError:
        return False
    return True


pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not _ambertools_available(), reason="AmberTools is not on PATH"),
]


def test_methanol_end_to_end(tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml", workdir=tmp_path / "wd")
    assert result.parameterized == ["LIG"]
    assert result.skipped == {}
    assert result.residue_xmls.keys() == {"LIG"}

    xml_text = Path(result.forcefield_xml).read_text()
    assert '<Residue name="LIG">' in xml_text

    # The combined XML must work alongside the base force field.
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    ff.createSystem(app.PDBFile(str(pdb)).topology)


def test_relative_workdir_and_output(tmp_path, monkeypatch):
    # Regression: relative paths used to break antechamber, whose cwd is
    # changed to the per-residue directory.
    monkeypatch.chdir(tmp_path)
    write_methanol_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml("in.pdb", "extras.xml", workdir="wd")
    assert result.parameterized == ["LIG"]
    assert Path(result.workdir).is_absolute()
    assert (tmp_path / "extras.xml").is_file()


def test_methanol_from_sdf_end_to_end(tmp_path):
    # residue_files: antechamber reads the drawn SDF instead of the extracted PDB.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml", workdir=tmp_path / "wd", residue_files={"LIG": sdf})
    assert result.parameterized == ["LIG"]
    assert not (tmp_path / "wd" / "LIG" / "LIG.pdb").exists()

    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    ff.createSystem(app.PDBFile(str(pdb)).topology)


def test_cleanup_removes_workdir_on_success(tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml", workdir=tmp_path / "wd", cleanup=True)
    assert result.workdir is None
    assert result.residue_xmls == {}
    assert not (tmp_path / "wd").exists()
    assert Path(result.forcefield_xml).is_file()


def test_cleanup_refuses_output_inside_workdir(tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    with pytest.raises(ValueError, match="inside the working directory"):
        build_forcefield_xml(pdb, tmp_path / "wd" / "extras.xml", workdir=tmp_path / "wd", cleanup=True)


def test_minimize_end_to_end(tmp_path):
    # The full path with real AmberTools parameters: antechamber charges and
    # parmchk2 constants have to produce a finite energy that a minimizer can
    # lower, not merely a System that builds.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml", workdir=tmp_path / "wd", minimize=True)

    ligand = result.minimizations["LIG"]
    assert ligand.n_atoms == 6
    assert math.isfinite(ligand.initial_energy)
    assert ligand.energy_change < 0

    assert result.full_minimization is not None
    assert math.isfinite(result.full_minimization.final_energy)


def test_skipped_residue_still_returns_result(tmp_path):
    # Regression: a skipped residue used to make default validation raise
    # after all the work was done, losing the result.
    pdb = write_methanol_pdb(tmp_path / "in.pdb", broken_gly=True)
    result = build_forcefield_xml(pdb, tmp_path / "extras.xml", workdir=tmp_path / "wd")
    assert result.parameterized == ["LIG"]
    assert "GLY" in result.skipped
    assert Path(result.forcefield_xml).is_file()
