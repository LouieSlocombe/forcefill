"""Tests for :mod:`forcefill.checks`: validating and minimizing a generated force field.

Hermetic: the force fields come from the committed AmberTools fixtures and from
hand-written XML, so no executable is needed.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from openmm import app

from forcefill import amber, checks
from tests.helpers import DATA, METHANOL_XYZ, methanol_positions, methanol_residue, write_methanol_pdb

_LIG_TEMPLATE_XML = """<ForceField>
 <AtomTypes>
  <Type name="XC" class="XC" element="C" mass="12.011"/>
  <Type name="XO" class="XO" element="O" mass="15.999"/>
  <Type name="XH" class="XH" element="H" mass="1.008"/>
 </AtomTypes>
 <Residues>
  <Residue name="LIG">
   <Atom name="C1" type="XC" charge="0.1"/>
   <Atom name="O1" type="XO" charge="-0.5"/>
   <Atom name="H1" type="XH" charge="0.1"/>
   <Atom name="H2" type="XH" charge="0.1"/>
   <Atom name="H3" type="XH" charge="0.1"/>
   {h4_atom}
   <Bond atomName1="C1" atomName2="O1"/>
   <Bond atomName1="C1" atomName2="H1"/>
   <Bond atomName1="C1" atomName2="H2"/>
   <Bond atomName1="C1" atomName2="H3"/>
   {h4_bond}
  </Residue>
 </Residues>
</ForceField>
"""


def _write_lig_xml(tmp_path: Path, complete: bool = True) -> tuple[Path, app.ForceField]:
    xml = tmp_path / "lig.xml"
    xml.write_text(
        _LIG_TEMPLATE_XML.format(
            h4_atom='<Atom name="H4" type="XH" charge="0.1"/>' if complete else "",
            h4_bond='<Bond atomName1="O1" atomName2="H4"/>' if complete else "",
        )
    )
    return xml, app.ForceField(str(xml))


def _methanol_ffxml(tmp_path: Path) -> str:
    """A real, fully parameterized LIG force field, assembled from the committed fixtures."""
    return amber.assemble_openmm_ffxml(
        {"LIG": DATA / "methanol.mol2"}, [DATA / "methanol.frcmod"], tmp_path / "methanol_ff.xml"
    )


def test_minimization_result_energy_change() -> None:
    result = checks.MinimizationResult(n_atoms=6, initial_energy=22.0, final_energy=19.0, max_force=6.0)
    assert result.energy_change == pytest.approx(-3.0)


# -- per-residue validation ------------------------------------------------


def test_validate_parameterized_residues_ok(tmp_path: Path) -> None:
    xml, ff = _write_lig_xml(tmp_path)
    checks._validate_parameterized_residues({"LIG": methanol_residue()}, ff, [str(xml)])


def test_validate_parameterized_residues_detects_mismatch(tmp_path: Path) -> None:
    # Template lacks the hydroxyl hydrogen -> graph mismatch -> no template
    # matches the actual residue.
    xml, ff = _write_lig_xml(tmp_path, complete=False)
    with pytest.raises(RuntimeError, match="residue LIG"):
        checks._validate_parameterized_residues({"LIG": methanol_residue()}, ff, [str(xml)])


def test_validate_forcefield_xml_accepts_prebuilt_forcefield(tmp_path: Path) -> None:
    xml, ff = _write_lig_xml(tmp_path)
    residue = methanol_residue()
    checks.validate_forcefield_xml(residue.chain.topology, xml, base_forcefield=(), forcefield=ff)


def test_validate_forcefield_xml_reports_failure(tmp_path: Path) -> None:
    xml, _ff = _write_lig_xml(tmp_path, complete=False)
    residue = methanol_residue()
    with pytest.raises(RuntimeError, match="Validation failed"):
        checks.validate_forcefield_xml(residue.chain.topology, xml, base_forcefield=())


# -- minimization ----------------------------------------------------------
#
# These use the committed AmberTools fixtures rather than _write_lig_xml: the
# hand-written template carries no bonded or nonbonded parameters, so every
# energy would be a trivially finite zero.


def test_minimize_reports_energies_and_forces(tmp_path: Path) -> None:
    residue = methanol_residue()
    result = checks.minimize_with_forcefield_xml(
        residue.chain.topology, methanol_positions(), _methanol_ffxml(tmp_path), base_forcefield=()
    )
    assert result.n_atoms == 6
    assert math.isfinite(result.initial_energy)
    assert math.isfinite(result.final_energy)
    # Relations only, never exact energies: the CPU platform is not reproducible
    # run to run once a system is large enough to split across threads.
    assert result.final_energy < result.initial_energy
    assert result.energy_change == pytest.approx(result.final_energy - result.initial_energy)
    assert result.max_force > 0


def test_minimize_detects_nonfinite_energy(tmp_path: Path) -> None:
    # Superpose O1 on C1: the H-C1-O1 angle term is then undefined and the
    # potential energy comes back NaN, which createSystem alone never notices.
    collapsed = list(METHANOL_XYZ)
    collapsed[1] = collapsed[0]
    residue = methanol_residue()
    with pytest.raises(RuntimeError, match="Minimization failed") as excinfo:
        checks.minimize_with_forcefield_xml(
            residue.chain.topology, methanol_positions(collapsed), _methanol_ffxml(tmp_path), base_forcefield=()
        )
    assert "before minimizing" in str(excinfo.value)


def test_minimize_reports_template_mismatch(tmp_path: Path) -> None:
    xml, _ff = _write_lig_xml(tmp_path, complete=False)  # template lacks the hydroxyl hydrogen
    residue = methanol_residue()
    with pytest.raises(RuntimeError, match="Minimization failed"):
        checks.minimize_with_forcefield_xml(residue.chain.topology, methanol_positions(), xml, base_forcefield=())


def test_minimize_names_the_topology_it_failed_on(tmp_path: Path) -> None:
    # A multi-residue topology is named by size rather than by residue, and the
    # untemplated GLY is what makes it fail.
    pdb = app.PDBFile(str(write_methanol_pdb(tmp_path / "in.pdb", broken_gly=True)))
    with pytest.raises(RuntimeError, match="Minimization failed") as excinfo:
        checks.minimize_with_forcefield_xml(pdb.topology, pdb.positions, _methanol_ffxml(tmp_path), base_forcefield=())
    assert "the topology (10 atoms, 2 residues)" in str(excinfo.value)


def test_minimize_reports_unknown_platform(tmp_path: Path) -> None:
    residue = methanol_residue()
    with pytest.raises(RuntimeError, match="Minimization failed") as excinfo:
        checks.minimize_with_forcefield_xml(
            residue.chain.topology,
            methanol_positions(),
            _methanol_ffxml(tmp_path),
            base_forcefield=(),
            platform_name="Bogus",
        )
    assert "Bogus" in str(excinfo.value)


def test_minimize_rejects_settings_that_would_silently_do_nothing(tmp_path: Path) -> None:
    # OpenMM accepts a negative tolerance and then minimizes nothing at all,
    # which would leave this check reporting success without having run.
    residue = methanol_residue()
    args = (residue.chain.topology, methanol_positions(), _methanol_ffxml(tmp_path))
    with pytest.raises(ValueError, match="tolerance"):
        checks.minimize_with_forcefield_xml(*args, base_forcefield=(), tolerance=-1.0)
    with pytest.raises(ValueError, match="max_iterations"):
        checks.minimize_with_forcefield_xml(*args, base_forcefield=(), max_iterations=-5)


def test_minimize_accepts_prebuilt_forcefield(tmp_path: Path) -> None:
    xml = _methanol_ffxml(tmp_path)
    residue = methanol_residue()
    result = checks.minimize_with_forcefield_xml(
        residue.chain.topology, methanol_positions(), xml, base_forcefield=(), forcefield=app.ForceField(xml)
    )
    assert math.isfinite(result.final_energy)
