"""Tests for :mod:`forcefill.nonstandard_ffxml`.

Hermetic: no AmberTools executables and no PDB fixtures are required.
openmm and parmed must be importable (they are package dependencies);
everything here is skipped cleanly when they are not.
"""

import logging
import math
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from openmm import Vec3, app, unit

import forcefill
from forcefill import nonstandard_ffxml
from tests.helpers import (
    METHANOL_ATOMS,
    METHANOL_XYZ,
    write_broken_gly_pdb,
    write_methanol_pdb,
    write_methanol_sdf,
    write_water_pdb,
)

DATA = Path(__file__).parent / "data"


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
    assert result.minimizations == {}
    assert result.full_minimization is None


def test_minimization_result_energy_change():
    result = forcefill.MinimizationResult(n_atoms=6, initial_energy=22.0, final_energy=19.0, max_force=6.0)
    assert result.energy_change == pytest.approx(-3.0)


def test_classify_ligand_is_parameterized():
    to_param, skipped = nonstandard_ffxml._classify_unmatched([StubResidue("LIG", n_atoms=10)])
    assert list(to_param) == ["LIG"]
    assert skipped == {}


def test_classify_standard_residue_is_skipped():
    to_param, skipped = nonstandard_ffxml._classify_unmatched([StubResidue("ALA", n_atoms=5)])
    assert to_param == {}
    assert "repair the structure" in skipped["ALA"]


def test_classify_monatomic_is_skipped():
    to_param, skipped = nonstandard_ffxml._classify_unmatched([StubResidue("ZN", n_atoms=1)])
    assert to_param == {}
    assert "monatomic" in skipped["ZN"]


def test_classify_covalently_linked_is_skipped():
    to_param, skipped = nonstandard_ffxml._classify_unmatched([StubResidue("PTM", n_atoms=8, external_bonds=1)])
    assert to_param == {}
    assert "covalently bonded" in skipped["PTM"]


def test_classify_picks_most_complete_copy():
    small = StubResidue("LIG", n_atoms=4)
    big = StubResidue("LIG", n_atoms=9)
    to_param, skipped = nonstandard_ffxml._classify_unmatched([small, big])
    assert to_param["LIG"] is big
    assert skipped == {}


def test_classify_skips_if_any_copy_linked():
    # The representative (most atoms) is free-standing, but another copy is
    # covalently linked: the name must still be skipped.
    free_big = StubResidue("SUG", n_atoms=12)
    linked_small = StubResidue("SUG", n_atoms=11, external_bonds=1)
    to_param, skipped = nonstandard_ffxml._classify_unmatched([free_big, linked_small])
    assert to_param == {}
    assert "covalently bonded" in skipped["SUG"]
    assert "1 of 2 copies" in skipped["SUG"]


def test_require_executable_raises_when_missing(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="AmberTools"):
        nonstandard_ffxml._require_executable("antechamber")


class _RunRecorder:
    """Stands in for _run: records calls and creates the '-o' output file."""

    def __init__(self):
        self.calls = []
        self.kwargs = []

    def __call__(self, cmd, cwd, **kwargs):
        argv = [str(c) for c in cmd]
        self.calls.append((argv, Path(cwd)))
        self.kwargs.append(kwargs)
        out = Path(argv[argv.index("-o") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.touch()


def test_run_antechamber_resolves_relative_paths(monkeypatch, tmp_path):
    recorder = _RunRecorder()
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "antechamber")
    monkeypatch.setattr(nonstandard_ffxml, "_run", recorder)
    monkeypatch.chdir(tmp_path)

    nonstandard_ffxml.run_antechamber("wd/LIG/LIG.pdb", "wd/LIG/LIG.mol2", "LIG")

    ((cmd, cwd),) = recorder.calls
    in_arg = Path(cmd[cmd.index("-i") + 1])
    out_arg = Path(cmd[cmd.index("-o") + 1])
    assert in_arg == (tmp_path / "wd/LIG/LIG.pdb").resolve()
    assert out_arg == (tmp_path / "wd/LIG/LIG.mol2").resolve()
    assert cwd == (tmp_path / "wd/LIG").resolve()


def test_run_parmchk2_resolves_relative_paths(monkeypatch, tmp_path):
    recorder = _RunRecorder()
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "parmchk2")
    monkeypatch.setattr(nonstandard_ffxml, "_run", recorder)
    monkeypatch.chdir(tmp_path)

    nonstandard_ffxml.run_parmchk2("wd/LIG/LIG.mol2", "wd/LIG/LIG.frcmod")

    ((cmd, cwd),) = recorder.calls
    in_arg = Path(cmd[cmd.index("-i") + 1])
    out_arg = Path(cmd[cmd.index("-o") + 1])
    assert in_arg == (tmp_path / "wd/LIG/LIG.mol2").resolve()
    assert out_arg == (tmp_path / "wd/LIG/LIG.frcmod").resolve()
    assert cwd == (tmp_path / "wd/LIG").resolve()


def test_run_antechamber_purge_scratch_flag(monkeypatch, tmp_path):
    recorder = _RunRecorder()
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "antechamber")
    monkeypatch.setattr(nonstandard_ffxml, "_run", recorder)

    nonstandard_ffxml.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG")
    nonstandard_ffxml.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG", purge_scratch=False)

    (cmd_default, _), (cmd_keep, _) = recorder.calls
    assert cmd_default[cmd_default.index("-pf") + 1] == "y"
    assert cmd_keep[cmd_keep.index("-pf") + 1] == "n"


def test_run_antechamber_infers_input_format(monkeypatch, tmp_path):
    recorder = _RunRecorder()
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "antechamber")
    monkeypatch.setattr(nonstandard_ffxml, "_run", recorder)

    nonstandard_ffxml.run_antechamber(tmp_path / "lig.sdf", tmp_path / "a.mol2", "LIG")
    nonstandard_ffxml.run_antechamber(tmp_path / "lig.mol2", tmp_path / "b.mol2", "LIG")
    nonstandard_ffxml.run_antechamber(tmp_path / "lig.xyz", tmp_path / "c.mol2", "LIG", input_format="mol2")

    formats = [cmd[cmd.index("-fi") + 1] for cmd, _ in recorder.calls]
    assert formats == ["sdf", "mol2", "mol2"]


def test_run_antechamber_unknown_suffix_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "antechamber")
    monkeypatch.setattr(nonstandard_ffxml, "_run", lambda *args, **kwargs: None)
    with pytest.raises(ValueError, match="input format"):
        nonstandard_ffxml.run_antechamber(tmp_path / "lig.xyz", tmp_path / "out.mol2", "LIG")


def test_run_antechamber_missing_output_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "antechamber")
    monkeypatch.setattr(nonstandard_ffxml, "_run", lambda cmd, cwd, **kw: None)  # writes nothing
    with pytest.raises(RuntimeError, match="did not write"):
        nonstandard_ffxml.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG")


def test_run_parmchk2_missing_output_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: "parmchk2")
    monkeypatch.setattr(nonstandard_ffxml, "_run", lambda cmd, cwd, **kw: None)
    with pytest.raises(RuntimeError, match="did not write"):
        nonstandard_ffxml.run_parmchk2(tmp_path / "in.mol2", tmp_path / "out.frcmod")


def test_run_nonzero_exit_includes_tails_and_hint(tmp_path):
    script = "import sys; print('OUT-MARKER'); print('ERR-MARKER', file=sys.stderr); sys.exit(3)"
    with pytest.raises(RuntimeError) as excinfo:
        nonstandard_ffxml._run([sys.executable, "-c", script], cwd=tmp_path, hint="HINT-TEXT")
    msg = str(excinfo.value)
    assert "exit code 3" in msg
    assert "OUT-MARKER" in msg
    assert "ERR-MARKER" in msg
    assert msg.endswith("HINT-TEXT")


def test_run_timeout(tmp_path):
    with pytest.raises(RuntimeError, match="timed out"):
        nonstandard_ffxml._run([sys.executable, "-c", "import time; time.sleep(30)"], cwd=tmp_path, timeout=0.2)


def test_bad_atom_type_rejected_early(tmp_path):
    with pytest.raises(ValueError, match="atom_type"):
        nonstandard_ffxml.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG", atom_type="gaff3")
    with pytest.raises(ValueError, match="atom_type"):
        nonstandard_ffxml.run_parmchk2(tmp_path / "in.mol2", tmp_path / "out.frcmod", atom_type="gaff3")
    with pytest.raises(ValueError, match="atom_type"):
        nonstandard_ffxml.locate_gaff_dat("gaff3")
    with pytest.raises(ValueError, match="atom_type"):
        forcefill.build_forcefield_xml(tmp_path / "absent.pdb", atom_type="gaff3")


def test_bad_charge_method_rejected_early(tmp_path):
    with pytest.raises(ValueError, match="charge_method"):
        nonstandard_ffxml.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG", charge_method="bbc")
    with pytest.raises(ValueError, match="charge_method"):
        forcefill.build_forcefield_xml(tmp_path / "absent.pdb", charge_method="bbc")


def test_warn_unused_overrides(caplog):
    to_param = {"LIG": StubResidue("LIG")}
    skipped = {"ZN": "monatomic species - ..."}
    with caplog.at_level(logging.WARNING):
        nonstandard_ffxml._warn_unused_overrides(
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


# -- per-residue validation ------------------------------------------------


def _methanol_residue():
    from openmm import app
    from openmm.app import element

    top = app.Topology()
    chain = top.addChain("A")
    res = top.addResidue("LIG", chain)
    atoms = {}
    for name, elem in [
        ("C1", element.carbon),
        ("O1", element.oxygen),
        ("H1", element.hydrogen),
        ("H2", element.hydrogen),
        ("H3", element.hydrogen),
        ("H4", element.hydrogen),
    ]:
        atoms[name] = top.addAtom(name, elem, res)
    for a, b in [("C1", "O1"), ("C1", "H1"), ("C1", "H2"), ("C1", "H3"), ("O1", "H4")]:
        top.addBond(atoms[a], atoms[b])
    return next(top.residues())


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


def _write_lig_xml(tmp_path, complete=True):
    from openmm import app

    xml = tmp_path / "lig.xml"
    xml.write_text(
        _LIG_TEMPLATE_XML.format(
            h4_atom='<Atom name="H4" type="XH" charge="0.1"/>' if complete else "",
            h4_bond='<Bond atomName1="O1" atomName2="H4"/>' if complete else "",
        )
    )
    return xml, app.ForceField(str(xml))


def test_validate_parameterized_residues_ok(tmp_path):
    xml, ff = _write_lig_xml(tmp_path)
    nonstandard_ffxml._validate_parameterized_residues({"LIG": _methanol_residue()}, ff, [str(xml)])


def test_validate_parameterized_residues_detects_mismatch(tmp_path):
    # Template lacks the hydroxyl hydrogen -> graph mismatch -> no template
    # matches the actual residue.
    xml, ff = _write_lig_xml(tmp_path, complete=False)
    with pytest.raises(RuntimeError, match="residue LIG"):
        nonstandard_ffxml._validate_parameterized_residues({"LIG": _methanol_residue()}, ff, [str(xml)])


def test_validate_forcefield_xml_accepts_prebuilt_forcefield(tmp_path):
    xml, ff = _write_lig_xml(tmp_path)
    residue = _methanol_residue()
    nonstandard_ffxml.validate_forcefield_xml(residue.chain.topology, xml, base_forcefield=(), forcefield=ff)


def test_validate_forcefield_xml_reports_failure(tmp_path):
    xml, _ff = _write_lig_xml(tmp_path, complete=False)
    residue = _methanol_residue()
    with pytest.raises(RuntimeError, match="Validation failed"):
        nonstandard_ffxml.validate_forcefield_xml(residue.chain.topology, xml, base_forcefield=())


# -- minimization ----------------------------------------------------------
#
# These use the committed AmberTools fixtures rather than _write_lig_xml: the
# hand-written template carries no bonded or nonbonded parameters, so every
# energy would be a trivially finite zero.


def _methanol_ffxml(tmp_path):
    """A real, fully parameterized LIG force field, assembled from the committed fixtures."""
    return nonstandard_ffxml.assemble_openmm_ffxml(
        {"LIG": DATA / "methanol.mol2"}, [DATA / "methanol.frcmod"], tmp_path / "methanol_ff.xml"
    )


def _methanol_positions(xyz=METHANOL_XYZ):
    return unit.Quantity([Vec3(*p) for p in xyz], unit.angstrom)


def test_minimize_reports_energies_and_forces(tmp_path):
    residue = _methanol_residue()
    result = nonstandard_ffxml.minimize_with_forcefield_xml(
        residue.chain.topology, _methanol_positions(), _methanol_ffxml(tmp_path), base_forcefield=()
    )
    assert result.n_atoms == 6
    assert math.isfinite(result.initial_energy)
    assert math.isfinite(result.final_energy)
    # Relations only, never exact energies: the CPU platform is not reproducible
    # run to run once a system is large enough to split across threads.
    assert result.final_energy < result.initial_energy
    assert result.energy_change == pytest.approx(result.final_energy - result.initial_energy)
    assert result.max_force > 0


def test_minimize_detects_nonfinite_energy(tmp_path):
    # Superpose O1 on C1: the H-C1-O1 angle term is then undefined and the
    # potential energy comes back NaN, which createSystem alone never notices.
    collapsed = list(METHANOL_XYZ)
    collapsed[1] = collapsed[0]
    residue = _methanol_residue()
    with pytest.raises(RuntimeError, match="Minimization failed") as excinfo:
        nonstandard_ffxml.minimize_with_forcefield_xml(
            residue.chain.topology, _methanol_positions(collapsed), _methanol_ffxml(tmp_path), base_forcefield=()
        )
    assert "before minimizing" in str(excinfo.value)


def test_minimize_reports_template_mismatch(tmp_path):
    xml, _ff = _write_lig_xml(tmp_path, complete=False)  # template lacks the hydroxyl hydrogen
    residue = _methanol_residue()
    with pytest.raises(RuntimeError, match="Minimization failed"):
        nonstandard_ffxml.minimize_with_forcefield_xml(
            residue.chain.topology, _methanol_positions(), xml, base_forcefield=()
        )


def test_minimize_names_the_topology_it_failed_on(tmp_path):
    # A multi-residue topology is named by size rather than by residue, and the
    # untemplated GLY is what makes it fail.
    pdb = app.PDBFile(str(write_methanol_pdb(tmp_path / "in.pdb", broken_gly=True)))
    with pytest.raises(RuntimeError, match="Minimization failed") as excinfo:
        nonstandard_ffxml.minimize_with_forcefield_xml(
            pdb.topology, pdb.positions, _methanol_ffxml(tmp_path), base_forcefield=()
        )
    assert "the topology (10 atoms, 2 residues)" in str(excinfo.value)


def test_minimize_reports_unknown_platform(tmp_path):
    residue = _methanol_residue()
    with pytest.raises(RuntimeError, match="Minimization failed") as excinfo:
        nonstandard_ffxml.minimize_with_forcefield_xml(
            residue.chain.topology,
            _methanol_positions(),
            _methanol_ffxml(tmp_path),
            base_forcefield=(),
            platform_name="Bogus",
        )
    assert "Bogus" in str(excinfo.value)


def test_minimize_rejects_settings_that_would_silently_do_nothing(tmp_path):
    # OpenMM accepts a negative tolerance and then minimizes nothing at all,
    # which would leave this check reporting success without having run.
    residue = _methanol_residue()
    args = (residue.chain.topology, _methanol_positions(), _methanol_ffxml(tmp_path))
    with pytest.raises(ValueError, match="tolerance"):
        nonstandard_ffxml.minimize_with_forcefield_xml(*args, base_forcefield=(), tolerance=-1.0)
    with pytest.raises(ValueError, match="max_iterations"):
        nonstandard_ffxml.minimize_with_forcefield_xml(*args, base_forcefield=(), max_iterations=-5)


def test_minimize_accepts_prebuilt_forcefield(tmp_path):
    xml = _methanol_ffxml(tmp_path)
    residue = _methanol_residue()
    result = nonstandard_ffxml.minimize_with_forcefield_xml(
        residue.chain.topology, _methanol_positions(), xml, base_forcefield=(), forcefield=app.ForceField(xml)
    )
    assert math.isfinite(result.final_energy)


# -- build_forcefield_xml orchestration (AmberTools faked) -----------------


@pytest.fixture
def fake_ambertools(monkeypatch):
    """Replace the AmberTools wrappers with fakes that install the committed fixtures."""
    calls = {"antechamber": [], "parmchk2": []}

    def fake_antechamber(input_pdb, output_mol2, residue_name, **kwargs):
        calls["antechamber"].append(
            {"input": str(input_pdb), "output": str(output_mol2), "residue": residue_name, **kwargs}
        )
        out = Path(output_mol2)
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(DATA / "methanol.mol2", out)
        return str(out)

    def fake_parmchk2(input_mol2, output_frcmod, atom_type="gaff2", timeout=None):
        calls["parmchk2"].append(
            {"input": str(input_mol2), "output": str(output_frcmod), "atom_type": atom_type, "timeout": timeout}
        )
        shutil.copyfile(DATA / "methanol.frcmod", output_frcmod)
        return str(output_frcmod)

    monkeypatch.setattr(nonstandard_ffxml, "_require_executable", lambda name: f"/fake/{name}")
    # The complete frcmod (parmchk2 -a Y) stands in for gaff2.dat.
    monkeypatch.setattr(nonstandard_ffxml, "locate_gaff_dat", lambda atom_type="gaff2": str(DATA / "methanol.frcmod"))
    monkeypatch.setattr(nonstandard_ffxml, "run_antechamber", fake_antechamber)
    monkeypatch.setattr(nonstandard_ffxml, "run_parmchk2", fake_parmchk2)
    return calls


def test_orchestration_end_to_end_with_fakes(fake_ambertools, tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    wd = tmp_path / "wd"
    result = forcefill.build_forcefield_xml(
        pdb,
        tmp_path / "extras.xml",
        base_forcefield=(),
        net_charges={"LIG": -1},
        multiplicities={"LIG": 3},
        antechamber_args=("-dr", "no"),
        workdir=wd,
        timeout=123,
    )

    assert result.parameterized == ["LIG"]
    assert result.skipped == {}
    assert result.forcefield_xml == str(tmp_path / "extras.xml")
    assert Path(result.forcefield_xml).is_file()
    assert Path(result.residue_xmls["LIG"]) == wd / "LIG" / "LIG.xml"
    assert (wd / "LIG" / "LIG.pdb").is_file()
    assert result.workdir == str(wd)
    # minimize is opt-in; the default run does no energy evaluation.
    assert result.minimizations == {}
    assert result.full_minimization is None

    (ante,) = fake_ambertools["antechamber"]
    assert ante["residue"] == "LIG"
    assert ante["net_charge"] == -1
    assert ante["multiplicity"] == 3
    assert ante["extra_args"] == ("-dr", "no")
    assert ante["timeout"] == 123
    (chk,) = fake_ambertools["parmchk2"]
    assert chk["atom_type"] == "gaff2"
    assert chk["timeout"] == 123


def test_orchestration_residue_files_bypass_extraction(fake_ambertools, tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    wd = tmp_path / "wd"
    result = forcefill.build_forcefield_xml(
        pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=wd, residue_files={"LIG": sdf}
    )
    assert result.parameterized == ["LIG"]
    (ante,) = fake_ambertools["antechamber"]
    assert ante["input"] == str(sdf)
    assert not (wd / "LIG" / "LIG.pdb").exists()  # extraction skipped


def test_orchestration_minimize(fake_ambertools, tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = forcefill.build_forcefield_xml(
        pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=tmp_path / "wd", minimize=True
    )
    assert result.minimizations.keys() == {"LIG"}
    assert result.minimizations["LIG"].n_atoms == 6
    assert math.isfinite(result.minimizations["LIG"].final_energy)
    # Nothing was skipped, so the whole input is minimized too. It is the same
    # six atoms here, but by a different route: the full topology, not a
    # per-residue subtopology.
    assert result.full_minimization is not None
    assert result.full_minimization.n_atoms == 6


def test_orchestration_minimize_without_validate(fake_ambertools, tmp_path):
    # minimize subsumes validate: it has to build the System to get an energy.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = forcefill.build_forcefield_xml(
        pdb, tmp_path / "extras.xml", base_forcefield=(), validate=False, minimize=True
    )
    try:
        assert result.minimizations.keys() == {"LIG"}
        assert result.full_minimization is not None
    finally:
        shutil.rmtree(result.workdir, ignore_errors=True)


def test_orchestration_minimize_skips_full_structure_when_residue_skipped(fake_ambertools, tmp_path, caplog):
    pdb = write_methanol_pdb(tmp_path / "in.pdb", broken_gly=True)
    with caplog.at_level(logging.WARNING):
        result = forcefill.build_forcefield_xml(
            pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=tmp_path / "wd", minimize=True
        )
    assert "GLY" in result.skipped
    # The parameterized residue is still minimized on its own; the full
    # structure cannot be, because GLY has no template.
    assert result.minimizations.keys() == {"LIG"}
    assert result.full_minimization is None
    assert "Skipping the full-structure checks" in caplog.text


def test_orchestration_cleanup_removes_workdir(fake_ambertools, tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    wd = tmp_path / "wd"
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=wd, cleanup=True)
    assert result.workdir is None
    assert result.residue_xmls == {}
    assert not wd.exists()
    assert Path(result.forcefield_xml).is_file()


def test_orchestration_failure_preserves_workdir(fake_ambertools, monkeypatch, tmp_path, caplog):
    def boom(*args, **kwargs):
        raise RuntimeError("antechamber exploded")

    monkeypatch.setattr(nonstandard_ffxml, "run_antechamber", boom)
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    wd = tmp_path / "wd"
    with caplog.at_level(logging.WARNING), pytest.raises(RuntimeError, match="exploded"):
        forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=wd, cleanup=True)
    assert wd.exists()
    assert "kept for debugging" in caplog.text


def test_orchestration_default_workdir_reported(fake_ambertools, tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), validate=False)
    try:
        assert result.workdir is not None
        assert Path(result.workdir).name.startswith("nonstandard_ff_")
        assert Path(result.residue_xmls["LIG"]).is_file()
    finally:
        shutil.rmtree(result.workdir, ignore_errors=True)


def test_nothing_to_parameterize_short_circuits(tmp_path):
    pdb = write_water_pdb(tmp_path / "w.pdb")
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=("amber14/tip3p.xml",))
    assert result.forcefield_xml is None
    assert result.parameterized == []
    assert not (tmp_path / "extras.xml").exists()


def test_everything_skipped_raises(tmp_path):
    pdb = write_broken_gly_pdb(tmp_path / "g.pdb")
    with pytest.raises(RuntimeError, match="none can be auto-parameterized"):
        forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=())


# -- ParmEd assembly against the committed fixtures ------------------------


def test_assemble_openmm_ffxml_semantic(tmp_path):
    out = tmp_path / "sub" / "lig.xml"  # exercises output-directory creation
    path = nonstandard_ffxml.assemble_openmm_ffxml({"LIG": DATA / "methanol.mol2"}, [DATA / "methanol.frcmod"], out)
    assert Path(path) == out

    tree = ET.parse(path)
    residue = tree.find(".//Residues/Residue[@name='LIG']")
    assert residue is not None
    atoms = residue.findall("Atom")
    assert len(atoms) == 6
    assert abs(sum(float(a.get("charge")) for a in atoms)) < 1e-6  # methanol is neutral
    assert len(residue.findall("Bond")) == 5
    for section in ("AtomTypes", "HarmonicBondForce", "HarmonicAngleForce", "NonbondedForce"):
        assert tree.find(f".//{section}") is not None, section


def test_load_residue_template_rejects_multi_residue_mol2(tmp_path):
    multi = tmp_path / "two.mol2"
    multi.write_text(
        "@<TRIPOS>MOLECULE\n"
        "TWO\n"
        "4 2 2 0 0\n"
        "SMALL\n"
        "NO_CHARGES\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        "      1 O1     0.0000  0.0000  0.0000 oh  1 R1  0.0000\n"
        "      2 H1     0.9600  0.0000  0.0000 ho  1 R1  0.0000\n"
        "      3 O2     5.0000  0.0000  0.0000 oh  2 R2  0.0000\n"
        "      4 H2     5.9600  0.0000  0.0000 ho  2 R2  0.0000\n"
        "@<TRIPOS>BOND\n"
        "     1     1     2 1\n"
        "     2     3     4 1\n"
        "@<TRIPOS>SUBSTRUCTURE\n"
        "     1 R1     1 TEMP  0 **** ****  0 ROOT\n"
        "     2 R2     3 TEMP  0 **** ****  0 ROOT\n"
    )
    with pytest.raises(ValueError, match="contains 2 residues"):
        nonstandard_ffxml._load_residue_template(multi, "X")


def test_load_residue_template_rejects_non_template(tmp_path):
    with pytest.raises(TypeError, match="ResidueTemplate"):
        nonstandard_ffxml._load_residue_template(DATA / "methanol.frcmod", "X")


# -- locate_gaff_dat -------------------------------------------------------


def test_locate_gaff_dat_search_order(monkeypatch, tmp_path):
    roots = {key: tmp_path / key for key in ("amberhome", "conda", "which")}
    dats = {}
    for key, root in roots.items():
        parm = root / "dat" / "leap" / "parm"
        parm.mkdir(parents=True)
        dats[key] = parm / "gaff2.dat"
        dats[key].write_text(key)
    monkeypatch.setenv("AMBERHOME", str(roots["amberhome"]))
    monkeypatch.setenv("CONDA_PREFIX", str(roots["conda"]))
    monkeypatch.setattr(shutil, "which", lambda name: str(roots["which"] / "bin" / "antechamber"))

    assert nonstandard_ffxml.locate_gaff_dat() == str(dats["amberhome"])
    dats["amberhome"].unlink()
    assert nonstandard_ffxml.locate_gaff_dat() == str(dats["conda"])
    dats["conda"].unlink()
    assert nonstandard_ffxml.locate_gaff_dat() == str(dats["which"])


def test_locate_gaff_dat_error_lists_candidates(monkeypatch, tmp_path):
    monkeypatch.delenv("AMBERHOME", raising=False)
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(FileNotFoundError, match="AMBERHOME") as excinfo:
        nonstandard_ffxml.locate_gaff_dat("gaff")
    assert str(tmp_path / "dat" / "leap" / "parm" / "gaff.dat") in str(excinfo.value)


# -- extract_residue_to_pdb ------------------------------------------------


def test_extract_residue_to_pdb_roundtrip(tmp_path):
    residue = _methanol_residue()
    positions = unit.Quantity([Vec3(*p) for p in METHANOL_XYZ], unit.angstrom)
    out = tmp_path / "LIG.pdb"
    path = nonstandard_ffxml.extract_residue_to_pdb(positions, residue, out)
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
        nonstandard_ffxml.extract_residue_to_pdb(positions, next(top.residues()), tmp_path / "unk.pdb")
    assert "no element" in caplog.text
