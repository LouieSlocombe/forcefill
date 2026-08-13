"""Tests for :mod:`forcefill.nonstandard_ffxml`.

Hermetic: no AmberTools executables and no PDB fixtures are required.
openmm and parmed must be importable (they are package dependencies);
everything here is skipped cleanly when they are not.
"""

import logging
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

import forcefill
from forcefill import nonstandard_ffxml


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
        )
    assert len(caplog.records) == 3  # 'lig', 'ZN', 'XYZ'; 'LIG' is fine
    text = caplog.text
    assert "'lig'" in text and "spelling" in text
    assert "'ZN'" in text and "skipped" in text
    assert "'XYZ'" in text


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


def test_validate_parameterized_residues_ok(tmp_path):
    xml = tmp_path / "lig.xml"
    xml.write_text(
        _LIG_TEMPLATE_XML.format(
            h4_atom='<Atom name="H4" type="XH" charge="0.1"/>',
            h4_bond='<Bond atomName1="O1" atomName2="H4"/>',
        )
    )
    nonstandard_ffxml._validate_parameterized_residues({"LIG": _methanol_residue()}, xml, base_forcefield=())


def test_validate_parameterized_residues_detects_mismatch(tmp_path):
    # Template lacks the hydroxyl hydrogen -> graph mismatch -> no template
    # matches the actual residue.
    xml = tmp_path / "lig.xml"
    xml.write_text(_LIG_TEMPLATE_XML.format(h4_atom="", h4_bond=""))
    with pytest.raises(RuntimeError, match="residue LIG"):
        nonstandard_ffxml._validate_parameterized_residues({"LIG": _methanol_residue()}, xml, base_forcefield=())
