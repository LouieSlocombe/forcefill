"""Tests for :mod:`forcefill.amber`: the AmberTools wrappers and the ParmEd assembly.

Hermetic: ``_run`` is stubbed out, so no real antechamber or parmchk2 is
invoked. The real-executable version lives in test_integration.py.
"""

from __future__ import annotations

import shutil
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from forcefill import amber
from forcefill._spec import PathLike
from tests.helpers import DATA


class _RunRecorder:
    """Stands in for _run: records calls and creates the '-o' output file."""

    def __init__(self) -> None:
        self.calls: list[tuple[list[str], Path]] = []
        self.kwargs: list[dict[str, float | str | None]] = []

    def __call__(self, cmd: Sequence[str], cwd: PathLike, **kwargs: float | str | None) -> None:
        argv = [str(c) for c in cmd]
        self.calls.append((argv, Path(cwd)))
        self.kwargs.append(kwargs)
        out = Path(argv[argv.index("-o") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.touch()


@pytest.fixture
def recorder(monkeypatch: pytest.MonkeyPatch) -> _RunRecorder:
    """Stub out ``_run`` and the executable lookup; hand back the recorder."""
    rec = _RunRecorder()
    monkeypatch.setattr(amber, "require_executable", lambda name: name)
    monkeypatch.setattr(amber, "_run", rec)
    return rec


# -- running the executables -----------------------------------------------


def test_require_executable_raises_when_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="AmberTools"):
        amber.require_executable("antechamber")


def test_run_antechamber_resolves_relative_paths(
    recorder: _RunRecorder, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.chdir(tmp_path)

    amber.run_antechamber("wd/LIG/LIG.pdb", "wd/LIG/LIG.mol2", "LIG")

    ((cmd, cwd),) = recorder.calls
    in_arg = Path(cmd[cmd.index("-i") + 1])
    out_arg = Path(cmd[cmd.index("-o") + 1])
    assert in_arg == (tmp_path / "wd/LIG/LIG.pdb").resolve()
    assert out_arg == (tmp_path / "wd/LIG/LIG.mol2").resolve()
    assert cwd == (tmp_path / "wd/LIG").resolve()


def test_run_parmchk2_resolves_relative_paths(
    recorder: _RunRecorder, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.chdir(tmp_path)

    amber.run_parmchk2("wd/LIG/LIG.mol2", "wd/LIG/LIG.frcmod")

    ((cmd, cwd),) = recorder.calls
    in_arg = Path(cmd[cmd.index("-i") + 1])
    out_arg = Path(cmd[cmd.index("-o") + 1])
    assert in_arg == (tmp_path / "wd/LIG/LIG.mol2").resolve()
    assert out_arg == (tmp_path / "wd/LIG/LIG.frcmod").resolve()
    assert cwd == (tmp_path / "wd/LIG").resolve()


def test_run_antechamber_purge_scratch_flag(recorder: _RunRecorder, tmp_path: Path) -> None:
    amber.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG")
    amber.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG", purge_scratch=False)

    (cmd_default, _), (cmd_keep, _) = recorder.calls
    assert cmd_default[cmd_default.index("-pf") + 1] == "y"
    assert cmd_keep[cmd_keep.index("-pf") + 1] == "n"


def test_run_antechamber_infers_input_format(recorder: _RunRecorder, tmp_path: Path) -> None:
    amber.run_antechamber(tmp_path / "lig.sdf", tmp_path / "a.mol2", "LIG")
    amber.run_antechamber(tmp_path / "lig.mol2", tmp_path / "b.mol2", "LIG")
    amber.run_antechamber(tmp_path / "lig.xyz", tmp_path / "c.mol2", "LIG", input_format="mol2")

    formats = [cmd[cmd.index("-fi") + 1] for cmd, _ in recorder.calls]
    assert formats == ["sdf", "mol2", "mol2"]


def test_run_antechamber_unknown_suffix_raises(recorder: _RunRecorder, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="input format"):
        amber.run_antechamber(tmp_path / "lig.xyz", tmp_path / "out.mol2", "LIG")


def test_run_antechamber_missing_output_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(amber, "require_executable", lambda name: name)
    monkeypatch.setattr(amber, "_run", lambda cmd, cwd, **kw: None)  # writes nothing
    with pytest.raises(RuntimeError, match="did not write"):
        amber.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG")


def test_run_parmchk2_missing_output_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(amber, "require_executable", lambda name: name)
    monkeypatch.setattr(amber, "_run", lambda cmd, cwd, **kw: None)
    with pytest.raises(RuntimeError, match="did not write"):
        amber.run_parmchk2(tmp_path / "in.mol2", tmp_path / "out.frcmod")


def test_run_nonzero_exit_includes_tails_and_hint(tmp_path: Path) -> None:
    script = "import sys; print('OUT-MARKER'); print('ERR-MARKER', file=sys.stderr); sys.exit(3)"
    with pytest.raises(RuntimeError) as excinfo:
        amber._run([sys.executable, "-c", script], cwd=tmp_path, hint="HINT-TEXT")
    msg = str(excinfo.value)
    assert "exit code 3" in msg
    assert "OUT-MARKER" in msg
    assert "ERR-MARKER" in msg
    assert msg.endswith("HINT-TEXT")


def test_run_timeout(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="timed out"):
        amber._run([sys.executable, "-c", "import time; time.sleep(30)"], cwd=tmp_path, timeout=0.2)


def test_bad_atom_type_rejected_early(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="atom_type"):
        amber.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG", atom_type="gaff3")
    with pytest.raises(ValueError, match="atom_type"):
        amber.run_parmchk2(tmp_path / "in.mol2", tmp_path / "out.frcmod", atom_type="gaff3")
    with pytest.raises(ValueError, match="atom_type"):
        amber.locate_gaff_dat("gaff3")


def test_bad_charge_method_rejected_early(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="charge_method"):
        amber.run_antechamber(tmp_path / "in.pdb", tmp_path / "out.mol2", "LIG", charge_method="bbc")


# -- locate_gaff_dat -------------------------------------------------------


def test_locate_gaff_dat_search_order(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
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

    assert amber.locate_gaff_dat() == str(dats["amberhome"])
    dats["amberhome"].unlink()
    assert amber.locate_gaff_dat() == str(dats["conda"])
    dats["conda"].unlink()
    assert amber.locate_gaff_dat() == str(dats["which"])


def test_locate_gaff_dat_error_lists_candidates(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv("AMBERHOME", raising=False)
    monkeypatch.setenv("CONDA_PREFIX", str(tmp_path))
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(FileNotFoundError, match="AMBERHOME") as excinfo:
        amber.locate_gaff_dat("gaff")
    assert str(tmp_path / "dat" / "leap" / "parm" / "gaff.dat") in str(excinfo.value)


# -- ParmEd assembly against the committed fixtures ------------------------


def test_assemble_openmm_ffxml_semantic(tmp_path: Path) -> None:
    out = tmp_path / "sub" / "lig.xml"  # exercises output-directory creation
    path = amber.assemble_openmm_ffxml({"LIG": DATA / "methanol.mol2"}, [DATA / "methanol.frcmod"], out)
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


def test_load_residue_template_rejects_multi_residue_mol2(tmp_path: Path) -> None:
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
        amber.load_residue_template(multi, "X")


def test_load_residue_template_rejects_non_template() -> None:
    with pytest.raises(TypeError, match="ResidueTemplate"):
        amber.load_residue_template(DATA / "methanol.frcmod", "X")
