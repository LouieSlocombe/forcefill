"""Tests for :mod:`forcefill.structure`: the build_forcefield_xml orchestration.

Hermetic: the AmberTools layer is replaced by the ``fake_ambertools`` fixture in
conftest.py, which installs the committed fixtures instead of running anything.
The real-executable version lives in test_integration.py.
"""

from __future__ import annotations

import logging
import math
import shutil
from pathlib import Path
from typing import Any, NoReturn

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

import forcefill
from forcefill import _pipeline, _spec, amber, checks, clean_structure, ligand, merge, structure, topology
from tests.helpers import (
    write_broken_gly_pdb,
    write_ligand_and_glycerol_pdb,
    write_ligand_and_water_pdb,
    write_methanol_pdb,
    write_methanol_sdf,
    write_water_pdb,
)


def test_public_api_resolves() -> None:
    for name in forcefill.__all__:
        assert hasattr(forcefill, name), name
    modules = (_pipeline, _spec, amber, checks, clean_structure, ligand, merge, structure, topology)
    # The package re-exports every module's __all__ and nothing else. Ordering
    # is ruff's business (RUF022), so compare contents, not sequence.
    assert set(forcefill.__all__) == {n for m in modules for n in m.__all__}
    assert len(forcefill.__all__) == len(set(forcefill.__all__))
    # No name may be exported by two modules, or one would shadow the other.
    exported = [n for m in modules for n in m.__all__]
    assert len(exported) == len(set(exported))
    assert forcefill.build_forcefield_xml is structure.build_forcefield_xml
    assert forcefill.clean_pdb is clean_structure.clean_pdb
    assert isinstance(forcefill.__version__, str)


def test_parameterization_result_defaults() -> None:
    result = forcefill.ParameterizationResult(forcefield_xml=None)
    assert result.residue_xmls == {}
    assert result.parameterized == []
    assert result.skipped == {}
    assert result.workdir is None
    assert result.minimizations == {}
    assert result.full_minimization is None
    assert result.cleaning is None


def test_bad_options_rejected_before_the_file_is_opened(tmp_path: Path) -> None:
    # Each raises on the option, not on the absent PDB, so a typo costs nothing.
    with pytest.raises(ValueError, match="atom_type"):
        forcefill.build_forcefield_xml(tmp_path / "absent.pdb", atom_type="gaff3")
    with pytest.raises(ValueError, match="charge_method"):
        forcefill.build_forcefield_xml(tmp_path / "absent.pdb", charge_method="bbc")


# -- orchestration (AmberTools faked) --------------------------------------


def test_orchestration_end_to_end_with_fakes(fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path) -> None:
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


def test_orchestration_residue_files_bypass_extraction(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
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


def _write_charged_methanol_sdf(path: Path) -> Path:
    """Methanol's atoms with an ``M  CHG`` record bolted on.

    Deliberately synthetic: real methoxide would be one hydrogen short, which the
    composition check would then (correctly) reject. Here the atoms must still
    match the PDB residue so that the *charge inference* is the only thing under
    test.
    """
    text = Path(write_methanol_sdf(path)).read_text()
    Path(path).write_text(text.replace("M  END", "M  CHG  1   2  -1\nM  END"))
    return path


def test_orchestration_ligands_spec_replaces_the_legacy_mappings(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    forcefill.build_forcefield_xml(
        pdb,
        tmp_path / "extras.xml",
        base_forcefield=(),
        workdir=tmp_path / "wd",
        ligands={"LIG": forcefill.LigandSpec(file=sdf, multiplicity=3, atom_type="gaff", charge_method="gas")},
    )
    (ante,) = fake_ambertools["antechamber"]
    assert ante["input"] == str(sdf)
    assert ante["multiplicity"] == 3
    assert ante["atom_type"] == "gaff"
    assert ante["charge_method"] == "gas"


def test_orchestration_infers_the_net_charge_from_the_ligand_file(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # The behaviour change: a supplied file that states a formal charge is no
    # longer silently treated as neutral.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    charged = _write_charged_methanol_sdf(tmp_path / "charged.sdf")
    with caplog.at_level(logging.WARNING):
        forcefill.build_forcefield_xml(
            pdb,
            tmp_path / "extras.xml",
            base_forcefield=(),
            workdir=tmp_path / "wd",
            residue_files={"LIG": charged},
            validate=False,
        )
    assert fake_ambertools["antechamber"][0]["net_charge"] == -1
    assert "Using net charge -1 for LIG" in caplog.text


def test_orchestration_explicit_net_charge_wins_over_a_silent_file(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    # A PDB states no formal charge, so the caller's value is used unchallenged.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    forcefill.build_forcefield_xml(
        pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=tmp_path / "wd", net_charges={"LIG": -1}
    )
    assert fake_ambertools["antechamber"][0]["net_charge"] == -1


def test_orchestration_mismatched_residue_file_fails_before_antechamber(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    # The point of the preflight: this used to surface only after antechamber,
    # as an opaque OpenMM template-match error.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    ben = Path(__file__).parent.parent / "examples" / "data" / "benzamidinium.sdf"
    with pytest.raises(ValueError, match="not the same molecule"):
        forcefill.build_forcefield_xml(
            pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=tmp_path / "wd", residue_files={"LIG": ben}
        )
    assert fake_ambertools["antechamber"] == []


def test_orchestration_strict_false_downgrades_the_mismatch(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    ben = Path(__file__).parent.parent / "examples" / "data" / "benzamidinium.sdf"
    with caplog.at_level(logging.WARNING):
        forcefill.build_forcefield_xml(
            pdb,
            tmp_path / "extras.xml",
            base_forcefield=(),
            workdir=tmp_path / "wd",
            residue_files={"LIG": ben},
            validate=False,
            strict=False,
        )
    assert "not the same molecule" in caplog.text
    assert len(fake_ambertools["antechamber"]) == 1


def test_orchestration_rejects_an_unknown_backend(tmp_path: Path) -> None:
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    with pytest.raises(ValueError, match="backend"):
        forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", backend="amoeba")


def test_orchestration_smirnoff_without_a_source_is_refused(tmp_path: Path) -> None:
    # Refused before any tool runs: a PDB residue has no bond orders for SMARTS
    # matching to work on.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    with pytest.raises(ValueError, match="no ligand source"):
        forcefill.build_forcefield_xml(
            pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=tmp_path / "wd", backend="smirnoff"
        )


def test_orchestration_minimize(fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path) -> None:
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


def test_orchestration_minimize_without_validate(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
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


def test_orchestration_minimize_skips_full_structure_when_residue_skipped(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
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


def test_orchestration_cleanup_removes_workdir(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    wd = tmp_path / "wd"
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=wd, cleanup=True)
    assert result.workdir is None
    assert result.residue_xmls == {}
    assert not wd.exists()
    assert Path(result.forcefield_xml).is_file()


def test_orchestration_cleanup_refuses_output_inside_the_workdir(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    wd = tmp_path / "wd"
    with pytest.raises(ValueError, match="inside the working directory"):
        forcefill.build_forcefield_xml(pdb, wd / "extras.xml", base_forcefield=(), workdir=wd, cleanup=True)


def test_orchestration_failure_preserves_workdir(
    fake_ambertools: dict[str, list[dict[str, Any]]],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    def boom(*args: object, **kwargs: object) -> NoReturn:
        raise RuntimeError("antechamber exploded")

    monkeypatch.setattr(amber, "run_antechamber", boom)
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    wd = tmp_path / "wd"
    with caplog.at_level(logging.WARNING), pytest.raises(RuntimeError, match="exploded"):
        forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=wd, cleanup=True)
    assert wd.exists()
    assert "kept for debugging" in caplog.text


def test_orchestration_default_workdir_reported(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), validate=False)
    try:
        assert result.workdir is not None
        assert Path(result.workdir).name.startswith("nonstandard_ff_")
        assert Path(result.residue_xmls["LIG"]).is_file()
    finally:
        shutil.rmtree(result.workdir, ignore_errors=True)


def test_nothing_to_parameterize_short_circuits(tmp_path: Path) -> None:
    pdb = write_water_pdb(tmp_path / "w.pdb")
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=("amber14/tip3p.xml",))
    assert result.forcefield_xml is None
    assert result.parameterized == []
    assert not (tmp_path / "extras.xml").exists()


def test_everything_skipped_raises(tmp_path: Path) -> None:
    pdb = write_broken_gly_pdb(tmp_path / "g.pdb")
    with pytest.raises(RuntimeError, match="none can be auto-parameterized"):
        forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=())


# -- clean_structure=True --------------------------------------------------


def test_clean_structure_removes_an_additive_before_parameterizing(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    # A free-standing glycerol is indistinguishable from a ligand, so without
    # cleaning it goes to antechamber - minutes of AM1-BCC on a cryoprotectant,
    # for meaningless charges, since X-ray additives carry no hydrogens.
    pdb = write_ligand_and_glycerol_pdb(tmp_path / "in.pdb")

    # validate=False: the fake antechamber returns methanol for every residue,
    # so GOL's template cannot match its own bond graph. What matters here is
    # only that GOL reached antechamber at all.
    dirty = forcefill.build_forcefield_xml(
        pdb, tmp_path / "dirty.xml", base_forcefield=(), workdir=tmp_path / "wd1", validate=False
    )
    assert dirty.parameterized == ["GOL", "LIG"]
    assert dirty.cleaning is None

    clean = forcefill.build_forcefield_xml(
        pdb, tmp_path / "clean.xml", base_forcefield=(), workdir=tmp_path / "wd2", clean_structure=True
    )
    assert clean.parameterized == ["LIG"]
    assert clean.cleaning.removed == {"GOL": ("additive", 1)}
    assert [c["residue"] for c in fake_ambertools["antechamber"]] == ["GOL", "LIG", "LIG"]


def test_clean_structure_lets_the_full_checks_run(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path
) -> None:
    # Crystallographic water that fails to match lands in `skipped`, and a
    # non-empty `skipped` suppresses the whole-structure validate/minimize.
    pdb = write_ligand_and_water_pdb(tmp_path / "in.pdb")

    dirty = forcefill.build_forcefield_xml(
        pdb, tmp_path / "dirty.xml", base_forcefield=(), workdir=tmp_path / "wd1", minimize=True
    )
    assert "HOH" in dirty.skipped
    assert dirty.full_minimization is None

    clean = forcefill.build_forcefield_xml(
        pdb,
        tmp_path / "clean.xml",
        base_forcefield=(),
        workdir=tmp_path / "wd2",
        clean_structure=True,
        minimize=True,
    )
    assert clean.skipped == {}
    assert clean.full_minimization is not None
    # The numbers describe the cleaned system: methanol alone, not the water too.
    assert clean.full_minimization.n_atoms == 6
    assert clean.cleaning.n_atoms_after == 6


def test_clean_structure_defaults_off(fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path) -> None:
    # Deleting atoms is never something the pipeline does unasked.
    pdb = write_ligand_and_water_pdb(tmp_path / "in.pdb")
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), workdir=tmp_path / "wd")
    assert result.cleaning is None
    assert "HOH" in result.skipped


def test_clean_structure_is_reported_on_the_early_return(tmp_path: Path) -> None:
    # A solvent-only PDB cleans down to nothing, so there is no XML to build -
    # but the caller still needs to hear what happened.
    pdb = write_water_pdb(tmp_path / "w.pdb")
    result = forcefill.build_forcefield_xml(pdb, tmp_path / "extras.xml", base_forcefield=(), clean_structure=True)
    assert result.forcefield_xml is None
    assert result.cleaning is not None
    assert result.cleaning.removed == {"HOH": ("water", 1)}
    assert result.cleaning.n_atoms_after == 0


def test_clean_structure_explains_an_override_aimed_at_a_removed_residue(
    fake_ambertools: dict[str, list[dict[str, Any]]], tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    pdb = write_ligand_and_glycerol_pdb(tmp_path / "in.pdb")
    with caplog.at_level(logging.WARNING):
        forcefill.build_forcefield_xml(
            pdb,
            tmp_path / "extras.xml",
            base_forcefield=(),
            workdir=tmp_path / "wd",
            clean_structure=True,
            net_charges={"GOL": 0},
        )
    assert "was removed from the structure by clean_structure=True" in caplog.text
