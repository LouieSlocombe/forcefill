"""Tests for :func:`forcefill.build_ligand_xml`: parameterizing ligands with no input structure.

Hermetic: the AmberTools layer is replaced by the ``fake_ambertools`` fixture in
conftest.py, shared with test_structure.py. The real-executable version lives in
test_integration.py.
"""

from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

import forcefill
from forcefill import LigandSpec, amber, build_ligand_xml
from tests.helpers import write_methanol_pdb, write_methanol_sdf


def test_single_file_derives_its_residue_name(fake_ambertools, tmp_path):
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    result = build_ligand_xml(sdf, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd")

    assert result.parameterized == ["LIG"]
    assert Path(result.forcefield_xml).is_file()
    assert Path(result.residue_xmls["LIG"]) == tmp_path / "wd" / "LIG" / "LIG.xml"
    # No input structure, so these three have nothing to describe.
    assert result.skipped == {}
    assert result.cleaning is None
    assert result.full_minimization is None
    (ante,) = fake_ambertools["antechamber"]
    assert ante["input"] == str(sdf)


def test_mapping_names_the_residue_explicitly(fake_ambertools, tmp_path):
    sdf = write_methanol_sdf(tmp_path / "whatever.sdf")
    result = build_ligand_xml({"ABC": sdf}, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd")
    assert result.parameterized == ["ABC"]
    assert fake_ambertools["antechamber"][0]["residue"] == "ABC"


def test_sequence_of_files(fake_ambertools, tmp_path):
    first = write_methanol_sdf(tmp_path / "aaa.sdf")
    second = write_methanol_sdf(tmp_path / "bbb.sdf")
    result = build_ligand_xml([first, second], tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd")
    assert result.parameterized == ["AAA", "BBB"]


def test_colliding_derived_names_are_refused(tmp_path):
    a = write_methanol_sdf(tmp_path / "ligand1.sdf")
    b = write_methanol_sdf(tmp_path / "ligand2.sdf")
    with pytest.raises(ValueError, match="both resolve to the residue name"):
        build_ligand_xml([a, b], tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd")


def test_per_ligand_settings_reach_antechamber(fake_ambertools, tmp_path):
    # A PDB, so net_charge is the caller's to state: a PDB carries no formal
    # charges, which is exactly when an explicit one is needed.
    pdb = write_methanol_pdb(tmp_path / "lig.pdb")
    build_ligand_xml(
        {"LIG": LigandSpec(file=pdb, net_charge=-1, multiplicity=3, atom_type="gaff", antechamber_args=("-s", "2"))},
        tmp_path / "out.xml",
        base_forcefield=(),
        workdir=tmp_path / "wd",
        antechamber_args=("-dr", "no"),
    )
    (ante,) = fake_ambertools["antechamber"]
    assert ante["net_charge"] == -1
    assert ante["multiplicity"] == 3
    assert ante["atom_type"] == "gaff"
    # Call-level args first, then the ligand's.
    assert ante["extra_args"] == ("-dr", "no", "-s", "2")
    assert fake_ambertools["parmchk2"][0]["atom_type"] == "gaff"


def test_net_charge_is_read_from_the_ligand_file(fake_ambertools, tmp_path, caplog):
    # The behaviour change: a file that says +1 is no longer silently treated as
    # neutral. examples/data/benzamidinium.sdf carries an M CHG record.
    ben = Path(__file__).parent.parent / "examples" / "data" / "benzamidinium.sdf"
    with caplog.at_level("WARNING"):
        build_ligand_xml(
            {"BEN": ben}, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd", validate=False
        )
    assert fake_ambertools["antechamber"][0]["net_charge"] == 1
    assert "Using net charge +1 for BEN" in caplog.text


def test_explicit_net_charge_contradicting_the_file_is_refused(fake_ambertools, tmp_path):
    ben = Path(__file__).parent.parent / "examples" / "data" / "benzamidinium.sdf"
    with pytest.raises(ValueError, match="formal charge"):
        build_ligand_xml(
            {"BEN": LigandSpec(file=ben, net_charge=0)},
            tmp_path / "out.xml",
            base_forcefield=(),
            workdir=tmp_path / "wd",
        )
    # And it cost nothing: the check runs before the expensive step.
    assert fake_ambertools["antechamber"] == []


def test_a_ligand_with_no_source_is_refused(tmp_path):
    with pytest.raises(ValueError, match="no source"):
        build_ligand_xml({"LIG": LigandSpec()}, tmp_path / "out.xml", workdir=tmp_path / "wd")


def test_a_bare_smiles_spec_needs_a_name(tmp_path):
    with pytest.raises(ValueError, match="no residue name"):
        build_ligand_xml(LigandSpec(smiles="CO"), tmp_path / "out.xml", workdir=tmp_path / "wd")


def test_no_ligands_is_refused(tmp_path):
    with pytest.raises(ValueError, match="no ligands"):
        build_ligand_xml([], tmp_path / "out.xml", workdir=tmp_path / "wd")


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [({"atom_type": "oplsaa"}, "atom_type"), ({"charge_method": "magic"}, "charge_method"), ({"backend": "x"}, "back")],
)
def test_bad_options_rejected_early(tmp_path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        build_ligand_xml("lig.sdf", tmp_path / "out.xml", **kwargs)


def test_validation_uses_the_molecule_as_its_own_topology(fake_ambertools, tmp_path):
    # With no structure there is nothing to match the template against, so the
    # mol2 antechamber wrote supplies the topology. It must still build a System.
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    result = build_ligand_xml(
        sdf, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd", validate=True, minimize=True
    )
    report = result.minimizations["LIG"]
    assert report.n_atoms == 6
    assert report.energy_change <= 0


def test_cleanup_removes_the_workdir(fake_ambertools, tmp_path):
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    result = build_ligand_xml(
        sdf, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd", cleanup=True, validate=False
    )
    assert result.workdir is None
    assert result.residue_xmls == {}
    assert not (tmp_path / "wd").exists()
    assert Path(result.forcefield_xml).is_file()


def test_cleanup_refuses_output_inside_the_workdir(fake_ambertools, tmp_path):
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")
    with pytest.raises(ValueError, match="inside the working directory"):
        build_ligand_xml(sdf, tmp_path / "wd" / "out.xml", workdir=tmp_path / "wd", cleanup=True)


def test_failure_keeps_the_workdir(fake_ambertools, monkeypatch, tmp_path, caplog):
    sdf = write_methanol_sdf(tmp_path / "lig.sdf")

    def boom(*args, **kwargs):
        raise RuntimeError("parmchk2 exploded")

    monkeypatch.setattr(amber, "run_parmchk2", boom)
    with caplog.at_level("WARNING"), pytest.raises(RuntimeError, match="exploded"):
        build_ligand_xml(sdf, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd", cleanup=True)
    assert (tmp_path / "wd").exists()
    assert "kept for debugging" in caplog.text


def test_a_pdb_ligand_file_still_works(fake_ambertools, tmp_path):
    # A PDB is a legitimate gaff input - it just cannot supply a formal charge.
    pdb = write_methanol_pdb(tmp_path / "lig.pdb")
    result = build_ligand_xml(pdb, tmp_path / "out.xml", base_forcefield=(), workdir=tmp_path / "wd")
    assert result.parameterized == ["LIG"]
    assert fake_ambertools["antechamber"][0]["net_charge"] == 0


def test_public_api_exposes_the_entry_point():
    assert forcefill.build_ligand_xml is build_ligand_xml
    assert "build_ligand_xml" in forcefill.__all__
    assert "LigandSpec" in forcefill.__all__
