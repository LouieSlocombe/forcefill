"""Tests for :mod:`forcefill.ligand_files`: reading ligand files and the preflight checks.

Hermetic - no AmberTools. Every reader is exercised twice, through RDKit and
through the text readers, because RDKit genuinely cannot parse some of what
forcefill is handed (the GAFF-typed mol2 antechamber writes, above all) and the
text path is what covers it.
"""

from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from openmm import app

from forcefill.ligand_files import (
    LigandFileInfo,
    check_geometry,
    check_matches_residue,
    inspect_ligand_file,
    residue_formula,
    residue_name_for,
    smiles_to_sdf,
    smiles_with_residue_geometry,
    split_multi_sdf,
)
from tests.helpers import write_methanol_pdb, write_methanol_sdf

DATA = Path(__file__).parent / "data"
EXAMPLES = Path(__file__).parent.parent / "examples" / "data"

#: Benzamidinium: 18 atoms, C7H9N2, formal charge +1 via an ``M  CHG`` record.
#: The example ligand, and the one whose charge people get wrong.
BENZAMIDINIUM = EXAMPLES / "benzamidinium.sdf"


#: Both reader paths. The text readers are not a second-class citizen: they are
#: what runs whenever RDKit cannot parse the file, so they are tested on every
#: format RDKit is.
BOTH_READERS = pytest.mark.parametrize("prefer_rdkit", [True, False], ids=["rdkit", "text"])


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------


@BOTH_READERS
def test_reads_methanol_sdf(tmp_path, prefer_rdkit):
    info = inspect_ligand_file(write_methanol_sdf(tmp_path / "lig.sdf"), prefer_rdkit=prefer_rdkit)
    assert info.formula == "CH4O"
    assert info.n_atoms == 6
    assert info.n_bonds == 5
    assert info.formal_charge == 0
    assert len(info.positions) == 6


@BOTH_READERS
def test_reads_the_antechamber_mol2(prefer_rdkit):
    # GAFF atom types (c3, oh, ho) are not SYBYL types and are ambiguous with
    # element symbols - ParmEd resolves them, a column parser cannot.
    info = inspect_ligand_file(DATA / "methanol.mol2", prefer_rdkit=prefer_rdkit)
    assert info.formula == "CH4O"
    assert info.n_bonds == 5


@BOTH_READERS
def test_reads_the_charged_example_ligand(prefer_rdkit):
    # The whole point of the inference: this file already says +1.
    info = inspect_ligand_file(BENZAMIDINIUM, prefer_rdkit=prefer_rdkit)
    assert info.formula == "C7H9N2"
    assert info.n_atoms == 18
    assert info.formal_charge == 1


@BOTH_READERS
def test_pdb_never_reports_a_confident_charge(tmp_path, prefer_rdkit):
    # A PDB carries no bond orders, so any formal charge read from one is a
    # guess. Reporting 0 as fact is how a charged ligand is silently
    # parameterized as neutral.
    info = inspect_ligand_file(write_methanol_pdb(tmp_path / "in.pdb"), prefer_rdkit=prefer_rdkit)
    assert info.formula == "CH4O"
    assert info.formal_charge is None


def test_reads_sdf_v3000(tmp_path):
    sdf = tmp_path / "v3000.sdf"
    sdf.write_text(
        "water\n  forcefill\n\n"
        "  0  0  0     0  0            999 V3000\n"
        "M  V30 BEGIN CTAB\n"
        "M  V30 COUNTS 3 2 0 0 0\n"
        "M  V30 BEGIN ATOM\n"
        "M  V30 1 O 0.0 0.0 0.0 0 CHG=-1\n"
        "M  V30 2 H 0.96 0.0 0.0 0\n"
        "M  V30 3 H -0.24 0.93 0.0 0\n"
        "M  V30 END ATOM\n"
        "M  V30 BEGIN BOND\n"
        "M  V30 1 1 1 2\n"
        "M  V30 2 1 1 3\n"
        "M  V30 END BOND\n"
        "M  V30 END CTAB\n"
        "M  END\n"
    )
    info = inspect_ligand_file(sdf, prefer_rdkit=False)
    assert info.formula == "H2O"
    assert info.n_bonds == 2
    assert info.formal_charge == -1


def test_v2000_legacy_charge_codes(tmp_path):
    # Pre-M-CHG files put the charge in atom-block column 37-39, where code 3 is +1.
    sdf = tmp_path / "legacy.sdf"
    sdf.write_text(
        "ion\n  forcefill\n\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
    )
    assert inspect_ligand_file(sdf, prefer_rdkit=False).formal_charge == 1


def test_unknown_suffix_is_rejected(tmp_path):
    bad = tmp_path / "lig.xyz"
    bad.write_text("1\n\nC 0 0 0\n")
    with pytest.raises(ValueError, match="unknown ligand file type"):
        inspect_ligand_file(bad)


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        inspect_ligand_file(tmp_path / "nope.sdf")


def test_malformed_sdf_is_rejected(tmp_path):
    bad = tmp_path / "bad.sdf"
    bad.write_text("title\nprog\ncomment\nnot a counts line\n")
    with pytest.raises(ValueError, match="counts line"):
        inspect_ligand_file(bad, prefer_rdkit=False)


# --------------------------------------------------------------------------
# Composition check
# --------------------------------------------------------------------------


def _methanol_residue(tmp_path):
    return next(app.PDBFile(str(write_methanol_pdb(tmp_path / "in.pdb"))).topology.residues())


def test_matching_file_passes(tmp_path):
    check_matches_residue(inspect_ligand_file(DATA / "methanol.mol2"), _methanol_residue(tmp_path), "LIG")


def test_mismatched_file_names_the_difference(tmp_path):
    with pytest.raises(ValueError) as exc:
        check_matches_residue(inspect_ligand_file(BENZAMIDINIUM), _methanol_residue(tmp_path), "LIG")
    message = str(exc.value)
    assert "C7H9N2" in message
    assert "CH4O" in message
    assert "C: file 7 vs PDB 1" in message


def test_mismatch_can_be_downgraded_to_a_warning(tmp_path, caplog):
    with caplog.at_level("WARNING"):
        check_matches_residue(inspect_ligand_file(BENZAMIDINIUM), _methanol_residue(tmp_path), "LIG", strict=False)
    assert "not the same molecule" in caplog.text


def test_residue_with_no_elements_is_reported_not_compared(tmp_path, caplog):
    top = app.Topology()
    residue = top.addResidue("LIG", top.addChain("A"))
    top.addAtom("C1", None, residue)
    with caplog.at_level("WARNING"):
        check_matches_residue(inspect_ligand_file(DATA / "methanol.mol2"), residue, "LIG")
    assert "no element assigned" in caplog.text


def test_residue_formula_counts_missing_elements():
    top = app.Topology()
    residue = top.addResidue("LIG", top.addChain("A"))
    top.addAtom("X1", None, residue)
    assert residue_formula(residue) == {"?": 1}


# --------------------------------------------------------------------------
# Geometry check
# --------------------------------------------------------------------------


def test_geometry_accepts_a_sane_molecule():
    check_geometry([(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)], "LIG")


def test_geometry_rejects_coincident_atoms():
    with pytest.raises(ValueError, match=r"0\.100 A apart"):
        check_geometry([(0.0, 0.0, 0.0), (0.1, 0.0, 0.0)], "LIG")


def test_geometry_rejects_a_missing_conformer():
    with pytest.raises(ValueError, match="no 3D conformer"):
        check_geometry([(0.0, 0.0, 0.0)] * 3, "LIG")


def test_geometry_rejects_non_finite_coordinates():
    with pytest.raises(ValueError, match="non-finite"):
        check_geometry([(0.0, 0.0, 0.0), (float("nan"), 0.0, 0.0)], "LIG")


def test_geometry_ignores_an_empty_position_list():
    check_geometry([], "LIG")


def test_geometry_can_be_downgraded_to_a_warning(caplog):
    with caplog.at_level("WARNING"):
        check_geometry([(0.0, 0.0, 0.0), (0.1, 0.0, 0.0)], "LIG", strict=False)
    assert "apart" in caplog.text


# --------------------------------------------------------------------------
# Naming and splitting
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("benzamidinium.sdf", "BEN"),
        ("lig.mol2", "LIG"),
        ("my-ligand_2.sdf", "MYL"),
        ("2ol.sdf", "L2O"),
    ],
)
def test_residue_name_from_file_name(filename, expected):
    assert residue_name_for(filename) == expected


def test_residue_name_needs_something_to_work_with():
    with pytest.raises(ValueError, match="no alphanumeric"):
        residue_name_for("---.sdf")


def _sdf_records(text, titles):
    """Concatenate *text* once per title, each terminated properly."""
    body = text.rstrip("\n").split("\n", 1)[1]
    return "".join(f"{title}\n{body}\n$$$$\n" for title in titles)


def test_split_multi_sdf_uses_the_title_line(tmp_path):
    source = tmp_path / "multi.sdf"
    source.write_text(_sdf_records(BENZAMIDINIUM.read_text(), ["AAA", "BBB"]))
    out = split_multi_sdf(source, tmp_path / "split")
    assert list(out) == ["AAA", "BBB"]
    for path in out.values():
        assert inspect_ligand_file(path, prefer_rdkit=False).formula == "C7H9N2"


def test_split_multi_sdf_falls_back_to_the_file_name(tmp_path):
    # An SDF's title line is line 1 and is very often blank; the program line
    # below it ("     RDKit          3D") must never become the residue name.
    source = tmp_path / "ligand.sdf"
    source.write_text(_sdf_records(BENZAMIDINIUM.read_text(), [""]))
    assert list(split_multi_sdf(source, tmp_path / "split")) == ["LIG"]


def test_split_multi_sdf_accepts_a_final_record_with_no_terminator(tmp_path):
    # examples/data/benzamidinium.sdf is exactly this shape.
    assert list(split_multi_sdf(BENZAMIDINIUM, tmp_path / "split")) == ["BEN"]


def test_split_multi_sdf_refuses_duplicate_names(tmp_path):
    source = tmp_path / "multi.sdf"
    source.write_text(_sdf_records(BENZAMIDINIUM.read_text(), ["SAME", "SAME"]))
    with pytest.raises(ValueError, match="already used"):
        split_multi_sdf(source, tmp_path / "split")


def test_split_multi_sdf_rejects_an_empty_file(tmp_path):
    empty = tmp_path / "empty.sdf"
    empty.write_text("\n\n")
    with pytest.raises(ValueError, match="no molecules"):
        split_multi_sdf(empty, tmp_path / "split")


# --------------------------------------------------------------------------
# SMILES conversion
# --------------------------------------------------------------------------


def test_smiles_to_sdf_embeds_hydrogens_and_charge(tmp_path):
    out = smiles_to_sdf("NC(=[NH2+])c1ccccc1", tmp_path / "ben.sdf", "BEN")
    info = inspect_ligand_file(out)
    assert info.formula == "C7H9N2"
    assert info.formal_charge == 1
    # A real conformer, not everything at the origin.
    assert len(set(info.positions)) == info.n_atoms


def test_smiles_to_sdf_is_reproducible(tmp_path):
    # The seed is fixed so a re-run gives the same conformer and so the same charges.
    first = Path(smiles_to_sdf("CCO", tmp_path / "a.sdf", "ETH")).read_text()
    second = Path(smiles_to_sdf("CCO", tmp_path / "b.sdf", "ETH")).read_text()
    assert first == second


def test_smiles_to_sdf_rejects_nonsense(tmp_path):
    with pytest.raises(RuntimeError, match="could not parse"):
        smiles_to_sdf("this is not a smiles", tmp_path / "x.sdf", "LIG")


def test_smiles_with_residue_geometry_keeps_the_structure_coordinates(tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    out = smiles_with_residue_geometry("CO", pdb, tmp_path / "lig.sdf", "LIG")
    info = inspect_ligand_file(out)
    assert info.formula == "CH4O"
    original = inspect_ligand_file(pdb, prefer_rdkit=False)
    # Flattened: pytest.approx does not compare sequences of tuples element-wise.
    flat = [c for xyz in info.positions for c in xyz]
    assert flat == pytest.approx([c for xyz in original.positions for c in xyz], abs=1e-3)


def test_smiles_with_residue_geometry_rejects_a_different_molecule(tmp_path):
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    with pytest.raises(RuntimeError, match="does not match the residue"):
        smiles_with_residue_geometry("CCO", pdb, tmp_path / "lig.sdf", "LIG")


def test_info_reports_where_it_came_from(tmp_path):
    info = inspect_ligand_file(write_methanol_sdf(tmp_path / "lig.sdf"), prefer_rdkit=False)
    assert info.source == "text"
    assert isinstance(info, LigandFileInfo)


def test_gaff_mol2_falls_back_to_the_text_reader(tmp_path):
    # Not a hypothetical: RDKit's mol2 parser expects SYBYL atom types and
    # returns None for the GAFF-typed file antechamber writes, so prefer_rdkit
    # silently lands on ParmEd. This is the reason the text readers exist now
    # that RDKit is a hard dependency.
    assert inspect_ligand_file(DATA / "methanol.mol2", prefer_rdkit=True).source == "text"
