"""Tests for :func:`forcefill.merge_ffxml`.

Hermetic: hand-written ffxml fragments, no AmberTools and no OpenFF. What is
easy to lose here is *why* two sections sometimes have to stay apart - see the
ordering test.
"""

import xml.etree.ElementTree as ET

import pytest

pytest.importorskip("openmm")

from openmm import app

from forcefill import merge_ffxml

#: Amber writes the 1-4 Coulomb scale as 5/6 in full precision; openmmforcefields
#: truncates it. The same number, and merging must not care.
AMBER_SCALE = "0.8333333333333334"
SMIRNOFF_SCALE = "0.8333333333"


def _write(tmp_path, name, body):
    path = tmp_path / name
    path.write_text(f"<ForceField>\n{body}\n</ForceField>\n")
    return path


def _types(*names, element="C", mass="12.011"):
    inner = "".join(f'<Type name="{n}" class="{n}" element="{element}" mass="{mass}"/>' for n in names)
    return f"<AtomTypes>{inner}</AtomTypes>"


def _residue(name, type_name):
    return f'<Residues><Residue name="{name}"><Atom name="C1" type="{type_name}" charge="0.0"/></Residue></Residues>'


def _nonbonded(scale, class_name):
    return (
        f'<NonbondedForce coulomb14scale="{scale}" lj14scale="0.5">'
        '<UseAttributeFromResidue name="charge"/>'
        f'<Atom class="{class_name}" sigma="0.3" epsilon="0.4"/>'
        "</NonbondedForce>"
    )


def _sections(path):
    return [(child.tag, child.attrib.get("ordering")) for child in ET.parse(str(path)).getroot()]


def test_merges_two_documents(tmp_path):
    a = _write(tmp_path, "a.xml", _types("ta") + _residue("AAA", "ta"))
    b = _write(tmp_path, "b.xml", _types("tb") + _residue("BBB", "tb"))
    out = merge_ffxml([a, b], tmp_path / "merged.xml")

    root = ET.parse(out).getroot()
    assert [t.get("name") for t in root.findall("./AtomTypes/Type")] == ["ta", "tb"]
    assert [r.get("name") for r in root.findall("./Residues/Residue")] == ["AAA", "BBB"]
    assert sorted(app.ForceField(out)._templates) == ["AAA", "BBB"]


def test_folds_sections_whose_scales_agree_to_within_tolerance(tmp_path):
    # Merging these into one section is what lets a GAFF and a SMIRNOFF ligand
    # share an output file at all.
    a = _write(tmp_path, "a.xml", _types("ta") + _residue("AAA", "ta") + _nonbonded(AMBER_SCALE, "ta"))
    b = _write(tmp_path, "b.xml", _types("tb") + _residue("BBB", "tb") + _nonbonded(SMIRNOFF_SCALE, "tb"))
    out = merge_ffxml([a, b], tmp_path / "merged.xml")

    nonbonded = ET.parse(out).getroot().findall("NonbondedForce")
    assert len(nonbonded) == 1
    assert len(nonbonded[0].findall("Atom")) == 2
    # The duplicated declaration is dropped: OpenMM rejects a second copy.
    assert len(nonbonded[0].findall("UseAttributeFromResidue")) == 1
    app.ForceField(out)


def test_keeps_sections_apart_when_their_scales_really_differ(tmp_path):
    a = _write(tmp_path, "a.xml", _types("ta") + _residue("AAA", "ta") + _nonbonded("0.5", "ta"))
    b = _write(tmp_path, "b.xml", _types("tb") + _residue("BBB", "tb") + _nonbonded("1.0", "tb"))
    out = merge_ffxml([a, b], tmp_path / "merged.xml")

    assert len(ET.parse(out).getroot().findall("NonbondedForce")) == 2
    # And OpenMM is the one that gets to complain about it, as it would for two
    # separate files.
    with pytest.raises(ValueError, match="1-4 scale"):
        app.ForceField(out)


def test_torsion_sections_with_different_ordering_stay_separate(tmp_path):
    # OpenMM reads `ordering` per <Improper> at parse time, so merging a GAFF
    # section (no ordering) with a SMIRNOFF one would silently reinterpret every
    # GAFF improper under the smirnoff convention.
    gaff = _write(
        tmp_path,
        "gaff.xml",
        _types("ta") + '<PeriodicTorsionForce><Improper class1="ta" periodicity1="2" phase1="3.14" k1="4.6"/>'
        "</PeriodicTorsionForce>",
    )
    sage = _write(
        tmp_path,
        "sage.xml",
        _types("tb")
        + '<PeriodicTorsionForce ordering="smirnoff"><Improper class1="tb" periodicity1="2" phase1="3.14" k1="1.1"/>'
        "</PeriodicTorsionForce>",
    )
    out = merge_ffxml([gaff, sage], tmp_path / "merged.xml")
    assert _sections(out).count(("PeriodicTorsionForce", None)) == 1
    assert _sections(out).count(("PeriodicTorsionForce", "smirnoff")) == 1


def test_identical_atom_types_are_deduplicated(tmp_path):
    a = _write(tmp_path, "a.xml", _types("shared") + _residue("AAA", "shared"))
    b = _write(tmp_path, "b.xml", _types("shared") + _residue("BBB", "shared"))
    out = merge_ffxml([a, b], tmp_path / "merged.xml")
    assert len(ET.parse(out).getroot().findall("./AtomTypes/Type")) == 1


def test_conflicting_atom_type_is_refused(tmp_path):
    a = _write(tmp_path, "a.xml", _types("shared", element="C", mass="12.011"))
    b = _write(tmp_path, "b.xml", _types("shared", element="N", mass="14.007"))
    with pytest.raises(ValueError, match="defines type 'shared'"):
        merge_ffxml([a, b], tmp_path / "merged.xml")


def test_duplicate_residue_name_is_refused(tmp_path):
    a = _write(tmp_path, "a.xml", _types("ta") + _residue("SAME", "ta"))
    b = _write(tmp_path, "b.xml", _types("tb") + _residue("SAME", "tb"))
    with pytest.raises(ValueError, match="already defines"):
        merge_ffxml([a, b], tmp_path / "merged.xml")


def test_rejects_a_document_that_is_not_a_forcefield(tmp_path):
    bad = tmp_path / "bad.xml"
    bad.write_text("<NotAForceField/>\n")
    with pytest.raises(ValueError, match="not an OpenMM force-field XML"):
        merge_ffxml([bad], tmp_path / "merged.xml")


def test_rejects_an_empty_file_list(tmp_path):
    with pytest.raises(ValueError, match="at least one"):
        merge_ffxml([], tmp_path / "merged.xml")


def test_creates_the_output_directory(tmp_path):
    a = _write(tmp_path, "a.xml", _types("ta") + _residue("AAA", "ta"))
    out = merge_ffxml([a], tmp_path / "nested" / "deeper" / "merged.xml")
    assert (tmp_path / "nested" / "deeper" / "merged.xml").is_file()
    app.ForceField(out)
