"""Tests for the SMIRNOFF backend.

These run openff-toolkit and openmmforcefields for real - there is no useful way
to fake the charge assignment - so they carry the ``smirnoff`` marker and are
deselected along with ``integration`` when you want only the fast tests.
"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from openff.toolkit import Molecule
from openff.units import unit
from openmm import app

from forcefill import LigandSpec, build_ligand_xml
from forcefill._spec import ResolvedSpec
from forcefill.smirnoff import installed_smirnoff_forcefields, smirnoff_residue_ffxml
from tests.helpers import write_methanol_pdb

EXAMPLES = Path(__file__).parent.parent / "examples" / "data"
BENZAMIDINIUM = EXAMPLES / "benzamidinium.sdf"

pytestmark = pytest.mark.smirnoff


def test_installed_forcefields_are_reported() -> None:
    installed = installed_smirnoff_forcefields()
    assert installed
    assert all(name.startswith("openff-") for name in installed)


def test_residue_template_is_renamed(tmp_path: Path) -> None:
    # openmmforcefields names the template with a mapped SMILES; a 60-character
    # SMILES in every error message is useless.
    xml = smirnoff_residue_ffxml(ResolvedSpec(name="MOL", smiles="CO"), tmp_path / "MOL.xml")
    names = [r.get("name") for r in ET.parse(xml).getroot().findall("./Residues/Residue")]
    assert names == ["MOL"]


def test_end_to_end_from_smiles(tmp_path: Path) -> None:
    result = build_ligand_xml(
        {"MOL": LigandSpec(smiles="CO")},
        tmp_path / "out.xml",
        backend="smirnoff",
        workdir=tmp_path / "wd",
        minimize=True,
    )
    assert result.parameterized == ["MOL"]
    report = result.minimizations["MOL"]
    assert report.n_atoms == 6
    assert report.energy_change <= 0
    # The product is an ordinary ffxml: it loads next to the standard set.
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)


def test_end_to_end_from_sdf_reads_the_charge(tmp_path: Path) -> None:
    result = build_ligand_xml(
        BENZAMIDINIUM, tmp_path / "out.xml", backend="smirnoff", workdir=tmp_path / "wd", minimize=True
    )
    assert result.parameterized == ["BEN"]
    assert result.minimizations["BEN"].n_atoms == 18


def test_explicit_charge_contradicting_the_molecule_is_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="formal charge"):
        build_ligand_xml(
            {"BEN": LigandSpec(file=BENZAMIDINIUM, net_charge=0)},
            tmp_path / "out.xml",
            backend="smirnoff",
            workdir=tmp_path / "wd",
        )


def test_mixed_backends_produce_one_loadable_xml(tmp_path: Path) -> None:
    # The merge is only worth having if the result actually loads: GAFF and
    # SMIRNOFF write the 1-4 scale to different precision and use different
    # improper ordering conventions.
    result = build_ligand_xml(
        {
            "BEN": LigandSpec(file=BENZAMIDINIUM, backend="smirnoff"),
            "MOL": LigandSpec(smiles="CO", backend="gaff"),
        },
        tmp_path / "mixed.xml",
        workdir=tmp_path / "wd",
        minimize=True,
    )
    assert result.parameterized == ["BEN", "MOL"]
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)
    assert {"BEN", "MOL"} <= set(forcefield._templates)

    root = ET.parse(result.forcefield_xml).getroot()
    # One NonbondedForce (the scales agree within tolerance) but two torsion
    # sections (the improper ordering conventions do not).
    assert len(root.findall("NonbondedForce")) == 1
    orderings = {section.get("ordering") for section in root.findall("PeriodicTorsionForce")}
    assert orderings == {None, "smirnoff"}


def test_smirnoff_refuses_a_ligand_with_no_bond_orders(tmp_path: Path) -> None:
    pdb = write_methanol_pdb(tmp_path / "lig.pdb")
    with pytest.raises(ValueError, match="bond orders"):
        build_ligand_xml(pdb, tmp_path / "out.xml", backend="smirnoff", workdir=tmp_path / "wd")


def test_a_named_forcefield_release_is_honoured(tmp_path: Path) -> None:
    release = installed_smirnoff_forcefields()[0]
    xml = smirnoff_residue_ffxml(ResolvedSpec(name="MOL", smiles="CO", forcefield=release), tmp_path / "MOL.xml")
    assert Path(xml).is_file()


# --------------------------------------------------------------------------
# Custom OFFXML files (BespokeFit output, and anything else shaped like it)
# --------------------------------------------------------------------------


def _bespoke_offxml(path: Path, smiles: str, *, release: str = "openff_unconstrained-2.2.1.offxml") -> Path:
    """Write an OFFXML shaped like BespokeFit's output: a release plus one bespoke torsion.

    Built from the installed release rather than committed, so the fixture never
    drifts from the toolkit the tests run against - and so half a megabyte of
    force field does not live in the repository.
    """
    from openff.toolkit import ForceField, Molecule

    forcefield = ForceField(release)
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    smirks = _bespoke_torsion_smirks(molecule)
    # Appended, so it wins: SMIRNOFF assignment is last match wins, which is how
    # BespokeFit's parameters override the general ones they sit on top of.
    forcefield.get_parameter_handler("ProperTorsions").add_parameter(
        {
            "smirks": smirks,
            "periodicity1": 1,
            "phase1": 0.0 * unit.degree,
            "k1": 1.234 * unit.kilocalorie_per_mole,
            "idivf1": 1.0,
        }
    )
    forcefield.to_file(str(path))
    return path


def _bespoke_torsion_smirks(molecule: Molecule) -> str:
    """A SMIRKS matching one torsion of *molecule* that no stock release contains.

    Every atom decorated with its element and connectivity - stock Sage torsions
    wildcard the outer two - which is the shape that makes a bespoke parameter
    specific to the molecule it was fitted for.
    """
    for bond in molecule.bonds:
        outer_first = [n for n in bond.atom1.bonded_atoms if n is not bond.atom2]
        outer_last = [n for n in bond.atom2.bonded_atoms if n is not bond.atom1]
        if not (outer_first and outer_last):
            continue
        atoms = [outer_first[0], bond.atom1, bond.atom2, outer_last[0]]
        # "~" (any bond), not "-": the stock patterns do the same, and a written
        # single bond would miss every torsion crossing a double or aromatic one.
        smirks = "~".join(
            f"[#{atom.atomic_number}X{len(list(atom.bonded_atoms))}:{i + 1}]" for i, atom in enumerate(atoms)
        )
        if molecule.chemical_environment_matches(smirks):
            return smirks
    raise AssertionError(f"no usable torsion in {molecule.to_smiles()}")


def test_a_bespoke_offxml_produces_a_loadable_ffxml(tmp_path: Path) -> None:
    offxml = _bespoke_offxml(tmp_path / "ben_bespoke.offxml", "NC(=[NH2+])c1ccccc1")
    result = build_ligand_xml(
        {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=offxml)},
        tmp_path / "out.xml",
        backend="smirnoff",
        workdir=tmp_path / "wd",
        minimize=True,
    )
    assert result.parameterized == ["BEN"]
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)


def test_a_bespoke_offxml_changes_the_torsions(tmp_path: Path) -> None:
    # The point of the exercise: the ffxml must differ from the stock one, or
    # the bespoke fit was paid for and thrown away.
    common = {"backend": "smirnoff", "workdir": tmp_path / "wd"}
    stock = build_ligand_xml({"BEN": LigandSpec(file=BENZAMIDINIUM)}, tmp_path / "stock.xml", **common)
    offxml = _bespoke_offxml(tmp_path / "ben_bespoke.offxml", "NC(=[NH2+])c1ccccc1")
    bespoke = build_ligand_xml(
        {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=offxml)}, tmp_path / "bespoke.xml", **common
    )
    assert Path(stock.forcefield_xml).read_text() != Path(bespoke.forcefield_xml).read_text()


def test_layering_a_release_and_a_file_is_accepted(tmp_path: Path) -> None:
    offxml = _bespoke_offxml(tmp_path / "ben_bespoke.offxml", "NC(=[NH2+])c1ccccc1")
    result = build_ligand_xml(
        {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=["openff_unconstrained-2.2.1.offxml", str(offxml)])},
        tmp_path / "out.xml",
        backend="smirnoff",
        workdir=tmp_path / "wd",
    )
    assert result.parameterized == ["BEN"]


def test_a_constrained_offxml_is_refused(tmp_path: Path) -> None:
    # openmmforcefields copies constraints into the residue template, so the
    # ligand's bonds would be rigid whatever createSystem was asked for.
    offxml = _bespoke_offxml(tmp_path / "constrained.offxml", "NC(=[NH2+])c1ccccc1", release="openff-2.2.1.offxml")
    with pytest.raises(ValueError, match="unconstrained"):
        build_ligand_xml(
            {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=offxml)},
            tmp_path / "out.xml",
            backend="smirnoff",
            workdir=tmp_path / "wd",
        )


def test_an_offxml_fitted_for_another_molecule_is_refused(tmp_path: Path) -> None:
    # Bespoke parameters match by SMIRKS alone, so the wrong file silently falls
    # back to stock parameters after the QC has already been paid for.
    offxml = _bespoke_offxml(tmp_path / "other.offxml", "CCCCCCO")
    with pytest.raises(ValueError, match=re.escape("nothing that openff-2.2.1 would not")):
        build_ligand_xml(
            {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=offxml)},
            tmp_path / "out.xml",
            backend="smirnoff",
            workdir=tmp_path / "wd",
        )


def test_the_wrong_molecule_is_only_a_warning_when_not_strict(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    offxml = _bespoke_offxml(tmp_path / "other.offxml", "CCCCCCO")
    with caplog.at_level(logging.WARNING, logger="forcefill.smirnoff"):
        result = build_ligand_xml(
            {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=offxml)},
            tmp_path / "out.xml",
            backend="smirnoff",
            workdir=tmp_path / "wd",
            strict=False,
        )
    assert result.parameterized == ["BEN"]
    assert "nothing that openff-2.2.1 would not" in caplog.text


def test_a_missing_offxml_fails_before_any_ligand_is_read(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # The whole point of the early gate: a typo in one path must not cost the
    # AM1-BCC charges of the ligands ahead of it.
    def fail(*args: object, **kwargs: object) -> None:
        raise AssertionError("a ligand was read before the force field was checked")

    monkeypatch.setattr("forcefill.smirnoff._load_molecule", fail)
    with pytest.raises(ValueError, match="No file exists at"):
        build_ligand_xml(
            {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=tmp_path / "typo.offxml")},
            tmp_path / "out.xml",
            backend="smirnoff",
            workdir=tmp_path / "wd",
        )
