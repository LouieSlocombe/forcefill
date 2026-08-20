"""Tests for the SMIRNOFF backend.

These run openff-toolkit and openmmforcefields for real - there is no useful way
to fake the charge assignment - so they carry the ``smirnoff`` marker and are
deselected along with ``integration`` when you want only the fast tests.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

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
