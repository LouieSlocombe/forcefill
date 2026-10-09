"""Shared backend preparation and output contracts."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest

from forcefill import _pipeline, smirnoff
from forcefill._pipeline import ResidueArtifacts
from forcefill._spec import PathLike, ResolvedSpec


def test_prepare_backends_preserves_preflight_order_and_profiles(monkeypatch: pytest.MonkeyPatch) -> None:
    events = []
    specs = {"LIG": ResolvedSpec(name="LIG")}
    profiles = {}

    def prepare_gaff(received: Mapping[str, ResolvedSpec], atom_type: str) -> str:
        assert received is specs
        assert atom_type == "gaff2"
        events.append("gaff")
        return "gaff2.dat"

    def prepare_smirnoff(received: Mapping[str, ResolvedSpec]) -> dict:
        assert received is specs
        events.append("smirnoff")
        return profiles

    def prepare_espaloma(received: Mapping[str, ResolvedSpec]) -> None:
        assert received is specs
        events.append("espaloma")

    def check_base(received: Mapping[str, ResolvedSpec], base: Sequence[str], measured: Mapping) -> None:
        assert received is specs
        assert base == ("base.xml",)
        assert measured is profiles
        events.append("base")

    monkeypatch.setattr(_pipeline, "prepare_gaff_backend", prepare_gaff)
    monkeypatch.setattr(_pipeline, "prepare_smirnoff_backend", prepare_smirnoff)
    monkeypatch.setattr(_pipeline, "prepare_espaloma_backend", prepare_espaloma)
    monkeypatch.setattr(_pipeline, "check_backends_match_base", check_base)

    gaff_dat, measured = _pipeline.prepare_backends(specs, "gaff2", ("base.xml",))
    assert gaff_dat == "gaff2.dat"
    assert measured is profiles
    assert events == ["gaff", "smirnoff", "espaloma", "base"]


@pytest.mark.parametrize(
    ("backends", "assembly_names", "merged_names"),
    [
        ((), (), None),
        ((("Z", "gaff"), ("A", "gaff")), ("A", "Z"), None),
        ((("Z", "smirnoff"),), None, ("Z.xml",)),
        ((("Z", "gaff"), ("B", "smirnoff"), ("A", "gaff")), ("A", "Z"), ("_gaff_combined.xml", "B.xml")),
    ],
)
def test_combining_assembles_gaff_once_at_the_right_destination(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    backends: tuple[tuple[str, str], ...],
    assembly_names: tuple[str, ...] | None,
    merged_names: tuple[str, ...] | None,
) -> None:
    specs = {name: ResolvedSpec(name=name, backend=backend) for name, backend in backends}
    artifacts = {
        name: ResidueArtifacts(xml=f"{name}.xml", mol2=f"{name}.mol2", frcmod=f"{name}.frcmod") for name in specs
    }
    assemblies = []
    merges = []
    output = tmp_path / "output.xml"
    workdir = tmp_path / "work"

    def assemble(templates: Mapping[str, PathLike], parameters: Sequence[PathLike], destination: PathLike) -> str:
        assemblies.append((templates, parameters, destination))
        return str(destination)

    def merge(files: Sequence[PathLike], destination: PathLike) -> str:
        merges.append((files, destination))
        return str(destination)

    monkeypatch.setattr(_pipeline.amber, "assemble_openmm_ffxml", assemble)
    monkeypatch.setattr(_pipeline, "merge_ffxml", merge)

    assert _pipeline.combine_residue_xmls(artifacts, specs, "gaff2.dat", output, workdir) == str(output)
    if assembly_names is not None:
        destination = workdir / "_gaff_combined.xml" if merged_names else output
        assert assemblies == [
            (
                {name: f"{name}.mol2" for name in assembly_names},
                ["gaff2.dat", *(f"{name}.frcmod" for name in assembly_names)],
                destination,
            )
        ]
    else:
        assert assemblies == []
    if merged_names is not None:
        expected_files = [str(workdir / name) if name == "_gaff_combined.xml" else name for name in merged_names]
        assert merges == [(expected_files, output)]
    else:
        assert merges == []


def test_generated_template_finalization_preserves_parameters_and_creates_parents(tmp_path: Path) -> None:
    ffxml = (
        '<ForceField><Residues><Residue name="mapped-smiles"><Atom name="C1" type="hash" charge="-0.123456789"/>'
        '<VirtualSite type="average2" index="1" atom1="0" atom2="0" weight1="0.5" weight2="0.5"/>'
        '</Residue></Residues><NonbondedForce coulomb14scale="0.8333333333" lj14scale="0.5"/></ForceField>'
    )
    output = tmp_path / "nested" / "ligand.xml"
    assert smirnoff._write_residue_template(ffxml, "LIG", output) == str(output)
    expected = ET.fromstring(ffxml)
    expected.find("./Residues/Residue").set("name", "LIG")
    assert ET.tostring(ET.parse(output).getroot()) == ET.tostring(expected)


@pytest.mark.parametrize("residues", ["", '<Residue name="one"/><Residue name="two"/>'])
def test_invalid_generated_templates_fail_before_creating_output(tmp_path: Path, residues: str) -> None:
    output = tmp_path / "nested" / "ligand.xml"
    with pytest.raises(RuntimeError, match="Expected exactly one residue template"):
        smirnoff._write_residue_template(f"<ForceField><Residues>{residues}</Residues></ForceField>", "LIG", output)
    assert not output.parent.exists()
