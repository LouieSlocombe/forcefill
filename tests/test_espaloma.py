"""Tests for the espaloma backend.

Mostly hermetic. ``espaloma`` itself is an optional dependency that pulls in
PyTorch, so the tests that would actually run a model are marked ``espaloma``
and skipped when it is absent - which is also the point of most of what is
tested here: everything reachable *without* espaloma installed must fail early
and say so, rather than after a ligand has been read and a model downloaded.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("openmmforcefields")

from forcefill import DEFAULT_ESPALOMA_FORCEFIELD, LigandSpec, build_ligand_xml
from forcefill import espaloma as espaloma_backend
from forcefill._pipeline import prepare_espaloma_backend
from forcefill._spec import ResolvedSpec

EXAMPLES = Path(__file__).parent.parent / "examples" / "data"
BENZAMIDINIUM = EXAMPLES / "benzamidinium.sdf"


def _has_espaloma() -> bool:
    try:
        import espaloma  # noqa: F401
    except ImportError:
        return False
    return True


needs_espaloma = pytest.mark.skipif(not _has_espaloma(), reason="espaloma is not installed")
needs_no_espaloma = pytest.mark.skipif(_has_espaloma(), reason="espaloma is installed")


def test_installed_models_are_reported() -> None:
    models = espaloma_backend.installed_espaloma_forcefields()
    assert models
    assert all(name.startswith("espaloma-") for name in models)
    assert DEFAULT_ESPALOMA_FORCEFIELD in models


def test_an_unknown_charge_method_is_refused_before_espaloma_is_needed() -> None:
    # Checked ahead of require_espaloma, so the message is about the typo rather
    # than about a missing optional dependency.
    with pytest.raises(ValueError, match="not one espaloma accepts"):
        espaloma_backend.espaloma_residue_ffxml(
            ResolvedSpec(name="LIG", smiles="CO", backend="espaloma", forcefield=DEFAULT_ESPALOMA_FORCEFIELD),
            "unused.xml",
            charge_method="am1bcc",  # the espaloma spelling is "am1-bcc"
        )


def test_the_charge_method_is_always_stated() -> None:
    # openmmforcefields picks a different default depending on whether
    # template_generator_kwargs was passed at all, so leaving it unset would make
    # the charges depend on that. forcefill's default has to be a real choice.
    assert espaloma_backend.DEFAULT_ESPALOMA_CHARGE_METHOD in espaloma_backend.ESPALOMA_CHARGE_METHODS
    assert espaloma_backend.DEFAULT_ESPALOMA_CHARGE_METHOD == "nn"


@needs_no_espaloma
def test_a_missing_espaloma_is_reported_by_name() -> None:
    with pytest.raises(RuntimeError, match="needs the 'espaloma' package"):
        espaloma_backend.require_espaloma()


@needs_no_espaloma
def test_the_backend_gate_fires_before_any_ligand_is_read(monkeypatch: pytest.MonkeyPatch) -> None:
    # The same contract as the gaff and smirnoff gates: a backend that cannot run
    # must say so before the first expensive step, not after.
    def fail(*args: object, **kwargs: object) -> None:
        raise AssertionError("a ligand was read before the backend was checked")

    monkeypatch.setattr("forcefill.smirnoff._load_molecule", fail)
    with pytest.raises(RuntimeError, match="needs the 'espaloma' package"):
        prepare_espaloma_backend(
            {"LIG": ResolvedSpec(name="LIG", smiles="CO", backend="espaloma", forcefield=DEFAULT_ESPALOMA_FORCEFIELD)}
        )


def test_the_gate_is_silent_when_no_ligand_uses_the_backend() -> None:
    # Nothing to check, and importing espaloma to find that out would make an
    # optional dependency mandatory for every build.
    prepare_espaloma_backend({"LIG": ResolvedSpec(name="LIG", smiles="CO", backend="smirnoff")})


def test_layering_models_is_refused_at_the_spec(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="single model"):
        build_ligand_xml(
            {"BEN": LigandSpec(file=BENZAMIDINIUM, forcefield=["espaloma-0.3.2", "extra.pt"])},
            tmp_path / "out.xml",
            backend="espaloma",
            workdir=tmp_path / "wd",
        )


def test_a_pdb_residue_is_refused_for_espaloma(tmp_path: Path) -> None:
    # Espaloma reads the chemical graph, which a PDB residue does not carry.
    with pytest.raises(ValueError, match=re.escape("espaloma backend but has no ligand source")):
        build_ligand_xml(
            {"BEN": LigandSpec()},
            tmp_path / "out.xml",
            backend="espaloma",
            workdir=tmp_path / "wd",
        )


@needs_espaloma
@pytest.mark.espaloma
def test_end_to_end_from_a_file(tmp_path: Path) -> None:
    from openmm import app

    result = build_ligand_xml(
        {"BEN": LigandSpec(file=BENZAMIDINIUM)},
        tmp_path / "out.xml",
        backend="espaloma",
        workdir=tmp_path / "wd",
        minimize=True,
    )
    assert result.parameterized == ["BEN"]
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)


@needs_espaloma
@pytest.mark.espaloma
def test_the_residue_template_is_renamed(tmp_path: Path) -> None:
    import xml.etree.ElementTree as ET

    xml = espaloma_backend.espaloma_residue_ffxml(
        ResolvedSpec(name="MOL", smiles="CO", backend="espaloma", forcefield=DEFAULT_ESPALOMA_FORCEFIELD),
        tmp_path / "MOL.xml",
    )
    names = [r.get("name") for r in ET.parse(xml).getroot().findall("./Residues/Residue")]
    assert names == ["MOL"]


@needs_espaloma
@pytest.mark.espaloma
def test_espaloma_and_smirnoff_ligands_merge_into_one_xml(tmp_path: Path) -> None:
    from openmm import app

    result = build_ligand_xml(
        {
            "BEN": LigandSpec(file=BENZAMIDINIUM, backend="espaloma"),
            "MOL": LigandSpec(smiles="CO", backend="smirnoff"),
        },
        tmp_path / "out.xml",
        workdir=tmp_path / "wd",
    )
    assert result.parameterized == ["BEN", "MOL"]
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)


def test_the_pipeline_wiring_works_without_espaloma(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Everything around ``EspalomaTemplateGenerator``, exercised without it.

    The generator is the only part that needs PyTorch, so stubbing it leaves the
    whole rest of the espaloma path under test on an ordinary machine: the
    backend dispatch, the mapped-SMILES rename, the topology the ligand is
    validated against, the merge, and the System build. The stub returns a real,
    complete ffxml - the committed AmberTools methanol fixtures - named the way
    openmmforcefields names one, so the rename has something to do.
    """
    import xml.etree.ElementTree as ET

    from openmm import app

    from forcefill import amber
    from tests.helpers import DATA

    source = amber.assemble_openmm_ffxml(
        {"[H][O][C]([H])([H])[H]": DATA / "methanol.mol2"}, [DATA / "methanol.frcmod"], tmp_path / "stub.xml"
    )
    ffxml = Path(source).read_text()
    seen: dict[str, object] = {}

    class StubGenerator:
        def __init__(self, molecules: list[object], forcefield: str, template_generator_kwargs: dict | None = None):
            seen["forcefield"] = forcefield
            seen["kwargs"] = template_generator_kwargs
            seen["n_molecules"] = len(molecules)

        def generate_residue_template(self, molecule: object) -> str:
            return ffxml

    monkeypatch.setattr(espaloma_backend, "EspalomaTemplateGenerator", StubGenerator)
    monkeypatch.setattr(espaloma_backend, "require_espaloma", lambda: None)

    result = build_ligand_xml(
        {"LIG": LigandSpec(smiles="CO")},
        tmp_path / "out.xml",
        backend="espaloma",
        espaloma_forcefield="espaloma-0.3.2",
        workdir=tmp_path / "wd",
        minimize=True,
    )

    # The model reached the generator as a bare string, not a tuple, and the
    # charge method was stated rather than left to the library's two defaults.
    assert seen["forcefield"] == "espaloma-0.3.2"
    assert seen["kwargs"] == {"charge_method": "nn"}
    assert seen["n_molecules"] == 1

    assert result.parameterized == ["LIG"]
    names = [r.get("name") for r in ET.parse(result.forcefield_xml).getroot().findall("./Residues/Residue")]
    assert names == ["LIG"], "the mapped-SMILES template name must be rewritten"
    assert result.minimizations["LIG"].n_atoms == 6
    app.ForceField("amber14-all.xml", "amber14/tip3p.xml", result.forcefield_xml)


def test_the_charge_method_reaches_the_generator(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen: dict[str, object] = {}

    class StubGenerator:
        def __init__(self, molecules: list[object], forcefield: str, template_generator_kwargs: dict | None = None):
            seen["kwargs"] = template_generator_kwargs
            raise RuntimeError("stop here: the constructor arguments are the whole point")

        def generate_residue_template(self, molecule: object) -> str:  # pragma: no cover - never reached
            raise AssertionError

    monkeypatch.setattr(espaloma_backend, "EspalomaTemplateGenerator", StubGenerator)
    monkeypatch.setattr(espaloma_backend, "require_espaloma", lambda: None)
    with pytest.raises(RuntimeError, match="stop here"):
        espaloma_backend.espaloma_residue_ffxml(
            ResolvedSpec(name="LIG", smiles="CO", backend="espaloma", forcefield=DEFAULT_ESPALOMA_FORCEFIELD),
            tmp_path / "LIG.xml",
            charge_method="gasteiger",
        )
    assert seen["kwargs"] == {"charge_method": "gasteiger"}
