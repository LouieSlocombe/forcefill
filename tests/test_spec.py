"""Tests for :mod:`forcefill._spec`: LigandSpec validation and defaulting.

Hermetic: pure data, no openmm and no AmberTools.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import pytest

from forcefill import LigandSpec
from forcefill._spec import DEFAULT_SMIRNOFF_FORCEFIELD, _Defaults, resolve_specs


def _defaults(names: Iterable[str] = ("LIG",), **kwargs: str | tuple[str, ...]) -> _Defaults:
    return _Defaults(names=frozenset(names), **kwargs)


# --------------------------------------------------------------------------
# LigandSpec construction
# --------------------------------------------------------------------------


def test_spec_defaults_are_all_inherit() -> None:
    spec = LigandSpec()
    assert spec.file is None
    assert spec.smiles is None
    assert spec.net_charge is None
    assert spec.atom_type is None
    assert spec.charge_method is None
    assert spec.backend is None
    assert spec.multiplicity == 1


def test_spec_rejects_file_and_smiles_together() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        LigandSpec(file="lig.sdf", smiles="CO")


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"atom_type": "oplsaa"}, "atom_type"),
        ({"charge_method": "magic"}, "charge_method"),
        ({"backend": "amoeba"}, "backend"),
        ({"multiplicity": 0}, "multiplicity"),
    ],
)
def test_spec_rejects_bad_values_at_construction(kwargs: Mapping[str, str | int], match: str) -> None:
    # A typo is caught where it was written, not several expensive steps later.
    with pytest.raises(ValueError, match=match):
        LigandSpec(**kwargs)


def test_spec_antechamber_args_are_a_tuple() -> None:
    # A caller's list must not be able to mutate underneath a frozen spec.
    args = ["-dr", "no"]
    spec = LigandSpec(antechamber_args=args)
    args.append("-x")
    assert spec.antechamber_args == ("-dr", "no")


def test_spec_charmm_files_are_a_tuple() -> None:
    files = ["lig.str"]
    spec = LigandSpec(charmm_files=files)
    files.append("more.prm")
    assert spec.charmm_files == ("lig.str",)


def test_spec_rejects_a_charmm_file_parmed_cannot_type() -> None:
    # ParmEd decides what a CHARMM file is from its suffix, so this one would
    # fail with a bare "Unrecognized file type" much later.
    with pytest.raises(ValueError, match="decides the file type from the suffix"):
        LigandSpec(charmm_files=["par_all36_cgenff.prm.txt"])


def test_spec_accepts_a_charmm_inp_whose_name_says_what_it_is() -> None:
    # ParmEd reads the kind of a .inp from the rest of the name, not the suffix.
    assert LigandSpec(charmm_files=["toppar/par_all36_cgenff.inp"]).charmm_files


def test_spec_rejects_a_charmm_inp_whose_name_does_not() -> None:
    with pytest.raises(ValueError, match="must contain 'par' or 'top'"):
        LigandSpec(charmm_files=["ligand.inp"])


# --------------------------------------------------------------------------
# resolve_specs
# --------------------------------------------------------------------------


def test_resolve_fills_in_call_level_defaults() -> None:
    (resolved,) = resolve_specs({}, defaults=_defaults(atom_type="gaff", backend="gaff")).values()
    assert resolved.name == "LIG"
    assert resolved.atom_type == "gaff"
    assert resolved.charge_method == "bcc"
    assert resolved.backend == "gaff"
    assert resolved.forcefield == DEFAULT_SMIRNOFF_FORCEFIELD


def test_resolve_spec_overrides_the_default() -> None:
    specs = resolve_specs({"LIG": LigandSpec(atom_type="gaff")}, defaults=_defaults(atom_type="gaff2"))
    assert specs["LIG"].atom_type == "gaff"


def test_resolve_appends_per_ligand_antechamber_args() -> None:
    # Call-level args come first, then the ligand's - both apply.
    specs = resolve_specs(
        {"LIG": LigandSpec(antechamber_args=("-s", "2"))},
        defaults=_defaults(antechamber_args=("-dr", "no")),
    )
    assert specs["LIG"].antechamber_args == ("-dr", "no", "-s", "2")


def test_resolve_folds_in_the_legacy_mappings() -> None:
    specs = resolve_specs(
        None,
        net_charges={"LIG": -1},
        multiplicities={"LIG": 3},
        residue_files={"LIG": "lig.sdf"},
        defaults=_defaults(),
    )
    assert (specs["LIG"].net_charge, specs["LIG"].multiplicity, specs["LIG"].file) == (-1, 3, "lig.sdf")


@pytest.mark.parametrize(
    ("legacy", "spec"),
    [
        ({"net_charges": {"LIG": -1}}, LigandSpec(net_charge=1)),
        ({"multiplicities": {"LIG": 3}}, LigandSpec(multiplicity=2)),
        ({"residue_files": {"LIG": "a.sdf"}}, LigandSpec(file="b.sdf")),
        ({"residue_files": {"LIG": "a.sdf"}}, LigandSpec(smiles="CO")),
    ],
)
def test_resolve_refuses_to_pick_a_winner(legacy: Mapping[str, Mapping[str, int | str]], spec: LigandSpec) -> None:
    # Silently choosing one of two contradictory settings is how a ligand gets
    # parameterized as the wrong molecule.
    with pytest.raises(ValueError, match="configured twice"):
        resolve_specs({"LIG": spec}, defaults=_defaults(), **legacy)


def test_resolve_allows_a_legacy_mapping_alongside_an_unrelated_spec_field() -> None:
    specs = resolve_specs(
        {"LIG": LigandSpec(atom_type="gaff")},
        net_charges={"LIG": -1},
        defaults=_defaults(),
    )
    assert (specs["LIG"].atom_type, specs["LIG"].net_charge) == ("gaff", -1)


def test_resolve_ignores_names_outside_the_requested_set() -> None:
    # The caller reports these: only it knows whether the residue was skipped,
    # cleaned away or simply misspelled.
    specs = resolve_specs({"XXX": LigandSpec(net_charge=2)}, defaults=_defaults(names=("LIG",)))
    assert set(specs) == {"LIG"}
    assert specs["LIG"].net_charge is None


def test_resolve_requires_a_source_for_the_smirnoff_backend() -> None:
    # SMIRNOFF matches SMARTS against the chemical graph, which a PDB residue
    # cannot supply.
    with pytest.raises(ValueError, match="no ligand source"):
        resolve_specs({}, defaults=_defaults(backend="smirnoff"))


def test_resolve_accepts_smirnoff_with_a_source() -> None:
    specs = resolve_specs({"LIG": LigandSpec(smiles="CO")}, defaults=_defaults(backend="smirnoff"))
    assert specs["LIG"].backend == "smirnoff"
    assert specs["LIG"].has_explicit_bonds


def test_resolve_lets_one_ligand_opt_out_of_the_default_backend() -> None:
    specs = resolve_specs(
        {"LIG": LigandSpec(backend="gaff"), "SMI": LigandSpec(smiles="CO")},
        defaults=_defaults(names=("LIG", "SMI"), backend="smirnoff"),
    )
    assert specs["LIG"].backend == "gaff"
    assert specs["SMI"].backend == "smirnoff"


# --------------------------------------------------------------------------
# The charmm backend
# --------------------------------------------------------------------------


def test_resolve_requires_charmm_files_for_the_charmm_backend() -> None:
    # forcefill converts CGenFF parameters; nothing in it can derive them.
    with pytest.raises(ValueError, match="no CHARMM files"):
        resolve_specs({}, defaults=_defaults(backend="charmm"))


def test_resolve_appends_per_ligand_charmm_files() -> None:
    # Shared toppar first, then the ligand's own stream file.
    specs = resolve_specs(
        {"LIG": LigandSpec(charmm_files=("lig.str",))},
        defaults=_defaults(backend="charmm", charmm_files=("par_all36_cgenff.prm",)),
    )
    assert specs["LIG"].charmm_files == ("par_all36_cgenff.prm", "lig.str")
    assert specs["LIG"].has_source


@pytest.mark.parametrize("source", [{"file": "lig.sdf"}, {"smiles": "CO"}])
def test_resolve_refuses_a_second_source_for_a_charmm_ligand(source: Mapping[str, str]) -> None:
    # The RESI block names its own atoms and charges, so a molecule file would
    # be read by nothing - silently, if it were allowed.
    with pytest.raises(ValueError, match="charmm backend and also sets"):
        resolve_specs(
            {"LIG": LigandSpec(charmm_files=("lig.str",), **source)},
            defaults=_defaults(backend="charmm"),
        )


def test_resolve_refuses_charmm_files_on_another_backend() -> None:
    with pytest.raises(ValueError, match="cannot read them"):
        resolve_specs({"LIG": LigandSpec(charmm_files=("lig.str",))}, defaults=_defaults(backend="gaff"))


def test_resolve_leaves_shared_charmm_files_off_another_backend() -> None:
    # A call-level toppar file is shared CHARMM input, not a per-ligand setting:
    # it simply does not apply to a gaff ligand, and does not make it look like
    # one that carries its own chemistry.
    specs = resolve_specs(
        {"LIG": LigandSpec(backend="gaff")},
        defaults=_defaults(backend="charmm", charmm_files=("par_all36_cgenff.prm",)),
    )
    assert specs["LIG"].charmm_files == ()
    assert not specs["LIG"].has_source
