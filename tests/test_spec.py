"""Tests for :mod:`forcefill._spec`: LigandSpec validation and defaulting.

Hermetic: pure data, no openmm and no AmberTools.
"""

import pytest

from forcefill import LigandSpec
from forcefill._spec import DEFAULT_SMIRNOFF_FORCEFIELD, _Defaults, resolve_specs


def _defaults(names=("LIG",), **kwargs):
    return _Defaults(names=frozenset(names), **kwargs)


# --------------------------------------------------------------------------
# LigandSpec construction
# --------------------------------------------------------------------------


def test_spec_defaults_are_all_inherit():
    spec = LigandSpec()
    assert spec.file is None
    assert spec.smiles is None
    assert spec.net_charge is None
    assert spec.atom_type is None
    assert spec.charge_method is None
    assert spec.backend is None
    assert spec.multiplicity == 1


def test_spec_rejects_file_and_smiles_together():
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
def test_spec_rejects_bad_values_at_construction(kwargs, match):
    # A typo is caught where it was written, not several expensive steps later.
    with pytest.raises(ValueError, match=match):
        LigandSpec(**kwargs)


def test_spec_antechamber_args_are_a_tuple():
    # A caller's list must not be able to mutate underneath a frozen spec.
    args = ["-dr", "no"]
    spec = LigandSpec(antechamber_args=args)
    args.append("-x")
    assert spec.antechamber_args == ("-dr", "no")


# --------------------------------------------------------------------------
# resolve_specs
# --------------------------------------------------------------------------


def test_resolve_fills_in_call_level_defaults():
    (resolved,) = resolve_specs({}, defaults=_defaults(atom_type="gaff", backend="gaff")).values()
    assert resolved.name == "LIG"
    assert resolved.atom_type == "gaff"
    assert resolved.charge_method == "bcc"
    assert resolved.backend == "gaff"
    assert resolved.forcefield == DEFAULT_SMIRNOFF_FORCEFIELD


def test_resolve_spec_overrides_the_default():
    specs = resolve_specs({"LIG": LigandSpec(atom_type="gaff")}, defaults=_defaults(atom_type="gaff2"))
    assert specs["LIG"].atom_type == "gaff"


def test_resolve_appends_per_ligand_antechamber_args():
    # Call-level args come first, then the ligand's - both apply.
    specs = resolve_specs(
        {"LIG": LigandSpec(antechamber_args=("-s", "2"))},
        defaults=_defaults(antechamber_args=("-dr", "no")),
    )
    assert specs["LIG"].antechamber_args == ("-dr", "no", "-s", "2")


def test_resolve_folds_in_the_legacy_mappings():
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
def test_resolve_refuses_to_pick_a_winner(legacy, spec):
    # Silently choosing one of two contradictory settings is how a ligand gets
    # parameterized as the wrong molecule.
    with pytest.raises(ValueError, match="configured twice"):
        resolve_specs({"LIG": spec}, defaults=_defaults(), **legacy)


def test_resolve_allows_a_legacy_mapping_alongside_an_unrelated_spec_field():
    specs = resolve_specs(
        {"LIG": LigandSpec(atom_type="gaff")},
        net_charges={"LIG": -1},
        defaults=_defaults(),
    )
    assert (specs["LIG"].atom_type, specs["LIG"].net_charge) == ("gaff", -1)


def test_resolve_ignores_names_outside_the_requested_set():
    # The caller reports these: only it knows whether the residue was skipped,
    # cleaned away or simply misspelled.
    specs = resolve_specs({"XXX": LigandSpec(net_charge=2)}, defaults=_defaults(names=("LIG",)))
    assert set(specs) == {"LIG"}
    assert specs["LIG"].net_charge is None


def test_resolve_requires_a_source_for_the_smirnoff_backend():
    # SMIRNOFF matches SMARTS against the chemical graph, which a PDB residue
    # cannot supply.
    with pytest.raises(ValueError, match="no ligand source"):
        resolve_specs({}, defaults=_defaults(backend="smirnoff"))


def test_resolve_accepts_smirnoff_with_a_source():
    specs = resolve_specs({"LIG": LigandSpec(smiles="CO")}, defaults=_defaults(backend="smirnoff"))
    assert specs["LIG"].backend == "smirnoff"
    assert specs["LIG"].has_explicit_bonds


def test_resolve_lets_one_ligand_opt_out_of_the_default_backend():
    specs = resolve_specs(
        {"LIG": LigandSpec(backend="gaff"), "SMI": LigandSpec(smiles="CO")},
        defaults=_defaults(names=("LIG", "SMI"), backend="smirnoff"),
    )
    assert specs["LIG"].backend == "gaff"
    assert specs["SMI"].backend == "smirnoff"
