"""Per-ligand settings: what to parameterize a residue from, and how.

One :class:`LigandSpec` describes one ligand. Every knob except *multiplicity*
defaults to ``None``, meaning "inherit the call-level default", so a spec states
only what it overrides::

    build_forcefield_xml(
        "complex.pdb",
        "extras.xml",
        ligands={"BEN": LigandSpec(file="ben.sdf", net_charge=1, atom_type="gaff")},
    )

The older ``net_charges`` / ``multiplicities`` / ``residue_files`` mappings are
still accepted and are folded into specs by :func:`resolve_specs`. A residue
named by both a legacy mapping and a spec raises rather than silently picking a
winner: a disagreement there is exactly the kind of thing that produces
plausible-but-wrong charges.

This module is pure data and validation. Like :mod:`forcefill._residue_names` it
imports nothing else from the package, so every other module can depend on it.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

#: Re-exported at package level. ``ResolvedSpec`` / ``resolve_specs`` are the
#: pipeline's internal form and stay importable from this module only.
__all__ = [
    "BACKENDS",
    "DEFAULT_SMIRNOFF_FORCEFIELD",
    "LigandSpec",
]

PathLike = str | os.PathLike

#: Parameterization backends. ``"gaff"`` is antechamber + parmchk2 + ParmEd;
#: ``"smirnoff"`` is openff-toolkit + openmmforcefields.
BACKENDS = ("gaff", "smirnoff")

#: SMIRNOFF release used when a spec does not name one. Pinned rather than
#: "latest": the force field version is part of the science, and a silent
#: upgrade underneath a saved XML would be invisible in the output.
DEFAULT_SMIRNOFF_FORCEFIELD = "openff-2.2.1"

#: Valid ``atom_type`` values: each names a GAFF parameter database ({atom_type}.dat).
ATOM_TYPES = ("gaff", "gaff2")

#: Charge methods antechamber accepts (see ``antechamber -L``). The list is
#: AmberTools-version-dependent (abcg2 needs >= 23); extend it rather than
#: bypassing the check if your antechamber knows more.
CHARGE_METHODS = ("bcc", "abcg2", "gas", "mul", "cm1", "cm2", "esp", "resp", "rc", "wc", "dc")


def check_choice(value: str, valid: Sequence[str], label: str) -> None:
    """Raise ValueError early for a typo'd option instead of a late, cryptic failure."""
    if value not in valid:
        raise ValueError(f"{label}={value!r} is not one of {list(valid)}")


@dataclass(frozen=True)
class LigandSpec:
    """How to parameterize one ligand.

    Args:
        file: Ligand file with explicit bonds (SDF/MOL2, or anything antechamber
            reads) to use instead of extracting the residue from the PDB. Bond
            orders and protonation as drawn spare antechamber having to
            re-perceive them from geometry, a classic source of silently wrong
            atom types. Mutually exclusive with *smiles*.
        smiles: SMILES for the ligand, converted to a 3D SDF before use. When
            the residue also exists in the input PDB, the PDB's coordinates are
            kept and only the bond orders come from the SMILES.
        net_charge: Formal net charge. ``None`` means "read it from *file* or
            *smiles*", falling back to 0 for a residue extracted from a PDB.
            Getting this right is essential for sensible AM1-BCC charges.
        multiplicity: Spin multiplicity passed to antechamber.
        atom_type: ``"gaff2"`` or ``"gaff"``; ``None`` inherits the call-level
            default. Ignored by the ``smirnoff`` backend.
        charge_method: antechamber charge method; ``None`` inherits. Ignored by
            the ``smirnoff`` backend.
        antechamber_args: Extra raw antechamber arguments for this ligand,
            appended after the call-level ones.
        backend: ``"gaff"`` or ``"smirnoff"``; ``None`` inherits. The
            ``smirnoff`` backend needs real bond orders, so a ligand using it
            must set *file* or *smiles*.
        forcefield: SMIRNOFF release for the ``smirnoff`` backend, e.g.
            ``"openff-2.2.1"``; ``None`` inherits. Ignored by ``gaff``.
    """

    file: PathLike | None = None
    smiles: str | None = None
    net_charge: int | None = None
    multiplicity: int = 1
    atom_type: str | None = None
    charge_method: str | None = None
    antechamber_args: Sequence[str] = ()
    backend: str | None = None
    forcefield: str | None = None

    def __post_init__(self) -> None:
        """Reject contradictory or typo'd settings at construction, not mid-pipeline."""
        if self.file is not None and self.smiles is not None:
            raise ValueError(
                f"LigandSpec has both file={self.file!r} and smiles={self.smiles!r}; "
                "they are two ways to say the same thing, so pass exactly one."
            )
        if self.atom_type is not None:
            check_choice(self.atom_type, ATOM_TYPES, "atom_type")
        if self.charge_method is not None:
            check_choice(self.charge_method, CHARGE_METHODS, "charge_method")
        if self.backend is not None:
            check_choice(self.backend, BACKENDS, "backend")
        if self.multiplicity < 1:
            raise ValueError(f"multiplicity={self.multiplicity!r} must be >= 1.")
        # Tuple, so a caller's list cannot mutate underneath a frozen dataclass.
        object.__setattr__(self, "antechamber_args", tuple(self.antechamber_args))


@dataclass(frozen=True)
class ResolvedSpec:
    """A :class:`LigandSpec` with every default filled in.

    What the pipeline actually consumes: no ``None`` to re-interpret at each
    step, and ``net_charge`` still optional only because inferring it needs the
    ligand file, which is read later.
    """

    name: str
    file: PathLike | None = None
    smiles: str | None = None
    net_charge: int | None = None
    multiplicity: int = 1
    atom_type: str = "gaff2"
    charge_method: str = "bcc"
    antechamber_args: tuple[str, ...] = ()
    backend: str = "gaff"
    forcefield: str = DEFAULT_SMIRNOFF_FORCEFIELD

    @property
    def has_explicit_bonds(self) -> bool:
        """True when the ligand comes from a source that carries bond orders."""
        return self.file is not None or self.smiles is not None

    def with_net_charge(self, net_charge: int) -> ResolvedSpec:
        """Copy of this spec with *net_charge* filled in."""
        return replace(self, net_charge=net_charge)

    def with_file(self, file: PathLike) -> ResolvedSpec:
        """Copy of this spec pointing at *file* (used after a SMILES is embedded to SDF)."""
        return replace(self, file=file, smiles=None)


@dataclass
class _Defaults:
    """Call-level settings a spec inherits when it leaves a field as None."""

    atom_type: str = "gaff2"
    charge_method: str = "bcc"
    antechamber_args: tuple[str, ...] = ()
    backend: str = "gaff"
    forcefield: str = DEFAULT_SMIRNOFF_FORCEFIELD
    #: Residue names to build specs for. Everything outside this set is reported
    #: by the caller as an override that matched nothing.
    names: frozenset[str] = field(default_factory=frozenset)


def _merge_legacy(
    name: str,
    spec: LigandSpec,
    net_charges: Mapping[str, int],
    multiplicities: Mapping[str, int],
    residue_files: Mapping[str, PathLike],
) -> LigandSpec:
    """Fold the legacy per-residue mappings into *spec*, raising on any disagreement."""
    updates: dict[str, object] = {}
    for label, mapping, attr, default in (
        ("net_charges", net_charges, "net_charge", None),
        ("multiplicities", multiplicities, "multiplicity", 1),
        ("residue_files", residue_files, "file", None),
    ):
        if name not in mapping:
            continue
        legacy = mapping[name]
        current = getattr(spec, attr)
        if current != default:
            raise ValueError(
                f"Residue {name} is configured twice: {label}[{name!r}]={legacy!r} "
                f"and LigandSpec.{attr}={current!r}. Set it in one place - the "
                "LigandSpec is the one that can express everything."
            )
        updates[attr] = legacy
    if "file" in updates and spec.smiles is not None:
        raise ValueError(
            f"Residue {name} is configured twice: residue_files[{name!r}]="
            f"{updates['file']!r} and LigandSpec.smiles={spec.smiles!r}. "
            "Pass exactly one source for the ligand."
        )
    return replace(spec, **updates) if updates else spec


def resolve_specs(
    ligands: Mapping[str, LigandSpec] | None,
    *,
    net_charges: Mapping[str, int] | None = None,
    multiplicities: Mapping[str, int] | None = None,
    residue_files: Mapping[str, PathLike] | None = None,
    defaults: _Defaults | None = None,
) -> dict[str, ResolvedSpec]:
    """Build one :class:`ResolvedSpec` per residue in ``defaults.names``.

    Merges *ligands* with the legacy mappings and fills every unset field from
    *defaults*. Keys naming a residue outside ``defaults.names`` are ignored
    here - the caller reports them, because only it knows whether the residue was
    skipped, cleaned away or simply misspelled.
    """
    ligands = dict(ligands or {})
    net_charges = dict(net_charges or {})
    multiplicities = dict(multiplicities or {})
    residue_files = dict(residue_files or {})
    defaults = defaults or _Defaults()

    resolved: dict[str, ResolvedSpec] = {}
    for name in sorted(defaults.names):
        spec = _merge_legacy(name, ligands.get(name, LigandSpec()), net_charges, multiplicities, residue_files)
        backend = spec.backend or defaults.backend
        if backend == "smirnoff" and not (spec.file or spec.smiles):
            raise ValueError(
                f"Residue {name} uses the smirnoff backend but has no ligand "
                "source. SMIRNOFF assigns parameters from the chemical graph, "
                "which a PDB residue does not carry (no bond orders), so pass "
                f"LigandSpec(file=...) or LigandSpec(smiles=...) for {name} - or "
                "use the gaff backend, which perceives bonds from geometry."
            )
        resolved[name] = ResolvedSpec(
            name=name,
            file=spec.file,
            smiles=spec.smiles,
            net_charge=spec.net_charge,
            multiplicity=spec.multiplicity,
            atom_type=spec.atom_type or defaults.atom_type,
            charge_method=spec.charge_method or defaults.charge_method,
            antechamber_args=(*defaults.antechamber_args, *spec.antechamber_args),
            backend=backend,
            forcefield=spec.forcefield or defaults.forcefield,
        )
    return resolved
