"""Per-ligand settings: what to parameterize a residue from, and how.

One :class:`LigandSpec` describes one ligand. Every field except *multiplicity*
defaults to ``None`` - inherit the call-level default - so a spec states only
what it overrides::

    build_forcefield_xml(
        "complex.pdb",
        "extras.xml",
        ligands={"BEN": LigandSpec(file="ben.sdf", net_charge=1, atom_type="gaff")},
    )

The legacy ``net_charges`` / ``multiplicities`` / ``residue_files`` mappings are
folded into specs by :func:`resolve_specs`; a residue named by both raises rather
than silently picking a winner.

Pure data and validation, importing nothing else from the package so every other
module can depend on it - which is why the shared :data:`PathLike` alias and the
cross-module settings live here too.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

#: Re-exported at package level. ``ResolvedSpec`` / ``resolve_specs`` are the
#: pipeline's internal form and stay importable from this module only.
__all__ = [
    "BACKENDS",
    "CHARMM_BASE_FORCEFIELD",
    "DEFAULT_BASE_FORCEFIELD",
    "DEFAULT_ESPALOMA_FORCEFIELD",
    "DEFAULT_SMIRNOFF_FORCEFIELD",
    "LigandSpec",
]

PathLike = str | os.PathLike

#: A SMIRNOFF force field: an installed release name, a path to an OFFXML file,
#: or a sequence of either layered left to right.
ForceFieldSelection = str | os.PathLike | Sequence[str | os.PathLike]

#: Base force field used to decide what counts as "non-standard", and loaded
#: underneath the generated XML for the validation and minimization checks.
DEFAULT_BASE_FORCEFIELD = ("amber14-all.xml", "amber14/tip3p.xml")

#: The CHARMM equivalent, for ``backend="charmm"``. Not interchangeable with
#: :data:`DEFAULT_BASE_FORCEFIELD`: Amber scales 1-4 interactions by 0.8333/0.5
#: and CHARMM by 1.0/1.0, and OpenMM will not load force fields that disagree.
#: ``charmm36.xml`` already carries every CGenFF atom type, which is what lets
#: the charmm backend emit a residue template alone - see :mod:`forcefill.charmm`.
CHARMM_BASE_FORCEFIELD = ("charmm36.xml", "charmm36/water.xml")

#: Parameterization backends. ``"gaff"`` is antechamber + parmchk2 + ParmEd;
#: ``"smirnoff"`` is openff-toolkit + openmmforcefields; ``"charmm"`` converts
#: CGenFF/ParamChem files with ParmEd; ``"espaloma"`` is the graph neural network
#: openmmforcefields wraps as ``EspalomaTemplateGenerator``.
BACKENDS = ("gaff", "smirnoff", "charmm", "espaloma")

#: File suffixes ``parmed.charmm.CharmmParameterSet`` dispatches on. It decides
#: what a file *is* from its name alone, so ``par_all36_cgenff.prm.txt`` raises
#: "Unrecognized file type" no matter what is inside it.
CHARMM_FILE_SUFFIXES = (".str", ".rtf", ".top", ".prm", ".par", ".inp")

#: SMIRNOFF release used when a spec does not name one. Pinned rather than
#: "latest": the version is part of the science, and a silent upgrade underneath
#: a saved XML would be invisible in the output. A spec may name a different
#: release, or one or more OFFXML files - see :func:`normalize_forcefield`.
DEFAULT_SMIRNOFF_FORCEFIELD = "openff-2.2.1"

#: Espaloma model used when a spec on the ``espaloma`` backend does not name one.
#: Pinned for the same reason as :data:`DEFAULT_SMIRNOFF_FORCEFIELD`. Unlike a
#: SMIRNOFF selection this is always a *single* entry - openmmforcefields'
#: ``EspalomaTemplateGenerator`` takes one model name or one ``.pt`` file, and
#: raises TypeError for anything else.
DEFAULT_ESPALOMA_FORCEFIELD = "espaloma-0.3.2"

#: Valid ``atom_type`` values: each names a GAFF parameter database ({atom_type}.dat).
ATOM_TYPES = ("gaff", "gaff2")

#: Charge methods antechamber accepts (``antechamber -L``). Version-dependent
#: (abcg2 needs AmberTools >= 23); extend it if your antechamber knows more.
CHARGE_METHODS = ("bcc", "abcg2", "gas", "mul", "cm1", "cm2", "esp", "resp", "rc", "wc", "dc")


def check_choice(value: str, valid: Sequence[str], label: str) -> None:
    """Raise ValueError early for a typo'd option instead of a late, cryptic failure."""
    if value not in valid:
        raise ValueError(f"{label}={value!r} is not one of {list(valid)}")


def normalize_forcefield(value: ForceFieldSelection, label: str) -> tuple[str, ...]:
    """Normalize a SMIRNOFF force-field selection to a non-empty tuple of strings.

    A selection is either an installed release name (``"openff-2.2.1"``) or a
    path to an OFFXML file - openmmforcefields takes both through the same
    argument, and a sequence layers them left to right, which is how a
    bespoke-torsions-only file is stacked on top of a stock release.
    """
    entries = (value,) if isinstance(value, str | os.PathLike) else tuple(value)
    if not entries:
        raise ValueError(
            f"{label} is empty. Name a SMIRNOFF release (e.g. "
            f"{DEFAULT_SMIRNOFF_FORCEFIELD!r}) or give the path to an OFFXML "
            "file; leave it unset to inherit the default."
        )
    normalized = tuple(str(entry) for entry in entries)
    for entry in normalized:
        check_not_inline_xml(entry, label)
    return normalized


def check_not_inline_xml(entry: str, label: str) -> None:
    """Refuse force-field *content* passed where a file name belongs.

    openmmforcefields accepts a force field as raw XML text, bytes or an open
    file object as well as a name or a path. forcefill does not, and the reason
    is traceability: every check it runs names the force field it checked, the
    generated XML records which file it came from, and half a megabyte of OFFXML
    inlined into an error message is unreadable. Caught here rather than inside
    :func:`forcefill.smirnoff.resolve_forcefield_files`, where it arrives as a
    filename that cannot be resolved and reads like a library bug.
    """
    if not entry.lstrip().startswith("<"):
        return
    raise ValueError(
        f"{label} was given OFFXML content rather than a file name. forcefill "
        "identifies a force field by the file it came from - the preflight "
        "checks name it, and so does the generated XML - so write the document "
        "to a '.offxml' file and pass that path instead. (openmmforcefields "
        "itself accepts inline XML; forcefill deliberately does not.)"
    )


def check_charmm_suffix(path: PathLike) -> None:
    """Raise unless *path* has a name ParmEd recognizes as a CHARMM file.

    ParmEd reads the file kind from the **name**, so a well-formed parameter file
    saved as ``par_all36_cgenff.prm.txt`` fails with a bare "Unrecognized file
    type". ``.inp`` is its exception: the kind comes from the rest of the name.
    """
    name = os.path.basename(str(path)).lower()
    if name.endswith(".inp") and ("par" in name or "top" in name):
        return
    if not name.endswith(CHARMM_FILE_SUFFIXES):
        raise ValueError(
            f"{path} is not a CHARMM file ParmEd will read: it decides the file "
            f"type from the suffix, which must be one of {list(CHARMM_FILE_SUFFIXES)}. "
            "Rename the file (a parameter set saved as '.prm.txt' is the usual "
            "case) rather than passing it as-is."
        )
    if name.endswith(".inp"):
        raise ValueError(
            f"{path} is a CHARMM '.inp' file, whose kind ParmEd reads from the "
            "rest of the name - it must contain 'par' or 'top'. Rename it, or "
            "give it a '.prm' or '.rtf' suffix instead."
        )


@dataclass(frozen=True)
class LigandSpec:
    """How to parameterize one ligand.

    Args:
        file: Ligand file with explicit bonds (SDF/MOL2, or anything antechamber
            reads), used instead of extracting the residue from the PDB. Drawn
            bond orders spare antechamber re-perceiving them from geometry, a
            classic source of silently wrong atom types. Mutually exclusive with
            *smiles*; unused by the ``charmm`` backend.
        smiles: SMILES for the ligand, converted to a 3D SDF before use. When the
            residue is also in the input PDB, its coordinates are kept and only
            the bond orders come from the SMILES.
        net_charge: Formal net charge. ``None`` reads it from *file* or *smiles*,
            falling back to 0 for a residue extracted from a PDB. Essential for
            sensible AM1-BCC charges.
        multiplicity: Spin multiplicity passed to antechamber.
        atom_type: ``"gaff2"`` or ``"gaff"``; ``None`` inherits. gaff only.
        charge_method: antechamber charge method; ``None`` inherits. gaff only.
        antechamber_args: Extra antechamber arguments for this ligand, appended
            after the call-level ones.
        backend: ``"gaff"``, ``"smirnoff"``, ``"charmm"`` or ``"espaloma"``;
            ``None`` inherits. smirnoff and espaloma need *file* or *smiles*,
            charmm needs *charmm_files*.
        forcefield: The force field for the ``smirnoff`` and ``espaloma``
            backends; ``None`` inherits. For smirnoff: an installed release name
            (``"openff-2.2.1"``), a path to an OFFXML file - a bespoke one from
            `BespokeFit <https://github.com/openforcefield/openff-bespokefit>`_,
            say - or a sequence of either, layered left to right. An OFFXML must
            be *unconstrained*; see :mod:`forcefill.smirnoff`. For espaloma: one
            model name (``"espaloma-0.3.2"``) or one ``.pt`` file, never a
            sequence - see :mod:`forcefill.espaloma`.
        charmm_files: CHARMM topology/parameter files for the ``charmm`` backend
            - typically one ``.str`` from ParamChem or the cgenff program, plus
            any extra ``.rtf``/``.prm``. Appended after the call-level ones.
    """

    file: PathLike | None = None
    smiles: str | None = None
    net_charge: int | None = None
    multiplicity: int = 1
    atom_type: str | None = None
    charge_method: str | None = None
    antechamber_args: Sequence[str] = ()
    backend: str | None = None
    forcefield: ForceFieldSelection | None = None
    charmm_files: Sequence[PathLike] = ()

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
        for path in self.charmm_files:
            check_charmm_suffix(path)
        if self.forcefield is not None:
            object.__setattr__(self, "forcefield", normalize_forcefield(self.forcefield, "LigandSpec.forcefield"))
        # Tuples, so a caller's list cannot mutate underneath a frozen dataclass.
        object.__setattr__(self, "antechamber_args", tuple(self.antechamber_args))
        object.__setattr__(self, "charmm_files", tuple(self.charmm_files))


@dataclass(frozen=True)
class ResolvedSpec:
    """A :class:`LigandSpec` with every default filled in.

    What the pipeline consumes: no ``None`` to re-interpret at each step.
    ``net_charge`` stays optional only because inferring it needs the ligand
    file, which is read later.
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
    forcefield: tuple[str, ...] = (DEFAULT_SMIRNOFF_FORCEFIELD,)
    charmm_files: tuple[PathLike, ...] = ()

    def __post_init__(self) -> None:
        """Normalize ``forcefield``, which a bare string satisfies without being one.

        ``"openff-2.2.1"`` is a perfectly good ``Sequence[str]``, so neither the
        annotation nor a type checker rejects it - and it would then be iterated
        one character at a time. Normalizing here means every construction site,
        including the ones in tests, gets the tuple the rest of the pipeline
        assumes.
        """
        object.__setattr__(self, "forcefield", normalize_forcefield(self.forcefield, "ResolvedSpec.forcefield"))

    @property
    def has_explicit_bonds(self) -> bool:
        """True when the ligand comes from a source that carries bond orders."""
        return self.file is not None or self.smiles is not None

    @property
    def has_source(self) -> bool:
        """True when the ligand carries its own chemistry, without an input structure."""
        return self.has_explicit_bonds or bool(self.charmm_files)

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
    forcefield: ForceFieldSelection = DEFAULT_SMIRNOFF_FORCEFIELD
    #: Kept apart from *forcefield* rather than sharing it: the two backends
    #: cannot read each other's selections, so one field would make
    #: ``backend="espaloma"`` silently inherit a SMIRNOFF release name.
    espaloma_forcefield: ForceFieldSelection = DEFAULT_ESPALOMA_FORCEFIELD
    charmm_files: tuple[PathLike, ...] = ()
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


def _forcefield_for(
    name: str,
    backend: str,
    spec: LigandSpec,
    default_smirnoff: tuple[str, ...],
    default_espaloma: tuple[str, ...],
) -> tuple[str, ...]:
    """The force-field selection for one resolved spec, defaulted by backend.

    The two backends that read this field cannot read each other's values, so
    the default has to follow the backend: a ligand switched to ``espaloma``
    without naming a model wants :data:`DEFAULT_ESPALOMA_FORCEFIELD`, not the
    SMIRNOFF release the call inherited.
    """
    if spec.forcefield:
        selection = normalize_forcefield(spec.forcefield, "LigandSpec.forcefield")
    else:
        selection = default_espaloma if backend == "espaloma" else default_smirnoff
    if backend == "espaloma" and len(selection) != 1:
        raise ValueError(
            f"Residue {name} uses the espaloma backend with "
            f"{list(selection)}, but an Espaloma force field is a single model: "
            "openmmforcefields' EspalomaTemplateGenerator takes one name "
            f"(e.g. {DEFAULT_ESPALOMA_FORCEFIELD!r}) or one '.pt' file and "
            "raises TypeError for anything else. Layering, which the smirnoff "
            "backend allows, has no meaning for a trained model."
        )
    return selection


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
    *defaults*. Keys outside ``defaults.names`` are ignored here - the caller
    reports them, since only it knows whether the residue was skipped, cleaned
    away or misspelled.
    """
    ligands = dict(ligands or {})
    net_charges = dict(net_charges or {})
    multiplicities = dict(multiplicities or {})
    residue_files = dict(residue_files or {})
    defaults = defaults or _Defaults()
    default_forcefield = normalize_forcefield(defaults.forcefield, "smirnoff_forcefield")
    default_espaloma = normalize_forcefield(defaults.espaloma_forcefield, "espaloma_forcefield")

    resolved: dict[str, ResolvedSpec] = {}
    for name in sorted(defaults.names):
        spec = _merge_legacy(name, ligands.get(name, LigandSpec()), net_charges, multiplicities, residue_files)
        backend = spec.backend or defaults.backend
        if backend in ("smirnoff", "espaloma") and not (spec.file or spec.smiles):
            what = "SMIRNOFF" if backend == "smirnoff" else "Espaloma"
            raise ValueError(
                f"Residue {name} uses the {backend} backend but has no ligand "
                f"source. {what} assigns parameters from the chemical graph, "
                "which a PDB residue does not carry (no bond orders), so pass "
                f"LigandSpec(file=...) or LigandSpec(smiles=...) for {name} - or "
                "use the gaff backend, which perceives bonds from geometry."
            )
        # Call-level charmm_files are shared CHARMM input and simply do not apply
        # to a ligand on another backend; a per-ligand one there is a mistake.
        charmm_files = (*defaults.charmm_files, *spec.charmm_files) if backend == "charmm" else ()
        if backend == "charmm":
            if not charmm_files:
                raise ValueError(
                    f"Residue {name} uses the charmm backend but has no CHARMM "
                    "files. forcefill converts CGenFF parameters, it does not "
                    "derive them - generate a stream file for the ligand with "
                    "ParamChem (cgenff.paramchem.org) or the cgenff program, "
                    f"then pass LigandSpec(charmm_files=['{name.lower()}.str']) "
                    "for it."
                )
            source = "smiles" if spec.smiles is not None else "file" if spec.file is not None else None
            if source is not None:
                raise ValueError(
                    f"Residue {name} uses the charmm backend and also sets "
                    f"{source}={getattr(spec, source)!r}. The charmm backend takes "
                    "the ligand's atoms, bonds and charges from its CHARMM files, "
                    f"so the {source} would be ignored - drop it, or switch this "
                    "ligand to the gaff or smirnoff backend."
                )
        elif spec.charmm_files:
            raise ValueError(
                f"Residue {name} was given charmm_files={list(spec.charmm_files)!r} "
                f"but uses the {backend} backend, which cannot read them - they "
                f"would be silently ignored. Set backend='charmm' for {name}, or "
                "drop the files."
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
            # Idempotent: LigandSpec.__post_init__ already normalized its own,
            # but doing it here too is what makes the tuple visible to a reader
            # (and a type checker) at the point ResolvedSpec is built.
            forcefield=_forcefield_for(name, backend, spec, default_forcefield, default_espaloma),
            charmm_files=charmm_files,
        )
    return resolved
