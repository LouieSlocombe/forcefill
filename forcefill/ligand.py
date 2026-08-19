"""Parameterize ligands on their own, with no input structure.

:func:`~forcefill.build_forcefield_xml` starts from a PDB and asks which of its
residues the base force field cannot match. This starts from the ligands
themselves::

    build_ligand_xml("benzamidinium.sdf", "ben.xml")
    build_ligand_xml(["lig1.sdf", "lig2.sdf"], "ligs.xml")
    build_ligand_xml({"BEN": LigandSpec(smiles="NC(=[NH2+])c1ccccc1")}, "ben.xml")
    build_ligand_xml("ben.str", "ben.xml", base_forcefield=CHARMM_BASE_FORCEFIELD)

Same backends, same preflight checks, same output: an ffxml you load next to the
standard force fields. What changes is the validation. With no structure there is
no bond graph to match the generated template against, so the molecule supplies
its own topology - narrower, but not weaker: it still proves the template covers
every atom and that no parameter is missing, and with ``minimize=True`` that the
numbers are physical. The charmm backend carries no Cartesian coordinates, so it
can be validated here but not minimized.

What this cannot tell you is whether the ligand matches the complex you
eventually load it with. If you have that structure, use
:func:`~forcefill.build_forcefield_xml`, which checks exactly that.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence

from openmm import app, unit

from . import _pipeline, charmm, ligand_files, smirnoff
from ._pipeline import ParameterizationResult, ResidueArtifacts
from ._spec import (
    ATOM_TYPES,
    BACKENDS,
    CHARGE_METHODS,
    CHARMM_FILE_SUFFIXES,
    DEFAULT_BASE_FORCEFIELD,
    DEFAULT_SMIRNOFF_FORCEFIELD,
    LigandSpec,
    PathLike,
    ResolvedSpec,
    _Defaults,
    check_choice,
    resolve_specs,
)
from .amber import DEFAULT_AMBERTOOLS_TIMEOUT, load_residue_template
from .checks import MinimizationResult, minimize_with_forcefield_xml, validate_forcefield_xml
from .preflight import preflight_specs

log = logging.getLogger(__name__)

__all__ = ["build_ligand_xml"]

#: What :func:`build_ligand_xml` accepts as its first argument.
LigandInput = PathLike | LigandSpec | Sequence[PathLike | LigandSpec] | Mapping[str, LigandSpec | PathLike]


def _coerce_spec(value: LigandSpec | PathLike) -> LigandSpec:
    """Turn a bare path into a spec; pass a spec through.

    A plain string is always a *path*, never a SMILES: the two are not reliably
    distinguishable (``C`` is both a filename and methane), and guessing wrong
    would silently parameterize the wrong molecule - so a SMILES must say so, as
    ``LigandSpec(smiles=...)``. A CHARMM suffix becomes ``charmm_files``: it is a
    parameter set, not a molecule file.
    """
    if isinstance(value, LigandSpec):
        return value
    if str(value).lower().endswith(CHARMM_FILE_SUFFIXES):
        return LigandSpec(backend="charmm", charmm_files=(value,))
    return LigandSpec(file=value)


def _name_for(spec: LigandSpec, index: int) -> str:
    """Residue name for a spec given without one: from its file name."""
    source = spec.file if spec.file is not None else next(iter(spec.charmm_files), None)
    if source is not None:
        return ligand_files.residue_name_for(source)
    raise ValueError(
        f"Ligand {index} was given as a SMILES with no residue name. A file name "
        "can supply one, a SMILES cannot - pass a mapping instead, as "
        "build_ligand_xml({'LIG': LigandSpec(smiles=...)}, ...)."
    )


def _normalize_ligands(ligands: LigandInput) -> dict[str, LigandSpec]:
    """Accept a path, a spec, a sequence of either, or a name -> spec mapping."""
    if isinstance(ligands, Mapping):
        return {str(name): _coerce_spec(value) for name, value in ligands.items()}

    items = (
        [ligands] if isinstance(ligands, (LigandSpec, str, os.PathLike)) else [_coerce_spec(item) for item in ligands]
    )
    out: dict[str, LigandSpec] = {}
    for index, item in enumerate(items, start=1):
        spec = _coerce_spec(item)
        name = _name_for(spec, index)
        if name in out:
            raise ValueError(
                f"Two ligands both resolve to the residue name {name!r} "
                f"({out[name].file} and {spec.file}). Residue names must be "
                "unique - pass a mapping to name them explicitly, as "
                "build_ligand_xml({'AAA': ..., 'BBB': ...}, ...)."
            )
        out[name] = spec
    if not out:
        raise ValueError("build_ligand_xml was given no ligands.")
    return out


def _ligand_topology(
    spec: ResolvedSpec,
    artifacts: ResidueArtifacts,
    base_forcefield: Sequence[str],
) -> tuple[app.Topology, unit.Quantity | None]:
    """Topology and coordinates for validating one ligand on its own.

    For gaff that comes from the mol2 antechamber wrote - the molecule as it was
    actually parameterized - and for smirnoff from the molecule itself. A charmm
    ligand has no Cartesian coordinates, so it returns None and only the graph
    can be checked.
    """
    if spec.backend == "charmm":
        return charmm.ligand_topology(spec, base_forcefield), None
    if spec.backend == "smirnoff":
        return smirnoff.ligand_topology(spec)
    structure = load_residue_template(artifacts.mol2, spec.name).to_structure()
    topology = structure.topology
    for residue in topology.residues():
        residue.name = spec.name
    return topology, structure.positions


def build_ligand_xml(
    ligands: LigandInput,
    output_xml: PathLike = "ligand_ff.xml",
    *,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    backend: str = "gaff",
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    smirnoff_forcefield: str = DEFAULT_SMIRNOFF_FORCEFIELD,
    charmm_files: Sequence[PathLike] = (),
    workdir: PathLike | None = None,
    cleanup: bool = False,
    validate: bool = True,
    minimize: bool = False,
    strict: bool = True,
    antechamber_args: Sequence[str] = (),
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
) -> ParameterizationResult:
    """Build an OpenMM force-field XML for one or more ligands, with no input structure.

    Args:
        ligands: What to parameterize. A path (``"lig.sdf"``), a
            :class:`~forcefill.LigandSpec`, a sequence of either, or a
            ``{residue_name: spec_or_path}`` mapping. Names not given explicitly
            are derived from the file name (``benzamidinium.sdf`` -> ``BEN``). A
            bare string is always a path; a SMILES must be given as
            ``LigandSpec(smiles=...)``. A CHARMM suffix (``ben.str``) is taken as
            that ligand's CGenFF parameters and puts it on the charmm backend.
        output_xml: Where to write the combined force-field XML.
        base_forcefield: ffxml files loaded underneath the generated one for the
            validation and minimization checks. Not used to decide what needs
            parameterizing - here that is the caller's list. Must be
            :data:`~forcefill.CHARMM_BASE_FORCEFIELD` for the charmm backend.
        backend: ``"gaff"`` (default), ``"smirnoff"`` or ``"charmm"``;
            per-ligand with ``LigandSpec(backend=...)``.
        atom_type: ``"gaff2"`` (default) or ``"gaff"``. gaff backend only.
        charge_method: antechamber charge method, default ``"bcc"``. gaff only.
        smirnoff_forcefield: SMIRNOFF release for the smirnoff backend.
        charmm_files: CHARMM topology/parameter files shared by every charmm
            ligand; per-ligand stream files go in
            ``LigandSpec(charmm_files=...)`` and are appended after these.
        workdir: Directory for intermediate files. A fresh temporary directory is
            created if not given, and kept unless *cleanup*.
        cleanup: Delete the working directory after a successful build. Refuses
            to run if *output_xml* resolves inside it.
        validate: Verify that ``base_forcefield + output_xml`` builds an
            ``openmm.System`` for each ligand on its own (default).
        minimize: Also energy-minimize each ligand in vacuum, which catches
            unphysical parameters that a System build accepts. Reported in
            ``minimizations``. Unavailable for the charmm backend, whose input
            carries no coordinates.
        strict: Raise on a ligand whose geometry is unusable (coincident atoms,
            no conformer) rather than warn.
        antechamber_args: Extra raw antechamber arguments for every ligand.
        timeout: Per-invocation ceiling in seconds for antechamber / parmchk2.
            Does not apply to the smirnoff backend, which runs in-process.

    Returns:
        ParameterizationResult. ``skipped``, ``cleaning`` and
        ``full_minimization`` are always empty here - there is no input
        structure to skip residues from, clean, or minimize as a whole.
    """
    check_choice(atom_type, ATOM_TYPES, "atom_type")
    check_choice(charge_method, CHARGE_METHODS, "charge_method")
    check_choice(backend, BACKENDS, "backend")

    requested = _normalize_ligands(ligands)
    specs = resolve_specs(
        requested,
        defaults=_Defaults(
            atom_type=atom_type,
            charge_method=charge_method,
            antechamber_args=tuple(antechamber_args),
            backend=backend,
            forcefield=smirnoff_forcefield,
            charmm_files=tuple(charmm_files),
            names=frozenset(requested),
        ),
    )
    for name, spec in specs.items():
        if not spec.has_source:
            raise ValueError(
                f"Ligand {name} has no source. With no input structure to "
                "extract it from, every ligand needs a file or a SMILES: "
                "build_ligand_xml({'"
                f"{name}"
                "': 'lig.sdf'}, ...) or LigandSpec(smiles=...)."
            )
        if minimize and spec.backend == "charmm":
            raise ValueError(
                f"minimize=True cannot be used for ligand {name}: its parameters "
                "come from CHARMM files, which record internal coordinates "
                "rather than Cartesian ones, so there is no geometry to "
                "minimize. Pass minimize=False - validate=True still proves the "
                "template covers every atom and that no parameter is missing - "
                "or use build_forcefield_xml() with a structure, where the "
                "coordinates come from the PDB."
            )
    log.info("Ligands to parameterize: %s", sorted(specs))

    gaff_dat = _pipeline.prepare_gaff_backend(specs, atom_type)
    _pipeline.check_backends_match_base(specs, base_forcefield)

    minimizations: dict[str, MinimizationResult] = {}
    with _pipeline.working_directory(workdir, output_xml, prefix="ligand_ff_", cleanup=cleanup) as wd:
        # No structure to compare against, but the net charge and the geometry
        # are still read and checked before the first antechamber run.
        specs = preflight_specs(specs, {}, None, wd, strict=strict, base_forcefield=base_forcefield)

        artifacts = _pipeline.parameterize_all(
            specs, {}, None, wd, gaff_dat=gaff_dat, timeout=timeout, base_forcefield=base_forcefield
        )
        residue_xmls = {name: art.xml for name, art in artifacts.items()}

        combined = _pipeline.combine_residue_xmls(artifacts, specs, gaff_dat, output_xml, wd)
        log.info("Wrote combined force-field XML: %s", combined)

        if validate or minimize:
            # Parse the (large) base force field + new XML once; every check uses it.
            files = [*base_forcefield, combined]
            forcefield = app.ForceField(*files)
            for name in sorted(specs):
                topology, positions = _ligand_topology(specs[name], artifacts[name], base_forcefield)
                if validate:
                    validate_forcefield_xml(topology, combined, base_forcefield, forcefield=forcefield)
                if minimize:
                    minimizations[name] = minimize_with_forcefield_xml(
                        topology, positions, combined, base_forcefield, forcefield=forcefield
                    )

    return ParameterizationResult(
        forcefield_xml=combined,
        residue_xmls={} if cleanup else residue_xmls,
        parameterized=sorted(artifacts),
        workdir=None if cleanup else str(wd),
        minimizations=minimizations,
    )
