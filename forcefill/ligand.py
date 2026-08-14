"""Parameterize ligands on their own, with no input structure.

:func:`~forcefill.build_forcefield_xml` starts from a PDB and asks which of its
residues the base force field cannot match. This starts from the ligands
themselves::

    build_ligand_xml("benzamidinium.sdf", "ben.xml")
    build_ligand_xml(["lig1.sdf", "lig2.sdf"], "ligs.xml")
    build_ligand_xml({"BEN": LigandSpec(smiles="NC(=[NH2+])c1ccccc1")}, "ben.xml")

Same backends, same preflight checks, same output: an ffxml you load next to the
standard force fields. What changes is what the validation can check. With no
structure there is no bond graph to match the generated template against, so the
molecule itself supplies the topology - which makes the check narrower but not
weaker: it still proves the template covers every atom and that no parameter is
missing, and with ``minimize=True`` that the numbers are physical.

The one thing this cannot tell you is whether the ligand matches the protein
complex you eventually load it with. If you have that structure, use
:func:`~forcefill.build_forcefield_xml` instead - it checks exactly that.
"""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

from openmm import app, unit

from . import ligand_files, smirnoff
from ._spec import (
    ATOM_TYPES,
    BACKENDS,
    CHARGE_METHODS,
    DEFAULT_SMIRNOFF_FORCEFIELD,
    LigandSpec,
    ResolvedSpec,
    _Defaults,
    resolve_specs,
)
from ._spec import check_choice as _check_choice
from .nonstandard_ffxml import (
    DEFAULT_AMBERTOOLS_TIMEOUT,
    DEFAULT_BASE_FORCEFIELD,
    MinimizationResult,
    ParameterizationResult,
    _combine_residue_xmls,
    _load_residue_template,
    _parameterize_one_residue,
    _require_executable,
    _ResidueArtifacts,
    locate_gaff_dat,
    minimize_with_forcefield_xml,
    preflight_specs,
    validate_forcefield_xml,
)

log = logging.getLogger(__name__)

__all__ = ["build_ligand_xml"]

PathLike = str | os.PathLike

#: What :func:`build_ligand_xml` accepts as its first argument.
LigandInput = PathLike | LigandSpec | Sequence[PathLike | LigandSpec] | Mapping[str, LigandSpec | PathLike]


def _coerce_spec(value: LigandSpec | PathLike) -> LigandSpec:
    """Turn a bare path into a spec; pass a spec through.

    A plain string is always a *path*, never a SMILES. The two are not reliably
    distinguishable (``C`` is both a valid filename and methane), and guessing
    wrong would silently parameterize the wrong molecule - so a SMILES must say
    so, as ``LigandSpec(smiles=...)``.
    """
    return value if isinstance(value, LigandSpec) else LigandSpec(file=value)


def _name_for(spec: LigandSpec, index: int) -> str:
    """Residue name for a spec given without one: from its file name."""
    if spec.file is not None:
        return ligand_files.residue_name_for(spec.file)
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


def _ligand_topology(spec: ResolvedSpec, artifacts: _ResidueArtifacts) -> tuple[app.Topology, unit.Quantity]:
    """Topology and coordinates for validating one ligand on its own.

    For gaff that comes from the mol2 antechamber wrote, which is the molecule as
    it was actually parameterized; for smirnoff, from the molecule itself.
    """
    if artifacts.mol2 is None:
        return smirnoff.ligand_topology(spec)
    structure = _load_residue_template(artifacts.mol2, spec.name).to_structure()
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
            ``{residue_name: spec_or_path}`` mapping. Residue names not given
            explicitly are derived from the file name
            (``benzamidinium.sdf`` -> ``BEN``). A bare string is always a file
            path; a SMILES must be given as ``LigandSpec(smiles=...)``.
        output_xml: Where to write the combined force-field XML.
        base_forcefield: ffxml files loaded underneath the generated one for the
            validation and minimization checks. Not used to decide what needs
            parameterizing - here that is the caller's list.
        backend: ``"gaff"`` (default) or ``"smirnoff"``; per-ligand with
            ``LigandSpec(backend=...)``.
        atom_type: ``"gaff2"`` (default) or ``"gaff"``. gaff backend only.
        charge_method: antechamber charge method, default ``"bcc"``. gaff only.
        smirnoff_forcefield: SMIRNOFF release for the smirnoff backend.
        workdir: Directory for intermediate files. A fresh temporary directory is
            created if not given, and kept unless *cleanup*.
        cleanup: Delete the working directory after a successful build. Refuses
            to run if *output_xml* resolves inside it.
        validate: Verify that ``base_forcefield + output_xml`` builds an
            ``openmm.System`` for each ligand on its own (default).
        minimize: Also energy-minimize each ligand in vacuum, which catches
            unphysical parameters that a System build accepts. Reported in
            ``minimizations``.
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
    _check_choice(atom_type, ATOM_TYPES, "atom_type")
    _check_choice(charge_method, CHARGE_METHODS, "charge_method")
    _check_choice(backend, BACKENDS, "backend")

    requested = _normalize_ligands(ligands)
    specs = resolve_specs(
        requested,
        defaults=_Defaults(
            atom_type=atom_type,
            charge_method=charge_method,
            antechamber_args=tuple(antechamber_args),
            backend=backend,
            forcefield=smirnoff_forcefield,
            names=frozenset(requested),
        ),
    )
    for name, spec in specs.items():
        if not spec.has_explicit_bonds:
            raise ValueError(
                f"Ligand {name} has no source. With no input structure to "
                "extract it from, every ligand needs a file or a SMILES: "
                "build_ligand_xml({'"
                f"{name}"
                "': 'lig.sdf'}, ...) or LigandSpec(smiles=...)."
            )
    log.info("Ligands to parameterize: %s", sorted(specs))

    gaff_dat: str | None = None
    if any(spec.backend == "gaff" for spec in specs.values()):
        _require_executable("antechamber")
        _require_executable("parmchk2")
        gaff_dat = locate_gaff_dat(atom_type)
        log.info("Using GAFF parameter database: %s", gaff_dat)

    workdir = Path(workdir).resolve() if workdir is not None else Path(tempfile.mkdtemp(prefix="ligand_ff_"))
    workdir.mkdir(parents=True, exist_ok=True)
    log.info("Intermediate files in %s", workdir)
    if cleanup and Path(output_xml).resolve().is_relative_to(workdir):
        raise ValueError(
            f"cleanup=True would delete the output XML: {output_xml} resolves "
            f"inside the working directory {workdir}. Write it elsewhere or "
            "pass cleanup=False."
        )

    try:
        # No structure, so there is nothing to compare against - but the net
        # charge and the geometry are still read and checked here, before the
        # first antechamber run.
        specs = preflight_specs(specs, {}, None, workdir, strict=strict)

        artifacts: dict[str, _ResidueArtifacts] = {}
        residue_xmls: dict[str, str] = {}
        minimizations: dict[str, MinimizationResult] = {}
        for name in sorted(specs):
            artifacts[name] = _parameterize_one_residue(
                specs[name], None, None, workdir / name, gaff_dat=gaff_dat, timeout=timeout
            )
            residue_xmls[name] = artifacts[name].xml

        combined = _combine_residue_xmls(artifacts, specs, gaff_dat, output_xml, workdir)
        log.info("Wrote combined force-field XML: %s", combined)

        if validate or minimize:
            # Parse the (large) base force field + new XML once; every check uses it.
            files = [*base_forcefield, combined]
            forcefield = app.ForceField(*files)
            for name in sorted(specs):
                topology, positions = _ligand_topology(specs[name], artifacts[name])
                if validate:
                    validate_forcefield_xml(topology, combined, base_forcefield, forcefield=forcefield)
                if minimize:
                    minimizations[name] = minimize_with_forcefield_xml(
                        topology, positions, combined, base_forcefield, forcefield=forcefield
                    )
    except Exception:
        # Never delete on failure: sqm.out and the intermediates are the post-mortem.
        log.warning("Intermediate files kept for debugging in %s", workdir)
        raise

    if cleanup:
        shutil.rmtree(workdir)
        log.info("Removed working directory %s", workdir)
        residue_xmls = {}

    return ParameterizationResult(
        forcefield_xml=combined,
        residue_xmls=residue_xmls,
        parameterized=sorted(artifacts),
        workdir=None if cleanup else str(workdir),
        minimizations=minimizations,
    )
