"""Identify the non-standard residues in a structure and build a force field for them.

The entry point that starts from a PDB; :func:`~forcefill.build_ligand_xml`
starts from the ligands themselves. The pipeline between the two is shared and
lives in :mod:`forcefill._pipeline`.

Pipeline:
    1. :func:`~forcefill.find_nonstandard_residues` asks the base force field
       which residues it cannot match.
    2. Those are classified, and only chemistry a stand-alone GAFF treatment is
       *valid* for is parameterized:

       * standard residues merely missing atoms -> skipped (repair the structure
         with PDBFixer / ``Modeller.addHydrogens`` instead)
       * monatomic species (ions) -> skipped (load an ion parameter file)
       * residues covalently bonded to neighbours -> skipped (a polymer-linked
         residue needs capping and a consistent charge derivation, pyRED-style)
       * free-standing hetero molecules (ligands, cofactors) -> parameterized

    3. Everything checkable is checked before anything expensive runs; see
       :mod:`forcefill.preflight`.
    4. Each unique residue goes through its backend: :mod:`forcefill.amber`,
       :mod:`forcefill.smirnoff` or :mod:`forcefill.charmm`.
    5. The per-residue XMLs are combined into one.
    6. That XML is validated, and optionally minimized, by
       :mod:`forcefill.checks`.

Requirements:
    ``openmm``, ``parmed``, ``rdkit``, ``openff-toolkit`` and
    ``openmmforcefields`` at import time, plus AmberTools (``antechamber``,
    ``parmchk2``) on ``PATH`` for the gaff backend. conda-forge is the
    recommended route for the whole stack; see ``environment.yml``.

Example:
    >>> from forcefill import build_forcefield_xml
    >>> result = build_forcefield_xml("complex.pdb", "extras.xml", net_charges={"LIG": -1})
    >>> result.parameterized
    ['LIG']

    then simulate with::

        ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
        system = ff.createSystem(pdb.topology, ...)

Notes:
    * A crystallization additive such as glycerol is a free-standing hetero
      molecule too, so step 2 parameterizes it. Strip those first with
      :func:`forcefill.clean_pdb`, or pass ``clean_structure=True`` to do it in
      memory.
    * Ligands must carry **all explicit hydrogens** with reasonable geometry;
      AM1-BCC charges are meaningless otherwise.
    * Load either the combined XML *or* the per-residue XMLs, never both -
      duplicate GAFF atom-type definitions would collide.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence

from openmm import app

from . import _pipeline
from ._pipeline import ParameterizationResult
from ._spec import (
    ATOM_TYPES,
    BACKENDS,
    CHARGE_METHODS,
    DEFAULT_BASE_FORCEFIELD,
    DEFAULT_SMIRNOFF_FORCEFIELD,
    LigandSpec,
    PathLike,
    _Defaults,
    check_choice,
    resolve_specs,
)
from .amber import DEFAULT_AMBERTOOLS_TIMEOUT
from .checks import (
    MinimizationResult,
    _minimize_parameterized_residues,
    _validate_parameterized_residues,
    minimize_with_forcefield_xml,
    validate_forcefield_xml,
)
from .clean_structure import CleaningResult, clean_topology
from .preflight import preflight_specs
from .topology import _classify_unmatched, _warn_unused_overrides, find_nonstandard_residues

log = logging.getLogger(__name__)

__all__ = ["build_forcefield_xml"]


def build_forcefield_xml(
    pdb_file: PathLike,
    output_xml: PathLike = "nonstandard_ff.xml",
    *,
    clean_structure: bool = False,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    ligands: Mapping[str, LigandSpec] | None = None,
    net_charges: Mapping[str, int] | None = None,
    multiplicities: Mapping[str, int] | None = None,
    residue_files: Mapping[str, PathLike] | None = None,
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
    """Identify non-standard residues in *pdb_file* and build an OpenMM force-field XML for them.

    Args:
        pdb_file: Input structure. Ligands must have explicit hydrogens; the
            element columns and (for hetero groups) CONECT records should be
            present.
        output_xml: Where to write the combined force-field XML.
        clean_structure: Remove crystallographic water, bulk counter-ions and
            crystallization additives *in memory* before anything else looks at
            the structure, via :func:`~forcefill.clean_topology` with its
            defaults - so structural metals are kept. Off by default because it
            deletes atoms; what went is reported in ``cleaning``. Everything
            downstream then describes the *cleaned* system, including
            ``full_minimization.n_atoms``. For anything beyond the defaults,
            call :func:`~forcefill.clean_pdb` yourself and pass the result in.
        base_forcefield: ffxml files defining what counts as "standard".
        ligands: ``{residue_name: LigandSpec}`` - per-ligand settings, each
            inheriting the call-level defaults below for whatever it does not
            set. The general form of *net_charges* / *multiplicities* /
            *residue_files*, which still work and are folded in; setting the
            same thing both ways raises.
        net_charges: ``{residue_name: net_charge}``. A ligand supplied as a file
            or SMILES states its own formal charge, so this is only needed for a
            residue extracted from the PDB - or to override the file, which
            raises if the two disagree.
        multiplicities: ``{residue_name: spin_multiplicity}``; defaults to 1.
        residue_files: ``{residue_name: ligand_file}`` - parameterize these
            residues from the given SDF/MOL2 (or other antechamber-readable)
            file instead of extracting them from the PDB, so antechamber need
            not re-perceive bond orders from geometry. The file must hold the
            *same* atoms and bonds (hydrogens included) as the PDB residue;
            *strict* checks that up front.
        backend: ``"gaff"`` (default), ``"smirnoff"`` or ``"charmm"``. SMIRNOFF
            assigns parameters from the chemical graph, so its ligands need a
            ``file`` or ``smiles``; charmm needs ``charmm_files``. Set it per
            ligand with ``LigandSpec(backend=...)`` to mix gaff and smirnoff.
            CHARMM mixes with neither - it scales 1-4 interactions differently -
            and needs *base_forcefield* set to
            :data:`~forcefill.CHARMM_BASE_FORCEFIELD`.
        atom_type: ``"gaff2"`` (default) or ``"gaff"``. gaff backend only.
        charge_method: antechamber charge method, default ``"bcc"`` (AM1-BCC).
            gaff backend only.
        smirnoff_forcefield: SMIRNOFF release for the smirnoff backend, default
            :data:`~forcefill.DEFAULT_SMIRNOFF_FORCEFIELD`. See
            ``forcefill.smirnoff.installed_smirnoff_forcefields()``.
        charmm_files: CHARMM topology/parameter files shared by every charmm
            ligand, e.g. an extra ``par_all36_cgenff.prm``. Per-ligand stream
            files go in ``LigandSpec(charmm_files=...)``, appended after these.
        workdir: Directory for intermediate files (per-residue PDB, mol2, frcmod
            and XML). A fresh temporary directory is created if not given, and
            kept unless ``cleanup=True``.
        cleanup: Delete the working directory after a *successful* build. The
            per-residue XMLs live there, so ``residue_xmls`` comes back empty
            and only the combined XML survives; on failure the directory is
            always kept for the post-mortem. Refuses to run if *output_xml*
            resolves inside it.
        validate: Verify that ``base_forcefield + output_xml`` builds a System
            for each parameterized residue on its own (default). When nothing
            was skipped, the full input topology is checked too; with skipped
            residues that check would always fail, so it is logged and omitted.
        minimize: Also energy-minimize each parameterized residue in vacuum (and
            the full topology, under the same condition as *validate*), reported
            in ``minimizations`` and ``full_minimization``. Catches unphysical
            parameters - a NaN charge, a zero force constant - that a System
            build accepts. Off by default: it costs an energy evaluation per
            residue. Subsumes *validate*, with a less specific error message
            when a template does not match.
        strict: Raise (default) rather than warn when a supplied ligand file is
            not the same molecule as the PDB residue, or a ligand has atoms on
            top of each other - before antechamber runs rather than an hour
            later. These are heuristics and the long tail is real.
        antechamber_args: Extra arguments appended to the antechamber command
            line (e.g. ``("-dr", "no")`` to relax the acdoctor checks).
            Per-ligand additions go in ``LigandSpec.antechamber_args``.
        timeout: Per-invocation ceiling in seconds for antechamber / parmchk2
            (None disables it). Default one hour. Does not apply to the smirnoff
            backend, which runs in-process.

    Returns:
        ParameterizationResult
    """
    check_choice(atom_type, ATOM_TYPES, "atom_type")
    check_choice(charge_method, CHARGE_METHODS, "charge_method")
    check_choice(backend, BACKENDS, "backend")
    ligands = dict(ligands or {})
    net_charges = dict(net_charges or {})
    multiplicities = dict(multiplicities or {})
    residue_files = dict(residue_files or {})

    pdb = app.PDBFile(str(pdb_file))
    topology, positions = pdb.topology, pdb.positions

    cleaning: CleaningResult | None = None
    if clean_structure:
        # Rebind both halves together, before anything holds a Residue: those
        # index into the topology's coordinate array, so swapping the topology
        # in later would silently address the wrong coordinates.
        topology, positions, cleaning = clean_topology(topology, positions)
        if cleaning.n_atoms_removed:
            log.warning(
                "Cleaned the structure in memory: removed %d atoms in %d "
                "residues (%s). The full-structure checks and "
                "full_minimization below describe the cleaned system, not %s.",
                cleaning.n_atoms_removed,
                cleaning.n_residues_removed,
                ", ".join(sorted(cleaning.removed)),
                pdb_file,
            )

    unmatched = find_nonstandard_residues(topology, base_forcefield)
    if not unmatched:
        log.info("All residues matched %s - nothing to parameterize.", list(base_forcefield))
        return ParameterizationResult(forcefield_xml=None, cleaning=cleaning)

    to_param, skipped = _classify_unmatched(unmatched)
    for name, reason in skipped.items():
        log.warning("Skipping %s: %s", name, reason)
    if not to_param:
        details = "\n".join(f"  {n}: {r}" for n, r in skipped.items())
        raise RuntimeError(f"Unmatched residues were found but none can be auto-parameterized:\n{details}")
    log.info("Residues to parameterize: %s", sorted(to_param))
    _warn_unused_overrides(
        to_param,
        skipped,
        net_charges,
        multiplicities,
        residue_files,
        removed=cleaning.removed if cleaning else (),
        ligands=ligands,
    )

    specs = resolve_specs(
        ligands,
        net_charges=net_charges,
        multiplicities=multiplicities,
        residue_files=residue_files,
        defaults=_Defaults(
            atom_type=atom_type,
            charge_method=charge_method,
            antechamber_args=tuple(antechamber_args),
            backend=backend,
            forcefield=smirnoff_forcefield,
            charmm_files=tuple(charmm_files),
            names=frozenset(to_param),
        ),
    )

    # Fail early if the tools a backend needs are absent, or if the base force
    # field is one its output could never be loaded with.
    gaff_dat = _pipeline.prepare_gaff_backend(specs, atom_type)
    _pipeline.check_backends_match_base(specs, base_forcefield)

    minimizations: dict[str, MinimizationResult] = {}
    full_minimization: MinimizationResult | None = None
    with _pipeline.working_directory(workdir, output_xml, prefix="nonstandard_ff_", cleanup=cleanup) as wd:
        # Check everything before the first expensive call: AM1-BCC can take an
        # hour per ligand, so a mistake in the last must not cost the first.
        specs = preflight_specs(specs, to_param, positions, wd, strict=strict, base_forcefield=base_forcefield)

        artifacts = _pipeline.parameterize_all(
            specs, to_param, positions, wd, gaff_dat=gaff_dat, timeout=timeout, base_forcefield=base_forcefield
        )
        residue_xmls = {name: art.xml for name, art in artifacts.items()}

        combined = _pipeline.combine_residue_xmls(artifacts, specs, gaff_dat, output_xml, wd)
        log.info("Wrote combined force-field XML: %s", combined)

        if validate or minimize:
            # Parse the (large) base force field + new XML once; every check uses it.
            files = [*base_forcefield, combined]
            forcefield = app.ForceField(*files)
            # Validate first: a template mismatch then reports itself as such,
            # not as the minimizer failing to build a System.
            if validate:
                _validate_parameterized_residues(to_param, forcefield, files)
            if minimize:
                minimizations = _minimize_parameterized_residues(
                    to_param, positions, combined, base_forcefield, forcefield
                )
            if skipped:
                log.warning(
                    "Skipping the full-structure checks: %d residue type(s) "
                    "were skipped (%s) and still have no template, so building "
                    "a System for the whole input cannot succeed. Repair or "
                    "parameterize those, then re-check with "
                    "validate_forcefield_xml() / minimize_with_forcefield_xml().",
                    len(skipped),
                    ", ".join(sorted(skipped)),
                )
            else:
                if validate:
                    validate_forcefield_xml(topology, combined, base_forcefield, forcefield=forcefield)
                if minimize:
                    full_minimization = minimize_with_forcefield_xml(
                        topology, positions, combined, base_forcefield, forcefield=forcefield
                    )

    return ParameterizationResult(
        forcefield_xml=combined,
        residue_xmls={} if cleanup else residue_xmls,
        parameterized=sorted(artifacts),
        skipped=skipped,
        workdir=None if cleanup else str(wd),
        minimizations=minimizations,
        full_minimization=full_minimization,
        cleaning=cleaning,
    )
