"""Identify the non-standard residues in a structure and build a force field for them.

The entry point that starts from a PDB. :func:`~forcefill.build_ligand_xml`
starts from the ligands themselves; the pipeline between the two is shared, and
lives in :mod:`forcefill._pipeline`.

Pipeline:
    1. :func:`~forcefill.find_nonstandard_residues` asks the base force field
       which residues it cannot match.
    2. Those are classified, and only chemistry a stand-alone GAFF treatment is
       actually *valid* for is parameterized:

       * standard residues that are merely missing atoms   -> reported, skipped
         (repair the structure with PDBFixer / ``Modeller.addHydrogens`` instead)
       * monatomic species (ions)                          -> reported, skipped
         (load an ion parameter file; never run antechamber on a bare ion)
       * residues covalently bonded to neighbours          -> reported, skipped
         (stand-alone GAFF is wrong for polymer-linked residues such as modified
         amino acids; those need capping + a consistent charge derivation, e.g.
         with pyRED or ffparam-style workflows)
       * free-standing hetero molecules (ligands, cofactors) -> parameterized

    3. Everything checkable is checked before anything expensive runs -
       see :mod:`forcefill.preflight`.
    4. Each unique residue goes through its backend: AmberTools for ``gaff``
       (:mod:`forcefill.amber`), openff-toolkit for ``smirnoff``
       (:mod:`forcefill.smirnoff`).
    5. The per-residue results are combined into one XML.
    6. That XML is validated, and optionally minimized, by
       :mod:`forcefill.checks`.

Requirements:
    * ``openmm``, ``parmed``, ``rdkit``, ``openff-toolkit`` and
      ``openmmforcefields`` at import time
    * AmberTools (``antechamber``, ``parmchk2``) on ``PATH`` for the gaff
      backend, e.g. ``conda install -c conda-forge ambertools``

    conda-forge is the recommended route for the whole stack; see
    ``environment.yml``.

Example:
    >>> from forcefill import build_forcefield_xml
    >>> result = build_forcefield_xml("complex.pdb", "extras.xml", net_charges={"LIG": -1})
    >>> result.parameterized
    ['LIG']

    then simulate with::

        ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
        system = ff.createSystem(pdb.topology, ...)

Notes:
    * Crystal structures carry water, buffer ions and crystallization
      additives, and a free-standing additive such as glycerol looks exactly
      like a ligand to step 2 - so it gets parameterized. Strip them first with
      :func:`forcefill.clean_pdb`, or pass ``clean_structure=True`` here to do
      it in memory.
    * Ligands must contain **all explicit hydrogens** with reasonable geometry;
      AM1-BCC charges are meaningless otherwise.
    * Load either the combined XML *or* the per-residue XMLs with a ForceField,
      never both at once (duplicate GAFF atom-type definitions would collide).
    * If you would rather not manage XML files at all, the same job can be done
      at runtime with ``openmmforcefields.generators.GAFFTemplateGenerator``.
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
        clean_structure: If True, remove crystallographic water, bulk
            counter-ions and crystallization additives *in memory* before
            anything else looks at the structure, using
            :func:`~forcefill.clean_topology` with its defaults - so structural
            metals are kept. Off by default: it deletes atoms, and a pipeline
            should never do that unasked. What went is reported in
            ``cleaning``. Note that everything downstream then describes the
            *cleaned* system rather than the file on disk: the full-structure
            ``validate`` / ``minimize`` checks, and therefore
            ``full_minimization.n_atoms``, refer to the stripped topology. For
            anything beyond the defaults - keeping a particular additive,
            stripping a metal - call :func:`~forcefill.clean_pdb` yourself and
            pass the cleaned file in.
        base_forcefield: ffxml files defining what counts as "standard".
        ligands: ``{residue_name: LigandSpec}`` - per-ligand settings. A spec
            says where the ligand comes from (``file``, ``smiles``) and how to
            treat it (``net_charge``, ``backend``, ``atom_type``, ...), and
            leaves everything it does not set to the call-level defaults below.
            This is the general form of *net_charges* / *multiplicities* /
            *residue_files*, which still work and are folded in; setting the
            same thing both ways raises.
        net_charges: ``{residue_name: net_charge}``. When a ligand is supplied
            as a file or SMILES its formal charge is read from there, so this is
            only needed for a residue extracted from the PDB - or to override
            what the file says, which raises if the two disagree.
        multiplicities: ``{residue_name: spin_multiplicity}``; defaults to 1.
        residue_files: ``{residue_name: ligand_file}`` - parameterize these
            residues from the given SDF/MOL2 (or other antechamber-readable)
            file instead of extracting them from the PDB. A file with explicit
            bond orders and protonation as drawn avoids antechamber having to
            re-perceive them from PDB geometry, a classic source of silently
            wrong atom types. The file must contain the *same* atoms and bonds
            (including hydrogens) as the residue in the PDB; *strict* checks
            that up front.
        backend: ``"gaff"`` (default) for GAFF/GAFF2 through antechamber, or
            ``"smirnoff"`` for OpenFF through openff-toolkit. SMIRNOFF assigns
            parameters from the chemical graph, so every ligand on it needs a
            ``file`` or ``smiles`` - a PDB residue carries no bond orders. Set
            it per ligand with ``LigandSpec(backend=...)`` to mix the two.
        atom_type: ``"gaff2"`` (default) or ``"gaff"``. gaff backend only.
        charge_method: antechamber charge method, default ``"bcc"`` (AM1-BCC).
            gaff backend only.
        smirnoff_forcefield: SMIRNOFF release for the smirnoff backend, default
            :data:`~forcefill.DEFAULT_SMIRNOFF_FORCEFIELD`. See
            ``forcefill.smirnoff.installed_smirnoff_forcefields()``.
        workdir: Directory for intermediate files (per-residue PDB, mol2,
            frcmod, per-residue XML). A fresh temporary directory is created
            if not given; its path is reported in the result and it is kept
            unless ``cleanup=True``.
        cleanup: If True, delete the working directory after a *successful*
            build. The per-residue XMLs live there, so ``residue_xmls`` comes
            back empty and ``workdir`` is None; only the combined XML
            survives. On failure the directory is always kept so ``sqm.out``
            and the intermediates can be inspected. Refuses to run if
            *output_xml* itself resolves inside the working directory.
        validate: If True (default), verify that ``base_forcefield +
            output_xml`` can build a System for each parameterized residue on
            its own. When no residues were skipped, additionally verify the
            full input topology (the cleaned one, with *clean_structure*);
            with skipped residues present that check would always fail (they
            still have no template), so it is logged and omitted instead.
        minimize: If True, additionally energy-minimize each parameterized
            residue in vacuum (and the full input topology, under the same
            no-residues-skipped condition as *validate*), reported back in
            ``minimizations`` and ``full_minimization``. This catches
            unphysical parameters - a NaN charge, a zero force constant -
            which building a System alone cannot. Off by default because it
            costs an energy evaluation per residue. It subsumes *validate*:
            minimizing implies building the System, just with a less specific
            error message when a template does not match.
        strict: If True (default), a supplied ligand file that is not the same
            molecule as the residue in the PDB, or a ligand with atoms on top of
            each other, is an error - raised before antechamber runs rather than
            discovered an hour later. Set False to downgrade both to warnings;
            these are heuristics and the long tail is real.
        antechamber_args: Extra raw arguments appended to the antechamber
            command line (e.g. ``("-dr", "no")`` to relax acdoctor structure
            checks). Per-ligand additions go in ``LigandSpec.antechamber_args``.
        timeout: Per-invocation ceiling in seconds for each antechamber /
            parmchk2 run (None disables it). Default one hour. Does not apply to
            the smirnoff backend, which runs in-process.

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
        # Both halves rebound together, before anything holds a Residue
        # reference: _classify_unmatched stores live Residue objects and
        # _residue_positions indexes positions by atom.index, so a topology
        # swapped in later would silently address the wrong coordinates.
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
            names=frozenset(to_param),
        ),
    )

    # Fail early if the tools a backend needs are absent.
    gaff_dat = _pipeline.prepare_gaff_backend(specs, atom_type)

    minimizations: dict[str, MinimizationResult] = {}
    full_minimization: MinimizationResult | None = None
    with _pipeline.working_directory(workdir, output_xml, prefix="nonstandard_ff_", cleanup=cleanup) as wd:
        # Everything checkable is checked before the first expensive call, so a
        # mistake in the last ligand does not cost the parameterization of the
        # first - antechamber's AM1-BCC can take an hour per ligand.
        specs = preflight_specs(specs, to_param, positions, wd, strict=strict)

        artifacts = _pipeline.parameterize_all(specs, to_param, positions, wd, gaff_dat=gaff_dat, timeout=timeout)
        residue_xmls = {name: art.xml for name, art in artifacts.items()}

        combined = _pipeline.combine_residue_xmls(artifacts, specs, gaff_dat, output_xml, wd)
        log.info("Wrote combined force-field XML: %s", combined)

        if validate or minimize:
            # Parse the (large) base force field + new XML once; every check uses it.
            files = [*base_forcefield, combined]
            forcefield = app.ForceField(*files)
            # Validate first: a template mismatch then reports itself as a bond-graph
            # problem rather than as the minimizer failing to build a System.
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
