"""The machinery :func:`~forcefill.build_forcefield_xml` and :func:`~forcefill.build_ligand_xml` share.

The two entry points differ only in where the ligands come from - a PDB's
unmatched residues, or the caller's own list. Everything after that is the same
work and lives here: preparing a backend, owning the working directory, running
one residue through to a per-residue XML, and combining those into one file.
Each entry point then just resolves its input into ``{name: ResolvedSpec}``,
calls in here, and assembles a :class:`ParameterizationResult`.

AmberTools is reached as ``amber.run_antechamber(...)`` rather than through a
from-import, deliberately: a from-import binds a *copy* at import time, which a
test stubbing ``forcefill.amber`` could never reach.
"""

from __future__ import annotations

import logging
import shutil
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

from openmm import app, unit

from . import amber, charmm, smirnoff
from ._spec import CHARMM_BASE_FORCEFIELD, DEFAULT_BASE_FORCEFIELD, PathLike, ResolvedSpec
from .checks import MinimizationResult
from .clean_structure import CleaningResult
from .merge import _SCALE_TOLERANCE, merge_ffxml
from .topology import extract_residue_to_pdb

log = logging.getLogger(__name__)

__all__ = ["ParameterizationResult"]

#: 1-4 scaling each backend's output declares, compared against the base force
#: field by :func:`check_backends_match_base`. The Amber pair is what ParmEd
#: writes from ``gaff*.dat`` and openmmforcefields for SMIRNOFF; the CHARMM pair
#: is what ParmEd writes from a ``CharmmParameterSet``.
_BACKEND_14_SCALES = {
    "gaff": (0.8333333333333334, 0.5),
    "smirnoff": (0.8333333333333334, 0.5),
    "charmm": (1.0, 1.0),
}


@dataclass
class ParameterizationResult:
    """What :func:`~forcefill.build_forcefield_xml` produced.

    Also the return type of :func:`~forcefill.build_ligand_xml`, where
    ``skipped``, ``cleaning`` and ``full_minimization`` are always empty: with no
    input structure there is nothing to skip, clean or minimize as a whole.
    """

    #: Path to the combined ffxml covering every parameterized residue
    #: (``None`` when nothing needed parameterizing).
    forcefield_xml: str | None
    #: Per-residue ffxml files, keyed by residue name (empty after
    #: ``cleanup=True`` - they live in the working directory).
    residue_xmls: dict[str, str] = field(default_factory=dict)
    #: Residue names that were successfully parameterized.
    parameterized: list[str] = field(default_factory=list)
    #: Residue names that were skipped, mapped to the reason.
    skipped: dict[str, str] = field(default_factory=dict)
    #: Directory holding intermediate files (per-residue PDB/mol2/frcmod);
    #: ``None`` when nothing was parameterized or after ``cleanup=True``.
    workdir: str | None = None
    #: Per-residue vacuum minimizations, keyed by residue name. Empty unless
    #: ``minimize=True``.
    minimizations: dict[str, MinimizationResult] = field(default_factory=dict)
    #: Minimization of the whole input topology - the *cleaned* topology when
    #: ``clean_structure=True``, so reconcile ``n_atoms`` against
    #: ``cleaning.n_atoms_after``. ``None`` unless ``minimize=True`` *and* no
    #: residue was skipped.
    full_minimization: MinimizationResult | None = None
    #: What ``clean_structure=True`` removed from the input. ``None`` when it
    #: was off, in which case nothing was deleted.
    cleaning: CleaningResult | None = None


@dataclass
class ResidueArtifacts:
    """What parameterizing one residue produced.

    ``mol2``/``frcmod`` are None for the smirnoff backend, which has no
    intermediate Amber files - only the finished XML is common to both.
    """

    xml: str
    mol2: str | None = None
    frcmod: str | None = None


def prepare_gaff_backend(specs: Mapping[str, ResolvedSpec], atom_type: str) -> str | None:
    """Check the gaff backend can run and return the GAFF database, or None if unused.

    Called before the working directory exists, so a missing AmberTools install
    fails immediately rather than after the first ligand has been read.
    """
    if not any(spec.backend == "gaff" for spec in specs.values()):
        return None
    amber.require_executable("antechamber")
    amber.require_executable("parmchk2")
    gaff_dat = amber.locate_gaff_dat(atom_type)
    log.info("Using GAFF parameter database: %s", gaff_dat)
    return gaff_dat


def check_backends_match_base(specs: Mapping[str, ResolvedSpec], base_forcefield: Sequence[str]) -> None:
    """Refuse a combination OpenMM could never load, before anything expensive runs.

    Amber-family force fields scale 1-4 interactions by 0.8333/0.5 and CHARMM by
    1.0/1.0, and OpenMM rejects a ``ForceField`` whose files disagree. Two
    combinations are therefore impossible rather than inadvisable, and are worth
    naming here instead of an hour into a build:

        * a charmm ligand mixed with a gaff or smirnoff one, whose merged XML
          could not be loaded at all;
        * a backend whose output does not match the base force field it is
          validated against.

    The base convention is read from the loaded force field, so a custom one is
    checked as accurately as the two presets.
    """
    backends = {spec.backend for spec in specs.values()}
    if not backends:
        return
    charmm_names = sorted(name for name, spec in specs.items() if spec.backend == "charmm")
    amber_names = sorted(name for name, spec in specs.items() if spec.backend != "charmm")
    if charmm_names and amber_names:
        raise ValueError(
            f"Cannot build one force field from both CHARMM and Amber-family "
            f"parameters: {charmm_names} use the charmm backend and {amber_names} "
            f"use {sorted(backends - {'charmm'})}. The two conventions scale 1-4 "
            "interactions differently (CHARMM 1.0/1.0, Amber 0.8333/0.5) and "
            "OpenMM will not load a force field that says both. Build them "
            "separately, against their own base force fields."
        )

    # Only one family is in play now, and gaff and smirnoff share a convention,
    # so any one backend answers for all of them.
    expected = _BACKEND_14_SCALES["charmm" if charmm_names else "gaff"]
    actual = charmm.base_14_scales(base_forcefield)
    # None: the base force field declares no non-bonded terms at all, so there is
    # nothing for the generated XML to contradict.
    if actual is None or _scales_agree(expected, actual):
        return
    wanted = CHARMM_BASE_FORCEFIELD if charmm_names else DEFAULT_BASE_FORCEFIELD
    raise ValueError(
        f"The {'/'.join(sorted(backends))} backend produces parameters with 1-4 "
        f"scaling {expected[0]:g}/{expected[1]:g} (coulomb/lj), but the base force "
        f"field {list(base_forcefield)} declares {actual[0]:g}/{actual[1]:g}. OpenMM "
        "cannot load the two together, so the generated XML would be unusable "
        f"even though it built. Pass base_forcefield={list(wanted)}, or switch "
        "backend to match the base force field you want."
    )


def _scales_agree(expected: tuple[float, float], actual: tuple[float, float]) -> bool:
    """Compare 1-4 scales with the tolerance OpenMM itself applies when merging them."""
    return all(abs(a - b) <= _SCALE_TOLERANCE for a, b in zip(expected, actual, strict=True))


@contextmanager
def working_directory(
    workdir: PathLike | None,
    output_xml: PathLike,
    *,
    prefix: str,
    cleanup: bool,
) -> Iterator[Path]:
    """Own the intermediate-file directory for one build.

    Creates *workdir* (or a fresh temporary directory named *prefix*), yields
    it, then keeps it on failure - ``sqm.out`` and the intermediates are the
    post-mortem - and removes it on success only with ``cleanup=True``.

    Refuses up front if ``cleanup`` would delete *output_xml* along with the
    directory, which is otherwise a silently empty result.
    """
    workdir = Path(workdir).resolve() if workdir is not None else Path(tempfile.mkdtemp(prefix=prefix))
    workdir.mkdir(parents=True, exist_ok=True)
    log.info("Intermediate files in %s", workdir)
    if cleanup and Path(output_xml).resolve().is_relative_to(workdir):
        raise ValueError(
            f"cleanup=True would delete the output XML: {output_xml} resolves "
            f"inside the working directory {workdir}. Write it elsewhere or "
            "pass cleanup=False."
        )

    try:
        yield workdir
    except Exception:
        log.warning("Intermediate files kept for debugging in %s", workdir)
        raise

    if cleanup:
        shutil.rmtree(workdir)
        log.info("Removed working directory %s", workdir)


def parameterize_one_residue(
    spec: ResolvedSpec,
    residue: app.topology.Residue | None,
    positions: unit.Quantity | None,
    res_dir: Path,
    *,
    gaff_dat: str | None = None,
    timeout: float | None = amber.DEFAULT_AMBERTOOLS_TIMEOUT,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
) -> ResidueArtifacts:
    """Run one residue through its backend to a per-residue XML.

    For ``gaff`` that is extract -> antechamber -> parmchk2 -> ParmEd, with a
    ``spec.file`` (SDF/MOL2 with explicit bonds) replacing the extraction step.
    For ``smirnoff`` it is one call into openmmforcefields; for ``charmm``, a
    conversion of the ligand's CGenFF files - the one output that is *not*
    self-contained, since it names atom types *base_forcefield* defines rather
    than redefining them. *residue* and *positions* are None in standalone mode,
    where there is no structure to extract from.
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    name = spec.name

    if spec.backend == "smirnoff":
        return ResidueArtifacts(xml=smirnoff.smirnoff_residue_ffxml(spec, res_dir / f"{name}.xml"))

    if spec.backend == "charmm":
        return ResidueArtifacts(xml=charmm.charmm_residue_ffxml(spec, res_dir / f"{name}.xml", base_forcefield))

    if spec.file is None:
        if residue is None or positions is None:
            raise ValueError(f"Residue {name} has no ligand file and no structure to extract it from.")
        antechamber_input: PathLike = extract_residue_to_pdb(positions, residue, res_dir / f"{name}.pdb")
    else:
        log.info("Using the supplied ligand file for %s: %s", name, spec.file)
        antechamber_input = spec.file

    # Explicitly against None: a net charge of 0 is a real, stated value, and
    # `or` would quietly conflate it with "not determined".
    net_charge = spec.net_charge if spec.net_charge is not None else 0
    log.info("antechamber: %s (net charge %+d, %s/%s)", name, net_charge, spec.atom_type, spec.charge_method)
    mol2 = amber.run_antechamber(
        antechamber_input,
        res_dir / f"{name}.mol2",
        name,
        net_charge=net_charge,
        multiplicity=spec.multiplicity,
        atom_type=spec.atom_type,
        charge_method=spec.charge_method,
        extra_args=spec.antechamber_args,
        timeout=timeout,
    )
    frcmod = amber.run_parmchk2(mol2, res_dir / f"{name}.frcmod", atom_type=spec.atom_type, timeout=timeout)
    # Per-residue template XML (self-contained).
    gaff_dat = gaff_dat or amber.locate_gaff_dat(spec.atom_type)
    xml = amber.assemble_openmm_ffxml({name: mol2}, [gaff_dat, frcmod], res_dir / f"{name}.xml")
    log.info("Wrote per-residue XML: %s", xml)
    return ResidueArtifacts(xml=xml, mol2=mol2, frcmod=frcmod)


def combine_residue_xmls(
    artifacts: Mapping[str, ResidueArtifacts],
    specs: Mapping[str, ResolvedSpec],
    gaff_dat: str | None,
    output_xml: PathLike,
    workdir: Path,
) -> str:
    """Write the one XML covering every parameterized residue.

    All-GAFF goes through ParmEd, which merges at the parameter-set level and
    writes only the atom types the templates actually use. Anything else is
    merged as finished XML instead - gaff, smirnoff and charmm share nothing
    upstream of that.
    """
    gaff = {name for name, spec in specs.items() if spec.backend == "gaff"}
    other_names = sorted(set(specs) - gaff)
    if not other_names:
        return amber.assemble_openmm_ffxml(
            {name: artifacts[name].mol2 for name in sorted(gaff)},
            [gaff_dat, *(artifacts[name].frcmod for name in sorted(gaff))],
            output_xml,
        )

    to_merge: list[PathLike] = []
    if gaff:
        # One ParmEd document for all the GAFF residues, so their shared atom
        # types are written once, then merged with the SMIRNOFF ones.
        to_merge.append(
            amber.assemble_openmm_ffxml(
                {name: artifacts[name].mol2 for name in sorted(gaff)},
                [gaff_dat, *(artifacts[name].frcmod for name in sorted(gaff))],
                workdir / "_gaff_combined.xml",
            )
        )
    to_merge += [artifacts[name].xml for name in other_names]
    return merge_ffxml(to_merge, output_xml)


def parameterize_all(
    specs: Mapping[str, ResolvedSpec],
    residues: Mapping[str, app.topology.Residue],
    positions: unit.Quantity | None,
    workdir: Path,
    *,
    gaff_dat: str | None,
    timeout: float | None,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
) -> dict[str, ResidueArtifacts]:
    """Run every spec through its backend, in a stable order."""
    return {
        name: parameterize_one_residue(
            specs[name],
            residues.get(name),
            positions,
            workdir / name,
            gaff_dat=gaff_dat,
            timeout=timeout,
            base_forcefield=base_forcefield,
        )
        for name in sorted(specs)
    }
