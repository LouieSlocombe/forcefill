"""Identify non-standard residues in a PDB with OpenMM and build a ready-to-use force-field XML for them.

Pipeline:
    1. ``ForceField.getUnmatchedResidues`` finds every residue the chosen base
       force field cannot match.
    2. Unmatched residues are classified:

       * standard residues that are merely missing atoms   -> reported, skipped
         (repair the structure with PDBFixer / ``Modeller.addHydrogens`` instead)
       * monatomic species (ions)                          -> reported, skipped
         (load an ion parameter file; never run antechamber on a bare ion)
       * residues covalently bonded to neighbours          -> skipped by default
         (stand-alone GAFF treatment is wrong for polymer-linked residues such
         as modified amino acids; those need capping + a consistent charge
         derivation, e.g. with pyRED or ffparam-style workflows)
       * free-standing hetero molecules (ligands, cofactors) -> parameterized

    3. Everything checkable is checked, before anything expensive runs: the net
       charge is read from any supplied ligand file, that file is confirmed to
       describe the same molecule as the residue, and the geometry is checked
       for the faults that produce NaN energies. See :func:`preflight_specs`.

    4. Each unique parameterizable residue goes through its backend. With
       ``backend="gaff"`` (the default) that is AmberTools:

       * ``antechamber``  assigns GAFF/GAFF2 atom types + AM1-BCC charges
         (-> ``<RES>.mol2`` residue template)
       * ``parmchk2``     generates any GAFF parameters that are missing
         (-> ``<RES>.frcmod``)

       With ``backend="smirnoff"`` it is openff-toolkit and openmmforcefields
       instead; see :mod:`forcefill.smirnoff`.

    5. ParmEd merges the GAFF parameter database, the frcmod files and the mol2
       templates into OpenMM ffxml: one XML per residue plus one combined XML
       containing every residue (atom types, residue templates with charges,
       bonds/angles/torsions and nonbonded parameters). Where the smirnoff
       backend is involved the combined file is produced by
       :func:`merge_ffxml`, which merges finished XML rather than parameter
       sets - the two backends share nothing upstream of that.

    6. The combined XML is validated by building an ``openmm.System`` from
       ``base force field + new XML`` for each parameterized residue on its
       own; when no residues were skipped, a System for the full input
       topology is built as well. With ``minimize=True`` each of those is also
       energy-minimized, which catches unphysical parameters that still build
       a perfectly valid System.

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

    To parameterize a ligand with no structure at all, use
    :func:`forcefill.build_ligand_xml`.

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
import math
import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import parmed
from openmm import Context, LocalEnergyMinimizer, Platform, VerletIntegrator, app, unit
from parmed.amber import AmberParameterSet
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.openmm import OpenMMParameterSet

from . import ligand_files, smirnoff
from ._residue_names import STANDARD_RESIDUES
from ._spec import ATOM_TYPES as _ATOM_TYPES
from ._spec import BACKENDS, DEFAULT_SMIRNOFF_FORCEFIELD, LigandSpec, ResolvedSpec, _Defaults, resolve_specs
from ._spec import CHARGE_METHODS as _CHARGE_METHODS
from ._spec import check_choice as _check_choice
from .clean_structure import CleaningResult, clean_topology

log = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_AMBERTOOLS_TIMEOUT",
    "DEFAULT_BASE_FORCEFIELD",
    "DEFAULT_MINIMIZATION_PLATFORM",
    "DEFAULT_MINIMIZATION_TOLERANCE",
    "MinimizationResult",
    "ParameterizationResult",
    "assemble_openmm_ffxml",
    "build_forcefield_xml",
    "extract_residue_to_pdb",
    "find_nonstandard_residues",
    "locate_gaff_dat",
    "merge_ffxml",
    "minimize_with_forcefield_xml",
    "run_antechamber",
    "run_parmchk2",
    "validate_forcefield_xml",
]

PathLike = str | os.PathLike

#: Base force field used to decide what counts as "non-standard".
DEFAULT_BASE_FORCEFIELD = ("amber14-all.xml", "amber14/tip3p.xml")

#: Ceiling for a single AmberTools invocation, in seconds. sqm's AM1-BCC on a
#: large ligand can legitimately take many minutes; nothing should take an hour.
DEFAULT_AMBERTOOLS_TIMEOUT: float = 3600.0

#: Platform used for the minimization checks. Pinned rather than left to OpenMM,
#: which picks the fastest available one - forcefill must never take a GPU the
#: caller wanted for something else.
DEFAULT_MINIMIZATION_PLATFORM = "CPU"

#: Minimizer convergence target on the RMS force, in kJ/mol/nm (OpenMM's own default).
DEFAULT_MINIMIZATION_TOLERANCE: float = 10.0

#: Unit the reported forces are converted to.
_FORCE_UNIT = unit.kilojoule_per_mole / unit.nanometer

#: File suffix -> antechamber -fi format, for run_antechamber's inference.
_ANTECHAMBER_FORMATS = {
    ".pdb": "pdb",
    ".mol2": "mol2",
    ".sdf": "sdf",
    ".sd": "sdf",
    ".mol": "mdl",
}


@dataclass
class MinimizationResult:
    """What one :func:`minimize_with_forcefield_xml` run measured.

    Plain floats rather than ``unit.Quantity``: energies in kJ/mol, forces in
    kJ/mol/nm. A Quantity cannot be ``math.isfinite``-checked or ``%``-formatted,
    which is all this object is ever used for.
    """

    #: Number of atoms in the minimized topology.
    n_atoms: int
    #: Potential energy of the input coordinates.
    initial_energy: float
    #: Potential energy after minimizing.
    final_energy: float
    #: Largest per-atom force magnitude after minimizing.
    max_force: float

    @property
    def energy_change(self) -> float:
        """Final minus initial energy; negative when the minimizer did its job."""
        return self.final_energy - self.initial_energy


@dataclass
class ParameterizationResult:
    """What :func:`build_forcefield_xml` produced.

    Also the return type of :func:`~forcefill.build_ligand_xml`, where
    ``skipped``, ``cleaning`` and ``full_minimization`` are always empty: with no
    input structure there is nothing to skip residues from, clean, or minimize
    as a whole.
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


# --------------------------------------------------------------------------
# Step 1: identification
# --------------------------------------------------------------------------


def find_nonstandard_residues(
    topology: app.Topology,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
) -> list[app.topology.Residue]:
    """Return every residue that *base_forcefield* has no template for.

    This is OpenMM's own definition of "non-standard": a residue whose
    element/bond graph matches no registered template. Note that standard
    residues with missing atoms (e.g. a protein without hydrogens) also fail
    to match; :func:`build_forcefield_xml` filters those out separately.
    """
    forcefield = app.ForceField(*base_forcefield)
    return forcefield.getUnmatchedResidues(topology)


def _classify_unmatched(
    unmatched: Sequence[app.topology.Residue],
) -> tuple[dict[str, app.topology.Residue], dict[str, str]]:
    """Split unmatched residues into those to parameterize and those to skip.

    Returns ``({name: representative_residue}, {name: skip_reason})``.
    """
    groups: dict[str, list[app.topology.Residue]] = defaultdict(list)
    for res in unmatched:
        groups[res.name.strip()].append(res)

    to_param: dict[str, app.topology.Residue] = {}
    skipped: dict[str, str] = {}
    for name, residues in groups.items():
        counts = sorted({sum(1 for _ in r.atoms()) for r in residues})
        rep = max(residues, key=lambda r: sum(1 for _ in r.atoms()))
        n_atoms = sum(1 for _ in rep.atoms())
        if len(counts) > 1:
            log.warning(
                "Copies of residue %s differ in atom count (%s); "
                "using the most complete copy (%d atoms) as the template.",
                name,
                counts,
                n_atoms,
            )
        if name in STANDARD_RESIDUES:
            skipped[name] = (
                f"standard residue that failed to match ({len(residues)} "
                "copies) - it is probably missing atoms or has non-standard "
                "atom names; repair the structure (e.g. with PDBFixer or "
                "Modeller.addHydrogens) instead of reparameterizing it"
            )
        elif n_atoms == 1:
            skipped[name] = (
                "monatomic species - use an ion parameter file for it (GAFF/antechamber cannot treat bare ions)"
            )
        elif n_linked := sum(1 for r in residues if any(True for _ in r.external_bonds())):
            skipped[name] = (
                f"covalently bonded to neighbouring residues ({n_linked} of "
                f"{len(residues)} copies) - a stand-alone GAFF "
                "parameterization is not valid for polymer-linked residues; "
                "cap the fragment and derive charges consistently with the "
                "backbone force field instead"
            )
        else:
            to_param[name] = rep
    return to_param, skipped


def _warn_unused_overrides(
    to_param: Mapping[str, app.topology.Residue],
    skipped: Mapping[str, str],
    net_charges: Mapping[str, int],
    multiplicities: Mapping[str, int],
    residue_files: Mapping[str, PathLike] | None = None,
    removed: Iterable[str] = (),
    ligands: Mapping[str, LigandSpec] | None = None,
) -> None:
    """Warn about ligands/net_charges/multiplicities/residue_files keys with no effect.

    A typo'd or case-mismatched key silently leaves the defaults (net
    charge 0, multiplicity 1, PDB extraction), which yields plausible but
    wrong AM1-BCC charges - the worst failure mode.

    *removed* names the residues cleaning deleted, so an override aimed at one
    of them reports the real cause rather than "matches no residue".
    """
    removed = set(removed)
    for label, mapping in (
        ("ligands", ligands or {}),
        ("net_charges", net_charges),
        ("multiplicities", multiplicities),
        ("residue_files", residue_files or {}),
    ):
        for key in mapping:
            if key in to_param:
                continue
            if key in removed:
                log.warning(
                    "%s[%r] has no effect: residue %s was removed from the structure by clean_structure=True.",
                    label,
                    key,
                    key,
                )
            elif key in skipped:
                log.warning(
                    "%s[%r] has no effect: residue %s is being skipped, not parameterized.",
                    label,
                    key,
                    key,
                )
            else:
                log.warning(
                    "%s[%r] does not match any residue selected for "
                    "parameterization %s; check the spelling and case of "
                    "the residue name.",
                    label,
                    key,
                    sorted(to_param),
                )


# --------------------------------------------------------------------------
# Step 2: extraction
# --------------------------------------------------------------------------


def _residue_subtopology(residue: app.topology.Residue) -> app.Topology:
    """Copy *residue* (atoms and internal bonds) into a fresh Topology."""
    sub_top = app.Topology()
    chain = sub_top.addChain("A")
    new_res = sub_top.addResidue(residue.name, chain)
    atom_map = {}
    for atom in residue.atoms():
        atom_map[atom] = sub_top.addAtom(atom.name, atom.element, new_res)
    for bond in residue.internal_bonds():
        sub_top.addBond(atom_map[bond.atom1], atom_map[bond.atom2])
    return sub_top


def _residue_positions(positions: unit.Quantity, residue: app.topology.Residue) -> unit.Quantity:
    """Slice *positions* (indexed by global atom index) down to *residue*, in its own atom order.

    The order must match :func:`_residue_subtopology`, which iterates the same
    ``residue.atoms()``.
    """
    return unit.Quantity(
        [positions[a.index].value_in_unit(unit.nanometer) for a in residue.atoms()],
        unit.nanometer,
    )


def extract_residue_to_pdb(
    positions: unit.Quantity,
    residue: app.topology.Residue,
    out_pdb: PathLike,
) -> str:
    """Write a single residue's atoms, internal bonds and coordinates to *out_pdb* and return the path."""
    for atom in residue.atoms():
        if atom.element is None:
            log.warning(
                "Atom %s in residue %s has no element assigned; antechamber "
                "may misread it. Check the element columns of the PDB.",
                atom.name,
                residue.name,
            )
    sub_top = _residue_subtopology(residue)

    out_pdb = Path(out_pdb)
    with open(out_pdb, "w") as fh:
        app.PDBFile.writeFile(sub_top, _residue_positions(positions, residue), fh)
    return str(out_pdb)


# --------------------------------------------------------------------------
# Step 3: AmberTools wrappers
# --------------------------------------------------------------------------


def _require_executable(name: str) -> str:
    exe = shutil.which(name)
    if exe is None:
        raise RuntimeError(
            f"'{name}' was not found on PATH. It is part of AmberTools; "
            "install it with e.g. 'conda install -c conda-forge ambertools' "
            "and try again."
        )
    return exe


def _tail(stream: str | bytes | None, limit: int = 2000) -> str:
    """Last *limit* characters of captured output; TimeoutExpired may carry bytes or None."""
    if stream is None:
        return ""
    if isinstance(stream, bytes):
        stream = stream.decode(errors="replace")
    return stream[-limit:]


def _run(
    cmd: Sequence[str],
    cwd: PathLike,
    *,
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
    hint: str = "",
) -> None:
    log.debug("Running: %s (cwd=%s)", " ".join(map(str, cmd)), cwd)
    argv = [str(c) for c in cmd]
    try:
        proc = subprocess.run(argv, cwd=str(cwd), capture_output=True, text=True, check=False, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Command timed out after {timeout:g} s:\n"
            f"  {' '.join(argv)}\n"
            f"--- stdout (tail) ---\n{_tail(exc.stdout)}\n"
            f"--- stderr (tail) ---\n{_tail(exc.stderr)}\n" + hint
        ) from exc
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}:\n"
            f"  {' '.join(argv)}\n"
            f"--- stdout (tail) ---\n{_tail(proc.stdout)}\n"
            f"--- stderr (tail) ---\n{_tail(proc.stderr)}\n" + hint
        )


def run_antechamber(
    input_file: PathLike,
    output_mol2: PathLike,
    residue_name: str,
    net_charge: int = 0,
    multiplicity: int = 1,
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    extra_args: Sequence[str] = (),
    purge_scratch: bool = True,
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
    input_format: str | None = None,
) -> str:
    """Assign atom types and partial charges with antechamber -> mol2.

    *input_file* may be a PDB or a ligand file with explicit bonds (SDF/MOL2);
    ``input_format`` (antechamber ``-fi``) is inferred from the suffix when not
    given. ``purge_scratch=False`` keeps antechamber's ANTECHAMBER_*/sqm
    scratch files after a successful run (they always survive a failed one),
    which is the way to audit suspicious charges. ``timeout`` is in seconds.
    """
    _check_choice(atom_type, _ATOM_TYPES, "atom_type")
    _check_choice(charge_method, _CHARGE_METHODS, "charge_method")
    exe = _require_executable("antechamber")
    # antechamber runs with cwd set to the output directory (it scatters
    # scratch files there); resolve both paths so relative ones survive.
    input_file = Path(input_file).resolve()
    output_mol2 = Path(output_mol2).resolve()
    if input_format is None:
        suffix = input_file.suffix.lower()
        input_format = _ANTECHAMBER_FORMATS.get(suffix)
        if input_format is None:
            raise ValueError(
                f"Cannot infer the antechamber input format from {input_file.name!r} "
                f"(known suffixes: {sorted(_ANTECHAMBER_FORMATS)}); pass input_format explicitly."
            )
    cmd = [
        exe,
        "-i",
        str(input_file),
        "-fi",
        input_format,
        "-o",
        str(output_mol2),
        "-fo",
        "mol2",
        "-c",
        charge_method,
        "-nc",
        str(net_charge),
        "-m",
        str(multiplicity),
        "-at",
        atom_type,
        "-rn",
        residue_name,
        "-pf",
        "y" if purge_scratch else "n",
        *extra_args,
    ]
    hint = (
        f"For antechamber AM1-BCC failures, inspect 'sqm.out' in {output_mol2.parent}; "
        "the most common causes are missing hydrogens or a wrong net "
        "charge (see the net_charges argument)."
    )
    _run(cmd, cwd=output_mol2.parent, timeout=timeout, hint=hint)
    if not output_mol2.is_file():
        raise RuntimeError(
            f"antechamber exited 0 but did not write {output_mol2}. "
            "acdoctor may have rejected the structure or the input may be "
            "malformed; re-run with purge_scratch=False and inspect the "
            f"scratch files in {output_mol2.parent} "
            "(extra_args=('-dr', 'no') relaxes the acdoctor checks)."
        )
    return str(output_mol2)


def run_parmchk2(
    input_mol2: PathLike,
    output_frcmod: PathLike,
    atom_type: str = "gaff2",
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
) -> str:
    """Generate missing GAFF parameters with parmchk2 -> frcmod."""
    _check_choice(atom_type, _ATOM_TYPES, "atom_type")
    exe = _require_executable("parmchk2")
    # Same cwd game as run_antechamber: resolve so relative paths survive.
    input_mol2 = Path(input_mol2).resolve()
    output_frcmod = Path(output_frcmod).resolve()
    cmd = [
        exe,
        "-i",
        str(input_mol2),
        "-f",
        "mol2",
        "-o",
        str(output_frcmod),
        "-s",
        atom_type,
    ]
    _run(cmd, cwd=output_frcmod.parent, timeout=timeout)
    if not output_frcmod.is_file():
        raise RuntimeError(f"parmchk2 exited 0 but did not write {output_frcmod}.")
    return str(output_frcmod)


def locate_gaff_dat(atom_type: str = "gaff2") -> str:
    """Find gaff.dat / gaff2.dat inside the AmberTools installation."""
    _check_choice(atom_type, _ATOM_TYPES, "atom_type")
    fname = f"{atom_type}.dat"
    candidates: list[Path] = []
    for env in ("AMBERHOME", "CONDA_PREFIX"):
        root = os.environ.get(env)
        if root:
            candidates.append(Path(root) / "dat" / "leap" / "parm" / fname)
    exe = shutil.which("antechamber")
    if exe:
        prefix = Path(exe).resolve().parent.parent
        candidates.append(prefix / "dat" / "leap" / "parm" / fname)
    for cand in candidates:
        if cand.is_file():
            return str(cand)
    raise FileNotFoundError(
        f"Could not locate {fname}. Set the AMBERHOME environment variable "
        f"to your AmberTools installation. Searched: "
        f"{[str(c) for c in candidates] or 'nowhere (no AMBERHOME/CONDA_PREFIX)'}"
    )


# --------------------------------------------------------------------------
# Step 4: ParmEd assembly (mol2 templates + Amber parameters -> OpenMM ffxml)
# --------------------------------------------------------------------------


def _load_residue_template(mol2_file: PathLike, name: str) -> ResidueTemplate:
    template = parmed.load_file(str(mol2_file))
    if isinstance(template, ResidueTemplateContainer):
        if len(template) != 1:
            raise ValueError(f"{mol2_file} contains {len(template)} residues; expected 1.")
        template = template[0]
    if not isinstance(template, ResidueTemplate):
        raise TypeError(f"{mol2_file} did not load as a ResidueTemplate (got {type(template).__name__}).")
    template.name = name
    return template


def assemble_openmm_ffxml(
    templates_mol2: Mapping[str, PathLike],
    parameter_files: Sequence[PathLike],
    output_xml: PathLike,
    write_unused: bool = False,
) -> str:
    """Merge mol2 residue templates and Amber parameter files into a single OpenMM force-field XML.

    ``parameter_files`` typically holds the GAFF database (gaff*.dat) plus the
    per-residue frcmod files. With ``write_unused=False`` only the atom types
    and parameters actually referenced by the residue templates are written,
    which keeps the XML small even though the full GAFF database is loaded.
    """
    params = AmberParameterSet(*[str(f) for f in parameter_files])
    omm_params = OpenMMParameterSet.from_parameterset(params)
    for name, mol2_file in templates_mol2.items():
        omm_params.residues[name] = _load_residue_template(mol2_file, name)

    output_xml = Path(output_xml)
    if output_xml.parent != Path(""):
        output_xml.parent.mkdir(parents=True, exist_ok=True)
    provenance = {"Info": "Generated by forcefill (antechamber/parmchk2 + ParmEd)"}
    omm_params.write(str(output_xml), provenance=provenance, write_unused=write_unused)
    return str(output_xml)


#: Tolerance for comparing numeric force-field attributes when merging. Matches
#: ``openmm.app.forcefield.NonbondedGenerator.SCALETOL``, so a merge that
#: succeeds here is one OpenMM would also accept across separate files.
_SCALE_TOLERANCE = 1e-5


def _attributes_compatible(a: Mapping[str, str], b: Mapping[str, str]) -> bool:
    """True when two force sections can be folded into one element.

    Numeric attributes are compared with a tolerance, because the same constant
    is written to different precision by different producers: Amber's 1-4
    Coulomb scale comes out of ParmEd as ``0.8333333333333334`` and out of
    openmmforcefields as ``0.8333333333``. Anything else must match exactly.
    """
    if set(a) != set(b):
        return False
    for key, left in a.items():
        right = b[key]
        if left == right:
            continue
        try:
            if abs(float(left) - float(right)) > _SCALE_TOLERANCE:
                return False
        except ValueError:
            return False
    return True


def _check_no_redefinition(
    section: ET.Element,
    incoming: ET.Element,
    source: PathLike,
    seen: dict[tuple[str, str], str],
) -> None:
    """Raise if *incoming* redefines an atom type or residue template already merged in."""
    for child in incoming:
        if child.tag not in ("Type", "Residue"):
            continue
        name = child.get("name")
        if name is None:
            continue
        key = (child.tag, name)
        if key not in seen:
            seen[key] = str(source)
            continue
        previous = seen[key]
        existing = next(
            (e for e in section if e.tag == child.tag and e.get("name") == name),
            None,
        )
        if child.tag == "Type" and existing is not None and existing.attrib == child.attrib:
            continue  # Identical redefinition: harmless, keep the one already there.
        fix = (
            "give one of them a different residue name"
            if child.tag == "Residue"
            else "regenerate one of them, since forcefill's backends name their atom types uniquely"
        )
        raise ValueError(
            f"Cannot merge {source}: it defines {child.tag.lower()} {name!r}, "
            f"which {previous} already defines differently. Two force fields "
            f"cannot share a name for different things - {fix}, or load the "
            "files separately with ForceField(...) instead of merging them."
        )


def _element_key(element: ET.Element) -> tuple[str, tuple[tuple[str, str], ...]] | None:
    """Identity of a leaf element for duplicate detection; None for one with children."""
    if len(element):
        return None
    return element.tag, tuple(sorted(element.attrib.items()))


def _extend_section(section: ET.Element, incoming: ET.Element) -> None:
    """Append *incoming*'s children to *section*, dropping exact duplicates of leaf elements.

    Every ffxml this merges declares ``<UseAttributeFromResidue name="charge"/>``
    - charges live on the residue templates, not the atom types - and OpenMM
    rejects a second copy of that declaration outright, because it removes the
    named attribute from the expected list and then cannot find it again. An
    identical leaf element carries no information the first copy did not, so
    dropping it is safe for every section, not just that one.
    """
    present = {key for child in section if (key := _element_key(child)) is not None}
    for child in incoming:
        key = _element_key(child)
        if key is not None:
            if key in present:
                continue
            present.add(key)
        section.append(child)


def merge_ffxml(xml_files: Sequence[PathLike], output_xml: PathLike) -> str:
    """Merge finished OpenMM force-field XML documents into one file.

    The counterpart of :func:`assemble_openmm_ffxml`, which merges at the ParmEd
    parameter-set level and so only works for Amber-style input. This merges the
    XML itself, which is what mixing backends needs: a GAFF ffxml and a SMIRNOFF
    ffxml have nothing in common upstream of the XML.

    Sections of the same kind are folded together when their attributes agree,
    and kept as separate siblings when they do not - which is exactly how OpenMM
    would see them as separate files. That distinction matters: GAFF and
    SMIRNOFF impropers use different ``ordering`` conventions, and OpenMM reads
    ``ordering`` per ``<Improper>`` at parse time, so keeping the two
    ``<PeriodicTorsionForce>`` sections apart is what preserves both.

    Args:
        xml_files: Documents to merge, in order. The first to define a name wins.
        output_xml: Where to write the merged document.

    Returns:
        The path written, as a string.

    Raises:
        ValueError: Two documents define the same atom type or residue template
            differently, or a file is not an OpenMM force-field XML.
    """
    xml_files = list(xml_files)
    if not xml_files:
        raise ValueError("merge_ffxml needs at least one XML file.")

    root = ET.Element("ForceField")
    sections: list[ET.Element] = []
    seen: dict[tuple[str, str], str] = {}
    for xml_file in xml_files:
        source_root = ET.parse(str(xml_file)).getroot()
        if source_root.tag != "ForceField":
            raise ValueError(f"{xml_file} is not an OpenMM force-field XML (root element is <{source_root.tag}>).")
        for incoming in source_root:
            section = next(
                (s for s in sections if s.tag == incoming.tag and _attributes_compatible(s.attrib, incoming.attrib)),
                None,
            )
            if section is None:
                section = ET.SubElement(root, incoming.tag, dict(incoming.attrib))
                sections.append(section)
            _check_no_redefinition(section, incoming, xml_file, seen)
            _extend_section(section, incoming)

    output_xml = Path(output_xml)
    if output_xml.parent != Path(""):
        output_xml.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root)
    ET.ElementTree(root).write(str(output_xml), encoding="utf-8", xml_declaration=True)
    log.info("Merged %d force-field XML files into %s", len(xml_files), output_xml)
    return str(output_xml)


# --------------------------------------------------------------------------
# Step 5: validation and minimization
# --------------------------------------------------------------------------


def _validate_parameterized_residues(
    residues: Mapping[str, app.topology.Residue],
    forcefield: app.ForceField,
    files: Sequence[str],
) -> None:
    """Check that *forcefield* (built from *files*) makes a System for each residue on its own.

    This checks exactly what was produced: that each generated template
    matches the residue's original bond graph and that no parameters are
    missing - independently of whether the rest of the input structure is
    complete. Raises RuntimeError on the first residue that fails. *files*
    only labels the error message.
    """
    files = list(files)
    name = None
    try:
        for name in sorted(residues):
            forcefield.createSystem(_residue_subtopology(residues[name]))
    except Exception as exc:
        raise RuntimeError(
            f"Validation failed: could not build an openmm.System for "
            f"residue {name} on its own from {files}.\n"
            f"OpenMM said: {exc}\n"
            "The generated template does not match the residue's bond "
            "graph, or parameters are missing. If you supplied this residue "
            "via residue_files, its atoms/bonds (including hydrogens) must "
            "match the residue in the PDB exactly."
        ) from exc
    log.info("Validation OK: per-residue Systems built for %s from %s", sorted(residues), files)


def validate_forcefield_xml(
    topology: app.Topology,
    xml_file: PathLike,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    *,
    forcefield: app.ForceField | None = None,
) -> None:
    """Raise RuntimeError unless ``base_forcefield + xml_file`` can build a System for *topology*.

    A pre-built *forcefield* - which must have been constructed from
    ``base_forcefield + xml_file`` - skips re-parsing the XML files.
    """
    files = [*base_forcefield, str(xml_file)]
    try:
        if forcefield is None:
            forcefield = app.ForceField(*files)
        forcefield.createSystem(topology)
    except Exception as exc:
        raise RuntimeError(
            f"Validation failed: could not build an openmm.System from "
            f"{files}.\nOpenMM said: {exc}\n"
            "Common causes: other parts of the structure are missing atoms "
            "or hydrogens (repair with PDBFixer first), or the system "
            "contains ions/waters that need an additional parameter file. "
            "Pass validate=False to skip this check."
        ) from exc
    log.info("Validation OK: System built from %s", files)


class _NonFiniteEnergyError(RuntimeError):
    """Raised by the finite-energy check so it passes through the OpenMM error wrapper unchanged."""


def _describe_topology(topology: app.Topology) -> str:
    """Name *topology* for an error message: a lone residue by name, anything else by size."""
    n_atoms = topology.getNumAtoms()
    n_residues = topology.getNumResidues()
    if n_residues == 1:
        return f"residue {next(topology.residues()).name} ({n_atoms} atoms)"
    return f"the topology ({n_atoms} atoms, {n_residues} residues)"


def minimize_with_forcefield_xml(
    topology: app.Topology,
    positions: unit.Quantity,
    xml_file: PathLike,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    *,
    forcefield: app.ForceField | None = None,
    max_iterations: int = 100,
    tolerance: float = DEFAULT_MINIMIZATION_TOLERANCE,
    nonbonded_method: Any = app.NoCutoff,
    constraints: Any = None,
    rigid_water: bool = False,
    platform_name: str = DEFAULT_MINIMIZATION_PLATFORM,
) -> MinimizationResult:
    """Energy-minimize *topology* with ``base_forcefield + xml_file`` and report what happened.

    This is the check :func:`validate_forcefield_xml` cannot make. Building a
    System only proves the templates match and no parameter is missing; a NaN
    charge, a zero force constant or a broken angle term all survive it and
    surface later as an exploding simulation. Computing an energy and taking a
    few minimizer steps catches them here.

    Raises RuntimeError if the XML will not load, the System will not build, or
    the potential energy is not finite before or after minimizing. The reported
    ``max_force`` and ``energy_change`` are diagnostics, not pass/fail criteria -
    what counts as converged depends on the system.

    Args:
        topology: Structure to minimize.
        positions: Coordinates for *topology*, indexed by global atom index and
            in the same order as ``topology.atoms()``.
        xml_file: The generated force-field XML to load on top of
            *base_forcefield*.
        base_forcefield: ffxml files loaded first, e.g. the standard Amber set.
        forcefield: A pre-built ForceField, which must have been constructed
            from ``base_forcefield + xml_file``; skips re-parsing the XML files.
        max_iterations: Minimizer iteration ceiling. The default is a sanity
            check, not convergence - pass 0 to run until *tolerance* is met.
        tolerance: Convergence target in kJ/mol/nm, applied to the RMS over all
            force *components*, so it is not comparable to the per-atom
            ``max_force`` reported back. Must be positive: OpenMM accepts a
            negative tolerance and then silently minimizes nothing.
        nonbonded_method: ``createSystem`` nonbonded method. The default,
            ``app.NoCutoff``, is exact but O(N^2); pass ``app.PME`` for a
            solvated periodic box.
        constraints: ``createSystem`` constraints. Deliberately ``None`` rather
            than ``app.HBonds``: constraining bonds is exactly what would hide
            a bad bond parameter.
        rigid_water: Deliberately False, unlike OpenMM's default. With
            constraints present ``getForces`` returns the unconstrained forces,
            which makes the reported ``max_force`` a meaningless residual.
        platform_name: OpenMM platform. Pinned to CPU by default so this check
            cannot take a GPU the caller wanted for something else.

    Returns:
        MinimizationResult
    """
    if tolerance <= 0:
        raise ValueError(f"tolerance={tolerance!r} must be positive; OpenMM silently skips minimization otherwise.")
    if max_iterations < 0:
        raise ValueError(f"max_iterations={max_iterations!r} must be >= 0 (0 means 'run until converged').")

    files = [*base_forcefield, str(xml_file)]
    subject = _describe_topology(topology)
    try:
        if forcefield is None:
            forcefield = app.ForceField(*files)
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=nonbonded_method,
            constraints=constraints,
            rigidWater=rigid_water,
        )
        # The integrator is never stepped; a Context just requires one.
        context = Context(system, VerletIntegrator(0.001 * unit.picoseconds), Platform.getPlatformByName(platform_name))
        context.setPositions(positions)

        initial = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        _require_finite_energy(initial, "before minimizing", subject, files)

        LocalEnergyMinimizer.minimize(context, tolerance, max_iterations)

        state = context.getState(getEnergy=True, getForces=True)
        final = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        _require_finite_energy(final, "after minimizing", subject, files)
        forces = state.getForces().value_in_unit(_FORCE_UNIT)
    except _NonFiniteEnergyError:
        raise
    except Exception as exc:
        raise RuntimeError(
            f"Minimization failed: could not minimize {subject} with {files}.\n"
            f"OpenMM said: {exc}\n"
            "Either the generated template does not match the topology, or a "
            "parameter file is missing for some other part of it."
        ) from exc

    result = MinimizationResult(
        n_atoms=topology.getNumAtoms(),
        initial_energy=initial,
        final_energy=final,
        max_force=max((math.sqrt(f[0] ** 2 + f[1] ** 2 + f[2] ** 2) for f in forces), default=0.0),
    )
    log.info(
        "Minimization OK: %s went %.1f -> %.1f kJ/mol (max force %.1f kJ/mol/nm) with %s",
        subject,
        result.initial_energy,
        result.final_energy,
        result.max_force,
        files,
    )
    return result


def _require_finite_energy(energy: float, when: str, subject: str, files: Sequence[str]) -> None:
    """Raise unless *energy* is finite.

    A non-finite energy is the failure this whole check exists to catch: it
    means the parameters are unphysical, which building a System cannot detect.
    """
    if math.isfinite(energy):
        return
    raise _NonFiniteEnergyError(
        f"Minimization failed: the potential energy of {subject} is {energy} "
        f"{when} with {list(files)}.\n"
        "The parameters are not physical. Look for NaN charges or zero-valued "
        "equilibrium bond lengths / angles in the generated XML, and for atoms "
        "sharing coordinates in the input structure."
    )


def _minimize_parameterized_residues(
    residues: Mapping[str, app.topology.Residue],
    positions: unit.Quantity,
    xml_file: PathLike,
    base_forcefield: Sequence[str],
    forcefield: app.ForceField,
    **kwargs: Any,
) -> dict[str, MinimizationResult]:
    """Minimize each residue on its own in vacuum; returns the reports keyed by residue name.

    The counterpart of :func:`_validate_parameterized_residues`: it checks the
    numbers rather than the graph, and for the same reason - each residue is
    tested independently of whether the rest of the input is complete.
    """
    return {
        name: minimize_with_forcefield_xml(
            _residue_subtopology(residues[name]),
            _residue_positions(positions, residues[name]),
            xml_file,
            base_forcefield,
            forcefield=forcefield,
            **kwargs,
        )
        for name in sorted(residues)
    }


# --------------------------------------------------------------------------
# Main entry point
# --------------------------------------------------------------------------


@dataclass
class _ResidueArtifacts:
    """What parameterizing one residue produced.

    ``mol2``/``frcmod`` are None for the smirnoff backend, which has no
    intermediate Amber files - only the finished XML is common to both.
    """

    xml: str
    mol2: str | None = None
    frcmod: str | None = None


def _parameterize_one_residue(
    spec: ResolvedSpec,
    residue: app.topology.Residue | None,
    positions: unit.Quantity | None,
    res_dir: Path,
    *,
    gaff_dat: str | None = None,
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
) -> _ResidueArtifacts:
    """Run one residue through its backend to a self-contained per-residue XML.

    For ``gaff`` that is extract -> antechamber -> parmchk2 -> ParmEd; a
    ``spec.file`` (SDF/MOL2 with explicit bonds) replaces the extraction step.
    For ``smirnoff`` it is a single call into openmmforcefields. *residue* and
    *positions* are None in standalone mode, where there is no input structure
    to extract from.
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    name = spec.name

    if spec.backend == "smirnoff":
        return _ResidueArtifacts(xml=smirnoff.smirnoff_residue_ffxml(spec, res_dir / f"{name}.xml"))

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
    mol2 = run_antechamber(
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
    frcmod = run_parmchk2(mol2, res_dir / f"{name}.frcmod", atom_type=spec.atom_type, timeout=timeout)
    # Per-residue template XML (self-contained).
    gaff_dat = gaff_dat or locate_gaff_dat(spec.atom_type)
    xml = assemble_openmm_ffxml({name: mol2}, [gaff_dat, frcmod], res_dir / f"{name}.xml")
    log.info("Wrote per-residue XML: %s", xml)
    return _ResidueArtifacts(xml=xml, mol2=mol2, frcmod=frcmod)


# --------------------------------------------------------------------------
# Preflight: everything worth knowing before an hour of sqm
# --------------------------------------------------------------------------


def _warn_if_no_hydrogens(name: str, residue: app.topology.Residue) -> None:
    """Warn about a ligand extracted from a PDB that carries no hydrogens at all."""
    n_hydrogens = sum(1 for a in residue.atoms() if a.element is not None and a.element.symbol == "H")
    n_heavy = sum(1 for _ in residue.atoms()) - n_hydrogens
    if n_hydrogens == 0 and n_heavy > 1:
        log.warning(
            "Residue %s contains no hydrogens; AM1-BCC charges will be wrong "
            "unless the molecule really has none. Add explicit hydrogens to the "
            "ligand before parameterizing.",
            name,
        )


def _resolve_smiles(spec: ResolvedSpec, residue: app.topology.Residue | None, positions, res_dir: Path) -> ResolvedSpec:
    """Embed a spec's SMILES to an SDF and return the spec pointing at it.

    When the residue is also present in the input structure, the structure's
    coordinates are kept and only the bond orders come from the SMILES - the
    geometry a crystal structure gives is better than anything embedding
    produces, and the atom count has to match it anyway.
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = res_dir / f"{spec.name}_smiles.sdf"
    if residue is None or positions is None:
        return spec.with_file(ligand_files.smiles_to_sdf(spec.smiles, out_sdf, spec.name))
    residue_pdb = extract_residue_to_pdb(positions, residue, res_dir / f"{spec.name}_extracted.pdb")
    return spec.with_file(ligand_files.smiles_with_residue_geometry(spec.smiles, residue_pdb, out_sdf, spec.name))


def preflight_specs(
    specs: Mapping[str, ResolvedSpec],
    residues: Mapping[str, app.topology.Residue],
    positions: unit.Quantity | None,
    workdir: Path,
    *,
    strict: bool = True,
) -> dict[str, ResolvedSpec]:
    """Check and complete every spec before anything expensive runs.

    Three things go wrong often enough to be worth a dedicated pass, and all
    three are visible in the input:

    * the net charge is left at 0 while the ligand file says otherwise, which
      produces plausible and wrong AM1-BCC charges;
    * the supplied ligand file is not the same molecule as the residue in the
      structure, which only surfaces at the very end as a template mismatch;
    * atoms sit on top of each other, which surfaces as a NaN energy.

    Running this first means a mistake in the last ligand does not cost the
    parameterization of the first. Returns the specs with SMILES resolved to
    files and inferred net charges filled in.

    Args:
        specs: Resolved specs, keyed by residue name.
        residues: The residues as they appear in the input structure; empty in
            standalone mode.
        positions: Coordinates for *residues*, or None in standalone mode.
        workdir: Where a SMILES-derived SDF is written.
        strict: Raise on a composition or geometry fault rather than warn.

    Returns:
        ``{residue_name: ResolvedSpec}``, ready to parameterize.
    """
    out: dict[str, ResolvedSpec] = {}
    for name in sorted(specs):
        spec = specs[name]
        residue = residues.get(name)
        if spec.smiles is not None:
            spec = _resolve_smiles(spec, residue, positions, workdir / name)

        if spec.file is None:
            if residue is not None:
                _warn_if_no_hydrogens(name, residue)
                if positions is not None:
                    ligand_files.check_geometry(
                        [tuple(p.value_in_unit(unit.angstrom)) for p in _residue_positions(positions, residue)],
                        name,
                        strict=strict,
                    )
            out[name] = spec
            continue

        info = ligand_files.inspect_ligand_file(spec.file)
        if residue is not None:
            ligand_files.check_matches_residue(info, residue, name, strict=strict)
        ligand_files.check_geometry(info.positions, name, strict=strict)
        out[name] = _apply_net_charge(spec, info)
    return out


def _apply_net_charge(spec: ResolvedSpec, info: ligand_files.LigandFileInfo) -> ResolvedSpec:
    """Fill in or cross-check the net charge against what the ligand file says.

    A file with real bond orders knows its own formal charge, so leaving
    ``net_charge`` unset is no longer a silent vote for 0. An explicit value that
    contradicts the file is refused outright: one of the two is wrong, and
    guessing which produces exactly the plausible-but-wrong charges this is meant
    to prevent.
    """
    if info.formal_charge is None:
        if spec.net_charge is None:
            log.warning(
                "Could not determine the net charge of %s from %s; assuming 0. "
                "Pass LigandSpec(net_charge=...) if that is wrong - it is the "
                "classic source of plausible but wrong AM1-BCC charges.",
                spec.name,
                Path(info.path).name,
            )
        return spec
    if spec.net_charge is None:
        log.warning(
            "Using net charge %+d for %s, read from %s. Pass LigandSpec(net_charge=...) to override it.",
            info.formal_charge,
            spec.name,
            Path(info.path).name,
        )
        return spec.with_net_charge(info.formal_charge)
    if spec.net_charge != info.formal_charge:
        raise ValueError(
            f"Residue {spec.name} was given net_charge={spec.net_charge:+d}, but "
            f"{Path(info.path).name} describes a molecule with a formal charge "
            f"of {info.formal_charge:+d}. One of them is wrong, and picking "
            "either would give charges that look reasonable and are not. Fix the "
            "protonation state in the ligand file, or correct net_charge."
        )
    return spec


def _combine_residue_xmls(
    artifacts: Mapping[str, _ResidueArtifacts],
    specs: Mapping[str, ResolvedSpec],
    gaff_dat: str | None,
    output_xml: PathLike,
    workdir: Path,
) -> str:
    """Write the one XML covering every parameterized residue.

    All-GAFF goes through ParmEd, which merges at the parameter-set level and
    writes only the atom types the templates actually use. Anything involving
    the smirnoff backend is merged as finished XML instead - the two backends
    share nothing upstream of that.
    """
    gaff = {name for name, spec in specs.items() if spec.backend == "gaff"}
    smirnoff_names = sorted(set(specs) - gaff)
    if not smirnoff_names:
        return assemble_openmm_ffxml(
            {name: artifacts[name].mol2 for name in sorted(gaff)},
            [gaff_dat, *(artifacts[name].frcmod for name in sorted(gaff))],
            output_xml,
        )

    to_merge: list[PathLike] = []
    if gaff:
        # One ParmEd document for all the GAFF residues, so their shared atom
        # types are written once, then merged with the SMIRNOFF ones.
        to_merge.append(
            assemble_openmm_ffxml(
                {name: artifacts[name].mol2 for name in sorted(gaff)},
                [gaff_dat, *(artifacts[name].frcmod for name in sorted(gaff))],
                workdir / "_gaff_combined.xml",
            )
        )
    to_merge += [artifacts[name].xml for name in smirnoff_names]
    return merge_ffxml(to_merge, output_xml)


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
    _check_choice(atom_type, _ATOM_TYPES, "atom_type")
    _check_choice(charge_method, _CHARGE_METHODS, "charge_method")
    _check_choice(backend, BACKENDS, "backend")
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
    gaff_dat: str | None = None
    if any(spec.backend == "gaff" for spec in specs.values()):
        _require_executable("antechamber")
        _require_executable("parmchk2")
        gaff_dat = locate_gaff_dat(atom_type)
        log.info("Using GAFF parameter database: %s", gaff_dat)

    workdir = Path(workdir).resolve() if workdir is not None else Path(tempfile.mkdtemp(prefix="nonstandard_ff_"))
    workdir.mkdir(parents=True, exist_ok=True)
    log.info("Intermediate files in %s", workdir)
    if cleanup and Path(output_xml).resolve().is_relative_to(workdir):
        raise ValueError(
            f"cleanup=True would delete the output XML: {output_xml} resolves "
            f"inside the working directory {workdir}. Write it elsewhere or "
            "pass cleanup=False."
        )

    try:
        # Everything checkable is checked before the first expensive call, so a
        # mistake in the last ligand does not cost the parameterization of the
        # first - antechamber's AM1-BCC can take an hour per ligand.
        specs = preflight_specs(specs, to_param, positions, workdir, strict=strict)

        artifacts: dict[str, _ResidueArtifacts] = {}
        residue_xmls: dict[str, str] = {}
        minimizations: dict[str, MinimizationResult] = {}
        full_minimization: MinimizationResult | None = None
        for name in sorted(specs):
            artifacts[name] = _parameterize_one_residue(
                specs[name],
                to_param[name],
                positions,
                workdir / name,
                gaff_dat=gaff_dat,
                timeout=timeout,
            )
            residue_xmls[name] = artifacts[name].xml

        combined = _combine_residue_xmls(artifacts, specs, gaff_dat, output_xml, workdir)
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
        skipped=skipped,
        workdir=None if cleanup else str(workdir),
        minimizations=minimizations,
        full_minimization=full_minimization,
        cleaning=cleaning,
    )
