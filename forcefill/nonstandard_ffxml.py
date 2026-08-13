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

    3. Each unique parameterizable residue is written to its own PDB and run
       through AmberTools:

       * ``antechamber``  assigns GAFF/GAFF2 atom types + AM1-BCC charges
         (-> ``<RES>.mol2`` residue template)
       * ``parmchk2``     generates any GAFF parameters that are missing
         (-> ``<RES>.frcmod``)

    4. ParmEd merges the GAFF parameter database, the frcmod files and the mol2
       templates into OpenMM ffxml: one XML per residue plus one combined XML
       containing every residue (atom types, residue templates with charges,
       bonds/angles/torsions and nonbonded parameters).

    5. The combined XML is validated by building an ``openmm.System`` from
       ``base force field + new XML`` for each parameterized residue on its
       own; when no residues were skipped, a System for the full input
       topology is built as well. With ``minimize=True`` each of those is also
       energy-minimized, which catches unphysical parameters that still build
       a perfectly valid System.

Requirements:
    * ``openmm >= 7.6`` and ``parmed >= 3.4`` (``pip install openmm parmed``)
    * AmberTools (``antechamber``, ``parmchk2``) on ``PATH`` for step 3, e.g.
      ``conda install -c conda-forge ambertools``

Example:
    >>> from forcefill import build_forcefield_xml
    >>> result = build_forcefield_xml("complex.pdb", "extras.xml", net_charges={"LIG": -1})
    >>> result.parameterized
    ['LIG']

    then simulate with::

        ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
        system = ff.createSystem(pdb.topology, ...)

Notes:
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
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import parmed
from openmm import Context, LocalEnergyMinimizer, Platform, VerletIntegrator, app, unit
from parmed.amber import AmberParameterSet
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.openmm import OpenMMParameterSet

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

#: Valid ``atom_type`` values: each names a parameter database ({atom_type}.dat).
_ATOM_TYPES = ("gaff", "gaff2")

#: Charge methods antechamber accepts (see ``antechamber -L``). The list is
#: AmberTools-version-dependent (abcg2 needs >= 23); extend it rather than
#: bypassing the check if your antechamber knows more.
_CHARGE_METHODS = ("bcc", "abcg2", "gas", "mul", "cm1", "cm2", "esp", "resp", "rc", "wc", "dc")

#: File suffix -> antechamber -fi format, for run_antechamber's inference.
_ANTECHAMBER_FORMATS = {
    ".pdb": "pdb",
    ".mol2": "mol2",
    ".sdf": "sdf",
    ".sd": "sdf",
    ".mol": "mdl",
}

#: Residue names the base force fields already know. Unmatched residues with
#: these names are almost always incomplete structures, not new chemistry.
_STANDARD_RESIDUES = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ASH",
        "GLH",
        "LYN",
        "CYX",
        "CYM",
        "HID",
        "HIE",
        "HIP",
        "ACE",
        "NME",
        "NMA",
        "DA",
        "DC",
        "DG",
        "DT",
        "DA3",
        "DA5",
        "DC3",
        "DC5",
        "DG3",
        "DG5",
        "DT3",
        "DT5",
        "A",
        "C",
        "G",
        "U",
        "A3",
        "A5",
        "C3",
        "C5",
        "G3",
        "G5",
        "U3",
        "U5",
        "HOH",
        "WAT",
        "H2O",
        "TIP",
        "TIP3",
        "TP3",
        "SPC",
        "SOL",
    ]
)


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
    """What :func:`build_forcefield_xml` produced."""

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
    #: Minimization of the whole input topology. ``None`` unless
    #: ``minimize=True`` *and* no residue was skipped.
    full_minimization: MinimizationResult | None = None


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
        if name in _STANDARD_RESIDUES:
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
) -> None:
    """Warn about net_charges/multiplicities/residue_files keys with no effect.

    A typo'd or case-mismatched key silently leaves the defaults (net
    charge 0, multiplicity 1, PDB extraction), which yields plausible but
    wrong AM1-BCC charges - the worst failure mode.
    """
    for label, mapping in (
        ("net_charges", net_charges),
        ("multiplicities", multiplicities),
        ("residue_files", residue_files or {}),
    ):
        for key in mapping:
            if key in to_param:
                continue
            if key in skipped:
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


def _check_choice(value: str, valid: Sequence[str], label: str) -> None:
    """Raise ValueError early for a typo'd option instead of a late, cryptic AmberTools failure."""
    if value not in valid:
        raise ValueError(f"{label}={value!r} is not one of {list(valid)}")


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


def _parameterize_one_residue(
    name: str,
    residue: app.topology.Residue,
    positions: unit.Quantity,
    res_dir: Path,
    *,
    gaff_dat: str,
    net_charge: int,
    multiplicity: int,
    atom_type: str,
    charge_method: str,
    antechamber_args: Sequence[str],
    timeout: float | None,
    input_file: PathLike | None = None,
) -> tuple[str, str, str]:
    """Run one residue through extract -> antechamber -> parmchk2 -> per-residue XML.

    A user-supplied *input_file* (SDF/MOL2 with explicit bonds) replaces the
    PDB extraction step. Returns ``(mol2_path, frcmod_path,
    per_residue_xml_path)``, all inside *res_dir*.
    """
    res_dir.mkdir(parents=True, exist_ok=True)

    if input_file is None:
        n_hydrogens = sum(1 for a in residue.atoms() if a.element is not None and a.element.symbol == "H")
        n_heavy = sum(1 for _ in residue.atoms()) - n_hydrogens
        if n_hydrogens == 0 and n_heavy > 1:
            log.warning(
                "Residue %s contains no hydrogens; AM1-BCC charges will be "
                "wrong unless the molecule really has none. Add explicit "
                "hydrogens to the ligand before parameterizing.",
                name,
            )
        antechamber_input: PathLike = extract_residue_to_pdb(positions, residue, res_dir / f"{name}.pdb")
    else:
        log.info("Using the supplied ligand file for %s: %s", name, input_file)
        antechamber_input = input_file

    log.info("antechamber: %s (net charge %+d, %s/%s)", name, net_charge, atom_type, charge_method)
    mol2 = run_antechamber(
        antechamber_input,
        res_dir / f"{name}.mol2",
        name,
        net_charge=net_charge,
        multiplicity=multiplicity,
        atom_type=atom_type,
        charge_method=charge_method,
        extra_args=antechamber_args,
        timeout=timeout,
    )
    frcmod = run_parmchk2(mol2, res_dir / f"{name}.frcmod", atom_type=atom_type, timeout=timeout)
    # Per-residue template XML (self-contained).
    xml = assemble_openmm_ffxml({name: mol2}, [gaff_dat, frcmod], res_dir / f"{name}.xml")
    log.info("Wrote per-residue XML: %s", xml)
    return mol2, frcmod, xml


def build_forcefield_xml(
    pdb_file: PathLike,
    output_xml: PathLike = "nonstandard_ff.xml",
    *,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    net_charges: Mapping[str, int] | None = None,
    multiplicities: Mapping[str, int] | None = None,
    residue_files: Mapping[str, PathLike] | None = None,
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    workdir: PathLike | None = None,
    cleanup: bool = False,
    validate: bool = True,
    minimize: bool = False,
    antechamber_args: Sequence[str] = (),
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
) -> ParameterizationResult:
    """Identify non-standard residues in *pdb_file* and build an OpenMM force-field XML for them.

    Args:
        pdb_file: Input structure. Ligands must have explicit hydrogens; the
            element columns and (for hetero groups) CONECT records should be
            present.
        output_xml: Where to write the combined force-field XML.
        base_forcefield: ffxml files defining what counts as "standard".
        net_charges: ``{residue_name: net_charge}``; defaults to 0 per
            residue. Getting this right is essential for sensible AM1-BCC
            charges.
        multiplicities: ``{residue_name: spin_multiplicity}``; defaults to 1.
        residue_files: ``{residue_name: ligand_file}`` - parameterize these
            residues from the given SDF/MOL2 (or other antechamber-readable)
            file instead of extracting them from the PDB. A file with explicit
            bond orders and protonation as drawn avoids antechamber having to
            re-perceive them from PDB geometry, a classic source of silently
            wrong atom types. The file must contain the *same* atoms and bonds
            (including hydrogens) as the residue in the PDB, or validation
            fails: the template is still matched against the PDB's bond graph.
        atom_type: ``"gaff2"`` (default) or ``"gaff"``.
        charge_method: antechamber charge method, default ``"bcc"`` (AM1-BCC).
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
            full input topology; with skipped residues present that check
            would always fail (they still have no template), so it is logged
            and omitted instead.
        minimize: If True, additionally energy-minimize each parameterized
            residue in vacuum (and the full input topology, under the same
            no-residues-skipped condition as *validate*), reported back in
            ``minimizations`` and ``full_minimization``. This catches
            unphysical parameters - a NaN charge, a zero force constant -
            which building a System alone cannot. Off by default because it
            costs an energy evaluation per residue. It subsumes *validate*:
            minimizing implies building the System, just with a less specific
            error message when a template does not match.
        antechamber_args: Extra raw arguments appended to the antechamber
            command line (e.g. ``("-dr", "no")`` to relax acdoctor structure
            checks).
        timeout: Per-invocation ceiling in seconds for each antechamber /
            parmchk2 run (None disables it). Default one hour.

    Returns:
        ParameterizationResult
    """
    _check_choice(atom_type, _ATOM_TYPES, "atom_type")
    _check_choice(charge_method, _CHARGE_METHODS, "charge_method")
    net_charges = dict(net_charges or {})
    multiplicities = dict(multiplicities or {})
    residue_files = dict(residue_files or {})

    pdb = app.PDBFile(str(pdb_file))
    topology, positions = pdb.topology, pdb.positions

    unmatched = find_nonstandard_residues(topology, base_forcefield)
    if not unmatched:
        log.info("All residues matched %s - nothing to parameterize.", list(base_forcefield))
        return ParameterizationResult(forcefield_xml=None)

    to_param, skipped = _classify_unmatched(unmatched)
    for name, reason in skipped.items():
        log.warning("Skipping %s: %s", name, reason)
    if not to_param:
        details = "\n".join(f"  {n}: {r}" for n, r in skipped.items())
        raise RuntimeError(f"Unmatched residues were found but none can be auto-parameterized:\n{details}")
    log.info("Residues to parameterize: %s", sorted(to_param))
    _warn_unused_overrides(to_param, skipped, net_charges, multiplicities, residue_files)

    # Fail early if AmberTools is absent.
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
        mol2_files: dict[str, str] = {}
        frcmod_files: dict[str, str] = {}
        residue_xmls: dict[str, str] = {}
        minimizations: dict[str, MinimizationResult] = {}
        full_minimization: MinimizationResult | None = None
        for name in sorted(to_param):
            mol2_files[name], frcmod_files[name], residue_xmls[name] = _parameterize_one_residue(
                name,
                to_param[name],
                positions,
                workdir / name,
                gaff_dat=gaff_dat,
                net_charge=net_charges.get(name, 0),
                multiplicity=multiplicities.get(name, 1),
                atom_type=atom_type,
                charge_method=charge_method,
                antechamber_args=antechamber_args,
                timeout=timeout,
                input_file=residue_files.get(name),
            )

        combined = assemble_openmm_ffxml(mol2_files, [gaff_dat, *frcmod_files.values()], output_xml)
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
        parameterized=sorted(mol2_files),
        skipped=skipped,
        workdir=None if cleanup else str(workdir),
        minimizations=minimizations,
        full_minimization=full_minimization,
    )
