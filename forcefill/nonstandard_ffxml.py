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
       topology is built as well.

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
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import parmed
from openmm import app, unit
from parmed.amber import AmberParameterSet
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.openmm import OpenMMParameterSet

log = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_AMBERTOOLS_TIMEOUT",
    "DEFAULT_BASE_FORCEFIELD",
    "ParameterizationResult",
    "assemble_openmm_ffxml",
    "build_forcefield_xml",
    "extract_residue_to_pdb",
    "find_nonstandard_residues",
    "locate_gaff_dat",
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

#: Valid ``atom_type`` values: each names a parameter database ({atom_type}.dat).
_ATOM_TYPES = ("gaff", "gaff2")

#: Charge methods antechamber accepts (see ``antechamber -L``). The list is
#: AmberTools-version-dependent (abcg2 needs >= 23); extend it rather than
#: bypassing the check if your antechamber knows more.
_CHARGE_METHODS = ("bcc", "abcg2", "gas", "mul", "cm1", "cm2", "esp", "resp", "rc", "wc", "dc")

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
) -> None:
    """Warn about net_charges/multiplicities keys that will have no effect.

    A typo'd or case-mismatched key silently leaves the defaults (net
    charge 0, multiplicity 1), which yields plausible but wrong AM1-BCC
    charges - the worst failure mode.
    """
    for label, mapping in (("net_charges", net_charges), ("multiplicities", multiplicities)):
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

    coords = unit.Quantity(
        [positions[a.index].value_in_unit(unit.nanometer) for a in residue.atoms()],
        unit.nanometer,
    )
    out_pdb = Path(out_pdb)
    with open(out_pdb, "w") as fh:
        app.PDBFile.writeFile(sub_top, coords, fh)
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
    input_pdb: PathLike,
    output_mol2: PathLike,
    residue_name: str,
    net_charge: int = 0,
    multiplicity: int = 1,
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    extra_args: Sequence[str] = (),
    purge_scratch: bool = True,
    timeout: float | None = DEFAULT_AMBERTOOLS_TIMEOUT,
) -> str:
    """Assign atom types and partial charges with antechamber -> mol2.

    ``purge_scratch=False`` keeps antechamber's ANTECHAMBER_*/sqm scratch files
    after a successful run (they always survive a failed one), which is the way
    to audit suspicious charges. ``timeout`` is in seconds.
    """
    _check_choice(atom_type, _ATOM_TYPES, "atom_type")
    _check_choice(charge_method, _CHARGE_METHODS, "charge_method")
    exe = _require_executable("antechamber")
    # antechamber runs with cwd set to the output directory (it scatters
    # scratch files there); resolve both paths so relative ones survive.
    input_pdb = Path(input_pdb).resolve()
    output_mol2 = Path(output_mol2).resolve()
    cmd = [
        exe,
        "-i",
        str(input_pdb),
        "-fi",
        "pdb",
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
# Step 5: validation
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
            "graph, or parameters are missing."
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
) -> tuple[str, str, str]:
    """Run one residue through extract -> antechamber -> parmchk2 -> per-residue XML.

    Returns ``(mol2_path, frcmod_path, per_residue_xml_path)``, all inside *res_dir*.
    """
    res_dir.mkdir(parents=True, exist_ok=True)

    n_hydrogens = sum(1 for a in residue.atoms() if a.element is not None and a.element.symbol == "H")
    n_heavy = sum(1 for _ in residue.atoms()) - n_hydrogens
    if n_hydrogens == 0 and n_heavy > 1:
        log.warning(
            "Residue %s contains no hydrogens; AM1-BCC charges will be "
            "wrong unless the molecule really has none. Add explicit "
            "hydrogens to the ligand before parameterizing.",
            name,
        )

    res_pdb = extract_residue_to_pdb(positions, residue, res_dir / f"{name}.pdb")
    log.info("antechamber: %s (net charge %+d, %s/%s)", name, net_charge, atom_type, charge_method)
    mol2 = run_antechamber(
        res_pdb,
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
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    workdir: PathLike | None = None,
    cleanup: bool = False,
    validate: bool = True,
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
    _warn_unused_overrides(to_param, skipped, net_charges, multiplicities)

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
            )

        combined = assemble_openmm_ffxml(mol2_files, [gaff_dat, *frcmod_files.values()], output_xml)
        log.info("Wrote combined force-field XML: %s", combined)

        if validate:
            # Parse the (large) base force field + new XML once; both checks use it.
            files = [*base_forcefield, combined]
            forcefield = app.ForceField(*files)
            _validate_parameterized_residues(to_param, forcefield, files)
            if skipped:
                log.warning(
                    "Skipping full-structure validation: %d residue type(s) "
                    "were skipped (%s) and still have no template, so building "
                    "a System for the whole input cannot succeed. Repair or "
                    "parameterize those, then check with "
                    "validate_forcefield_xml().",
                    len(skipped),
                    ", ".join(sorted(skipped)),
                )
            else:
                validate_forcefield_xml(topology, combined, base_forcefield, forcefield=forcefield)
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
    )
