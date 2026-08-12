#!/usr/bin/env python3
"""
nonstandard_ffxml.py
====================

Identify non-standard residues in a PDB file with OpenMM and produce a
ready-to-use OpenMM force-field XML for them.

Pipeline
--------
1.  ``ForceField.getUnmatchedResidues`` finds every residue the chosen base
    force field cannot match.
2.  Unmatched residues are classified:

    * standard residues that are merely missing atoms   -> reported, skipped
      (repair the structure with PDBFixer / ``Modeller.addHydrogens`` instead)
    * monatomic species (ions)                          -> reported, skipped
      (load an ion parameter file; never run antechamber on a bare ion)
    * residues covalently bonded to neighbours          -> skipped by default
      (stand-alone GAFF treatment is wrong for polymer-linked residues such
      as modified amino acids; those need capping + a consistent charge
      derivation, e.g. with pyRED or ffparam-style workflows)
    * free-standing hetero molecules (ligands, cofactors) -> parameterized

3.  Each unique parameterizable residue is written to its own PDB and run
    through AmberTools:

    * ``antechamber``  assigns GAFF/GAFF2 atom types + AM1-BCC charges
      (-> ``<RES>.mol2`` residue template)
    * ``parmchk2``     generates any GAFF parameters that are missing
      (-> ``<RES>.frcmod``)

4.  ParmEd merges the GAFF parameter database, the frcmod files and the mol2
    templates into OpenMM ffxml: one XML per residue plus one combined XML
    containing every residue (atom types, residue templates with charges,
    bonds/angles/torsions and nonbonded parameters).

5.  The combined XML is validated by building an ``openmm.System`` from
    ``base force field + new XML``.

Requirements
------------
* ``openmm >= 7.6`` and ``parmed >= 3.4`` (``pip install openmm parmed``)
* AmberTools (``antechamber``, ``parmchk2``) on ``PATH`` for step 3, e.g.
  ``conda install -c conda-forge ambertools``

Example
-------
>>> from nonstandard_ffxml import build_forcefield_xml
>>> result = build_forcefield_xml("complex.pdb", "extras.xml",
...                               net_charges={"LIG": -1})
>>> result.parameterized
['LIG']

then simulate with::

    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", "extras.xml")
    system = ff.createSystem(pdb.topology, ...)

Command line::

    python nonstandard_ffxml.py complex.pdb -o extras.xml --charge LIG=-1

Notes
-----
* Ligands must contain **all explicit hydrogens** with reasonable geometry;
  AM1-BCC charges are meaningless otherwise.
* Load either the combined XML *or* the per-residue XMLs with a ForceField,
  never both at once (duplicate GAFF atom-type definitions would collide).
* If you would rather not manage XML files at all, the same job can be done
  at runtime with ``openmmforcefields.generators.GAFFTemplateGenerator``.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Union

import parmed
from parmed.amber import AmberParameterSet
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.openmm import OpenMMParameterSet

from openmm import app, unit

log = logging.getLogger("nonstandard_ffxml")

PathLike = Union[str, os.PathLike]

#: Base force field used to decide what counts as "non-standard".
DEFAULT_BASE_FORCEFIELD = ("amber14-all.xml", "amber14/tip3p.xml")

#: Residue names the base force fields already know. Unmatched residues with
#: these names are almost always incomplete structures, not new chemistry.
_STANDARD_RESIDUES = frozenset(
    """
    ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP
    TYR VAL ASH GLH LYN CYX CYM HID HIE HIP ACE NME NMA
    DA DC DG DT DA3 DA5 DC3 DC5 DG3 DG5 DT3 DT5
    A C G U A3 A5 C3 C5 G3 G5 U3 U5
    HOH WAT H2O TIP TIP3 TP3 SPC SOL
    """.split()
)


@dataclass
class ParameterizationResult:
    """What :func:`build_forcefield_xml` produced."""

    #: Path to the combined ffxml covering every parameterized residue
    #: (``None`` when nothing needed parameterizing).
    forcefield_xml: Optional[str]
    #: Per-residue ffxml files, keyed by residue name.
    residue_xmls: Dict[str, str] = field(default_factory=dict)
    #: Residue names that were successfully parameterized.
    parameterized: List[str] = field(default_factory=list)
    #: Residue names that were skipped, mapped to the reason.
    skipped: Dict[str, str] = field(default_factory=dict)
    #: Directory holding intermediate files (per-residue PDB/mol2/frcmod).
    workdir: Optional[str] = None


# --------------------------------------------------------------------------
# Step 1: identification
# --------------------------------------------------------------------------

def find_nonstandard_residues(
    topology: app.Topology,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
) -> List[app.topology.Residue]:
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
) -> tuple[Dict[str, app.topology.Residue], Dict[str, str]]:
    """Split unmatched residues into {name: representative} to parameterize
    and {name: reason} to skip."""
    groups: Dict[str, List[app.topology.Residue]] = defaultdict(list)
    for res in unmatched:
        groups[res.name.strip()].append(res)

    to_param: Dict[str, app.topology.Residue] = {}
    skipped: Dict[str, str] = {}
    for name, residues in groups.items():
        counts = sorted({sum(1 for _ in r.atoms()) for r in residues})
        rep = max(residues, key=lambda r: sum(1 for _ in r.atoms()))
        n_atoms = sum(1 for _ in rep.atoms())
        if len(counts) > 1:
            log.warning(
                "Copies of residue %s differ in atom count (%s); "
                "using the most complete copy (%d atoms) as the template.",
                name, counts, n_atoms,
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
                "monatomic species - use an ion parameter file for it "
                "(GAFF/antechamber cannot treat bare ions)"
            )
        elif any(True for _ in rep.external_bonds()):
            skipped[name] = (
                "covalently bonded to neighbouring residues - a stand-alone "
                "GAFF parameterization is not valid for polymer-linked "
                "residues; cap the fragment and derive charges consistently "
                "with the backbone force field instead"
            )
        else:
            to_param[name] = rep
    return to_param, skipped


# --------------------------------------------------------------------------
# Step 2: extraction
# --------------------------------------------------------------------------

def extract_residue_to_pdb(
    topology: app.Topology,
    positions: unit.Quantity,
    residue: app.topology.Residue,
    out_pdb: PathLike,
) -> str:
    """Write a single residue (its atoms, internal bonds and coordinates)
    to *out_pdb* and return the path."""
    sub_top = app.Topology()
    chain = sub_top.addChain("A")
    new_res = sub_top.addResidue(residue.name, chain)
    atom_map = {}
    for atom in residue.atoms():
        if atom.element is None:
            log.warning(
                "Atom %s in residue %s has no element assigned; antechamber "
                "may misread it. Check the element columns of the PDB.",
                atom.name, residue.name,
            )
        atom_map[atom] = sub_top.addAtom(atom.name, atom.element, new_res)
    for bond in residue.internal_bonds():
        sub_top.addBond(atom_map[bond.atom1], atom_map[bond.atom2])

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

def _require_executable(name: str) -> str:
    exe = shutil.which(name)
    if exe is None:
        raise RuntimeError(
            f"'{name}' was not found on PATH. It is part of AmberTools; "
            "install it with e.g. 'conda install -c conda-forge ambertools' "
            "and try again."
        )
    return exe


def _run(cmd: Sequence[str], cwd: PathLike) -> None:
    log.debug("Running: %s (cwd=%s)", " ".join(map(str, cmd)), cwd)
    proc = subprocess.run(
        [str(c) for c in cmd], cwd=str(cwd), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}:\n"
            f"  {' '.join(map(str, cmd))}\n"
            f"--- stdout (tail) ---\n{proc.stdout[-2000:]}\n"
            f"--- stderr (tail) ---\n{proc.stderr[-2000:]}\n"
            f"For antechamber AM1-BCC failures, inspect 'sqm.out' in {cwd}; "
            "the most common causes are missing hydrogens or a wrong net "
            "charge (see the net_charges argument)."
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
) -> str:
    """Assign atom types and partial charges with antechamber -> mol2."""
    exe = _require_executable("antechamber")
    cmd = [
        exe,
        "-i", str(input_pdb), "-fi", "pdb",
        "-o", str(output_mol2), "-fo", "mol2",
        "-c", charge_method,
        "-nc", str(net_charge),
        "-m", str(multiplicity),
        "-at", atom_type,
        "-rn", residue_name,
        "-pf", "y",
        *extra_args,
    ]
    _run(cmd, cwd=Path(output_mol2).parent)
    return str(output_mol2)


def run_parmchk2(
    input_mol2: PathLike,
    output_frcmod: PathLike,
    atom_type: str = "gaff2",
) -> str:
    """Generate missing GAFF parameters with parmchk2 -> frcmod."""
    exe = _require_executable("parmchk2")
    cmd = [
        exe,
        "-i", str(input_mol2), "-f", "mol2",
        "-o", str(output_frcmod),
        "-s", atom_type,
    ]
    _run(cmd, cwd=Path(output_frcmod).parent)
    return str(output_frcmod)


def locate_gaff_dat(atom_type: str = "gaff2") -> str:
    """Find gaff.dat / gaff2.dat inside the AmberTools installation."""
    fname = f"{atom_type}.dat"
    candidates: List[Path] = []
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

def _load_amber_parameters(parameter_files: Sequence[PathLike]) -> AmberParameterSet:
    files = [str(f) for f in parameter_files]
    try:
        return AmberParameterSet(*files)
    except TypeError:  # very old ParmEd: load one at a time
        params = AmberParameterSet()
        for f in files:
            params.load_parameters(f)
        return params


def _load_residue_template(mol2_file: PathLike, name: str) -> ResidueTemplate:
    template = parmed.load_file(str(mol2_file))
    if isinstance(template, ResidueTemplateContainer):
        if len(template) != 1:
            raise ValueError(
                f"{mol2_file} contains {len(template)} residues; expected 1."
            )
        template = template[0]
    if not isinstance(template, ResidueTemplate):
        raise TypeError(
            f"{mol2_file} did not load as a ResidueTemplate "
            f"(got {type(template).__name__})."
        )
    template.name = name
    return template


def assemble_openmm_ffxml(
    templates_mol2: Mapping[str, PathLike],
    parameter_files: Sequence[PathLike],
    output_xml: PathLike,
    write_unused: bool = False,
) -> str:
    """Merge mol2 residue templates and Amber parameter files (gaff*.dat,
    frcmod) into a single OpenMM force-field XML.

    With ``write_unused=False`` only the atom types and parameters actually
    referenced by the residue templates are written, which keeps the XML
    small even though the full GAFF database is loaded.
    """
    params = _load_amber_parameters(parameter_files)
    omm_params = OpenMMParameterSet.from_parameterset(params)
    for name, mol2_file in templates_mol2.items():
        omm_params.residues[name] = _load_residue_template(mol2_file, name)

    output_xml = Path(output_xml)
    if output_xml.parent != Path(""):
        output_xml.parent.mkdir(parents=True, exist_ok=True)
    provenance = {
        "Info": "Generated by nonstandard_ffxml.py "
                "(antechamber/parmchk2 + ParmEd)"
    }
    try:
        omm_params.write(str(output_xml), provenance=provenance,
                         write_unused=write_unused)
    except TypeError:  # older ParmEd without one of the keywords
        omm_params.write(str(output_xml))
    return str(output_xml)


# --------------------------------------------------------------------------
# Step 5: validation
# --------------------------------------------------------------------------

def validate_forcefield_xml(
    topology: app.Topology,
    xml_file: PathLike,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
) -> None:
    """Raise RuntimeError unless base_forcefield + xml_file can build a
    System for *topology*."""
    files = [*base_forcefield, str(xml_file)]
    try:
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

def build_forcefield_xml(
    pdb_file: PathLike,
    output_xml: PathLike = "nonstandard_ff.xml",
    *,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    net_charges: Optional[Mapping[str, int]] = None,
    multiplicities: Optional[Mapping[str, int]] = None,
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    workdir: Optional[PathLike] = None,
    validate: bool = True,
    antechamber_args: Sequence[str] = (),
) -> ParameterizationResult:
    """Identify non-standard residues in *pdb_file* and build an OpenMM
    force-field XML for them.

    Parameters
    ----------
    pdb_file
        Input structure. Ligands must have explicit hydrogens; the element
        columns and (for hetero groups) CONECT records should be present.
    output_xml
        Where to write the combined force-field XML.
    base_forcefield
        ffxml files defining what counts as "standard".
    net_charges
        ``{residue_name: net_charge}``; defaults to 0 per residue. Getting
        this right is essential for sensible AM1-BCC charges.
    multiplicities
        ``{residue_name: spin_multiplicity}``; defaults to 1.
    atom_type
        ``"gaff2"`` (default) or ``"gaff"``.
    charge_method
        antechamber charge method, default ``"bcc"`` (AM1-BCC).
    workdir
        Directory for intermediate files (per-residue PDB, mol2, frcmod,
        per-residue XML). A fresh temporary directory is created and kept
        if not given; its path is reported in the result.
    validate
        If True (default), verify that ``base_forcefield + output_xml`` can
        build a System for the full input topology.
    antechamber_args
        Extra raw arguments appended to the antechamber command line
        (e.g. ``("-dr", "no")`` to relax acdoctor structure checks).

    Returns
    -------
    ParameterizationResult
    """
    net_charges = dict(net_charges or {})
    multiplicities = dict(multiplicities or {})

    pdb = app.PDBFile(str(pdb_file))
    topology, positions = pdb.topology, pdb.positions

    unmatched = find_nonstandard_residues(topology, base_forcefield)
    if not unmatched:
        log.info("All residues matched %s - nothing to parameterize.",
                 list(base_forcefield))
        return ParameterizationResult(forcefield_xml=None)

    to_param, skipped = _classify_unmatched(unmatched)
    for name, reason in skipped.items():
        log.warning("Skipping %s: %s", name, reason)
    if not to_param:
        details = "\n".join(f"  {n}: {r}" for n, r in skipped.items())
        raise RuntimeError(
            "Unmatched residues were found but none can be "
            f"auto-parameterized:\n{details}"
        )
    log.info("Residues to parameterize: %s", sorted(to_param))

    # Fail early if AmberTools is absent.
    _require_executable("antechamber")
    _require_executable("parmchk2")
    gaff_dat = locate_gaff_dat(atom_type)
    log.info("Using GAFF parameter database: %s", gaff_dat)

    workdir = Path(workdir) if workdir is not None else Path(
        tempfile.mkdtemp(prefix="nonstandard_ff_")
    )
    workdir.mkdir(parents=True, exist_ok=True)
    log.info("Intermediate files in %s", workdir)

    mol2_files: Dict[str, str] = {}
    frcmod_files: Dict[str, str] = {}
    residue_xmls: Dict[str, str] = {}
    for name in sorted(to_param):
        residue = to_param[name]
        res_dir = workdir / name
        res_dir.mkdir(exist_ok=True)

        n_hydrogens = sum(
            1 for a in residue.atoms()
            if a.element is not None and a.element.symbol == "H"
        )
        n_heavy = sum(1 for _ in residue.atoms()) - n_hydrogens
        if n_hydrogens == 0 and n_heavy > 1:
            log.warning(
                "Residue %s contains no hydrogens; AM1-BCC charges will be "
                "wrong unless the molecule really has none. Add explicit "
                "hydrogens to the ligand before parameterizing.", name,
            )

        res_pdb = extract_residue_to_pdb(
            topology, positions, residue, res_dir / f"{name}.pdb"
        )
        log.info("antechamber: %s (net charge %+d, %s/%s)",
                 name, net_charges.get(name, 0), atom_type, charge_method)
        mol2_files[name] = run_antechamber(
            res_pdb, res_dir / f"{name}.mol2", name,
            net_charge=net_charges.get(name, 0),
            multiplicity=multiplicities.get(name, 1),
            atom_type=atom_type, charge_method=charge_method,
            extra_args=antechamber_args,
        )
        frcmod_files[name] = run_parmchk2(
            mol2_files[name], res_dir / f"{name}.frcmod", atom_type=atom_type
        )
        # Per-residue template XML (self-contained).
        residue_xmls[name] = assemble_openmm_ffxml(
            {name: mol2_files[name]},
            [gaff_dat, frcmod_files[name]],
            res_dir / f"{name}.xml",
        )
        log.info("Wrote per-residue XML: %s", residue_xmls[name])

    combined = assemble_openmm_ffxml(
        mol2_files, [gaff_dat, *frcmod_files.values()], output_xml
    )
    log.info("Wrote combined force-field XML: %s", combined)

    if validate:
        validate_forcefield_xml(topology, combined, base_forcefield)

    return ParameterizationResult(
        forcefield_xml=combined,
        residue_xmls=residue_xmls,
        parameterized=sorted(mol2_files),
        skipped=skipped,
        workdir=str(workdir),
    )


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def _parse_charge_specs(specs: Sequence[str]) -> Dict[str, int]:
    charges: Dict[str, int] = {}
    for spec in specs:
        try:
            name, value = spec.split("=")
            charges[name.strip()] = int(value)
        except ValueError as exc:
            raise SystemExit(
                f"Bad --charge specification '{spec}' (expected RES=INT)"
            ) from exc
    return charges


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Parameterize non-standard residues in a PDB file and "
                    "write an OpenMM force-field XML for them.",
    )
    parser.add_argument("pdb", help="input PDB file")
    parser.add_argument("-o", "--output", default="nonstandard_ff.xml",
                        help="combined output ffxml (default: %(default)s)")
    parser.add_argument("--charge", action="append", default=[],
                        metavar="RES=INT",
                        help="net charge of a residue, repeatable "
                             "(e.g. --charge LIG=-1)")
    parser.add_argument("--atom-type", choices=["gaff", "gaff2"],
                        default="gaff2")
    parser.add_argument("--charge-method", default="bcc",
                        help="antechamber charge method (default: bcc)")
    parser.add_argument("--base-ff", nargs="+",
                        default=list(DEFAULT_BASE_FORCEFIELD),
                        help="base ffxml files (default: %(default)s)")
    parser.add_argument("--workdir", default=None,
                        help="directory for intermediate files")
    parser.add_argument("--no-validate", action="store_true",
                        help="skip building a System as a final check")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    result = build_forcefield_xml(
        args.pdb, args.output,
        base_forcefield=tuple(args.base_ff),
        net_charges=_parse_charge_specs(args.charge),
        atom_type=args.atom_type,
        charge_method=args.charge_method,
        workdir=args.workdir,
        validate=not args.no_validate,
    )
    print(f"\nParameterized residues : {', '.join(result.parameterized) or '-'}")
    for name, reason in result.skipped.items():
        print(f"Skipped {name:<12}: {reason}")
    print(f"Combined force field   : {result.forcefield_xml}")
    print(f"Intermediates          : {result.workdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
