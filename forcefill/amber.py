"""The AmberTools layer: run antechamber and parmchk2, then assemble their output with ParmEd.

Two halves, in pipeline order.

**The executables.** ``antechamber`` assigns GAFF/GAFF2 atom types and AM1-BCC
charges (-> mol2), ``parmchk2`` fills in whatever GAFF parameters are missing
(-> frcmod). Neither is a Python package, so both must be on ``PATH``:
``conda install -c conda-forge ambertools``. :func:`require_executable` is the
single place that says so, and :func:`locate_gaff_dat` finds the parameter
database inside the same installation.

**The assembly.** :func:`assemble_openmm_ffxml` hands ParmEd the GAFF database,
the frcmods and the mol2 templates and gets one OpenMM ffxml back. That works
only for Amber-style input; merging a finished SMIRNOFF document into the same
file is :mod:`forcefill.merge`'s job.

Nothing here knows about residues, specs or the pipeline - it is the thinnest
wrapper over AmberTools that still gives a useful error message.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path

import parmed
from parmed.amber import AmberParameterSet
from parmed.modeller import ResidueTemplate, ResidueTemplateContainer
from parmed.openmm import OpenMMParameterSet

from ._spec import ATOM_TYPES, CHARGE_METHODS, PathLike, check_choice

log = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_AMBERTOOLS_TIMEOUT",
    "assemble_openmm_ffxml",
    "locate_gaff_dat",
    "run_antechamber",
    "run_parmchk2",
]

#: Ceiling for a single AmberTools invocation, in seconds. sqm's AM1-BCC on a
#: large ligand can legitimately take many minutes; nothing should take an hour.
DEFAULT_AMBERTOOLS_TIMEOUT: float = 3600.0

#: File suffix -> antechamber -fi format, for run_antechamber's inference.
_ANTECHAMBER_FORMATS = {
    ".pdb": "pdb",
    ".mol2": "mol2",
    ".sdf": "sdf",
    ".sd": "sdf",
    ".mol": "mdl",
}


# --------------------------------------------------------------------------
# Running the executables
# --------------------------------------------------------------------------


def require_executable(name: str) -> str:
    """Return the path to *name* on PATH, or raise saying how to install it."""
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
    check_choice(atom_type, ATOM_TYPES, "atom_type")
    check_choice(charge_method, CHARGE_METHODS, "charge_method")
    exe = require_executable("antechamber")
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
    check_choice(atom_type, ATOM_TYPES, "atom_type")
    exe = require_executable("parmchk2")
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
    check_choice(atom_type, ATOM_TYPES, "atom_type")
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
# ParmEd assembly (mol2 templates + Amber parameters -> OpenMM ffxml)
# --------------------------------------------------------------------------


def load_residue_template(mol2_file: PathLike, name: str) -> ResidueTemplate:
    """Load a single-residue mol2 as a ParmEd template, renamed to *name*."""
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
        omm_params.residues[name] = load_residue_template(mol2_file, name)

    output_xml = Path(output_xml)
    if output_xml.parent != Path(""):
        output_xml.parent.mkdir(parents=True, exist_ok=True)
    provenance = {"Info": "Generated by forcefill (antechamber/parmchk2 + ParmEd)"}
    omm_params.write(str(output_xml), provenance=provenance, write_unused=write_unused)
    return str(output_xml)
