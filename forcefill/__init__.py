from importlib import metadata as _metadata

from .nonstandard_ffxml import (
    DEFAULT_AMBERTOOLS_TIMEOUT,
    DEFAULT_BASE_FORCEFIELD,
    ParameterizationResult,
    assemble_openmm_ffxml,
    build_forcefield_xml,
    extract_residue_to_pdb,
    find_nonstandard_residues,
    locate_gaff_dat,
    run_antechamber,
    run_parmchk2,
    validate_forcefield_xml,
)

# Single source of truth for the version is pyproject.toml; installed metadata
# carries it. The fallback covers running from an uninstalled checkout.
try:
    __version__ = _metadata.version("forcefill")
except _metadata.PackageNotFoundError:  # uninstalled checkout
    __version__ = "0.0.0+unknown"

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
