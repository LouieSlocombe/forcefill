from importlib import metadata as _metadata

from .clean_structure import (
    ADDITIVE_RESIDUES,
    BULK_ION_RESIDUES,
    STRUCTURAL_METAL_RESIDUES,
    WATER_RESIDUES,
    CleaningResult,
    clean_pdb,
    clean_topology,
)
from .nonstandard_ffxml import (
    DEFAULT_AMBERTOOLS_TIMEOUT,
    DEFAULT_BASE_FORCEFIELD,
    DEFAULT_MINIMIZATION_PLATFORM,
    DEFAULT_MINIMIZATION_TOLERANCE,
    MinimizationResult,
    ParameterizationResult,
    assemble_openmm_ffxml,
    build_forcefield_xml,
    extract_residue_to_pdb,
    find_nonstandard_residues,
    locate_gaff_dat,
    minimize_with_forcefield_xml,
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
    "ADDITIVE_RESIDUES",
    "BULK_ION_RESIDUES",
    "DEFAULT_AMBERTOOLS_TIMEOUT",
    "DEFAULT_BASE_FORCEFIELD",
    "DEFAULT_MINIMIZATION_PLATFORM",
    "DEFAULT_MINIMIZATION_TOLERANCE",
    "STRUCTURAL_METAL_RESIDUES",
    "WATER_RESIDUES",
    "CleaningResult",
    "MinimizationResult",
    "ParameterizationResult",
    "assemble_openmm_ffxml",
    "build_forcefield_xml",
    "clean_pdb",
    "clean_topology",
    "extract_residue_to_pdb",
    "find_nonstandard_residues",
    "locate_gaff_dat",
    "minimize_with_forcefield_xml",
    "run_antechamber",
    "run_parmchk2",
    "validate_forcefield_xml",
]
