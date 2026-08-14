from importlib import metadata as _metadata

from ._pipeline import ParameterizationResult
from ._spec import (
    BACKENDS,
    DEFAULT_BASE_FORCEFIELD,
    DEFAULT_SMIRNOFF_FORCEFIELD,
    LigandSpec,
)
from .amber import (
    DEFAULT_AMBERTOOLS_TIMEOUT,
    assemble_openmm_ffxml,
    locate_gaff_dat,
    run_antechamber,
    run_parmchk2,
)
from .checks import (
    DEFAULT_MINIMIZATION_PLATFORM,
    DEFAULT_MINIMIZATION_TOLERANCE,
    MinimizationResult,
    minimize_with_forcefield_xml,
    validate_forcefield_xml,
)
from .clean_structure import (
    ADDITIVE_RESIDUES,
    BULK_ION_RESIDUES,
    STRUCTURAL_METAL_RESIDUES,
    WATER_RESIDUES,
    CleaningResult,
    clean_pdb,
    clean_topology,
)
from .ligand import build_ligand_xml
from .merge import merge_ffxml
from .structure import build_forcefield_xml
from .topology import extract_residue_to_pdb, find_nonstandard_residues

# Single source of truth for the version is pyproject.toml; installed metadata
# carries it. The fallback covers running from an uninstalled checkout.
try:
    __version__ = _metadata.version("forcefill")
except _metadata.PackageNotFoundError:  # uninstalled checkout
    __version__ = "0.0.0+unknown"

# The reading and conversion helpers stay behind their modules rather than being
# re-exported: `forcefill.ligand_files.inspect_ligand_file(...)` and
# `forcefill.smirnoff.installed_smirnoff_forcefields()` say where they belong,
# and the top level stays about the pipeline.
__all__ = [
    "ADDITIVE_RESIDUES",
    "BACKENDS",
    "BULK_ION_RESIDUES",
    "DEFAULT_AMBERTOOLS_TIMEOUT",
    "DEFAULT_BASE_FORCEFIELD",
    "DEFAULT_MINIMIZATION_PLATFORM",
    "DEFAULT_MINIMIZATION_TOLERANCE",
    "DEFAULT_SMIRNOFF_FORCEFIELD",
    "STRUCTURAL_METAL_RESIDUES",
    "WATER_RESIDUES",
    "CleaningResult",
    "LigandSpec",
    "MinimizationResult",
    "ParameterizationResult",
    "assemble_openmm_ffxml",
    "build_forcefield_xml",
    "build_ligand_xml",
    "clean_pdb",
    "clean_topology",
    "extract_residue_to_pdb",
    "find_nonstandard_residues",
    "locate_gaff_dat",
    "merge_ffxml",
    "minimize_with_forcefield_xml",
    "run_antechamber",
    "run_parmchk2",
    "validate_forcefield_xml",
]
