from importlib import metadata as _metadata

from ._pipeline import ParameterizationResult
from ._spec import (
    BACKENDS,
    CHARMM_BASE_FORCEFIELD,
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

# The version lives in pyproject.toml and reaches here through the installed
# metadata; the fallback covers an uninstalled checkout.
try:
    __version__ = _metadata.version("forcefill")
except _metadata.PackageNotFoundError:  # uninstalled checkout
    __version__ = "0.0.0+unknown"

# The reading and conversion helpers are not re-exported: they say where they
# belong (`forcefill.ligand_files.inspect_ligand_file(...)`,
# `forcefill.charmm.read_charmm_files(...)`), and the top level stays about the
# pipeline.
# `clean_structure` is also the name of a *module* here; it stays out of
# __all__ because at the top level that name reads as build_forcefield_xml's
# `clean_structure=` switch. Its public names are re-exported above.
__all__ = [
    "ADDITIVE_RESIDUES",
    "BACKENDS",
    "BULK_ION_RESIDUES",
    "CHARMM_BASE_FORCEFIELD",
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
