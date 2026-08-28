# API reference

Thirty-one names are importable straight from `forcefill`. Everything else is
reached through its module — `forcefill.ligand_files.inspect_ligand_file(...)`,
`forcefill.smirnoff.installed_smirnoff_forcefields()` — which is deliberate: the
top level stays about the pipeline, and the module name says where a helper
belongs.

## Where each top-level name lives

{doc}`entrypoints`
: `build_forcefield_xml`, `build_ligand_xml`, `LigandSpec`,
  `ParameterizationResult`

{doc}`constants`
: `BACKENDS`, `DEFAULT_BASE_FORCEFIELD`, `CHARMM_BASE_FORCEFIELD`,
  `DEFAULT_SMIRNOFF_FORCEFIELD`, `DEFAULT_ESPALOMA_FORCEFIELD`,
  `WATER_RESIDUES`, `BULK_ION_RESIDUES`, `STRUCTURAL_METAL_RESIDUES`,
  `ADDITIVE_RESIDUES`

{doc}`topology`
: `find_nonstandard_residues`, `extract_residue_to_pdb`

{doc}`clean_structure`
: `clean_pdb`, `clean_topology`, `CleaningResult`

{doc}`checks`
: `validate_forcefield_xml`, `minimize_with_forcefield_xml`,
  `MinimizationResult`, `add_extra_particles`,
  `residue_templates_with_virtual_sites`, `DEFAULT_MINIMIZATION_PLATFORM`,
  `DEFAULT_MINIMIZATION_TOLERANCE`

{doc}`amber`
: `run_antechamber`, `run_parmchk2`, `locate_gaff_dat`,
  `assemble_openmm_ffxml`, `DEFAULT_AMBERTOOLS_TIMEOUT`

{doc}`merge`
: `merge_ffxml`

## Modules with a public surface of their own

{doc}`ligand_files`, {doc}`preflight`, {doc}`smirnoff`, {doc}`espaloma` and
{doc}`charmm` are not re-exported at the top level, but each has documented
functions worth calling directly — reading a ligand file, listing the installed
SMIRNOFF releases or Espaloma models, or converting a CGenFF stream file on its
own.

```{toctree}
:maxdepth: 2

entrypoints
constants
topology
clean_structure
preflight
checks
ligand_files
merge
amber
smirnoff
espaloma
charmm
```
