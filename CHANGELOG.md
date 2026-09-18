# Changelog

All notable changes to forcefill are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the version numbers
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html): the public API
is everything exported from `forcefill/__init__.py`, so a breaking change to any
of it requires a major bump.

## [Unreleased]

## [1.0.0] - 2026-09-18

First release. The API below is now covered by the compatibility promise above.

### Added

- **`build_forcefield_xml(pdb_file, output_xml, ...)`** — takes a PDB, finds
  every residue the base force field has no template for
  (`ForceField.getUnmatchedResidues`), parameterizes the ones a stand-alone
  small-molecule treatment is valid for, and writes one OpenMM ffxml you load
  alongside the standard force fields. Returns a `ParameterizationResult`
  carrying the per-residue XMLs, the skipped residues with their reasons, and
  the working directory.
- **`build_ligand_xml(ligands, output_xml, ...)`** — the same pipeline for
  ligands with no structure at all, from SDF, MOL2, SMILES or a `LigandSpec`,
  singly or as a batch.
- **Four backends**, selected with `backend=` globally or per ligand:
  - `"gaff"` (default) — `antechamber` GAFF2 atom types with AM1-BCC charges,
    plus `parmchk2` for the missing parameters.
  - `"smirnoff"` — OpenFF, default `openff-2.2.1`, through
    `openmmforcefields`. Accepts a released force field by name, a local
    OFFXML, or a layered list of both.
  - `"espaloma"` — the Espaloma graph network, default `espaloma-0.3.2`.
    Optional dependency; it pulls in PyTorch.
  - `"charmm"` — converts the CGenFF stream file ParamChem returns, against
    `CHARMM_BASE_FORCEFIELD`. Handles the LJ conversion that ParmEd's own
    CHARMM writer gets wrong.
  Backends may be mixed across residues in one run; the resulting XMLs merge
  into a single loadable file.
- **`LigandSpec`** — per-ligand settings (input file or SMILES, net charge,
  multiplicity, backend, force field, CHARMM files, atom type, charge method,
  extra `antechamber` arguments), overriding the run-wide defaults.
- **BespokeFit support** — a BespokeFit OFFXML is accepted wherever a SMIRNOFF
  force field is, with checks that the fit actually belongs to the molecule it
  is being applied to, so bought QC is not silently applied to the wrong ligand.
- **Preflight checks**, before anything expensive: the net charge is read from
  the ligand file, a supplied file is confirmed to be the same molecule as the
  residue it stands in for, and the geometry is screened for the faults that
  produce NaN energies.
- **Post-parameterization validation** — an `openmm.System` is built from
  `base force field + new XML` for every parameterized residue on its own, and
  for the whole input when nothing was skipped, so a template that does not
  match its residue fails at build time instead of at simulation time.
  `minimize=True` additionally catches parameters that are unphysical rather
  than merely absent. `validate_forcefield_xml`, `minimize_with_forcefield_xml`
  and `add_extra_particles` are public for use on their own.
- **`clean_pdb` / `clean_topology`** — drop water, bulk buffer ions and
  crystallization additives while keeping structural metals, either as a
  pre-step (`clean_structure=True`) or on their own. The residue sets
  (`WATER_RESIDUES`, `BULK_ION_RESIDUES`, `ADDITIVE_RESIDUES`,
  `STRUCTURAL_METAL_RESIDUES`) are public and inspectable.
- **Refusal with a reason** — chemistry a stand-alone treatment is not valid
  for (covalently bound residues, monatomic species, polymer fragments) is
  skipped and reported in `result.skipped`, rather than silently given wrong
  parameters.
- **Supporting public API** — `merge_ffxml`, `find_nonstandard_residues`,
  `extract_residue_to_pdb`, `assemble_openmm_ffxml`, `run_antechamber`,
  `run_parmchk2`, `locate_gaff_dat`, `residue_templates_with_virtual_sites`,
  and the `DEFAULT_*` constants.
- **Typing** — the package is fully annotated and ships a PEP 561 `py.typed`
  marker, so annotations are visible to type checkers in consuming projects.
- **Documentation** at [forcefill.readthedocs.io](https://forcefill.readthedocs.io),
  covering installation, a quickstart, the four backends, CHARMM/CGenFF,
  BespokeFit, cleaning, what gets skipped and why, and an API reference.

### Known limitations

- Covalently bound ligands are skipped deliberately: a stand-alone treatment of
  a polymer-linked residue is wrong whichever backend produces it.
- There is no command-line interface yet.
- Re-running into an existing working directory repeats finished `antechamber`
  jobs; nothing is cached between runs.

[Unreleased]: https://github.com/LouieSlocombe/forcefill/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/LouieSlocombe/forcefill/releases/tag/v1.0.0
