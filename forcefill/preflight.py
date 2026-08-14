"""Check and complete every ligand spec before anything expensive runs.

antechamber's AM1-BCC can take an hour on a drug-sized ligand. Three mistakes
waste it, and all three are visible in the input:

    * the **net charge** is left at 0 while the ligand file says the molecule is
      an amidinium - the charges come out plausible and wrong, and nothing
      downstream notices;
    * the supplied file is **not the same molecule** as the residue in the
      structure, which otherwise surfaces at the very end as an opaque OpenMM
      "no template matched";
    * atoms sit **on top of each other**, which surfaces as a NaN energy.

Running this as one pass up front also means a mistake in the last ligand does
not cost the parameterization of the first.

The reading and the checks themselves live in :mod:`forcefill.ligand_files`;
this module is what applies them to a set of specs. The checks that run *after*
parameterization, on the generated force field, are :mod:`forcefill.checks`.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

from openmm import app, unit

from . import ligand_files
from ._spec import ResolvedSpec
from .topology import _residue_positions, extract_residue_to_pdb

log = logging.getLogger(__name__)

__all__ = ["preflight_specs"]


def _warn_if_no_hydrogens(name: str, residue: app.topology.Residue) -> None:
    """Warn about a ligand extracted from a PDB that carries no hydrogens at all."""
    n_hydrogens = sum(1 for a in residue.atoms() if a.element is not None and a.element.symbol == "H")
    n_heavy = sum(1 for _ in residue.atoms()) - n_hydrogens
    if n_hydrogens == 0 and n_heavy > 1:
        log.warning(
            "Residue %s contains no hydrogens; AM1-BCC charges will be wrong "
            "unless the molecule really has none. Add explicit hydrogens to the "
            "ligand before parameterizing.",
            name,
        )


def _resolve_smiles(spec: ResolvedSpec, residue: app.topology.Residue | None, positions, res_dir: Path) -> ResolvedSpec:
    """Embed a spec's SMILES to an SDF and return the spec pointing at it.

    When the residue is also present in the input structure, the structure's
    coordinates are kept and only the bond orders come from the SMILES - the
    geometry a crystal structure gives is better than anything embedding
    produces, and the atom count has to match it anyway.
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = res_dir / f"{spec.name}_smiles.sdf"
    if residue is None or positions is None:
        return spec.with_file(ligand_files.smiles_to_sdf(spec.smiles, out_sdf, spec.name))
    residue_pdb = extract_residue_to_pdb(positions, residue, res_dir / f"{spec.name}_extracted.pdb")
    return spec.with_file(ligand_files.smiles_with_residue_geometry(spec.smiles, residue_pdb, out_sdf, spec.name))


def _apply_net_charge(spec: ResolvedSpec, info: ligand_files.LigandFileInfo) -> ResolvedSpec:
    """Fill in or cross-check the net charge against what the ligand file says.

    A file with real bond orders knows its own formal charge, so leaving
    ``net_charge`` unset is no longer a silent vote for 0. An explicit value that
    contradicts the file is refused outright: one of the two is wrong, and
    guessing which produces exactly the plausible-but-wrong charges this is meant
    to prevent.
    """
    if info.formal_charge is None:
        if spec.net_charge is None:
            log.warning(
                "Could not determine the net charge of %s from %s; assuming 0. "
                "Pass LigandSpec(net_charge=...) if that is wrong - it is the "
                "classic source of plausible but wrong AM1-BCC charges.",
                spec.name,
                Path(info.path).name,
            )
        return spec
    if spec.net_charge is None:
        log.warning(
            "Using net charge %+d for %s, read from %s. Pass LigandSpec(net_charge=...) to override it.",
            info.formal_charge,
            spec.name,
            Path(info.path).name,
        )
        return spec.with_net_charge(info.formal_charge)
    if spec.net_charge != info.formal_charge:
        raise ValueError(
            f"Residue {spec.name} was given net_charge={spec.net_charge:+d}, but "
            f"{Path(info.path).name} describes a molecule with a formal charge "
            f"of {info.formal_charge:+d}. One of them is wrong, and picking "
            "either would give charges that look reasonable and are not. Fix the "
            "protonation state in the ligand file, or correct net_charge."
        )
    return spec


def preflight_specs(
    specs: Mapping[str, ResolvedSpec],
    residues: Mapping[str, app.topology.Residue],
    positions: unit.Quantity | None,
    workdir: Path,
    *,
    strict: bool = True,
) -> dict[str, ResolvedSpec]:
    """Check and complete every spec before anything expensive runs.

    Returns the specs with SMILES resolved to files and inferred net charges
    filled in. See the module docstring for what is checked and why.

    Args:
        specs: Resolved specs, keyed by residue name.
        residues: The residues as they appear in the input structure; empty in
            standalone mode.
        positions: Coordinates for *residues*, or None in standalone mode.
        workdir: Where a SMILES-derived SDF is written.
        strict: Raise on a composition or geometry fault rather than warn.

    Returns:
        ``{residue_name: ResolvedSpec}``, ready to parameterize.
    """
    out: dict[str, ResolvedSpec] = {}
    for name in sorted(specs):
        spec = specs[name]
        residue = residues.get(name)
        if spec.smiles is not None:
            spec = _resolve_smiles(spec, residue, positions, workdir / name)

        if spec.file is None:
            if residue is not None:
                _warn_if_no_hydrogens(name, residue)
                if positions is not None:
                    ligand_files.check_geometry(
                        [tuple(p.value_in_unit(unit.angstrom)) for p in _residue_positions(positions, residue)],
                        name,
                        strict=strict,
                    )
            out[name] = spec
            continue

        info = ligand_files.inspect_ligand_file(spec.file)
        if residue is not None:
            ligand_files.check_matches_residue(info, residue, name, strict=strict)
        ligand_files.check_geometry(info.positions, name, strict=strict)
        out[name] = _apply_net_charge(spec, info)
    return out
