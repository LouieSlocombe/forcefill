"""Parameterize a ligand with a SMIRNOFF force field (OpenFF Sage) instead of GAFF.

The alternative to the AmberTools path in :mod:`forcefill.amber`. SMIRNOFF
assigns parameters by matching SMARTS patterns against the chemical graph, so
there are no atom types to get wrong - but also no way to guess the graph from
coordinates. Hence the one requirement this backend adds: **a ligand must arrive
as a file with bond orders (SDF/MOL2) or as a SMILES**, never as a bare PDB
residue.

The work is done by :class:`openmmforcefields.generators.SMIRNOFFTemplateGenerator`,
whose ``generate_residue_template`` returns a self-contained ffxml string for one
molecule. Two properties of that output shape this module:

    * its atom types are named by a hash of the molecule, so two SMIRNOFF ffxmls
      - or a SMIRNOFF and a GAFF one - never collide and can be merged with
      :func:`~forcefill.merge_ffxml`;
    * it names the residue template with a mapped SMILES, so
      :func:`smirnoff_residue_ffxml` rewrites it to the residue name the rest of
      forcefill uses.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from pathlib import Path

from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from ._spec import PathLike, ResolvedSpec

log = logging.getLogger(__name__)

__all__ = ["installed_smirnoff_forcefields", "ligand_topology", "smirnoff_residue_ffxml"]

#: File suffixes ``openff.toolkit.Molecule.from_file`` reads reliably. MOL2 is
#: accepted but discouraged: the toolkit reads it through RDKit, whose MOL2
#: parser rejects the GAFF-typed files antechamber writes.
_MOLECULE_FORMATS = {".sdf", ".sd", ".mol", ".mol2"}


def installed_smirnoff_forcefields() -> list[str]:
    """Names of the SMIRNOFF releases available locally, e.g. ``['openff-2.2.1', ...]``."""
    return list(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)


def _load_molecule(spec: ResolvedSpec):  # -> an openff Molecule
    """Build an OpenFF Molecule from the spec's file or SMILES, with a 3D conformer."""
    if spec.smiles is not None:
        molecule = Molecule.from_smiles(spec.smiles, allow_undefined_stereo=True)
    else:
        path = Path(spec.file)
        if path.suffix.lower() not in _MOLECULE_FORMATS:
            raise ValueError(
                f"The smirnoff backend cannot read {path.name!r} for residue "
                f"{spec.name}: it needs a file carrying bond orders "
                f"(one of {sorted(_MOLECULE_FORMATS)}), and {path.suffix or 'no suffix'} "
                "does not. A PDB records no bond orders at all - supply an SDF, "
                "a SMILES, or use backend='gaff'."
            )
        molecule = Molecule.from_file(str(path), allow_undefined_stereo=True)
        if isinstance(molecule, list):
            if len(molecule) != 1:
                raise ValueError(
                    f"{path} holds {len(molecule)} molecules but residue "
                    f"{spec.name} needs exactly one. Split it first with "
                    "forcefill.ligand_files.split_multi_sdf()."
                )
            molecule = molecule[0]

    if spec.net_charge is not None:
        total = round(molecule.total_charge.m)
        if total != spec.net_charge:
            raise ValueError(
                f"Residue {spec.name} was given net_charge={spec.net_charge:+d}, "
                f"but its structure has a formal charge of {total:+d}. SMIRNOFF "
                "takes the charge from the chemical graph, so the two must "
                "agree - fix the protonation in the ligand file/SMILES, or drop "
                "net_charge and let it be read from there."
            )

    if not molecule.n_conformers:
        # Charges are conformer-dependent; letting the toolkit pick silently
        # would make the output depend on an invisible default.
        log.info("Generating a conformer for %s: the input carried none.", spec.name)
        molecule.generate_conformers(n_conformers=1)
    molecule.name = spec.name
    return molecule


def _rename_residue_template(ffxml: str, name: str) -> str:
    """Rewrite the generated template's residue name to *name*.

    openmmforcefields names it with a mapped SMILES. OpenMM matches templates by
    graph, so that would work - but the name is what appears in every error
    message, in the merged XML and in any ``registerTemplate`` override, and a
    60-character SMILES there is useless.
    """
    root = ET.fromstring(ffxml)
    residues = root.findall("./Residues/Residue")
    if len(residues) != 1:
        raise RuntimeError(
            f"Expected exactly one residue template for {name} from "
            f"openmmforcefields, got {len(residues)}. This is a change in that "
            "library's output; report it against forcefill."
        )
    log.debug("Renaming SMIRNOFF template %r to %s", residues[0].get("name"), name)
    residues[0].set("name", name)
    return ET.tostring(root, encoding="unicode")


def ligand_topology(spec: ResolvedSpec):  # -> (openmm Topology, Quantity)
    """Return ``(topology, positions)`` for the spec's molecule, for validating it on its own.

    Standalone mode has no input structure to check the generated template
    against, so the molecule supplies one. Rebuilding it here rather than
    carrying it out of :func:`smirnoff_residue_ffxml` costs a file read and no
    charge assignment, which is the expensive half.
    """
    molecule = _load_molecule(spec)
    topology = molecule.to_topology().to_openmm()
    # OpenFF names the residue after the molecule or leaves it UNK; matching is
    # by graph either way, but the name is what error messages print.
    for residue in topology.residues():
        residue.name = spec.name
    return topology, molecule.conformers[0].to_openmm()


def smirnoff_residue_ffxml(spec: ResolvedSpec, output_xml: PathLike) -> str:
    """Write a SMIRNOFF force-field XML for one ligand and return the path.

    Args:
        spec: The ligand. Must carry ``file`` or ``smiles``; ``atom_type`` and
            ``charge_method`` do not apply and are ignored (SMIRNOFF has neither
            atom types nor a choice of charge model at this level).
        output_xml: Where to write the per-residue XML.

    Returns:
        The path written, as a string.

    Raises:
        RuntimeError: openmmforcefields produced something unexpected.
        ValueError: The ligand source is unusable, or its formal charge
            contradicts an explicit ``net_charge``.
    """
    molecule = _load_molecule(spec)
    log.info(
        "smirnoff: %s (net charge %+d, %s)",
        spec.name,
        round(molecule.total_charge.m),
        spec.forcefield,
    )
    generator = SMIRNOFFTemplateGenerator(molecules=[molecule], forcefield=spec.forcefield)
    ffxml = _rename_residue_template(generator.generate_residue_template(molecule), spec.name)

    output_xml = Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    output_xml.write_text(ffxml)
    log.info("Wrote per-residue XML: %s", output_xml)
    return str(output_xml)
