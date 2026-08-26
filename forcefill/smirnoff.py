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

``forcefield`` need not be an installed release. openmmforcefields takes a path
to an OFFXML through the same argument, which is how a bespoke force field from
`BespokeFit <https://github.com/openforcefield/openff-bespokefit>`_ - Sage plus
torsions fitted to QC data for one molecule - is used here. forcefill neither
runs nor requires BespokeFit; it consumes the file, and checks the two things
that otherwise go wrong silently:

    * **constraints.** openmmforcefields copies every constraint the force field
      assigns into the residue template, so a *constrained* OFFXML produces a
      ligand whose X-H bonds are rigid no matter what ``createSystem`` was asked
      for. :func:`check_unconstrained` refuses one - BespokeFit's own default
      starting point is ``openff_unconstrained-*``, so this only catches the
      override.
    * **the wrong file for this ligand.** Bespoke parameters are identified by
      nothing but their SMIRKS, so pairing molecule A's force field with
      molecule B silently falls back to stock parameters after the QC has
      already been paid for. :func:`check_forcefield_applies` says so.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING

from openff.toolkit import ForceField, Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from ._spec import DEFAULT_SMIRNOFF_FORCEFIELD, PathLike, ResolvedSpec

if TYPE_CHECKING:
    from openmm import app, unit

log = logging.getLogger(__name__)

__all__ = [
    "check_custom_forcefield",
    "check_forcefield_applies",
    "check_unconstrained",
    "forcefield_14_scales",
    "installed_smirnoff_forcefields",
    "is_custom_forcefield",
    "ligand_topology",
    "load_forcefield",
    "resolve_forcefield_files",
    "smirnoff_residue_ffxml",
]

#: File suffixes ``openff.toolkit.Molecule.from_file`` reads reliably. MOL2 is
#: accepted but discouraged: the toolkit reads it through RDKit, whose MOL2
#: parser rejects the GAFF-typed files antechamber writes.
_MOLECULE_FORMATS = {".sdf", ".sd", ".mol", ".mol2"}


def installed_smirnoff_forcefields() -> list[str]:
    """Names of the SMIRNOFF releases available locally, e.g. ``['openff-2.2.1', ...]``."""
    return list(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)


def is_custom_forcefield(forcefield: Sequence[str]) -> bool:
    """True unless *forcefield* is a single installed SMIRNOFF release name.

    The same test openmmforcefields makes when it decides whether to resolve an
    entry against its own table or hand it to the toolkit as a path. One
    recognized name is a release forcefill can trust; anything else - a path, or
    several entries layered - is a file it has to check.
    """
    return len(forcefield) != 1 or forcefield[0] not in installed_smirnoff_forcefields()


def resolve_forcefield_files(forcefield: Sequence[str]) -> list[str]:
    """The OFFXML files openmmforcefields will actually read for this selection.

    A release name is an openmmforcefields *alias*: ``"openff-2.2.1"`` resolves
    to ``openff_unconstrained-2.2.1.offxml``, a filename the toolkit knows and
    the alias one it does not. Going through the generator rather than
    reimplementing the table means forcefill checks the same bytes that will be
    parameterized from.
    """
    generator = SMIRNOFFTemplateGenerator(forcefield=list(forcefield))
    paths = generator.smirnoff_filenames
    if any(path is None for path in paths):
        unresolved = [entry for entry, path in zip(forcefield, paths, strict=True) if path is None]
        raise RuntimeError(
            f"openmmforcefields loaded {list(forcefield)} but could not say which "
            f"file {unresolved} came from. This is a change in that library's "
            "output; report it against forcefill."
        )
    return list(paths)


def load_forcefield(forcefield: Sequence[str]) -> ForceField:
    """Load a SMIRNOFF force field, naming the entry at fault if it will not load.

    openmmforcefields raises for the whole selection at once, and only once a
    molecule has been read. Loading here instead means a mistyped path in the
    fifth ligand does not cost the AM1-BCC charges of the first four.
    """
    try:
        return ForceField(*resolve_forcefield_files(forcefield))
    except Exception as exc:
        missing = [entry for entry in forcefield if entry.endswith(".offxml") and not Path(entry).is_file()]
        hint = (
            f" No file exists at {missing[0]!r}."
            if missing
            else " It is neither an installed release name nor a readable OFFXML document."
        )
        raise ValueError(
            f"Could not load the SMIRNOFF force field {list(forcefield)}.{hint} "
            f"Name an installed release (..., {installed_smirnoff_forcefields()[-3:]}) "
            "or give the path to an OFFXML file - a bespoke one from BespokeFit, "
            f"for instance.\n  {type(exc).__name__}: {exc}"
        ) from exc


@cache
def _stock_torsion_smirks(release: str) -> frozenset[str]:
    """Every proper-torsion SMIRKS in an installed release.

    Cached like :func:`forcefill.charmm._base_profile`, and for the same reason:
    a Sage OFFXML is half a megabyte and the answer is the same for every ligand
    in a build.
    """
    stock = ForceField(*resolve_forcefield_files((release,)))
    return frozenset(parameter.smirks for parameter in stock.get_parameter_handler("ProperTorsions").parameters)


def forcefield_14_scales(forcefield: ForceField) -> tuple[float, float]:
    """The ``(coulomb, lj)`` 1-4 scaling *forcefield* declares.

    What openmmforcefields writes into the generated ffxml, and so what has to
    agree with the base force field. Every stock Sage release says 0.8333/0.5,
    but an arbitrary OFFXML is not obliged to.
    """
    return (
        float(forcefield.get_parameter_handler("Electrostatics").scale14),
        float(forcefield.get_parameter_handler("vdW").scale14),
    )


def check_custom_forcefield(spec: ResolvedSpec, forcefield: ForceField, *, strict: bool = True) -> None:
    """Run both custom-OFFXML checks for one ligand.

    What :mod:`forcefill.preflight` calls. Reads the molecule again rather than
    threading it out of :func:`smirnoff_residue_ffxml`: that costs a file read
    and no charge assignment, which is the expensive half, and keeps every
    preflight check in one pass over the specs.
    """
    molecule = _load_molecule(spec)
    check_unconstrained(forcefield, molecule, spec.name, spec.forcefield)
    check_forcefield_applies(forcefield, molecule, spec.name, spec.forcefield, strict=strict)


def check_unconstrained(forcefield: ForceField, molecule: Molecule, name: str, selection: Sequence[str]) -> None:
    """Raise if *forcefield* constrains any bond in *molecule*.

    openmmforcefields copies every constraint the force field assigns into the
    residue template (as ``<Constraint>`` elements), so a constrained OFFXML
    yields a ligand whose X-H bonds are rigid regardless of the ``constraints``
    argument the caller later passes to ``createSystem`` - and OpenMM reports
    nothing, because a template is entitled to declare constraints.

    The test has to be per-molecule, not per-file: the *unconstrained* releases
    still carry the two rigid-TIP3P patterns, which no ligand matches. Only the
    constrained ones add ``[#1:1]-[*:2]``, which every X-H does.
    """
    constrained = forcefield.label_molecules(molecule.to_topology())[0]["Constraints"]
    if not constrained:
        return
    raise ValueError(
        f"The SMIRNOFF force field {list(selection)} constrains {len(constrained)} "
        f"bond(s) of residue {name}. openmmforcefields writes those into the "
        "residue template, so the ligand's bonds would stay rigid whatever "
        "constraints=... createSystem is given, and nothing would report the "
        "disagreement. Use the unconstrained build of the same force field - the "
        "OpenFF releases ship as 'openff_unconstrained-<version>.offxml' "
        "alongside 'openff-<version>.offxml', and BespokeFit fits against the "
        "unconstrained one by default."
    )


def check_forcefield_applies(
    forcefield: ForceField, molecule: Molecule, name: str, selection: Sequence[str], *, strict: bool = True
) -> None:
    """Check a custom force field actually contributes something to *molecule*.

    A bespoke parameter is identified by nothing but its SMIRKS, so a force field
    fitted for one molecule and applied to another simply fails to match and
    falls back to the stock parameters underneath it - after the quantum
    chemistry has been paid for, and with nothing anywhere reporting it.

    The signal is a proper torsion assigned from a pattern the stock release does
    not contain: a bespoke SMIRKS never appears in a stock force field, so a
    genuine bespoke file always has at least one, even when it was fitted against
    a different Sage version. Version skew can only make this check *quieter*,
    never make it fire wrongly.
    """
    labels = forcefield.label_molecules(molecule.to_topology())[0]["ProperTorsions"]
    assigned = {parameter.smirks for parameter in labels.values()}
    extra = assigned - _stock_torsion_smirks(DEFAULT_SMIRNOFF_FORCEFIELD)
    if extra:
        log.info("%s: %d of %d proper torsions come from %s", name, len(extra), len(assigned), list(selection))
        return
    message = (
        f"The SMIRNOFF force field {list(selection)} gives residue {name} nothing "
        f"that {DEFAULT_SMIRNOFF_FORCEFIELD} would not: all {len(assigned)} of its "
        "proper torsions were assigned from stock patterns. A bespoke force field "
        "matches by SMIRKS alone, so one fitted for a different molecule falls "
        "back to stock parameters silently. Check this is the force field fitted "
        f"for {name}, or drop it and use the release directly."
    )
    if strict:
        raise ValueError(message)
    log.warning("%s", message)


def _load_molecule(spec: ResolvedSpec) -> Molecule:
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


def ligand_topology(spec: ResolvedSpec) -> tuple[app.Topology, unit.Quantity]:
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
        ", ".join(spec.forcefield),
    )
    generator = SMIRNOFFTemplateGenerator(molecules=[molecule], forcefield=list(spec.forcefield))
    ffxml = _rename_residue_template(generator.generate_residue_template(molecule), spec.name)

    output_xml = Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    output_xml.write_text(ffxml)
    log.info("Wrote per-residue XML: %s", output_xml)
    return str(output_xml)
