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
      already been paid for. :func:`check_forcefield_applies` says so, comparing
      every handler's assignments - not only the torsions - against the release
      the file layers on, and standing down entirely for a selection that is
      released chemistry named by path.

Requires openmmforcefields >= 0.16, which is where ``smirnoff_filenames``, a
multi-file ``forcefield=`` selection, and constraints and virtual sites in the
generated template all arrive. Earlier releases strip constraints and take a
single force field only, so both the check above and the layering below would be
silently or loudly wrong against them.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING

from openff.toolkit import ForceField, Molecule
from openff.toolkit.typing.engines.smirnoff import get_available_force_fields
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
    "is_stock_forcefield",
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
        raise ValueError(
            f"openmmforcefields parsed {unresolved} but it corresponds to no file "
            "on disk, so forcefill cannot record which force field the generated "
            "XML came from. That happens when a force field is given as content "
            "rather than as a name or a path - write it to a '.offxml' file and "
            "pass that instead."
        )
    return list(paths)


@cache
def _stock_forcefield_paths() -> frozenset[str]:
    """Resolved paths of every OFFXML the installed OpenFF packages ship.

    Read from the toolkit rather than from a list of release names: the released
    force fields live in an installed data package, so "is this file a release?"
    is answered exactly by where it sits, with nothing parsed.
    """
    return frozenset(str(Path(path).resolve()) for path in get_available_force_fields(full_paths=True))


@cache
def _stock_forcefield_names() -> frozenset[str]:
    """File names of the released OFFXMLs, for an entry the toolkit resolves itself."""
    return frozenset(Path(path).name for path in _stock_forcefield_paths())


def _is_stock_entry(entry: str) -> bool:
    """True when *entry* names a released OpenFF force field rather than a local file.

    Decided by resolved path wherever there is one, not by name alone: a file in
    the working directory called ``openff_unconstrained-2.2.1.offxml`` is a
    bespoke force field with a confusing name, and calling it stock would skip
    the very check it needs.
    """
    if entry in installed_smirnoff_forcefields():
        return True
    path = Path(entry)
    if path.is_file():
        return str(path.resolve()) in _stock_forcefield_paths()
    # Nothing here by that name, so the toolkit will resolve it from its own
    # search path - which holds released force fields and nothing else.
    return path.name in _stock_forcefield_names()


def is_stock_forcefield(forcefield: Sequence[str]) -> bool:
    """True when every entry in *forcefield* is a released OpenFF force field.

    Distinct from the negation of :func:`is_custom_forcefield`, which asks
    whether openmmforcefields will recognize the *name*. A release can also be
    named by path - ``openff_unconstrained-2.3.0.offxml``, or any release newer
    than the installed openmmforcefields has an alias for - and that is still
    stock chemistry, with no bespoke parameters for
    :func:`check_forcefield_applies` to look for.
    """
    return bool(forcefield) and all(_is_stock_entry(entry) for entry in forcefield)


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
def _stock_smirks(release: str) -> frozenset[tuple[str, str]]:
    """Every ``(handler, SMIRKS)`` pair an installed release defines.

    Every handler, not just ``ProperTorsions``: a force field can be customized
    in its charges, its vdW, its bonded terms or its virtual sites, and reading
    only the torsions would call all of those "nothing a stock release would not
    give you".

    Cached like :func:`forcefill.charmm._base_profile`, and for the same reason:
    a Sage OFFXML is half a megabyte and the answer is the same for every ligand
    in a build.
    """
    stock = ForceField(*resolve_forcefield_files((release,)))
    return frozenset(
        (handler, parameter.smirks)
        for handler in stock.registered_parameter_handlers
        for parameter in getattr(stock.get_parameter_handler(handler), "parameters", ())
    )


def _assigned_smirks(forcefield: ForceField, molecule: Molecule) -> frozenset[tuple[str, str]]:
    """Every ``(handler, SMIRKS)`` pair *forcefield* actually assigns to *molecule*.

    What was assigned, not what the file contains: a parameter that matches
    nothing in this molecule is exactly the case being looked for.

    Most handlers map each atom tuple to one parameter, but ``VirtualSites``
    maps it to a *list* of them - one atom can carry several extra sites - so
    every value is flattened before its SMIRKS is read. Reading ``.smirks`` off
    the list instead would silently drop the whole handler.
    """
    labels = forcefield.label_molecules(molecule.to_topology())[0]
    return frozenset(
        (handler, parameter.smirks)
        for handler, assignments in labels.items()
        for value in assignments.values()
        for parameter in (value if isinstance(value, (list, tuple)) else [value])
        if hasattr(parameter, "smirks")
    )


def _reference_release(selection: Sequence[str]) -> str:
    """The released force field a custom *selection* is measured against.

    A bespoke file is nearly always layered on a stock release
    (``["openff-2.2.1", "bespoke.offxml"]``), and that release is the honest
    baseline - measuring against :data:`DEFAULT_SMIRNOFF_FORCEFIELD` instead
    would credit the bespoke file with every parameter the two releases happen to
    differ by. The last stock entry wins, matching how the toolkit layers them.
    """
    stock = [entry for entry in selection if _is_stock_entry(entry)]
    return stock[-1] if stock else DEFAULT_SMIRNOFF_FORCEFIELD


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

    The signal is a parameter assigned from a pattern the stock release does not
    contain: a bespoke SMIRKS never appears in a stock force field, so a genuine
    bespoke file always has at least one, even when it was fitted against a
    different Sage version. Version skew can only make this check *quieter*,
    never make it fire wrongly.

    Two things keep it from firing on a force field that is merely not the
    default. A selection that is entirely released chemistry - a release named by
    path, or one newer than the installed openmmforcefields has an alias for - has
    no bespoke parameters to look for and is not examined at all. And the
    comparison covers every handler, so a file that customizes charges, vdW,
    bonded terms or virtual sites rather than torsions is recognized as
    contributing something.
    """
    if is_stock_forcefield(selection):
        log.info(
            "%s: %s is released OpenFF chemistry, so there are no bespoke parameters to check for.",
            name,
            list(selection),
        )
        return
    reference = _reference_release(selection)
    assigned = _assigned_smirks(forcefield, molecule)
    extra = assigned - _stock_smirks(reference)
    if extra:
        handlers = ", ".join(sorted({handler for handler, _ in extra}))
        log.info(
            "%s: %d of %d assigned parameters come from %s (%s)",
            name,
            len(extra),
            len(assigned),
            list(selection),
            handlers,
        )
        return
    message = (
        f"The SMIRNOFF force field {list(selection)} gives residue {name} nothing "
        f"that {reference} would not: all {len(assigned)} of the parameters it "
        "assigns come from stock patterns. A bespoke force field matches by "
        "SMIRKS alone, so one fitted for a different molecule falls back to stock "
        "parameters silently. Check this is the force field fitted for "
        f"{name}, or drop it and use the release directly."
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
