"""Parameterize a ligand from CHARMM/CGenFF files instead of GAFF or SMIRNOFF.

The third backend, and the only one that does not *derive* parameters. GAFF runs
antechamber, SMIRNOFF matches SMARTS patterns; CGenFF parameters come from the
``cgenff`` program or the ParamChem web service, neither of which is
redistributable. What this module does is convert what they emit - a CHARMM
stream file (``.str``) holding one ``RESI`` block and whatever parameters had to
be assigned by analogy - into an OpenMM ffxml the rest of forcefill can validate,
merge and minimize like any other.

**CHARMM is not interchangeable with Amber.** The two conventions disagree about
1-4 scaling (Amber 0.8333/0.5, CHARMM 1.0/1.0) and OpenMM refuses to load force
fields that disagree, so a CHARMM ligand needs
:data:`~forcefill.CHARMM_BASE_FORCEFIELD` underneath it and cannot share a build
with a gaff or smirnoff one. :func:`forcefill._pipeline.check_backends_match_base`
is what says so up front.

**The output is a residue template, not a self-contained force field.**
``charmm36.xml`` already carries every CGenFF atom type and its parameters, so
the generated XML names those types rather than redefining them. That is what
makes the backend usable with nothing but a ParamChem ``.str``: no CHARMM toppar
distribution to download, and no chance of a locally-supplied parameter file
silently overriding the base force field.

Three properties of ParmEd's writer shape :func:`charmm_residue_ffxml`, and each
one is a silent wrong answer rather than an error if it is not handled:

    * it writes ``sigma="1.0" epsilon="0.0"`` for every atom type whose
      non-bonded parameters it does not have, which *overrides* the real ones
      when the file is loaded after ``charmm36.xml``;
    * with ``separate_ljforce=False`` the ligand's Lennard-Jones energy is
      counted twice, once in ``NonbondedForce`` and once in the
      ``CustomNonbondedForce`` that ``charmm36.xml`` builds for its NBFIX pairs;
    * given a stream file that names atom types without defining their masses -
      which is every ParamChem stream file - it drops the residue template with
      only a warning and writes an empty document.

So the conversion injects the missing masses from the base force field, writes
with ``separate_ljforce=True``, strips back everything the base force field
already owns, and refuses to return a document with no residue in it.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType

from openmm import app
from parmed.charmm import CharmmParameterSet
from parmed.modeller import ResidueTemplate
from parmed.openmm import OpenMMParameterSet
from parmed.topologyobjects import AtomType

from ._spec import CHARMM_BASE_FORCEFIELD, PathLike, ResolvedSpec, check_charmm_suffix

log = logging.getLogger(__name__)

__all__ = [
    "atom_elements",
    "base_14_scales",
    "base_atom_types",
    "base_residue_names",
    "charmm_residue_ffxml",
    "ligand_topology",
    "read_charmm_files",
    "residue_template",
]

#: Sections of the generated XML whose per-type entries are pruned when the base
#: force field already defines the type. ``AtomTypes`` keys its children by
#: ``name``, the two force sections by ``class``; for a CHARMM parameter set
#: ParmEd writes the two identically, which is why one table covers both.
_PRUNED_SECTIONS = {"AtomTypes": "name", "NonbondedForce": "class", "LennardJonesForce": "class"}


@lru_cache(maxsize=4)
def _base_profile(
    base_forcefield: tuple[str, ...],
) -> tuple[dict[str, tuple[float, int]], tuple[float, float] | None, frozenset[str]]:
    """``(atom types, 1-4 scales, residue names)`` for a base force field.

    One parse answers every question the pipeline asks of the base force field,
    and it is cached: reading ``charmm36.xml`` is 6 MB of XML, and every residue
    in a build would otherwise pay for it again. These are static data files, so
    the cache cannot go stale within a process. The scales are None when the
    force field declares no non-bonded terms at all.
    """
    forcefield = app.ForceField(*base_forcefield)
    types = {}
    for name, atom_type in forcefield._atomTypes.items():
        element = atom_type.element
        # Lone pairs and Drude particles have no element; ParmEd writes those as
        # atomic number 0, which is also what an extra point is.
        types[name] = (atom_type.mass, element.atomic_number if element is not None else 0)
    nonbonded = next((g for g in forcefield.getGenerators() if type(g).__name__ == "NonbondedGenerator"), None)
    scales = None if nonbonded is None else (nonbonded.coulomb14scale, nonbonded.lj14scale)
    return types, scales, frozenset(forcefield._templates)


def base_atom_types(base_forcefield: Sequence[str] = CHARMM_BASE_FORCEFIELD) -> Mapping[str, tuple[float, int]]:
    """``{type name: (mass, atomic number)}`` for every atom type *base_forcefield* defines.

    The CGenFF atom types live in ``charmm36.xml``, so this is what lets a bare
    ParamChem stream file - which names types like ``CG331`` without ever saying
    what they weigh - be converted at all. Read-only: the table behind it is
    cached, and a caller editing it would change what every later build sees.
    """
    return MappingProxyType(_base_profile(tuple(base_forcefield))[0])


def base_14_scales(base_forcefield: Sequence[str]) -> tuple[float, float] | None:
    """``(coulomb14scale, lj14scale)`` a base force field declares, or None if it declares none.

    Read from the loaded force field rather than guessed from the file names,
    which is what makes the compatibility check in
    :func:`forcefill._pipeline.check_backends_match_base` work for a custom base
    force field as well as for the two named presets. ``None`` means there is
    nothing to be incompatible *with* - an empty base force field, or one with no
    non-bonded terms at all.
    """
    return _base_profile(tuple(base_forcefield))[1]


def base_residue_names(base_forcefield: Sequence[str] = CHARMM_BASE_FORCEFIELD) -> frozenset[str]:
    """Residue template names *base_forcefield* already defines.

    ``charmm36.xml`` carries 814 of them - every amino acid, nucleotide, lipid
    and CGenFF model compound - so a ligand named after one is much easier to
    hit here than under an Amber base force field, and worth catching by name.
    """
    return _base_profile(tuple(base_forcefield))[2]


def read_charmm_files(files: Sequence[PathLike]) -> CharmmParameterSet:
    """Read CHARMM topology/parameter/stream files into one parameter set.

    Args:
        files: ``.str``/``.rtf``/``.prm``/``.par``/``.top``/``.inp`` paths.
            ParmEd reads topologies first, then parameters, then stream files,
            whatever order they are given in.

    Returns:
        The combined :class:`parmed.charmm.CharmmParameterSet`.

    Raises:
        ValueError: No files were given, or one has a suffix ParmEd will not
            recognize.
    """
    files = list(files)
    if not files:
        raise ValueError("read_charmm_files needs at least one CHARMM file.")
    for path in files:
        check_charmm_suffix(path)
    # CharmmParameterSet dispatches with str.endswith, so a PathLike must be
    # converted before it gets there or every file is "unrecognized".
    return CharmmParameterSet(*[str(f) for f in files])


def residue_template(params: CharmmParameterSet, name: str, files: Sequence[PathLike]) -> ResidueTemplate:
    """Pick the residue *name* out of *params*, renaming a lone unnamed match.

    A stream file written for one ligand holds one ``RESI``, whose name is
    whatever was typed into ParamChem and need not be the residue name in the
    structure - so a single template is accepted and renamed, exactly as
    :func:`forcefill.amber.load_residue_template` renames a mol2. Anything
    ambiguous is refused: picking one of several would be a guess about which
    molecule is being parameterized.
    """
    residues = params.residues
    if not residues:
        raise ValueError(
            f"{[str(f) for f in files]} defines no residue template, so there is "
            f"nothing to parameterize {name} with. A CGenFF stream file needs a "
            "'read rtf card' section containing the RESI block for the ligand."
        )
    if name in residues:
        return residues[name]
    if len(residues) == 1:
        (found,) = residues.values()
        log.info("Renaming CHARMM residue template %r to %s.", found.name, name)
        found.name = name
        return found
    raise ValueError(
        f"{[str(f) for f in files]} defines {len(residues)} residue templates "
        f"({sorted(residues)}) and none of them is {name}. Name the ligand's RESI "
        f"{name} in the stream file, or use one of those names as the residue "
        "name, so which molecule is being parameterized is not a guess."
    )


def _inject_atom_types(params: CharmmParameterSet, template: ResidueTemplate, known: Mapping[str, tuple[float, int]]):
    """Give *params* an AtomType for every type *template* uses but does not define.

    Without this ParmEd's writer drops the residue template and produces an empty
    document, warning but not failing - and a ParamChem stream file never carries
    the ``MASS`` records that would avoid it, because the types it uses are
    CGenFF's own.
    """
    missing = sorted({atom.type for atom in template.atoms} - set(params.atom_types))
    unknown = [t for t in missing if t not in known]
    if unknown:
        raise ValueError(
            f"Residue {template.name} uses atom type(s) {unknown} that neither "
            "the CHARMM files nor the base force field define. Either the type "
            "names are misspelled, or they are new types whose masses and "
            "non-bonded parameters must be supplied - add the parameter file "
            "that defines them to charmm_files."
        )
    for index, name in enumerate(missing, start=1):
        mass, atomic_number = known[name]
        params.atom_types[name] = AtomType(name, index, mass, atomic_number=atomic_number)
    if missing:
        log.debug("Took masses for %s from the base force field: %s", template.name, missing)
    return missing


def _prune_base_definitions(xml_file: Path, known: Mapping[str, tuple[float, int]]) -> list[str]:
    """Strip everything the base force field already owns; returns what went.

    ParmEd writes an ``<AtomTypes>`` entry and a zero-valued non-bonded entry for
    every type the residue uses, including the ones it only knows the mass of.
    Left in place, loading the file after ``charmm36.xml`` redefines those types
    and replaces their real Lennard-Jones parameters with ``epsilon=0`` - silently,
    since OpenMM treats a later definition as an override rather than a clash.
    """
    root = ET.parse(str(xml_file)).getroot()
    dropped: list[str] = []
    for section in list(root):
        key = _PRUNED_SECTIONS.get(section.tag)
        if key is not None:
            for child in list(section):
                if child.get(key) in known:
                    section.remove(child)
                    dropped.append(f"{section.tag}/{child.get(key)}")
        # ParmEd emits empty sections (an <AmoebaUreyBradleyForce/> for a ligand
        # with no Urey-Bradley terms); pruning creates more. Neither carries
        # anything, and <NonbondedForce> survives on its <UseAttributeFromResidue>.
        if len(section) == 0 and section.tag != "Info":
            root.remove(section)
    ET.indent(root)
    ET.ElementTree(root).write(str(xml_file), encoding="utf-8", xml_declaration=True)
    return dropped


def _require_residue(xml_file: Path, name: str, files: Sequence[PathLike]) -> None:
    """Raise unless the written document still contains the residue template."""
    root = ET.parse(str(xml_file)).getroot()
    if any(residue.get("name") == name for residue in root.findall("./Residues/Residue")):
        return
    raise RuntimeError(
        f"The force field written for {name} from {[str(f) for f in files]} "
        "contains no residue template. ParmEd drops a template whose atom types "
        "it cannot resolve, warning rather than failing, so this is what an "
        "unusable stream file looks like. Check that the RESI block's atom types "
        "are real CGenFF types."
    )


def charmm_residue_ffxml(
    spec: ResolvedSpec,
    output_xml: PathLike,
    base_forcefield: Sequence[str] = CHARMM_BASE_FORCEFIELD,
) -> str:
    """Convert one ligand's CHARMM files into an OpenMM force-field XML and return the path.

    The result names CGenFF's own atom types rather than redefining them, so it
    is meaningful only when loaded on top of *base_forcefield*. It carries the
    residue template plus any parameters the CHARMM files themselves define -
    for a ParamChem stream file, the terms it had to assign by analogy.

    Args:
        spec: The ligand. Must carry ``charmm_files``; ``atom_type``,
            ``charge_method`` and ``forcefield`` do not apply and are ignored.
        output_xml: Where to write the per-residue XML.
        base_forcefield: The CHARMM force field this will be loaded with, whose
            atom types the output is allowed to refer to rather than define.

    Returns:
        The path written, as a string.

    Raises:
        ValueError: The CHARMM files are unreadable, define no usable residue
            template, or use an atom type nothing defines.
        RuntimeError: The conversion produced a document with no residue in it.
    """
    params = read_charmm_files(spec.charmm_files)
    template = residue_template(params, spec.name, spec.charmm_files)
    if spec.name in base_residue_names(base_forcefield):
        raise ValueError(
            f"The base force field {list(base_forcefield)} already defines a "
            f"residue template named {spec.name}, so adding another would make "
            "OpenMM refuse to load them together. charmm36.xml carries every "
            "amino acid, nucleotide and CGenFF model compound, so this is easy "
            f"to hit by accident - give the ligand a different residue name in "
            "both the structure and the stream file."
        )
    known = base_atom_types(base_forcefield)
    _inject_atom_types(params, template, known)

    net_charge = round(sum(atom.charge for atom in template.atoms))
    log.info(
        "charmm: %s (net charge %+d, %d atoms, from %s)",
        spec.name,
        net_charge,
        len(template.atoms),
        ", ".join(Path(f).name for f in spec.charmm_files),
    )

    omm_params = OpenMMParameterSet.from_parameterset(params)
    # A toppar file carries hundreds of residues; only the ligand's belongs here,
    # and writing the rest collides with the templates charmm36.xml already has.
    omm_params.residues = {spec.name: template}

    output_xml = Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    omm_params.write(
        str(output_xml),
        provenance={"Info": "Generated by forcefill (CHARMM/CGenFF via ParmEd)"},
        # Not optional: charmm36.xml puts Lennard-Jones interactions in a
        # separate force for its NBFIX pairs, and a NonbondedForce that also
        # carries them makes every ligand LJ term count twice.
        separate_ljforce=True,
        write_unused=False,
    )
    dropped = _prune_base_definitions(output_xml, known)
    log.debug("Dropped %d definition(s) the base force field already owns: %s", len(dropped), dropped)
    _require_residue(output_xml, spec.name, spec.charmm_files)
    log.info("Wrote per-residue XML: %s", output_xml)
    return str(output_xml)


def atom_elements(
    template: ResidueTemplate,
    base_forcefield: Sequence[str] = CHARMM_BASE_FORCEFIELD,
) -> list[app.Element | None]:
    """Elements for a residue template's atoms, in order; None where undeterminable.

    A ``RESI`` block names atom *types*, not elements, and ParmEd leaves
    ``Atom.atomic_number`` at 0 unless the same files also carried the ``MASS``
    records - which a ParamChem stream file does not. So the base force field is
    what turns ``CG331`` back into carbon.
    """
    known = base_atom_types(base_forcefield)
    elements = []
    for atom in template.atoms:
        atomic_number = atom.atomic_number or known.get(atom.type, (0.0, 0))[1]
        elements.append(app.Element.getByAtomicNumber(atomic_number) if atomic_number else None)
    return elements


def ligand_topology(
    spec: ResolvedSpec,
    base_forcefield: Sequence[str] = CHARMM_BASE_FORCEFIELD,
) -> app.Topology:
    """Return a Topology for the spec's ligand, for validating it on its own.

    Coordinates deliberately absent: a CHARMM stream file records internal
    coordinates rather than Cartesian ones, and ParmEd leaves the template's
    positions at the origin. Validation only needs the bond graph; minimization
    needs real geometry, which means starting from a structure with
    :func:`~forcefill.build_forcefield_xml`.
    """
    params = read_charmm_files(spec.charmm_files)
    template = residue_template(params, spec.name, spec.charmm_files)
    topology = app.Topology()
    residue = topology.addResidue(spec.name, topology.addChain("A"))
    atoms = {
        atom.name: topology.addAtom(atom.name, element, residue)
        for atom, element in zip(template.atoms, atom_elements(template, base_forcefield), strict=True)
    }
    for bond in template.bonds:
        topology.addBond(atoms[bond.atom1.name], atoms[bond.atom2.name])
    return topology
