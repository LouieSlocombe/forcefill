"""Merge finished OpenMM force-field XML documents into one file.

The counterpart of :func:`forcefill.amber.assemble_openmm_ffxml`, which merges
at the ParmEd parameter-set level and so handles Amber-style input only. Mixing
backends needs this instead: a GAFF ffxml and a SMIRNOFF ffxml have nothing in
common upstream of the XML.

Pure ``xml.etree`` - it knows nothing about GAFF, SMIRNOFF or the pipeline, only
what OpenMM will accept when it reads the result back.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from pathlib import Path

from ._spec import PathLike

log = logging.getLogger(__name__)

__all__ = ["merge_ffxml"]

#: Tolerance for comparing numeric force-field attributes when merging. Matches
#: ``openmm.app.forcefield.NonbondedGenerator.SCALETOL``, so a merge that
#: succeeds here is one OpenMM would also accept across separate files.
_SCALE_TOLERANCE = 1e-5


def _attributes_compatible(a: Mapping[str, str], b: Mapping[str, str]) -> bool:
    """True when two force sections can be folded into one element.

    Numeric attributes are compared with a tolerance: the same constant is
    written to different precision by different producers, e.g. Amber's 1-4
    Coulomb scale is ``0.8333333333333334`` from ParmEd and ``0.8333333333``
    from openmmforcefields. Anything else must match exactly.
    """
    if set(a) != set(b):
        return False
    for key, left in a.items():
        right = b[key]
        if left == right:
            continue
        try:
            if abs(float(left) - float(right)) > _SCALE_TOLERANCE:
                return False
        except ValueError:
            return False
    return True


def _check_no_redefinition(
    section: ET.Element,
    incoming: ET.Element,
    source: PathLike,
    seen: dict[tuple[str, str], str],
) -> None:
    """Raise if *incoming* redefines an atom type or residue template already merged in."""
    for child in incoming:
        if child.tag not in ("Type", "Residue"):
            continue
        name = child.get("name")
        if name is None:
            continue
        key = (child.tag, name)
        if key not in seen:
            seen[key] = str(source)
            continue
        previous = seen[key]
        existing = next(
            (e for e in section if e.tag == child.tag and e.get("name") == name),
            None,
        )
        if child.tag == "Type" and existing is not None and existing.attrib == child.attrib:
            continue  # Identical redefinition: harmless, keep the one already there.
        fix = (
            "give one of them a different residue name"
            if child.tag == "Residue"
            else "regenerate one of them, since forcefill's backends name their atom types uniquely"
        )
        raise ValueError(
            f"Cannot merge {source}: it defines {child.tag.lower()} {name!r}, "
            f"which {previous} already defines differently. Two force fields "
            f"cannot share a name for different things - {fix}, or load the "
            "files separately with ForceField(...) instead of merging them."
        )


def _element_key(element: ET.Element) -> tuple[str, tuple[tuple[str, str], ...]] | None:
    """Identity of a leaf element for duplicate detection; None for one with children."""
    if len(element):
        return None
    return element.tag, tuple(sorted(element.attrib.items()))


def _extend_section(section: ET.Element, incoming: ET.Element) -> None:
    """Append *incoming*'s children to *section*, dropping exact duplicates of leaf elements.

    Every ffxml this merges declares ``<UseAttributeFromResidue name="charge"/>``
    and OpenMM rejects a second copy outright: it removes the named attribute
    from the expected list and then cannot find it again. An identical leaf
    carries nothing the first copy did not, so dropping it is safe generally.
    """
    present = {key for child in section if (key := _element_key(child)) is not None}
    for child in incoming:
        key = _element_key(child)
        if key is not None:
            if key in present:
                continue
            present.add(key)
        section.append(child)


def merge_ffxml(xml_files: Sequence[PathLike], output_xml: PathLike) -> str:
    """Merge finished OpenMM force-field XML documents into one file.

    Sections of the same kind are folded together when their attributes agree,
    and kept as separate siblings when they do not - which is exactly how OpenMM
    would see them as separate files. That distinction matters: GAFF and
    SMIRNOFF impropers use different ``ordering`` conventions, and OpenMM reads
    ``ordering`` per ``<Improper>`` at parse time, so keeping the two
    ``<PeriodicTorsionForce>`` sections apart is what preserves both.

    Args:
        xml_files: Documents to merge, in order. The first to define a name wins.
        output_xml: Where to write the merged document.

    Returns:
        The path written, as a string.

    Raises:
        ValueError: Two documents define the same atom type or residue template
            differently, or a file is not an OpenMM force-field XML.
    """
    xml_files = list(xml_files)
    if not xml_files:
        raise ValueError("merge_ffxml needs at least one XML file.")

    root = ET.Element("ForceField")
    sections: list[ET.Element] = []
    seen: dict[tuple[str, str], str] = {}
    for xml_file in xml_files:
        source_root = ET.parse(str(xml_file)).getroot()
        if source_root.tag != "ForceField":
            raise ValueError(f"{xml_file} is not an OpenMM force-field XML (root element is <{source_root.tag}>).")
        for incoming in source_root:
            section = next(
                (s for s in sections if s.tag == incoming.tag and _attributes_compatible(s.attrib, incoming.attrib)),
                None,
            )
            if section is None:
                section = ET.SubElement(root, incoming.tag, dict(incoming.attrib))
                sections.append(section)
            _check_no_redefinition(section, incoming, xml_file, seen)
            _extend_section(section, incoming)

    output_xml = Path(output_xml)
    if output_xml.parent != Path(""):
        output_xml.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root)
    ET.ElementTree(root).write(str(output_xml), encoding="utf-8", xml_declaration=True)
    log.info("Merged %d force-field XML files into %s", len(xml_files), output_xml)
    return str(output_xml)
