"""Read a ligand file well enough to check it before the expensive step runs.

:func:`inspect_ligand_file` reads elements, bonds, coordinates and formal charge
out of SDF (V2000 and V3000), MOL2 and PDB; :func:`check_matches_residue` and
:func:`check_geometry` are the checks built on it, applied to a whole run by
:mod:`forcefill.preflight`.

RDKit does the reading where it can, and small text readers cover what it cannot
- chiefly the GAFF-typed mol2 antechamber writes, whose atom types are not SYBYL
ones. ``formal_charge`` is the one field the text readers cannot always recover,
which is why ``prefer_rdkit`` exists and defaults to True.
"""

from __future__ import annotations

import logging
import math
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import parmed
from openmm import app, unit
from parmed.modeller import ResidueTemplateContainer
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from ._spec import PathLike

log = logging.getLogger(__name__)

__all__ = [
    "LigandFileInfo",
    "check_geometry",
    "check_matches_residue",
    "inspect_ligand_file",
    "residue_formula",
    "residue_name_for",
    "smiles_to_sdf",
    "smiles_with_residue_geometry",
    "split_multi_sdf",
]

#: Suffixes :func:`inspect_ligand_file` knows how to read.
_READABLE = {".sdf", ".sd", ".mol", ".mol2", ".pdb"}

#: Closer than this (angstrom) and two atoms are on top of each other: the
#: standard cause of an infinite Coulomb energy that only shows up as a NaN.
MIN_ATOM_SEPARATION = 0.5

#: Legacy SDF V2000 atom-block charge codes (column 37-39). Superseded by
#: ``M  CHG`` records, which win whenever any are present.
_V2000_CHARGE_CODES = {1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3}


@dataclass
class LigandFileInfo:
    """What one ligand file says about the molecule in it."""

    #: The file this was read from.
    path: str
    #: Element symbols with their counts. ``"?"`` for an atom whose element
    #: could not be determined.
    elements: Counter[str] = field(default_factory=Counter)
    #: Bond count, or None when the format does not carry one.
    n_bonds: int | None = None
    #: Formal net charge, or None when it could not be determined. MOL2 carries
    #: no formal charges, so this is the rounded sum of its partial charges and
    #: is None when those are all zero.
    formal_charge: int | None = None
    #: Coordinates in angstrom, in file order.
    positions: list[tuple[float, float, float]] = field(default_factory=list)
    #: How it was read: ``"rdkit"`` or ``"text"``.
    source: str = "text"

    @property
    def n_atoms(self) -> int:
        """Total atom count."""
        return sum(self.elements.values())

    @property
    def formula(self) -> str:
        """Molecular formula in Hill notation (C, then H, then alphabetical)."""
        return _hill_formula(self.elements)


def _hill_formula(elements: Counter[str]) -> str:
    """Format an element count as ``C7H8N2``: carbon, hydrogen, then the rest alphabetically."""
    rest = sorted(k for k in elements if k not in ("C", "H"))
    order = [k for k in ("C", "H") if elements[k]] + rest
    return "".join(f"{k}{elements[k]}" if elements[k] > 1 else k for k in order)


def residue_formula(residue: app.topology.Residue) -> Counter[str]:
    """Element counts for an OpenMM residue; atoms with no element count as ``"?"``."""
    return Counter(a.element.symbol if a.element is not None else "?" for a in residue.atoms())


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------


def _rdkit_mol_from_file(path: Path) -> Chem.Mol | None:
    """Load *path* with RDKit, or return None if RDKit cannot read it.

    None means the format defeated RDKit - chiefly the GAFF-typed mol2
    antechamber writes, whose atom types are not SYBYL ones - not that RDKit is
    missing; it is a hard dependency. The text readers below cover those.
    """
    # RDKit narrates every valence complaint to stderr; the caller gets a proper
    # message from the checks below instead.
    RDLogger.DisableLog("rdApp.*")
    try:
        suffix = path.suffix.lower()
        if suffix in (".sdf", ".sd", ".mol"):
            supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
            mols = [m for m in supplier if m is not None]
            if len(mols) > 1:
                log.warning(
                    "%s holds %d molecules; reading the first. Use split_multi_sdf() to parameterize them all.",
                    path,
                    len(mols),
                )
            return mols[0] if mols else None
        if suffix == ".mol2":
            return Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
        if suffix == ".pdb":
            return Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
    finally:
        RDLogger.EnableLog("rdApp.*")
    return None


def _info_from_rdkit(mol: Chem.Mol, path: Path) -> LigandFileInfo:
    """Build a LigandFileInfo from an RDKit molecule."""
    positions: list[tuple[float, float, float]] = []
    if mol.GetNumConformers():
        conf = mol.GetConformer()
        positions = [(p.x, p.y, p.z) for p in (conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms()))]
    return LigandFileInfo(
        path=str(path),
        elements=Counter(a.GetSymbol() for a in mol.GetAtoms()),
        n_bonds=mol.GetNumBonds(),
        # A PDB carries no bond orders, so RDKit's formal charge for one is 0 by
        # construction rather than by evidence. Reporting that as a fact is how
        # a charged ligand silently gets parameterized as neutral - exactly what
        # the inference is meant to prevent.
        formal_charge=None if path.suffix.lower() == ".pdb" else Chem.GetFormalCharge(mol),
        positions=positions,
        source="rdkit",
    )


def _read_sdf_v3000(lines: list[str], path: Path) -> LigandFileInfo:
    """Parse the V3000 connection table: ``M  V30`` records rather than fixed columns."""
    elements: Counter[str] = Counter()
    positions: list[tuple[float, float, float]] = []
    charge = 0
    n_bonds = 0
    section = ""
    for raw in lines:
        if not raw.startswith("M  V30 "):
            continue
        body = raw[7:].strip()
        upper = body.upper()
        if upper.startswith("BEGIN "):
            section = upper[6:].strip()
            continue
        if upper.startswith("END "):
            section = ""
            continue
        if section == "ATOM":
            fields = body.split()
            # index symbol x y z aamap [KEY=value ...]
            if len(fields) < 6:
                continue
            elements[fields[1]] += 1
            positions.append((float(fields[2]), float(fields[3]), float(fields[4])))
            for extra in fields[6:]:
                if extra.upper().startswith("CHG="):
                    charge += int(extra[4:])
        elif section == "BOND":
            n_bonds += 1
    return LigandFileInfo(path=str(path), elements=elements, n_bonds=n_bonds, formal_charge=charge, positions=positions)


def _read_sdf(path: Path) -> LigandFileInfo:
    """Parse the first molecule of an SDF/MOL file without RDKit."""
    lines = path.read_text().splitlines()
    if len(lines) < 4:
        raise ValueError(f"{path} is too short to be an SDF/MOL file.")
    counts = lines[3]
    if "V3000" in counts.upper():
        return _read_sdf_v3000(lines, path)

    try:
        n_atoms = int(counts[0:3])
        n_bonds = int(counts[3:6])
    except ValueError as exc:
        raise ValueError(f"{path} line 4 is not an SDF counts line: {counts!r}") from exc

    elements: Counter[str] = Counter()
    positions: list[tuple[float, float, float]] = []
    legacy_charge = 0
    for line in lines[4 : 4 + n_atoms]:
        elements[line[31:34].strip()] += 1
        positions.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
        code = line[36:39].strip()
        if code:
            legacy_charge += _V2000_CHARGE_CODES.get(int(code), 0)

    # M CHG records supersede the atom-block codes entirely: the format says a
    # file carrying any of them has zeroed the legacy fields.
    m_charge = 0
    saw_m_charge = False
    for line in lines[4 + n_atoms + n_bonds :]:
        if not line.startswith("M  CHG"):
            continue
        saw_m_charge = True
        # M  CHG  <count>  <atom> <charge> [<atom> <charge> ...]
        pairs = line.split()[3:]
        m_charge += sum(int(pairs[i + 1]) for i in range(0, len(pairs) - 1, 2))
    formal_charge = m_charge if saw_m_charge else legacy_charge
    return LigandFileInfo(
        path=str(path), elements=elements, n_bonds=n_bonds, formal_charge=formal_charge, positions=positions
    )


def _read_mol2(path: Path) -> LigandFileInfo:
    """Read a MOL2 with ParmEd, which knows the atom types RDKit's parser rejects.

    Hand-parsing the type column does not work: antechamber writes GAFF types
    (``c3``, ``ho``, ``oh``), not SYBYL ones (``C.3``), and they are ambiguous
    with element symbols - GAFF's ``ca`` is an aromatic carbon, not calcium.

    MOL2 has no formal-charge field, so the net charge is the rounded sum of the
    partial charges, reported only when they are not all zero and the sum is
    close enough to an integer to be real.
    """
    template = parmed.load_file(str(path))
    if isinstance(template, ResidueTemplateContainer):
        if len(template) > 1:
            log.warning(
                "%s holds %d residues; reading the first. Use split_multi_sdf() or "
                "separate files to parameterize them all.",
                path,
                len(template),
            )
        template = template[0]

    atoms = list(template.atoms)
    elements = Counter(
        app.element.Element.getByAtomicNumber(a.atomic_number).symbol if a.atomic_number > 0 else "?" for a in atoms
    )
    positions = [(a.xx, a.xy, a.xz) for a in atoms] if all(hasattr(a, "xx") for a in atoms) else []

    partial_sum = sum(a.charge for a in atoms)
    formal_charge = None
    if any(a.charge for a in atoms):
        rounded = round(partial_sum)
        if abs(partial_sum - rounded) < 0.05:
            formal_charge = int(rounded)
        else:
            log.warning(
                "%s: the partial charges sum to %+.3f, which is not close to an "
                "integer; not inferring a net charge from it.",
                path,
                partial_sum,
            )
    return LigandFileInfo(
        path=str(path),
        elements=elements,
        n_bonds=len(template.bonds),
        formal_charge=formal_charge,
        positions=positions,
    )


def _read_pdb(path: Path) -> LigandFileInfo:
    """Read a ligand PDB with OpenMM's own parser; formal charge is not recoverable from PDB."""
    pdb = app.PDBFile(str(path))
    top = pdb.topology
    positions = [tuple(v.value_in_unit(unit.angstrom)) for v in pdb.positions]
    return LigandFileInfo(
        path=str(path),
        elements=Counter(a.element.symbol if a.element is not None else "?" for a in top.atoms()),
        n_bonds=top.getNumBonds(),
        formal_charge=None,
        positions=positions,
    )


def inspect_ligand_file(path: PathLike, *, prefer_rdkit: bool = True) -> LigandFileInfo:
    """Read elements, bonds, coordinates and formal charge out of a ligand file.

    Uses RDKit when it is installed - it perceives formal charges reliably and
    handles the format corners - and falls back to a small text parser for SDF
    (V2000/V3000) and MOL2, and to OpenMM's reader for PDB. Only the first
    molecule of a multi-molecule file is read; see :func:`split_multi_sdf`.

    Args:
        path: The ligand file. Suffix decides the format.
        prefer_rdkit: Set False to force the text parser, which is how the
            fallback gets exercised in tests.

    Returns:
        LigandFileInfo

    Raises:
        FileNotFoundError: The file does not exist.
        ValueError: The suffix is not one this can read, or the contents do not
            parse as that format.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Ligand file not found: {path}")
    suffix = path.suffix.lower()
    if suffix not in _READABLE:
        raise ValueError(
            f"Cannot inspect {path.name!r}: unknown ligand file type {suffix!r} (known: {sorted(_READABLE)})."
        )

    if prefer_rdkit:
        mol = _rdkit_mol_from_file(path)
        if mol is not None:
            return _info_from_rdkit(mol, path)
        log.debug("RDKit could not read %s; falling back to the text parser.", path)

    if suffix == ".mol2":
        return _read_mol2(path)
    if suffix == ".pdb":
        return _read_pdb(path)
    return _read_sdf(path)


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def _report(message: str, strict: bool) -> None:
    """Raise when *strict*, warn otherwise - the shape every check below shares."""
    if strict:
        raise ValueError(message)
    log.warning("%s (strict=False, continuing anyway)", message)


def check_matches_residue(
    info: LigandFileInfo,
    residue: app.topology.Residue,
    name: str,
    *,
    strict: bool = True,
) -> None:
    """Check that a supplied ligand file is the same molecule as the PDB residue.

    The generated template is matched against the *PDB's* bond graph, so a file
    differing by even one hydrogen produces a template that cannot match - and
    only after antechamber has run. This says so in a second, and names the
    difference.

    Args:
        info: What :func:`inspect_ligand_file` read from the supplied file.
        residue: The residue as it appears in the input structure.
        name: Residue name, for the message.
        strict: Raise on a mismatch (default) rather than warn.

    Raises:
        ValueError: The compositions differ and *strict* is True.
    """
    expected = residue_formula(residue)
    if expected == info.elements:
        return
    if "?" in expected:
        log.warning(
            "Residue %s has %d atom(s) with no element assigned, so it cannot be "
            "compared with %s. Fill in the element columns (77-78) of the PDB.",
            name,
            expected["?"],
            info.path,
        )
        return

    differences = [
        f"{symbol}: file {info.elements.get(symbol, 0)} vs PDB {expected.get(symbol, 0)}"
        for symbol in sorted(set(expected) | set(info.elements))
        if expected.get(symbol, 0) != info.elements.get(symbol, 0)
    ]
    _report(
        f"The ligand file for {name} is not the same molecule as the residue in "
        f"the structure: {Path(info.path).name} has {info.formula} "
        f"({info.n_atoms} atoms), residue {name} has {_hill_formula(expected)} "
        f"({sum(expected.values())} atoms). Differing elements - "
        f"{'; '.join(differences)}. The generated template is matched against "
        "the structure's bond graph, so this cannot produce a working force "
        "field. Supply the file the structure was built from, including its "
        "hydrogens and protonation state.",
        strict,
    )


def check_geometry(
    positions: list[tuple[float, float, float]],
    name: str,
    *,
    strict: bool = True,
    min_separation: float = MIN_ATOM_SEPARATION,
) -> None:
    """Check *positions* (angstrom) for the geometry faults that produce NaN energies.

    Catches coincident atoms, non-finite coordinates and an all-zero conformer
    (a molecule written without 3D coordinates) - the standard causes of the
    infinite energies ``minimize=True`` otherwise only finds at the very end.

    Args:
        positions: Coordinates in angstrom.
        name: Residue name, for the message.
        strict: Raise (default) rather than warn.
        min_separation: Atoms closer than this are considered coincident.

    Raises:
        ValueError: A fault was found and *strict* is True.
    """
    if not positions:
        return
    bad = [i for i, p in enumerate(positions) if not all(math.isfinite(c) for c in p)]
    if bad:
        _report(f"Residue {name} has non-finite coordinates for atom index(es) {bad[:10]}.", strict)
        return
    if len(positions) > 1 and all(p == (0.0, 0.0, 0.0) for p in positions):
        _report(
            f"Every atom of residue {name} sits at the origin: the file carries "
            "no 3D conformer. Generate coordinates before parameterizing - "
            "AM1-BCC charges are derived from the geometry.",
            strict,
        )
        return
    for i, j, distance in _close_pairs(positions, min_separation):
        _report(
            f"Atoms {i} and {j} of residue {name} are {distance:.3f} A apart, "
            f"closer than {min_separation} A. Coincident atoms give an infinite "
            "Coulomb energy; check the structure for duplicated or unresolved "
            "atoms (a leftover altLoc is the usual cause).",
            strict,
        )
        return


def _close_pairs(
    positions: list[tuple[float, float, float]], min_separation: float
) -> Iterator[tuple[int, int, float]]:
    """Yield ``(i, j, distance)`` for every atom pair closer than *min_separation* angstrom."""
    cutoff = min_separation**2
    for i, a in enumerate(positions):
        for j in range(i + 1, len(positions)):
            b = positions[j]
            squared = (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
            if squared < cutoff:
                yield i, j, math.sqrt(squared)


# --------------------------------------------------------------------------
# Conversion
# --------------------------------------------------------------------------


def _write_sdf(mol: Chem.Mol, out_sdf: PathLike) -> str:
    """Write one molecule to *out_sdf*, creating its parent directory."""
    out_sdf = Path(out_sdf)
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(out_sdf))
    try:
        writer.write(mol)
    finally:
        writer.close()
    return str(out_sdf)


def smiles_to_sdf(smiles: str, out_sdf: PathLike, name: str = "LIG", *, random_seed: int = 0xF0) -> str:
    """Embed *smiles* as a 3D SDF with explicit hydrogens and return the path.

    Add hydrogens, embed a conformer, relax it with MMFF - so the geometry the
    charges are derived from is reasonable, not merely valid. The seed is fixed
    so a re-run reproduces the same conformer and therefore the same charges.

    Raises:
        RuntimeError: RDKit could not parse the SMILES or embed a conformer.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse the SMILES for {name}: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError(
            f"RDKit could not generate a 3D conformer for {name} from {smiles!r}. "
            "Supply a prepared SDF/MOL2 with LigandSpec(file=...) instead."
        )
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("_Name", name)

    out_sdf = _write_sdf(mol, out_sdf)
    log.info("Embedded %s from SMILES %r: %s", name, smiles, out_sdf)
    return out_sdf


def smiles_with_residue_geometry(
    smiles: str,
    residue_pdb: PathLike,
    out_sdf: PathLike,
    name: str = "LIG",
) -> str:
    """Apply *smiles*' bond orders to the geometry in *residue_pdb* and write an SDF.

    For a ligand present in the structure whose bonds the structure cannot state:
    coordinates and atom count come from the PDB, bond orders and formal charges
    from the SMILES. The two must describe the same molecule, hydrogens included
    - ``AssignBondOrdersFromTemplate`` enforces that.

    Raises:
        RuntimeError: RDKit could not read either input, or the SMILES does not
            describe the same molecule as the residue.
    """
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise RuntimeError(f"RDKit could not parse the SMILES for {name}: {smiles!r}")
    template = Chem.AddHs(template)
    extracted = Chem.MolFromPDBFile(str(residue_pdb), removeHs=False, sanitize=False)
    if extracted is None:
        raise RuntimeError(f"RDKit could not read the extracted residue {name} from {residue_pdb}.")
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, extracted)
    except ValueError as exc:
        raise RuntimeError(
            f"The SMILES for {name} does not match the residue in the structure: {exc}\n"
            f"  SMILES: {smiles}\n"
            "Both must be the same molecule with the same hydrogens. Check the "
            "protonation state - it is the usual difference."
        ) from exc
    mol.SetProp("_Name", name)

    out_sdf = _write_sdf(mol, out_sdf)
    log.info("Applied the SMILES bond orders for %s to the structure's geometry: %s", name, out_sdf)
    return out_sdf


def split_multi_sdf(path: PathLike, outdir: PathLike) -> dict[str, str]:
    """Split a multi-molecule SDF into one file per molecule, keyed by residue name.

    Names come from each record's title line, uppercased and cleaned; records
    with no usable title are named ``<stem>1``, ``<stem>2``, ... A duplicate name
    raises rather than overwriting - two molecules cannot share a template.

    Returns:
        ``{residue_name: sdf_path}``, in file order.

    Raises:
        ValueError: The file holds no molecules, or two records claim one name.
    """
    path = Path(path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = _sdf_records(path)
    if not records:
        raise ValueError(f"{path} contains no molecules.")

    out: dict[str, str] = {}
    for index, record in enumerate(records, start=1):
        # record[0] is the title line, often blank - hence the line-exact split
        # above. Trimming blank lines instead would promote line 2, the program
        # line ("     RDKit          3D"), into the residue name.
        title = record[0].strip()
        name = _residue_name_from(title) or _residue_name_from(f"{path.stem}{index}") or f"L{index:02d}"
        if name in out:
            raise ValueError(
                f"{path} record {index} is named {name!r}, which record "
                f"{list(out).index(name) + 1} already used. Residue names must be "
                "unique; rename the titles or split the file yourself."
            )
        out[name] = str(outdir / f"{name}.sdf")
        Path(out[name]).write_text("\n".join([*record, "$$$$", ""]))
    log.info("Split %s into %d molecules: %s", path, len(out), sorted(out))
    return out


def _sdf_records(path: Path) -> list[list[str]]:
    """Split an SDF into per-molecule line lists on the ``$$$$`` terminator."""
    records: list[list[str]] = []
    current: list[str] = []
    for line in path.read_text().splitlines():
        if line.strip() == "$$$$":
            if any(entry.strip() for entry in current):
                records.append(current)
            current = []
            continue
        current.append(line)
    # A final molecule with no terminator is malformed but common enough to accept.
    if any(entry.strip() for entry in current):
        records.append(current)
    return records


def _residue_name_from(text: str) -> str:
    """Turn arbitrary text into a PDB-style residue name: alphanumeric, uppercase, <= 3 characters."""
    cleaned = "".join(c for c in text if c.isalnum()).upper()
    if cleaned and cleaned[0].isdigit():
        cleaned = f"L{cleaned}"
    return cleaned[:3]


def residue_name_for(path: PathLike) -> str:
    """Derive a residue name from a ligand file's stem (``benzamidinium.sdf`` -> ``BEN``).

    Only the file name is consulted - it is what the caller can see and change.
    Pass an explicit mapping to :func:`~forcefill.build_ligand_xml` when the
    derived name is not the one you want.
    """
    stem = Path(path).stem
    name = _residue_name_from(stem)
    if not name:
        raise ValueError(
            f"Cannot derive a residue name from {Path(path).name!r} - it has no "
            "alphanumeric characters. Pass a name explicitly, as "
            "build_ligand_xml({'LIG': '<file>'}, ...)."
        )
    if name != stem.upper():
        log.info("Using residue name %s for %s.", name, path)
    return name
