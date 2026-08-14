"""Shared in-memory structure builders for the tests (no AmberTools required)."""

from pathlib import Path

from openmm import Vec3, app, unit
from openmm.app import element

#: Committed fixtures: byte-for-byte AmberTools output the assembly tests read.
DATA = Path(__file__).parent / "data"

#: Idealized methanol geometry, matching tests/data/methanol.mol2 (angstrom).
METHANOL_XYZ = [
    (0.000, 0.000, 0.000),
    (1.410, 0.000, 0.000),
    (-0.360, 1.030, 0.000),
    (-0.360, -0.515, 0.892),
    (-0.360, -0.515, -0.892),
    (1.730, 0.890, 0.000),
]

METHANOL_ATOMS = [
    ("C1", element.carbon),
    ("O1", element.oxygen),
    ("H1", element.hydrogen),
    ("H2", element.hydrogen),
    ("H3", element.hydrogen),
    ("H4", element.hydrogen),
]

METHANOL_BONDS = [("C1", "O1"), ("C1", "H1"), ("C1", "H2"), ("C1", "H3"), ("O1", "H4")]


def add_methanol_residue(top, chain_id="A", name="LIG"):
    """Append a methanol residue to *top*; returns its coordinate list (angstrom tuples)."""
    chain = top.addChain(chain_id)
    res = top.addResidue(name, chain)
    atoms = {}
    for atom_name, elem in METHANOL_ATOMS:
        atoms[atom_name] = top.addAtom(atom_name, elem, res)
    for a, b in METHANOL_BONDS:
        top.addBond(atoms[a], atoms[b])
    return list(METHANOL_XYZ)


def methanol_residue(name="LIG"):
    """A lone methanol residue in a topology of its own, reachable as ``residue.chain.topology``."""
    top = app.Topology()
    add_methanol_residue(top, name=name)
    return next(top.residues())


def methanol_positions(xyz=METHANOL_XYZ):
    """Methanol's coordinates as an OpenMM Quantity, in angstrom."""
    return unit.Quantity([Vec3(*p) for p in xyz], unit.angstrom)


def add_broken_gly_residue(top, chain_id="B"):
    """Append a hydrogen-stripped free glycine (standard name, missing atoms -> gets skipped)."""
    chain = top.addChain(chain_id)
    gly = top.addResidue("GLY", chain)
    atoms = {}
    for name, elem in [
        ("N", element.nitrogen),
        ("CA", element.carbon),
        ("C", element.carbon),
        ("O", element.oxygen),
    ]:
        atoms[name] = top.addAtom(name, elem, gly)
    top.addBond(atoms["N"], atoms["CA"])
    top.addBond(atoms["CA"], atoms["C"])
    top.addBond(atoms["C"], atoms["O"])
    return [
        (8.000, 8.000, 8.000),
        (9.450, 8.000, 8.000),
        (10.050, 9.350, 8.000),
        (11.280, 9.400, 8.000),
    ]


#: Idealized glycerol geometry (angstrom): the archetypal crystallization
#: additive, and the residue name most often mistaken for a ligand.
GLYCEROL_ATOMS = [
    ("C1", element.carbon),
    ("C2", element.carbon),
    ("C3", element.carbon),
    ("O1", element.oxygen),
    ("O2", element.oxygen),
    ("O3", element.oxygen),
]

GLYCEROL_XYZ = [
    (0.000, 0.000, 0.000),
    (1.520, 0.000, 0.000),
    (2.100, 1.400, 0.000),
    (-0.480, 1.330, 0.000),
    (2.000, -0.700, 1.150),
    (3.520, 1.400, 0.000),
]

GLYCEROL_BONDS = [("C1", "C2"), ("C2", "C3"), ("C1", "O1"), ("C2", "O2"), ("C3", "O3")]


def add_water_residue(top, chain_id="W", name="HOH", origin=(0.0, 0.0, 0.0), virtual_site=False):
    """Append one water to *top*; returns its coordinate list (angstrom tuples).

    With *virtual_site*, adds the massless ``M`` particle of a 4-site model -
    element None, which is why the cleaner classifies water before it looks at
    elements.
    """
    chain = top.addChain(chain_id)
    res = top.addResidue(name, chain)
    o = top.addAtom("O", element.oxygen, res)
    h1 = top.addAtom("H1", element.hydrogen, res)
    h2 = top.addAtom("H2", element.hydrogen, res)
    top.addBond(o, h1)
    top.addBond(o, h2)
    dx, dy, dz = origin
    xyz = [(dx, dy, dz), (dx + 0.96, dy, dz), (dx - 0.24, dy + 0.93, dz)]
    if virtual_site:
        top.addAtom("M", None, res)
        xyz.append((dx + 0.15, dy + 0.06, dz))
    return xyz


def add_ion_residue(top, name, elem, chain_id="I", origin=(0.0, 0.0, 0.0)):
    """Append a one-atom ion residue (name and atom name both *name*)."""
    chain = top.addChain(chain_id)
    res = top.addResidue(name, chain)
    top.addAtom(name, elem, res)
    return [origin]


def add_glycerol_residue(top, chain_id="G", name="GOL", origin=(0.0, 0.0, 0.0)):
    """Append a glycerol (heavy atoms only, as X-ray additives come); returns its coordinates."""
    chain = top.addChain(chain_id)
    res = top.addResidue(name, chain)
    atoms = {}
    for atom_name, elem in GLYCEROL_ATOMS:
        atoms[atom_name] = top.addAtom(atom_name, elem, res)
    for a, b in GLYCEROL_BONDS:
        top.addBond(atoms[a], atoms[b])
    dx, dy, dz = origin
    return [(x + dx, y + dy, z + dz) for x, y, z in GLYCEROL_XYZ]


def bond_across_residues(top, res_name_a, res_name_b):
    """Bond the first atom of the first *res_name_a* to the first atom of the first *res_name_b*.

    Makes an otherwise free-standing residue covalently linked, which is the
    one condition under which the cleaner refuses to delete it.
    """
    by_name = {}
    for res in top.residues():
        by_name.setdefault(res.name, res)
    top.addBond(next(by_name[res_name_a].atoms()), next(by_name[res_name_b].atoms()))


def _write_pdb(path, top, xyz):
    positions = unit.Quantity([Vec3(*p) for p in xyz], unit.angstrom)
    with open(path, "w") as fh:
        app.PDBFile.writeFile(top, positions, fh)
    return path


def write_methanol_pdb(path, broken_gly=False):
    """Write methanol as hetero residue LIG, optionally plus a hydrogen-stripped GLY that gets skipped."""
    top = app.Topology()
    xyz = add_methanol_residue(top)
    if broken_gly:
        xyz += add_broken_gly_residue(top)
    return _write_pdb(path, top, xyz)


def write_broken_gly_pdb(path):
    """Write only the hydrogen-stripped glycine: everything unmatched gets skipped."""
    top = app.Topology()
    xyz = add_broken_gly_residue(top)
    return _write_pdb(path, top, xyz)


def write_methanol_sdf(path):
    """Write methanol as a V2000 SDF with explicit bonds, matching METHANOL_XYZ."""
    symbols = [elem.symbol for _, elem in METHANOL_ATOMS]
    lines = [
        "methanol",
        "  forcefill",
        "",
        f"{len(symbols):3d}{len(METHANOL_BONDS):3d}  0  0  0  0  0  0  0  0999 V2000",
    ]
    for (x, y, z), sym in zip(METHANOL_XYZ, symbols, strict=True):
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {sym:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    index = {name: i + 1 for i, (name, _) in enumerate(METHANOL_ATOMS)}
    for a, b in METHANOL_BONDS:
        lines.append(f"{index[a]:3d}{index[b]:3d}  1  0")
    lines += ["M  END", "$$$$", ""]
    Path(path).write_text("\n".join(lines))
    return path


def write_water_pdb(path):
    """Write a single TIP3P-matchable water (residue HOH with O/H1/H2)."""
    top = app.Topology()
    return _write_pdb(path, top, add_water_residue(top))


def write_ligand_and_water_pdb(path):
    """Write methanol plus one hydrogen-bearing water: only the water is solvent."""
    top = app.Topology()
    xyz = add_methanol_residue(top)
    xyz += add_water_residue(top, origin=(5.0, 0.0, 0.0))
    return _write_pdb(path, top, xyz)


def write_ligand_and_glycerol_pdb(path):
    """Write methanol plus a glycerol: a free-standing additive is indistinguishable from a ligand."""
    top = app.Topology()
    xyz = add_methanol_residue(top)
    xyz += add_glycerol_residue(top, origin=(0.0, -8.0, 0.0))
    return _write_pdb(path, top, xyz)


def write_dirty_pdb(path, waters=3):
    """Write what comes off the PDB: a ligand plus water, a counter-ion, a structural metal and glycerol.

    Everything but LIG is in its own chain, so deleting a category empties a
    chain and exercises Modeller dropping it.
    """
    top = app.Topology()
    xyz = add_methanol_residue(top)
    for i in range(waters):
        xyz += add_water_residue(top, chain_id=f"W{i}", origin=(5.0 + 3.0 * i, 0.0, 0.0))
    xyz += add_ion_residue(top, "NA", element.sodium, chain_id="N", origin=(-5.0, 0.0, 0.0))
    xyz += add_ion_residue(top, "CA", element.calcium, chain_id="C", origin=(-8.0, 0.0, 0.0))
    xyz += add_glycerol_residue(top, origin=(0.0, -8.0, 0.0))
    return _write_pdb(path, top, xyz)
