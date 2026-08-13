"""Shared in-memory structure builders for the tests (no AmberTools required)."""

from openmm import Vec3, app, unit
from openmm.app import element

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


def write_water_pdb(path):
    """Write a single TIP3P-matchable water (residue HOH with O/H1/H2)."""
    top = app.Topology()
    chain = top.addChain("W")
    res = top.addResidue("HOH", chain)
    o = top.addAtom("O", element.oxygen, res)
    h1 = top.addAtom("H1", element.hydrogen, res)
    h2 = top.addAtom("H2", element.hydrogen, res)
    top.addBond(o, h1)
    top.addBond(o, h2)
    return _write_pdb(path, top, [(0.0, 0.0, 0.0), (0.96, 0.0, 0.0), (-0.24, 0.93, 0.0)])
