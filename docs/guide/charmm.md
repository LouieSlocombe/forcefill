# CHARMM and CGenFF

```python
from forcefill import build_forcefield_xml, CHARMM_BASE_FORCEFIELD, LigandSpec

build_forcefield_xml(
    "complex.pdb",
    "extras.xml",
    base_forcefield=CHARMM_BASE_FORCEFIELD,  # charmm36.xml + charmm36/water.xml
    backend="charmm",
    ligands={"LIG": LigandSpec(charmm_files=["lig.str"])},
)
```

**forcefill converts CGenFF parameters; it cannot derive them.** There is no
CHARMM equivalent of antechamber to call: parameters come from
[ParamChem](https://cgenff.paramchem.org) or the licensed `cgenff` program, both
of which emit a CHARMM stream file. Give forcefill that file and it produces a
validated ffxml — the part that is fiddly enough to get quietly wrong.

**CHARMM is not interchangeable with Amber.** Amber scales 1-4 interactions by
0.8333/0.5 and CHARMM by 1.0/1.0, and OpenMM refuses to load force fields that
disagree:

```
ValueError: Found multiple NonbondedForce tags with different 1-4 scales
```

So a CHARMM ligand needs `base_forcefield=CHARMM_BASE_FORCEFIELD`, and cannot
share a build with a `gaff` or `smirnoff` one. Both mistakes are refused up
front, by name, rather than left to OpenMM at the end of the run.

The generated XML is a **residue template**, not a self-contained force field:
`charmm36.xml` already carries all 412 CGenFF atom types and their parameters,
so forcefill names them rather than redefining them, and adds only the terms the
stream file itself supplies (what ParamChem assigned by analogy). So no CHARMM
toppar download is needed — a ParamChem `.str` and OpenMM's own `charmm36.xml`
are enough — but loading the result *without* `charmm36.xml` underneath it will
not work, since it refers to atom types it does not define.

Three ways this goes silently wrong if done by hand with ParmEd, all of which
forcefill handles — they are why the backend is a module and not a three-line
script:

| What ParmEd does | What it costs |
|---|---|
| Writes `sigma="1.0" epsilon="0.0"` for atom types it only knows the mass of | Loaded after `charmm36.xml`, those **override the real Lennard-Jones parameters**. OpenMM treats a repeated atom type as an override, not a clash, so there is no error |
| Defaults to `separate_ljforce=False` | The ligand's LJ energy is **counted twice** — once in `NonbondedForce`, once in the `CustomNonbondedForce` `charmm36.xml` builds for its NBFIX pairs |
| Drops a residue template whose atom types it cannot resolve, with a `ParameterWarning` | A ParamChem `.str` never carries `MASS` records, so the default outcome is an **empty ffxml that loads fine and parameterizes nothing** |

Two smaller things follow from `charmm36.xml` shipping 814 residue templates —
every amino acid, nucleotide, lipid and CGenFF model compound:

* a fragment-sized ligand may be matched by the base force field already, in
  which case there is nothing to parameterize and forcefill says so;
* a ligand named after one of them (`MET`, `PC`, `CA`) is refused, because
  OpenMM will not load two templates with the same name.

Finally, `build_ligand_xml` can validate a CHARMM ligand but not minimize it: a
stream file records internal coordinates, not Cartesian ones, so there is no
geometry to minimize. Use `build_forcefield_xml`, where the coordinates come
from the structure.

