# Test data

Reference input the hermetic tests read, so the ParmEd assembly paths and the
`build_forcefield_xml` orchestration can be tested without antechamber
installed. Do not hand-edit or let whitespace hooks touch these files
(`.pre-commit-config.yaml` excludes this directory).

## AmberTools output (the gaff backend)

Byte-for-byte as the executables wrote them.

| File | What it is |
|---|---|
| `methanol.mol2` | antechamber output for the test methanol (residue `LIG`): GAFF2 atom types + AM1-BCC charges |
| `methanol.frcmod` | parmchk2 output with `-a Y`, i.e. a *complete* parameter set (MASS/BOND/ANGLE/DIHE/NONBON) for the mol2 — it stands in for `gaff2.dat` in tests |

Regenerate (AmberTools required; the input PDB is written by
`tests/helpers.py:write_methanol_pdb`):

```bash
antechamber -i methanol.pdb -fi pdb -o methanol.mol2 -fo mol2 \
    -c bcc -nc 0 -m 1 -at gaff2 -rn LIG -pf y
parmchk2 -i methanol.mol2 -f mol2 -o methanol.frcmod -s gaff2 -a Y
```

## CGenFF stream files (the charmm backend)

| File | What it is |
|---|---|
| `methanol_cgenff.str` | CGenFF parameters for the test methanol (residue `LIG`), with a `read param` section so the analogy-parameter path is exercised |
| `chloroethanol_cgenff.str` | The same for 2-chloroethanol (residue `CET`), used by the tests that go through a structure |

These two are **hand-written stand-ins for ParamChem output**, not real
`cgenff` runs — the program is licensed and its output is not redistributable.
They are written in the format it emits, and the atom types and the parameters
in `methanol_cgenff.str` are the real CGenFF ones, taken from
`top_all36_cgenff.rtf` / `par_all36_cgenff.prm`:

* `methanol_cgenff.str` reproduces `RESI MEOH` with the test methanol's atom
  names (`C1`/`O1`/`H1`–`H4`, matching `tests/helpers.py:METHANOL_ATOMS`);
* `chloroethanol_cgenff.str` combines the chloro group of `RESI CLET` with the
  hydroxyl group of `RESI ETOH`, with the charge on `C1` adjusted so each group
  is neutral. That charge is an analogy, not a published value.

The tests assert mechanics — that a residue template is written, that nothing
`charmm36.xml` already defines is redefined, and that its Lennard-Jones
parameters survive — so nothing depends on those charges being the ones
ParamChem would produce. Do not copy either file into real work: run the ligand
through ParamChem (<https://cgenff.paramchem.org>) instead.

Two properties are load-bearing and easy to break by editing:

* **Both stream files' atom names must match the structure helpers**, since
  those are what the residue template is checked against.
* **2-chloroethanol, not methanol, is the ligand for the structure tests.**
  `charmm36.xml` ships 814 residue templates including every CGenFF model
  compound, so methanol is matched by the base force field and never reaches a
  backend at all.
