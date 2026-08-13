# Test data

Byte-for-byte AmberTools reference output used by the hermetic tests, so the
ParmEd assembly path and the `build_forcefield_xml` orchestration can be tested
without antechamber installed. Do not hand-edit or let whitespace hooks touch
these files (`.pre-commit-config.yaml` excludes this directory).

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
