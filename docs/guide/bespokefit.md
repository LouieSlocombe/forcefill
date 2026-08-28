# Bespoke torsions from BespokeFit

[BespokeFit](https://github.com/openforcefield/openff-bespokefit) fits torsion
parameters to quantum-chemical torsion drives for one specific molecule and
writes an OFFXML — a stock OpenFF release with the bespoke torsions layered on
top. That file is not something OpenMM can load. Hand it to the `smirnoff`
backend and forcefill does the last mile:

```python
build_ligand_xml(
    {"BEN": LigandSpec(file="ben.sdf", forcefield="ben_bespoke.offxml")},
    "ben.xml",
    backend="smirnoff",
)
```

`forcefield` takes an installed release name (the default, `openff-2.2.1`), a
path to an OFFXML, or a list of either layered left to right — the last is for a
file carrying only the bespoke torsions, stacked on the release they were fitted
against. It is a per-ligand setting, so one ligand in a series can have a bespoke
force field while the rest use the release.

forcefill neither runs nor requires BespokeFit; producing the input is a separate
job needing psi4 or xtb, torsiondrive and ForceBalance, and hours of compute:

```
openff-bespoke executor run --file ben.sdf --output-force-field ben_bespoke.offxml
```

What forcefill adds is the checking. Three mistakes here cost you the quantum
chemistry you already paid for, and none of them report themselves:

| Mistake | What it costs |
|---|---|
| A **constrained** OFFXML — `openff-2.2.1.offxml` rather than `openff_unconstrained-2.2.1.offxml` | openmmforcefields copies every constraint into the residue template, so the ligand's X–H bonds stay rigid whatever `constraints=` you pass `createSystem`. OpenMM says nothing: a template is entitled to declare constraints ([openmmforcefields#428](https://github.com/openmm/openmmforcefields/issues/428)). BespokeFit starts from the unconstrained build by default, so this catches the override |
| The OFFXML **fitted for a different molecule** | Bespoke parameters match by SMIRKS alone, so the wrong file simply falls back to the stock parameters underneath. You get plain Sage and no warning. forcefill compares what the file assigns against what the release would and says so (`strict=False` downgrades it to a warning) |
| An OFFXML whose **1-4 scaling** differs from the base force field | The same clash as CHARMM-vs-Amber (see {doc}`charmm`), and refused the same way — but measured from the file rather than assumed, since only a stock release is guaranteed to say 0.8333/0.5 |

The "wrong molecule" check compares every handler's assignments, not just the
torsions, so an OFFXML that customizes charges, vdW, bonded terms or virtual
sites is recognized as contributing something. And a selection that is entirely
*released* chemistry is not examined at all — a release named by its path
(`openff_unconstrained-2.3.0.offxml`), which is the only way to reach one newer
than your openmmforcefields has an alias for, is stock chemistry with no bespoke
parameters to look for. The constraint and 1-4 checks still run on it.

All three are refused before the first ligand is read, so a typo in the fifth
ligand's path does not cost the AM1-BCC charges of the first four.

[`examples/parameterize_ligand_bespoke.py`](https://github.com/LouieSlocombe/forcefill/blob/main/examples/parameterize_ligand_bespoke.py) runs the whole thing — including the
three refusals — against a stand-in OFFXML synthesized from the installed
release, so it needs no BespokeFit and no QC.

