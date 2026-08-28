# Checks

forcefill checks twice: once on the input, before anything expensive runs,
and once on the output, against the force field it just produced.

## Before the expensive step

antechamber's AM1-BCC can take an hour on a drug-sized ligand. Three mistakes
that used to cost that hour — or worse, silently produce wrong numbers — are
caught in the first second:

- **Net charge read from the file.** An SDF or MOL2 states its own formal
  charge, so `net_charge` no longer defaults to a silent 0 for those. This
  applies only to a ligand supplied as a file or SMILES, never to a residue
  extracted from a PDB, and an explicit `net_charge` always wins — but one that
  contradicts the file raises rather than picking a side.
- **The ligand file must be the residue in the PDB.** A mismatch used to appear
  only at the end, as an opaque "no template matched". Now it names the
  difference: `ben.sdf has C7H9N2 (18 atoms), residue BEN has C7H8N2 (17)`.
- **Geometry sanity.** Coincident atoms, non-finite coordinates or a molecule
  written with no conformer — the standard causes of the NaN energies that
  `minimize=True` otherwise only finds at the very end.

`strict=False` downgrades the last two to warnings; these are heuristics and the
long tail is real.

## After: the parameters, not just the templates

Building a `System` says nothing about whether the numbers in it are physical: a
NaN charge or a zero force constant survives it and only shows up later as an
exploding simulation. `minimize=True` adds an energy evaluation and a short
minimization — of each parameterized residue in vacuum, and of the whole input
when nothing was skipped — and raises if the potential energy is not finite at
either end:

```python
result = build_forcefield_xml("complex.pdb", "extras.xml", minimize=True)
lig = result.minimizations["LIG"]
print(f"{lig.initial_energy:.0f} -> {lig.final_energy:.0f} kJ/mol")
print(result.full_minimization.max_force)  # kJ/mol/nm
```

`max_force` and `energy_change` are reported for inspection, not enforced —
what counts as converged depends on the system. The same check is available on
its own as `minimize_with_forcefield_xml(topology, positions, xml)`, which takes
the OpenMM knobs (`nonbonded_method`, `max_iterations`, `platform_name`) that
the pipeline leaves at their defaults.

