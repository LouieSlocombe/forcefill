# What gets skipped, and why

Most of the value is in what forcefill *refuses* to parameterize. Running
antechamber on the wrong kind of residue produces plausible-looking but
physically wrong parameters, so these are reported and skipped:

| Unmatched residue | Action | Do this instead |
|---|---|---|
| Standard residue (e.g. `ALA`, `HOH`) that failed to match | skip | It is missing atoms or has non-standard atom names — repair the structure with [PDBFixer](https://github.com/openmm/pdbfixer) or `Modeller.addHydrogens` |
| Monatomic species (ions such as `ZN`, `NA`) | skip | GAFF/antechamber cannot treat bare ions. The standard base force fields already define the common ones, so one that *still* failed to match is usually a residue- or atom-name mismatch, or a charge state no template carries — rename it, or load a set that covers it (note openmmforcefields' `amber/ions/*.xml` **replace** the ions in the bundled water file rather than adding to them, so they collide with `amber14-all.xml`) |
| Residue covalently bonded to its neighbours (modified amino acids, glycans) | skip | Stand-alone GAFF is not valid for polymer-linked residues; cap the fragment and derive charges consistently with the backbone force field (pyRED- or ffparam-style workflows) |
| Free-standing hetero molecule (ligand, cofactor) | **parameterize** | — |

The last row cuts both ways: a glycerol or a sulfate left over from the
crystallization drop *is* a free-standing hetero molecule, so it gets
parameterized too. Strip those first — see {doc}`cleaning`.

