# Things to get right

- **Explicit hydrogens.** Ligands must carry all of them, with reasonable
  geometry; AM1-BCC charges are meaningless otherwise. forcefill warns when a
  ligand has none.
- **Net charge for a PDB-extracted ligand.** There is nothing to read it from,
  so pass `net_charges={"RES": q}` yourself. forcefill warns about keys that
  match no residue (typos, case mismatches).
- **Connectivity.** Element columns and (for hetero groups) CONECT records
  should be present in the PDB.
- **One XML at a time.** Load either the combined XML *or* the per-residue
  XMLs, never both — the duplicated GAFF atom-type definitions would collide.
- **Cleaning changes what the checks describe.** With `clean_structure=True` the
  full-structure `validate`/`minimize` results, and so
  `full_minimization.n_atoms`, refer to the *cleaned* topology, not the file on
  disk. Reconcile against `cleaning.n_atoms_after`.
- **The periodic box survives the strip.** A de-solvated structure keeps the box
  vectors of the solvated one, so a later PME run would use a mostly-empty box.
  Reset them yourself before simulating it directly.
- **Virtual sites are extra particles, and they are yours to add.** A SMIRNOFF
  force field with a `VirtualSites` handler — a sigma-hole model on a halogen,
  say — makes openmmforcefields write `<VirtualSite>` into the residue template,
  so the template describes more particles than your topology has atoms and
  OpenMM matches *nothing*: "the residue is missing 1 extra site". forcefill's
  own checks handle it, and `MinimizationResult.n_atoms` then counts particles
  rather than atoms. When you load the XML yourself, add them first:

  ```python
  modeller = Modeller(pdb.topology, pdb.positions)
  modeller.addExtraParticles(ff)  # or forcefill.add_extra_particles(...)
  system = ff.createSystem(modeller.topology)
  ```

  and let OpenMM place them — `context.computeVirtualSites()` after
  `setPositions`. Modeller cannot reconstruct a two-particle bond-charge frame
  and leaves the site somewhere else entirely, or at NaN.

