# Defaults and residue tables

Every name on this page is importable from the top level — `from forcefill import
DEFAULT_BASE_FORCEFIELD` — even though the tables that define them live in
internal modules.

The two remaining public constants are documented alongside the code that uses
them: `DEFAULT_AMBERTOOLS_TIMEOUT` on {doc}`amber`, and
`DEFAULT_MINIMIZATION_PLATFORM` / `DEFAULT_MINIMIZATION_TOLERANCE` on
{doc}`checks`.

## Backends and force fields

```{eval-rst}
.. currentmodule:: forcefill._spec

.. autodata:: BACKENDS
.. autodata:: DEFAULT_BASE_FORCEFIELD
.. autodata:: CHARMM_BASE_FORCEFIELD
.. autodata:: DEFAULT_SMIRNOFF_FORCEFIELD
.. autodata:: DEFAULT_ESPALOMA_FORCEFIELD
```

## Residue-name tables

How {func}`~forcefill.clean_topology` classifies each residue it meets. They are
plain `frozenset`s of PDB chemical component IDs, so you can test membership,
extend one for a `keep=` or `extra_remove=` argument, or read them to find out
what a default run would have deleted.

```{eval-rst}
.. currentmodule:: forcefill._residue_names

.. autodata:: WATER_RESIDUES
.. autodata:: BULK_ION_RESIDUES
.. autodata:: STRUCTURAL_METAL_RESIDUES
.. autodata:: ADDITIVE_RESIDUES
```
