# Relation to `openmmforcefields`

[`openmmforcefields`](https://github.com/openmm/openmmforcefields) ships three
template generators, which do the same parameterization on the fly at
`createSystem` time. forcefill uses two of them —
`SMIRNOFFTemplateGenerator` for `backend="smirnoff"` and
`EspalomaTemplateGenerator` for `backend="espaloma"` — through
`generate_residue_template`, which hands back a finished ffxml rather than
registering a callback.

`GAFFTemplateGenerator` is the one forcefill does not use. It shells out to
`antechamber` and `parmchk2` exactly as forcefill's own `gaff` backend does, so
going through it would save no dependency and cost the control over charge
method, net charge, multiplicity and raw antechamber arguments that backend
exposes.

The trade-off forcefill is for is the opposite one: explicit, inspectable,
versionable XML artifacts, produced once, with the skip-classification and
preflight checks telling you which residues need a different treatment entirely.

Two consequences of taking the `generate_residue_template` route are worth
knowing. Its output carries **constraints and virtual sites** (both new in
openmmforcefields 0.16), so forcefill refuses a constrained OFFXML and adds the
extra particles before checking — see {doc}`bespokefit` and {doc}`gotchas`. And the
**TinyDB cache** the generators accept as `cache=` is consulted only on the
callback path, never in `generate_residue_template`, so passing it here would do
nothing; caching the per-residue XMLs is forcefill's own job and is on the
roadmap.

