# Four backends

| | `backend="gaff"` (default) | `backend="smirnoff"` | `backend="espaloma"` | `backend="charmm"` |
|---|---|---|---|---|
| Parameters | GAFF/GAFF2 atom types, AM1-BCC charges | OpenFF Sage, SMARTS-matched | Espaloma, **predicted** by a graph network | CGenFF, **converted, not derived** |
| Needs | AmberTools on `PATH` | nothing beyond the install | the optional `espaloma` package (PyTorch); the model downloads on first use | a CGenFF stream file for the ligand |
| Ligand source | PDB residue, SDF, MOL2 or SMILES | **SDF, MOL2 or SMILES only** | **SDF, MOL2 or SMILES only** | **CHARMM `.str`/`.rtf`/`.prm` only** |
| Base force field | `amber14` (default) | `amber14` (default) | `amber14` (default) | **`CHARMM_BASE_FORCEFIELD`** |

SMIRNOFF matches SMARTS against the chemical graph and Espaloma reads that graph
directly, and a PDB records no bond orders — hence the `file`-or-`smiles`
requirement, which both state up front.

The three Amber-family backends can be mixed in one call: forcefill writes one
combined XML and OpenMM loads it. That works because SMIRNOFF and Espaloma name
their atom types by a hash of the molecule, so nothing collides, and because the
merge keeps the `<PeriodicTorsionForce>` sections apart — GAFF and SMIRNOFF
impropers use different `ordering` conventions. CHARMM cannot join them; see {doc}`charmm`.

Espaloma takes its model where the others take a force field:

```python
build_ligand_xml("ben.sdf", "ben.xml", backend="espaloma")  # espaloma-0.3.2
build_ligand_xml(
    "ben.sdf", "ben.xml", backend="espaloma", espaloma_forcefield="my_model.pt", espaloma_charge_method="am1-bcc"
)
```

The charge model is always stated explicitly, never left to the library:
openmmforcefields picks a *different* default depending on whether
`template_generator_kwargs` was passed at all, so leaving it unset would make the
charges depend on how the generator happened to be constructed. forcefill asks
for `"nn"` — espaloma's own prediction — unless you say otherwise.

