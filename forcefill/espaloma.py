"""Parameterize a ligand with Espaloma, the graph neural network, instead of GAFF or Sage.

The third of openmmforcefields' template generators, and the third way into the
same output. Where :mod:`forcefill.amber` types atoms and looks parameters up in
``gaff2.dat``, and :mod:`forcefill.smirnoff` matches SMARTS patterns, Espaloma
reads the molecular graph and *predicts* every valence parameter and partial
charge at once. It shares SMIRNOFF's one hard requirement, and for the same
reason: the graph is the input, so **a ligand must arrive as a file with bond
orders (SDF/MOL2) or as a SMILES**, never as a bare PDB residue.

The work is done by
:class:`openmmforcefields.generators.EspalomaTemplateGenerator`, whose
``generate_residue_template`` has the same signature and the same output shape as
the SMIRNOFF one - hashed atom types, a mapped-SMILES residue name this module
rewrites - so everything downstream of it is shared with
:mod:`forcefill.smirnoff` rather than reimplemented.

Three things about that generator shape this module:

    * **espaloma is an optional dependency.** openmmforcefields imports it only
      when the generator is constructed, and it pulls in PyTorch, so it is not
      part of forcefill's install. :func:`require_espaloma` says so before a
      ligand is read rather than after.
    * **the model is downloaded, not installed.** A name like
      ``"espaloma-0.3.2"`` is resolved against ``~/.espaloma`` and fetched from
      GitHub if it is not there, so the first run on a machine needs network
      access. A local ``.pt`` file avoids that.
    * **the charge method depends on how it was called.** openmmforcefields
      defaults ``charge_method`` to ``"from-molecule"`` when no
      ``template_generator_kwargs`` are passed and to ``"nn"`` when an empty dict
      is - two different charge models from the same call. forcefill always
      passes one explicitly, so the charges do not depend on that.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from openmmforcefields.generators import EspalomaTemplateGenerator

from ._spec import PathLike, ResolvedSpec
from .smirnoff import _load_molecule, _rename_residue_template
from .smirnoff import ligand_topology as _smirnoff_ligand_topology

if TYPE_CHECKING:
    from openmm import app, unit

log = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_ESPALOMA_CHARGE_METHOD",
    "ESPALOMA_CHARGE_METHODS",
    "espaloma_residue_ffxml",
    "installed_espaloma_forcefields",
    "ligand_topology",
    "require_espaloma",
]

#: Charge models espaloma's ``create_openmm_system`` accepts.
ESPALOMA_CHARGE_METHODS = ("nn", "am1-bcc", "gasteiger", "from-molecule")

#: Charge model forcefill asks for, always explicitly. ``"nn"`` is espaloma's
#: own prediction, which is the reason to use this backend at all;
#: ``"from-molecule"``, openmmforcefields' default on one of its two code paths,
#: silently falls back to ``"nn"`` when the molecule carries no charges - which
#: is every molecule forcefill builds.
DEFAULT_ESPALOMA_CHARGE_METHOD = "nn"


def installed_espaloma_forcefields() -> list[str]:
    """Names of the Espaloma models openmmforcefields knows, e.g. ``['espaloma-0.3.2']``.

    A *name*, not an installed file: the model itself is fetched on first use.
    ``forcefield`` also accepts a path to a ``.pt`` file or a URL, neither of
    which appears here.
    """
    return list(EspalomaTemplateGenerator.INSTALLED_FORCEFIELDS)


def require_espaloma() -> None:
    """Raise unless the ``espaloma`` package can be imported.

    openmmforcefields raises for this too, but only from inside the generator
    constructor - after the ligand has been read and, for the first ligand on a
    machine, after the model has been downloaded. Checked here, an environment
    without espaloma fails immediately and says how to fix it.
    """
    try:
        import espaloma  # type: ignore[import-not-found]  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "The espaloma backend needs the 'espaloma' package, which is not "
            "installed. It is an optional dependency of openmmforcefields (it "
            "pulls in PyTorch), so forcefill does not install it either:\n"
            "    conda install -c conda-forge espaloma\n"
            "Use backend='smirnoff' or backend='gaff' instead if you do not "
            f"want it.\n  {type(exc).__name__}: {exc}"
        ) from exc


def _model(forcefield: tuple[str, ...]) -> str:
    """The single model name or ``.pt`` path from a resolved spec's selection.

    ``ResolvedSpec.forcefield`` is always a tuple because a SMIRNOFF selection
    can layer several files. An Espaloma model cannot, and openmmforcefields
    rejects anything but a bare string with a bald ``TypeError``, so the arity is
    checked where the spec is built (:func:`forcefill._spec.resolve_specs`) and
    only unwrapped here.
    """
    (model,) = forcefield
    return model


def ligand_topology(spec: ResolvedSpec) -> tuple[app.Topology, unit.Quantity]:
    """Return ``(topology, positions)`` for the spec's molecule, for validating it on its own.

    Identical to :func:`forcefill.smirnoff.ligand_topology`, and deliberately the
    same code: both backends take the molecule as their input, so the topology to
    check the generated template against is the same object built the same way.
    Only the parameters differ.
    """
    return _smirnoff_ligand_topology(spec)


def espaloma_residue_ffxml(
    spec: ResolvedSpec,
    output_xml: PathLike,
    *,
    charge_method: str = DEFAULT_ESPALOMA_CHARGE_METHOD,
) -> str:
    """Write an Espaloma force-field XML for one ligand and return the path.

    Args:
        spec: The ligand. Must carry ``file`` or ``smiles``; ``atom_type`` and
            ``charge_method`` do not apply and are ignored (Espaloma predicts
            both the types and the charges). ``forcefield`` names the model.
        output_xml: Where to write the per-residue XML.
        charge_method: One of :data:`ESPALOMA_CHARGE_METHODS`. Passed explicitly
            on every call - see the module docstring for why.

    Returns:
        The path written, as a string.

    Raises:
        RuntimeError: espaloma is not installed, or openmmforcefields produced
            something unexpected.
        ValueError: The ligand source is unusable, its formal charge contradicts
            an explicit ``net_charge``, or *charge_method* is not one espaloma
            knows.
    """
    if charge_method not in ESPALOMA_CHARGE_METHODS:
        raise ValueError(
            f"charge_method={charge_method!r} is not one espaloma accepts; choose from {list(ESPALOMA_CHARGE_METHODS)}."
        )
    require_espaloma()
    model = _model(spec.forcefield)
    molecule = _load_molecule(spec)
    log.info(
        "espaloma: %s (net charge %+d, %s, charges=%s)",
        spec.name,
        round(molecule.total_charge.m),
        model,
        charge_method,
    )
    generator = EspalomaTemplateGenerator(
        molecules=[molecule],
        forcefield=model,
        template_generator_kwargs={"charge_method": charge_method},
    )
    ffxml = _rename_residue_template(generator.generate_residue_template(molecule), spec.name)

    output_xml = Path(output_xml)
    output_xml.parent.mkdir(parents=True, exist_ok=True)
    output_xml.write_text(ffxml)
    log.info("Wrote per-residue XML: %s", output_xml)
    return str(output_xml)
