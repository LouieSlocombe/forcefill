"""Prove the generated force field works, by building an ``openmm.System`` from it.

The checks that run *after* parameterization; the input-side ones are
:mod:`forcefill.preflight`. There are two:

    * :func:`validate_forcefield_xml` builds a System, which proves the generated
      template matches the topology's bond graph and that no parameter is
      missing - and nothing else.
    * :func:`minimize_with_forcefield_xml` also computes an energy and takes a
      few minimizer steps. A NaN charge, a zero force constant or a broken angle
      term all survive a System build and surface later as an exploding
      simulation; this is what catches them.

Both run per residue as well as on the whole structure, so a template that does
not match reports itself as that rather than as an incomplete structure.

Both also have to account for **virtual sites**. openmmforcefields >= 0.16 writes
a ``<VirtualSite>`` into the residue template for every extra site a SMIRNOFF
force field assigns, so a force field carrying a ``VirtualSites`` handler - a
sigma-hole model on a halogen, say - describes more particles than the molecule
has atoms. Neither an input PDB nor a topology built from the molecule carries
those particles, and OpenMM then matches nothing at all: *"the residue is missing
1 extra site"*. :func:`add_extra_particles` puts them in, and
:func:`residue_templates_with_virtual_sites` says when that is needed.
"""

from __future__ import annotations

import logging
import math
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from openmm import Context, LocalEnergyMinimizer, Platform, VerletIntegrator, app, unit

from ._spec import DEFAULT_BASE_FORCEFIELD, PathLike
from .topology import _describe_topology, _residue_positions, _residue_subtopology

log = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_MINIMIZATION_PLATFORM",
    "DEFAULT_MINIMIZATION_TOLERANCE",
    "MinimizationResult",
    "add_extra_particles",
    "minimize_with_forcefield_xml",
    "residue_templates_with_virtual_sites",
    "validate_forcefield_xml",
]

#: Platform for the minimization checks. Pinned rather than left to OpenMM's
#: fastest-available pick, which could take a GPU the caller wanted elsewhere.
DEFAULT_MINIMIZATION_PLATFORM = "CPU"

#: Minimizer convergence target on the RMS force, in kJ/mol/nm (OpenMM's own default).
DEFAULT_MINIMIZATION_TOLERANCE: float = 10.0

#: Unit the reported forces are converted to.
_FORCE_UNIT = unit.kilojoule_per_mole / unit.nanometer


@dataclass
class MinimizationResult:
    """What one :func:`minimize_with_forcefield_xml` run measured.

    Plain floats rather than ``unit.Quantity``: energies in kJ/mol, forces in
    kJ/mol/nm. A Quantity cannot be ``math.isfinite``-checked or ``%``-formatted,
    which is all this object is ever used for.
    """

    #: Number of atoms in the minimized topology.
    n_atoms: int
    #: Potential energy of the input coordinates.
    initial_energy: float
    #: Potential energy after minimizing.
    final_energy: float
    #: Largest per-atom force magnitude after minimizing.
    max_force: float

    @property
    def energy_change(self) -> float:
        """Final minus initial energy; negative when the minimizer did its job."""
        return self.final_energy - self.initial_energy


# --------------------------------------------------------------------------
# Virtual sites: extra particles the templates declare but no topology carries
# --------------------------------------------------------------------------


def residue_templates_with_virtual_sites(xml_file: PathLike) -> list[str]:
    """Names of the residue templates in *xml_file* that declare virtual sites.

    Read from the XML rather than from a loaded ``ForceField``, whose templates
    are private: this is the same document OpenMM will read, and the question is
    only whether any ``<Residue>`` in it has a ``<VirtualSite>`` child.

    Returns:
        Residue names, in document order; empty when the file declares none,
        which is the case for every gaff build and for SMIRNOFF force fields
        without a ``VirtualSites`` handler.
    """
    root = ET.parse(str(xml_file)).getroot()
    return [
        residue.get("name", "")
        for residue in root.findall("./Residues/Residue")
        if residue.find("VirtualSite") is not None
    ]


def add_extra_particles(
    topology: app.Topology,
    positions: unit.Quantity,
    forcefield: app.ForceField,
) -> tuple[app.Topology, unit.Quantity]:
    """Return *topology* and *positions* with *forcefield*'s virtual sites added.

    A thin wrapper over ``Modeller.addExtraParticles``, which is what OpenMM
    itself points at when a template describes more particles than the residue
    has atoms.

    The returned coordinates for the *extra sites* are placeholders: Modeller
    reconstructs a site from the atoms defining its frame, and for the
    two-particle ``localCoords`` frame openmmforcefields writes for a SMIRNOFF
    bond-charge site that frame is degenerate - the result comes back somewhere
    else entirely, or as NaN. The ``System`` is the authority on where an extra
    site belongs, so a caller must let OpenMM place them:
    ``context.setPositions(positions)`` then ``context.computeVirtualSites()``,
    which is required for any system with virtual sites and is what
    :func:`minimize_with_forcefield_xml` does. Positions of the *real* atoms come
    back untouched.

    Args:
        topology: Structure the generated templates will be matched against.
        positions: Coordinates for *topology*.
        forcefield: The force field whose templates declare the extra sites -
            base force field plus the generated XML.

    Returns:
        ``(topology, positions)``, both extended. The atom count grows, so a
        :class:`MinimizationResult` measured afterwards reports the particle
        count rather than the number of real atoms.
    """
    modeller = app.Modeller(topology, positions)
    with warnings.catch_warnings():
        # A degenerate frame surfaces as a numpy RuntimeWarning from openmm.vec3.
        # It is expected here, and whatever it produces is overwritten by
        # computeVirtualSites, so it would only read as a fault.
        warnings.simplefilter("ignore", RuntimeWarning)
        modeller.addExtraParticles(forcefield)
    added = modeller.topology.getNumAtoms() - topology.getNumAtoms()
    if added:
        log.info("Added %d virtual-site particle(s) before checking %s", added, _describe_topology(topology))
    return modeller.topology, modeller.positions


def _expand_virtual_sites(
    topology: app.Topology,
    positions: unit.Quantity | None,
    forcefield: app.ForceField,
    *,
    needed: bool,
) -> tuple[app.Topology, unit.Quantity | None]:
    """Add extra particles when the generated templates declare them, else pass through.

    Kept separate from :func:`add_extra_particles` so the callers stay readable:
    the common case is that nothing declares a virtual site and the topology must
    reach OpenMM untouched, and the charmm backend has no Cartesian coordinates
    to place a site from.
    """
    if not needed or positions is None:
        return topology, positions
    return add_extra_particles(topology, positions, forcefield)


# --------------------------------------------------------------------------
# Validation: does a System build at all
# --------------------------------------------------------------------------


def _validate_parameterized_residues(
    residues: Mapping[str, app.topology.Residue],
    forcefield: app.ForceField,
    files: Sequence[str],
    positions: unit.Quantity | None = None,
    *,
    virtual_sites: bool = False,
) -> None:
    """Check that *forcefield* (built from *files*) makes a System for each residue on its own.

    Proves each generated template matches its residue's bond graph and that no
    parameter is missing, independently of whether the rest of the structure is
    complete. Raises on the first residue that fails; *files* only labels the
    error message.

    *positions* and *virtual_sites* are needed together, and only when a
    generated template declares extra sites: placing one takes the coordinates
    of the atoms defining its frame.
    """
    files = list(files)
    name = None
    try:
        for name in sorted(residues):
            topology, _ = _expand_virtual_sites(
                _residue_subtopology(residues[name]),
                _residue_positions(positions, residues[name]) if positions is not None else None,
                forcefield,
                needed=virtual_sites,
            )
            forcefield.createSystem(topology)
    except Exception as exc:
        raise RuntimeError(
            f"Validation failed: could not build an openmm.System for "
            f"residue {name} on its own from {files}.\n"
            f"OpenMM said: {exc}\n"
            f"{_diagnose(exc)}"
        ) from exc
    log.info("Validation OK: per-residue Systems built for %s from %s", sorted(residues), files)


#: What OpenMM says when a template describes particles the residue does not have.
_EXTRA_SITE_MARKERS = ("extra site", "extra particle")


def _diagnose(exc: Exception) -> str:
    """The advice to append to a failed System build, chosen from what OpenMM said.

    A missing extra site reads like a missing atom but is neither the structure's
    fault nor repairable with PDBFixer, so it gets its own answer: the particle
    belongs to the force field and :func:`add_extra_particles` is what adds it.
    """
    if any(marker in str(exc).lower() for marker in _EXTRA_SITE_MARKERS):
        return (
            "The generated template declares virtual sites that the topology "
            "does not carry - a SMIRNOFF force field with a VirtualSites "
            "handler, whose extra particles openmmforcefields writes into the "
            "residue template. This is not a fault in the structure: add the "
            "particles with forcefill.checks.add_extra_particles(topology, "
            "positions, forcefield), or OpenMM's Modeller.addExtraParticles(), "
            "before building the System."
        )
    return (
        "The generated template does not match the residue's bond graph, or "
        "parameters are missing. If you supplied this residue via "
        "residue_files, its atoms/bonds (including hydrogens) must match the "
        "residue in the PDB exactly."
    )


def validate_forcefield_xml(
    topology: app.Topology,
    xml_file: PathLike,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    *,
    forcefield: app.ForceField | None = None,
) -> None:
    """Raise RuntimeError unless ``base_forcefield + xml_file`` can build a System for *topology*.

    A pre-built *forcefield* - which must have been constructed from
    ``base_forcefield + xml_file`` - skips re-parsing the XML files.
    """
    files = [*base_forcefield, str(xml_file)]
    try:
        if forcefield is None:
            forcefield = app.ForceField(*files)
        forcefield.createSystem(topology)
    except Exception as exc:
        hint = (
            _diagnose(exc)
            if any(marker in str(exc).lower() for marker in _EXTRA_SITE_MARKERS)
            else (
                "Common causes: other parts of the structure are missing atoms "
                "or hydrogens (repair with PDBFixer first), or the system "
                "contains ions/waters that need an additional parameter file."
            )
        )
        raise RuntimeError(
            f"Validation failed: could not build an openmm.System from "
            f"{files}.\nOpenMM said: {exc}\n"
            f"{hint} Pass validate=False to skip this check."
        ) from exc
    log.info("Validation OK: System built from %s", files)


# --------------------------------------------------------------------------
# Minimization: are the numbers in it physical
# --------------------------------------------------------------------------


class _NonFiniteEnergyError(RuntimeError):
    """Raised by the finite-energy check so it passes through the OpenMM error wrapper unchanged."""


def _require_finite_energy(energy: float, when: str, subject: str, files: Sequence[str]) -> None:
    """Raise unless *energy* is finite.

    A non-finite energy is the failure this whole check exists to catch: it
    means the parameters are unphysical, which building a System cannot detect.
    """
    if math.isfinite(energy):
        return
    raise _NonFiniteEnergyError(
        f"Minimization failed: the potential energy of {subject} is {energy} "
        f"{when} with {list(files)}.\n"
        "The parameters are not physical. Look for NaN charges or zero-valued "
        "equilibrium bond lengths / angles in the generated XML, and for atoms "
        "sharing coordinates in the input structure."
    )


def minimize_with_forcefield_xml(
    topology: app.Topology,
    positions: unit.Quantity,
    xml_file: PathLike,
    base_forcefield: Sequence[str] = DEFAULT_BASE_FORCEFIELD,
    *,
    forcefield: app.ForceField | None = None,
    max_iterations: int = 100,
    tolerance: float = DEFAULT_MINIMIZATION_TOLERANCE,
    nonbonded_method: Any = app.NoCutoff,
    constraints: Any = None,
    rigid_water: bool = False,
    platform_name: str = DEFAULT_MINIMIZATION_PLATFORM,
) -> MinimizationResult:
    """Energy-minimize *topology* with ``base_forcefield + xml_file`` and report what happened.

    The check :func:`validate_forcefield_xml` cannot make: a NaN charge, a zero
    force constant or a broken angle term all survive a System build and surface
    later as an exploding simulation.

    Raises RuntimeError if the XML will not load, the System will not build, or
    the potential energy is not finite before or after minimizing. The reported
    ``max_force`` and ``energy_change`` are diagnostics, not pass/fail criteria -
    what counts as converged depends on the system.

    Args:
        topology: Structure to minimize.
        positions: Coordinates for *topology*, indexed by global atom index and
            in the same order as ``topology.atoms()``.
        xml_file: The generated force-field XML to load on top of
            *base_forcefield*.
        base_forcefield: ffxml files loaded first, e.g. the standard Amber set.
        forcefield: A pre-built ForceField, which must have been constructed
            from ``base_forcefield + xml_file``; skips re-parsing the XML files.
        max_iterations: Minimizer iteration ceiling - a sanity check, not
            convergence. Pass 0 to run until *tolerance* is met.
        tolerance: Convergence target in kJ/mol/nm, applied to the RMS over all
            force *components*, so not comparable to the per-atom ``max_force``
            reported back. Must be positive: OpenMM accepts a negative tolerance
            and then silently minimizes nothing.
        nonbonded_method: ``createSystem`` nonbonded method. The default,
            ``app.NoCutoff``, is exact but O(N^2); pass ``app.PME`` for a
            solvated periodic box.
        constraints: ``createSystem`` constraints. Deliberately ``None`` rather
            than ``app.HBonds``: constraining bonds would hide a bad bond
            parameter.
        rigid_water: Deliberately False, unlike OpenMM's default. With
            constraints present ``getForces`` returns the unconstrained forces,
            making the reported ``max_force`` a meaningless residual.
        platform_name: OpenMM platform, pinned to CPU so this check cannot take
            a GPU the caller wanted for something else.

    Returns:
        MinimizationResult
    """
    if tolerance <= 0:
        raise ValueError(f"tolerance={tolerance!r} must be positive; OpenMM silently skips minimization otherwise.")
    if max_iterations < 0:
        raise ValueError(f"max_iterations={max_iterations!r} must be >= 0 (0 means 'run until converged').")

    files = [*base_forcefield, str(xml_file)]
    subject = _describe_topology(topology)
    try:
        if forcefield is None:
            forcefield = app.ForceField(*files)
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=nonbonded_method,
            constraints=constraints,
            rigidWater=rigid_water,
        )
        # The integrator is never stepped; a Context just requires one.
        context = Context(system, VerletIntegrator(0.001 * unit.picoseconds), Platform.getPlatformByName(platform_name))
        context.setPositions(positions)
        # Required for any system with virtual sites, and a no-op without them.
        # Extra sites arrive from add_extra_particles as placeholders - misplaced
        # or NaN, depending on the frame - and this is what puts them where the
        # System says they belong. Without it the energy measured below is of a
        # geometry the force field never describes.
        context.computeVirtualSites()

        initial = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        _require_finite_energy(initial, "before minimizing", subject, files)

        LocalEnergyMinimizer.minimize(context, tolerance, max_iterations)

        state = context.getState(getEnergy=True, getForces=True)
        final = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        _require_finite_energy(final, "after minimizing", subject, files)
        forces = state.getForces().value_in_unit(_FORCE_UNIT)
    except _NonFiniteEnergyError:
        raise
    except Exception as exc:
        raise RuntimeError(
            f"Minimization failed: could not minimize {subject} with {files}.\n"
            f"OpenMM said: {exc}\n"
            "Either the generated template does not match the topology, or a "
            "parameter file is missing for some other part of it."
        ) from exc

    result = MinimizationResult(
        n_atoms=topology.getNumAtoms(),
        initial_energy=initial,
        final_energy=final,
        max_force=max((math.sqrt(f[0] ** 2 + f[1] ** 2 + f[2] ** 2) for f in forces), default=0.0),
    )
    log.info(
        "Minimization OK: %s went %.1f -> %.1f kJ/mol (max force %.1f kJ/mol/nm) with %s",
        subject,
        result.initial_energy,
        result.final_energy,
        result.max_force,
        files,
    )
    return result


def _minimize_parameterized_residues(
    residues: Mapping[str, app.topology.Residue],
    positions: unit.Quantity,
    xml_file: PathLike,
    base_forcefield: Sequence[str],
    forcefield: app.ForceField,
    *,
    virtual_sites: bool = False,
    **kwargs: Any,
) -> dict[str, MinimizationResult]:
    """Minimize each residue on its own in vacuum; returns the reports keyed by residue name.

    The counterpart of :func:`_validate_parameterized_residues`: it checks the
    numbers rather than the graph, and for the same reason - each residue is
    tested independently of whether the rest of the input is complete.
    """
    out: dict[str, MinimizationResult] = {}
    for name in sorted(residues):
        topology, residue_positions = _expand_virtual_sites(
            _residue_subtopology(residues[name]),
            _residue_positions(positions, residues[name]),
            forcefield,
            needed=virtual_sites,
        )
        out[name] = minimize_with_forcefield_xml(
            topology,
            residue_positions,
            xml_file,
            base_forcefield,
            forcefield=forcefield,
            **kwargs,
        )
    return out
