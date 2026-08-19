"""Tests for the charmm backend: converting CGenFF files into an OpenMM force field.

Hermetic and fast, unlike the gaff tests: the conversion needs no external
executable and no CHARMM toppar download - ``charmm36.xml`` ships with OpenMM and
already carries every CGenFF atom type.

Most of what is checked is *absence*: that the generated XML does not redefine
what the base force field already owns. ParmEd writes those definitions by
default, with zero-valued Lennard-Jones parameters for the types it only knows
the mass of, and loading them silently replaces the real ones - a failure with no
error message, which makes
:func:`test_base_lennard_jones_survives_the_conversion` the test that matters
most here.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

import openmm
from openmm import app, unit

from forcefill import (
    CHARMM_BASE_FORCEFIELD,
    LigandSpec,
    build_forcefield_xml,
    build_ligand_xml,
    charmm,
)
from forcefill._pipeline import check_backends_match_base
from forcefill._spec import ResolvedSpec
from forcefill.preflight import preflight_specs
from tests.helpers import (
    CHLOROETHANOL_ATOMS,
    DATA,
    METHANOL_ATOMS,
    METHANOL_BONDS,
    METHANOL_XYZ,
    write_chloroethanol_pdb,
    write_methanol_pdb,
)

#: The committed CGenFF stream file for the test methanol (residue ``LIG``).
STREAM = DATA / "methanol_cgenff.str"

#: 2-chloroethanol, for the tests that go through a structure. charmm36.xml
#: already has a template matching methanol, so that one never reaches a backend.
CHLOROETHANOL_STREAM = DATA / "chloroethanol_cgenff.str"

#: An ammonium ion: a residue whose charges sum to something other than zero, so
#: "the net charge was read" and "the net charge was assumed" can be told apart.
AMMONIUM = (
    "RESI AMM 1.000\nGROUP\nATOM N1 NG3P3 -0.320\nATOM H1 HGP2 0.330\n"
    "ATOM H2 HGP2 0.330\nATOM H3 HGP2 0.330\nATOM H4 HGP2 0.330\n\n"
    "BOND N1 H1\nBOND N1 H2\nBOND N1 H3\nBOND N1 H4\n"
)

#: What the conversion must produce: the residue template alone, naming CGenFF's
#: own atom types. Hand-written, so the energy comparison below has a reference
#: sharing no code with the thing it checks.
REFERENCE_XML = """<ForceField>
 <Residues>
  <Residue name="LIG">
   <Atom name="C1" type="CG331" charge="-0.04"/>
   <Atom name="O1" type="OG311" charge="-0.65"/>
   <Atom name="H1" type="HGA3" charge="0.09"/>
   <Atom name="H2" type="HGA3" charge="0.09"/>
   <Atom name="H3" type="HGA3" charge="0.09"/>
   <Atom name="H4" type="HGP1" charge="0.42"/>
   <Bond atomName1="C1" atomName2="O1"/>
   <Bond atomName1="C1" atomName2="H1"/>
   <Bond atomName1="C1" atomName2="H2"/>
   <Bond atomName1="C1" atomName2="H3"/>
   <Bond atomName1="O1" atomName2="H4"/>
  </Residue>
 </Residues>
</ForceField>"""

#: What ParmEd writes if the pruning is skipped: the same template, plus atom
#: types and non-bonded entries redefining charmm36's with epsilon=0. OpenMM
#: takes the later definition without complaint - the quiet failure the backend
#: exists to prevent.
UNPRUNED_XML = REFERENCE_XML.replace(
    "</ForceField>",
    """ <AtomTypes>
  <Type element="C" name="CG331" class="CG331" mass="12.011"/>
  <Type element="O" name="OG311" class="OG311" mass="15.9994"/>
  <Type element="H" name="HGA3" class="HGA3" mass="1.008"/>
  <Type element="H" name="HGP1" class="HGP1" mass="1.008"/>
 </AtomTypes>
 <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">
  <UseAttributeFromResidue name="charge"/>
  <Atom class="CG331" sigma="1.0" epsilon="0.0"/>
  <Atom class="OG311" sigma="1.0" epsilon="0.0"/>
  <Atom class="HGA3" sigma="1.0" epsilon="0.0"/>
  <Atom class="HGP1" sigma="1.0" epsilon="0.0"/>
 </NonbondedForce>
 <LennardJonesForce lj14scale="1.0">
  <Atom class="CG331" sigma="1.0" epsilon="0.0"/>
  <Atom class="OG311" sigma="1.0" epsilon="0.0"/>
  <Atom class="HGA3" sigma="1.0" epsilon="0.0"/>
  <Atom class="HGP1" sigma="1.0" epsilon="0.0"/>
 </LennardJonesForce>
</ForceField>""",
)


def convert(tmp_path, files=(STREAM,), name="LIG"):
    """Run one spec through the backend and return the path written."""
    spec = ResolvedSpec(name=name, backend="charmm", charmm_files=tuple(files))
    return charmm.charmm_residue_ffxml(spec, tmp_path / f"{name}.xml", CHARMM_BASE_FORCEFIELD)


def write_stream(path, body):
    """Write a minimal CHARMM stream file with *body* as its rtf section."""
    Path(path).write_text(f"* test\n*\n\nread rtf card append\n* test\n*\n36 1\n\n{body}\n\nEND\nRETURN\n")
    return path


def methanol_pair(name="LIG", separation=8.0):
    """Two copies of the test methanol, far enough apart to have real intermolecular LJ."""
    top = app.Topology()
    chain = top.addChain("A")
    xyz = []
    for copy in range(2):
        residue = top.addResidue(name, chain)
        atoms = {n: top.addAtom(n, element, residue) for n, element in METHANOL_ATOMS}
        for a, b in METHANOL_BONDS:
            top.addBond(atoms[a], atoms[b])
        xyz += [(x + copy * separation, y, z) for x, y, z in METHANOL_XYZ]
    return top, unit.Quantity([openmm.Vec3(*p) for p in xyz], unit.angstrom)


def energy(files, topology, positions, *, uncharged=False):
    """Potential energy of *topology* under *files*; with *uncharged*, what is left is the LJ."""
    system = app.ForceField(*files).createSystem(topology, nonbondedMethod=app.NoCutoff)
    if uncharged:
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    _, sigma, epsilon = force.getParticleParameters(i)
                    force.setParticleParameters(i, 0.0, sigma, epsilon)
                for i in range(force.getNumExceptions()):
                    p, q, _, sigma, epsilon = force.getExceptionParameters(i)
                    force.setExceptionParameters(i, p, q, 0.0, sigma, epsilon)
    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CPU"))
    context.setPositions(positions)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


# --------------------------------------------------------------------------
# What the conversion produces
# --------------------------------------------------------------------------


def test_the_residue_template_is_written(tmp_path):
    root = ET.parse(convert(tmp_path)).getroot()
    (residue,) = root.findall("./Residues/Residue")
    assert residue.get("name") == "LIG"
    assert [(a.get("name"), a.get("type")) for a in residue.findall("Atom")] == [
        ("C1", "CG331"),
        ("O1", "OG311"),
        ("H1", "HGA3"),
        ("H2", "HGA3"),
        ("H3", "HGA3"),
        ("H4", "HGP1"),
    ]
    assert len(residue.findall("Bond")) == len(METHANOL_BONDS)
    assert round(sum(float(a.get("charge")) for a in residue.findall("Atom")), 6) == 0.0


def test_the_base_force_fields_definitions_are_not_repeated(tmp_path):
    root = ET.parse(convert(tmp_path)).getroot()
    # CG331 and friends belong to charmm36.xml. Redefining them here would
    # override it - with epsilon=0, since a stream file carries no LJ terms.
    assert root.find("./AtomTypes") is None
    assert root.findall("./NonbondedForce/Atom") == []
    assert root.find("./LennardJonesForce") is None
    # ...but the charges do have to be declared as living on the residue.
    assert root.findall("./NonbondedForce/UseAttributeFromResidue")


def test_the_stream_files_own_parameters_are_kept(tmp_path):
    root = ET.parse(convert(tmp_path)).getroot()
    # The three terms the stream file's `read param` section defines: these are
    # what ParamChem assigns by analogy, and are the reason it exists.
    assert [b.get("class1") for b in root.findall("./HarmonicBondForce/Bond")] == ["CG331"]
    assert [a.get("class3") for a in root.findall("./HarmonicAngleForce/Angle")] == ["HGP1"]
    assert [t.get("class1") for t in root.findall("./PeriodicTorsionForce/Proper")] == ["HGA3"]


def test_charmm_1_4_scaling_is_declared(tmp_path):
    nonbonded = ET.parse(convert(tmp_path)).getroot().find("./NonbondedForce")
    assert (nonbonded.get("coulomb14scale"), nonbonded.get("lj14scale")) == ("1.0", "1.0")


# --------------------------------------------------------------------------
# What it means once OpenMM reads it back
# --------------------------------------------------------------------------


def test_a_system_builds_on_the_charmm_base(tmp_path):
    topology, positions = methanol_pair()
    result = energy([*CHARMM_BASE_FORCEFIELD, convert(tmp_path)], topology, positions)
    assert result == pytest.approx(144.8, abs=1.0)


def test_base_lennard_jones_survives_the_conversion(tmp_path):
    """The generated XML must leave charmm36's own Lennard-Jones parameters alone.

    ParmEd writes ``sigma=1.0 epsilon=0.0`` for every type whose non-bonded
    parameters it lacks, and OpenMM lets a later file redefine an atom type
    without a word, so leaving those in place quietly changes the physics. Two
    separated copies of the ligand put the intermolecular LJ on the table, where
    the difference shows up.
    """
    topology, positions = methanol_pair()
    reference = tmp_path / "reference.xml"
    reference.write_text(REFERENCE_XML)
    unpruned = tmp_path / "unpruned.xml"
    unpruned.write_text(UNPRUNED_XML)

    def lj(xml_file):
        return energy([*CHARMM_BASE_FORCEFIELD, str(xml_file)], topology, positions, uncharged=True)

    assert lj(convert(tmp_path)) == pytest.approx(lj(reference), rel=1e-9)
    # Teeth: without the pruning the answer is different - not by much, and
    # never loudly, which is the whole problem.
    assert lj(unpruned) != pytest.approx(lj(reference), rel=1e-6)


def charged_spec(tmp_path, **kwargs):
    stream = write_stream(tmp_path / "amm.str", AMMONIUM)
    return ResolvedSpec(name="AMM", backend="charmm", charmm_files=(stream,), **kwargs)


def preflight(spec, tmp_path):
    return preflight_specs({spec.name: spec}, {}, None, tmp_path, base_forcefield=CHARMM_BASE_FORCEFIELD)[spec.name]


def test_the_net_charge_comes_from_the_stream_file(tmp_path):
    # The RESI block's charges sum to the formal charge, the same way an SDF's
    # M CHG record states it - so it never has to be given by hand.
    assert preflight(charged_spec(tmp_path), tmp_path).net_charge == 1


def test_an_agreeing_net_charge_is_left_alone(tmp_path):
    assert preflight(charged_spec(tmp_path, net_charge=1), tmp_path).net_charge == 1


def test_a_contradictory_net_charge_is_refused(tmp_path):
    with pytest.raises(ValueError, match="describes a molecule with a formal charge"):
        preflight(charged_spec(tmp_path, net_charge=0), tmp_path)


def test_ligand_topology_reproduces_the_bond_graph():
    spec = ResolvedSpec(name="LIG", backend="charmm", charmm_files=(STREAM,))
    topology = charmm.ligand_topology(spec)
    assert [r.name for r in topology.residues()] == ["LIG"]
    assert [a.name for a in topology.atoms()] == [name for name, _ in METHANOL_ATOMS]
    assert topology.getNumBonds() == len(METHANOL_BONDS)


# --------------------------------------------------------------------------
# What it refuses
# --------------------------------------------------------------------------


def test_a_file_parmed_cannot_type_is_refused(tmp_path):
    # ParmEd decides what a file is from its name, so this is a real trap: the
    # contents are fine and it still fails.
    renamed = tmp_path / "methanol_cgenff.str.txt"
    renamed.write_text(STREAM.read_text())
    with pytest.raises(ValueError, match="decides the file type from the suffix"):
        convert(tmp_path, files=(renamed,))


def test_no_files_at_all_is_refused():
    with pytest.raises(ValueError, match="at least one CHARMM file"):
        charmm.read_charmm_files([])


def test_an_unknown_atom_type_is_refused(tmp_path):
    stream = write_stream(
        tmp_path / "bad.str",
        "RESI LIG 0.000\nGROUP\nATOM C1 XX999 0.000\nATOM O1 OG311 0.000\n\nBOND C1 O1\n",
    )
    with pytest.raises(ValueError, match=r"atom type\(s\) \['XX999'\]"):
        convert(tmp_path, files=(stream,))


def test_a_file_with_no_residue_is_refused(tmp_path):
    empty = tmp_path / "params_only.str"
    empty.write_text("* test\n*\n\nread param card flex append\n* test\n*\n\nEND\nRETURN\n")
    with pytest.raises(ValueError, match="defines no residue template"):
        convert(tmp_path, files=(empty,))


def test_an_empty_document_is_refused(tmp_path, monkeypatch):
    """The backstop for ParmEd dropping a residue template with only a warning.

    Nothing forcefill can pass it triggers that today - unresolvable atom types
    are caught first - so the writer is stubbed. The guard exists because the
    failure is silent: an ffxml that loads fine and parameterizes nothing.
    """
    from parmed.openmm import OpenMMParameterSet

    def write_nothing(self, dest, **kwargs):
        Path(dest).write_text("<ForceField/>")

    monkeypatch.setattr(OpenMMParameterSet, "write", write_nothing)
    with pytest.raises(RuntimeError, match="contains no residue template"):
        convert(tmp_path)


def test_a_lone_residue_is_renamed_to_the_one_asked_for(tmp_path):
    # The RESI name is whatever was typed into ParamChem; the structure decides
    # what the residue is actually called.
    root = ET.parse(convert(tmp_path, name="ABC")).getroot()
    assert [r.get("name") for r in root.findall("./Residues/Residue")] == ["ABC"]


def test_several_residues_with_no_match_is_refused(tmp_path):
    stream = write_stream(
        tmp_path / "two.str",
        "RESI AAA 0.000\nGROUP\nATOM C1 CG331 0.000\nATOM O1 OG311 0.000\n\nBOND C1 O1\n\n"
        "RESI BBB 0.000\nGROUP\nATOM C1 CG331 0.000\nATOM O1 OG311 0.000\n\nBOND C1 O1\n",
    )
    with pytest.raises(ValueError, match=r"defines 2 residue templates .*'AAA', 'BBB'"):
        convert(tmp_path, files=(stream,), name="LIG")


# --------------------------------------------------------------------------
# The combinations OpenMM could never load
# --------------------------------------------------------------------------


def charmm_spec(name="LIG"):
    return ResolvedSpec(name=name, backend="charmm", charmm_files=(STREAM,))


def test_charmm_with_the_amber_base_is_refused():
    with pytest.raises(ValueError, match="1-4 scaling 1/1"):
        check_backends_match_base({"LIG": charmm_spec()}, ("amber14-all.xml", "amber14/tip3p.xml"))


def test_gaff_with_the_charmm_base_is_refused():
    with pytest.raises(ValueError, match=r"gaff backend produces parameters with 1-4 scaling 0\.8333"):
        check_backends_match_base({"LIG": ResolvedSpec(name="LIG")}, CHARMM_BASE_FORCEFIELD)


def test_charmm_mixed_with_gaff_is_refused():
    specs = {"LIG": charmm_spec(), "BEN": ResolvedSpec(name="BEN", backend="gaff")}
    with pytest.raises(ValueError, match="both CHARMM and Amber-family"):
        check_backends_match_base(specs, CHARMM_BASE_FORCEFIELD)


def test_matching_combinations_pass():
    check_backends_match_base({"LIG": charmm_spec()}, CHARMM_BASE_FORCEFIELD)
    check_backends_match_base({"LIG": ResolvedSpec(name="LIG")}, ("amber14-all.xml", "amber14/tip3p.xml"))
    # An empty base force field declares nothing, so nothing can contradict it.
    check_backends_match_base({"LIG": charmm_spec()}, ())
    check_backends_match_base({}, CHARMM_BASE_FORCEFIELD)


# --------------------------------------------------------------------------
# Through the two entry points
# --------------------------------------------------------------------------


def test_build_ligand_xml_takes_a_stream_file_directly(tmp_path):
    # A CHARMM suffix is unambiguous, so the backend and the residue name can
    # both be read off the path - as they are for an SDF.
    stream = tmp_path / "lig.str"
    stream.write_text(STREAM.read_text())
    result = build_ligand_xml(
        stream,
        tmp_path / "out.xml",
        base_forcefield=CHARMM_BASE_FORCEFIELD,
        workdir=tmp_path / "wd",
    )
    assert result.parameterized == ["LIG"]
    assert result.minimizations == {}
    assert '<Residue name="LIG">' in Path(result.forcefield_xml).read_text()


def test_build_ligand_xml_names_the_ligand_explicitly(tmp_path):
    result = build_ligand_xml(
        {"ABC": LigandSpec(backend="charmm", charmm_files=(STREAM,))},
        tmp_path / "out.xml",
        base_forcefield=CHARMM_BASE_FORCEFIELD,
        workdir=tmp_path / "wd",
    )
    assert result.parameterized == ["ABC"]


def test_a_name_the_base_force_field_already_uses_is_refused(tmp_path):
    # charmm36.xml has 814 templates, so this is easy to hit: MET is methionine.
    with pytest.raises(ValueError, match="already defines a residue template named MET"):
        build_ligand_xml(
            {"MET": LigandSpec(backend="charmm", charmm_files=(STREAM,))},
            tmp_path / "out.xml",
            base_forcefield=CHARMM_BASE_FORCEFIELD,
            workdir=tmp_path / "wd",
        )


def test_two_charmm_ligands_merge_into_one_file(tmp_path):
    result = build_ligand_xml(
        {
            "LIG": LigandSpec(backend="charmm", charmm_files=(STREAM,)),
            "CET": LigandSpec(backend="charmm", charmm_files=(CHLOROETHANOL_STREAM,)),
        },
        tmp_path / "out.xml",
        base_forcefield=CHARMM_BASE_FORCEFIELD,
        workdir=tmp_path / "wd",
    )
    assert result.parameterized == ["CET", "LIG"]
    root = ET.parse(result.forcefield_xml).getroot()
    assert sorted(r.get("name") for r in root.findall("./Residues/Residue")) == ["CET", "LIG"]
    # Still nothing the base force field already owns, after merging.
    assert root.find("./AtomTypes") is None
    # One NonbondedForce, not two: merge_ffxml folds the compatible sections and
    # drops the second <UseAttributeFromResidue>, which OpenMM would reject.
    (nonbonded,) = root.findall("./NonbondedForce")
    assert len(nonbonded.findall("UseAttributeFromResidue")) == 1


def test_build_ligand_xml_refuses_to_minimize_a_charmm_ligand(tmp_path):
    # A stream file records internal coordinates, so every atom would start at
    # the origin - an infinite Coulomb energy dressed up as a parameter problem.
    with pytest.raises(ValueError, match="no geometry to minimize"):
        build_ligand_xml(
            {"LIG": LigandSpec(backend="charmm", charmm_files=(STREAM,))},
            tmp_path / "out.xml",
            base_forcefield=CHARMM_BASE_FORCEFIELD,
            workdir=tmp_path / "wd",
            minimize=True,
        )


def test_build_forcefield_xml_parameterizes_a_pdb_residue(tmp_path):
    pdb = write_chloroethanol_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml(
        pdb,
        tmp_path / "extras.xml",
        base_forcefield=CHARMM_BASE_FORCEFIELD,
        backend="charmm",
        ligands={"CET": LigandSpec(charmm_files=(CHLOROETHANOL_STREAM,))},
        workdir=tmp_path / "wd",
        minimize=True,
    )
    assert result.parameterized == ["CET"]
    assert result.skipped == {}
    # The structure supplies the coordinates the stream file does not, so the
    # numbers can be checked here even though build_ligand_xml cannot.
    assert result.minimizations["CET"].energy_change < 0
    assert result.full_minimization.n_atoms == len(CHLOROETHANOL_ATOMS)


def test_charmm36_matches_the_small_ligands_it_already_ships(tmp_path):
    # charmm36.xml carries every CGenFF model compound, so a fragment-sized
    # ligand may need nothing at all - and then no backend runs.
    pdb = write_methanol_pdb(tmp_path / "in.pdb")
    result = build_forcefield_xml(
        pdb, tmp_path / "extras.xml", base_forcefield=CHARMM_BASE_FORCEFIELD, backend="charmm"
    )
    assert result.forcefield_xml is None
    assert result.parameterized == []


def test_a_stream_file_for_a_different_molecule_is_caught_first(tmp_path):
    pdb = write_chloroethanol_pdb(tmp_path / "in.pdb")
    stream = write_stream(
        tmp_path / "wrong.str",
        "RESI CET 0.000\nGROUP\nATOM C1 CG331 0.000\nATOM O1 OG311 0.000\n\nBOND C1 O1\n",
    )
    with pytest.raises(ValueError, match="not the same molecule as the residue"):
        build_forcefield_xml(
            pdb,
            tmp_path / "extras.xml",
            base_forcefield=CHARMM_BASE_FORCEFIELD,
            backend="charmm",
            ligands={"CET": LigandSpec(charmm_files=(stream,))},
            workdir=tmp_path / "wd",
        )
