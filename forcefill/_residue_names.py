"""PDB residue-name tables shared by the parameterization and cleaning steps.

Pure data: this module imports nothing from the rest of the package, so both
:mod:`forcefill.nonstandard_ffxml` and :mod:`forcefill.clean_structure` can
depend on it without a cycle.

The names are PDB chemical component IDs (the residue-name column, 18-20), not
element symbols and not atom names. That distinction matters most for ``CA``:
as a *residue* name it is a calcium ion, as an *atom* name it is the alpha
carbon of every amino acid. Everything here is matched against residue names
only.

Deliberately absent from :data:`ADDITIVE_RESIDUES`, because removing them by
default would be wrong more often than right:

* ``BEN`` - benzamidine, a genuine inhibitor (and this package's own example
  ligand).
* ``NAG``, ``BMA``, ``MAN``, ``FUC``, ``GLC`` and the other glycans - usually
  covalently linked to the protein and part of the structure.
* ``HEM``, ``NAD``, ``FAD``, ``ATP``, ``ADP``, ``SAM`` and other cofactors.
* ``UNL``, ``UNX``, ``UNK`` - unidentified density; deleting it silently is too
  aggressive a default.
* ``EDT`` (EDTA) - a chelator, which may be a real feature of the structure.

Each is one ``extra_remove=(...)`` away. In the other direction, ``IMD``
(imidazole), ``AZI`` (azide) and ``SO4``/``PO4`` *are* removed by default but
can occupy genuine binding sites; ``keep=(...)`` covers that case.
"""

from __future__ import annotations

#: Water, including the model-specific aliases the base force fields use.
#: ``DOD`` is deuterated water, which appears in neutron structures.
WATER_RESIDUES = frozenset(
    [
        "HOH",
        "WAT",
        "H2O",
        "TIP",
        "TIP3",
        "TP3",
        "SPC",
        "SOL",
        "DOD",
    ]
)

#: The 20 amino acids plus the protonation variants and caps the Amber force
#: fields name separately.
_AMINO_ACIDS = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ASH",
        "GLH",
        "LYN",
        "CYX",
        "CYM",
        "HID",
        "HIE",
        "HIP",
        "ACE",
        "NME",
        "NMA",
    ]
)

#: DNA and RNA residues, including the 3'/5' terminal variants.
_NUCLEIC_ACIDS = frozenset(
    [
        "DA",
        "DC",
        "DG",
        "DT",
        "DA3",
        "DA5",
        "DC3",
        "DC5",
        "DG3",
        "DG5",
        "DT3",
        "DT5",
        "A",
        "C",
        "G",
        "U",
        "A3",
        "A5",
        "C3",
        "C5",
        "G3",
        "G5",
        "U3",
        "U5",
    ]
)

#: Residue names the base force fields already know. Unmatched residues with
#: these names are almost always incomplete structures, not new chemistry.
STANDARD_RESIDUES = _AMINO_ACIDS | _NUCLEIC_ACIDS | WATER_RESIDUES

#: Monatomic salt ions: group-1 cations and group-17 anions, plus the CHARMM
#: and Amber aliases. These come from the buffer or from neutralizing the
#: simulation box - they occupy no defined site, coordinate nothing
#: directionally, and you will re-add them with ``Modeller.addSolvent``
#: anyway - so the cleaner removes them by default.
BULK_ION_RESIDUES = frozenset(
    [
        "NA",
        "NA+",
        "SOD",
        "K",
        "K+",
        "POT",
        "LI",
        "LI+",
        "RB",
        "RB+",
        "CS",
        "CS+",
        "CES",
        "CL",
        "CL-",
        "CLA",
        "BR",
        "BR-",
        "IOD",
        "I",
        "F",
    ]
)

#: Metals that are frequently structural or catalytic: buried, directionally
#: coordinated by specific side chains, and often required for the fold or the
#: chemistry. Trypsin's Ca2+ (PDB 3PTB, residue CA 480) is the worked example.
#: Kept by default and reported in ``CleaningResult.retained``, because
#: deleting a needed metal is silent and wrong while keeping an unwanted one is
#: visible and reversible. The phasing heavy atoms (Hg, Pt, Au, Pb and the
#: lanthanides) sit here for the same reason, even though stripping them is
#: usually right.
STRUCTURAL_METAL_RESIDUES = frozenset(
    [
        "CA",
        "CAL",
        "MG",
        "ZN",
        "ZN2",
        "MN",
        "MN3",
        "FE",
        "FE2",
        "CU",
        "CU1",
        "NI",
        "NI2",
        "CO",
        "3CO",
        "CD",
        "HG",
        "PB",
        "PT",
        "AU",
        "AG",
        "SR",
        "BA",
        "MO",
        "W",
        "V",
        "CR",
        "SM",
        "EU",
        "GD",
        "YB",
        "TB",
        "LA",
        "CE",
    ]
)

#: Ion residue name -> element symbol. Every name in
#: :data:`BULK_ION_RESIDUES` and :data:`STRUCTURAL_METAL_RESIDUES` appears here,
#: so an ion is never deleted on the strength of its name alone: the residue's
#: single atom must carry the matching element too.
ION_ELEMENTS = {
    "NA": "Na",
    "NA+": "Na",
    "SOD": "Na",
    "K": "K",
    "K+": "K",
    "POT": "K",
    "LI": "Li",
    "LI+": "Li",
    "RB": "Rb",
    "RB+": "Rb",
    "CS": "Cs",
    "CS+": "Cs",
    "CES": "Cs",
    "CL": "Cl",
    "CL-": "Cl",
    "CLA": "Cl",
    "BR": "Br",
    "BR-": "Br",
    "IOD": "I",
    "I": "I",
    "F": "F",
    "CA": "Ca",
    "CAL": "Ca",
    "MG": "Mg",
    "ZN": "Zn",
    "ZN2": "Zn",
    "MN": "Mn",
    "MN3": "Mn",
    "FE": "Fe",
    "FE2": "Fe",
    "CU": "Cu",
    "CU1": "Cu",
    "NI": "Ni",
    "NI2": "Ni",
    "CO": "Co",
    "3CO": "Co",
    "CD": "Cd",
    "HG": "Hg",
    "PB": "Pb",
    "PT": "Pt",
    "AU": "Au",
    "AG": "Ag",
    "SR": "Sr",
    "BA": "Ba",
    "MO": "Mo",
    "W": "W",
    "V": "V",
    "CR": "Cr",
    "SM": "Sm",
    "EU": "Eu",
    "GD": "Gd",
    "YB": "Yb",
    "TB": "Tb",
    "LA": "La",
    "CE": "Ce",
}

#: Cryoprotectants and polyols.
_CRYOPROTECTANTS = frozenset(
    [
        "GOL",
        "EDO",
        "MPD",
        "MRD",
        "PGO",
        "PGR",
        "PDO",
        "BU3",
        "BUD",
        "PEG",
        "PGE",
        "PG4",
        "1PE",
        "2PE",
        "P6G",
        "PE4",
        "PE8",
        "12P",
        "15P",
        "XPE",
        "7PE",
        "SUC",
        "TRE",
        "DIO",
        "DOX",
    ]
)

#: Organic solvents from the crystallization cocktail.
_SOLVENTS = frozenset(["DMS", "DMF", "MOH", "EOH", "IPA", "ACY", "TFA"])

#: Buffer components.
_BUFFERS = frozenset(["EPE", "TRS", "MES", "BTB", "BIS", "CAC", "MPO", "NHE", "CXS", "TAM"])

#: Precipitants and cocktail ions. All polyatomic, so the monatomic gate that
#: guards :data:`BULK_ION_RESIDUES` cannot apply to them.
_PRECIPITANTS = frozenset(
    [
        "SO4",
        "PO4",
        "POP",
        "NO3",
        "CO3",
        "BCT",
        "ACT",
        "FMT",
        "CIT",
        "FLC",
        "TLA",
        "MLA",
        "MLI",
        "OXL",
        "SIN",
        "SCN",
        "AZI",
        "NH4",
        "TMA",
    ]
)

#: Reducing agents and His-tag elution reagents.
_REDUCTANTS = frozenset(["BME", "DTT", "DTU", "IMD"])

#: Crystallization additives: cryoprotectants, solvents, buffers, precipitants
#: and reductants. These are artefacts of how the crystal was grown, not of the
#: biology, and left in place they are exactly what antechamber wastes AM1-BCC
#: cycles on - with meaningless results, since X-ray additives carry no
#: hydrogens.
ADDITIVE_RESIDUES = _CRYOPROTECTANTS | _SOLVENTS | _BUFFERS | _PRECIPITANTS | _REDUCTANTS
