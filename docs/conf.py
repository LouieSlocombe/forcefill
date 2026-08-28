"""Sphinx configuration for the forcefill documentation.

The docs build is pip-only and installs forcefill with ``--no-deps``: two of its
runtime dependencies cannot be installed from PyPI at all (``openmmforcefields``
stops at 0.15.1 there, below the 0.16 floor, and ``parmed`` ships no wheels), so
autodoc imports the package against mocks instead. See ``autodoc_mock_imports``
below, and .readthedocs.yaml for the other half of the arrangement.
"""

from __future__ import annotations

from importlib import metadata

# -- Project information -----------------------------------------------------

project = "forcefill"
author = "Louie Slocombe"
copyright = "2026, Louie Slocombe"

# The version is static in pyproject.toml and read back through the installed
# metadata, exactly as forcefill/__init__.py does it - never duplicated here.
try:
    release = metadata.version("forcefill")
except metadata.PackageNotFoundError:  # uninstalled checkout
    release = "0.0.0+unknown"
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    # The docstrings are Google-style (ruff enforces `convention = "google"`),
    # with reST roles inside the prose. napoleon turns the sections into field
    # lists; the roles are already reST and need nothing.
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Signatures here run to 22 parameters; unqualified names keep them readable, and
# api/constants.md relies on this to document the public constants through the
# private modules that define them without printing the private path.
add_module_names = False

# -- Autodoc -----------------------------------------------------------------

# openmm is deliberately absent: forcefill/checks.py evaluates
# `unit.kilojoule_per_mole / unit.nanometer` at import time, and Sphinx's mock
# objects implement attribute access and calls but not division, so a mocked
# openmm turns that line into a TypeError. It is installed for real instead
# (docs/requirements.txt). Mocking the top-level `openff` covers openff.toolkit
# and everything under it.
autodoc_mock_imports = ["parmed", "rdkit", "openff", "openmmforcefields"]

# Kept in the signature rather than moved into the parameter descriptions.
# "description" makes Sphinx evaluate the annotations through
# typing.get_type_hints(), and forcefill.smirnoff, forcefill.espaloma and
# forcefill._pipeline import their annotation types only under `if
# TYPE_CHECKING:` - so evaluation raises NameError, which fail_on_warning turns
# into a failed build. Every module has `from __future__ import annotations`, so
# signature mode stringifies the source text without evaluating anything.
autodoc_typehints = "signature"

# Render defaults as the names the docstrings already cross-reference
# (`base_forcefield=DEFAULT_BASE_FORCEFIELD`), not as expanded values.
autodoc_preserve_defaults = True

autodoc_member_order = "bysource"

napoleon_google_docstring = True
napoleon_numpy_docstring = False

# Left off on purpose. Mocked third-party annotations (parmed.modeller.
# ResidueTemplate, openff.toolkit.Molecule) have no resolvable targets, so
# nitpicky mode plus fail_on_warning would fail every build.
nitpicky = False

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "openmm": ("https://docs.openmm.org/latest/api-python/", None),
    "openff.toolkit": ("https://docs.openforcefield.org/projects/toolkit/en/stable/", None),
    "parmed": ("https://parmed.github.io/ParmEd/html/", None),
}

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist"]
myst_heading_anchors = 3

# -- HTML output -------------------------------------------------------------

html_theme = "furo"
html_title = f"forcefill {release}"
html_theme_options = {
    "source_repository": "https://github.com/LouieSlocombe/forcefill/",
    "source_branch": "main",
    "source_directory": "docs/",
}


# -- Public-path cross-references --------------------------------------------


def _resolve_public_alias(app, env, node, contnode):
    """Resolve ``forcefill.<name>`` to wherever that name is actually documented.

    The docstrings reference the public import path - ``:func:`~forcefill.
    build_forcefield_xml```, ``:data:`~forcefill.CHARMM_BASE_FORCEFIELD``` - because
    that is how the package is used. autodoc registers each object under the
    module that *defines* it (``forcefill.structure.build_forcefield_xml``),
    because that is where the source link and the ``#:`` attribute comments come
    from. Without this hook the ~28 public-path references in the docstrings
    would silently render as plain text instead of links.

    Only an unambiguous single match is redirected, so a name documented in two
    places is left alone rather than linked to the wrong one.
    """
    from sphinx.util.nodes import make_refnode

    if node.get("refdomain") != "py":
        return None
    prefix, _, name = node.get("reftarget", "").rpartition(".")
    if prefix != "forcefill" or not name:
        return None

    matches = [
        (fullname, entry)
        for fullname, entry in env.domains["py"].objects.items()
        if fullname.startswith("forcefill.") and fullname.rpartition(".")[2] == name
    ]
    if len(matches) != 1:
        return None

    fullname, entry = matches[0]
    return make_refnode(app.builder, node["refdoc"], entry.docname, entry.node_id, contnode, fullname)


def setup(app):
    app.connect("missing-reference", _resolve_public_alias)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
