"""Check that the tag, ``pyproject.toml`` and ``CHANGELOG.md`` agree on a version.

Run by the ``guard`` job in ``.github/workflows/release.yml`` before anything is
built, because a version uploaded to PyPI can never be re-uploaded: the wrong
number published is the wrong number forever.

Prints the version on stdout for the workflow to capture; exits non-zero with the
disagreement on stderr. ``GITHUB_REF`` is only checked when it names a tag, so a
manual dry run validates the changelog without needing one.

Needs Python >= 3.11 for ``tomllib`` - below the package's own 3.10 floor on
purpose, because this only ever runs on the 3.12 the release workflow pins.

    python .github/check_release_version.py
"""

from __future__ import annotations

import os
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]

    ref = os.environ.get("GITHUB_REF", "")
    if ref.startswith("refs/tags/"):
        tag = ref.removeprefix("refs/tags/")
        if tag != f"v{version}":
            sys.exit(f"tag {tag} does not match the pyproject.toml version {version} (expected v{version})")

    changelog = (ROOT / "CHANGELOG.md").read_text()
    if not re.search(rf"^## \[{re.escape(version)}\]", changelog, re.MULTILINE):
        sys.exit(f"CHANGELOG.md has no '## [{version}]' section - add one before tagging")

    print(version)


if __name__ == "__main__":
    main()
