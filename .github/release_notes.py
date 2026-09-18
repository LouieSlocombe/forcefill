"""Extract one version's section from ``CHANGELOG.md`` as GitHub release notes.

Used by the ``github-release`` job in ``.github/workflows/release.yml``. The notes
are the changelog section rather than a generated commit list because the section
is already written, and saying what changed is what the changelog is for.

    python .github/release_notes.py 1.0.0 > release-notes.md
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def section_for(version: str, changelog: str) -> str:
    """Return the body of the ``## [version]`` section, without its heading.

    The section ends at the next ``## [`` heading, at the block of link-reference
    definitions the file ends with, or at end of file - whichever comes first.
    Without the second of those, the last section in the file would carry the
    link definitions into the release notes.
    """
    match = re.search(
        rf"^## \[{re.escape(version)}\].*?$(.*?)(?=^## \[|^\[[^\]]+\]:\s|\Z)",
        changelog,
        re.MULTILINE | re.DOTALL,
    )
    if match is None:
        sys.exit(f"CHANGELOG.md has no '## [{version}]' section")
    return match.group(1).strip()


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} <version>")
    print(section_for(sys.argv[1], (ROOT / "CHANGELOG.md").read_text()))


if __name__ == "__main__":
    main()
