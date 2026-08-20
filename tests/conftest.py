"""Fixtures shared across the hermetic test modules."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("openmm")
pytest.importorskip("parmed")

from forcefill import amber
from forcefill._spec import PathLike
from tests.helpers import DATA


@pytest.fixture
def fake_ambertools(monkeypatch: pytest.MonkeyPatch) -> dict[str, list[dict[str, Any]]]:
    """Replace the AmberTools layer with fakes that install the committed fixtures.

    Patches :mod:`forcefill.amber` rather than any one caller: the pipeline
    reaches AmberTools through the module (``amber.run_antechamber(...)``, never
    a from-import), so one patch covers both entry points.

    ``shutil.which`` is stubbed for the two executables too, so a code path that
    quietly goes looking for a real one fails here rather than passing on
    whichever machine happens to have AmberTools installed.
    """
    calls: dict[str, list[dict[str, Any]]] = {"antechamber": [], "parmchk2": []}

    def fake_antechamber(input_file: PathLike, output_mol2: PathLike, residue_name: str, **kwargs: object) -> str:
        calls["antechamber"].append(
            {"input": str(input_file), "output": str(output_mol2), "residue": residue_name, **kwargs}
        )
        out = Path(output_mol2)
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(DATA / "methanol.mol2", out)
        return str(out)

    def fake_parmchk2(
        input_mol2: PathLike, output_frcmod: PathLike, atom_type: str = "gaff2", timeout: float | None = None
    ) -> str:
        calls["parmchk2"].append(
            {"input": str(input_mol2), "output": str(output_frcmod), "atom_type": atom_type, "timeout": timeout}
        )
        shutil.copyfile(DATA / "methanol.frcmod", output_frcmod)
        return str(output_frcmod)

    real_which = shutil.which

    def which_without_ambertools(name: str, *args: Any, **kwargs: Any) -> str | None:
        return None if name in ("antechamber", "parmchk2") else real_which(name, *args, **kwargs)

    monkeypatch.setattr(shutil, "which", which_without_ambertools)
    monkeypatch.setattr(amber, "require_executable", lambda name: f"/fake/{name}")
    # The complete frcmod (parmchk2 -a Y) stands in for gaff2.dat.
    monkeypatch.setattr(amber, "locate_gaff_dat", lambda atom_type="gaff2": str(DATA / "methanol.frcmod"))
    monkeypatch.setattr(amber, "run_antechamber", fake_antechamber)
    monkeypatch.setattr(amber, "run_parmchk2", fake_parmchk2)
    return calls
