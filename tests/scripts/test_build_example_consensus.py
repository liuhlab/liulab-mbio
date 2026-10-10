"""What the example consensus comes to against the product it was made from.

The sequence verification page explains three outcomes from this one result. The build itself
runs under `pixi run examples-check`; this is the other half, that the result still gives them.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from mbio.io import read_record
from mbio.verification import SequencingResult, regions, verify

if TYPE_CHECKING:
    from types import ModuleType

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "build_example_consensus", REPO / "scripts/build_example_consensus.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_junctions_pass_gfp_fails_and_the_ampr_change_carries_no_verdict() -> None:
    script = _script()
    product = read_record(script.PRODUCT)
    bases, _ = script.planted(product)

    made = verify(product, [SequencingResult("consensus", bases)], regions(product))

    assert [(one.name, one.status) for one in made.checks] == [
        ("ATGA junction", "pass"),
        ("GFP insert", "fail"),
        ("TGGC junction", "pass"),
    ]
    (outside,) = [one for one in made.disagreements if not one.regions]
    assert outside.features == ("AmpR",)
    assert not made.verified
