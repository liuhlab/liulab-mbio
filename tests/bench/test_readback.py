"""Judging a construct read back well by well, and the clean-colony curve picking is sized on."""

import pytest

from mbio.bench.readback import CLEAN_COLONY_CURVE, WellVerdict, clean_colony_chance, reformat
from mbio.checks import Check
from mbio.protocol.model import Well


def test_reformatting_compacts_out_failures_and_keeps_the_uncalled():
    """Compacting out a well nobody read throws away a design that may be clean."""
    thin = Check("reads_per_well", None, 4.0, "too thin to call")
    deep = Check("reads_per_well", "pass", 400.0, "deep enough")
    kept = WellVerdict(Well("picked", "A1"), (thin,))
    wrong = Check("designed region", "fail", 7.0, "substitution at 6")
    right = Check("designed region", "pass", 7.0, "every base read, and none disagrees")
    gone = WellVerdict(Well("picked", "A2"), (deep, wrong))
    passed = WellVerdict(Well("picked", "A3"), (deep, right))
    assert (kept.status, gone.status, passed.status) == ("pass", "fail", "pass")
    assert kept.called is False
    assert [one.well.well for one in reformat((kept, gone, passed))] == ["A1", "A3"]


def test_the_clean_colony_curve_returns_its_measured_anchors():
    """Every anchor is Lund's own; between them the curve is interpolated and says so."""
    for fragments, chance in CLEAN_COLONY_CURVE:
        assert clean_colony_chance(fragments) == pytest.approx(chance)
    assert clean_colony_chance(1) == 1.0
    assert clean_colony_chance(20) == 0.0
    assert 0.846 < clean_colony_chance(4) < 0.938
