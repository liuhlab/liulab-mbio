"""What the working vector's build reads off its parent, and what it refuses.

The build itself runs under `pixi run examples-check`, which is a step of the `docs` CI job
rather than of the gate. These are the parent's side of that promise — the eight sites the
build takes out are where it says they are — and the two guards that stand between a stated
base and the wrong record.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from mbio.enzymes import get_enzyme
from mbio.sequence import SequenceRecord
from mbio.sites import find_sites
from mbio.translate import translate

if TYPE_CHECKING:
    from types import ModuleType

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "build_working_vector", REPO / "scripts/build_working_vector.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_parent_carries_the_eight_sites_the_build_takes_out(plvx: SequenceRecord) -> None:
    """The six BsaI and two BsmBI sites, where section 7 of the research note counts them."""
    found = {
        name: sorted(one.span.start for one in find_sites(plvx, get_enzyme(name)))
        for name in ("BsaI", "BsmBI")
    }

    assert found == {"BsaI": [454, 3724, 5730, 6472, 7112, 8720], "BsmBI": [3876, 5636]}


def test_the_build_leaves_no_site_and_marks_every_base_it_changed(plvx: SequenceRecord) -> None:
    """Eight bases go, eight marks arrive, and both proteins read as the parent reads them."""
    record, report = _script().rebuild(plvx)

    assert (len(record), record.topology) == (len(plvx), "circular")
    assert not find_sites(record, [get_enzyme("BsaI"), get_enzyme("BsmBI")])
    marks = [one for one in record.features if one.name.endswith("site removed")]
    assert len(marks) == len(report.changes) + len(_script().EDITS) == 8
    for name in ("PuroR", "AmpR"):
        before = next(one for one in plvx.features if one.name == name)
        after = next(one for one in record.features if one.name == name)
        assert translate(record.extract(after)) == translate(plvx.extract(before))


def test_a_record_too_short_for_a_stated_position_is_refused(plvx: SequenceRecord) -> None:
    """The containment check: a position is held to the record's length, never to `fits`.

    `SequenceRecord.fits` admits an end past the length, because that is how
    `docs/adr/0001-coordinates.md` spells a span across the origin of a circular record. A
    stated base is not such a span, so the build compares against the length itself.
    """
    short = SequenceRecord(plvx.sequence[:1000], topology="circular", name="short")

    with pytest.raises(ValueError, match="BsaI at 3725: the record is 1000 bases long"):
        _script().rebuild(short)


def test_a_parent_reading_another_base_where_one_is_stated_is_refused(
    plvx: SequenceRecord,
) -> None:
    """A different deposit is refused rather than edited: the build measured this one."""
    other = SequenceRecord(
        f"{plvx.sequence[:457]}T{plvx.sequence[458:]}", topology="circular", name="other"
    )

    with pytest.raises(ValueError, match="position 457 reads 'T' and this build states 'C'"):
        _script().rebuild(other)
