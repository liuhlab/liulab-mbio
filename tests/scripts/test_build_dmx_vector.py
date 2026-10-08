"""What the destination's rebuild reads off its two parents.

The rebuild itself runs under `pixi run examples-check`, which is a step of the `docs` CI job
rather than of the gate. These are the parents' side of that promise: the two tracked GenBank
files carry the bases, the topology and the one coding sequence each that
`scripts/build_dmx_vector.py` takes, so a re-conversion that loses any of it fails here rather
than in the example.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from types import ModuleType

    import pytest

    from liulab_mbio.sequence import Feature, SequenceRecord

REPO = Path(__file__).resolve().parents[2]


def _script() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "build_dmx_vector", REPO / "scripts/build_dmx_vector.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Registered before it runs: `dataclasses` resolves annotations through `sys.modules`.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _cds(record: SequenceRecord, name: str) -> Feature:
    found = [one for one in record.features if one.name == name and one.type == "CDS"]
    assert len(found) == 1, f"{name}: {len(found)} coding sequences"
    return found[0]


def test_dmx0001_carries_the_ampr_the_rebuild_replaces(dmx0001: SequenceRecord) -> None:
    """The parent, and the marker the rebuild takes out: a joined CDS on the reverse strand."""
    from liulab_mbio.sequence import Strand

    assert (len(dmx0001), dmx0001.topology) == (5839, "circular")
    ampr = _cds(dmx0001, "AmpR")
    assert ampr.strand is Strand.REVERSE
    assert [(one.start, one.end) for one in ampr.segments] == [(1524, 2316), (2316, 2385)]


def test_pcr_blunt_ii_topo_carries_the_marker_that_goes_in(
    pcr_blunt_ii_topo: SequenceRecord,
) -> None:
    """The carrier, and the coding sequence the rebuild lifts out of it."""
    from liulab_mbio.sequence import Strand

    assert (len(pcr_blunt_ii_topo), pcr_blunt_ii_topo.topology) == (3519, "circular")
    marker = _cds(pcr_blunt_ii_topo, "NeoR/KanR")
    assert marker.strand is Strand.FORWARD
    assert [(one.start, one.end) for one in marker.segments] == [(1236, 2031)]


def test_the_genbank_written_carries_each_feature_colour(tmp_path: Path) -> None:
    """A coloured record reads back coloured: the writer spells the note the reader takes."""
    from liulab_mbio.io import read_record
    from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand

    features = (
        Feature("ori", "rep_origin", (Segment(0, 12),), color="#ffff00"),
        Feature("lacI", "CDS", (Segment(16, 28),), strand=Strand.REVERSE, color="#993366"),
        Feature("joined", "misc_feature", (Segment(32, 40), Segment(48, 56)), color="#ccffcc"),
        Feature("plain", "misc_feature", (Segment(44, 48),)),
    )
    record = SequenceRecord("ACGT" * 16, name="coloured", features=features)
    path = tmp_path / "coloured.gb"
    _script().write_genbank(record, path)

    again = read_record(path)
    assert [one.color for one in again.features] == [one.color for one in features]
    assert [one.name for one in again.features] == [one.name for one in features]
    assert [[(s.start, s.end) for s in one.segments] for one in again.features] == [
        [(s.start, s.end) for s in one.segments] for one in features
    ]


def test_the_genbank_written_loses_a_segment_name_and_says_which(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """The lossy edge the writer documents, named on stderr rather than lost in silence."""
    from liulab_mbio.io import read_record
    from liulab_mbio.sequence import Feature, Segment, SequenceRecord

    promoter = Feature(
        "lac promoter",
        "promoter",
        (Segment(0, 6, name="-35"), Segment(6, 24), Segment(24, 31, name="-10")),
    )
    record = SequenceRecord("ACGT" * 16, name="named", features=(promoter,))
    path = tmp_path / "named.gb"
    _script().write_genbank(record, path)

    (again,) = read_record(path).features
    assert [(one.start, one.end) for one in again.segments] == [(0, 6), (6, 24), (24, 31)]
    assert [one.name for one in again.segments] == ["", "", ""]
    said = capsys.readouterr().err
    assert [name for name in ("-35", "-10") if f"{name!r}" in said] == ["-35", "-10"]


def test_the_definition_line_ends_in_one_full_stop_however_the_description_ends(
    tmp_path: Path,
) -> None:
    """Biopython terminates ``DEFINITION`` itself, so a description's own stop would double it."""
    from liulab_mbio.sequence import SequenceRecord

    for index, described in enumerate(("a stuffer.", "a stuffer")):
        record = SequenceRecord("ACGT" * 16, name="stuffer", notes={"Description": described})
        path = tmp_path / f"stuffer-{index}.gb"
        _script().write_genbank(record, path)

        (line,) = [one for one in path.read_text().splitlines() if one.startswith("DEFINITION")]
        assert line == "DEFINITION  a stuffer."
