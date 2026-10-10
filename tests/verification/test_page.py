"""The page sequence verification writes for one clone: its verdict, a table row per region, the
map, and a close-up of each disagreement with each Sanger result's trace under its bases."""

import random
import re
import socket
from pathlib import Path

import pytest

from mbio.io import read_record
from mbio.plot import sequence_view
from mbio.sequence import Feature, Segment, SequenceRecord, reverse_complement
from mbio.verification import page
from mbio.verification.cli import read_result
from mbio.verification.judge import CONSENSUS_CAVEAT, Verification, regions, verify
from mbio.verification.result import SequencingResult
from mbio.verification.trace import Channels, read_channels, read_trace

from ..html import parse

#: The worked Golden Gate example, and the consensus made from its product with two changes.
EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "pUC19-GFP"

#: Where a substitution is planted, counted into the trusted bases.
PLANTED = 300

_OTHER = {"A": "C", "C": "G", "G": "T", "T": "A"}


@pytest.fixture(autouse=True)
def offline(monkeypatch: pytest.MonkeyPatch) -> None:
    """Nothing reaches the network: opening a connection fails the test."""

    def refuse(*_: object, **__: object) -> None:
        raise AssertionError("the page reached for the network")

    monkeypatch.setattr(socket.socket, "connect", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)


@pytest.fixture(scope="module")
def trace(data_dir: Path) -> tuple[SequencingResult, Channels]:
    """A 3730xl read and the channels its bases were called from."""
    path = data_dir / "3730.ab1"
    return read_trace(path), read_channels(path)


@pytest.fixture(scope="module")
def worked() -> tuple[page.Page, Verification, tuple[Feature, ...]]:
    """The worked clone's page: its consensus changed once in GFP and once in AmpR."""
    product = read_record(EXAMPLE / "product.dna")
    consensus = read_result(EXAMPLE / "verification" / "consensus.fasta")
    judged = regions(product)
    made = verify(product, (consensus,), judged)
    return page.draw_page(product, (consensus,), made), made, judged


def test_the_page_holds_the_verdict_a_row_per_region_and_a_close_up_per_disagreement(
    worked: tuple[page.Page, Verification, tuple[Feature, ...]],
) -> None:
    laid, made, judged = worked
    [inside, outside] = made.disagreements
    [failing] = inside.regions

    assert not laid.verified
    assert [row.where for row in laid.rows] == [
        *(one.name for one in judged),
        "AmpR, outside every region",
        "consensus.fasta",
    ]
    verdicts = {row.where: row.verdict for row in laid.rows}
    assert verdicts[failing] == "fail"
    assert verdicts["AmpR, outside every region"] is None
    assert laid.rows[-1].detail.endswith(CONSENSUS_CAVEAT)
    assert [one.disagreement for one in laid.close_ups] == [inside, outside]
    # On the map, each disagreement is lit, coloured by its region's verdict, or grey outside.
    assert laid.map is not None
    names = [one.name for one in laid.close_ups]
    assert laid.map.highlight == tuple(names)
    colours = {
        arrow.item.name: arrow.span.color
        for arrow in laid.map.layout.arrows
        if arrow.item.name in names
    }
    assert colours == {names[0]: "#cc3311", names[1]: "#bbbbbb"}
    written = parse(laid.html)
    assert [heading.text for heading in written.find_all("h1")] == ["Not verified"]
    assert len(written.find_all("tbody")[0].find_all("tr")) == len(laid.rows)
    assert len(written.find_all("figure", cls="close-up")) == 2


def _own_bases(result: SequencingResult) -> list[str]:
    first, last = result.trusted_span
    return list(result.bases[first:last])


def _laid(
    record: SequenceRecord, result: SequencingResult, channels: Channels | None = None
) -> page.Page:
    whole = Feature("read", "misc_feature", (Segment(0, len(record)),))
    verification = verify(record, (result,), (whole,))
    return page.draw_page(record, (result,), verification, channels=(channels,) if channels else ())


def test_a_sanger_close_up_has_its_trace_under_it_each_peak_over_its_base(
    trace: tuple[SequencingResult, Channels],
) -> None:
    result, channels = trace
    bases = _own_bases(result)
    bases[PLANTED] = _OTHER[bases[PLANTED]]
    laid = _laid(SequenceRecord("".join(bases), name="own bases"), result, channels)

    assert not laid.verified
    [close] = laid.close_ups
    assert (close.disagreement.kind, close.disagreement.start) == ("substitution", PLANTED)
    view = close.drawing.sequence_view
    assert view is not None
    [row] = view.rows
    assert row.start <= PLANTED < row.end
    [strip] = row.strips
    first = strip.track.curves[0]
    drawn = {
        sample: point
        for stroke in strip.strokes
        if stroke.curve is first
        for sample, point in zip(stroke.samples, stroke.points, strict=True)
    }
    # The record is the read's own trusted bases, so the read's base `at` is the record's.
    offset = result.trusted_span[0]
    for at in range(row.start, row.end):
        x = drawn[channels.peaks[offset + at]].x
        assert x == pytest.approx((at - row.start + 0.5) * sequence_view.CELL)
    # The read's arrow is marked at the planted base, and the page writes the strip.
    [mark] = row.mismatches
    assert mark.box.x // sequence_view.CELL == PLANTED - row.start
    assert len(parse(laid.html).find_all("g", cls="track")) == 1


def test_two_traces_of_one_name_each_draw_their_own_channels(
    trace: tuple[SequencingResult, Channels],
) -> None:
    """Channels pair with results by position, as checks and placements do, not by name."""
    result, channels = trace
    halved = Channels(
        {base: tuple(value // 2 for value in scan) for base, scan in channels.scans.items()},
        channels.peaks,
    )
    bases = _own_bases(result)
    bases[PLANTED] = _OTHER[bases[PLANTED]]
    record = SequenceRecord("".join(bases), name="own bases")
    both = (result, result)
    made = verify(record, both, ())
    laid = page.draw_page(record, both, made, channels=(channels, halved))
    [close] = laid.close_ups
    view = close.drawing.sequence_view
    assert view is not None
    [row] = view.rows
    assert [strip.track.curves[0].values for strip in row.strips] == [
        channels.scans["A"],
        halved.scans["A"],
    ]
    with pytest.raises(ValueError, match="one for each"):
        page.draw_page(record, both, made, channels=(channels,))


@pytest.mark.parametrize("reverse", [False, True])
def test_each_peak_s_tallest_curve_is_the_one_for_the_top_strand_base_over_it(
    trace: tuple[SequencingResult, Channels], reverse: bool
) -> None:
    """A read along the bottom strand draws each channel as the base it pairs with on the top."""
    result, channels = trace
    bases = _own_bases(result)
    bases[PLANTED] = _OTHER[bases[PLANTED]]
    own = "".join(bases)
    record = SequenceRecord(reverse_complement(own) if reverse else own, name="own bases")
    [close] = _laid(record, result, channels).close_ups
    view = close.drawing.sequence_view
    assert view is not None
    [row] = view.rows
    [strip] = row.strips
    curves = strip.track.curves
    for sample, base in strip.track.peaks:
        if row.start <= base < row.end and base != close.disagreement.start:
            tallest = max(range(4), key=lambda one: curves[one].values[sample])
            assert "ACGT"[tallest] == record.sequence[base]


def test_an_insertion_draws_on_the_map_and_in_its_close_up(
    trace: tuple[SequencingResult, Channels],
) -> None:
    result, _ = trace
    bases = _own_bases(result)
    del bases[PLANTED]
    laid = _laid(SequenceRecord("".join(bases), name="one base short"), result)

    [close] = laid.close_ups
    one = close.disagreement
    assert (one.kind, len(one.bases)) == ("insertion", 1)
    assert laid.map is not None
    assert [
        label.item.label for label in laid.map.layout.labels if label.item.kind == "insertion"
    ] == ["+1 bp"]
    view = close.drawing.sequence_view
    assert view is not None
    marks = [mark for row in view.rows for mark in row.insertions]
    assert [mark.item.name for mark in marks] == [close.name]


def test_a_page_where_every_result_is_gated_holds_only_the_gate_rows_and_not_verified(
    trace: tuple[SequencingResult, Channels],
) -> None:
    result, _ = trace
    rng = random.Random(1)
    unrelated = SequenceRecord("".join(rng.choice("ACGT") for _ in range(1200)), name="other")
    laid = _laid(unrelated, result)

    assert [(row.where, row.verdict) for row in laid.rows] == [(result.name, "fail")]
    assert laid.rows[0].detail.startswith("does not read as this product")
    assert laid.map is None
    assert laid.close_ups == ()
    written = parse(laid.html)
    assert [heading.text for heading in written.find_all("h1")] == ["Not verified"]
    assert len(written.find_all("tr")) == 2
    assert not written.find_all("svg")


def test_the_page_is_one_file_that_loads_nothing(
    trace: tuple[SequencingResult, Channels], tmp_path: Path
) -> None:
    result, channels = trace
    bases = _own_bases(result)
    bases[PLANTED] = _OTHER[bases[PLANTED]]
    laid = _laid(SequenceRecord("".join(bases), name="own bases"), result, channels)

    written = laid.write(tmp_path / "made" / "verification.html")
    html = written.read_text(encoding="utf-8")
    assert written.parent.name == "made"
    assert not re.search(r"""(src|href)\s*=\s*["']?(https?:)?//""", html)
    assert not re.search(r"url\((?!data:)", html)
    document = parse(html)
    assert not document.find_all("script")
    assert not document.find_all("link")
    assert "@font-face" in html
