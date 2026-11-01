"""The page sequence verification writes for one clone: its verdict, a table row per region, and
molviz's map and close-up of each disagreement, still, in one file that loads nothing."""

import random
import re
import socket
from pathlib import Path

import pytest

from mbio.io import read_record
from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.verification import page
from mbio.verification.cli import read_result
from mbio.verification.judge import CONSENSUS_CAVEAT, regions, verify
from mbio.verification.result import SequencingResult
from mbio.verification.trace import Channels, read_channels, read_trace

from ..html import parse
from .conftest import Handed

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


def test_the_page_holds_the_verdict_a_row_per_region_and_each_drawing_lit_in_one_offline_file(
    handed: Handed, tmp_path: Path
) -> None:
    product = read_record(EXAMPLE / "product.dna")
    consensus = read_result(EXAMPLE / "verification" / "consensus.fasta")
    judged = regions(product)
    made = verify(product, (consensus,), judged)
    laid = page.draw_page(product, (consensus,), made)
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
    assert [(one.disagreement, one.name) for one in laid.close_ups] == [
        (inside, "696 A>G"),
        (outside, "2687 C>T"),
    ]
    # The map without cut sites, each disagreement and the result lit and the product's own
    # primers not; then each close-up, 30 bases either side in one row, as molviz's how-to draws
    # this clone.
    assert [options for _, options in handed] == [
        {"kinds_off": ["cuts"], "highlight": ["696 A>G", "2687 C>T", "consensus.fasta"]},
        *(
            {
                "span": {"start": start, "end": start + 61},
                "map": "line",
                "sequence_view": True,
                "bases_per_row": 61,
                "kinds_off": ["cuts"],
                "highlight": [name, "consensus.fasta"],
            }
            for start, name in ((665, "696 A>G"), (2656, "2687 C>T"))
        ),
    ]

    written = laid.write(tmp_path / "made" / "verification.html")
    assert written.parent.name == "made"
    html = written.read_text(encoding="utf-8")
    document = parse(html)
    assert [heading.text for heading in document.find_all("h1")] == ["Not verified"]
    assert len(document.find_all("tbody")[0].find_all("tr")) == len(laid.rows)
    assert len(document.find_all("figure", cls="close-up")) == 2
    assert len(document.find_all("figure")) == 3
    assert all(len(one.find_all("svg")) == 1 for one in document.find_all("figure"))
    assert not re.search(r"""(src|href)\s*=\s*["']?(https?:)?//""", html)
    assert not re.search(r"url\((?!data:)", html)
    assert not document.find_all("script")
    assert not document.find_all("link")
    assert "@font-face" in html


def _own_bases(result: SequencingResult) -> list[str]:
    first, last = result.trusted_span
    return list(result.bases[first:last])


def _laid(record: SequenceRecord, result: SequencingResult) -> page.Page:
    whole = Feature("read", "misc_feature", (Segment(0, len(record)),))
    return page.draw_page(record, (result,), verify(record, (result,), (whole,)))


def test_two_traces_of_one_name_draw_under_the_close_up_as_two_reads_each_lit(
    trace: tuple[SequencingResult, Channels], handed: Handed
) -> None:
    """molviz refuses two reads of one name, so the page draws and lights them by suffix."""
    result, channels = trace
    bases = _own_bases(result)
    bases[PLANTED] = _OTHER[bases[PLANTED]]
    record = SequenceRecord("".join(bases), name="own bases")
    both = (result, result)
    laid = page.draw_page(record, both, verify(record, both, ()), channels=(channels, channels))

    [close] = laid.close_ups
    assert (close.disagreement.kind, close.disagreement.start) == ("substitution", PLANTED)
    names = [result.name, f"{result.name} (2)"]
    [_, (drawn, options)] = handed
    assert [track["name"] for track in drawn["tracks"]] == names
    assert options["span"] == {"start": PLANTED - 30, "end": PLANTED + 31}
    assert options["highlight"] == [close.name, *names]


def test_an_insertion_is_drawn_at_its_point_and_named_as_the_command_names_it(
    trace: tuple[SequencingResult, Channels], handed: Handed
) -> None:
    result, _ = trace
    bases = _own_bases(result)
    del bases[PLANTED]
    laid = _laid(SequenceRecord("".join(bases), name="one base short"), result)

    [close] = laid.close_ups
    one = close.disagreement
    assert (one.kind, len(one.bases)) == ("insertion", 1)
    assert close.name == one.said(len(bases)) == f"insertion after {PLANTED}"
    assert close.caption.startswith(f"insertion of {one.bases} in ")
    [disagreement] = handed[0][0]["disagreements"]
    assert (disagreement["name"], disagreement["at"]) == (close.name, PLANTED)


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
