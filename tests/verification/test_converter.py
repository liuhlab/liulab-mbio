"""A clone's verification handed to molviz: its reads, disagreements and traces, held to its check."""

import dataclasses
from pathlib import Path

import molviz
import pytest

from mbio.io import read_record
from mbio.plot.converter import molviz_record
from mbio.sequence import Feature, Segment, SequenceRecord, reverse_complement
from mbio.verification import Channels, SequencingResult, read_channels, read_trace, regions
from mbio.verification import verify as verified
from mbio.verification.cli import read_result
from mbio.verification.converter import molviz_verification

EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "pUC19-GFP"

#: Where a substitution is planted in a trace's own bases, counted into its trusted bases.
PLANTED = 300

_OTHER = {"A": "C", "C": "G", "G": "T", "T": "A"}


def test_the_docs_clone_converts_to_the_reads_and_disagreements_molviz_shows() -> None:
    product = read_record(EXAMPLE / "product.dna")
    consensus = read_result(EXAMPLE / "verification" / "consensus.fasta")
    made = verified(product, (consensus,), regions(product))
    converted = molviz_verification(product, (consensus,), made)
    molviz.check(molviz_record(product) | converted)
    # As molviz's record reference draws this clone: one change in a failing region, red, and one
    # outside every region, grey.
    assert converted == {
        "reads": [{"name": "consensus.fasta", "start": 0, "end": 3347, "strand": 1}],
        "disagreements": [
            {
                "name": "696 A>G",
                "start": 695,
                "end": 696,
                "bases": "G",
                "reads": ["consensus.fasta"],
                "color": "#cc3311",
            },
            {
                "name": "2687 C>T",
                "start": 2686,
                "end": 2687,
                "bases": "T",
                "reads": ["consensus.fasta"],
                "color": "#bbbbbb",
            },
        ],
    }


@pytest.fixture(scope="module")
def trace(data_dir: Path) -> tuple[SequencingResult, Channels, str]:
    """A 3730xl read, its channels, and its own trusted bases with one base changed."""
    path = data_dir / "3730.ab1"
    result = read_trace(path)
    first, last = result.trusted_span
    bases = list(result.bases[first:last])
    bases[PLANTED] = _OTHER[bases[PLANTED]]
    return result, read_channels(path), "".join(bases)


def test_two_results_of_one_name_become_two_reads_told_apart_by_a_suffix(
    trace: tuple[SequencingResult, Channels, str],
) -> None:
    """Channels and disagreements pair with results by position, as checks and placements do."""
    result, channels, own = trace
    halved = Channels(
        {base: tuple(value // 2 for value in scan) for base, scan in channels.scans.items()},
        channels.peaks,
    )
    record = SequenceRecord(own, name="own bases")
    # The second reads the planted base as the record has it, so only the first shows it.
    planted = result.trusted_span[0] + PLANTED
    agreeing = own[PLANTED].join((result.bases[:planted], result.bases[planted + 1 :]))
    both = (result, dataclasses.replace(result, bases=agreeing))
    made = verified(record, both, ())
    converted = molviz_verification(record, both, made, channels=(channels, halved))
    molviz.check(molviz_record(record) | converted)
    names = [result.name, f"{result.name} (2)"]
    assert [read["name"] for read in converted["reads"]] == names
    assert [one["reads"] for one in converted["disagreements"]] == [names[:1]]
    assert [(track["name"], track["curves"]["A"]) for track in converted["tracks"]] == [
        (names[0], list(channels.scans["A"])),
        (names[1], list(halved.scans["A"])),
    ]
    with pytest.raises(ValueError, match="one for each"):
        molviz_verification(record, both, made, channels=(channels,))


@pytest.mark.parametrize("reverse", [False, True])
def test_a_traces_tallest_curve_at_each_peak_is_the_top_strand_base_under_it(
    trace: tuple[SequencingResult, Channels, str], reverse: bool
) -> None:
    result, channels, own = trace
    record = SequenceRecord(reverse_complement(own) if reverse else own, name="own bases")
    made = verified(record, (result,), (Feature("read", "misc_feature", (Segment(0, len(own)),)),))
    converted = molviz_verification(record, (result,), made, channels=(channels,))
    molviz.check(molviz_record(record) | converted)
    [track] = converted["tracks"]
    assert track["name"] == result.name
    assert set(track["curves"]) == set("ACGT")
    planted = len(own) - 1 - PLANTED if reverse else PLANTED
    for sample, at in track["peaks"]:
        if at != planted:
            tallest = max("ACGT", key=lambda base: track["curves"][base][sample])
            assert tallest == record.sequence[at]
