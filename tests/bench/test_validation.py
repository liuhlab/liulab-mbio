"""Colony PCR and Sanger reads pinned on the pUC19 and GFP fixtures.

The source is `docs/research/primer-design-and-pcr.md` for the bands and the read geometry.
"""

from itertools import pairwise

import pytest

from liulab_mbio import edits
from liulab_mbio.bench.gels import LADDER_100_BP
from liulab_mbio.bench.validation import (
    COLONY_ALLOWANCE,
    COLONY_FLANK,
    CORRECT_CLONE,
    EMPTY_CLONE,
    JUNCTION_OFFSET,
    REVERSED_CLONE,
    SANGER_ALLOWANCE,
    SANGER_FLANK,
    ColonyCheck,
    colony_pcr_check,
    sanger_primers,
)
from liulab_mbio.sequence import Primer, SequenceRecord, reverse_complement

#: The multiple cloning site of the fixture, 0-based and half-open, and the GFP that replaces it.
MCS = (395, 452)
#: Addgene's 23-mer M13/pUC pair, which the note recommends over the 17-mers for colony PCR.
M13_FORWARD = Primer("M13/pUC Forward", "CCCAGTCACGACGTTGTAAAACG")
M13_REVERSE = Primer("M13/pUC Reverse", "AGCGGATAACAATTTCACACAGG")


@pytest.fixture(scope="module")
def product(puc19: SequenceRecord, gfp: SequenceRecord) -> SequenceRecord:
    """pUC19 with GFP in place of the whole multiple cloning site."""
    edited, _ = edits.replace(puc19, *MCS, gfp.sequence)
    return edited


@pytest.fixture(scope="module")
def junctions(gfp: SequenceRecord) -> tuple[int, int]:
    return (MCS[0], MCS[0] + len(gfp))


@pytest.fixture(scope="module")
def three_junctions(gfp: SequenceRecord) -> tuple[int, int, int]:
    """The same product read as two inserts, split inside GFP."""
    return (MCS[0], MCS[0] + 350, MCS[0] + len(gfp))


def reversed_plasmids(
    product: SequenceRecord, junctions: tuple[int, ...]
) -> tuple[SequenceRecord, ...]:
    """The product with each insert, in turn, the other way round.

    These inserts were pasted in with no overhang, so turning the top strand in place is exact.
    """
    made = []
    for first, last in pairwise(junctions):
        flipped, _ = edits.replace(
            product, first, last, reverse_complement(product.sequence[first:last])
        )
        made.append(flipped)
    return tuple(made)


@pytest.fixture(scope="module")
def m13_check(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> ColonyCheck:
    """The Addgene pair across one insert, and the plasmid holding that insert reversed."""
    return colony_pcr_check(
        product,
        junctions,
        vector=puc19,
        primers=(M13_FORWARD, M13_REVERSE),
        reversed_inserts=reversed_plasmids(product, junctions),
    )


@pytest.fixture(scope="module")
def directional_check(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> ColonyCheck:
    """The same pair where the ends are directional, so no plasmid holds the insert reversed."""
    return colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))


@pytest.fixture(scope="module")
def two_insert_check(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> ColonyCheck:
    """Two inserts, each with a junction primer and a plasmid holding it the other way round."""
    return colony_pcr_check(
        product,
        three_junctions,
        vector=puc19,
        flank=60,
        insert_primer=True,
        reversed_inserts=reversed_plasmids(product, three_junctions),
    )


def bands(check: ColonyCheck, name: str) -> tuple[int, ...]:
    return next(clone.bands_bp for clone in check.clones if clone.name == name)


def test_the_product_band_is_the_empty_band_less_what_the_insert_replaced(
    directional_check: ColonyCheck,
) -> None:
    # 137 bp empty, less the 57 bp multiple cloning site, plus 717 bp of GFP.
    assert bands(directional_check, "Correct clone") == (797,)
    assert bands(directional_check, "Empty vector") == (137,)


def test_two_flanking_primers_cannot_tell_the_insert_round_the_other_way(
    m13_check: ColonyCheck,
) -> None:
    assert bands(m13_check, "Reversed insert") == bands(m13_check, "Correct clone")
    assert not m13_check.tells_orientation


def test_a_clone_that_cannot_turn_round_simulates_no_reversed_insert(
    m13_check: ColonyCheck, directional_check: ColonyCheck
) -> None:
    every, directional = m13_check, directional_check
    assert [clone.name for clone in every.clones] == [CORRECT_CLONE, EMPTY_CLONE, REVERSED_CLONE]
    assert directional.clones == every.clones[:2]
    assert not directional.reversed_clones
    assert [lane.label for lane in directional.gel.lanes] == [CORRECT_CLONE, EMPTY_CLONE]


def test_a_reverse_distance_of_its_own_gives_a_reversed_insert_bands_of_its_own(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(
        product,
        junctions,
        vector=puc19,
        flank=60,
        reverse_flank=120,
        insert_primer=True,
        reversed_inserts=reversed_plasmids(product, junctions),
    )
    # 60 bases of vector before the first junction and 120 past the last, each primer free to
    # move a little either way, so the junction primer reaches the near flank in a correct clone
    # and the far one in a reversed clone.
    assert bands(check, "Correct clone") == (160, 879)
    assert bands(check, "Reversed insert") == (220, 879)
    assert bands(check, "Empty vector") == (219,)
    assert check.tells_orientation


def test_the_gel_carries_one_lane_per_clone_and_a_ladder_for_the_range(
    m13_check: ColonyCheck,
) -> None:
    gel = m13_check.gel
    assert [lane.label for lane in gel.lanes] == [
        "Correct clone",
        "Empty vector",
        "Reversed insert",
    ]
    assert gel.ladder is LADDER_100_BP
    assert m13_check.agarose_percent == 2.0


def test_the_annealing_temperature_and_extension_come_from_the_polymerase(
    directional_check: ColonyCheck,
) -> None:
    # OneTaq: the 23-mer M13/pUC pair anneals at 51.5 C, and 797 bp extends for a minute.
    assert directional_check.annealing_temperature == pytest.approx(51.5, abs=0.1)
    assert directional_check.extension_seconds == 60


def test_sanger_primers_read_from_outside_each_junction(
    product: SequenceRecord, junctions: tuple[int, int], gfp: SequenceRecord
) -> None:
    forward, reverse = sanger_primers(product, junctions)
    # Genewiz's 100 bases at the closest, and no further out than the allowance past it.
    for read in (forward, reverse):
        assert SANGER_FLANK <= read.distance_bp <= SANGER_FLANK + SANGER_ALLOWANCE
    assert forward.read_bp == forward.distance_bp + len(gfp)
    assert reverse.read_bp == reverse.distance_bp + len(gfp)


def test_colony_pcr_check_refuses_a_junction_pair_it_cannot_place(
    product: SequenceRecord, puc19: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="junction"):
        colony_pcr_check(product, (395,), vector=puc19)


def test_a_junction_primer_is_designed_for_every_insert(two_insert_check: ColonyCheck) -> None:
    assert len(two_insert_check.primers) == 4
    assert all(report.status != "fail" for report in two_insert_check.reports)
    assert [clone.name for clone in two_insert_check.clones] == [
        "Correct clone",
        "Empty vector",
        "Reversed insert 1",
        "Reversed insert 2",
    ]


def test_the_correct_clone_shows_one_band_reading_each_junction(
    two_insert_check: ColonyCheck,
) -> None:
    # 60 bases of vector then 100 into each insert, each within its allowance of that, and the
    # flanking pair across the whole span.
    assert bands(two_insert_check, "Correct clone") == (160, 492, 837)
    assert bands(two_insert_check, "Empty vector") == (177,)


def test_each_reversed_plasmid_gets_a_lane_of_its_own(two_insert_check: ColonyCheck) -> None:
    correct = bands(two_insert_check, "Correct clone")
    assert bands(two_insert_check, "Reversed insert 1") != correct
    assert bands(two_insert_check, "Reversed insert 2") != correct
    assert two_insert_check.tells_orientation


def test_sanger_primers_read_from_outside_the_whole_inserted_span(
    product: SequenceRecord, three_junctions: tuple[int, ...], gfp: SequenceRecord
) -> None:
    forward, reverse = sanger_primers(product, three_junctions)
    assert forward.read_bp == forward.distance_bp + len(gfp)
    assert reverse.read_bp == reverse.distance_bp + len(gfp)


def test_an_insert_too_short_for_a_junction_primer_loses_it_rather_than_refusing_the_check(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    # A linker or a tag has no room to anneal a primer 100 bases inside it, whatever made it, so
    # the middle insert here is shorter than the offset and the allowance together.
    first, split, last = three_junctions
    short = (first, split - (JUNCTION_OFFSET - COLONY_ALLOWANCE), split, last)
    check = colony_pcr_check(
        product,
        short,
        vector=puc19,
        flank=60,
        insert_primer=True,
        reversed_inserts=reversed_plasmids(product, short),
    )
    # Numbered for the insert each reads out of, so the one with no room leaves a gap.
    assert [primer.name for primer in check.primers] == [
        "Colony PCR forward",
        "Colony PCR reverse",
        "Junction reverse 1",
        "Junction reverse 3",
    ]
    # The flanking pair still reads across it, and the gel says plainly that it cannot tell it
    # turned round rather than claiming a lane it does not have.
    correct = bands(check, "Correct clone")
    assert bands(check, "Reversed insert 2") == correct
    assert bands(check, "Reversed insert 1") != correct
    assert bands(check, "Reversed insert 3") != correct
    assert not check.tells_orientation


def test_an_insert_of_exactly_the_offset_and_the_allowance_keeps_its_junction_primer(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    # Fewer bases than the two together is what has no room, so an insert of exactly that many
    # has some: its whole 5' window lies inside it, short of the far junction.
    first, last = junctions
    room = JUNCTION_OFFSET + COLONY_ALLOWANCE
    check = colony_pcr_check(
        product, (first, first + room, last), vector=puc19, flank=60, insert_primer=True
    )
    assert [primer.name for primer in check.primers] == [
        "Colony PCR forward",
        "Colony PCR reverse",
        "Junction reverse 1",
        "Junction reverse 2",
    ]
    one_shorter = colony_pcr_check(
        product, (first, first + room - 1, last), vector=puc19, flank=60, insert_primer=True
    )
    assert [primer.name for primer in one_shorter.primers] == [
        "Colony PCR forward",
        "Colony PCR reverse",
        "Junction reverse 2",
    ]


def test_a_junction_set_whose_last_position_passes_the_length_reads_across_the_origin(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    # The same two junctions read the other way round, so the span between them is the vector
    # and the last passes the length: what a set in insert order does across the origin.
    wrapping = (junctions[1], junctions[0] + len(product))
    span = wrapping[1] - wrapping[0]
    check = colony_pcr_check(product, wrapping, vector=puc19, flank=COLONY_FLANK)
    forward, reverse = sanger_primers(product, wrapping)
    # One band across the origin, a flank either side of it within the allowance, and not the
    # GFP a sorted set would have read instead.
    (band,) = bands(check, CORRECT_CLONE)
    assert span + 2 * (COLONY_FLANK - COLONY_ALLOWANCE) <= band
    assert band <= span + 2 * (COLONY_FLANK + COLONY_ALLOWANCE)
    assert forward.read_bp == forward.distance_bp + span
    assert reverse.read_bp == reverse.distance_bp + span
