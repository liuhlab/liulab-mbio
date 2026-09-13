"""Colony PCR and Sanger reads pinned on the pUC19 and GFP fixtures.

The source is `docs/research/primer-design-and-pcr.md` for the bands and the read geometry.
"""

from pathlib import Path

import pytest

from liulab_mbio import edits
from liulab_mbio.bench.gels import LADDER_100_BP
from liulab_mbio.bench.validation import ColonyCheck, colony_pcr_check, sanger_primers
from liulab_mbio.io import read_record
from liulab_mbio.primers import amplicon_sizes
from liulab_mbio.sequence import Primer, SequenceRecord

DATA = Path(__file__).parent / "data"

#: The multiple cloning site of the fixture, 0-based and half-open, and the GFP that replaces it.
MCS = (395, 452)
#: Addgene's 23-mer M13/pUC pair, which the note recommends over the 17-mers for colony PCR.
M13_FORWARD = Primer("M13/pUC Forward", "CCCAGTCACGACGTTGTAAAACG")
M13_REVERSE = Primer("M13/pUC Reverse", "AGCGGATAACAATTTCACACAGG")


@pytest.fixture(scope="module")
def puc19() -> SequenceRecord:
    return read_record(DATA / "pUC19.dna")


@pytest.fixture(scope="module")
def gfp() -> SequenceRecord:
    return read_record(DATA / "GFP.dna")


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


def bands(check: ColonyCheck, name: str) -> tuple[int, ...]:
    return next(clone.bands_bp for clone in check.clones if clone.name == name)


def test_the_m13_pair_gives_nebs_empty_vector_band(puc19: SequenceRecord) -> None:
    assert amplicon_sizes(M13_FORWARD, M13_REVERSE, puc19) == (137,)


def test_the_product_band_is_the_empty_band_less_what_the_insert_replaced(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    # 137 bp empty, less the 57 bp multiple cloning site, plus 717 bp of GFP.
    assert bands(check, "Correct clone") == (797,)
    assert bands(check, "Empty vector") == (137,)


def test_two_flanking_primers_cannot_tell_the_insert_round_the_other_way(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    assert bands(check, "Reversed insert") == bands(check, "Correct clone")
    assert not check.tells_orientation


def test_a_junction_primer_tells_orientation_when_the_flanks_differ(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(
        product,
        junctions,
        vector=puc19,
        primers=(M13_FORWARD, M13_REVERSE),
        insert_primer=True,
    )
    assert len(check.primers) == 3
    assert bands(check, "Reversed insert") != bands(check, "Correct clone")
    assert check.tells_orientation
    assert bands(check, "Empty vector") == (137,)


def test_a_designed_pair_flanks_both_junctions(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, flank=60)
    assert bands(check, "Correct clone") == (837,)
    assert bands(check, "Empty vector") == (177,)
    assert all(report.status != "fail" for report in check.reports)


def test_the_gel_carries_one_lane_per_clone_and_a_ladder_for_the_range(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    gel = check.gel
    assert [lane.label for lane in gel.lanes] == [
        "Correct clone",
        "Empty vector",
        "Reversed insert",
    ]
    assert gel.ladder is LADDER_100_BP
    assert check.agarose_percent == 2.0


def test_the_annealing_temperature_and_extension_come_from_the_polymerase(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    # OneTaq: the 23-mer M13/pUC pair anneals at 51.5 C, and 797 bp extends for a minute.
    assert check.annealing_temperature == pytest.approx(51.5, abs=0.1)
    assert check.extension_seconds == 60


def test_sanger_primers_read_from_outside_each_junction(
    product: SequenceRecord, junctions: tuple[int, int], gfp: SequenceRecord
) -> None:
    forward, reverse = sanger_primers(product, junctions)
    assert forward.distance_bp >= 100
    assert reverse.distance_bp >= 100
    assert forward.read_bp == forward.distance_bp + len(gfp)
    assert reverse.read_bp == reverse.distance_bp + len(gfp)


def test_colony_pcr_check_refuses_a_junction_pair_it_cannot_place(
    product: SequenceRecord, puc19: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="junction"):
        colony_pcr_check(product, (395,), vector=puc19)


def test_a_junction_primer_is_designed_for_every_insert(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    check = colony_pcr_check(product, three_junctions, vector=puc19, flank=60, insert_primer=True)
    assert len(check.primers) == 4
    assert [clone.name for clone in check.clones] == [
        "Correct clone",
        "Empty vector",
        "Reversed insert 1",
        "Reversed insert 2",
    ]


def test_the_correct_clone_shows_one_band_reading_each_junction(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    check = colony_pcr_check(product, three_junctions, vector=puc19, flank=60, insert_primer=True)
    # 60 bases of vector then 100 into each insert, and the flanking pair across the whole span.
    assert bands(check, "Correct clone") == (160, 510, 837)
    assert bands(check, "Empty vector") == (177,)


def test_each_insert_gets_its_own_reversed_lane(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    check = colony_pcr_check(product, three_junctions, vector=puc19, flank=60, insert_primer=True)
    correct = bands(check, "Correct clone")
    assert bands(check, "Reversed insert 1") != correct
    assert bands(check, "Reversed insert 2") != correct
    assert check.tells_orientation


def test_sanger_primers_read_from_outside_the_whole_inserted_span(
    product: SequenceRecord, three_junctions: tuple[int, ...], gfp: SequenceRecord
) -> None:
    forward, reverse = sanger_primers(product, three_junctions)
    assert forward.read_bp == forward.distance_bp + len(gfp)
    assert reverse.read_bp == reverse.distance_bp + len(gfp)
