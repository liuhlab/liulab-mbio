"""A clone's sequencing results held against its record: what each region and each result comes to.

Every record here is a short synthetic product, a vector holding one insert between two
junctions, so each alignment runs in milliseconds.
"""

import random

import pytest

from mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement
from mbio.verification.judge import Verification, verify
from mbio.verification.result import SequencingResult

_RNG = random.Random(565)
LEFT, INSERT, RIGHT = ("".join(_RNG.choice("ACGT") for _ in range(n)) for n in (200, 700, 200))
LENGTH = len(LEFT) + len(INSERT) + len(RIGHT)
ORI = Feature("ori", "rep_origin", (Segment(40, 120),))
JUNCTION_1 = Feature("junction 1", "misc_feature", (Segment(195, 205),))
THE_INSERT = Feature("insert", "misc_feature", (Segment(205, 895),))
JUNCTION_2 = Feature("junction 2", "misc_feature", (Segment(895, 905),))
REGIONS = (JUNCTION_1, THE_INSERT, JUNCTION_2)
PRODUCT = SequenceRecord(LEFT + INSERT + RIGHT, topology="circular", features=(ORI,))


def bases(start: int, end: int, changes: dict[int, str] | None = None) -> str:
    """Return the product's bases over a span, with a base at each named position replaced."""
    read = list(PRODUCT.bases(start, end))
    for position, base in (changes or {}).items():
        read[position - start] = base
    return "".join(read)


def other(position: int) -> str:
    """Return a base unlike the product's at `position` and both its neighbours."""
    near = set(PRODUCT.bases(position - 1, position + 2))
    return next(one for one in "ACGT" if one not in near)


def statuses(verification: Verification) -> dict[str, str | None]:
    """Return each region's verdict by name."""
    return {one.name: one.status for one in verification.checks}


def test_planted_changes_are_placed_on_the_insert_and_fail_it():
    """Three point changes in an 800-base read are listed where they lie, and nothing gates it."""
    read = bases(150, 500, {300: other(300)}) + other(500) + bases(500, 600) + bases(603, 953)
    found = verify(PRODUCT, (SequencingResult("r", read),), REGIONS)
    assert [(one.kind, one.start, one.end) for one in found.disagreements] == [
        ("substitution", 300, 301),
        ("insertion", 500, 500),
        ("deletion", 600, 603),
    ]
    assert all(one.regions == ("insert",) for one in found.disagreements)
    assert statuses(found) == {"junction 1": "pass", "insert": "fail", "junction 2": "pass"}
    assert found.result_checks[0].status == "pass"
    assert found.status == "fail"
    assert found.verified is False


def test_a_read_across_the_origin_lands_in_one_piece():
    across = Feature("backbone", "misc_feature", (Segment(1080, 1120),))
    found = verify(PRODUCT, (SequencingResult("r", bases(1000, 1200)),), (across,))
    assert found.placements[0].span == Segment(1000, 1200)
    assert statuses(found) == {"backbone": "pass"}
    assert found.verified is True


def test_a_read_on_the_reverse_strand_reports_the_top_strand():
    changed = other(400)
    read = reverse_complement(bases(150, 950, {400: changed}))
    found = verify(PRODUCT, (SequencingResult("r", read),), REGIONS)
    assert found.placements[0].strand is Strand.REVERSE
    assert [(one.kind, one.start, one.bases) for one in found.disagreements] == [
        ("substitution", 400, changed)
    ]


def test_a_junction_left_unread_carries_no_verdict_and_the_clone_is_not_verified():
    found = verify(PRODUCT, (SequencingResult("r", bases(150, 880)),), REGIONS)
    assert statuses(found) == {"junction 1": "pass", "insert": None, "junction 2": None}
    assert "unread: 881 .. 895" in found.checks[1].detail
    assert found.verified is False


def test_two_results_contradicting_each_other_at_a_base_warn():
    one = SequencingResult("forward", bases(150, 950))
    two = SequencingResult("reverse", reverse_complement(bases(150, 950, {400: other(400)})))
    found = verify(PRODUCT, (one, two), REGIONS)
    assert statuses(found)["insert"] == "warn"
    assert found.disagreements[0].results == ("reverse",)


def test_a_result_carrying_a_second_consensus_is_mixed_and_fails():
    mixed = SequencingResult("r", bases(150, 950), others=(bases(150, 950, {400: other(400)}),))
    found = verify(PRODUCT, (mixed,), REGIONS)
    assert found.result_checks[0].status == "fail"
    assert set(statuses(found).values()) == {"fail"}


def test_a_disagreement_outside_every_region_is_named_by_its_feature_and_judges_nothing():
    found = verify(PRODUCT, (SequencingResult("r", bases(0, 200, {80: other(80)})),), REGIONS)
    (outside,) = found.disagreements
    assert (outside.regions, outside.features) == ((), ("ori",))
    assert "fail" not in statuses(found).values()


def test_an_unrelated_read_is_gated_and_shows_nothing():
    unrelated = "".join(random.Random(1).choice("ACGT") for _ in range(600))
    found = verify(PRODUCT, (SequencingResult("r", unrelated),), REGIONS)
    check = found.result_checks[0]
    assert check.status == "fail"
    assert check.detail.startswith("does not read as this product: ")
    assert found.disagreements == ()
    assert found.placements[0].span is None


def test_a_different_version_of_the_insert_is_gated():
    """About three bases in four agree, as two versions of one domain might."""
    changes = {position: other(position) for position in range(205, 895, 4)}
    found = verify(PRODUCT, (SequencingResult("r", bases(150, 950, changes)),), REGIONS)
    assert found.result_checks[0].status == "fail"
    assert found.disagreements == ()


def test_a_gated_result_leaves_the_regions_to_the_good_one():
    unrelated = "".join(random.Random(2).choice("ACGT") for _ in range(600))
    good = SequencingResult("good", bases(150, 950))
    found = verify(PRODUCT, (SequencingResult("bad", unrelated), good), REGIONS)
    assert set(statuses(found).values()) == {"pass"}
    assert [one.status for one in found.result_checks] == ["fail", "pass"]
    assert found.verified is False


def test_an_empty_vector_read_fails_the_missing_insert():
    """A short insert lost whole reads as one deletion, wherever its equal ends let it slide."""
    insert = Feature("insert", "misc_feature", (Segment(200, 224),))
    product = SequenceRecord(LEFT + INSERT[:24] + RIGHT, topology="circular")
    found = verify(product, (SequencingResult("r", LEFT[60:] + RIGHT[:140]),), (insert,))
    assert [(one.kind, one.end - one.start) for one in found.disagreements] == [("deletion", 24)]
    assert statuses(found) == {"insert": "fail"}


@pytest.mark.parametrize(
    ("result", "detail"),
    [
        (SequencingResult("r", "ACGT", quality=(0, 0, 0, 0), trusted=(0, 0)), "no quality values"),
        (SequencingResult("r", "ACGT", trusted=(2, 2)), "the trimmed span is empty"),
        (SequencingResult("r", "ACGT", quality=(40, 12, 15, 40), trusted=(1, 3)), "quality 20"),
    ],
)
def test_a_result_trusting_no_base_carries_no_verdict(result, detail):
    check = verify(PRODUCT, (result,), REGIONS).result_checks[0]
    assert check.status is None
    assert detail in check.detail


def test_a_base_counts_only_where_its_quality_or_depth_trusts_it():
    """Below Q20 a base reads nothing; a depth from 11 to 19 warns; 10 or fewer reads nothing."""
    read = bases(150, 950, {400: other(400)})
    quality = tuple(10 if index == 250 else 40 for index in range(800))
    low = verify(PRODUCT, (SequencingResult("r", read, quality=quality),), REGIONS)
    assert low.disagreements == ()
    assert "unread: 401" in low.checks[1].detail
    depth = tuple(15 if index < 100 else 30 for index in range(800))
    thin = verify(PRODUCT, (SequencingResult("r", bases(150, 950), depth=depth),), REGIONS)
    assert statuses(thin) == {"junction 1": "warn", "insert": "warn", "junction 2": "pass"}
    shallow = tuple(5 if index == 250 else 30 for index in range(800))
    unread = verify(PRODUCT, (SequencingResult("r", bases(150, 950), depth=shallow),), REGIONS)
    assert statuses(unread)["insert"] is None


def test_a_result_refuses_what_does_not_fit_its_bases():
    with pytest.raises(ValueError, match="quality"):
        SequencingResult("r", "ACGT", quality=(40, 40))
    with pytest.raises(ValueError, match="depth"):
        SequencingResult("r", "ACGT", depth=(30,) * 5)
    with pytest.raises(ValueError, match="outside"):
        SequencingResult("r", "ACGT", trusted=(2, 6))


def test_bases_past_a_linear_record_are_not_listed_and_count_on_neither_side():
    """The record says nothing of what lies beyond it, so more flank than region still passes."""
    region = Feature("designed region", "misc_feature", (Segment(0, 60),))
    consensus = LEFT[:100] + INSERT[:60] + RIGHT[:100]
    found = verify(SequenceRecord(INSERT[:60]), (SequencingResult("r", consensus),), (region,))
    assert found.disagreements == ()
    assert found.verified is True
