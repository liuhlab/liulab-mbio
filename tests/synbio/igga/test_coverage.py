"""Library coverage: the construct count, and the colonies one round needs for it."""

import pytest

from liulab_synbio.igga.coverage import (
    REPRESENTATION_MARKS,
    RoundCoverage,
    absent_probability,
    colonies_for_completeness,
    colonies_for_coverage,
    constructs,
    plan_coverage,
    reads_for_representation,
    skew_ratio,
)


def test_the_construct_count_is_the_product_of_the_part_list_sizes() -> None:
    assert constructs((24, 24, 24)) == 13824


def test_a_library_needs_a_part_list() -> None:
    with pytest.raises(ValueError, match="at least one part list"):
        constructs(())


def test_a_part_list_cannot_be_empty() -> None:
    with pytest.raises(ValueError, match="at least one part"):
        constructs((24, 0))


def test_colonies_are_the_coverage_asked_for_times_the_products() -> None:
    assert colonies_for_coverage(576, 10) == 5760


def test_a_fractional_colony_rounds_up() -> None:
    assert colonies_for_coverage(7, 1.5) == 11


def test_missing_a_product_gets_less_likely_with_more_colonies() -> None:
    assert absent_probability(100, 500) == pytest.approx(0.0066, abs=0.0001)
    assert absent_probability(100, 1000) < absent_probability(100, 500)


def test_no_colonies_leaves_every_product_absent() -> None:
    assert absent_probability(100, 0) == 1.0


def test_completeness_follows_clarke_and_carbon() -> None:
    # N = ln(1 - P) / ln(1 - f), f = 1/100, asked of all 100 members at once.
    assert colonies_for_completeness(100, 0.95) == 754


def test_a_lone_product_needs_one_colony() -> None:
    assert colonies_for_completeness(1, 0.95) == 1


def test_completeness_refuses_a_certainty() -> None:
    with pytest.raises(ValueError, match="between 0 and 1"):
        colonies_for_completeness(100, 1.0)


def test_every_round_is_counted_on_its_own() -> None:
    rows = plan_coverage((4, 6), completeness=0.99)

    assert [row.number for row in rows] == [1, 2]
    assert [row.part_list_size for row in rows] == [4, 6]
    assert [row.products for row in rows] == [4, 24]
    assert [row.colonies for row in rows] == [21, 183]


def test_the_multiple_a_round_works_out_at_is_derived() -> None:
    """The project states the completeness; the ratio is what it costs, and it moves."""
    rows = plan_coverage((4, 6), completeness=0.99)

    assert [round(row.coverage, 2) for row in rows] == [5.25, 7.62]


def test_a_round_reports_what_its_colonies_would_miss() -> None:
    (row,) = plan_coverage((100,), completeness=0.99)

    assert isinstance(row, RoundCoverage)
    assert row.completeness == 0.99
    assert row.absent_probability == pytest.approx(0.0001, abs=0.00001)


def test_completeness_has_to_be_asked_for() -> None:
    with pytest.raises(ValueError, match="between 0 and 1"):
        plan_coverage((4,), completeness=0)


def test_the_skew_ratio_is_the_ninetieth_over_the_tenth_percentile():
    """Joung's measure: an even library answers 1, and a spread one answers more."""
    assert skew_ratio([7, 7, 7, 7]) == 1.0
    assert skew_ratio([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) == 5.0


def test_a_library_a_tenth_of_which_went_unread_carries_no_skew_ratio():
    """The ratio would divide by zero, and the share seen is the number that says it."""
    with pytest.raises(ValueError, match="10th percentile"):
        skew_ratio([0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9])


def test_the_read_depth_follows_from_the_library_width():
    """Joung judges at 100 reads a member, so the depth is computed and never stated."""
    assert reads_for_representation(13824) == 1_382_400
    assert reads_for_representation(10, per_member=150) == 1500


def test_the_sourced_marks_are_joungs():
    """What the representation read is held to, until a project tightens it."""
    assert REPRESENTATION_MARKS.seen == 0.995
    assert REPRESENTATION_MARKS.skew == 10.0
    assert REPRESENTATION_MARKS.reads_per_member == 100
