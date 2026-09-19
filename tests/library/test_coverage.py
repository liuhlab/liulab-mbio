"""Library coverage: the construct count, and the colonies one round needs for it."""

import pytest

from liulab_mbio.library.coverage import (
    RoundCoverage,
    absent_probability,
    colonies_for_completeness,
    colonies_for_coverage,
    constructs,
    plan_coverage,
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
    rows = plan_coverage((4, 6), coverage=10)

    assert [row.number for row in rows] == [1, 2]
    assert [row.part_list_size for row in rows] == [4, 6]
    assert [row.products for row in rows] == [4, 24]
    assert [row.colonies for row in rows] == [40, 240]


def test_a_round_reports_what_its_colonies_would_miss() -> None:
    (row,) = plan_coverage((100,), coverage=1)

    assert isinstance(row, RoundCoverage)
    assert row.absent_probability == pytest.approx(0.366, abs=0.001)


def test_coverage_has_to_be_asked_for() -> None:
    with pytest.raises(ValueError, match="coverage must be positive"):
        plan_coverage((4,), coverage=0)
