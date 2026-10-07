"""Seating: one part a well, and no library at the end of it."""

import pytest

from liulab_synbio import seating


def test_seating_ends_in_one_plasmid_a_well_and_no_pool():
    """A round hands back a library; this hands back parts that are still told apart."""
    seated = seating.seat_parts(["FLAG", "GS linker", "NLS", "degron"])
    assert seated.products == 4
    assert sorted(seated.plate.seating) == ["A1", "A2", "A3", "A4"]
    assert not hasattr(seated, "library")


def test_the_step_says_what_comes_out_and_never_names_a_library():
    """The reaction the round model had no room for is expressible, and reads as itself."""
    step = seating.seating_step(seating.seat_parts(["FLAG", "NLS"]))
    said = " ".join((*step.instructions, *step.expected, *step.notes)).lower()
    assert seating.ENZYME.lower() in said
    assert "no pool and no library" in said
    assert "2 carrier plasmids" in " ".join(step.expected)


def test_the_plate_is_the_smallest_format_that_holds_the_parts():
    """Format is one parameter, so seating picks rather than hard-coding a plate."""
    assert seating.seat_parts([f"p{n}" for n in range(11)]).plate.wells == 12
    assert seating.seat_parts([str(n) for n in range(13)]).plate.wells == 24


def test_two_parts_cannot_share_a_name():
    """A well is named by the part sitting in it, so a repeat would lose one of them."""
    with pytest.raises(ValueError, match="share a name"):
        seating.seat_parts(["NLS", "NLS"])
    with pytest.raises(ValueError, match="at least one"):
        seating.seat_parts([])
