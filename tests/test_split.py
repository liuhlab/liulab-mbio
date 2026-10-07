"""Splitting a cargo into oligo-sized fragments that assemble back into it."""

import pytest

from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.split import Budget, fewest_pieces, split_cargo

#: A cargo with no repeat, so every window offers real choices rather than one overhang.
CARGO = SequenceRecord(
    "".join(
        "ACGTTGCAAGCTTCGAGGATCCTAGCATGCTAGGCCATTAACGGTTACGATCGATTACGCAT"[i % 62] for i in range(900)
    )
)


def test_the_fragments_assemble_back_into_the_cargo():
    split = split_cargo(CARGO, "BsaI", budget=Budget(350, 74))
    assert split.reassembled() == str(CARGO.sequence)
    assert split.pieces == fewest_pieces(len(CARGO.sequence), Budget(350, 74))


def test_every_fragment_fits_the_oligo_the_last_one_included():
    """The bound is checked on the final fragment, not only on each new cut.

    A length bound checked only on the way forward lets the last fragment overshoot, which is
    what emits an oligo the vendor will not make.
    """
    budget = Budget(350, 74)
    split = split_cargo(CARGO, "BsaI", budget=budget)
    assert [one.length for one in split.fragments][-1] <= budget.span
    assert all(one.length <= budget.span for one in split.fragments)


def test_a_fragment_carries_an_overhang_at_each_end_and_the_budget_counts_both():
    split = split_cargo(CARGO, "BsaI", budget=Budget(350, 74))
    spans = sum(one.length for one in split.fragments)
    assert spans == len(CARGO.sequence) + 4 * (split.pieces - 1)
    for one in split.fragments:
        assert len(one.overhangs[0]) == 4
        assert len(one.overhangs[1]) == 4


def test_each_junction_carries_the_window_the_program_proved_and_not_a_constant():
    wide = split_cargo(CARGO, "BsaI", budget=Budget(350, 74))
    narrow = split_cargo(CARGO, "BsaI", budget=Budget(250, 74))
    assert len({one.window for one in wide.junctions}) > 1
    assert max(one.window for one in wide.junctions) > max(one.window for one in narrow.junctions)


def test_an_overhang_held_out_by_name_is_never_taken():
    split = split_cargo(CARGO, "BsaI", budget=Budget(200, 74))
    held = split.overhangs[0]
    again = split_cargo(CARGO, "BsaI", budget=Budget(200, 74), reserved=[held])
    assert held not in again.overhangs


def test_a_cargo_fitting_one_oligo_is_one_fragment_with_no_junction():
    short = SequenceRecord(str(CARGO.sequence)[:200])
    split = split_cargo(short, "BsaI", budget=Budget(350, 74))
    assert split.pieces == 1
    assert split.junctions == ()
    assert split.reassembled() == str(short.sequence)


def test_it_refuses_a_cargo_past_the_ceiling_and_says_where_the_ceiling_is():
    with pytest.raises(ValueError, match="supply"):
        split_cargo(SequenceRecord("AT" * 4000), "BsaI", budget=Budget(120, 74))


def test_a_circular_record_has_no_first_fragment():
    with pytest.raises(ValueError, match="circular"):
        split_cargo(
            SequenceRecord("ACGT" * 200, topology="circular"), "BsaI", budget=Budget(350, 74)
        )
