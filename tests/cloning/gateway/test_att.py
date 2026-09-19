"""The att arithmetic, as units: the sequences, the pairings, and what a junction spells."""

from __future__ import annotations

import pytest

from liulab_mbio.cloning.gateway.att import (
    CROSSOVER,
    MISMATCHES,
    OVERLAP_BP,
    PARTNERS,
    REGION_BP,
    REGIONS,
    att_pair,
    core,
    find_att_sites,
    identify,
    joined,
    overlap,
    writes,
)
from liulab_mbio.sequence import SequenceRecord, Strand, reverse_complement

NUMBERS = (1, 2)


def test_every_shipped_site_is_one_region_long_and_shares_its_core_with_its_partner() -> None:
    assert {len(bases) for bases in REGIONS.values()} == {REGION_BP}
    for name, partner in PARTNERS.items():
        assert core(REGIONS[name]) == core(REGIONS[partner])


@pytest.mark.parametrize("number", NUMBERS)
def test_bp_writes_attl_and_attr_and_lr_writes_them_back_to_attb(number: int) -> None:
    attb, attp = REGIONS[f"attB{number}"], REGIONS[f"attP{number}"]

    attl, attr = joined(attp, attb), joined(attb, attp)

    assert (attl, attr) == (REGIONS[f"attL{number}"], REGIONS[f"attR{number}"])
    # The round trip is on attB: the by-product's attP is shorter than the donor vector's,
    # because attR was built with 33 bp of the arm deleted.
    assert joined(attr, attl) == attb


@pytest.mark.parametrize("number", NUMBERS)
def test_writes_names_the_site_each_junction_spells(number: int) -> None:
    assert writes(f"attR{number}", f"attL{number}") == f"attB{number}"
    assert writes(f"attL{number}", f"attR{number}") == f"attP{number}"
    assert writes(f"attP{number}", f"attB{number}") == f"attL{number}"
    assert writes(f"attB{number}", f"attP{number}") == f"attR{number}"


def test_writes_refuses_two_sites_that_do_not_react() -> None:
    with pytest.raises(KeyError, match="attB1 reacts with attP1"):
        writes("attB1", "attP2")


def test_identify_names_every_shipped_site() -> None:
    assert {name: identify(bases) for name, bases in REGIONS.items()} == {
        name: name for name in REGIONS
    }


def test_a_drifted_flank_still_reads_as_its_own_site_and_never_as_another() -> None:
    outside = [at for at in range(REGION_BP) if not CROSSOVER <= at < CROSSOVER + OVERLAP_BP]

    for name, bases in REGIONS.items():
        for at in outside:
            for base in "ACGT":
                drifted = identify(bases[:at] + base + bases[at + 1 :])
                assert drifted in (name, None)


def test_a_change_inside_the_overlap_is_never_read_as_that_site() -> None:
    for name, bases in REGIONS.items():
        for at in range(CROSSOVER, CROSSOVER + OVERLAP_BP):
            for base in "ACGT":
                if base != bases[at]:
                    assert identify(bases[:at] + base + bases[at + 1 :]) != name


def test_identify_refuses_anything_that_is_not_a_region() -> None:
    with pytest.raises(ValueError, match="25 bases, got 10"):
        identify("ACGTACGTAC")


def test_a_site_is_found_on_either_strand_and_across_the_origin() -> None:
    filler = "TTGACCGCTTAAGGCTAACGTCAGT"
    forward = SequenceRecord(filler + REGIONS["attL1"] + filler)
    reverse = SequenceRecord(filler + reverse_complement(REGIONS["attL1"]) + filler)
    wrapped = SequenceRecord(
        REGIONS["attR2"][10:] + filler + REGIONS["attR2"][:10], topology="circular"
    )

    assert [(one.name, one.start, one.strand) for one in find_att_sites(forward)] == [
        ("attL1", 25, Strand.FORWARD)
    ]
    assert [(one.name, one.strand) for one in find_att_sites(reverse)] == [
        ("attL1", Strand.REVERSE)
    ]
    found = find_att_sites(wrapped)
    assert [(one.name, one.start, one.end) for one in found] == [("attR2", 40, 65)]
    assert wrapped.extract(found[0].span) == REGIONS["attR2"]


def test_att_pair_refuses_a_record_carrying_none_and_says_what_it_looked_for() -> None:
    with pytest.raises(ValueError, match="looked for attL1 and attL2 on either strand"):
        att_pair(SequenceRecord("ACGT" * 20, name="pUC19"), "attL")


def test_att_pair_refuses_a_record_carrying_the_same_site_twice() -> None:
    filler = "TTGACCGCTTAAGGCTAACGTCAGT"
    twice = SequenceRecord(
        REGIONS["attL1"]
        + filler
        + reverse_complement(REGIONS["attL2"])
        + filler
        + REGIONS["attL1"],
        topology="circular",
        name="pENTR-GFP",
    )

    with pytest.raises(ValueError, match="carries 2 attL1 sites, needing one"):
        att_pair(twice, "attL")


def test_att_pair_refuses_two_sites_that_do_not_point_at_each_other() -> None:
    filler = "TTGACCGCTTAAGGCTAACGTCAGT"
    both = SequenceRecord(REGIONS["attL1"] + filler + REGIONS["attL2"], name="pBAD")

    with pytest.raises(ValueError, match="attL1 and attL2 on the same strand"):
        att_pair(both, "attL")


def test_att_pair_refuses_a_kind_this_module_holds_no_pair_for() -> None:
    with pytest.raises(KeyError, match="no att site called 'attX1'"):
        att_pair(SequenceRecord("ACGT" * 20), "attX")


def test_the_shipped_tolerance_is_as_far_as_a_window_may_drift_and_stay_unambiguous() -> None:
    apart = min(
        sum(one != other for one, other in zip(REGIONS[first], REGIONS[second], strict=True))
        for first in REGIONS
        for second in REGIONS
        if first != second and overlap(REGIONS[first]) == overlap(REGIONS[second])
    )

    assert apart > MISMATCHES
