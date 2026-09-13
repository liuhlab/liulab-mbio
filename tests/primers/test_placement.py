import pytest

from liulab_mbio.primers import Placement, find_binding_sites, find_priming_sites
from liulab_mbio.sequence import BindingSite, Segment, SequenceRecord, Strand

from .sequences import M13_FWD, M13_REV, MCS_FWD


def test_a_primer_is_placed_on_a_template_by_its_three_prime_match(puc19) -> None:
    assert find_binding_sites(M13_FWD, puc19) == (BindingSite(378, 395, Strand.FORWARD),)
    assert find_binding_sites(M13_REV, puc19) == (BindingSite(464, 481, Strand.REVERSE),)
    tailed = "TTGAAGACAA" + MCS_FWD
    assert find_binding_sites(tailed, puc19) == (BindingSite(452, 472, Strand.FORWARD),)


def test_a_binding_site_may_run_across_the_origin(puc19) -> None:
    across = puc19.sequence[-10:] + puc19.sequence[:10]
    assert find_binding_sites(across, puc19) == (BindingSite(2676, 2696, Strand.FORWARD),)


def test_a_placement_bounds_at_least_one_end() -> None:
    with pytest.raises(ValueError, match="at least one end"):
        Placement()


def test_a_placement_reads_each_end_of_a_site_by_its_strand(puc19) -> None:
    placement = Placement(five_prime=Segment(452, 463))
    assert placement.allows(BindingSite(452, 472, Strand.FORWARD), puc19)
    assert not placement.allows(BindingSite(451, 472, Strand.FORWARD), puc19)
    # A reverse primer's 5' end is the end of its binding site, not its start.
    assert placement.allows(BindingSite(432, 452, Strand.REVERSE), puc19)
    assert not placement.allows(BindingSite(452, 472, Strand.REVERSE), puc19)


def test_a_placement_across_the_origin_holds_both_sides_of_it(puc19) -> None:
    across = Placement(five_prime=Segment(len(puc19) - 5, len(puc19) + 5))
    assert across.allows(BindingSite(len(puc19) - 3, len(puc19) + 17, Strand.FORWARD), puc19)
    assert across.allows(BindingSite(2, 22, Strand.FORWARD), puc19)
    assert not across.allows(BindingSite(10, 30, Strand.FORWARD), puc19)


def test_a_placement_on_a_linear_template_does_not_count_round_it() -> None:
    template = SequenceRecord("ACGT" * 10)
    placement = Placement(five_prime=Segment(0, 3))
    assert placement.allows(BindingSite(1, 21, Strand.FORWARD), template)
    # The far end of a linear template is not its origin, so a 5' end there lies outside.
    assert not placement.allows(BindingSite(20, 40, Strand.REVERSE), template)


def test_priming_sites_carry_their_strand_mismatches_and_tm(puc19) -> None:
    sites = find_priming_sites(M13_FWD, puc19)
    assert [site.site for site in sites] == [BindingSite(378, 395, Strand.FORWARD)]
    assert sites[0].mismatches == 0
    assert sites[0].tm == pytest.approx(54.61, abs=0.01)
