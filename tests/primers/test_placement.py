import pytest

from liulab_mbio.primers import find_binding_sites, find_priming_sites
from liulab_mbio.sequence import BindingSite, Strand

from .sequences import M13_FWD, M13_REV, MCS_FWD


def test_a_primer_is_placed_on_a_template_by_its_three_prime_match(puc19) -> None:
    assert find_binding_sites(M13_FWD, puc19) == (BindingSite(378, 395, Strand.FORWARD),)
    assert find_binding_sites(M13_REV, puc19) == (BindingSite(464, 481, Strand.REVERSE),)
    tailed = "TTGAAGACAA" + MCS_FWD
    assert find_binding_sites(tailed, puc19) == (BindingSite(452, 472, Strand.FORWARD),)


def test_a_binding_site_may_run_across_the_origin(puc19) -> None:
    across = puc19.sequence[-10:] + puc19.sequence[:10]
    assert find_binding_sites(across, puc19) == (BindingSite(2676, 2696, Strand.FORWARD),)


def test_priming_sites_carry_their_strand_mismatches_and_tm(puc19) -> None:
    sites = find_priming_sites(M13_FWD, puc19)
    assert [site.site for site in sites] == [BindingSite(378, 395, Strand.FORWARD)]
    assert sites[0].mismatches == 0
    assert sites[0].tm == pytest.approx(54.61, abs=0.01)
