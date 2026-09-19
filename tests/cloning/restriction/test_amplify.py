"""The PCR that puts a site on each end of an insert no plasmid holds.

The insert is the GFP record the rest of the suite designs on, amplified into the backbone the
vector's own digest leaves. `test_plan.py` is where the protocol the PCR joins is checked.
"""

import pytest

from liulab_mbio.cloning.restriction.amplify import Amplicon, amplified
from liulab_mbio.cloning.restriction.digest import Piece, opened, resolve
from liulab_mbio.sequence import SequenceRecord, reverse_complement
from liulab_mbio.sites import SPACER_LENGTH, find_sites


@pytest.fixture(scope="module")
def backbone(puc19: SequenceRecord) -> Piece:
    """pUC19 opened with EcoRI and BamHI, which is what names the enzyme for each end."""
    return opened(puc19, resolve(["EcoRI", "BamHI"]))[0]


@pytest.fixture(scope="module")
def amplicon(gfp: SequenceRecord, backbone: Piece) -> Amplicon:
    """GFP amplified with an EcoRI tail at one end and a BamHI tail at the other."""
    return amplified(gfp, into=backbone)


def test_each_tail_is_a_spacer_and_the_site_of_the_end_it_has_to_anneal_to(amplicon, gfp):
    forward, reverse = amplicon.primers
    assert (forward.name, reverse.name) == ("GFP forward", "GFP reverse")
    assert forward.sequence.startswith(amplicon.left_tail)
    assert reverse.sequence.startswith(reverse_complement(amplicon.right_tail))
    assert amplicon.record.sequence == amplicon.left_tail + gfp.sequence + amplicon.right_tail
    # Each tail is the spacer the note asks for and then the site, and the plan says which
    # enzyme each end is for.
    for _, tail, enzyme in amplicon.ends:
        assert enzyme.site in tail
        assert len(tail) - len(enzyme.site) == SPACER_LENGTH
    assert (amplicon.left_enzyme.name, amplicon.right_enzyme.name) == ("EcoRI", "BamHI")
    # The two tails spell one site of each enzyme and no more, so the amplicon cuts cleanly.
    found = [site for site in find_sites(amplicon.record, ["EcoRI", "BamHI"]) if site.cuts]
    assert sorted(site.enzyme.name for site in found) == ["BamHI", "EcoRI"]


def test_one_enzyme_puts_its_own_site_on_both_tails(puc19, gfp):
    # A tail spells its own enzyme's site by design, so that enzyme is not its own avoid.
    one = amplified(gfp, into=opened(puc19, resolve(["SmaI"]))[0])
    assert one.left_enzyme.name == one.right_enzyme.name == "SmaI"
    assert len([site for site in find_sites(one.record, "SmaI") if site.cuts]) == 2
