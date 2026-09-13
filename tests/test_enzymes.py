import pytest

from liulab_mbio.enzymes import Enzyme
from liulab_mbio.sequence import SequenceRecord, Strand

# A BsaI site at index 4, its spacer base, and the four bases it leaves single-stranded.
FORWARD = "AAAAGGTCTCGTTTTCCCC"
# The same stretch read on the other strand: the site now sits at index 9.
REVERSE = "GGGGAAAACGAGACCTTTT"
# The same stretch rotated so that both cuts fall past the end of the record.
ROTATED = "GTTTTCCCCAAAAGGTCTC"

BSAI = Enzyme("BsaI", "GGTCTC", top_cut=7, bottom_cut=11)


def test_a_type_iis_enzyme_cuts_outside_its_site_and_leaves_a_5_prime_overhang() -> None:
    # REBASE writes BsaI as GGTCTC(1/5): one base past the site on the top strand, five on
    # the bottom.
    bsai = Enzyme("BsaI", "GGTCTC", top_cut=7, bottom_cut=11)
    assert bsai.type == "IIS"
    assert bsai.end == "5'"
    assert bsai.overhang_length == 4


def test_sapi_leaves_a_three_base_overhang() -> None:
    sapi = Enzyme("SapI", "GCTCTTC", top_cut=8, bottom_cut=11)  # GCTCTTC(1/4)
    assert (sapi.type, sapi.end, sapi.overhang_length) == ("IIS", "5'", 3)


@pytest.mark.parametrize(
    ("name", "site", "top_cut", "bottom_cut", "end", "length"),
    [
        ("EcoRI", "GAATTC", 1, 5, "5'", 4),  # G^AATTC
        ("KpnI", "GGTACC", 5, 1, "3'", 4),  # GGTAC^C
        ("SmaI", "CCCGGG", 3, 3, "blunt", 0),  # CCC^GGG
    ],
)
def test_a_type_ii_enzyme_cuts_inside_its_site(
    name: str, site: str, top_cut: int, bottom_cut: int, end: str, length: int
) -> None:
    enzyme = Enzyme(name, site, top_cut=top_cut, bottom_cut=bottom_cut)
    assert (enzyme.type, enzyme.end, enzyme.overhang_length) == ("II", end, length)


def test_a_forward_site_cuts_downstream_of_the_site() -> None:
    assert BSAI.cut_positions(4, Strand.FORWARD) == (11, 15)


def test_a_reverse_site_cuts_upstream_of_the_site_on_the_top_strand() -> None:
    # The site occupies 9-15, so the top strand is cut five bases before it.
    assert BSAI.cut_positions(9, Strand.REVERSE) == (4, 8)


def test_the_overhang_is_the_top_strand_between_the_two_cuts() -> None:
    assert BSAI.overhang(SequenceRecord(FORWARD), 4, Strand.FORWARD) == "TTTT"
    assert BSAI.overhang(SequenceRecord(REVERSE), 9, Strand.REVERSE) == "AAAA"


def test_an_overhang_across_the_origin_of_a_circular_record_reads_through_it() -> None:
    record = SequenceRecord(ROTATED, topology="circular")
    assert BSAI.cut_positions(13, Strand.FORWARD) == (20, 24)
    assert BSAI.overhang(record, 13, Strand.FORWARD) == "TTTT"


def test_a_blunt_cutter_leaves_no_overhang() -> None:
    smai = Enzyme("SmaI", "CCCGGG", top_cut=3, bottom_cut=3)
    assert smai.overhang(SequenceRecord("AAACCCGGGTTT"), 3, Strand.FORWARD) == ""


def test_a_site_too_close_to_the_end_of_a_linear_record_cannot_be_cut() -> None:
    with pytest.raises(ValueError, match="cut"):
        BSAI.overhang(SequenceRecord(ROTATED), 13, Strand.FORWARD)
