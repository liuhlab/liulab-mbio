from pathlib import Path

import pytest

from liulab_mbio.edits import EditReport
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.io import read_record
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.sites import (
    digest,
    find_sites,
    free_enzymes,
    has_site,
    insert_site,
    primer_tail,
    site_counts,
)

DATA = Path(__file__).parent / "data"

# A BsaI site at index 4, its spacer base, and the four bases it leaves single-stranded.
FORWARD = "AAAAGGTCTCGTTTTCCCC"
# The same stretch read on the other strand: the site now sits at index 9.
REVERSE = "GGGGAAAACGAGACCTTTT"
# The same stretch rotated so that both cuts fall past the end of the record.
ROTATED = "GTTTTCCCCAAAAGGTCTC"
# Twenty bases whose BsaI site begins at 17 and finishes in the first three.
ACROSS = "CTC" + "A" * 14 + "GGT"

BSAI = get_enzyme("BsaI")
# An invented enzyme whose site carries an IUPAC code and is not its own reverse complement.
FOOI = Enzyme("FooI", "GGWCA", top_cut=1, bottom_cut=5)


@pytest.fixture(scope="module")
def puc19() -> SequenceRecord:
    return read_record(DATA / "pUC19.dna")


@pytest.fixture(scope="module")
def gfp() -> SequenceRecord:
    return read_record(DATA / "GFP.dna")


def test_a_forward_site_reports_its_span_its_cuts_and_its_overhang() -> None:
    (site,) = find_sites(SequenceRecord(FORWARD), "BsaI")
    assert site.enzyme is BSAI
    assert (site.start, site.end, site.strand) == (4, 10, Strand.FORWARD)
    assert site.span == Segment(4, 10)
    assert (site.top_cut, site.bottom_cut) == (11, 15)
    assert site.overhang == "TTTT"
    assert site.cuts


def test_a_reverse_site_reports_the_top_strand_span_and_cuts_upstream_of_it() -> None:
    (site,) = find_sites(SequenceRecord(REVERSE), BSAI)
    assert (site.start, site.end, site.strand) == (9, 15, Strand.REVERSE)
    assert (site.top_cut, site.bottom_cut) == (4, 8)
    assert site.overhang == "AAAA"


def test_a_site_across_the_origin_is_found_and_its_cuts_reduce_into_the_record() -> None:
    (site,) = find_sites(SequenceRecord(ACROSS, topology="circular"), "BsaI")
    assert site.span == Segment(17, 23)
    assert (site.top_cut, site.bottom_cut) == (4, 8)
    assert site.overhang == "AAAA"


def test_a_linear_record_ending_inside_the_cut_reports_the_site_but_no_overhang() -> None:
    (site,) = find_sites(SequenceRecord(ROTATED), "BsaI")
    assert site.start == 13
    assert site.overhang is None
    assert not site.cuts


def test_a_palindromic_site_is_counted_once_and_reported_on_the_forward_strand() -> None:
    sites = find_sites(SequenceRecord("AAAGAATTCAAA"), "EcoRI")
    assert [(s.start, s.strand) for s in sites] == [(3, Strand.FORWARD)]


def test_an_iupac_code_in_the_site_matches_every_base_it_admits() -> None:
    (site,) = find_sites(SequenceRecord("AAGGACAAA"), FOOI)
    assert site.start == 2
    assert site.certain


def test_an_iupac_code_in_the_template_is_reported_as_a_site_that_may_be_there() -> None:
    (site,) = find_sites(SequenceRecord("AAAAGGNCTCGTTTT"), "BsaI")
    assert (site.start, site.overhang) == (4, "TTTT")
    assert not site.certain


def test_a_template_base_the_site_cannot_admit_is_not_a_site() -> None:
    assert find_sites(SequenceRecord("AAGGGCAAA"), FOOI) == ()


def test_several_enzymes_are_searched_at_once_and_the_hits_come_back_in_order() -> None:
    record = SequenceRecord("AAAGAATTCAAAGGTCTCAAA")
    assert [(s.enzyme.name, s.start) for s in find_sites(record, ["EcoRI", "BsaI"])] == [
        ("EcoRI", 3),
        ("BsaI", 12),
    ]


def test_the_fixture_tables_hold_the_sites_issue_1_lists(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    # 1-based starts, as issue #1 writes them, with the strand each site lies on.
    assert _found(puc19, "BsaI") == [(1766, Strand.REVERSE)]
    assert _found(puc19, "BsmBI") == [(51, Strand.REVERSE), (2683, Strand.FORWARD)]
    assert _found(puc19, "BbsI") == []
    assert _found(puc19, "SapI") == [(683, Strand.FORWARD)]
    assert _found(puc19, "PaqCI") == []
    assert _found(gfp, "BsaI") == [(644, Strand.REVERSE)]
    for name in ("BsmBI", "BbsI", "SapI", "PaqCI"):
        assert _found(gfp, name) == []


def test_the_bsmbi_site_at_2683_runs_across_the_origin_of_puc19(puc19: SequenceRecord) -> None:
    across = next(s for s in find_sites(puc19, "BsmBI") if s.start == 2682)
    assert across.span == Segment(2682, 2688)
    assert across.end > len(puc19)
    assert puc19.extract(across.span) == "CGTCTC"


def test_has_site_answers_whether_one_enzyme_still_cuts(puc19: SequenceRecord) -> None:
    assert has_site(puc19, "BsaI")
    assert not has_site(puc19, "BbsI")


def test_site_counts_and_free_enzymes_rank_the_type_iis_set(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    names = ["BsaI", "BsmBI", "BbsI", "SapI", "PaqCI"]
    assert site_counts([puc19, gfp], names) == {
        "BsaI": 2,
        "BsmBI": 2,
        "BbsI": 0,
        "SapI": 1,
        "PaqCI": 0,
    }
    assert [e.name for e in free_enzymes([puc19, gfp], names)] == ["BbsI", "PaqCI"]


def test_an_unknown_enzyme_name_is_refused() -> None:
    with pytest.raises(KeyError, match="NoSuchI"):
        find_sites(SequenceRecord("AAAA"), "NoSuchI")


def test_one_cut_splits_a_linear_record_into_two_fragments() -> None:
    record = SequenceRecord(FORWARD)
    left, right = digest(record, "BsaI")
    # The top strand is severed at 11, so the two top strands are 11 and 8 bases of the 19.
    assert (left.start, left.end, left.length) == (0, 11, 11)
    assert (right.start, right.end, right.length) == (11, 19, 8)
    # The four single-stranded bases sit on the top strand of one fragment and the bottom
    # strand of the other, so both ends spell the same thing and they anneal.
    assert (left.left_overhang, left.right_overhang) == ("", "TTTT")
    assert (right.left_overhang, right.right_overhang) == ("TTTT", "")


def test_a_circular_record_gives_one_fragment_per_cut(puc19: SequenceRecord) -> None:
    # BsmBI cuts pUC19's top strand at 3 and at 45.
    first, second = digest(puc19, "BsmBI")
    assert [fragment.length for fragment in (first, second)] == [42, 2644]
    assert sum(fragment.length for fragment in (first, second)) == len(puc19)
    # The second fragment runs across the origin, so it ends past the record's length.
    assert (second.start, second.end) == (45, 2689)
    assert first.right_overhang == second.left_overhang
    assert second.right_overhang == first.left_overhang
    assert len(first.left_overhang) == 4


def test_an_uncut_linear_record_is_one_blunt_fragment() -> None:
    (whole,) = digest(SequenceRecord("AAAACCCC"), "BsaI")
    assert (whole.start, whole.end, whole.length) == (0, 8, 8)
    assert (whole.left_overhang, whole.right_overhang) == ("", "")


def test_an_uncut_circular_record_is_not_digested_at_all() -> None:
    assert digest(SequenceRecord("AAAACCCC", topology="circular"), "BsaI") == ()


def test_a_site_whose_cut_falls_off_a_linear_end_does_not_split_it() -> None:
    (whole,) = digest(SequenceRecord(ROTATED), "BsaI")
    assert whole.length == len(ROTATED)


def test_a_blunt_cutter_leaves_fragments_with_no_overhang() -> None:
    left, right = digest(SequenceRecord("AAACCCGGGTTT"), "SmaI")
    assert (left.length, right.length) == (6, 6)
    assert {left.right_overhang, right.left_overhang} == {""}


def test_a_site_put_into_a_record_is_found_there() -> None:
    edited, report = insert_site(SequenceRecord("AAAACCCC"), "BsaI", 4)
    assert edited.sequence == "AAAAGGTCTCCCCC"
    assert [(s.start, s.strand) for s in find_sites(edited, "BsaI")] == [(4, Strand.FORWARD)]
    assert report == EditReport()


def test_a_site_may_be_put_in_pointing_the_other_way() -> None:
    edited, _ = insert_site(SequenceRecord("AAAACCCC"), "BsaI", 4, strand=Strand.REVERSE)
    assert edited.sequence == "AAAAGAGACCCCCC"
    assert [(s.start, s.strand) for s in find_sites(edited, "BsaI")] == [(4, Strand.REVERSE)]


def test_putting_a_site_in_shifts_what_the_record_annotates() -> None:
    record = SequenceRecord(
        "AAAACCCC", features=(Feature("tag", "misc_feature", (Segment(4, 8),)),)
    )
    edited, _ = insert_site(record, "BsaI", 4)
    assert edited.features[0].segments == (Segment(10, 14),)


def test_a_site_holding_an_iupac_code_cannot_be_put_in_as_written() -> None:
    with pytest.raises(ValueError, match="FooI"):
        insert_site(SequenceRecord("AAAACCCC"), FOOI, 4)


@pytest.mark.parametrize(
    ("name", "overhang"),
    [("BsaI", "AATG"), ("BbsI", "AATG"), ("SapI", "ATG"), ("PaqCI", "AATG"), ("BtgZI", "AATG")],
)
def test_a_primer_tail_is_cut_to_leave_exactly_the_overhang_it_was_asked_for(
    name: str, overhang: str
) -> None:
    enzyme = get_enzyme(name)
    tail = primer_tail(enzyme, overhang)
    (site,) = find_sites(SequenceRecord(tail), enzyme)
    assert site.overhang == overhang
    assert tail.endswith(overhang)
    # Six flanking bases, then the site, then whatever the enzyme reaches over to cut.
    assert len(tail) == 6 + enzyme.bottom_cut


def test_a_primer_tail_holds_the_recognition_site_once() -> None:
    tail = primer_tail(get_enzyme("BsaI"), "AATG")
    assert tail.count("GGTCTC") == 1
    assert len(find_sites(SequenceRecord(tail), "BsaI")) == 1


def test_a_flank_that_spells_a_second_site_is_refused() -> None:
    with pytest.raises(ValueError, match="BsaI"):
        primer_tail(get_enzyme("BsaI"), "AATG", flank="GGTCTC")


def test_a_dcm_site_is_refused_only_where_the_supplier_says_dcm_impairs_the_enzyme() -> None:
    # NEB: BsaI is impaired by overlapping Dcm methylation, and BsmBI is not sensitive to it.
    with pytest.raises(ValueError, match=r"[Dd]cm"):
        primer_tail(get_enzyme("BsaI"), "AATG", flank="ACCAGG")
    assert primer_tail(get_enzyme("BsmBI"), "AATG", flank="ACCAGG").startswith("ACCAGG")


def test_an_overhang_the_enzyme_would_not_leave_is_refused() -> None:
    with pytest.raises(ValueError, match="4"):
        primer_tail(get_enzyme("BsaI"), "AAT")


def test_a_type_ii_enzyme_gets_its_site_and_no_overhang_to_choose() -> None:
    tail = primer_tail(get_enzyme("EcoRI"))
    assert tail.endswith("GAATTC")
    assert len(tail) == 6 + 6
    with pytest.raises(ValueError, match="overhang"):
        primer_tail(get_enzyme("EcoRI"), "AATG")


def _found(record: SequenceRecord, name: str) -> list[tuple[int, Strand]]:
    return [(site.start + 1, site.strand) for site in find_sites(record, name)]
