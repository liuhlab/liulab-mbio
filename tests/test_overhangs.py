"""Whether two cut ends anneal, and the fidelity the shipped data scores a set of overhangs at.

Against Pryor 2020's own worked examples: the numbers the paper reports are the specification.
"""

import pytest

from liulab_mbio.enzymes import EndType
from liulab_mbio.overhangs import STRONG_LIGATION, End, compatible, fidelity, ligation_matrix

# Pryor 2020's worked example: the eleven overhangs the plant synthetic biology community
# standardised on, which the paper scores at 81% with BsmBI-v2 and 42 C / 16 C cycling.
PLANT = ("GGAG", "TGAC", "TCCC", "TACT", "CCAT", "AATG", "AGCC", "TTCG", "GCTT", "GGTA", "CGCT")

# Potapov 2018, Table 1, set 1: the MoClo standard overhangs extended to fifteen, chosen for
# cycled assembly and reported there at 98.5%.
HIGH_FIDELITY = (
    "TGCC", "GCAA", "ACTA", "TTAC", "CAGA", "TGTG", "GAGC", "AGGA",
    "ATTC", "CGAA", "ATAG", "AAGG", "AACT", "AAAA", "ACCG",
)  # fmt: skip


def test_two_enzymes_leaving_the_same_overhang_anneal() -> None:
    # SalI and XhoI both leave a 5' TCGA, which is what lets either cut into the other's site.
    assert compatible(End.cut_by("SalI", "TCGA"), End.cut_by("XhoI", "TCGA"))


def test_the_same_bases_on_opposite_strands_do_not_anneal() -> None:
    # HindIII leaves AGCT on the top strand and SacI leaves AGCT on the bottom.
    assert not compatible(End.cut_by("HindIII", "AGCT"), End.cut_by("SacI", "AGCT"))


def test_two_unequal_overhangs_do_not_anneal() -> None:
    assert not compatible(End.cut_by("EcoRI", "AATT"), End.cut_by("BamHI", "GATC"))


def test_two_blunt_ends_anneal() -> None:
    assert compatible(End.cut_by("SmaI", ""), End.cut_by("PmeI", ""))


def test_a_blunt_end_does_not_anneal_to_an_overhang() -> None:
    # SmaI and XmaI read one site and cut it in different places.
    assert not compatible(End.cut_by("SmaI", ""), End.cut_by("XmaI", "CCGG"))


def test_an_overhang_the_enzyme_would_not_leave_is_refused() -> None:
    with pytest.raises(ValueError, match="4-base"):
        End.cut_by("EcoRI", "AAT")


@pytest.mark.parametrize(("overhang", "end_type"), [("AATT", "blunt"), ("", "5'")])
def test_an_end_whose_bases_contradict_its_end_type_is_refused(
    overhang: str, end_type: EndType
) -> None:
    with pytest.raises(ValueError, match="cannot leave"):
        End(overhang, end_type)


def test_a_watson_crick_pair_is_a_row_against_the_reverse_complement_of_a_column() -> None:
    matrix = ligation_matrix("BsaI")
    assert matrix is not None

    assert matrix.count("TTTT", "AAAA") == matrix.count("AAAA", "TTTT")
    assert matrix.count("TTTT", "AAAA") > matrix.count("TTTT", "TTTT")
    assert (matrix.overhang_length, len(matrix.overhangs)) == (4, 256)


def test_the_published_plant_overhang_set_scores_as_the_paper_reports_it() -> None:
    report = fidelity(PLANT, "BsmBI")

    assert report.measured
    assert "Pryor" in report.source
    # Pryor 2020, Fig 4A: 81% for this set with BsmBI-v2 and 42 C / 16 C cycling.
    assert round(report.value, 2) == 0.81


def test_the_mispair_the_paper_blames_is_the_one_reported_worst() -> None:
    report = fidelity(PLANT, "BsmBI")

    worst = {report.mismatches[0][:2], report.mismatches[1][:2]}
    assert worst == {("GGTA", "TACT"), ("TACT", "GGTA")}
    # Pryor 2020, Fig 4A: dropping that mispair takes the same set to 92%.
    assert round(fidelity([one for one in PLANT if one != "GGTA"], "BsmBI").value, 2) == 0.92


def test_the_set_getset_extended_to_twenty_scores_as_the_paper_reports_it() -> None:
    # Pryor 2020, Fig 4B: the nine GetSet added take the plant set from 81% to 80%.
    extended = (*PLANT, "ACCT", "CCGC", "ACAA", "AACA", "GAAA", "CAAG", "GCAC", "TAGA", "AAAT")

    assert round(fidelity(extended, "BsmBI").value, 2) == 0.80


def test_a_high_fidelity_set_scores_above_a_set_of_near_duplicates() -> None:
    near = ("AAAA", "AAAC", "AAAG", "AAAT", "AACA")

    assert fidelity(HIGH_FIDELITY, "BsaI").value > fidelity(near, "BsaI").value


def test_sapi_is_scored_on_the_three_base_matrix() -> None:
    report = fidelity(("AAA", "CGT", "TGC"), "SapI")

    assert report.measured
    assert 0.0 < report.value <= 1.0


def test_an_overhang_the_enzyme_could_not_leave_is_refused_by_the_scorer() -> None:
    with pytest.raises(ValueError, match="3-base"):
        fidelity(("AAAA", "CCGT"), "SapI")


def test_an_enzyme_with_no_matrix_falls_back_to_the_rules_and_says_so() -> None:
    report = fidelity(("AATG", "GCTT", "TACA"), "PaqCI")

    assert not report.measured
    assert report.enzyme_specific
    assert report.label == "rule-based estimate"
    assert "rule" in report.source
    assert report.ligations == ()


def test_the_rule_based_fallback_marks_a_set_of_near_duplicates_down() -> None:
    spread = fidelity(("AATG", "GCTT", "TACA"), "PaqCI").value
    crowded = fidelity(("AATG", "AATC", "TACA"), "PaqCI").value

    assert crowded < spread <= 1.0


@pytest.mark.parametrize("name", ["BsaI", "BsmBI", "Esp3I", "BbsI", "SapI"])
def test_every_watson_crick_pair_the_shipped_data_covers_ligates_strongly(name: str) -> None:
    # So `FidelityReport.weak` is empty for any set scored on shipped data, and a weak pair
    # would be news rather than routine.
    matrix = ligation_matrix(name)
    assert matrix is not None

    lowest = min(matrix.normalised(one, _reverse(one)) for one in matrix.overhangs)

    assert lowest >= STRONG_LIGATION
    assert fidelity(matrix.overhangs[:4], name).weak == ()


def _reverse(overhang: str) -> str:
    """The other strand of an overhang, which is the column its row pairs with."""
    return overhang.translate(str.maketrans("ACGT", "TGCA"))[::-1]
