"""Whether two cut ends anneal, and the fidelity the shipped data scores a set of overhangs at.

Against Pryor 2020's own worked examples: the numbers the paper reports are the specification.
"""

from pathlib import Path

import pytest

from mbio.enzymes import EndType, Enzyme
from mbio.ligase import LigaseProfile
from mbio.overhangs import (
    MIN_DISTANCE,
    MODEST_MISMATCH,
    STRONG_LIGATION,
    End,
    compatible,
    fidelity,
    ligation_matrix,
    on_target,
    refusal,
)

# Pryor 2020's worked example: the eleven overhangs the plant synthetic biology community
# standardised on, which the paper scores at 81% with BsmBI-v2 and 42 C / 16 C cycling.
PLANT = ("GGAG", "TGAC", "TCCC", "TACT", "CCAT", "AATG", "AGCC", "TTCG", "GCTT", "GGTA", "CGCT")

# Potapov 2018, Table 1, set 1: the MoClo standard overhangs extended to fifteen, chosen for
# cycled assembly and reported there at 98.5%.
HIGH_FIDELITY = (
    "TGCC", "GCAA", "ACTA", "TTAC", "CAGA", "TGTG", "GAGC", "AGGA",
    "ATTC", "CGAA", "ATAG", "AAGG", "AACT", "AAAA", "ACCG",
)  # fmt: skip

# Every Type IIS enzyme the package ships leaves three or four bases, and shipped matrices cover
# both, so a stand-in always exists for one of them. This record is what an enzyme outside that
# range would be, and it is the only way left to reach the rules.
WIDE = Enzyme("five-base cutter", "CTGGAG", top_cut=20, bottom_cut=25)


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


def test_a_reserved_overhang_binds_as_a_taken_one_does_and_is_refused_as_reserved() -> None:
    # A reserved overhang is held out by its own name, its reverse complement with it, and it
    # keeps the distance rule the way a junction already taken would.
    held = refusal("AGGT", "BsaI", reserved=["ACCT"])
    near = refusal("AGGA", "BsaI", reserved=["AGGT"])

    assert held is not None
    assert near is not None
    assert (held.rule, near.rule) == ("reserved", "near-duplicate")


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
    # Named in passing, by author, year and table: a check detail is read at the bench, and the
    # reference in full belongs on the page that carries the number.
    assert report.source.startswith("Pryor 2020 S")
    assert "PLoS One" not in report.source
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


def test_an_enzyme_nobody_measured_is_scored_on_a_matrix_of_its_own_overhang_length() -> None:
    report = fidelity(("AATG", "GCTT", "TACA"), "PaqCI")

    assert report.measured
    assert not report.enzyme_specific
    # Esp3I is the four-base matrix with the most ligations behind it.
    assert report.label == "measured with Esp3I, not specific to PaqCI"
    assert "PaqCI" in report.source
    assert len(report.ligations) == 3


@pytest.mark.parametrize(("enzyme", "standing_in"), [("BspQI", "SapI"), ("BpiI", "BbsI-HF")])
def test_an_enzyme_reading_a_shipped_enzyme_s_site_takes_that_enzyme_s_matrix(
    enzyme: str, standing_in: str
) -> None:
    # An isoschizomer cuts the same site the same way, so the measurement is of both of them.
    overhangs = ("AAT", "GCT") if enzyme == "BspQI" else ("AATG", "GCTT")

    report = fidelity(overhangs, enzyme)

    assert report.measured
    assert not report.enzyme_specific
    assert report.label == f"measured with {standing_in}, not specific to {enzyme}"


def test_the_rules_score_an_enzyme_no_shipped_matrix_shares_an_overhang_length_with() -> None:
    report = fidelity(("AATGC", "GCTTA", "TACAG"), WIDE)

    assert not report.measured
    assert report.enzyme_specific
    assert report.label == "rule-based estimate"
    assert "rule" in report.source
    assert report.ligations == ()


def test_the_rule_based_fallback_marks_a_set_of_near_duplicates_down() -> None:
    spread = fidelity(("AATGC", "GCTTA", "TACAG"), WIDE).value
    crowded = fidelity(("AATGC", "AATCC", "TACAG"), WIDE).value

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


@pytest.mark.parametrize("name", ["BsaI", "BsmBI", "Esp3I", "BbsI", "SapI"])
def test_the_distance_rule_is_read_off_the_shipped_data_and_not_off_a_standard(name: str) -> None:
    """`MIN_DISTANCE` is where the measured mis-ligations stop, so the two may not drift apart.

    Two is the smallest separation at which every shipped matrix holds its cross-ligations
    under `MODEST_MISMATCH`. One base apart leaves hundreds of pairs at or above it. BsaI's one
    pair is the only one the rule lets through, and it sits far under what one base apart
    reaches; the matrix records it in both directions.
    """
    matrix = ligation_matrix(name)
    assert matrix is not None

    near: list[float] = []
    far: list[float] = []
    for row, columns in matrix.counts.items():
        for column in columns:
            seen = matrix.normalised(row, column)
            if column == _reverse(row) or seen < MODEST_MISMATCH:
                continue
            crowded = _apart(row, _reverse(column)) < MIN_DISTANCE
            (near if crowded else far).append(seen)

    assert len(far) == (2 if name == "BsaI" else 0)
    assert len(near) > 100
    assert max(far, default=0.0) < max(near) / 5


# Bilotti 2022's `File S7. T7 PEG`, correct Watson-Crick pairs per 100,000 events: the best of
# the overhangs the AP-1 example chose, and the worst T7 punishes that the rules still allow.
# Nothing of that workbook is here -- two numbers from a CC BY paper, cited, and the profile
# holding them is written by this file.
AGAT_ON_T7_PEG = 175.8
TAGA_ON_T7_PEG = 60.8


def t7_peg(**rates: float) -> LigaseProfile:
    """A sheet of a workbook the user holds, joining each overhang at the rate asked for.

    A rate is per 100,000 events, which the filler row makes exact.
    """
    counts = {one: {_reverse(one): round(rate * 10)} for one, rate in rates.items()}
    counts["AAAA"] = {"TTTT": 1_000_000 - sum(row[one] for row in counts.values() for one in row)}
    return LigaseProfile(
        Path("File S1_NAR.xlsx"),
        conditions="File S7. T7 PEG",
        overhang_length=4,
        observations=1_000_000,
        counts=counts,
    )


def test_an_overhang_above_the_floor_is_silent_and_one_below_it_is_named() -> None:
    profile = t7_peg(AGAT=AGAT_ON_T7_PEG, TAGA=TAGA_ON_T7_PEG)

    report = on_target(("AGAT", "TAGA"), profile)

    assert [round(one.rate, 1) for one in report.rates] == [AGAT_ON_T7_PEG, TAGA_ON_T7_PEG]
    assert report.weak == ("TAGA",)
    assert report.floor == STRONG_LIGATION


def test_the_report_says_which_profile_and_which_sheet_the_rates_came_from() -> None:
    report = on_target(("AGAT",), t7_peg(AGAT=AGAT_ON_T7_PEG))

    assert report.label == "measured on File S7. T7 PEG, read from File S1_NAR.xlsx"


def test_a_caller_may_hold_a_set_to_a_floor_of_its_own() -> None:
    profile = t7_peg(AGAT=AGAT_ON_T7_PEG)

    assert on_target(("AGAT",), profile, floor=200).weak == ("AGAT",)


def test_an_overhang_the_profile_does_not_measure_is_refused_by_length() -> None:
    with pytest.raises(ValueError, match="4-base overhangs"):
        on_target(("AGATC",), t7_peg(AGAT=AGAT_ON_T7_PEG))


def _apart(one: str, other: str) -> int:
    """How far two overhangs of a set stand, which counts each one's reverse complement too."""
    return min(
        sum(a != b for a, b in zip(one, partner, strict=True))
        for partner in (other, _reverse(other))
    )


def _reverse(overhang: str) -> str:
    """The other strand of an overhang, which is the column its row pairs with."""
    return overhang.translate(str.maketrans("ACGT", "TGCA"))[::-1]
