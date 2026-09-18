"""Choosing the enzyme, and designing the overhangs its junctions are cut to."""

import pytest

from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.goldengate.design import (
    GOLDEN_GATE_ENZYMES,
    LAST_RESORT,
    choose_enzyme,
    design_overhangs,
)
from liulab_mbio.overhangs import Junction, ligation_matrix
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.sites import find_sites, primer_tail


def test_the_fixtures_leave_bbsi_and_paqci_free(puc19: SequenceRecord, gfp: SequenceRecord) -> None:
    free = [choice.enzyme.name for choice in choose_enzyme([puc19, gfp]) if choice.free]

    assert set(free) >= {"BbsI", "PaqCI"}
    assert "BsaI" not in free


def test_the_chosen_enzyme_is_free_of_sites_in_the_parts(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    chosen = choose_enzyme([puc19, gfp])[0]

    assert chosen.enzyme.name == "BbsI"
    assert (chosen.sites, chosen.free) == (0, True)
    assert find_sites(puc19, chosen.enzyme) == ()
    assert find_sites(gfp, chosen.enzyme) == ()


def test_a_tie_between_free_enzymes_goes_to_the_one_whose_overhangs_can_be_scored(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    ranked = [choice.enzyme.name for choice in choose_enzyme([puc19, gfp]) if choice.free]

    assert ranked.index("BbsI") < ranked.index("PaqCI")
    assert ligation_matrix("BbsI") is not None
    assert ligation_matrix("PaqCI") is None


def test_every_free_enzyme_ranks_ahead_of_every_enzyme_with_a_site(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    ranked = choose_enzyme([puc19, gfp])

    sites = [choice.sites for choice in ranked if choice.enzyme.name not in LAST_RESORT]
    assert sites == sorted(sites)


def test_no_domestication_is_proposed_while_an_enzyme_is_free(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    ranked = choose_enzyme([puc19, gfp])

    assert not any(choice.changes or choice.outside_cds or choice.unchanged for choice in ranked)


def test_domestication_is_proposed_when_no_enzyme_is_free(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    (choice,) = choose_enzyme([puc19, gfp], enzymes=["BsaI"])

    assert (choice.sites, choice.free) == (2, False)
    # One codon in GFP and one in AmpR, both synonymous, so the site goes and the protein stays.
    assert {change.amino_acid for change in choice.changes} == {"D", "G"}
    assert (choice.outside_cds, choice.unchanged) == ((), ())
    assert choice.clean


def test_a_site_lying_in_no_coding_sequence_is_reported_rather_than_edited(
    puc19: SequenceRecord,
) -> None:
    (choice,) = choose_enzyme([puc19], enzymes=["BsmBI"])

    assert choice.changes == ()
    assert len(choice.outside_cds) == 2
    assert not choice.clean


def test_btgzi_ranks_last_because_neb_publishes_no_golden_gate_protocol_for_it() -> None:
    ranked = choose_enzyme([SequenceRecord("ACGT" * 40)])

    assert all(choice.free for choice in ranked)
    assert ranked[-1].enzyme.name == "BtgZI"
    assert "BtgZI" in LAST_RESORT
    assert "BtgZI" in GOLDEN_GATE_ENZYMES


def test_only_an_enzyme_cutting_outside_its_site_joins_the_golden_gate_set() -> None:
    # SrfI and PmeI ship for the library scheme, and cut inside their own site.
    assert {"SrfI", "PmeI"}.isdisjoint(GOLDEN_GATE_ENZYMES)
    assert all(get_enzyme(name).type == "IIS" for name in GOLDEN_GATE_ENZYMES)


def test_a_scarless_junction_takes_the_overhang_the_record_already_spells() -> None:
    record = SequenceRecord("AAAACCTGAGGGGTTTT")
    junction = Junction("insert", record=record, position=4, scarless=True)

    chosen = design_overhangs([junction], "BsaI")

    assert chosen.overhangs == ("CCTG",)
    assert chosen.choices[0].offset == 0


def test_a_scarless_junction_slides_within_its_window_to_escape_a_rule() -> None:
    # `GGCC` at 4 is a palindrome, so the junction moves by one base rather than failing.
    record = SequenceRecord("AAAAGGCCTTTTAAAA")
    junction = Junction("insert", record=record, position=4, scarless=True, window=2)

    chosen = design_overhangs([junction], "BsaI")

    assert chosen.overhangs != ("GGCC",)
    assert chosen.choices[0].offset != 0
    assert [rejection.rule for rejection in chosen.choices[0].rejected] == ["palindrome"]


def test_a_palindrome_is_refused_because_a_fragment_would_ligate_to_itself() -> None:
    junction = Junction("fixed", overhang="AATT")

    with pytest.raises(ValueError, match="palindrome"):
        design_overhangs([junction], "BsaI")


def test_an_overhang_already_in_the_set_is_refused_as_a_repeat() -> None:
    with pytest.raises(ValueError, match="repeat"):
        design_overhangs(
            [Junction("one", overhang="AATG"), Junction("two", overhang="AATG")], "BsaI"
        )


def test_the_reverse_complement_of_a_chosen_overhang_is_the_same_junction_twice() -> None:
    with pytest.raises(ValueError, match="repeat"):
        design_overhangs(
            [Junction("one", overhang="AATG"), Junction("two", overhang="CATT")], "BsaI"
        )


def test_a_near_duplicate_is_refused_and_the_distance_rule_is_a_dial() -> None:
    junctions = [Junction("one", overhang="AATG"), Junction("two", overhang="AATC")]

    with pytest.raises(ValueError, match="near-duplicate"):
        design_overhangs(junctions, "BsaI")

    assert design_overhangs(junctions, "BsaI", min_distance=1).overhangs == ("AATG", "AATC")


def test_an_overhang_of_one_base_kind_is_refused_unless_the_caller_allows_it() -> None:
    junction = Junction("fixed", overhang="GCCG")

    with pytest.raises(ValueError, match="uniform"):
        design_overhangs([junction], "BsaI")

    assert design_overhangs([junction], "BsaI", allow_uniform=True).overhangs == ("GCCG",)


def test_an_overhang_of_the_wrong_length_for_the_enzyme_is_refused() -> None:
    with pytest.raises(ValueError, match="length"):
        design_overhangs([Junction("fixed", overhang="AAT")], "BsaI")


def test_sapi_designs_the_three_base_overhangs_it_leaves() -> None:
    junctions = [Junction("one"), Junction("two")]

    chosen = design_overhangs(junctions, "SapI")

    assert {len(overhang) for overhang in chosen.overhangs} == {3}
    assert len(set(chosen.overhangs)) == 2


def test_a_free_junction_is_given_an_overhang_that_passes_every_rule() -> None:
    chosen = design_overhangs([Junction("one"), Junction("two"), Junction("three")], "BbsI")

    assert len(set(chosen.overhangs)) == 3
    for overhang in chosen.overhangs:
        assert len(overhang) == 4
        assert primer_tail("BbsI", overhang)


def test_an_in_frame_junction_moves_only_by_whole_codons() -> None:
    # A coding sequence whose junction bases at the codon boundary spell a palindrome, so the
    # junction has to move three bases rather than one to stay in frame.
    record = SequenceRecord(
        "ATGAAAGGCCTTGAAGATCCGTATAAA",
        features=(Feature("cds", "CDS", (Segment(0, 27),), strand=Strand.FORWARD),),
    )
    junction = Junction("fusion", record=record, position=6, in_frame=True, window=6)

    chosen = design_overhangs([junction], "BsaI")

    assert chosen.choices[0].offset % 3 == 0
    assert chosen.choices[0].offset != 0
    assert chosen.overhangs[0] == record.sequence[6 + chosen.choices[0].offset :][:4]


def test_an_in_frame_junction_off_a_codon_boundary_is_refused() -> None:
    record = SequenceRecord(
        "ATGAAAGAACTTGAAGATCCGTATAAA",
        features=(Feature("cds", "CDS", (Segment(0, 27),), strand=Strand.FORWARD),),
    )

    with pytest.raises(ValueError, match="codon boundary"):
        design_overhangs([Junction("fusion", record=record, position=7, in_frame=True)], "BsaI")


def test_a_junction_with_no_acceptable_overhang_names_the_rules_that_refused_it() -> None:
    record = SequenceRecord("AAAAGGCCTTTT")
    junction = Junction("stuck", record=record, position=4, scarless=True)

    with pytest.raises(ValueError, match="palindrome"):
        design_overhangs([junction], "BsaI")


def test_a_rejection_names_the_candidate_the_rule_and_why() -> None:
    # The junction sits where the record spells the overhang the other junction already took.
    record = SequenceRecord("AAAAAGGCTTTTAAAA")
    junctions = [
        Junction("one", overhang="AGGC"),
        Junction("two", record=record, position=4, scarless=True, window=1),
    ]

    chosen = design_overhangs(junctions, "BsaI")

    (rejection,) = chosen.choices[1].rejected
    assert (rejection.overhang, rejection.rule) == ("AGGC", "repeat")
    assert "AGGC" in rejection.detail
    assert chosen.choices[1].offset == -1


def test_the_designed_set_carries_its_own_fidelity(
    puc19: SequenceRecord, gfp: SequenceRecord
) -> None:
    enzyme = choose_enzyme([puc19, gfp])[0].enzyme

    chosen = design_overhangs([Junction("left"), Junction("right")], enzyme)

    assert chosen.fidelity.measured
    assert chosen.fidelity.value > 0.9
    assert chosen.enzyme is get_enzyme("BbsI")
