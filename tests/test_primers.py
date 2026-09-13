"""Pinned Tm and Ta values come from NEB's own calculator and API, through
``docs/research/primer-design-and-pcr.md``.
"""

import dataclasses
from pathlib import Path

import pytest

from liulab_mbio.primers import (
    ONETAQ,
    PHUSION,
    Q5,
    TAQ,
    THRESHOLDS,
    Band,
    design_pair,
    design_primer,
    evaluate_pair,
    evaluate_primer,
    find_binding_sites,
    find_priming_sites,
    melting_temperature,
)
from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)


@pytest.fixture(scope="module")
def puc19() -> SequenceRecord:
    from Bio import SeqIO

    record = SeqIO.read(Path(__file__).parent / "data" / "pUC19.dna", "snapgene")
    return SequenceRecord(str(record.seq), topology="circular", name="pUC19")


M13_FWD = "GTAAAACGACGGCCAGT"
M13_REV = "CAGGAAACAGCTATGAC"
#: The 23-mer M13/pUC pair, which colony PCR uses instead of the 17-mers.
PUC_FWD = "CCCAGTCACGACGTTGTAAAACG"
PUC_REV = "AGCGGATAACAATTTCACACAGG"
#: Anneals to pUC19 just after the MCS, reading along the top strand.
MCS_FWD = "GGCGTAATCATGGTCATAGC"


@pytest.mark.parametrize(
    ("polymerase", "forward_tm", "reverse_tm"),
    [
        (Q5, 62.27, 56.31),
        (PHUSION, 56.57, 50.89),
        (TAQ, 53.96, 47.99),
        (ONETAQ, 53.82, 47.85),
    ],
    ids=lambda value: getattr(value, "name", value),
)
def test_tm_matches_nebs_calculator_for_the_polymerase(polymerase, forward_tm, reverse_tm) -> None:
    assert melting_temperature(M13_FWD, polymerase) == pytest.approx(forward_tm, abs=0.01)
    assert melting_temperature(M13_REV, polymerase) == pytest.approx(reverse_tm, abs=0.01)


def test_tm_defaults_to_q5() -> None:
    assert melting_temperature(M13_REV) == melting_temperature(M13_REV, Q5)


@pytest.mark.parametrize(
    ("polymerase", "annealing_temperature"),
    [(Q5, 57.3), (PHUSION, 54.8), (TAQ, 43.0), (ONETAQ, 42.8)],
    ids=lambda value: getattr(value, "name", value),
)
def test_annealing_temperature_follows_nebs_rule_for_the_polymerase(
    polymerase, annealing_temperature
) -> None:
    tms = (melting_temperature(M13_FWD, polymerase), melting_temperature(M13_REV, polymerase))
    assert polymerase.annealing_temperature(*tms) == pytest.approx(annealing_temperature, abs=0.1)


def test_the_annealing_temperature_is_capped() -> None:
    long_pair = ("CGCCAGGGTTTTCCCAGTCACGACGTTG", "GCGGATAACAATTTCACACAGGAAACAGCTATGAC")
    tms = tuple(melting_temperature(primer, Q5) for primer in long_pair)
    assert Q5.annealing_temperature(*tms) == 72.0
    assert TAQ.annealing_temperature(80.0, 75.0) == 68.0


def test_the_twenty_three_mer_pair_anneals_above_nebs_floor_in_onetaq() -> None:
    tms = (melting_temperature(PUC_FWD, ONETAQ), melting_temperature(PUC_REV, ONETAQ))
    assert ONETAQ.annealing_temperature(*tms) == pytest.approx(51.5, abs=0.1)
    assert ONETAQ.annealing_temperature(*tms) > ONETAQ.annealing_min


def test_extension_time_rounds_the_amplicon_up_to_whole_kilobases() -> None:
    assert Q5.extension_seconds(2686) == 60
    assert Q5.extension_seconds(1000) == 20
    assert PHUSION.extension_seconds(1000) == 15
    assert TAQ.extension_seconds(103) == 60
    assert ONETAQ.extension_seconds(1001) == 120


def test_a_primer_without_a_tail_anneals_along_its_whole_length() -> None:
    report = evaluate_primer(Primer("M13 fwd", M13_FWD))
    assert report["length"].value == 17
    assert report["length"].status == "warn"
    assert report["gc_percent"].value == pytest.approx(52.9, abs=0.1)
    assert report["gc_percent"].status == "pass"
    assert report["gc_clamp"].value == 3
    assert report["gc_clamp"].status == "pass"
    assert report["tm"].value == pytest.approx(62.27, abs=0.01)
    assert report["tm"].status == "pass"
    assert report["tm_full"].value == pytest.approx(report["tm"].value)
    # The worst of its checks: the annealing region is shorter than primer3's minimum.
    assert report.status == "warn"


def test_a_tail_raises_the_full_primer_tm_and_leaves_the_annealing_region_alone() -> None:
    tailed = Primer(
        "MCS fwd",
        "TTGGTCTCA" + MCS_FWD,
        binding_sites=(BindingSite(452, 452 + len(MCS_FWD), Strand.FORWARD),),
    )
    report = evaluate_primer(tailed)
    assert report["length"].value == 20
    assert report["length"].status == "pass"
    assert report["tm"].value == pytest.approx(melting_temperature(MCS_FWD))
    assert report["tm_full"].value > report["tm"].value + 5


def test_end_stability_is_reported_but_not_judged() -> None:
    # primer3's own value for a 3' end of CCAGT, as a delta G.
    check = evaluate_primer(Primer("M13 fwd", M13_FWD))["end_stability"]
    assert check.value == pytest.approx(-4.0, abs=0.01)
    assert check.status == "pass"
    assert "not judged" in check.detail


def test_runs_and_repeats_are_counted_over_the_whole_primer() -> None:
    runs = evaluate_primer(Primer("run", "GCTAAAAAAAGCTGCTGATCG"))
    assert runs["mononucleotide_run"].value == 7
    assert runs["mononucleotide_run"].status == "fail"
    guanines = evaluate_primer(Primer("g run", "GCTGGGGATCGATCGATCGA"))
    assert guanines["mononucleotide_run"].value == 4
    assert guanines["mononucleotide_run"].status == "warn"
    repeats = evaluate_primer(Primer("repeat", "GCGATATATATATATGCG"))
    assert repeats["dinucleotide_repeat"].value == 6
    assert repeats["dinucleotide_repeat"].status == "warn"
    assert repeats["mononucleotide_run"].value == 1


def test_a_hairpin_above_the_threshold_warns() -> None:
    report = evaluate_primer(Primer("hairpin", "GCGCGCGAAAAACGCGCGC"))
    assert report["hairpin"].value == pytest.approx(88.65, abs=0.01)
    assert report["hairpin"].status == "warn"
    assert evaluate_primer(Primer("M13 fwd", M13_FWD))["hairpin"].status == "pass"


def test_a_three_prime_anchored_self_dimer_fails_where_an_internal_one_warns() -> None:
    anchored = evaluate_primer(Primer("anchored", "ATCGATCGATCAGCGCGCGCGC"))
    assert anchored["self_dimer_3prime"].value == pytest.approx(57.0, abs=0.01)
    assert anchored["self_dimer_3prime"].status == "fail"
    assert anchored.status == "fail"
    internal = evaluate_primer(Primer("internal", "GCGCGCGCGCTTTTTTTTTTTTTTT"))
    assert internal["self_dimer_3prime"].status == "pass"
    assert internal["self_dimer"].status == "warn"


def test_a_primer_longer_than_primer3_allows_is_judged_on_its_three_prime_end() -> None:
    long_tail = "GCTAGCTGACTGACTGATCGATCGATCGTAGCTAGCTGATCGATCGATGCTAGCTGA"
    report = evaluate_primer(Primer("long", long_tail + MCS_FWD))
    assert "60" in report["hairpin"].detail


def test_a_primer_is_placed_on_a_template_by_its_three_prime_match(puc19) -> None:
    assert find_binding_sites(M13_FWD, puc19) == (BindingSite(378, 395, Strand.FORWARD),)
    assert find_binding_sites(M13_REV, puc19) == (BindingSite(464, 481, Strand.REVERSE),)
    tailed = "TTGAAGACAA" + MCS_FWD
    assert find_binding_sites(tailed, puc19) == (BindingSite(452, 472, Strand.FORWARD),)


def test_a_binding_site_may_run_across_the_origin(puc19) -> None:
    across = puc19.sequence[-10:] + puc19.sequence[:10]
    assert find_binding_sites(across, puc19) == (BindingSite(2676, 2696, Strand.FORWARD),)


def test_a_primer_that_binds_one_place_has_no_off_target(puc19) -> None:
    report = evaluate_primer(Primer("M13 fwd", M13_FWD), puc19)
    assert report["binding_sites"].value == 1
    assert report["binding_sites"].status == "pass"
    assert report["off_target"].value == 0
    assert report["off_target"].status == "pass"
    assert report.status == "warn"  # its length, as without a template


def test_a_second_copy_of_the_annealing_region_is_an_off_target_site() -> None:
    template = SequenceRecord(
        "A" * 30 + MCS_FWD + "C" * 30 + reverse_complement(MCS_FWD) + "T" * 30,
        topology="circular",
    )
    primer = Primer("MCS fwd", MCS_FWD, binding_sites=(BindingSite(30, 50, Strand.FORWARD),))
    report = evaluate_primer(primer, template)
    off_target = report["off_target"]
    assert off_target.value == 1
    assert off_target.status == "warn"
    assert "80" in off_target.detail


def test_priming_sites_carry_their_strand_mismatches_and_tm(puc19) -> None:
    sites = find_priming_sites(M13_FWD, puc19)
    assert [site.site for site in sites] == [BindingSite(378, 395, Strand.FORWARD)]
    assert sites[0].mismatches == 0
    assert sites[0].tm == pytest.approx(54.61, abs=0.01)


def test_a_primer_that_binds_nowhere_fails(puc19) -> None:
    report = evaluate_primer(Primer("elsewhere", "ATGAGTAAAGGAGAAGAACTTTTC"), puc19)
    assert report["binding_sites"].value == 0
    assert report["binding_sites"].status == "fail"


def test_a_forward_primer_is_designed_from_a_position_to_a_target_tm(puc19) -> None:
    primer = design_primer(puc19, 452, Strand.FORWARD, name="MCS fwd")
    assert primer.name == "MCS fwd"
    assert primer.sequence == MCS_FWD
    assert primer.binding_sites == (BindingSite(452, 472, Strand.FORWARD),)
    assert melting_temperature(primer.sequence) == pytest.approx(62.0, abs=2.0)


def test_a_reverse_primer_reads_back_from_its_position(puc19) -> None:
    primer = design_primer(puc19, 396, Strand.REVERSE)
    site = primer.binding_sites[0]
    assert site.strand is Strand.REVERSE
    assert site.end == 396
    assert primer.sequence == reverse_complement(puc19.extract(Segment(site.start, site.end)))
    assert melting_temperature(primer.sequence) == pytest.approx(62.0, abs=2.0)


def test_a_tail_stays_outside_the_binding_site(puc19) -> None:
    tail = "TTGGTCTCAAATG"
    primer = design_primer(puc19, 452, Strand.FORWARD, tail=tail)
    site = primer.binding_sites[0]
    assert primer.sequence == tail + MCS_FWD
    assert site.end - site.start == len(primer.sequence) - len(tail)


def test_a_primer_may_be_designed_across_the_origin(puc19) -> None:
    primer = design_primer(puc19, len(puc19) - 8, Strand.FORWARD)
    site = primer.binding_sites[0]
    assert site.start == len(puc19) - 8
    assert site.end > len(puc19)
    assert primer.sequence == puc19.extract(Segment(site.start, site.end))


def test_a_pair_is_designed_with_matched_tms(puc19) -> None:
    forward, reverse = design_pair(puc19, 378, 481)
    assert forward.binding_sites[0].start == 378
    assert reverse.binding_sites[0].end == 481
    tms = [melting_temperature(primer.sequence) for primer in (forward, reverse)]
    assert abs(tms[0] - tms[1]) <= 2.0


def test_a_pair_may_amplify_across_the_origin(puc19) -> None:
    forward, reverse = design_pair(puc19, 452, 396 + len(puc19))
    assert forward.binding_sites[0].start == 452
    assert reverse.binding_sites[0].end % len(puc19) == 396


def test_design_refuses_a_position_a_linear_template_cannot_hold() -> None:
    template = SequenceRecord("ACGT" * 10)
    with pytest.raises(ValueError, match="fit"):
        design_primer(template, 38, Strand.FORWARD)


def test_a_pair_is_judged_with_its_amplicon(puc19) -> None:
    report = evaluate_pair(Primer("M13 fwd", M13_FWD), Primer("M13 rev", M13_REV), puc19)
    assert report["amplicon_size"].value == 103
    assert report["products"].value == 1
    assert report["products"].status == "pass"
    assert report["tm_difference"].value == pytest.approx(62.27 - 56.31, abs=0.02)
    assert report["tm_difference"].status == "warn"
    assert report["heterodimer"].status == "pass"
    assert report.amplicon_length == 103
    assert report.annealing_temperature == pytest.approx(57.3, abs=0.1)
    assert report.extension_seconds == 20


def test_the_colony_pcr_pair_gives_its_empty_vector_band(puc19) -> None:
    report = evaluate_pair(
        Primer("M13/pUC fwd", PUC_FWD), Primer("M13/pUC rev", PUC_REV), puc19, polymerase=ONETAQ
    )
    assert report["amplicon_size"].value == 137
    assert report["tm_difference"].status == "pass"
    assert report.annealing_temperature == pytest.approx(51.5, abs=0.1)
    assert report.extension_seconds == 60


def test_an_amplicon_may_run_across_the_origin(puc19) -> None:
    forward, reverse = design_pair(puc19, 452, 396 + len(puc19))
    report = evaluate_pair(forward, reverse, puc19)
    assert report.amplicon_length == 396 + len(puc19) - 452
    assert report["products"].status == "pass"


def test_tails_count_towards_the_amplicon(puc19) -> None:
    forward, reverse = design_pair(
        puc19, 378, 481, forward_tail="TTGGTCTCA", reverse_tail="TTGGTCTCA"
    )
    report = evaluate_pair(forward, reverse, puc19)
    assert report.amplicon_length == 103 + 9 + 9


def test_a_three_prime_anchored_heterodimer_fails(puc19) -> None:
    primer = Primer("dimer", "ATCGATCGATCAGCGCGCGCGC")
    report = evaluate_pair(primer, primer, puc19)
    assert report["heterodimer_3prime"].value == pytest.approx(57.0, abs=0.01)
    assert report["heterodimer_3prime"].status == "fail"


def test_a_pair_that_cannot_face_each_other_makes_no_product() -> None:
    template = SequenceRecord(MCS_FWD + "TTTT" + PUC_FWD)
    report = evaluate_pair(Primer("one", MCS_FWD), Primer("two", PUC_FWD), template)
    assert report["products"].value == 0
    assert report["products"].status == "fail"
    assert report.amplicon_length is None
    assert report.extension_seconds is None
    assert report.status == "fail"


def test_every_threshold_comes_from_one_place() -> None:
    strict = dataclasses.replace(THRESHOLDS, tm=Band(63.0, 64.0, 62.5, 65.0))
    report = evaluate_primer(Primer("M13 fwd", M13_FWD), thresholds=strict)
    assert report["tm"].status == "fail"
