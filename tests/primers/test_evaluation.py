"""Pinned Tm and Ta values come from NEB's own calculator and API, through
``docs/research/primer-design-and-pcr.md``.
"""

import pytest

from liulab_mbio.primers import (
    ONETAQ,
    design_pair,
    evaluate_pair,
    evaluate_primer,
    melting_temperature,
)
from liulab_mbio.sequence import BindingSite, Primer, SequenceRecord, Strand, reverse_complement

from .sequences import M13_FWD, M13_REV, MCS_FWD, PUC_FWD, PUC_REV


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


def test_a_check_no_sourced_threshold_judges_carries_no_verdict() -> None:
    # primer3's own value for a 3' end of CCAGT, as a delta G.
    report = evaluate_primer(Primer("M13 fwd", M13_FWD))
    assert report["end_stability"].value == pytest.approx(-4.0, abs=0.01)
    assert report["end_stability"].status is None
    assert report["tm_full"].status is None
    # The report's verdict is the worst of the checks something judged.
    assert report.status == "warn"


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


def test_a_primer_that_binds_nowhere_fails(puc19) -> None:
    report = evaluate_primer(Primer("elsewhere", "ATGAGTAAAGGAGAAGAACTTTTC"), puc19)
    assert report["binding_sites"].value == 0
    assert report["binding_sites"].status == "fail"


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
