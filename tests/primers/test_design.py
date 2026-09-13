import pytest

from liulab_mbio.primers import (
    THRESHOLDS,
    THRESHOLDS_FOR,
    design_pair,
    design_primer,
    evaluate_primer,
    melting_temperature,
)
from liulab_mbio.sequence import BindingSite, Segment, SequenceRecord, Strand, reverse_complement

from .sequences import MCS_FWD


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


def test_a_design_keeps_inside_the_length_band_it_is_judged_by(puc19) -> None:
    # No region 18 bases or longer reaches the Tm band here; a 16-mer lands on the target.
    for thresholds in (THRESHOLDS, THRESHOLDS_FOR["sequencing"]):
        primer = design_primer(puc19, 32, Strand.FORWARD, thresholds=thresholds)
        assert evaluate_primer(primer, thresholds=thresholds)["length"].status == "pass"
    forward, reverse = design_pair(puc19, 32, 232)
    assert [evaluate_primer(one)["length"].status for one in (forward, reverse)] == ["pass"] * 2
    # A sequencing primer may be 16 bases, so it keeps that one.
    sequencing = design_primer(puc19, 32, Strand.FORWARD, thresholds=THRESHOLDS_FOR["sequencing"])
    assert sequencing.sequence == "TGACACATGCAGCTCC"


def test_a_design_keeps_inside_the_bands_its_length_decides_where_a_length_can(puc19) -> None:
    decided = ("length", "gc_percent", "gc_clamp", "tm")
    # The 20-mer nearest the target Tm ends on four Gs and Cs; the 19-mer passes.
    primer = evaluate_primer(design_primer(puc19, 210, Strand.FORWARD))
    assert [primer[name].status for name in decided] == ["pass"] * 4
    # The reverse primer warns at every length here, and does not drag its partner with it.
    forward, reverse = (evaluate_primer(one) for one in design_pair(puc19, 455, 395 + len(puc19)))
    assert [forward[name].status for name in decided] == ["pass"] * 4
    assert reverse.status == "warn"


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
