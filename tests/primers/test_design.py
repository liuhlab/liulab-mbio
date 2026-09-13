import dataclasses

import pytest

from liulab_mbio.checks import STATUSES
from liulab_mbio.primers import (
    THRESHOLDS,
    THRESHOLDS_FOR,
    Band,
    PairReport,
    PrimerReport,
    design_pair,
    design_primer,
    evaluate_pair,
    evaluate_primer,
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


def test_a_design_gives_up_the_length_band_to_warn_on_fewer_checks(puc19) -> None:
    # Every region 18 bases or longer here warns on three checks or fails on Tm. The 16-mer
    # warns on two, and the band a length passes on ranks below the count of warnings.
    primer = design_primer(puc19, 32, Strand.FORWARD)
    assert primer.sequence == "TGACACATGCAGCTCC"
    assert warned(evaluate_primer(primer, puc19)) == ["length", "gc_clamp"]
    # A sequencing primer may be 16 bases, so the same region is inside its band.
    thresholds = THRESHOLDS_FOR["sequencing"]
    sequencing = design_primer(puc19, 32, Strand.FORWARD, thresholds=thresholds)
    report = evaluate_primer(sequencing, puc19, thresholds=thresholds)
    assert (sequencing.sequence, report["length"].status) == (primer.sequence, "pass")


def test_a_design_keeps_inside_the_bands_its_length_decides_where_a_length_can(puc19) -> None:
    decided = ("length", "gc_percent", "gc_clamp", "tm")
    # The 20-mer nearest the target Tm ends on four Gs and Cs; the 19-mer passes.
    primer = evaluate_primer(design_primer(puc19, 210, Strand.FORWARD))
    assert [primer[name].status for name in decided] == ["pass"] * 4
    # The reverse primer warns at every length here, and does not drag its partner with it.
    forward, reverse = (evaluate_primer(one) for one in design_pair(puc19, 455, 395 + len(puc19)))
    assert [forward[name].status for name in decided] == ["pass"] * 4
    assert reverse.status == "warn"


def test_a_design_weighs_a_check_that_a_length_does_not_decide(puc19) -> None:
    # The 24-mer nearest the target Tm folds into a hairpin at 54 °C, which the four checks a
    # length decides cannot see. A 22-mer at the same position passes every check.
    assert warned(evaluate_primer(at(puc19, 2116, Strand.REVERSE, 24), puc19)) == ["hairpin"]
    primer = design_primer(puc19, 2116, Strand.REVERSE)
    assert primer.sequence == "CCATAACCATGAGTGATAACAC"
    assert warned(evaluate_primer(primer, puc19)) == []


def test_a_passing_primer_is_chosen_wherever_a_length_at_that_position_passes(puc19) -> None:
    for position, strand in (
        (452, Strand.FORWARD),
        (2116, Strand.REVERSE),
        (32, Strand.FORWARD),
        (len(puc19) - 8, Strand.FORWARD),
    ):
        every = [evaluate_primer(at(puc19, position, strand, size), puc19) for size in SIZES]
        chosen = evaluate_primer(design_primer(puc19, position, strand), puc19)
        assert shape(chosen) == min(shape(one) for one in every), position


def test_a_passing_pair_is_chosen_wherever_a_pair_of_lengths_passes(puc19) -> None:
    # A narrow length band, so every pair of lengths can be judged here.
    thresholds = dataclasses.replace(THRESHOLDS, length=Band(18, 24, 18, 24))
    sizes = range(18, 25)
    forward, reverse = design_pair(puc19, 378, 481, thresholds=thresholds)
    chosen = evaluate_pair(forward, reverse, puc19, thresholds=thresholds)
    every = [
        evaluate_pair(
            at(puc19, 378, Strand.FORWARD, one),
            at(puc19, 481, Strand.REVERSE, other),
            puc19,
            thresholds=thresholds,
        )
        for one in sizes
        for other in sizes
    ]
    assert shape(chosen) == min(shape(one) for one in every)


def test_the_same_inputs_choose_the_same_primer(puc19) -> None:
    assert design_primer(puc19, 452, Strand.FORWARD) == design_primer(puc19, 452, Strand.FORWARD)
    assert design_pair(puc19, 378, 481) == design_pair(puc19, 378, 481)


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
    report = evaluate_pair(forward, reverse, puc19)
    assert report["tm_difference"].status == "pass"
    assert report.status == "pass"


def test_a_pair_may_amplify_across_the_origin(puc19) -> None:
    forward, reverse = design_pair(puc19, 452, 396 + len(puc19))
    assert forward.binding_sites[0].start == 452
    assert reverse.binding_sites[0].end % len(puc19) == 396


def test_design_refuses_a_position_a_linear_template_cannot_hold() -> None:
    template = SequenceRecord("ACGT" * 10)
    with pytest.raises(ValueError, match="fit"):
        design_primer(template, 38, Strand.FORWARD)


#: Every length a design considers, which is the band `Thresholds.length` does not fail on.
SIZES = range(15, 36)


def at(template: SequenceRecord, position: int, strand: Strand, size: int) -> Primer:
    """The primer a design considers at this position and length."""
    start = (position if strand is Strand.FORWARD else position - size) % len(template)
    bases = template.extract(Segment(start, start + size))
    sequence = bases if strand is Strand.FORWARD else reverse_complement(bases)
    return Primer("", sequence, binding_sites=(BindingSite(start, start + size, strand),))


def warned(report: PrimerReport | PairReport) -> list[str]:
    """The checks of a report that carried a verdict and did not pass."""
    return [check.name for check in report.checks if check.status not in (None, "pass")]


def shape(report: PrimerReport | PairReport) -> tuple[int, int]:
    """How well a primer or a pair came out: its worst status, then what did not pass."""
    if isinstance(report, PairReport):
        left = warned(report.forward) + warned(report.reverse) + warned(report)
    else:
        left = warned(report)
    return STATUSES.index(report.status), len(left)
