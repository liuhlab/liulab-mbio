import dataclasses
from collections import Counter
from collections.abc import Container, Iterator
from itertools import islice

import pytest

from liulab_mbio.checks import STATUSES
from liulab_mbio.primers import (
    TARGET_TM,
    THRESHOLDS,
    THRESHOLDS_FOR,
    Band,
    PairReport,
    Placement,
    PrimerReport,
    design_pair,
    design_primer,
    evaluate_pair,
    evaluate_primer,
    melting_temperature,
    ranked_pairs,
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


def test_an_anchored_placement_chooses_what_a_bare_position_does(puc19) -> None:
    anchored = Placement(five_prime=Segment(452, 453))
    primer = design_primer(puc19, 452, Strand.FORWARD, placement=anchored)
    assert primer == design_primer(puc19, 452, Strand.FORWARD)
    assert primer.binding_sites[0].start == 452


def test_a_placement_near_a_target_chooses_the_best_primer_it_allows(puc19) -> None:
    # A sequencing primer's 3' end, 100 to 130 bases outside a junction.
    junction = 480
    near = Placement(three_prime=Segment(junction - 130, junction - 99))
    chosen = design_primer(puc19, junction - 130, Strand.FORWARD, placement=near, thresholds=NARROW)
    assert 100 <= junction - chosen.binding_sites[0].end <= 130
    assert best(chosen, puc19, near, Strand.FORWARD)


def test_a_placement_may_cross_the_origin(puc19) -> None:
    across = Placement(five_prime=Segment(len(puc19) - 6, len(puc19) + 7))
    chosen = design_primer(puc19, len(puc19), Strand.FORWARD, placement=across, thresholds=NARROW)
    every = list(allowed(puc19, across, Strand.FORWARD))
    assert any(primer.binding_sites[0].end > len(puc19) for primer in every)
    assert best(chosen, puc19, across, Strand.FORWARD)


def test_the_option_nearest_the_position_asked_for_wins(puc19) -> None:
    # Three copies of one stretch of pUC19, and a placement two copies wide, so every candidate
    # has a twin 50 bases away spelling the same bases and scoring the same on every check.
    # Only nearness to the position asked for tells the two apart.
    unit = puc19.extract(Segment(452, 502))
    template = SequenceRecord(unit * 3, topology="circular")
    thresholds = dataclasses.replace(THRESHOLDS, length=Band(20, 20, 20, 20))
    across = Placement(five_prime=Segment(0, 101))
    chosen = design_primer(template, 40, Strand.FORWARD, placement=across, thresholds=thresholds)
    twins = [
        start
        for start in range(101)
        if template.extract(Segment(start, start + 20)) == chosen.sequence
    ]
    assert len(twins) >= 2
    assert chosen.binding_sites[0].start == min(twins, key=lambda start: abs(start - 40))
    # Equidistant from the two, and the lower position settles it.
    middle = (twins[0] + twins[1]) // 2
    tied = design_primer(template, middle, Strand.FORWARD, placement=across, thresholds=thresholds)
    assert tied.binding_sites[0].start == twins[0]


def test_a_placement_holding_no_annealing_region_is_refused(puc19) -> None:
    # Its two ends lie 5 bases apart, shorter than any length the band allows.
    placement = Placement(five_prime=Segment(452, 453), three_prime=Segment(457, 458))
    with pytest.raises(ValueError, match="placement"):
        design_primer(puc19, 452, Strand.FORWARD, placement=placement)


def test_the_same_inputs_choose_the_same_primer_inside_a_placement(puc19) -> None:
    placement = Placement(five_prime=Segment(448, 457))
    first = design_primer(puc19, 452, Strand.FORWARD, placement=placement, thresholds=NARROW)
    second = design_primer(puc19, 452, Strand.FORWARD, placement=placement, thresholds=NARROW)
    assert first == second


def test_pairs_come_back_in_the_order_judging_every_pair_gives(puc19, every_pair) -> None:
    assert {pair.status for pair in every_pair} == {"pass", "warn", "fail"}
    ranked = list(pairs_in_regions(puc19))
    assert Counter(ranked) == Counter(every_pair)
    assert [rank(pair) for pair in ranked] == sorted(rank(pair) for pair in every_pair)


def test_the_first_pair_in_rank_order_is_the_pair_a_design_chooses(puc19) -> None:
    first = next(pairs_in_regions(puc19))
    assert (first.forward.primer, first.reverse.primer) == pair_in_regions(puc19)


def test_an_excluded_binding_site_is_never_part_of_a_pair(puc19, every_pair) -> None:
    excluded = {primer.binding_sites[0] for primer in pair_in_regions(puc19)}
    left = [pair for pair in every_pair if not excluded & sites(pair)]
    assert Counter(pairs_in_regions(puc19, excluded)) == Counter(left)
    chosen = evaluate_pair(*pair_in_regions(puc19, excluded), puc19, thresholds=NARROW)
    assert rank(chosen) == min(rank(pair) for pair in left)


def test_a_site_excluded_while_pairs_are_taken_stays_out_of_every_later_pair(puc19) -> None:
    order = list(pairs_in_regions(puc19))
    excluded: set[BindingSite] = set()
    pairs = pairs_in_regions(puc19, excluded)
    # Past the passing pairs, where pairs already judged wait for their turn.
    taken = list(islice(pairs, 70))
    excluded.add(order[70].forward.primer.binding_sites[0])
    kept = [pair for pair in taken if not excluded & sites(pair)]
    assert kept + list(pairs) == list(pairs_in_regions(puc19, excluded))


def test_excluding_every_binding_site_at_one_end_is_refused(puc19, every_pair) -> None:
    forwards = {pair.forward.primer.binding_sites[0] for pair in every_pair}
    reverses = {pair.reverse.primer.binding_sites[0] for pair in every_pair}
    with pytest.raises(ValueError, match="every forward binding site"):
        pairs_in_regions(puc19, forwards)
    with pytest.raises(ValueError, match="every reverse binding site"):
        pair_in_regions(puc19, reverses)


#: Every length a design considers, which is the band `Thresholds.length` does not fail on.
SIZES = range(15, 36)

#: A narrow band, so a placement's options stay few enough to judge every one of them.
NARROW = dataclasses.replace(THRESHOLDS, length=Band(18, 24, 18, 24))

#: A pair's two regions on pUC19, where pairs pass, warn and fail.
START, END = 1160, 1310
FORWARD_REGION = Placement(five_prime=Segment(1160, 1163), three_prime=Segment(1180, 1184))
REVERSE_REGION = Placement(five_prime=Segment(1310, 1313), three_prime=Segment(1286, 1290))


def pairs_in_regions(
    template: SequenceRecord, exclude: Container[BindingSite] = ()
) -> Iterator[PairReport]:
    """Every pair the two regions allow, in rank order."""
    return ranked_pairs(
        template,
        START,
        END,
        forward_placement=FORWARD_REGION,
        reverse_placement=REVERSE_REGION,
        thresholds=NARROW,
        exclude=exclude,
    )


def pair_in_regions(
    template: SequenceRecord, exclude: Container[BindingSite] = ()
) -> tuple[Primer, Primer]:
    """The pair a design chooses in the two regions."""
    return design_pair(
        template,
        START,
        END,
        forward_placement=FORWARD_REGION,
        reverse_placement=REVERSE_REGION,
        thresholds=NARROW,
        exclude=exclude,
    )


def sites(report: PairReport) -> set[BindingSite]:
    """The binding sites of a pair's two primers."""
    return {report.forward.primer.binding_sites[0], report.reverse.primer.binding_sites[0]}


@pytest.fixture(scope="module")
def every_pair(puc19) -> list[PairReport]:
    """Every pair the two regions allow, each judged by `evaluate_pair`."""
    reverses = list(allowed(puc19, REVERSE_REGION, Strand.REVERSE))
    return [
        evaluate_pair(one, other, puc19, thresholds=NARROW)
        for one in allowed(puc19, FORWARD_REGION, Strand.FORWARD)
        for other in reverses
    ]


def rank(report: PairReport) -> tuple[float, ...]:
    """How a design ranks a pair in the regions, as far as its docstring says."""
    primers = (report.forward, report.reverse)
    fives = (
        primers[0].primer.binding_sites[0].start - START,
        primers[1].primer.binding_sites[0].end - END,
    )
    return (
        *shape(report),
        sum(one["length"].status != "pass" for one in primers),
        sum(abs(one["tm"].value - TARGET_TM) for one in primers),
        sum(abs(five) for five in fives),
    )


def allowed(
    template: SequenceRecord, placement: Placement, strand: Strand, sizes: range = range(18, 25)
):
    """Every primer a placement allows, which is what a design is choosing between."""
    for start in range(len(template)):
        for size in sizes:
            site = BindingSite(start, start + size, strand)
            if template.topology != "circular" and site.end > len(template):
                continue
            if placement.allows(site, template):
                bases = template.extract(Segment(site.start, site.end))
                sequence = bases if strand is Strand.FORWARD else reverse_complement(bases)
                yield Primer("", sequence, binding_sites=(site,))


def best(chosen: Primer, template: SequenceRecord, placement: Placement, strand: Strand) -> bool:
    """Whether no primer the placement allows came out better than the one chosen."""
    every = [
        evaluate_primer(primer, template, thresholds=NARROW)
        for primer in allowed(template, placement, strand)
    ]
    return shape(evaluate_primer(chosen, template, thresholds=NARROW)) == min(
        shape(one) for one in every
    )


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
