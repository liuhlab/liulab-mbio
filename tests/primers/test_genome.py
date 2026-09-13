"""The genome check runs the real `ipcr` on FASTA files each test writes, a few kilobases at most."""

import dataclasses
import random
from itertools import islice
from pathlib import Path

import pytest

from liulab_mbio.primers import THRESHOLDS, PairReport, Placement, amplicon_sizes, ranked_pairs
from liulab_mbio.primers.genome import (
    Locus,
    design_pair_on_genome,
    evaluate_on_genome,
    evaluate_pair_on_genome,
)
from liulab_mbio.sequence import BindingSite, Primer, Segment, SequenceRecord, reverse_complement

from ..fasta import write_fasta
from .sequences import M13_FWD, M13_REV, PUC_FWD, PUC_REV

#: Pairs whose primers differ from each other in eight bases or more, so none primes another's site.
ONE = ("AGATCATACACGAAGGAACC", "ACATCCCTTCATACGCTTTG")
MASKED = ("AACGTCTCTGGACTAACGAC", "TCCCTTTGCTTGTGATATCC")
ALONE = ("TGTAATCTCAAAGCTACCGC", "GACACAGAAGATGATCTGTC")
NEAR = ("CAGATTGATGATCATCGTCG", "CGACAAATACCGTTGGCGAC")


def _pair(sequences: tuple[str, str]) -> tuple[Primer, Primer]:
    return Primer("fwd", sequences[0]), Primer("rev", sequences[1])


def _mismatched(sequence: str, *positions: int) -> str:
    bases = list(sequence)
    for position in positions:
        bases[position] = reverse_complement(bases[position])
    return "".join(bases)


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """Six sequences, each planting what one test looks for between random bases."""
    rng = random.Random(2686)

    def gap(size: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(size))

    records = {
        # ONE's amplicon at 100-460, and a second forward site at 220 inside it.
        "chrI__ce11": gap(100)
        + ONE[0]
        + gap(100)
        + ONE[0]
        + gap(200)
        + reverse_complement(ONE[1])
        + gap(100),
        "chrII": (
            gap(100) + MASKED[0] + gap(250) + reverse_complement(MASKED[1]) + gap(100)
        ).lower(),
        "chrIII": gap(100) + ALONE[0] + gap(300) + reverse_complement(ALONE[0]) + gap(100),
        "chrIV": gap(100) + ALONE[1] + gap(200) + reverse_complement(ALONE[1]) + gap(100),
        # Two mismatches at the 5' end, and three mid-primer, which melt too far below a match.
        "chrV": gap(100) + _mismatched(NEAR[0], 0, 1) + gap(200) + reverse_complement(NEAR[1]),
        "chrVI": gap(100)
        + _mismatched(NEAR[0], 8, 10, 12)
        + gap(200)
        + reverse_complement(NEAR[1]),
    }
    path = tmp_path_factory.mktemp("genome") / "planted.fa"
    path.write_text("".join(f">{name} planted\n{bases}\n" for name, bases in records.items()))
    return path


@pytest.fixture(scope="module")
def reports(genome):
    """One run over every pair, ONE twice: named at its amplicon, and at a place it makes none."""
    return evaluate_on_genome(
        [_pair(ONE), _pair(MASKED), _pair(ALONE), _pair(NEAR), _pair(ONE)],
        genome,
        "planted",
        intended=[Locus("chrI__ce11", 100, 460), None, None, None, Locus("chrI__ce11", 0, 460)],
    )


def _found(report):
    return [
        (amplicon.sequence_name, amplicon.start, amplicon.end, amplicon.length, amplicon.made_by)
        for amplicon in report.amplicons
    ]


def test_on_a_plasmid_the_amplicons_are_the_template_checks(puc19, tmp_path) -> None:
    fasta = tmp_path / "pUC19.fa"
    fasta.write_text(f">pUC19\n{puc19.sequence}\n")
    across = puc19.sequence[-10:] + puc19.sequence[:10]
    pairs = [
        _pair(sequences)
        for sequences in (
            (M13_FWD, M13_REV),
            (PUC_FWD, PUC_REV),
            (across, M13_REV),
            (M13_REV, M13_FWD),
        )
    ]
    reports = evaluate_on_genome(pairs, fasta, "pUC19", circular=True)
    for (forward, reverse), report in zip(pairs, reports, strict=True):
        assert [one.length for one in report.amplicons] == list(
            amplicon_sizes(forward, reverse, puc19)
        )
    # Where the template check places each pair's binding sites, across the origin for the third.
    assert [
        (one.start, one.end, one.made_by) for report in reports for one in report.amplicons
    ] == [
        (378, 481, "pair"),
        (363, 500, "pair"),
        (2676, 481 + len(puc19), "pair"),
        (378, 481, "pair"),
    ]


def test_a_second_binding_site_gives_an_off_target_amplicon(reports) -> None:
    report = reports[0]
    assert report.assembly == "planted"
    assert _found(report) == [
        ("chrI__ce11", 100, 460, 360, "pair"),
        ("chrI__ce11", 220, 460, 240, "pair"),
    ]
    assert [one.start for one in report.off_target] == [220]
    off_target = report["off_target_amplicons"]
    assert (off_target.status, off_target.value) == ("warn", 1)
    assert "chrI__ce11:220-460" in off_target.detail
    assert "240 bp" in off_target.detail


def test_an_intended_amplicon_is_marked_and_one_that_is_absent_warns(reports) -> None:
    present, absent = reports[0], reports[4]
    assert [one.intended for one in present.amplicons] == [True, False]
    assert present["intended_amplicon"].status == "pass"
    assert not any(one.intended for one in absent.amplicons)
    assert absent["intended_amplicon"].status == "warn"
    assert "chrI__ce11:0-460" in absent["intended_amplicon"].detail


def test_a_primer_facing_itself_makes_an_amplicon_alone(reports) -> None:
    report = reports[2]
    assert _found(report) == [
        ("chrIII", 100, 440, 340, "forward"),
        ("chrIV", 100, 340, 240, "reverse"),
    ]
    assert report.status == "warn"


def test_a_soft_masked_sequence_still_yields_its_amplicon(reports) -> None:
    assert _found(reports[1]) == [("chrII", 100, 390, 290, "pair")]


def test_near_matches_are_searched_only_up_to_the_size_limit(genome, reports) -> None:
    # The default limit lies above this genome's size.
    near = reports[3]
    assert near.near_matches_checked
    assert _found(near) == [("chrV", 100, 340, 240, "pair")]
    assert near.amplicons[0].mismatches == (2, 0)
    small = dataclasses.replace(THRESHOLDS, genome_size_limit=genome.stat().st_size - 1)
    perfect = evaluate_pair_on_genome(*_pair(NEAR), genome, "planted", thresholds=small)
    assert not perfect.near_matches_checked
    assert perfect.amplicons == ()
    assert "near-match sites were not checked" in perfect["off_target_amplicons"].detail


def test_a_missing_tool_or_a_run_past_its_time_limit_raises_a_clear_error(
    genome, tmp_path, monkeypatch
) -> None:
    with pytest.raises(TimeoutError, match="killed"):
        evaluate_pair_on_genome(*_pair(ONE), genome, "planted", timeout=0)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(RuntimeError, match="pixi"):
        evaluate_pair_on_genome(*_pair(ONE), genome, "planted")


#: A region to amplify, in a record named as a chimera genome spells one, and its flank width.
REGION = Locus("I__ce11", 200, 300)
FLANK = 25
#: Where a test plants a second place the design's primers could prime.
ELSEWHERE = "chrII__ecHT115"
#: Bases padding what is planted there, so an off-target amplicon does not start at zero.
PAD = "".join(random.Random(58).choice("ACGT") for _ in range(60))


@pytest.fixture(scope="module")
def sequence() -> str:
    """The record the region lies in: random bases, so nothing primes twice by accident."""
    rng = random.Random(1901)
    return "".join(rng.choice("ACGT") for _ in range(600))


@pytest.fixture(scope="module")
def template(sequence) -> str:
    """What a design works on: the region and a flank each side."""
    return sequence[REGION.start - FLANK : REGION.end + FLANK]


@pytest.fixture(scope="module")
def order(template) -> list[PairReport]:
    """The pairs a design chooses between, best first, as its docstring places them."""
    record = SequenceRecord(template)
    left, right = FLANK, len(template) - FLANK
    return list(
        islice(
            ranked_pairs(
                record,
                left,
                right,
                forward_placement=Placement(
                    five_prime=Segment(0, left), three_prime=Segment(0, left + 1)
                ),
                reverse_placement=Placement(
                    five_prime=Segment(right + 1, len(record) + 1),
                    three_prime=Segment(right, len(record)),
                ),
            ),
            20,
        )
    )


@pytest.fixture(scope="module")
def together(sequence, template, order, tmp_path_factory) -> Path:
    """A genome where the best-ranked pair's two primers amplify somewhere else as well."""
    pair = order[0]
    planted = [_copied(template, site) for site in _sites(pair)]
    return _written(tmp_path_factory, "together.fa", sequence, PAD.join(["", *planted, ""]))


@pytest.fixture(scope="module")
def alone(sequence, template, order, tmp_path_factory) -> Path:
    """A genome where the best-ranked pair's reverse primer amplifies on its own."""
    copy = reverse_complement(_copied(template, _sites(order[0])[1]))
    planted = PAD.join(["", copy, reverse_complement(copy), ""])
    return _written(tmp_path_factory, "alone.fa", sequence, planted)


@pytest.fixture(scope="module")
def everywhere(sequence, template, tmp_path_factory) -> Path:
    """A genome carrying the whole template twice, so no pair at all is specific."""
    return _written(tmp_path_factory, "everywhere.fa", sequence, PAD + template + PAD)


def test_the_best_ranked_pair_with_no_off_target_amplicon_is_the_one_chosen(
    together, order
) -> None:
    design = design_pair_on_genome(together, "planted", REGION, flank=FLANK)
    assert (design.specific, design.rounds, design.exhausted) == (True, 1, False)
    assert design.template == Locus(REGION.sequence_name, REGION.start - FLANK, REGION.end + FLANK)
    assert design.genome.off_target == ()
    assert [one.sequence_name for one in design.genome.amplicons] == [REGION.sequence_name]
    # Every pair ranked above the one chosen makes an off-target amplicon, and it makes none.
    chosen = [_primers(pair) for pair in order].index((design.forward, design.reverse))
    above = order[:chosen]
    assert above
    reports = evaluate_on_genome(
        [_primers(pair) for pair in above],
        together,
        "planted",
        intended=[_amplicon(pair) for pair in above],
    )
    assert all(report.off_target for report in reports)


def test_a_primer_that_amplifies_on_its_own_is_kept_out_of_every_later_pair(alone, order) -> None:
    # One pair a search, so a pair carrying that primer is checked again unless it is kept out.
    thresholds = dataclasses.replace(THRESHOLDS, genome_pairs_per_search=1, genome_rounds=2)
    ruled_out = _sites(order[0])[1]
    assert _sites(order[1])[1].start == ruled_out.start
    design = design_pair_on_genome(alone, "planted", REGION, flank=FLANK, thresholds=thresholds)
    assert (design.specific, design.rounds) == (True, 2)
    assert design.reverse.binding_sites[0].start != ruled_out.start


def test_a_pair_that_amplifies_only_together_leaves_each_of_its_primers_free(
    together, order
) -> None:
    thresholds = dataclasses.replace(THRESHOLDS, genome_pairs_per_search=1)
    design = design_pair_on_genome(together, "planted", REGION, flank=FLANK, thresholds=thresholds)
    assert (design.specific, design.rounds) == (True, 2)
    assert {design.forward, design.reverse} & set(_primers(order[0]))


def test_when_no_pair_is_specific_the_best_ranked_pair_checked_comes_back(
    everywhere, order
) -> None:
    design = design_pair_on_genome(everywhere, "planted", REGION, flank=FLANK)
    assert (design.specific, design.rounds, design.exhausted) == (False, 3, False)
    assert (design.forward, design.reverse) == _primers(order[0])
    assert [one.sequence_name for one in design.genome.off_target] == [ELSEWHERE]
    assert design.genome.off_target[0].length == design.pair.amplicon_length
    assert design.genome["off_target_amplicons"].status == "warn"


def test_a_search_that_runs_out_of_pairs_stops_and_says_so(everywhere) -> None:
    # Flanks this narrow allow nine pairs in all, fewer than one search checks.
    design = design_pair_on_genome(everywhere, "planted", REGION, flank=16)
    assert (design.specific, design.rounds, design.exhausted) == (False, 1, True)


def test_a_fasta_with_no_index_is_refused(sequence, tmp_path) -> None:
    path = tmp_path / "bare.fa"
    path.write_text(f">{REGION.sequence_name}\n{sequence}\n")
    with pytest.raises(FileNotFoundError, match="genome assembly register"):
        design_pair_on_genome(path, "planted", REGION, flank=FLANK)


def test_a_region_with_no_flank_beside_it_is_refused(together) -> None:
    with pytest.raises(ValueError, match="flank"):
        design_pair_on_genome(together, "planted", Locus(REGION.sequence_name, 0, 100), flank=FLANK)


def _written(factory, name: str, sequence: str, planted: str) -> Path:
    path = factory.mktemp("region") / name
    write_fasta(path, {REGION.sequence_name: sequence, ELSEWHERE: planted})
    return path


def _copied(template: str, site: BindingSite) -> str:
    """A binding site's bases, between neighbours complementing the ones it has on the template.

    A primer with another 3' end mismatches the copy where it reaches a neighbour, so only the
    primers of the site planted prime here.
    """
    flip = str.maketrans("ACGT", "TGCA")
    return (
        template[max(site.start - 4, 0) : site.start].translate(flip)
        + template[site.start : site.end]
        + template[site.end : site.end + 4].translate(flip)
    )


def _sites(pair: PairReport) -> tuple[BindingSite, BindingSite]:
    """Where a pair's two primers anneal on the template."""
    return pair.forward.primer.binding_sites[0], pair.reverse.primer.binding_sites[0]


def _primers(pair: PairReport) -> tuple[Primer, Primer]:
    return pair.forward.primer, pair.reverse.primer


def _amplicon(pair: PairReport) -> Locus:
    """Where a pair's amplicon lies on the genome, counted from the template's own start."""
    forward, reverse = _sites(pair)
    start = REGION.start - FLANK
    return Locus(REGION.sequence_name, start + forward.start, start + reverse.end)
