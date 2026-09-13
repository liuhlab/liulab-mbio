"""The genome check runs the real `ipcr` on FASTA files each test writes, a few kilobases at most."""

import dataclasses
import random

import pytest

from liulab_mbio.primers import THRESHOLDS, amplicon_sizes
from liulab_mbio.primers.genome import Locus, evaluate_on_genome, evaluate_pair_on_genome
from liulab_mbio.sequence import Primer, reverse_complement

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
