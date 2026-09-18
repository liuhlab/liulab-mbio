"""Each verdict on its own, over records small enough to read.

The plan wires these together; what is pinned here is the branch each one takes, and above all
the two the ticket turns on: a pair nothing sourced can judge carries no verdict rather than a
pass, and a methylation answer that depends on the flanking bases is computed from them.
"""

import dataclasses

import pytest

from liulab_mbio.bench.amounts import dna_amount
from liulab_mbio.cloning.restriction.digest import Diagnostic
from liulab_mbio.cloning.restriction.ligation import Junction
from liulab_mbio.cloning.restriction.verdicts import (
    buffer_check,
    cleanup_check,
    diagnostic_check,
    frame_check,
    methylation_check,
    ratio_check,
    temperature_check,
)
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord

#: Bases spelling no site of any enzyme named below, for padding a record written in code.
FILLER = "ACGT" * 6


def plasmid(sequence: str, name: str) -> SequenceRecord:
    """A circular record written in code."""
    return SequenceRecord(sequence, topology="circular", name=name)


def test_a_pair_supplied_in_one_buffer_passes_and_a_mixed_pair_carries_no_verdict():
    together = buffer_check([get_enzyme("EcoRI"), get_enzyme("BamHI")])
    assert (together.status, together.value) == ("pass", 1)
    assert "rCutSmart Buffer" in together.detail
    mixed = buffer_check([get_enzyme("EcoRI"), get_enzyme("BsmBI")])
    assert mixed.status is None
    assert mixed.value == 2
    assert "no verdict rather than a pass" in mixed.detail
    assert "Double Digest Finder" in mixed.detail


def test_two_enzymes_wanting_different_temperatures_are_reported():
    agreed = temperature_check([get_enzyme("EcoRI"), get_enzyme("BamHI")])
    assert (agreed.status, agreed.detail) == ("pass", "every one incubates at 37 °C")
    apart = temperature_check([get_enzyme("EcoRI"), get_enzyme("BsmBI")])
    assert (apart.status, apart.value) == ("warn", 2)
    assert "EcoRI at 37 °C and BsmBI at 55 °C" in apart.detail
    assert "two digests" in apart.detail


def test_an_enzyme_heat_cannot_stop_says_so_and_asks_for_the_clean_up():
    stopped = cleanup_check([get_enzyme("EcoRI")])
    assert (stopped.status, stopped.value) == ("pass", 0)
    assert "65 °C for 20 minutes" in stopped.detail
    unstopped = cleanup_check([get_enzyme("EcoRI"), get_enzyme("BamHI")])
    assert unstopped.value == 1
    assert "no heat inactivation for BamHI-HF" in unstopped.detail
    assert "gel purification" in unstopped.detail


def test_methylation_blocked_by_overlapping_is_read_off_the_flanks_and_not_the_enzyme_name():
    # XbaI reads TCTAGA and its supplier calls it blocked by an OVERLAPPING Dam site, which is
    # GATC. Only the first record spells one across it.
    over = plasmid(FILLER + "GATCTAGA" + FILLER, "pDam")
    clear = plasmid(FILLER + "ATCTAGAA" + FILLER, "pClear")
    xbai = [get_enzyme("XbaI")]
    fired = methylation_check(xbai, [over])
    assert (fired.status, fired.value) == ("warn", 1)
    assert "pDam carries a Dam site" in fired.detail
    assert "C2925" in fired.detail
    assert methylation_check(xbai, [clear]).status == "pass"
    assert methylation_check([get_enzyme("EcoRI")], [over]).status == "pass"


def coding(name: str, *spans: tuple[int, int]) -> Feature:
    """One coding sequence over those spans."""
    return Feature(name, "CDS", tuple(Segment(start, end) for start, end in spans))


def test_a_junction_inside_a_coding_sequence_reports_the_frame_the_added_bases_leave():
    junctions = (Junction(10, "AATT", "backbone", "insert"),)
    record = SequenceRecord("A" * 60, topology="circular", name="p")
    kept = frame_check(
        dataclasses.replace(record, features=(coding("tag", (0, 10), (19, 30)),)), junctions
    )
    assert (kept.status, kept.value) == ("pass", 1)
    assert "9 bases" in kept.detail
    shifted = frame_check(
        dataclasses.replace(record, features=(coding("tag", (0, 10), (18, 30)),)), junctions
    )
    assert shifted.status == "warn"
    assert "2 past a whole number of codons" in shifted.detail
    assert frame_check(record, junctions).value == 0


def test_the_ligation_ratio_is_judged_in_picomoles_and_shown_in_nanograms_too():
    backbone = dna_amount("backbone", 2665, pmol=0.02)
    inside = ratio_check(backbone, dna_amount("insert", 723, pmol=0.06))
    assert (inside.status, inside.value) == ("pass", 3.0)
    assert "0.06 pmol (26.72 ng)" in inside.detail
    outside = ratio_check(backbone, dna_amount("insert", 723, pmol=0.4))
    assert outside.status == "warn"


def test_a_diagnostic_digest_that_gives_the_vector_s_own_bands_fails():
    pair = (get_enzyme("EcoRI"), get_enzyme("BamHI"))
    told = diagnostic_check(Diagnostic(pair, (2665, 723), (2665, 21), ("clone", "pUC19")))
    assert told.status == "pass"
    assert "2665 bp and 723 bp" in told.detail
    same = diagnostic_check(Diagnostic(pair, (2665, 21), (2665, 21), ("clone", "pUC19")))
    assert same.status == "fail"
    assert "cannot tell one from the other" in same.detail


@pytest.mark.parametrize("name", ["EcoRI", "XbaI"])
def test_every_shipped_multiple_cloning_site_enzyme_is_judged_on_its_own(name: str):
    # One enzyme is the shape #132 makes possible, so no check may assume a pair.
    assert buffer_check([get_enzyme(name)]).status == "pass"
    assert temperature_check([get_enzyme(name)]).status == "pass"
    assert cleanup_check([get_enzyme(name)]).status == "pass"
