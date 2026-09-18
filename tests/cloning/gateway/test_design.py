"""The attB tail an insert is amplified with, and the PCR that puts it there."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from liulab_mbio.cloning.gateway.design import (
    C_TERMINAL_FRAME,
    FUSIONS,
    N_TERMINAL_FRAME,
    SPACER,
    TAIL_RUN,
    amplify_attb,
    attb_tails,
    attb_thresholds,
)
from liulab_mbio.primers.polymerase import Q5, melting_temperature
from liulab_mbio.primers.thresholds import THRESHOLDS
from liulab_mbio.sequence import SequenceRecord, reverse_complement

if TYPE_CHECKING:
    from liulab_mbio.cloning.gateway import Amplicon

#: The whole tail MAN0000470 page 13 asks for, as the note's own table writes it (note §7).
FORWARD_TAIL = "GGGGACAAGTTTGTACAAAAAAGCAGGCT"
REVERSE_TAIL = "GGGGACCACTTTGTACAAGAAAGCTGGGT"


def test_the_tail_is_the_four_g_residues_and_the_whole_att_site() -> None:
    forward, reverse = attb_tails()

    assert (forward, reverse) == (FORWARD_TAIL, REVERSE_TAIL)
    assert [len(forward), len(reverse)] == [29, 29]
    assert forward.startswith(SPACER)
    assert reverse.startswith(SPACER)


def test_a_fusion_adds_the_frame_bases_the_note_counts_and_nothing_else() -> None:
    lengths = {fusion: [len(tail) for tail in attb_tails(fusion)] for fusion in FUSIONS}

    assert lengths == {
        "none": [29, 29],
        "N-terminal": [31, 29],
        "C-terminal": [29, 30],
        "both": [31, 30],
    }
    assert attb_tails("both") == (FORWARD_TAIL + N_TERMINAL_FRAME, REVERSE_TAIL + C_TERMINAL_FRAME)


def test_the_n_terminal_frame_bases_make_no_stop_codon_with_the_att_site() -> None:
    forward, _ = attb_tails("N-terminal")

    assert forward[-3:] not in {"TAA", "TAG", "TGA"}
    assert N_TERMINAL_FRAME not in {"AA", "AG", "GA"}
    assert (len(N_TERMINAL_FRAME), len(C_TERMINAL_FRAME)) == (2, 1)


def test_the_run_band_is_widened_only_to_the_run_the_tail_itself_spells() -> None:
    widened = attb_thresholds()

    assert "A" * TAIL_RUN in attb_tails()[0]
    assert widened.mononucleotide_run.high == THRESHOLDS.mononucleotide_run.high
    assert widened.mononucleotide_run.grade(TAIL_RUN) == "warn"
    assert widened.mononucleotide_run.grade(TAIL_RUN + 1) == "fail"


def test_the_amplicon_is_the_forward_tail_the_insert_and_the_reverse_tail_turned_round(
    amplicon: Amplicon, gfp: SequenceRecord
) -> None:
    forward, reverse = amplicon.tails

    assert amplicon.record.sequence == forward + gfp.sequence + reverse_complement(reverse)
    assert amplicon.length == len(gfp) + len(forward) + len(reverse)
    assert amplicon.primers[0].sequence.startswith(forward)
    assert amplicon.primers[1].sequence.startswith(reverse)


def test_each_primer_reads_its_tm_off_the_annealing_region_alone(
    amplicon: Amplicon, gfp: SequenceRecord
) -> None:
    for report, tail in zip(amplicon.reports, amplicon.tails, strict=True):
        annealing = report.primer.sequence[len(tail) :]

        assert gfp.sequence.startswith(annealing) or gfp.sequence.endswith(
            reverse_complement(annealing)
        )
        assert report["tm"].value == pytest.approx(melting_temperature(annealing, Q5))
        assert report["tm"].value < report["tm_full"].value
        assert report.status != "fail"


def test_an_insert_no_annealing_region_fits_is_refused_naming_the_record() -> None:
    with pytest.raises(ValueError, match=r"tiny takes no attB primer.*12 bases long"):
        amplify_attb(SequenceRecord("ATGCATGCATGC", name="tiny"))
