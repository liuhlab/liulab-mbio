import pytest

from mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    across_the_origin,
    counted_round,
    position_text,
    reverse_complement,
)


def _feature(*segments: tuple[int, int], strand: Strand = Strand.FORWARD) -> Feature:
    return Feature(
        name="f",
        type="misc_feature",
        segments=tuple(Segment(start, end) for start, end in segments),
        strand=strand,
    )


def test_a_record_holds_its_sequence_in_upper_case() -> None:
    record = SequenceRecord("acgTN", topology="circular")
    assert record.sequence == "ACGTN"
    assert record.topology == "circular"
    assert len(record) == 5


def test_a_record_refuses_letters_outside_the_iupac_dna_alphabet() -> None:
    with pytest.raises(ValueError, match="X"):
        SequenceRecord("ACGX")


def test_a_record_refuses_an_unknown_topology() -> None:
    with pytest.raises(ValueError, match="topology"):
        SequenceRecord("ACGT", topology="cirular")  # pyright: ignore[reportArgumentType]


def test_a_segment_is_half_open_and_must_hold_at_least_one_base() -> None:
    assert Segment(0, 1).end == 1
    with pytest.raises(ValueError, match="start"):
        Segment(3, 3)
    with pytest.raises(ValueError, match="start"):
        Segment(-1, 2)


def test_a_feature_needs_at_least_one_segment() -> None:
    with pytest.raises(ValueError, match="segment"):
        _feature()


def test_a_feature_may_end_at_the_last_base_of_a_linear_record() -> None:
    record = SequenceRecord("ACGTACGT", features=(_feature((0, 8)),))
    assert record.features[0].segments == (Segment(0, 8),)


def test_a_feature_may_not_run_past_the_end_of_a_linear_record() -> None:
    with pytest.raises(ValueError, match="linear"):
        SequenceRecord("ACGTACGT", features=(_feature((6, 9)),))


def test_a_segment_across_the_origin_ends_past_the_length_of_a_circular_record() -> None:
    record = SequenceRecord("ACGTACGT", topology="circular", features=(_feature((6, 10)),))
    assert len(record.features) == 1


def test_a_segment_may_not_start_past_the_origin_or_wrap_more_than_once() -> None:
    with pytest.raises(ValueError, match="circular"):
        SequenceRecord("ACGTACGT", topology="circular", features=(_feature((8, 10)),))
    with pytest.raises(ValueError, match="circular"):
        SequenceRecord("ACGTACGT", topology="circular", features=(_feature((2, 11)),))


def test_a_feature_reads_the_same_bases_with_its_later_segments_past_the_length(
    across_origin: SequenceRecord,
) -> None:
    across, site = across_origin.features
    top = across_origin.sequence
    assert across.segments == (Segment(90, 98), Segment(101, 108))
    assert across_origin.extract(across) == top[90:98] + top[1:8]
    assert across_origin.extract(site) == reverse_complement(top[96:] + top[:4])


def test_a_later_segment_across_the_origin_starts_past_the_length_within_one_turn() -> None:
    def circular(*segments: tuple[int, int]) -> SequenceRecord:
        return SequenceRecord("ACGTACGT", topology="circular", features=(_feature(*segments),))

    assert circular((6, 8), (9, 10)).extract(circular((6, 8), (9, 10)).features[0]) == "GTC"
    with pytest.raises(ValueError, match="circular"):
        circular((2, 6), (9, 11))
    with pytest.raises(ValueError, match="circular"):
        circular((6, 8), (1, 2))


def test_sorting_a_feature_across_the_origin_leaves_its_segments_as_they_are(
    across_origin: SequenceRecord,
) -> None:
    for feature in across_origin.features:
        assert tuple(sorted(feature.segments, key=lambda one: one.start)) == feature.segments


def test_a_feature_at_both_ends_of_a_linear_record_lists_its_higher_start_first() -> None:
    record = SequenceRecord("ACGTACGT", features=(_feature((6, 8), (0, 2)),))
    assert record.features[0].segments == (Segment(6, 8), Segment(0, 2))
    assert record.extract(record.features[0]) == "GTAC"


def test_reverse_complement_pairs_iupac_codes_in_either_case() -> None:
    assert reverse_complement("AACRYK") == "MRYGTT"
    assert reverse_complement("aacg") == "cgtt"


LINEAR = SequenceRecord("AACCGGTTAC")
CIRCULAR = SequenceRecord("AACCGGTTAC", topology="circular")


def test_extract_reads_a_forward_feature_off_the_top_strand() -> None:
    assert LINEAR.extract(_feature((2, 6))) == "CCGG"


def test_extract_reads_a_reverse_feature_as_its_reverse_complement() -> None:
    assert LINEAR.extract(_feature((0, 4), strand=Strand.REVERSE)) == "GGTT"


def test_extract_joins_segments_in_top_strand_order_before_complementing() -> None:
    assert LINEAR.extract(_feature((0, 3), (7, 10))) == "AACTAC"
    assert LINEAR.extract(_feature((0, 3), (7, 10), strand=Strand.REVERSE)) == "GTAGTT"


def test_extract_reads_a_segment_across_the_origin() -> None:
    assert CIRCULAR.extract(Segment(8, 12)) == "ACAA"
    assert CIRCULAR.extract(_feature((8, 12), strand=Strand.REVERSE)) == "TTGT"
    assert CIRCULAR.extract(Segment(3, 13)) == "CGGTTACAAC"


def test_extract_refuses_a_span_that_does_not_fit_the_record() -> None:
    with pytest.raises(ValueError, match="circular"):
        CIRCULAR.extract(Segment(10, 14))
    with pytest.raises(ValueError, match="circular"):
        CIRCULAR.extract(_feature((2, 13)))
    with pytest.raises(ValueError, match="linear"):
        LINEAR.extract(Segment(8, 12))


def test_a_primer_holds_its_sequence_in_upper_case_and_refuses_other_letters() -> None:
    assert Primer("M13 fwd", "gtaaaacgacggccagt").sequence == "GTAAAACGACGGCCAGT"
    with pytest.raises(ValueError, match="U"):
        Primer("bad", "ACGU")


def test_a_binding_site_lies_on_one_strand() -> None:
    with pytest.raises(ValueError, match="strand"):
        BindingSite(0, 4, Strand.BOTH)


def test_a_binding_site_across_the_origin_fits_a_circular_record_only() -> None:
    primer = Primer("p", "TTGT", binding_sites=(BindingSite(8, 12, Strand.REVERSE),))
    assert SequenceRecord("AACCGGTTAC", topology="circular", primers=(primer,)).primers == (primer,)
    with pytest.raises(ValueError, match="primer 'p'"):
        SequenceRecord("AACCGGTTAC", primers=(primer,))


def test_extras_a_reader_keeps_for_its_writer_do_not_affect_equality() -> None:
    notes = {"Description": "cloning vector"}
    kept = SequenceRecord("ACGT", name="pUC19", notes=notes, extras={"cache": b"\x00"})
    assert kept == SequenceRecord("ACGT", name="pUC19", notes=notes)
    assert kept != SequenceRecord("ACGT", name="pUC18", notes=notes)


def test_a_span_runs_across_the_origin_when_it_ends_past_the_length() -> None:
    assert across_the_origin(Segment(8, 12), 10)
    assert across_the_origin(BindingSite(8, 12, Strand.REVERSE), 10)
    assert not across_the_origin(Segment(6, 10), 10)


def test_a_feature_runs_across_the_origin_when_its_segments_count_round_past_the_length() -> None:
    assert across_the_origin(_feature((8, 10), (11, 12)), 10)
    assert across_the_origin(_feature((8, 10), (0, 2)), 10)
    assert not across_the_origin(_feature((0, 2), (8, 10)), 10)
    assert not across_the_origin(_feature((0, 4), (3, 6)), 10)


def test_a_span_across_the_origin_covers_the_positions_either_side_of_it() -> None:
    assert CIRCULAR.covers(Segment(8, 12), 1)
    assert CIRCULAR.covers(Segment(8, 12), 11)
    assert not CIRCULAR.covers(Segment(8, 12), 2)
    assert CIRCULAR.covers(_feature((2, 4), (8, 12)), 0)


def test_a_span_covers_only_its_own_positions_on_a_linear_record() -> None:
    assert LINEAR.covers(Segment(2, 6), 5)
    assert not LINEAR.covers(Segment(2, 6), 6)
    assert not LINEAR.covers(Segment(2, 6), 12)


def test_a_span_covers_a_shorter_one_lying_inside_it_across_the_origin() -> None:
    assert CIRCULAR.covers(Segment(7, 13), Segment(9, 12))
    assert CIRCULAR.covers(Segment(7, 13), BindingSite(0, 2, Strand.FORWARD))
    assert not CIRCULAR.covers(Segment(7, 13), Segment(2, 4))
    assert not LINEAR.covers(Segment(2, 6), Segment(4, 7))


def test_bases_read_across_the_origin_counting_any_position_round_the_circle() -> None:
    assert CIRCULAR.bases(8, 12) == "ACAA"
    assert CIRCULAR.bases(-2, 2) == CIRCULAR.bases(18, 22) == "ACAA"
    with pytest.raises(ValueError, match="circular"):
        CIRCULAR.bases(3, 14)


def test_bases_read_a_linear_record_only_inside_it() -> None:
    assert LINEAR.bases(2, 6) == "CCGG"
    assert LINEAR.bases(4, 4) == ""
    with pytest.raises(ValueError, match="linear"):
        LINEAR.bases(8, 12)


def test_a_span_fits_a_circular_record_starting_inside_it_and_at_most_one_turn_long() -> None:
    assert CIRCULAR.fits(8, 12)
    assert CIRCULAR.fits(0, 10)
    assert not CIRCULAR.fits(10, 12)
    assert not CIRCULAR.fits(2, 13)


def test_a_span_fits_a_linear_record_only_inside_it() -> None:
    assert LINEAR.fits(2, 10)
    assert LINEAR.fits(10, 10)
    assert not LINEAR.fits(8, 12)
    assert not LINEAR.fits(-1, 2)


def test_positions_counted_round_in_reading_order_pass_the_length_after_the_origin() -> None:
    assert counted_round((91, 2, 5), 100) == (91, 102, 105)
    assert counted_round((2, 5, 91), 100) == (2, 5, 91)


def test_a_position_a_person_reads_counts_from_one_and_comes_round_the_origin() -> None:
    assert position_text(0, 2686) == "1"
    assert position_text(2683, 2686) == "2684"
    assert position_text(2688, 2686) == "3"
