import pytest

from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement


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
