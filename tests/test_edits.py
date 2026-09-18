import pytest

from liulab_mbio.edits import (
    EditReport,
    annealed,
    carried,
    delete,
    flipped,
    insert,
    replace,
    rotate,
)
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand


def _feature(name: str, *spans: tuple[int, int], strand: Strand = Strand.FORWARD) -> Feature:
    segments = tuple(Segment(start, end) for start, end in spans)
    return Feature(name, "misc_feature", segments, strand=strand)


A, C, G = _feature("a", (0, 4)), _feature("c", (4, 8)), _feature("g", (8, 12))
BASES = "AAAACCCCGGGG"
RECORD = SequenceRecord(BASES, features=(A, C, G))


def test_insert_shifts_features_after_it_and_leaves_those_before() -> None:
    record = SequenceRecord(
        "AAAACCCCGGGG", features=(_feature("a", (0, 4)), _feature("g", (8, 12)))
    )
    edited, report = insert(record, 6, "TT")
    assert edited.sequence == "AAAACCTTCCGGGG"
    assert edited.features == (_feature("a", (0, 4)), _feature("g", (10, 14)))
    assert report == EditReport()


def test_delete_trims_the_features_it_cuts_and_reports_them() -> None:
    edited, report = delete(RECORD, 2, 6)
    assert edited.sequence == "AACCGGGG"
    assert edited.features == (_feature("a", (0, 2)), _feature("c", (2, 4)), _feature("g", (4, 8)))
    assert report == EditReport(trimmed=(A, C))


def test_delete_drops_a_feature_it_removes_entirely() -> None:
    edited, report = delete(RECORD, 3, 9)
    assert edited.sequence == "AAAGGG"
    assert edited.features == (_feature("a", (0, 3)), _feature("g", (3, 6)))
    assert report == EditReport(trimmed=(A, G), dropped=(C,))


def test_an_edit_inside_a_feature_keeps_it_spanning_the_new_bases() -> None:
    edited, report = replace(RECORD, 5, 6, "TTT")
    assert edited.sequence == "AAAACTTTCCGGGG"
    assert edited.features == (A, _feature("c", (4, 10)), _feature("g", (10, 14)))
    assert report == EditReport(changed=(C,))


def test_an_edit_across_the_origin_keeps_the_remaining_bases_and_moves_the_origin() -> None:
    record = SequenceRecord(BASES, topology="circular", features=(G,))
    edited, report = delete(record, 10, 14)
    assert edited.sequence == "AACCCCGG"
    assert edited.features == (_feature("g", (6, 8)),)
    assert report == EditReport(trimmed=(G,))


def test_a_feature_across_the_origin_is_shifted_by_an_edit_elsewhere() -> None:
    across = _feature("across", (10, 14))
    record = SequenceRecord(BASES, topology="circular", features=(across,))
    edited, report = delete(record, 4, 6)
    assert edited.sequence == "AAAACCGGGG"
    assert edited.features == (_feature("across", (8, 12)),)
    assert report == EditReport()


def test_rotate_moves_the_origin_and_carries_spans_across_it() -> None:
    primer = Primer("p", "CCCC", binding_sites=(BindingSite(4, 8, Strand.FORWARD),))
    record = SequenceRecord(BASES, topology="circular", features=(C,), primers=(primer,))
    turned = rotate(record, 6)
    assert turned.sequence == "CCGGGGAAAACC"
    assert turned.features == (_feature("c", (10, 14)),)
    assert turned.primers[0].binding_sites == (BindingSite(10, 14, Strand.FORWARD),)
    assert rotate(turned, 6) == record


def test_rotate_refuses_a_linear_record() -> None:
    with pytest.raises(ValueError, match="circular"):
        rotate(RECORD, 3)


def test_an_edit_drops_the_binding_sites_it_overlaps_and_shifts_the_rest() -> None:
    primer = Primer(
        "p",
        "CCCC",
        binding_sites=(BindingSite(4, 8, Strand.FORWARD), BindingSite(8, 12, Strand.REVERSE)),
    )
    record = SequenceRecord(BASES, primers=(primer,))
    edited, report = delete(record, 5, 7)
    assert edited.primers[0].binding_sites == (BindingSite(6, 10, Strand.REVERSE),)
    assert report == EditReport(dropped_sites=((primer, BindingSite(4, 8, Strand.FORWARD)),))


def test_an_edit_refuses_a_span_that_does_not_fit() -> None:
    with pytest.raises(ValueError, match="linear"):
        delete(RECORD, 8, 13)


def test_an_edit_drops_only_the_segment_it_removes_from_a_multi_segment_feature() -> None:
    record = SequenceRecord(BASES, features=(_feature("m", (0, 4), (8, 12)),))
    edited, report = delete(record, 8, 12)
    assert edited.sequence == "AAAACCCC"
    assert edited.features == (_feature("m", (0, 4)),)
    assert report == EditReport(trimmed=(record.features[0],))


def test_flipped_turns_every_span_and_strand_end_for_end() -> None:
    primer = Primer("p", "CCCC", binding_sites=(BindingSite(4, 8, Strand.FORWARD),))
    record = SequenceRecord(
        BASES, features=(_feature("m", (0, 4), (8, 12), strand=Strand.REVERSE),), primers=(primer,)
    )
    turned = flipped(record)
    assert turned.sequence == "CCCCGGGGTTTT"
    assert turned.features == (_feature("m", (0, 4), (8, 12), strand=Strand.FORWARD),)
    assert turned.primers[0].binding_sites == (BindingSite(4, 8, Strand.REVERSE),)
    assert flipped(turned) == record


def test_flipped_keeps_a_span_across_the_origin_across_it() -> None:
    primer = Primer("p", "GAAA", binding_sites=(BindingSite(11, 15, Strand.FORWARD),))
    record = SequenceRecord(
        BASES, topology="circular", features=(_feature("across", (9, 14)),), primers=(primer,)
    )
    turned = flipped(record)
    # The feature read GGGAA and the primer site GAAA; turned over, each reads the other strand.
    assert turned.features == (_feature("across", (10, 15), strand=Strand.REVERSE),)
    assert turned.extract(Segment(10, 15)) == "TTCCC"
    assert turned.primers[0].binding_sites == (BindingSite(9, 13, Strand.REVERSE),)
    assert turned.extract(Segment(9, 13)) == "TTTC"
    assert flipped(turned) == record


def test_carried_cuts_a_feature_down_to_the_span_and_shifts_it_by_the_offset() -> None:
    features, primers = carried(RECORD, 2, 10, offset=-2)
    assert features == (_feature("a", (0, 2)), _feature("c", (2, 6)), _feature("g", (6, 8)))
    assert primers == ()


def test_a_feature_meeting_a_span_across_the_origin_twice_keeps_a_segment_for_each() -> None:
    record = SequenceRecord(BASES, topology="circular", features=(_feature("split", (4, 10)),))
    features, _ = carried(record, 8, 18, offset=-8)
    assert features == (_feature("split", (0, 2), (8, 10)),)


def test_carried_keeps_a_primer_only_where_a_whole_binding_site_survives() -> None:
    whole = Primer("whole", "AAAA", binding_sites=(BindingSite(0, 4, Strand.FORWARD),))
    cut = Primer("cut", "CCCC", binding_sites=(BindingSite(6, 10, Strand.REVERSE),))
    record = SequenceRecord(BASES, topology="circular", primers=(whole, cut))
    _, primers = carried(record, 8, 18, offset=-8)
    # "cut" annealed across the edge of the span, so it has nowhere left to sit.
    assert [one.name for one in primers] == ["whole"]
    assert primers[0].binding_sites == (BindingSite(4, 8, Strand.FORWARD),)


def test_annealed_puts_a_primer_at_the_edge_its_strand_reads_from() -> None:
    primer = Primer("p", "GGGGAACCGG", binding_sites=(BindingSite(0, 6, Strand.FORWARD),))
    forward = annealed(primer, 4, Strand.FORWARD)
    reverse = annealed(primer, 20, Strand.REVERSE)
    assert forward.binding_sites == (BindingSite(4, 10, Strand.FORWARD),)
    assert reverse.binding_sites == (BindingSite(14, 20, Strand.REVERSE),)
