from liulab_mbio.edits import EditReport, delete, insert, replace
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand


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


def test_an_edit_drops_only_the_segment_it_removes_from_a_multi_segment_feature() -> None:
    record = SequenceRecord(BASES, features=(_feature("m", (0, 4), (8, 12)),))
    edited, report = delete(record, 8, 12)
    assert edited.sequence == "AAAACCCC"
    assert edited.features == (_feature("m", (0, 4)),)
    assert report == EditReport(trimmed=(record.features[0],))
