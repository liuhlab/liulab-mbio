"""What a product says about itself, read across the origin of a circular record.

Each record is 100 bases written in code, so every span below is read off it by eye.
"""

from liulab_mbio.bench.phenotype import read_phenotype, selection_marker
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand


def plasmid(*features: Feature) -> SequenceRecord:
    """A circular record of 100 bases carrying these features."""
    return SequenceRecord("ACGT" * 25, topology="circular", features=features)


def feature(name: str, kind: str, *spans: tuple[int, int]) -> Feature:
    """A forward feature over those spans."""
    return Feature(
        name, kind, tuple(Segment(start, end) for start, end in spans), strand=Strand.FORWARD
    )


def test_a_coding_sequence_past_the_origin_inside_an_insert_across_it_is_the_insert_s():
    product = plasmid(feature("GFP", "CDS", (2, 8)))
    read = read_phenotype(product, (90, 110), vector=plasmid(), span=(90, 91))
    assert read.coding is not None
    assert read.coding.name == "GFP"


def test_a_promoter_across_the_origin_is_measured_from_its_reading_end():
    # It reads 90..98 then 100..108, so it ends at base 8 and the insert at 20 is 12 bases on.
    product = plasmid(
        feature("P", "promoter", (90, 98), (100, 108)), feature("RBS", "RBS", (10, 15))
    )
    read = read_phenotype(product, (20, 40), vector=plasmid(), span=(20, 21))
    assert (read.gap_bp, read.ribosome_binding_site) == (12, True)


def test_a_coding_sequence_meeting_a_replaced_span_past_the_origin_is_interrupted():
    # The span runs 95..105, so bases 0..4 of the vector are replaced too.
    for spans in [((0, 50),), ((2, 8), (20, 30))]:
        vector = plasmid(feature("lacZ", "CDS", *spans))
        read = read_phenotype(plasmid(), (20, 40), vector=vector, span=(95, 105))
        assert read.reporter is not None
        assert read.reporter.name == "lacZ"


def test_a_marker_inside_a_replaced_span_past_the_origin_is_passed_over():
    vector = plasmid(feature("AmpR", "CDS", (0, 50)), feature("KanR", "CDS", (60, 80)))
    marker = selection_marker(vector, outside=(95, 105))
    assert marker is not None
    assert marker.name == "KanR"


def test_an_insertion_that_replaces_nothing_breaks_a_coding_sequence_only_inside_it():
    # KanR runs across the origin; an insertion at AmpR's first base leaves AmpR whole.
    vector = plasmid(feature("AmpR", "CDS", (10, 60)), feature("KanR", "CDS", (90, 110)))
    for at, broken, marker in ((30, "AmpR", "KanR"), (0, "KanR", "AmpR"), (10, None, "AmpR")):
        read = read_phenotype(vector, (0, 0), vector=vector, span=(at, at))
        assert (_name(read.reporter), _name(read.marker)) == (broken, marker), at


def _name(one: Feature | None) -> str | None:
    return None if one is None else one.name
