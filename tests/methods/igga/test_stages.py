"""What a round is plated on, read off the vector's own marker rather than the paper's."""

from mbio.sequence import Feature, Segment, SequenceRecord
from synbio.igga import stages


def marked(marker: str | None) -> SequenceRecord:
    """A vector annotating `marker` as its selection marker, or nothing where it is `None`."""
    features = (Feature(marker, "CDS", (Segment(10, 100),)),) if marker else ()
    return SequenceRecord("TA" * 200, topology="circular", name="vector", features=features)


def test_a_round_is_plated_on_the_vectors_own_marker_and_not_the_papers() -> None:
    """The method rebuilt its DMX vector KanR, so the published carbenicillin does not carry."""
    assert stages.selection_for(marked("KanR")) == "50 µg/mL kanamycin"
    assert stages.selection_for(marked("AmpR")) == "100 µg/mL carbenicillin"


def test_a_vector_naming_no_marker_names_no_drug_and_keeps_the_hole() -> None:
    """H22 stands only where nothing in the record answers it."""
    assert stages.selection_for(marked(None)) == ""
    assert stages.ROUND_SELECTION in stages.holes_for(marked(None))
    assert stages.ROUND_SELECTION not in stages.holes_for(marked("KanR"))
    assert next(one.id for one in stages.holes_for(marked(None))) == "H22"
