from pathlib import Path

from liulab_mbio.sequence import Feature, Segment, Strand
from liulab_mbio.snapgene import read_dna

DATA = Path(__file__).parent / "data"
PUC19 = DATA / "pUC19.dna"
GFP = DATA / "GFP.dna"


def _feature(record_features: tuple[Feature, ...], name: str) -> Feature:
    (feature,) = (f for f in record_features if f.name == name)
    return feature


def test_read_dna_gives_the_sequence_and_topology() -> None:
    puc19 = read_dna(PUC19)
    gfp = read_dna(GFP)
    assert (len(puc19), puc19.topology) == (2686, "circular")
    assert puc19.sequence.startswith("TCGCGCGTTTCGGTGATGACGG")
    assert (len(gfp), gfp.topology) == (717, "linear")
    assert gfp.sequence.startswith("ATGAGTAAAGGAGAAGAACTTTTCACTGG")


def test_read_dna_converts_feature_ranges_to_half_open_segments() -> None:
    features = read_dna(PUC19).features
    assert [f.name for f in features] == [
        "lac operator",
        "M13 rev",
        "M13 fwd",
        "lac promoter",
        "MCS",
        "AmpR promoter",
        "lacZα",  # noqa: RUF001 - the fixture's feature name
        "ori",
        "AmpR",
    ]
    m13_fwd = _feature(features, "M13 fwd")
    assert (m13_fwd.type, m13_fwd.strand, m13_fwd.segments) == (
        "primer_bind",
        Strand.FORWARD,
        (Segment(378, 395),),
    )
    assert _feature(features, "lacZα").segments == (Segment(145, 469),)  # noqa: RUF001
    assert _feature(read_dna(GFP).features, "GFP").segments == (Segment(0, 717),)


def test_read_dna_maps_directionality_to_strand() -> None:
    features = read_dna(PUC19).features
    assert _feature(features, "M13 rev").strand == Strand.REVERSE
    assert _feature(features, "MCS").strand == Strand.NONE


def test_read_dna_keeps_named_segments_and_colours() -> None:
    features = read_dna(PUC19).features
    promoter = _feature(features, "lac promoter")
    assert promoter.color == "#ffffff"
    assert promoter.segments == (
        Segment(512, 519, name="-10"),
        Segment(519, 537),
        Segment(537, 543, name="-35"),
    )
    amp = _feature(features, "AmpR")
    assert (amp.color, amp.segments) == (
        "#ccffcc",
        (Segment(1625, 2417), Segment(2417, 2486, name="signal sequence")),
    )


def test_read_dna_keeps_qualifiers_with_their_integer_and_text_values() -> None:
    amp = _feature(read_dna(PUC19).features, "AmpR")
    assert amp.qualifiers["codon_start"] == (1,)
    assert amp.qualifiers["product"] == ("<html><body>β-lactamase</body></html>",)
    assert set(amp.qualifiers) == {
        "codon_start",
        "gene",
        "note",
        "product",
        "transl_table",
        "translation",
    }
