import struct
from pathlib import Path

from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, Strand
from liulab_mbio.snapgene import read_dna

DATA = Path(__file__).parent / "data"
PUC19 = DATA / "pUC19.dna"
GFP = DATA / "GFP.dna"

# Shaped like the primers packet SnapGene writes, with binding sites on the GFP fixture.
PRIMERS_XML = (
    '<?xml version="1.0"?><Primers nextValidID="2"><HybridizationParams'
    ' minContinuousMatchLen="10" allowMismatch="1" minMeltingTemperature="40"'
    ' showAdditionalFivePrimeMatches="1" minimumFivePrimeAnnealing="15"/>'
    '<Primer recentID="0" name="GFP fwd" sequence="CAGTCAGGATCCATGAGTAAAGGAGAAGAACT"'
    ' description="&lt;html&gt;&lt;body&gt;forward&lt;/body&gt;&lt;/html&gt;"'
    ' dateAdded="2023-05-18T20:25:03Z">'
    '<BindingSite location="0-19" boundStrand="0" annealedBases="ATGAGTAAAGGAGAAGAACT"'
    ' meltingTemperature="57"><Component bases="CAGTCAGGATCC"/>'
    '<Component hybridizedRange="0-19" bases="ATGAGTAAAGGAGAAGAACT"/></BindingSite>'
    '<BindingSite simplified="1" location="0-19" boundStrand="0"'
    ' annealedBases="ATGAGTAAAGGAGAAGAACT" meltingTemperature="57">'
    '<Component hybridizedRange="0-19" bases="ATGAGTAAAGGAGAAGAACT"/></BindingSite></Primer>'
    '<Primer recentID="1" name="GFP rev" sequence="CTGACAGAATTCTTTGTATAGTTCATCCATGG"'
    ' description=""><BindingSite location="697-716" boundStrand="1"'
    ' annealedBases="CCATGGATGAACTATACAAA" meltingTemperature="55">'
    '<Component bases="CTGACAGAATTC"/>'
    '<Component hybridizedRange="697-716" bases="CCATGGATGAACTATACAAA"/>'
    "</BindingSite></Primer></Primers>"
).encode()


def _rewritten(path: Path, kind: int, payload: bytes) -> bytes:
    """The file's packets with every packet of `kind` replaced by one holding `payload`."""
    data, out, offset = path.read_bytes(), bytearray(), 0
    while offset < len(data):
        packet_kind, length = struct.unpack_from(">BI", data, offset)
        packet = data[offset : offset + 5 + length]
        out += struct.pack(">BI", kind, len(payload)) + payload if packet_kind == kind else packet
        offset += 5 + length
    return bytes(out)


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


def test_read_dna_converts_binding_sites_and_ignores_simplified_duplicates(tmp_path: Path) -> None:
    path = tmp_path / "primers.dna"
    path.write_bytes(_rewritten(GFP, 0x05, PRIMERS_XML))
    assert read_dna(path).primers == (
        Primer(
            "GFP fwd",
            "CAGTCAGGATCCATGAGTAAAGGAGAAGAACT",
            binding_sites=(BindingSite(0, 20, Strand.FORWARD),),
            description="<html><body>forward</body></html>",
        ),
        Primer(
            "GFP rev",
            "CTGACAGAATTCTTTGTATAGTTCATCCATGG",
            binding_sites=(BindingSite(697, 717, Strand.REVERSE),),
        ),
    )


def test_read_dna_takes_the_notes_and_reads_the_map_label_as_the_name() -> None:
    record = read_dna(PUC19)
    assert record.name == "pUC19"
    assert read_dna(GFP).name == "GFP"
    assert record.notes["Type"] == "Synthetic"
    assert record.notes["CreatedBy"] == "New England Biolabs"
    assert record.notes["Comments"] == "See also GenBank accession L09137."
    assert record.notes["Description"].startswith("<html><body>Standard <i>E. coli</i> vector")
    assert "CustomMapLabel" not in record.notes


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
