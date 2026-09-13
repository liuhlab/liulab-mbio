import dataclasses
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

from liulab_mbio.edits import delete
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand
from liulab_mbio.snapgene import read_dna, write_dna

DATA = Path(__file__).parent / "data"
PUC19 = DATA / "pUC19.dna"
GFP = DATA / "GFP.dna"

#: The packets the model holds, which the writer builds rather than keeping.
MODELLED = (0x00, 0x05, 0x06, 0x0A)


def _packets(data: bytes) -> list[tuple[int, bytes]]:
    packets, offset = [], 0
    while offset < len(data):
        kind, length = struct.unpack_from(">BI", data, offset)
        packets.append((kind, data[offset + 5 : offset + 5 + length]))
        offset += 5 + length
    return packets


def _flags(data: bytes) -> int:
    return next(payload[0] for kind, payload in _packets(data) if kind == 0x00)


# Shaped like the primers packet SnapGene writes, with binding sites on the GFP fixture.
PRIMERS_XML = (
    b'<?xml version="1.0"?><Primers nextValidID="2"><HybridizationParams'
    b' minContinuousMatchLen="10" allowMismatch="1" minMeltingTemperature="40"'
    b' showAdditionalFivePrimeMatches="1" minimumFivePrimeAnnealing="15"/>'
    b'<Primer recentID="0" name="GFP fwd" sequence="CAGTCAGGATCCATGAGTAAAGGAGAAGAACT"'
    b' description="&lt;html&gt;&lt;body&gt;forward&lt;/body&gt;&lt;/html&gt;"'
    b' dateAdded="2023-05-18T20:25:03Z">'
    b'<BindingSite location="0-19" boundStrand="0" annealedBases="ATGAGTAAAGGAGAAGAACT"'
    b' meltingTemperature="57"><Component bases="CAGTCAGGATCC"/>'
    b'<Component hybridizedRange="0-19" bases="ATGAGTAAAGGAGAAGAACT"/></BindingSite>'
    b'<BindingSite simplified="1" location="0-19" boundStrand="0"'
    b' annealedBases="ATGAGTAAAGGAGAAGAACT" meltingTemperature="57">'
    b'<Component hybridizedRange="0-19" bases="ATGAGTAAAGGAGAAGAACT"/></BindingSite></Primer>'
    b'<Primer recentID="1" name="GFP rev" sequence="CTGACAGAATTCTTTGTATAGTTCATCCATGG"'
    b' description=""><BindingSite location="697-716" boundStrand="1"'
    b' annealedBases="CCATGGATGAACTATACAAA" meltingTemperature="55">'
    b'<Component bases="CTGACAGAATTC"/>'
    b'<Component hybridizedRange="697-716" bases="CCATGGATGAACTATACAAA"/>'
    b"</BindingSite></Primer></Primers>"
)


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


def test_write_dna_round_trips_a_record_through_a_file(tmp_path: Path) -> None:
    for fixture in (PUC19, GFP):
        record = read_dna(fixture)
        path = tmp_path / fixture.name
        write_dna(record, path)
        assert read_dna(path) == record


def test_write_dna_keeps_the_packets_the_model_does_not_hold(tmp_path: Path) -> None:
    path = tmp_path / "same.dna"
    write_dna(read_dna(PUC19), path)
    before, after = _packets(PUC19.read_bytes()), _packets(path.read_bytes())
    assert [p for p in after if p[0] not in MODELLED] == [p for p in before if p[0] not in MODELLED]
    assert dict(after)[0x00] == dict(before)[0x00]


def test_write_dna_drops_the_cut_site_cache_once_the_sequence_changes(tmp_path: Path) -> None:
    edited, _ = delete(read_dna(PUC19), 100, 200)
    path = tmp_path / "edited.dna"
    write_dna(edited, path)
    assert [kind for kind, _ in _packets(path.read_bytes())] == [
        0x09,  # cookie
        0x00,  # sequence
        0x08,  # end properties
        0x0A,  # features
        0x05,  # primers
        0x06,  # notes
        0x0D,  # display settings
        0x1C,  # enzyme visibilities
        0x0E,  # custom enzyme sets
    ]


def test_the_sequence_packet_flags_hold_topology_strandedness_and_methylation(
    tmp_path: Path,
) -> None:
    assert _flags(PUC19.read_bytes()) == 0x1F  # circular, double, Dam, Dcm and EcoKI methylated
    assert _flags(GFP.read_bytes()) == 0x02  # linear and double-stranded
    linear = dataclasses.replace(read_dna(PUC19), topology="linear")
    write_dna(linear, tmp_path / "linear.dna")
    assert _flags((tmp_path / "linear.dna").read_bytes()) == 0x1E
    write_dna(SequenceRecord("ACGT"), tmp_path / "new.dna")
    assert _flags((tmp_path / "new.dna").read_bytes()) == 0x02


def test_write_dna_writes_a_span_across_the_origin_as_a_wrapped_range(tmp_path: Path) -> None:
    site = Feature(
        "BsmBI",
        "misc_feature",
        (Segment(2682, 2688),),
        strand=Strand.FORWARD,
        color="#ff0000",
    )
    record = dataclasses.replace(read_dna(PUC19), features=(site,))
    path = tmp_path / "wrapped.dna"
    write_dna(record, path)
    assert b'range="2683-2"' in path.read_bytes()
    assert read_dna(path).features == (site,)


def test_write_dna_fills_a_hole_between_segments_with_a_gap_segment(tmp_path: Path) -> None:
    feature = Feature(
        "exons", "mRNA", (Segment(10, 20), Segment(30, 40)), strand=Strand.FORWARD, color="#ff0000"
    )
    record = dataclasses.replace(read_dna(GFP), features=(feature,))
    path = tmp_path / "gap.dna"
    write_dna(record, path)
    assert b'range="21-30" color="noColor" type="gap"' in path.read_bytes()
    assert read_dna(path).features == (feature,)


def test_write_dna_writes_primers_that_were_built_in_code(tmp_path: Path) -> None:
    primer = Primer(
        "GFP fwd",
        "CAGTCAGGATCCATGAGTAAAGGAGAAGAACT",
        binding_sites=(BindingSite(0, 20, Strand.FORWARD),),
        description="forward",
    )
    record = dataclasses.replace(read_dna(GFP), primers=(primer,))
    path = tmp_path / "primers.dna"
    write_dna(record, path)
    assert b'annealedBases="ATGAGTAAAGGAGAAGAACT"' in path.read_bytes()
    assert read_dna(path).primers == (primer,)


def test_write_dna_pairs_each_binding_site_with_a_simplified_copy(tmp_path: Path) -> None:
    primer = Primer(
        "GFP fwd",
        "CAGTCAGGATCCATGAGTAAAGGAGAAGAACT",
        binding_sites=(BindingSite(0, 20, Strand.FORWARD),),
    )
    path = tmp_path / "primers.dna"
    write_dna(dataclasses.replace(read_dna(GFP), primers=(primer,)), path)
    (packet,) = (payload for kind, payload in _packets(path.read_bytes()) if kind == 0x05)
    site, simplified = ET.fromstring(packet).iter("BindingSite")
    assert simplified.attrib == {"simplified": "1", **site.attrib}
    assert [component.attrib for component in simplified] == [
        component.attrib for component in site
    ]


def test_write_dna_has_snapgene_label_the_map_with_the_name(tmp_path: Path) -> None:
    path = tmp_path / "file-name.dna"
    for record in (
        SequenceRecord("ACGT", name="probe"),
        dataclasses.replace(read_dna(PUC19), name="renamed"),
    ):
        write_dna(record, path)
        notes = ET.fromstring(dict(_packets(path.read_bytes()))[0x06])
        assert [node.text for node in notes.iter("UseCustomMapLabel")] == ["1"]


def test_the_reader_agrees_with_biopython_on_both_fixtures() -> None:
    from Bio import SeqIO

    for fixture in (PUC19, GFP):
        record, reference = read_dna(fixture), SeqIO.read(fixture, "snapgene")
        assert record.sequence == str(reference.seq).upper()
        assert record.topology == reference.annotations["topology"]
        assert {(f.name, f.segments[0].start, f.segments[-1].end) for f in record.features} == {
            (
                f.qualifiers["label"][0],
                min(int(part.start) for part in f.location.parts),
                max(int(part.end) for part in f.location.parts),
            )
            for f in reference.features
        }


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
