import random
from pathlib import Path

import pytest

from liulab_mbio.io import read_record, read_region
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.snapgene import write_dna

from .fasta import write_fasta

GENBANK = """LOCUS       mini                      12 bp    DNA     circular SYN 01-JAN-2020
DEFINITION  a tiny circular record.
FEATURES             Location/Qualifiers
     misc_feature    join(11..12,1..2)
                     /label="across"
     CDS             complement(5..8)
                     /label="rev"
                     /codon_start=1
ORIGIN
        1 aaaaccccgggg
//
"""


def test_read_record_reads_a_snapgene_file(gfp_file: Path) -> None:
    record = read_record(gfp_file)
    assert (len(record), record.topology, record.name) == (717, "linear", "GFP")
    assert record.features[0].name == "GFP"


def test_read_record_reads_genbank_topology_and_features(tmp_path: Path) -> None:
    path = tmp_path / "mini.gb"
    path.write_text(GENBANK)
    record = read_record(path)
    assert (record.sequence, record.topology, record.name) == ("AAAACCCCGGGG", "circular", "mini")
    across, reverse = record.features
    assert (across.name, across.segments, across.type) == (
        "across",
        (Segment(10, 14),),
        "misc_feature",
    )
    assert (reverse.name, reverse.segments, reverse.strand) == (
        "rev",
        (Segment(4, 8),),
        Strand.REVERSE,
    )
    assert reverse.qualifiers == {"codon_start": (1,)}
    assert record.notes["Description"] == "a tiny circular record"


#: 100 bp, circular, numbered as GenBank does: bases 2 to 8 are C, 91 to 98 are G, the rest A.
GAPPED = """LOCUS       gapped                   100 bp    DNA     circular SYN 01-JAN-2020
FEATURES             Location/Qualifiers
     misc_feature    {location}
                     /label="g"
{qualifiers}ORIGIN
        1 acccccccaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
       61 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa ggggggggaa
//
"""


def _gapped(tmp_path: Path, location: str, qualifiers: str = "") -> SequenceRecord:
    """Read `GAPPED` with its feature at `location`, and these qualifier lines after its label."""
    path = tmp_path / "gapped.gb"
    path.write_text(GAPPED.format(location=location, qualifiers=qualifiers))
    return read_record(path)


def test_read_record_keeps_a_join_across_the_origin_in_the_order_it_reads(
    tmp_path: Path,
) -> None:
    record = _gapped(tmp_path, "join(91..98,2..8)")
    (feature,) = record.features
    assert feature.segments == (Segment(90, 98), Segment(1, 8))
    assert record.extract(feature) == "GGGGGGGGCCCCCCC"
    record = _gapped(tmp_path, "complement(join(91..98,2..8))")
    (feature,) = record.features
    assert (feature.segments, feature.strand) == ((Segment(90, 98), Segment(1, 8)), Strand.REVERSE)
    assert record.extract(feature) == "GGGGGGGCCCCCCCC"


def test_read_record_joins_two_parts_only_where_the_feature_reads_across_the_origin(
    tmp_path: Path,
) -> None:
    record = _gapped(tmp_path, "complement(join(91..100,1..8))")
    (feature,) = record.features
    assert feature.segments == (Segment(90, 108),)
    assert record.extract(feature) == "GGGGGGGTTTCCCCCCCC"
    # These two meet at the origin only if read the other way round.
    assert _gapped(tmp_path, "join(1..8,91..100)").features[0].segments == (
        Segment(0, 8),
        Segment(90, 100),
    )


def test_read_record_colours_a_segment_list_across_the_origin(tmp_path: Path) -> None:
    record = _gapped(
        tmp_path,
        "join(91..98,2..8)",
        """                     /note="This feature has 2 segments:
                        1: 91 .. 98 / #ff0000
                        2: 2 .. 8 / #00ff00"
""",
    )
    (feature,) = record.features
    assert (feature.color, feature.segments, feature.qualifiers) == (
        "#ff0000",
        (Segment(90, 98), Segment(1, 8, color="#00ff00")),
        {},
    )


def _named(record: SequenceRecord, name: str) -> Feature:
    (feature,) = (f for f in record.features if f.name == name)
    return feature


def test_read_record_reads_the_colours_snapgene_exports_to_genbank(data_dir: Path) -> None:
    record = read_record(data_dir / "colour-test-snapgene.gbk")
    colours = {feature.name: (feature.color, feature.segments) for feature in record.features}
    assert colours["split"] == (
        "#ff0000",
        (Segment(100, 160), Segment(200, 260, color="#00ff00"), Segment(300, 360, color="#0000ff")),
    )
    assert colours["split-same-first"] == (
        "#ffcc00",
        (Segment(400, 420), Segment(440, 460, color="#123abc")),
    )
    assert colours["lac promoter"] == (
        "#ffffff",
        (Segment(512, 519, name="-10"), Segment(519, 537), Segment(537, 543, name="-35")),
    )
    assert colours["across"] == ("#993366", (Segment(2600, 2750),))
    assert (colours["white"][0], colours["M13 fwd"][0]) == ("#ffffff", "#a020f0")
    assert record.features[0].color is None  # the source feature, which has no colour note
    amp = _named(record, "AmpR")
    assert (amp.color, amp.segments) == (
        "#ccffcc",
        (Segment(1625, 2417), Segment(2417, 2486, name="signal sequence")),
    )
    assert amp.qualifiers["note"] == (
        "confers resistance to ampicillin, carbenicillin, and related antibiotics",
        "Cleavage site after base 2417",
    )
    assert _named(record, "ori").qualifiers["direction"] == ("LEFT",)
    assert not [
        note
        for feature in record.features
        for note in feature.qualifiers.get("note", ())
        if "#" in str(note)
    ]


def test_a_dna_file_written_from_the_export_keeps_its_colours_and_no_colour_note(
    data_dir: Path, tmp_path: Path
) -> None:
    record = read_record(data_dir / "colour-test-snapgene.gbk")
    path = tmp_path / "colour-test.dna"
    write_dna(record, path)
    # All but the source feature, which has no colour until the writer gives it SnapGene's grey.
    assert read_record(path).features[1:] == record.features[1:]
    assert b"color:" not in path.read_bytes()
    assert b"segments:" not in path.read_bytes()


#: Colour notes SnapGene would not write, beside the colour qualifiers of other tools.
NOTED = """LOCUS       noted                     12 bp    DNA     circular SYN 01-JAN-2020
FEATURES             Location/Qualifiers
     misc_feature    1..4
                     /note="color: #FFCC00"
                     /note="color: black"
                     /note="its color: #123456 was chosen by hand"
                     /note="color: #12ab3C"
     misc_feature    join(5..6,9..10)
                     /note="This feature has 2 segments:
                        1: 5 .. 6 / #ff0000
                        2: 9 .. 10 / black"
                     /ApEinfo_fwdcolor="#ff0000"
                     /ApEinfo_revcolor="#00ff00"
                     /color="#0000ff"
ORIGIN
        1 aaaaccccgggg
//
"""


def test_read_record_reads_only_a_whole_colour_note_of_hex_colours(tmp_path: Path) -> None:
    path = tmp_path / "noted.gb"
    path.write_text(NOTED)
    noted, listed = read_record(path).features
    assert noted.color == "#12ab3C"
    assert noted.qualifiers == {"note": ("color: black", "its color: #123456 was chosen by hand")}
    assert listed.color is None
    assert listed.qualifiers == {
        "note": ("This feature has 2 segments: 1: 5 .. 6 / #ff0000 2: 9 .. 10 / black",),
        "ApEinfo_fwdcolor": ("#ff0000",),
        "ApEinfo_revcolor": ("#00ff00",),
        "color": ("#0000ff",),
    }


def test_read_record_reads_fasta_as_a_linear_record(tmp_path: Path) -> None:
    path = tmp_path / "gfp.fa"
    path.write_text(">gfp\nAAAACCCCGGGG\n")
    record = read_record(path)
    assert (record.sequence, record.topology, record.name) == ("AAAACCCCGGGG", "linear", "gfp")
    assert record.features == ()
    assert dict(record.notes) == {}


def test_read_record_refuses_a_suffix_it_does_not_know(tmp_path: Path) -> None:
    path = tmp_path / "mini.xyz"
    path.write_text("")
    with pytest.raises(ValueError, match="xyz"):
        read_record(path)


#: Two records, one soft-masked and one not, written ten bases to a line.
GENOME = {
    "chrI__ce11": "".join(random.Random(7).choice("acgt") for _ in range(95)),
    "chrII": "ACGT" * 10,
}
WIDTH = 10


@pytest.fixture(scope="module")
def indexed(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The two records as a FASTA, with the index beside them."""
    path = tmp_path_factory.mktemp("genome") / "planted.fa"
    write_fasta(path, GENOME, width=WIDTH)
    return path


def test_read_region_reads_one_span_of_one_sequence(indexed: Path) -> None:
    record = read_region(indexed, "chrI__ce11", 23, 71)
    assert record.sequence == GENOME["chrI__ce11"][23:71].upper()
    assert (record.name, record.topology) == ("chrI__ce11", "linear")
    assert read_region(indexed, "chrII", 0, 40).sequence == GENOME["chrII"]


def test_read_region_clips_a_span_to_the_sequence(indexed: Path) -> None:
    assert read_region(indexed, "chrII", -5, 12).sequence == GENOME["chrII"][:12]
    assert read_region(indexed, "chrII", 32, 500).sequence == GENOME["chrII"][32:]
    with pytest.raises(ValueError, match="none of it"):
        read_region(indexed, "chrII", 40, 60)


def test_read_region_reads_no_byte_outside_the_span(tmp_path: Path) -> None:
    path = tmp_path / "planted.fa"
    write_fasta(path, GENOME, width=WIDTH)
    # Every line but the ones bases 20 to 79 lie on is nonsense the model would refuse, so a
    # reader of the whole sequence could not answer at all.
    lines = path.read_text().splitlines()
    path.write_text(
        "".join(
            f"{line}\n" if line.startswith(">") or index in range(3, 9) else "?" * WIDTH + "\n"
            for index, line in enumerate(lines)
        )
    )
    assert read_region(path, "chrI__ce11", 23, 71).sequence == GENOME["chrI__ce11"][23:71].upper()


def test_read_region_refuses_a_fasta_with_no_index(tmp_path: Path) -> None:
    path = tmp_path / "planted.fa"
    write_fasta(path, GENOME, width=WIDTH)
    path.with_name(path.name + ".fai").unlink()
    with pytest.raises(FileNotFoundError, match="genome assembly register"):
        read_region(path, "chrII", 0, 10)


def test_read_region_refuses_a_name_the_index_does_not_hold(indexed: Path) -> None:
    with pytest.raises(ValueError, match="chrIII"):
        read_region(indexed, "chrIII", 0, 10)
