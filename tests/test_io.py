import random
from pathlib import Path

import pytest

from liulab_mbio.io import read_record, read_region
from liulab_mbio.sequence import Segment, Strand

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
