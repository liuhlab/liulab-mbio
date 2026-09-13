from pathlib import Path

import pytest

from liulab_mbio.io import read_record
from liulab_mbio.sequence import Segment, Strand

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
