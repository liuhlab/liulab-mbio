"""The regions a verification judges, read off the junctions a plan tagged."""

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from mbio.io import read_record
from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.snapgene import write_dna
from mbio.verification.judge import JUNCTION_TAG, regions

#: The Golden Gate plan's product, as `scripts/check_examples.py` keeps it current.
PRODUCT = Path(__file__).parents[2] / "docs" / "examples" / "pUC19-GFP" / "product.dna"


@pytest.fixture(scope="module")
def product() -> SequenceRecord:
    return read_record(PRODUCT)


def _seen(record: SequenceRecord) -> list[tuple[str, list[tuple[int, int]], tuple]]:
    """Each region's name, spans and tag, which is what a file must carry back."""
    return [
        (one.name, [(s.start, s.end) for s in one.segments], one.qualifiers.get(JUNCTION_TAG, ()))
        for one in regions(record)
    ]


def test_a_product_gives_its_junctions_and_the_insert_between_them(product):
    assert _seen(product) == [
        ("ATGA junction", [(395, 399)], ("GFP",)),
        ("GFP insert", [(399, 1112)], ()),
        ("TGGC junction", [(1112, 1116)], ("pUC19 backbone",)),
    ]


def test_the_tag_survives_a_dna_file(product, tmp_path):
    write_dna(product, tmp_path / "again.dna")
    assert _seen(read_record(tmp_path / "again.dna")) == _seen(product)


def test_the_tag_survives_a_genbank_file(product, tmp_path):
    features = []
    for one in product.features:
        strand = int(one.strand) or None
        parts = [SimpleLocation(s.start, s.end, strand=strand) for s in one.segments]
        if strand == -1:
            parts.reverse()
        location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        qualifiers = {"label": [one.name], **{k: list(v) for k, v in one.qualifiers.items()}}
        features.append(SeqFeature(location, type=one.type, qualifiers=qualifiers))
    written = SeqRecord(
        Seq(product.sequence),
        id="product",
        features=features,
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    SeqIO.write(written, tmp_path / "product.gb", "genbank")
    assert _seen(read_record(tmp_path / "product.gb")) == _seen(product)


def test_names_given_are_judged_instead(product):
    assert [one.name for one in regions(product, ["GFP", "ori"])] == ["GFP", "ori"]
    with pytest.raises(ValueError, match="no feature named 'mCherry'"):
        regions(product, ["GFP", "mCherry"])


def _junction(start: int, end: int, part: str) -> Feature:
    return Feature("j", "misc_feature", (Segment(start, end),), qualifiers={JUNCTION_TAG: (part,)})


def test_a_repeated_part_takes_an_ordinal_and_an_insert_past_the_origin_comes_back():
    # The junction at 95 runs across the origin, so no stretch holds base 0 and the one after
    # it starts a turn on.
    record = SequenceRecord(
        "A" * 100,
        topology="circular",
        features=(_junction(20, 30, "GFP"), _junction(60, 70, "GFP"), _junction(95, 105, "ori")),
    )
    inserts = [(one.name, one.segments[0]) for one in regions(record) if one.name != "j"]
    assert inserts == [
        ("ori insert", Segment(5, 20)),
        ("GFP insert", Segment(30, 60)),
        ("GFP insert 2", Segment(70, 95)),
    ]


def test_a_linear_record_has_no_stretch_round_from_its_last_junction():
    record = SequenceRecord("A" * 100, features=(_junction(20, 30, "GFP"), _junction(60, 70, "x")))
    assert [one.name for one in regions(record)] == ["j", "GFP insert", "j"]
