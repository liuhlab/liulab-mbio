"""The regions a verification judges, read off the junctions a plan tagged."""

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from mbio.cloning.restriction import plan_restriction
from mbio.edits import rotate
from mbio.io import read_record
from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.snapgene import write_dna
from mbio.verification.judge import BACKBONE, JUNCTION_TAG, regions


@pytest.fixture(scope="module")
def product(plan) -> SequenceRecord:
    """The shared Golden Gate plan's product: GFP in the pUC19 cloning site."""
    return plan.product


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
        ("TGGC junction", [(1112, 1116)], (BACKBONE,)),
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


def test_a_repeated_part_takes_an_ordinal():
    record = SequenceRecord(
        "A" * 100,
        topology="circular",
        features=(
            _junction(20, 30, "GFP"),
            _junction(40, 50, "GFP"),
            _junction(60, 70, BACKBONE),
        ),
    )
    inserts = [(one.name, one.segments[0]) for one in regions(record) if one.name != "j"]
    assert inserts == [("GFP insert", Segment(30, 40)), ("GFP insert 2", Segment(50, 60))]


@pytest.mark.parametrize("turn", [444, 450])
def test_a_product_whose_origin_lies_in_a_junction_judges_the_insert_and_not_the_vector(
    puc19, gfp, turn
):
    # pUC19 turned so its own first base lands inside a junction's bases: at 444 inside the
    # junction the insert follows, at 450 inside the one the vector follows. The product keeps
    # that origin, so where base 0 falls says nothing about which side the vector lies on.
    product = plan_restriction(rotate(puc19, turn), gfp).product
    judged = sorted((one.name, one.qualifiers.get(JUNCTION_TAG, ())) for one in regions(product))
    assert judged == [
        ("AAGCTT junction", (BACKBONE,)),
        ("GCATGC junction", ("GFP amplicon",)),
        ("GFP amplicon insert", ()),
    ]
    # The one stretch judged is GFP, and not the 2.6 kb of vector on the other side of it.
    insert = next(one for one in regions(product) if one.name.endswith("insert"))
    assert len(product.extract(insert)) == len(gfp)


def test_a_linear_record_has_no_stretch_round_from_its_last_junction():
    record = SequenceRecord("A" * 100, features=(_junction(20, 30, "GFP"), _junction(60, 70, "x")))
    assert [one.name for one in regions(record)] == ["j", "GFP insert", "j"]
