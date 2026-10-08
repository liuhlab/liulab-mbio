"""Dressing a split cargo as an orderable oligo pool."""

from decimal import Decimal

import pytest

from liulab_mbio.bench.pools import (
    OligoLayout,
    Pool,
    PrimerSite,
    build_oligo,
    headroom,
    oligo_record,
    pad_bases,
    pool_item,
    pool_sheet,
    primer_inventory,
)
from liulab_mbio.bench.prices import Band
from liulab_mbio.sequence import SequenceRecord, Strand, reverse_complement
from liulab_mbio.sites import find_sites
from liulab_mbio.split import split_cargo

CARGO = SequenceRecord(
    "".join(
        "ACGTTGCAAGCTTCGAGGATCCTAGCATGCTAGGCCATTAACGGTTACGATCGATTACGCAT"[i % 62] for i in range(700)
    ),
    name="demo",
)
FORWARD = (PrimerSite("P1 batch forward", "OP1", "AAACACGTGGCAAACATTCC"),)
REVERSE = (
    PrimerSite("P2 gene reverse", "OP2", "AAACCGGAGCCATACAGTAC"),
    PrimerSite("P3 batch outer", "OP3", "AAAGCACTCTTAGGCCTCTG"),
)


@pytest.fixture
def pool():
    layout = OligoLayout(350, enzyme="BsmBI")
    split = split_cargo(CARGO, layout.cutter, budget=layout.budget)
    oligos = tuple(
        build_oligo(
            split,
            one,
            name=f"demo_f{one.index + 1}",
            source="demo",
            layout=layout,
            forward=FORWARD,
            reverse=REVERSE,
            avoid=["BsaI", "BbsI"],
        )
        for one in split.fragments
    )
    return Pool("demo", layout=layout, oligos=oligos, primers=(*FORWARD, *REVERSE))


def test_the_layout_leaves_the_span_the_split_is_budgeted_by():
    layout = OligoLayout(350, enzyme="BsmBI")
    assert layout.overhead == 74
    assert layout.budget.span == 276


def test_every_oligo_is_one_length_so_the_pool_sits_inside_the_uniformity_band(pool):
    assert {len(one) for one in pool.oligos} == {350}
    assert pool.spread == 0.0


def test_an_oligo_carries_two_designed_cuts_and_no_other(pool):
    for one in pool.oligos:
        assert len(find_sites(SequenceRecord(one.sequence), "BsmBI")) == 2


def test_an_oligo_says_which_cargo_and_which_fragment_it_came_from(pool):
    assert [one.fragment for one in pool.oligos] == list(range(1, pool.count + 1))
    assert {one.source for one in pool.oligos} == {"demo"}
    assert all(one.fragments == pool.count for one in pool.oligos)
    rows = pool_sheet(pool).splitlines()
    assert rows[0].split("\t")[:4] == ["name", "sequence", "length", "source"]
    assert len(rows) == pool.count + 1


def test_the_inventory_names_every_primer_with_its_role_and_its_sequence(pool):
    rows = [line.split("\t") for line in primer_inventory(pool).splitlines()[1:]]
    assert [row[0] for row in rows] == ["OP1", "OP2", "OP3"]
    assert [row[1] for row in rows] == [one.role for one in (*FORWARD, *REVERSE)]
    assert [row[2] for row in rows] == [one.sequence for one in (*FORWARD, *REVERSE)]


def test_the_filler_spells_no_forbidden_site_in_the_context_it_sits_in():
    filler = pad_bases("GAGACG", "TTTTTTTT", 80, avoid=["BsaI", "BsmBI"], seed=3)
    whole = "GAGACG" + filler + "TTTTTTTT"
    assert len(filler) == 80
    assert not find_sites(SequenceRecord(whole[6:]), ["BsaI", "BsmBI"])


def test_the_filler_is_the_same_bases_on_a_second_run():
    assert pad_bases("AAAA", "TTTT", 40, avoid=["BsaI"], seed=7) == pad_bases(
        "AAAA", "TTTT", 40, avoid=["BsaI"], seed=7
    )


def test_a_quantity_at_the_top_of_its_band_says_it_has_no_slack():
    assert str(headroom("length", 350, [Band("length", Decimal(301), Decimal(350))])) == (
        "350 length, no slack above it at all"
    )
    assert headroom("length", 400, [Band("length", Decimal(301), Decimal(350))]) is None


def test_the_bill_line_carries_the_count_and_the_band_without_any_price(pool):
    item = pool_item(
        pool,
        key="oligo-pool",
        bands={"length": (Band("length", Decimal(301), Decimal(350)),)},
    )
    assert item.quantity == pool.count
    assert item.unit == "oligos"
    assert [str(one) for one in item.headroom] == ["350 length, no slack above it at all"]


def test_an_oligo_draws_as_a_record_with_every_primer_where_it_reads() -> None:
    """The three-primer figure is drawn over this: two nested pairs on one bar."""
    cargo = SequenceRecord("ATG" + "ACGTTGCA" * 8, name="block")
    layout = OligoLayout(160, enzyme="BsaI")
    sites = (
        PrimerSite("forward", "F1", "ACACGACGCTCTTCCGATCT"[:20]),
        PrimerSite("inner", "I1", "TGTGCTGCGAGAAGGCTAGA"),
        PrimerSite("outer", "O1", "GTGACTGGAGTTCAGACGTG"),
    )
    split = split_cargo(cargo, layout.cutter, budget=layout.budget)
    one = build_oligo(
        split,
        split.fragments[0],
        name="block_f1",
        source="block",
        layout=layout,
        forward=sites[:1],
        reverse=sites[1:],
    )

    record = oligo_record(one, layout=layout, primers=sites)

    assert record.sequence == one.sequence
    placed = {
        p.name: (p.binding_sites[0].start, p.binding_sites[0].end, p.binding_sites[0].strand)
        for p in record.primers
    }
    assert placed["F1"] == (0, 20, Strand.FORWARD)
    assert placed["I1"][2] is Strand.REVERSE
    assert placed["O1"][1] == len(one.sequence)
    # Each primer lies where the oligo spells it, forward or reverse-complemented.
    for name, (start, end, strand) in placed.items():
        spelt = next(s.sequence for s in sites if s.name == name)
        here = one.sequence[start:end]
        assert here == (spelt if strand is Strand.FORWARD else reverse_complement(spelt))
    [fragment] = [f for f in record.features if f.name.startswith("block fragment")]
    assert one.sequence[fragment.segments[0].start : fragment.segments[0].end].endswith(
        split.sequence(split.fragments[0])[-10:]
    )


def test_an_oligo_naming_a_primer_the_pool_does_not_carry_is_refused() -> None:
    cargo = SequenceRecord("ATG" + "ACGTTGCA" * 8, name="block")
    layout = OligoLayout(160, enzyme="BsaI")
    sites = (
        PrimerSite("forward", "F1", "ACACGACGCTCTTCCGATCT"),
        PrimerSite("inner", "I1", "TGTGCTGCGAGAAGGCTAGA"),
        PrimerSite("outer", "O1", "GTGACTGGAGTTCAGACGTG"),
    )
    split = split_cargo(cargo, layout.cutter, budget=layout.budget)
    one = build_oligo(
        split,
        split.fragments[0],
        name="block_f1",
        source="block",
        layout=layout,
        forward=sites[:1],
        reverse=sites[1:],
    )
    with pytest.raises(ValueError, match="which the pool does not carry"):
        oligo_record(one, layout=layout, primers=sites[:2])
