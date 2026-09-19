"""`recombine`, the simulation below a plan: the product, its pieces and its two junctions."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from liulab_mbio.cloning.gateway.att import CROSSOVER, REGIONS, find_att_sites
from liulab_mbio.cloning.gateway.recombination import recombine
from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)

from .records import att_site, destination_vector

if TYPE_CHECKING:
    from liulab_mbio.cloning.gateway.recombination import Recombination


@pytest.fixture(scope="module")
def lr(entry: SequenceRecord, destination: SequenceRecord) -> Recombination:
    """The LR reaction a plan runs: that entry clone into that destination vector."""
    return recombine(entry, destination, reaction="LR")


@pytest.fixture(scope="module")
def bp(insert: SequenceRecord, donor: SequenceRecord) -> Recombination:
    """The BP reaction before it: that attB insert into that donor vector."""
    return recombine(insert, donor, reaction="BP")


def test_each_junction_spells_the_att_site_the_arithmetic_names(lr: Recombination) -> None:
    product = lr.product
    first, second = lr.junctions

    assert (first.name, second.name) == ("attB1", "attB2")
    assert first.bases == REGIONS["attB1"]
    assert second.bases == REGIONS["attB2"]
    assert product.extract(first.span) == REGIONS["attB1"]
    assert reverse_complement(product.extract(second.span)) == REGIONS["attB2"]


def test_the_product_carries_the_insert_between_the_two_att_sites(
    lr: Recombination, gfp: SequenceRecord
) -> None:
    product = lr.product
    first, second = lr.junctions

    assert product.topology == "circular"
    assert product.extract(Segment(first.end, second.start)) == gfp.sequence
    assert [one.name for one in find_att_sites(product)] == ["attB1", "attB2"]


def test_the_expression_clone_is_the_backbone_and_the_insert_and_nothing_else(
    lr: Recombination, destination: SequenceRecord
) -> None:
    dropped = destination.extract(Segment(*lr.cassette))

    assert len(lr.product) == lr.backbone.length + lr.moved.length
    assert lr.backbone.length + len(dropped) == len(destination)
    assert dropped not in lr.product.sequence
    assert "ccdB" not in {feature.name for feature in lr.product.features}


def test_every_feature_of_both_records_is_carried_to_its_new_coordinates(
    lr: Recombination,
) -> None:
    product = lr.product
    named = {feature.name: feature for feature in product.features}

    assert set(named) == {"AmpR", "T7 promoter", "RBS", "GFP", "attB1", "attB2"}
    assert product.extract(named["GFP"]).startswith("ATG")
    starts = [feature.segments[0].start for feature in product.features]
    assert starts == sorted(starts)


def test_a_primer_of_the_source_records_is_annotated_where_it_still_binds(
    entry: SequenceRecord, destination: SequenceRecord, gfp: SequenceRecord
) -> None:
    at = entry.sequence.find(gfp.sequence)
    kept = Primer(
        "reads the insert",
        gfp.sequence[:20],
        binding_sites=(BindingSite(at, at + 20, Strand.FORWARD),),
    )
    dropped = Primer(
        "reads the backbone",
        entry.sequence[:20],
        binding_sites=(BindingSite(0, 20, Strand.FORWARD),),
    )
    carried = recombine(
        SequenceRecord(
            entry.sequence,
            topology="circular",
            name=entry.name,
            features=entry.features,
            primers=(kept, dropped),
        ),
        destination,
        reaction="LR",
    ).product

    annotated = {primer.name: primer for primer in carried.primers}
    assert set(annotated) == {"reads the insert"}
    site = annotated["reads the insert"].binding_sites[0]
    assert carried.sequence[site.start : site.end] == kept.sequence


def test_an_entry_clone_read_on_the_other_strand_makes_the_same_product(
    entry: SequenceRecord, destination: SequenceRecord, lr: Recombination
) -> None:
    turned = SequenceRecord(
        reverse_complement(entry.sequence), topology="circular", name=entry.name
    )

    made = recombine(turned, destination, reaction="LR")

    assert made.product.sequence == lr.product.sequence


def test_a_destination_vector_with_a_drifted_flank_is_still_recognised(
    entry: SequenceRecord, gfp: SequenceRecord
) -> None:
    drifted = destination_vector(name="pDEST-drifted")
    at = drifted.sequence.find(att_site("attR1"))
    changed = SequenceRecord(
        drifted.sequence[:at] + "T" + drifted.sequence[at + 1 :],
        topology="circular",
        name=drifted.name,
        features=drifted.features,
    )

    made = recombine(entry, changed, reaction="LR")

    assert made.junctions[0].bases == "T" + REGIONS["attB1"][1:]
    assert made.status == "pass"
    assert gfp.sequence in made.product.sequence


def test_the_segment_that_moved_is_bounded_inside_each_att_region_not_at_its_edge(
    lr: Recombination,
) -> None:
    first, second = lr.junctions
    start, end = lr.boundaries

    assert (start, end) == (first.start + CROSSOVER, second.end - CROSSOVER)
    assert lr.product.extract(Segment(start, end)) == lr.moved.bases


def test_bp_writes_an_entry_clone_carrying_the_insert_between_its_attl_sites(
    bp: Recombination, gfp: SequenceRecord
) -> None:
    first, second = bp.junctions

    assert (first.name, second.name) == ("attL1", "attL2")
    assert first.bases == REGIONS["attL1"]
    assert second.bases == REGIONS["attL2"]
    assert bp.product.topology == "circular"
    assert bp.product.extract(Segment(first.end, second.start)) == gfp.sequence
    assert [one.name for one in find_att_sites(bp.product)] == ["attL1", "attL2"]
    assert "GGGG" + REGIONS["attB1"] not in bp.product.sequence
