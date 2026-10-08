"""The method's two figures: the assembly rows, and the oligo each PCR reads by its own pair."""

import random
from pathlib import Path

import pytest

from liulab_mbio.bench.pools import PrimerSite
from liulab_mbio.cloning.plan import PRODUCT_FILE
from liulab_mbio.edits import rotate
from liulab_synbio.igga.cargo import PRIMER_LENGTH, design_pool
from liulab_synbio.igga.figures import (
    OLIGO_FILE,
    SOURCE_KEY,
    assembly_rows,
    pool_pcr_figure,
)
from liulab_synbio.igga.method import ORTHOGONAL_SPLIT
from liulab_synbio.igga.project import Build
from liulab_synbio.igga.rounds import ROUND_FILE, assemble_rounds, representative

from .test_cargo import bases, block
from .test_rounds import POSITIONS, built, carrier, scheme

#: How much vector the test carrier spells either side of its cassette, which is room enough for
#: the margin a round's map draws around one.
FLANK = 800


@pytest.fixture(scope="module")
def made():
    return scheme()


@pytest.fixture(scope="module")
def parts(made):
    return representative(built(made)[1], POSITIONS)


@pytest.fixture(scope="module")
def rounds(made, parts):
    """Rounds on a vector with room either side of its cassette for a span to be drawn."""
    return assemble_rounds(carrier(made, flank=FLANK), parts, made, POSITIONS, name="library")


@pytest.fixture(scope="module")
def pooled():
    rng = random.Random(2)
    roles = [role for role, count in ORTHOGONAL_SPLIT for _ in range(count)]
    primers = tuple(
        PrimerSite(role, f"P{index}", bases(rng, PRIMER_LENGTH)) for index, role in enumerate(roles)
    )
    project = Build(
        "tiny",
        positions=("N",),
        parts=Path("parts.fasta"),
        vector=Path("vector.dna"),
        host="human",
        oligo_length=350,
        batch_size=2,
        completeness=0.99,
    )
    return design_pool([block("a", bases(rng, 200), 0)], project, primers=primers)


def test_a_round_draws_what_it_opens_and_what_it_makes(rounds):
    """The two rows are the files the plan writes, so the figure and the directory agree."""
    figure = assembly_rows(rounds, lit=2, at="../")

    assert figure.records == (
        "../" + ROUND_FILE.format(number=1),
        "../" + ROUND_FILE.format(number=2),
    )
    assert figure.linear
    assert figure.enzymes is not None
    assert set(figure.enzymes) == {
        rounds[0].scheme.internal_enzyme,
        rounds[0].scheme.external_enzyme,
    }
    assert figure.citation is not None
    assert figure.citation.source == SOURCE_KEY


def test_the_last_round_makes_the_product(rounds):
    figure = assembly_rows(rounds, lit=len(rounds))

    assert figure.records == (ROUND_FILE.format(number=len(rounds) - 1), PRODUCT_FILE)


def test_round_one_opens_the_vector_record_the_build_names(rounds):
    """Round 1 opens no round, so its first row is the vector the plan wrote."""
    figure = assembly_rows(rounds, lit=1, vector="block-vector-1.dna", at="../")

    assert figure.records == ("../block-vector-1.dna", "../" + ROUND_FILE.format(number=1))
    assert "what it opens above, what it makes below" in figure.caption


def test_round_one_draws_one_row_where_no_vector_record_is_written(rounds):
    figure = assembly_rows(rounds, lit=1)

    assert figure.records == (ROUND_FILE.format(number=1),)
    assert "above" not in figure.caption


def test_both_rows_are_drawn_over_the_cassette_and_not_the_backbone(rounds):
    """A span the two records share: what the round changed, and a margin either side."""
    one = rounds[1]
    figure = assembly_rows(rounds, lit=2)

    assert figure.span is not None
    start, end = figure.span
    assert start < one.entry.start
    assert end > one.scar.end
    assert end - start < len(one.destination)
    assert one.destination.fits(start, end)


def test_a_cassette_across_the_origin_draws_the_records_whole(made, parts):
    """The edit moves the origin of what it made, so no one span falls on both records."""
    turned = rotate(carrier(made, flank=FLANK), FLANK + 10)
    across = assemble_rounds(turned, parts, made, POSITIONS)

    assert assembly_rows(across, lit=1).span is None


def test_a_cassette_at_the_end_of_the_opened_record_draws_the_records_whole(made, parts):
    """A span past the opened record's last base runs across its origin and not the other's."""
    near = assemble_rounds(rotate(carrier(made, flank=FLANK), FLANK + 50), parts, made, POSITIONS)
    one = near[0]
    assert one.excised.end <= len(one.destination)
    assert one.scar.end > len(one.destination)

    assert assembly_rows(near, lit=1).span is None


def test_a_lit_round_lights_the_part_that_round_joined(rounds):
    """Which round the step is at is the parameter the published figure cannot offer."""
    figure = assembly_rows(rounds, lit=2)

    assert figure.highlight == (rounds[1].part.name,)
    assert rounds[1].part.name in figure.caption
    assert rounds[1].entry_overhang in figure.caption


@pytest.mark.parametrize("lit", [-1, 0, 99])
def test_a_round_this_build_does_not_run_is_refused(rounds, lit):
    with pytest.raises(ValueError, match="not one"):
        assembly_rows(rounds, lit=lit)


def test_an_assembly_figure_of_no_rounds_is_refused():
    with pytest.raises(ValueError, match="at least one round"):
        assembly_rows((), lit=1)


def test_each_pcr_lights_its_own_pair_of_the_three_roles(pooled):
    """PCR1 is the outer pair and PCR2 the inner, sharing one forward primer."""
    forward, inner, outer = pooled.pool.oligos[0].primers

    first = pool_pcr_figure(pooled, stage="PCR1")
    second = pool_pcr_figure(pooled, stage="PCR2")

    assert first.records == (OLIGO_FILE,)
    assert first.highlight == (forward, outer)
    assert second.highlight == (forward, inner)
    assert first.enzymes == (pooled.pool.layout.cutter.name,)
    assert "PCR1" in first.caption
    assert "batch" in first.caption
    assert "block" in second.caption
