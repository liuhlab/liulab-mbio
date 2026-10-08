"""The method's two figures: the assembly rows, and the oligo each PCR reads by its own pair."""

import random
from pathlib import Path

import pytest

from liulab_mbio.bench.pools import PrimerSite
from liulab_mbio.cloning.plan import PRODUCT_FILE
from liulab_synbio.igga.cargo import PRIMER_LENGTH, design_pool
from liulab_synbio.igga.figures import (
    OLIGO_FILE,
    SOURCE_KEY,
    assembly_rows,
    pool_pcr_figure,
)
from liulab_synbio.igga.method import ORTHOGONAL_SPLIT
from liulab_synbio.igga.project import Build
from liulab_synbio.igga.rounds import ROUND_FILE, assemble_rounds

from .test_cargo import bases, block
from .test_rounds import POSITIONS, built, carrier, scheme


@pytest.fixture(scope="module")
def made():
    return scheme()


@pytest.fixture(scope="module")
def rounds(made):
    from liulab_synbio.igga.rounds import representative

    parts = representative(built(made)[1], POSITIONS)
    return assemble_rounds(carrier(made), parts, made, POSITIONS, name="library")


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


def test_the_assembly_figure_names_the_record_every_round_writes(rounds):
    """The rows are the files the plan writes, so the figure and the directory agree."""
    figure = assembly_rows(rounds)

    assert figure.records == (
        *(ROUND_FILE.format(number=one.number) for one in rounds[:-1]),
        PRODUCT_FILE,
    )
    assert figure.linear
    assert figure.enzymes is not None
    assert set(figure.enzymes) == {
        rounds[0].scheme.internal_enzyme,
        rounds[0].scheme.external_enzyme,
    }
    assert figure.citation is not None
    assert figure.citation.source == SOURCE_KEY


def test_a_lit_round_lights_the_part_that_round_joined(rounds):
    """Which round is lit is the parameter the published figure cannot offer."""
    figure = assembly_rows(rounds, lit=2)

    assert figure.highlight == (rounds[1].part.name,)
    assert rounds[1].part.name in figure.caption
    assert rounds[1].entry_overhang in figure.caption


def test_no_round_lit_leaves_every_row_in_colour(rounds):
    figure = assembly_rows(rounds)
    assert figure.highlight == ()
    assert "in order" in figure.caption


@pytest.mark.parametrize("lit", [-1, 99])
def test_a_round_this_build_does_not_run_is_refused(rounds, lit):
    with pytest.raises(ValueError, match="not one"):
        assembly_rows(rounds, lit=lit)


def test_an_assembly_figure_of_no_rounds_is_refused():
    with pytest.raises(ValueError, match="at least one round"):
        assembly_rows(())


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
