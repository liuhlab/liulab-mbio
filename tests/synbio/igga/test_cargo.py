"""The oligo pool's own rules, on blocks and primer sheets small enough to read.

`test_ap1_demo.py` keeps what the pool does at AP-1's scale; here is what each function decides
alone: which primer a block takes, what a sheet or a band may say, and what is refused.
"""

import random
from dataclasses import replace
from decimal import Decimal
from pathlib import Path

import pytest

from liulab_mbio.bench.pools import PrimerSite
from liulab_mbio.sequence import Segment
from liulab_synbio.igga.cargo import (
    PRIMER_LENGTH,
    design_pool,
    read_bands,
    read_primers,
)
from liulab_synbio.igga.method import ORTHOGONAL_SPLIT
from liulab_synbio.igga.parts import Part
from liulab_synbio.igga.project import Project

INNER, FORWARD, OUTER = (count for _, count in ORTHOGONAL_SPLIT)
TOTAL = INNER + FORWARD + OUTER


def bases(rng, length):
    """Draw a block with no T, so it spells no site of an enzyme the method reserves."""
    return "".join(rng.choice("ACG") for _ in range(length))


def block(name, sequence, index=0):
    return Part(
        name,
        "N",
        index=index,
        sequence=sequence,
        barcode="",
        protein="",
        coding=Segment(0, len(sequence)),
    )


def sheet(tmp_path, rows, header="name\tsequence"):
    path = tmp_path / "primers.tsv"
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return path


def rows(count):
    return [f"p{index}\t{'ACG' * 7}" for index in range(count)]


@pytest.fixture(scope="module")
def project():
    return Project(
        "tiny",
        positions=("N",),
        parts=Path("parts.fasta"),
        vector=Path("vector.dna"),
        host="human",
        oligo_length=350,
        batch_size=2,
        coverage=1.0,
    )


@pytest.fixture(scope="module")
def primers():
    rng = random.Random(2)
    roles = [role for role, count in ORTHOGONAL_SPLIT for _ in range(count)]
    return tuple(
        PrimerSite(role, f"P{index}", bases(rng, PRIMER_LENGTH)) for index, role in enumerate(roles)
    )


@pytest.fixture(scope="module")
def parts():
    """A one-oligo block, a three-oligo block and a short one: three blocks, two batches."""
    rng = random.Random(2)
    return [
        block("a", bases(rng, 200), 0),
        block("b", bases(rng, 600), 1),
        block("c", bases(rng, 90), 2),
    ]


@pytest.fixture(scope="module")
def pooled(parts, project, primers):
    return design_pool(parts, project, primers=primers)


def test_a_block_becomes_one_named_oligo_a_fragment(pooled):
    assert [(one.name, one.source, one.fragment, one.fragments) for one in pooled.pool.oligos] == [
        ("a_f1", "a", 1, 1),
        ("b_f1", "b", 1, 3),
        ("b_f2", "b", 2, 3),
        ("b_f3", "b", 3, 3),
        ("c_f1", "c", 1, 1),
    ]
    assert [str(one.cargo.name) for one in pooled.splits] == ["a", "b", "c"]


def test_the_count_equals_the_floor_where_no_block_forces_one_more(pooled):
    assert pooled.floor == pooled.pool.count == 5
    assert pooled.over_floor == ()


def test_a_batch_shares_its_outer_and_forward_primers_and_a_gene_its_own_inner(pooled):
    used = {one.source: one.primers for one in pooled.pool.oligos}
    assert used["a"] == ("P96", "P0", "P131")
    assert used["b"] == ("P96", "P1", "P131")
    assert used["c"] == ("P97", "P0", "P132")
    assert [one.name for one in pooled.pool.primers] == ["P96", "P0", "P131", "P1", "P97", "P132"]


def test_no_internal_junction_spells_the_vector_s_two_overhangs(pooled, project):
    held = {project.scheme.entry_overhang, project.scheme.scar_overhang}
    for oligo in pooled.pool.oligos:
        assert not held & set(oligo.overhangs[1:-1])


def test_the_same_seed_writes_the_same_oligos_and_another_draws_other_filler(
    parts, project, primers, pooled
):
    again = design_pool(parts, project, primers=primers)
    other = design_pool(parts, project, primers=primers, seed=5)
    assert [one.sequence for one in again.pool.oligos] == [
        one.sequence for one in pooled.pool.oligos
    ]
    assert other.pool.oligos[0].sequence != pooled.pool.oligos[0].sequence


def test_the_bill_line_carries_its_key_and_the_headroom_of_the_bands_it_is_given(
    parts, project, primers, pooled
):
    assert pooled.item.headroom == ()
    banded = design_pool(
        parts,
        project,
        primers=primers,
        bands=read_bands({"count": ["1-4", "5-"], "length": ["301-350"]}),
        key="vendor-pool",
    )
    assert banded.item.key == "vendor-pool"
    assert {(one.quantity, one.band.low) for one in banded.item.headroom} == {
        ("count", 5),
        ("length", 301),
    }


def test_a_block_that_spells_no_legal_overhang_set_is_refused_and_says_so(project, primers):
    with pytest.raises(ValueError, match=r"block 'flat'.*no legal overhang set exists"):
        design_pool([block("flat", "A" * 600)], project, primers=primers)


def test_a_primer_set_short_of_the_allotment_is_refused_and_says_what_it_allots(
    parts, project, primers
):
    with pytest.raises(ValueError, match=rf"allots {TOTAL} primers .* {TOTAL - 1} were given"):
        design_pool(parts, project, primers=primers[:-1])


def test_an_empty_block_is_refused_by_name(parts, project, primers):
    hollow = Part("hollow", "N", index=3, sequence="", barcode="", protein="", coding=Segment(0, 1))
    with pytest.raises(ValueError, match="block 'hollow' is empty"):
        design_pool([*parts, hollow], project, primers=primers)


def test_a_block_carrying_the_synthesis_site_names_the_block_and_the_offset(project, primers):
    """A hand-built block spelling BsmBI on its bottom strand: the split refuses it."""
    carrier = block("carrier", bases(random.Random(3), 97) + "GAGACG" + bases(random.Random(4), 97))
    with pytest.raises(ValueError, match=r"block 'carrier'.*BsmBI site\(s\), the first at 97 "):
        design_pool([carrier], project, primers=primers, seed=7)


def test_a_batch_wider_than_the_plate_of_inner_primers_is_refused(parts, project, primers):
    wide = replace(project, batch_size=INNER + 1)
    with pytest.raises(ValueError, match=rf"batch_size is {INNER + 1}.*{INNER} .*one plate"):
        design_pool(parts, wide, primers=primers)


def test_a_sheet_allots_its_primers_to_the_roles_in_order(tmp_path):
    found = read_primers(sheet(tmp_path, rows(TOTAL)))
    assert [one.name for one in found] == [f"p{index}" for index in range(TOTAL)]
    roles = [one.role for one in found]
    assert [roles.count(role) for role, _ in ORTHOGONAL_SPLIT] == [INNER, FORWARD, OUTER]
    assert roles == sorted(roles, key=[role for role, _ in ORTHOGONAL_SPLIT].index)


def test_a_sheet_may_carry_a_byte_order_mark_blank_rows_and_primers_past_the_set(tmp_path):
    path = sheet(tmp_path, ["", *rows(TOTAL + 3)])
    path.write_text("﻿" + path.read_text(encoding="utf-8"), encoding="utf-8")
    assert len(read_primers(path)) == TOTAL


@pytest.mark.parametrize("header", ["sequence\tname", "name"])
def test_a_sheet_without_the_name_and_sequence_columns_is_not_a_primer_sheet(tmp_path, header):
    with pytest.raises(ValueError, match="is not a primer sheet"):
        read_primers(sheet(tmp_path, rows(TOTAL), header=header))


def test_an_empty_file_is_not_a_primer_sheet(tmp_path):
    path = tmp_path / "empty.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="is not a primer sheet"):
        read_primers(path)


def test_a_primer_with_no_sequence_is_named_by_its_line(tmp_path):
    with pytest.raises(ValueError, match=r"line 3 .* no sequence"):
        read_primers(sheet(tmp_path, ["p0\tACGT", "p1\t", *rows(TOTAL)]))


def test_a_sheet_short_of_the_set_says_how_many_it_holds(tmp_path):
    with pytest.raises(ValueError, match=rf"holds {TOTAL - 1} primers"):
        read_primers(sheet(tmp_path, rows(TOTAL - 1)))


def test_a_band_with_no_top_is_open_and_a_stated_one_holds_both_ends():
    bands = read_bands({"count": ["1-96", "97-"], "length": ["301-350"]})
    assert [one.high for one in bands["count"]] == [Decimal(96), None]
    assert bands["count"][1].holds(Decimal(10_000))
    assert bands["length"][0].holds(Decimal(301))
    assert not bands["length"][0].holds(Decimal(351))


def test_nothing_stated_reads_no_bands():
    assert read_bands({}) == {}


def test_a_quantity_the_pool_is_not_banded_by_is_refused():
    with pytest.raises(ValueError, match="not by 'price'"):
        read_bands({"price": ["1-2"]})


@pytest.mark.parametrize("written", ["350", "-350"])
def test_a_band_not_written_low_dash_high_is_refused(written):
    with pytest.raises(ValueError, match="not written 'low-high'"):
        read_bands({"length": [written]})
