"""Planning a whole library: the way in, the files it writes, and the vector-standard seam.

The method is `IGGA` and the project is written into a temporary directory, as a user writes one.
The part lists are three of two members, so a whole build stays small enough for the gate.
"""

import dataclasses
from typing import cast

import pytest

from liulab_mbio.checks import worst
from liulab_mbio.cloning.plan import PRODUCT_FILE, PROTOCOL_DATA_FILE, PROTOCOL_FILE
from liulab_mbio.protocol import read_protocol
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_mbio.sites import digest
from liulab_mbio.snapgene import write_dna
from liulab_mbio.translate import reverse_translate
from liulab_synbio.library.method import IGGA
from liulab_synbio.library.plan import (
    BARCODE_FILE,
    CHANGE_FILE,
    PARTS_FILE,
    Kind,
    plan_library,
    read_part_lists,
)
from liulab_synbio.library.project import Project
from liulab_synbio.library.rounds import ROUND_FILE

HOST = "e-coli-k12"
COVERAGE = 10.0

#: The positions one build fills, named so a record's own name says which it belongs to.
POSITIONS = ("N", "bZIP", "C")

#: Three part lists of two members, each name saying which position it fills.
LISTS = (
    {"N_a": "MKTAEK", "N_b": "MKTCEK"},
    {"bZIP_a": "WQAFAK", "bZIP_b": "WQAYAK"},
    {"C_a": "MKTHGK", "C_b": "MKTWGK"},
)


def pad(length: int) -> str:
    """`length` bases spelling no site any of the method's enzymes reads."""
    return ("TA" * length)[:length]


def fasta(lists=LISTS) -> str:
    """These part lists as one FASTA, each record named for the position it fills."""
    return "".join(
        f">{name} a description\n{sequence}\n" for one in lists for name, sequence in one.items()
    )


@pytest.fixture(scope="module")
def scheme():
    return IGGA


@pytest.fixture(scope="module")
def inputs(tmp_path_factory, scheme):
    """A parts FASTA, a vector already carrying a stuffer, and a bare one, all on disk."""
    out = tmp_path_factory.mktemp("inputs")
    (out / "parts.fasta").write_text(fasta(), encoding="utf-8")
    write_dna(
        SequenceRecord(
            pad(80) + scheme.internal_stuffer + pad(80), topology="circular", name="carrier"
        ),
        out / "carrier.dna",
    )
    write_dna(
        SequenceRecord(pad(200), topology="circular", name="bare"),
        out / "bare.dna",
    )
    return out


def project(inputs, *, vector: str = "carrier.dna") -> Project:
    """The project a test plans, naming the inputs written into `inputs`."""
    return Project(
        "library",
        positions=POSITIONS,
        parts=inputs / "parts.fasta",
        vector=inputs / vector,
        host=HOST,
        oligo_length=350,
        batch_size=96,
        coverage=COVERAGE,
    )


@pytest.fixture(scope="module")
def carrier(inputs):
    """The vector the plan was made from, as the plan read it."""
    from liulab_mbio.io import read_record

    return read_record(inputs / "carrier.dna")


@pytest.fixture(scope="module")
def plan(inputs):
    return plan_library(project(inputs), parts=LISTS)


@pytest.fixture(scope="module")
def written(plan, tmp_path_factory):
    return plan, plan.write(tmp_path_factory.mktemp("library"))


@pytest.fixture(scope="module")
def protocol(plan):
    return plan.protocol()


def test_one_call_plans_every_round_and_every_part(plan):
    assert len(plan.rounds) == len(POSITIONS) == 3
    assert len(plan.parts) == sum(len(one) for one in LISTS) == 6
    assert plan.constructs == 8
    assert [one.position for one in plan.rounds] == ["N", "bZIP", "C"]
    assert plan.status == "pass"


def test_write_puts_every_file_in_one_directory_and_the_protocol_reads_back(written):
    _plan, files = written

    names = {path.name for path in (files.parts, files.barcodes, files.changes)}
    assert names == {PARTS_FILE, BARCODE_FILE, CHANGE_FILE}
    assert [path.name for path in files.records] == [
        ROUND_FILE.format(number=1),
        ROUND_FILE.format(number=2),
        PRODUCT_FILE,
    ]
    assert files.protocol_data.name == PROTOCOL_DATA_FILE
    assert files.protocol.name == PROTOCOL_FILE
    for path in (files.parts, files.barcodes, files.changes, *files.records):
        assert path.read_bytes()
    assert len({path.parent for path in (files.parts, files.protocol, *files.records)}) == 1
    back = read_protocol(files.protocol_data)
    assert back.title
    assert back.steps
    assert files.protocol.read_text(encoding="utf-8")


def test_the_same_inputs_write_the_same_bytes(plan, tmp_path):
    once = plan.write(tmp_path / "once")
    twice = plan.write(tmp_path / "twice")

    for one, other in (
        (once.parts, twice.parts),
        (once.barcodes, twice.barcodes),
        (once.changes, twice.changes),
        (once.protocol_data, twice.protocol_data),
        (once.protocol, twice.protocol),
        *zip(once.records, twice.records, strict=True),
    ):
        assert one.read_bytes() == other.read_bytes(), one.name


def test_every_step_says_what_a_good_result_looks_like_and_carries_its_mixes(plan, protocol):
    tables = [table for step in protocol.steps for table in step.tables]
    programs = [program for step in protocol.steps for program in step.programs]

    for step in protocol.steps:
        assert step.expected, step.title
        assert step.instructions, step.title
    # Two digests and one ligation a round, and a program for each digest and for the growth.
    assert len(tables) == 3 * len(plan.rounds)
    assert len(programs) == 3 * len(plan.rounds)
    for program in programs:
        assert program.stages


def test_the_protocol_carries_the_traps_this_method_has(protocol):
    said = " ".join(note for step in protocol.steps for note in step.notes)
    confirm = " ".join(protocol.steps[-1].notes)

    assert "2 volumes here and 1 after the ligation" in said
    # The designed set is indel-aware, so the share is nil — and it is printed rather than implied,
    # because a set designed on mismatches alone leaves a share that is not.
    assert "0.0% of the single-base deletions" in confirm
    assert "keep the barcodes away from where a primer anneals" in confirm


def test_a_compatible_vector_pins_position_one_to_the_overhang_its_stuffer_spells(
    plan, scheme, carrier
):
    assert plan.destination.record.sequence == carrier.sequence
    assert plan.destination.edit is None
    assert plan.standard.entry_overhangs[0] == scheme.entry_overhang
    excised = [
        piece
        for piece in digest(carrier, scheme.internal)
        if (piece.left_overhang, piece.right_overhang)
        == (plan.standard.entry_overhangs[0], plan.standard.scar_overhang)
    ]
    assert len(excised) == 1


def test_a_retrofitted_vector_carries_the_overhang_the_standard_chose(scheme, inputs):
    made = plan_library(project(inputs, vector="bare.dna"), parts=LISTS, site=(100, 140))

    assert made.destination.edit is not None
    opened = [
        piece
        for piece in digest(made.destination.record, scheme.internal)
        if (piece.start, piece.end)
        == (made.destination.stuffer.start, made.destination.stuffer.end)
    ]
    assert len(opened) == 1
    assert opened[0].left_overhang == made.standard.entry_overhangs[0]
    assert opened[0].right_overhang == made.standard.scar_overhang
    assert made.status == "pass"


def test_the_plan_s_status_is_the_worst_over_every_round_s_checks(plan):
    named = [check.name for check in plan.checks]

    for one in plan.rounds:
        for check in one.checks:
            assert f"round {one.number} {check.name}" in named
    assert "round 3 reading frame" in named
    assert plan.status == worst(check.status for check in plan.checks)
    # A stuffer nothing excises is what a surviving site or a lost cut looks like to a round.
    broken = dataclasses.replace(
        plan,
        rounds=(*plan.rounds[:-1], dataclasses.replace(plan.rounds[-1], stuffer=Segment(0, 4))),
    )
    assert broken.status == "fail"
    assert next(one for one in broken.checks if one.name == "round 3 opens").status == "fail"


def test_dna_input_is_read_for_its_protein_and_kept_rather_than_re_coded(inputs, plan):
    # Coded for another host, so every codon differs from the one this host would have written.
    supplied = {
        name: reverse_translate(protein, host="human")
        for one in LISTS
        for name, protein in one.items()
    }
    lists = tuple({name: supplied[name] for name in one} for one in LISTS)

    made = plan_library(project(inputs), parts=lists, kind="dna")

    assert [one.protein for one in made.parts] == [one.protein for one in plan.parts]
    # The bases handed over are the ones kept; which codons survive is `design_parts`'s promise.
    assert [one.coding_sequence for one in made.parts] != [
        one.coding_sequence for one in plan.parts
    ]
    for part in made.parts:
        assert part.coding_sequence in supplied[part.name] or part.changes


def test_dna_that_is_not_a_coding_sequence_is_refused(inputs):
    lists = ({"N_a": "ATGAA", "N_b": "ATGAAA"}, LISTS[1], LISTS[2])

    with pytest.raises(ValueError, match="part 'N_a' is not a coding sequence"):
        plan_library(project(inputs), parts=lists, kind="dna")


def test_a_stop_inside_a_coded_part_is_refused(inputs):
    lists = ({"N_a": "ATGTAAAAA", "N_b": "ATGAAAAAA"}, LISTS[1], LISTS[2])

    with pytest.raises(ValueError, match="part 'N_a' spells a stop"):
        plan_library(project(inputs), parts=lists, kind="dna")


def test_a_fasta_is_sorted_into_one_part_list_a_position(tmp_path):
    path = tmp_path / "parts.fasta"
    path.write_text(fasta())

    lists = read_part_lists(path, POSITIONS)

    assert [sorted(one) for one in lists] == [sorted(one) for one in LISTS]


def test_a_name_that_says_no_position_is_refused_naming_it(tmp_path):
    path = tmp_path / "parts.fasta"
    path.write_text(">Q_zero\nMKTAEK\n")

    with pytest.raises(ValueError, match="the name 'Q_zero' says no position"):
        read_part_lists(path, POSITIONS)


def test_a_name_that_says_two_positions_is_refused_as_ambiguous(tmp_path):
    path = tmp_path / "parts.fasta"
    path.write_text(">N_C_both\nMKTAEK\n")

    with pytest.raises(ValueError, match="says 2 positions"):
        read_part_lists(path, POSITIONS)


def test_a_position_no_record_names_is_refused(tmp_path):
    path = tmp_path / "parts.fasta"
    path.write_text(">N_a\nMKTAEK\n>bZIP_a\nWQAFAK\n")

    with pytest.raises(ValueError, match="no record names position"):
        read_part_lists(path, POSITIONS)


def test_one_name_cannot_fill_two_part_lists(inputs):
    lists = (LISTS[0], LISTS[1], {"N_a": "MKTHGK", "C_b": "MKTWGK"})

    with pytest.raises(ValueError, match="names a part in two part lists"):
        plan_library(project(inputs), parts=lists)


def test_part_lists_must_be_one_a_position(inputs):
    with pytest.raises(ValueError, match="part list"):
        plan_library(project(inputs), parts=LISTS[:2])


def test_a_kind_that_is_neither_is_refused(inputs):
    # Cast deliberately: the annotation already forbids this, and the runtime guard is what a
    # caller reaching the function from the command line or from JSON actually meets.
    with pytest.raises(ValueError, match="kind is 'protein' or 'dna'"):
        plan_library(project(inputs), parts=LISTS, kind=cast(Kind, "rna"))


def test_each_round_is_sized_for_the_coverage_asked_for(plan):
    assert [row.products for row in plan.coverage] == [2, 4, 8]
    assert [row.colonies for row in plan.coverage] == [20, 40, 80]
    assert [one.coverage.colonies for one in plan.bench] == [20, 40, 80]
    assert [one.number for one in plan.bench] == [1, 2, 3]
