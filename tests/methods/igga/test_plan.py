"""Planning a whole library: the way in, the files it writes, and the vector-standard seam.

The method is `IGGA` and the build is written into a temporary directory, as a user writes one.
The part lists are three of two members, so a whole build stays small enough for the gate.
"""

import dataclasses
import math
from pathlib import Path
from typing import cast

import pytest

from mbio.bench.amounts import DNA_VOLUME_UL
from mbio.bench.prices import read_prices
from mbio.checks import worst
from mbio.cloning.plan import PRODUCT_FILE
from mbio.io import read_record
from mbio.protocol.model import read_project
from mbio.protocol.render import INDEX_FILE, PROJECT_DATA_FILE
from mbio.sequence import Feature, Segment, SequenceRecord
from mbio.sites import digest, find_sites
from mbio.snapgene import write_dna
from mbio.translate import reverse_translate
from synbio.igga.bench import (
    DIGEST_NG,
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    pool_floor_ng_ul,
)
from synbio.igga.gate import check_product
from synbio.igga.method import IGGA
from synbio.igga.plan import (
    BARCODE_FILE,
    BLOCK_VECTOR_FILE,
    CHANGE_FILE,
    PARTS_FILE,
    WORKING_VECTOR_FILE,
    Kind,
    plan_igga,
    read_part_lists,
)
from synbio.igga.project import Build, FinalAssembly
from synbio.igga.protocols import ASSEMBLY, CREATION, FINAL, ORDERING
from synbio.igga.protocols.run import ROUND_EQUIPMENT
from synbio.igga.rounds import ROUND_FILE
from synbio.igga.vector import cargo_enzyme

from ...chains import whole

HOST = "e-coli-k12"
COMPLETENESS = 0.99

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


def build(
    inputs, *, vector: str = "carrier.dna", working: str | None = None, primers: Path | None = None
) -> Build:
    """The build a test plans, naming the inputs written into `inputs`."""
    return Build(
        "library",
        positions=POSITIONS,
        parts=inputs / "parts.fasta",
        vector=inputs / vector,
        working_vector=None if working is None else inputs / working,
        host=HOST,
        oligo_length=350,
        batch_size=96,
        completeness=COMPLETENESS,
        primers=primers,
    )


@pytest.fixture(scope="module")
def carrier(inputs):
    """The vector the plan was made from, as the plan read it."""
    from mbio.io import read_record

    return read_record(inputs / "carrier.dna")


@pytest.fixture(scope="module")
def plan(inputs):
    return plan_igga(build(inputs), parts=LISTS)


@pytest.fixture(scope="module")
def written(plan, tmp_path_factory):
    return plan, plan.write(tmp_path_factory.mktemp("library"))


#: The orthogonal primer set the demo's pool is amplified by, which is the one set this repo
#: ships. A pool is what makes the blocks need a vector to supply their stuffers.
PRIMERS = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library" / "primers.tsv"


@pytest.fixture(scope="module")
def pooled(inputs, tmp_path_factory):
    """The same build with a pool designed, and the files it writes."""
    made = plan_igga(build(inputs, primers=PRIMERS), parts=LISTS)
    return made, made.write(tmp_path_factory.mktemp("pooled"))


@pytest.fixture(scope="module")
def protocol(plan):
    """Every protocol of the run as one, which is what most of these tests ask about."""
    return whole(plan.chain())


def test_a_plain_run_orders_its_blocks_assembles_and_moves_the_library(plan):
    assert [one.title for one in plan.chain().protocols] == [ORDERING, ASSEMBLY, FINAL]


def test_a_pool_splits_ordering_from_making_the_cargo(pooled):
    made, _files = pooled

    # This build states no primer plates, so the pool's primers are ordered with the pool.
    assert [one.title for one in made.chain().protocols] == [ORDERING, CREATION, ASSEMBLY, FINAL]


def test_a_run_without_a_pool_takes_its_cargo_from_the_vendors_tube(plan, pooled):
    """A run that archives nothing takes the blocks as the vendor shipped them."""
    made, _files = pooled
    plain, with_pool = plan.chain(), made.chain()
    rounds = [
        next(one for one in chain.protocols if one.title == ASSEMBLY)
        for chain in (plain, with_pool)
    ]

    assert [one.name for one in rounds[0].consumes] == ["synthesised blocks", "block vector 1"]
    assert "cargo archive plate" in {one.name for one in rounds[1].consumes}
    for chain in (plain, with_pool):
        assert chain.audit()[0].status == "pass", chain.audit()[0].detail


def test_each_protocol_of_the_chain_answers_for_its_own_page(plan):
    run = plan.chain()
    ordering, assembly = run.protocols[0], run.protocols[1]

    assert ordering.summary == "Order every block this library is built from."
    # Only the protocols running a round's chemistry buy the reagents the rounds share.
    assert not ordering.equipment
    assert assembly.equipment[: len(ROUND_EQUIPMENT)] == ROUND_EQUIPMENT


def test_each_page_names_the_sequencer_its_own_reads_ask_for(plan):
    """The instruments follow from the reads a page runs, never from a sentence beside them."""
    chain = plan.chain()
    assembly, final = chain.protocols[1], chain.protocols[2]
    reads = plan.reads

    assert reads is not None
    taken = {
        assembly: (reads.linkage, reads.representation),
        final: (reads.final_representation,),
    }
    for page, pairs in taken.items():
        pairs = tuple(one for one in pairs if one is not None)
        listed = [one for one in page.equipment if "sequencer" in one]
        assert len(listed) == len({one.platform for one in pairs}), page.title
        for one in pairs:
            wanted = f"{one.platform.replace(' ', '-').capitalize()} sequencer"
            assert one.name in next(item for item in listed if item.startswith(wanted))

    # The page that never runs the linkage read does not ask for a long-read instrument.
    assert not [one for one in final.equipment if one.startswith("Long-read")]


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
    assert files.protocol_data.name == PROJECT_DATA_FILE
    assert files.protocol[0].name == INDEX_FILE
    for path in (files.parts, files.barcodes, files.changes, *files.records):
        assert path.read_bytes()
    assert len({path.parent for path in (files.parts, *files.records)}) == 1
    assert {path.parent for path in files.protocol} == {files.protocol_data.parent}
    back = read_project(files.protocol_data)
    assert back.title
    assert [one.title for one in back.protocols]
    assert all(path.read_text(encoding="utf-8") for path in files.protocol)


def test_a_pool_writes_the_block_vector_each_position_closes_into(pooled):
    """Only the cargo is synthesised, so a vector has to supply the stuffers either side of it."""
    made, files = pooled

    assert [path.name for path in files.block_vectors] == [
        BLOCK_VECTOR_FILE.format(number=number) for number in (1, 2, 3)
    ]
    assert [one.record.name for one in made.block_vectors] == [
        f"carrier {position}" for position in POSITIONS
    ]
    for one, path in zip(made.block_vectors, files.block_vectors, strict=True):
        assert read_record(path).sequence == one.record.sequence
    assert set(files.block_vectors) <= set(files.paths)


def test_a_build_naming_no_primer_set_writes_no_block_vector(written):
    """A block ordered whole carries its own stuffers, so nothing has to supply them."""
    made, files = written

    assert made.pool is None
    assert made.block_vectors == ()
    assert files.block_vectors == ()


def test_the_same_inputs_write_the_same_bytes(plan, tmp_path):
    once = plan.write(tmp_path / "once")
    twice = plan.write(tmp_path / "twice")

    for one, other in (
        (once.parts, twice.parts),
        (once.barcodes, twice.barcodes),
        (once.changes, twice.changes),
        (once.protocol_data, twice.protocol_data),
        *zip(once.protocol, twice.protocol, strict=True),
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
    # The final assembly adds one more growth; without a working vector it names no other.
    assert len(tables) == 3 * len(plan.rounds)
    assert len(programs) == 3 * len(plan.rounds) + 1
    for program in programs:
        assert program.stages


def test_the_protocol_carries_the_traps_this_method_has(protocol):
    said = " ".join(n.text for step in protocol.steps for n in step.noted)
    linkage, representation = (
        " ".join(n.text for n in step.noted) for step in protocol.steps[-7:-5]
    )

    assert "2 volumes here and 1 after the ligation" in said
    # The designed set is indel-aware, so the share is nil — and it is printed rather than implied,
    # because a set designed on mismatches alone leaves a share that is not. Exactly none reads
    # as none; a share merely rounding to 0.0% still prints the figure.
    assert "None of the single-base deletions" in linkage
    assert "Keep the barcodes away from where a primer anneals" in representation


def test_the_finished_library_is_read_for_linkage_and_for_representation(protocol):
    """The method's last two steps, and its rule that only one of them repeats."""
    linkage, representation = protocol.steps[-7:-5]

    assert [linkage.title, representation.title] == ["Read linkage", "Read representation"]
    # Linkage keeps H28 because this build states no mark of its own: no source sets one for
    # barcode-to-part fidelity. Representation is held to Joung's bar, so it carries none.
    assert [hole.id for hole in linkage.holes] == ["H28"]
    assert representation.holes == ()
    said = " ".join(n.text for step in protocol.steps for n in step.noted)
    assert said.count("after every later bottleneck") == 1


def test_a_build_may_state_the_four_numbers_the_method_leaves_open(pooled):
    """Stated, each prints as this run's own and its hole goes; stated nowhere, the hole stands."""
    made, _files = pooled
    stated = dataclasses.replace(
        made,
        build=dataclasses.replace(
            made.build,
            linkage_fidelity=0.9,
            final_assembly=FinalAssembly(75.0),
            pcr1_cycles=16,
            pcr2_cycles=18,
        ),
    )
    before, after = whole(made.chain()), whole(stated.chain())
    said = " ".join(
        text
        for step in after.steps
        for text in (*step.instructions, *(n.text for n in step.noted), *step.expected)
    )

    assert [hole.id for step in before.steps for hole in step.holes] == [
        "H29",
        "H30",
        "H28",
        "H31",
        "H32",
        "H31",
        "H24",
    ]
    # What is left is what no number closes: this build names no working vector, and the
    # carrier it ran in presents no cut that frees the cargo.
    assert [hole.id for step in after.steps for hole in step.holes] == ["H31", "H32", "H31"]
    assert "16 cycles is what this run measured for its own polymerase" in said
    assert "18 cycles is what this run measured" in said
    assert "75 ng of working vector" in said
    assert "75 ng of vector is what this run measured" in said
    assert "this run's own mark" in said


def test_both_read_steps_name_their_pair_and_its_amplicon(plan, protocol):
    """The plan designs the pairs, so neither step sends anyone to the bench without one."""
    linkage, representation = protocol.steps[-7:-5]
    pairs = plan.reads

    for step, pair in ((linkage, pairs.linkage), (representation, pairs.representation)):
        said = " ".join(step.instructions)
        # The oligos are named, never spelled out: the sheet is where a sequence belongs.
        assert pair.forward.name in said
        assert pair.reverse.name in said
        assert pair.forward.sequence not in said
        assert f"{pair.amplicon_length} bp amplicon" in said
        # And where to find the pair, which is the sheet the plan wrote it to.
        assert "library-read-primers.tsv" in said
    assert "long read" in " ".join(linkage.instructions)


def test_every_sheet_the_run_writes_is_one_its_pages_can_link(plan, pooled):
    """A page links a filename only where the run states the file, so the run states them all."""
    run = plan.chain()
    assert run.files == (
        "../parts.tsv",
        "../barcodes.tsv",
        "../changes.tsv",
        "../library-read-primers.tsv",
    )
    assert [one.files for one in run.protocols] == [run.files] * len(run.protocols)
    # A build with a pool writes the two pool sheets and a block vector a position besides.
    made, _files = pooled
    assert set(made.chain().files) - set(run.files) == {
        "../pool.tsv",
        "../pool-primers.tsv",
        "../block-vector-1.dna",
        "../block-vector-2.dna",
        "../block-vector-3.dna",
    }


def test_the_representation_step_states_the_marks_and_the_depth_they_take(protocol):
    """Joung's three, and the read depth that follows from the library's own width."""
    said = " ".join(protocol.steps[-6].expected)

    assert "99.5%" in said
    assert "skew ratio below 10" in said
    assert "100 or more reads a member" in said


def test_no_emitted_protocol_carries_the_read_primer_hole(protocol):
    """H27 is closed: the plan designs both pairs rather than naming none."""
    assert "H27" not in [hole.id for step in protocol.steps for hole in step.holes]


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
    made = plan_igga(build(inputs, vector="bare.dna"), parts=LISTS, site=(100, 140))

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


#: A coding sequence spelling PaqCI's site, which is the enzyme the bare vector leaves free.
SPELLS_CARGO = "ATGCACCTGCAAGAAAAA"


def coded() -> tuple[dict[str, str], ...]:
    """The part lists as DNA, the first member spelling the cargo enzyme's site."""
    first, *rest = LISTS
    return (
        {
            name: SPELLS_CARGO if name == "N_a" else reverse_translate(protein, host="human")
            for name, protein in first.items()
        },
        *(
            {name: reverse_translate(protein, host="human") for name, protein in one.items()}
            for one in rest
        ),
    )


def test_a_working_vector_fixes_the_cargo_enzyme_before_a_block_is_designed(inputs):
    """The vector is an input and the blocks are an output, so the pipeline does the ordering."""
    lists = coded()
    bare = SequenceRecord(pad(200), topology="circular", name="bare")
    alone = cargo_enzyme([bare], scheme=IGGA).enzyme
    assert alone is not None
    assert find_sites(SequenceRecord(lists[0]["N_a"]), alone)

    made = plan_igga(
        build(inputs, working="bare.dna"), parts=lists, kind="dna", working_site=(100, 140)
    )

    assert made.working is not None
    assert made.working.enzyme == alone
    # Reserved from the vector alone, so the one part that spelled it gives the site up.
    assert not any(find_sites(SequenceRecord(one.sequence), alone) for one in made.parts)
    changed = next(one for one in made.parts if one.name == "N_a")
    assert [one.site.enzyme for one in changed.changes] == [alone]


def test_a_working_vector_is_written_out_with_its_cassette_in(inputs, tmp_path):
    """The backbone the build names does not open; the record the plan writes does."""
    made = plan_igga(
        build(inputs, working="bare.dna"), parts=coded(), kind="dna", working_site=(100, 140)
    )

    files = made.write(tmp_path)

    assert made.working is not None
    assert files.working_vector is not None
    assert files.working_vector.name == WORKING_VECTOR_FILE
    assert read_record(files.working_vector).sequence == made.working.record.sequence
    assert files.working_vector in files.paths


def test_a_build_naming_no_working_vector_reserves_nothing_of_its_own(inputs):
    """Nothing is added to the reserved set where there is no working vector to read one off."""
    lists = coded()
    alone = cargo_enzyme([SequenceRecord(pad(200))], scheme=IGGA).enzyme
    assert alone is not None

    made = plan_igga(build(inputs), parts=lists, kind="dna")

    assert made.working is None
    assert any(find_sites(SequenceRecord(one.sequence), alone) for one in made.parts)


def test_the_plan_is_judged_by_the_gate_and_by_nothing_of_its_own(plan, inputs):
    named = {check.name for check in plan.checks}

    assert plan.checks == plan.verdict.checks
    assert plan.status == worst(check.status for check in plan.checks) == "pass"
    assert {"destination opens", "donor releases", "cargo frame", "terminal stop"} <= named
    # The protocol prints one badge a check name, not one a molecule judged.
    assert {check.name for check in plan.verdict.summary} == named
    assert len(plan.verdict.summary) < len(plan.checks)
    # A product the internal enzyme no longer opens is what a lost cut looks like to the gate.
    site = IGGA.internal.site
    broken = SequenceRecord(
        str(plan.product.sequence).replace(site, "A" * len(site), 1), topology="circular"
    )
    judged = check_product(
        broken,
        build=plan.build,
        barcodes={one.position: [one.barcode] for one in plan.representative_parts},
    )
    opens = next(one for one in judged if one.name == "product opens")
    assert opens.status == "fail"
    assert opens.check.value == 1


def test_dna_input_is_read_for_its_protein_and_kept_rather_than_re_coded(inputs, plan):
    # Coded for another host, so every codon differs from the one this host would have written.
    supplied = {
        name: reverse_translate(protein, host="human")
        for one in LISTS
        for name, protein in one.items()
    }
    lists = tuple({name: supplied[name] for name in one} for one in LISTS)

    made = plan_igga(build(inputs), parts=lists, kind="dna")

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
        plan_igga(build(inputs), parts=lists, kind="dna")


def test_a_stop_inside_a_coded_part_is_refused(inputs):
    lists = ({"N_a": "ATGTAAAAA", "N_b": "ATGAAAAAA"}, LISTS[1], LISTS[2])

    with pytest.raises(ValueError, match="part 'N_a' spells a stop"):
        plan_igga(build(inputs), parts=lists, kind="dna")


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
        plan_igga(build(inputs), parts=lists)


def test_part_lists_must_be_one_a_position(inputs):
    with pytest.raises(ValueError, match="part list"):
        plan_igga(build(inputs), parts=LISTS[:2])


def test_a_kind_that_is_neither_is_refused(inputs):
    # Cast deliberately: the annotation already forbids this, and the runtime guard is what a
    # caller reaching the function from the command line or from JSON actually meets.
    with pytest.raises(ValueError, match="kind is 'protein' or 'dna'"):
        plan_igga(build(inputs), parts=LISTS, kind=cast(Kind, "rna"))


def test_each_round_is_sized_for_the_completeness_asked_for(plan):
    assert [row.products for row in plan.coverage] == [2, 4, 8]
    assert [row.colonies for row in plan.coverage] == [8, 21, 51]
    assert [one.coverage.colonies for one in plan.bench] == [8, 21, 51]
    assert [one.number for one in plan.bench] == [1, 2, 3]


def test_the_protocol_reads_the_colony_count_as_a_floor_and_not_a_multiple(protocol):
    """Every page the count reaches states the chance asked for; the multiple is what it costs."""
    assert (
        protocol.overview["Completeness"]
        == "0.99 chance nothing is missing, 51 colonies at the end"
    )
    assert "51 colonies for a 0.99 chance that none of its 8 products is missing" in " ".join(
        protocol.highlights
    )
    growth = next(one for one in protocol.steps if one.title.startswith("Round 3: recover"))
    assert "At least 51 net colonies" in growth.expected[0]
    assert "floor for the 0.99 chance this design asked for" in growth.expected[0]
    assert "it works out at 6x this round's products" in " ".join(n.text for n in growth.noted)
    final = next(one for one in protocol.steps if one.title.startswith("Clean the assembly up"))
    assert (
        "At least 51 net colonies: the floor for the 0.99 chance this design asked for"
        in final.expected[0]
    )
    # The count the reader works out is read against the same floor, on both steps.
    assert [one.calculator and one.calculator.floor for one in (growth, final)] == [51, 51]
    assert growth.calculator is not None
    assert growth.calculator.counted == "round 3 titre"


#: A price record as a user writes one: the synthesis order banded by count and by length, and
#: nothing else priced.
PRICES = """key,item,bands,charge,basis,currency
synthesised blocks,gene fragments,count 1-100; length_nt 1-2000,1200.00,per order,USD
"""


def test_the_bill_computes_its_quantities_and_holes_the_money_with_no_record(plan, protocol):
    bill = protocol.bill

    assert bill is not None
    blocks = bill.rows[0]
    assert (blocks.item, blocks.quantity) == ("Synthesised blocks", len(plan.parts))
    assert blocks.charge == ""
    assert blocks.hole is not None
    assert blocks.hole.kind == "price"
    assert bill.total == ""


@pytest.fixture(scope="module")
def priced(plan, tmp_path_factory):
    """The same run with a price record loaded, which is what makes its bill cite one."""
    record = tmp_path_factory.mktemp("prices") / "prices.csv"
    record.write_text(PRICES, encoding="utf-8")
    return dataclasses.replace(plan, prices=read_prices(record)).chain()


def test_a_price_record_prices_the_bill_and_reports_its_headroom(priced):
    bill = priced.bill

    assert bill is not None
    blocks = bill.rows[0]
    assert blocks.charge == "1200.00"
    assert "94 below the next band" in blocks.headroom
    assert blocks.citation is not None
    assert (bill.currency, bill.total) == ("USD", "1200.00")
    # Everything else the run buys is still a hole, and no figure is estimated for one.
    assert all(row.hole is not None for row in bill.rows[1:])


def test_a_run_holding_a_price_record_names_it_as_a_source_of_the_run(priced):
    """The bill is the run's, so the record it cites is named there and on no protocol."""
    assert priced.bill is not None
    assert priced.bill.cited == frozenset({"prices"})
    assert "prices" in priced.sources
    assert not [one for one in priced.protocols if "prices" in one.sources]


def said_by(protocol) -> str:
    """Every line of the protocol that could name a drug."""
    return " ".join(
        [one.note or "" for one in protocol.materials]
        + [line for step in protocol.steps for line in step.instructions]
    )


def test_a_vector_naming_no_marker_keeps_the_hole_rather_than_naming_a_drug(protocol):
    """The toy vector annotates none, so nothing invents the published carbenicillin for it."""
    assert "carbenicillin" not in said_by(protocol)
    assert "H22" in [one.id for one in protocol.holes]


def test_a_vector_that_names_its_marker_plates_every_round_on_it(plan, carrier):
    """The marker is a fact in the record, so the drug follows from it and the hole is answered."""
    kanr = dataclasses.replace(carrier, features=(Feature("KanR", "CDS", (Segment(10, 100),)),))

    protocol = whole(dataclasses.replace(plan, vector=kanr).chain())

    said = said_by(protocol)
    assert "LB with 50 µg/mL kanamycin, the destination vector's own marker (KanR)" in said
    assert "Plate a measured dilution of the recovery on 50 µg/mL kanamycin" in said
    assert "H22" not in [one.id for one in protocol.holes]


def test_the_pooling_step_states_the_mass_the_round_then_digests(plan, protocol):
    """One number, said twice: a pool sized against anything else would not fit the digest."""
    step = next(one for one in protocol.steps if one.title == "Pool each part list")
    said = " ".join((*step.instructions, *step.expected))

    for row in plan.bench:
        assert f"{row.donor_digest.nanograms:,.0f} ng" in said


def test_the_pooling_floor_is_computed_from_the_digest_and_never_typed(protocol):
    """The floor is what the digest leaves room for, so the two can never disagree."""
    step = next(one for one in protocol.steps if one.title == "Pool each part list")

    limit = DIGEST_NG / (DIGEST_VOLUME_UL - 2 * ENZYME_UL)
    floor = pool_floor_ng_ul()

    # Rounded up to the tenth a page states, which is also what `fits` accepts: at the bare
    # limit the DNA fills the tube exactly and nothing is left to make the volume up with.
    assert floor == math.ceil(limit * 10) / 10
    assert floor > limit
    assert f"at least {floor:g} ng/µL" in " ".join(step.instructions)


def test_the_digest_table_pipettes_the_volume_the_stated_floor_implies(plan, protocol):
    """The floor a page states and the volume its digest table asks for are one number."""
    floor = pool_floor_ng_ul()
    step = next(one for one in protocol.steps if one.key == "round-1-release")
    (table,) = step.tables
    pool = next(one for one in table.components if not one.master_mix)

    # The stand-in volume an unmeasured DNA carries would say the pool is at 1,000 ng/µL.
    assert pool.volume_ul != DNA_VOLUME_UL
    assert pool.volume_ul == round(plan.bench[0].donor_digest.nanograms / floor, 2)
    assert sum(one.volume_ul for one in table.components) == DIGEST_VOLUME_UL


def test_the_two_stuffer_figures_are_one_stuffer_counted_two_ways(plan, protocol, scheme):
    """30 and 34 are the whole stuffer and its excised core, and no page leaves them at odds."""
    kept = plan.rounds[-1].stuffer
    whole_stuffer = len(scheme.internal_stuffer)
    core = len(scheme.internal_stuffer_core)

    opened = next(one for one in protocol.steps if one.key == "round-1-open")
    read = next(one for one in protocol.steps if one.key == "read-representation")
    said = " ".join(read.instructions)

    # Each round excises the core; what the last round leaves is the span the read anchors in.
    assert plan.rounds[0].excised.length == core
    assert f"The {core} bp internal stuffer comes out" in " ".join(opened.expected)
    assert kept.end - kept.start == core
    assert f"{core} bp internal stuffer every member keeps" in said
    assert f"{whole_stuffer} bp internal stuffer" not in said

    # The method's own provenance is where the two numbers meet, so it carries both.
    assert f"{whole_stuffer}-base internal stuffer" in scheme.source
    assert f"{core}-base core" in scheme.source
