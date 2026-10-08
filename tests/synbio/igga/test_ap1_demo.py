"""The AP-1 demo planned end to end, from the amino acid sequences and nothing else.

The inputs are `docs/examples/ap1-library`, specified by `docs/research/ap1-demo-project.md`:
72 proteins, a project file and a destination. Every assertion here is on the planned product and
the blocks that make it, not on how the design reached them.
"""

import re
from collections import Counter
from dataclasses import replace
from itertools import groupby
from pathlib import Path

import pytest

from liulab_mbio.barcodes import MAX_HOMOPOLYMER
from liulab_mbio.bench.amounts import dna_amount
from liulab_mbio.bench.materials import CUVETTE_ON_ICE, POLYMERASE_ON_ICE
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.io import read_record
from liulab_mbio.protocol.model import Citation, write_project, write_protocol
from liulab_mbio.sequence import SequenceRecord, span_text
from liulab_mbio.sites import digest, find_sites
from liulab_mbio.translate import translate
from liulab_synbio.igga import plan_igga
from liulab_synbio.igga.bench import CUTSMART, SPRI_BEADS, STRAIN
from liulab_synbio.igga.cargo import cargo_record
from liulab_synbio.igga.method import IGGA
from liulab_synbio.igga.protocols import ASSEMBLY, CREATION, FINAL
from liulab_synbio.igga.protocols.run import CUVETTES, FINAL_SELECTIVE, PREP_KIT, SELECTIVE
from liulab_synbio.igga.reads import ALLOWANCE, FLANK
from liulab_synbio.igga.vector import released_cargo

from ...chains import whole

DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library"

#: The enzymes the method reserves: the three a round uses, and BsmBI, which seats a part in its
#: carrier and so must not cut a designed block either.
RESERVED = ("BsaI", "BsmBI", "BbsI", "SrfI", "PmeI")

#: What each block's own stuffers put there, and nothing else: BsaI releases the part, BbsI
#: opens the destination it becomes, SrfI and PmeI shred what a round throws away.
EXPECTED_SITES = {"BsaI": 2, "BbsI": 2, "SrfI": 1, "PmeI": 2}


@pytest.fixture(scope="module")
def plan():
    """The demo as its own page runs it: the domesticated backbone, and the prices it ships."""
    return plan_igga(DEMO / "project.json", working_site="EGFP", prices=DEMO / "prices.csv")


@pytest.fixture(scope="module")
def protocol(plan):
    """Every protocol of the run as one, which is what most of these tests ask about."""
    return whole(plan.chain())


def test_the_demo_plans_every_part_and_the_whole_library(plan):
    assert [len(one) for one in plan.part_lists] == [24, 24, 24]
    assert len(plan.parts) == 72
    assert plan.constructs == 24**3 == 13824
    assert [one.position for one in plan.rounds] == ["N", "DBD", "C"]
    assert plan.status == "pass"


def test_no_reserved_enzyme_reads_a_part_outside_its_stuffers(plan):
    enzymes = tuple(get_enzyme(one) for one in RESERVED)
    for part in plan.parts:
        assert not find_sites(SequenceRecord(part.coding_sequence), enzymes), part.name
        found = find_sites(SequenceRecord(part.sequence), enzymes)
        assert Counter(site.enzyme.name for site in found) == EXPECTED_SITES, part.name


def test_the_product_reads_in_frame_from_its_first_part_to_its_last_barcode(plan):
    bases = str(plan.product.sequence)
    start = bases.index(plan.representative_parts[0].coding_sequence)
    stretch = bases[start : plan.rounds[-1].retained.end]
    assert len(stretch) % 3 == 0
    assert "*" not in translate(stretch)


def test_the_product_keeps_one_barcode_a_round_in_reverse_order(plan):
    bases = str(plan.product.sequence)
    scar = plan.scheme.cloning_scar
    codes = [one.barcode for one in plan.representative_parts]
    assert scar.join(reversed(codes)) in bases
    retained = bases[plan.rounds[-1].retained.start : plan.rounds[-1].retained.end]
    assert len(retained) == plan.project.retained_length == 75
    assert "*" not in translate(retained)


def test_the_block_meets_the_homopolymer_cap_at_a_scar_junction_and_never_passes_it(plan):
    # The block reads barcode-scar-barcode, so a run grows across a junction no barcode holds on
    # its own. Counted here plainly, over the whole context, so it fails if the designer's own
    # junction arithmetic ever admits a run this count can see.
    scar = plan.scheme.cloning_scar
    runs = [
        max(len(tuple(same)) for _, same in groupby(scar + one.barcode + scar))
        for one in plan.parts
    ]
    assert max(runs) <= MAX_HOMOPOLYMER


def test_the_pool_is_one_oligo_a_fragment_every_one_at_the_project_length(plan):
    pool = plan.pool.pool
    assert {len(one) for one in pool.oligos} == {plan.project.oligo_length}
    assert pool.spread == 0.0
    assert pool.count == len(pool.oligos)


def test_the_design_lands_on_the_arithmetic_floor_of_131(plan):
    """Splitting the cargo rather than the block leaves no block over its own floor.

    N_ATF7 overran while the stuffers were synthesised with it: they forced its cut positions,
    and the first forced one spelled GCCG, which is refused as one base kind.
    """
    assert plan.pool.floor == 131
    assert plan.pool.pool.count == 131
    assert plan.pool.over_floor == ()


def test_every_oligo_traces_to_its_protein_and_its_fragment(plan):
    names = {one.name for one in plan.parts}
    for oligo in plan.pool.pool.oligos:
        assert oligo.source in names
        assert 1 <= oligo.fragment <= oligo.fragments


def test_every_cargo_reassembles_from_its_own_fragments(plan):
    """What is synthesised is the cargo; the stuffers either side are the destination's."""
    for split, part in zip(plan.pool.splits, plan.parts, strict=True):
        assert split.reassembled() == str(cargo_record(part, plan.scheme).sequence)
        assert split.reassembled() in part.sequence


def test_every_position_has_a_block_vector_its_whole_part_list_closes_into(plan):
    """A part enters on its own position's overhang, and one backbone presents one pair."""
    assert [one.record.name for one in plan.block_vectors] == [
        "DMX-iGGA N",
        "DMX-iGGA DBD",
        "DMX-iGGA C",
    ]
    for index, one in enumerate(plan.block_vectors):
        backbone = next(
            piece
            for piece in digest(one.record, plan.scheme.internal)
            if (piece.start, piece.end) != (one.stuffer.start, one.stuffer.end)
        )
        cargo = [cargo_record(part, plan.scheme) for part in plan.parts if part.index == index]
        assert len(cargo) == 24
        assert {(str(each.sequence)[:4], str(each.sequence)[-4:]) for each in cargo} == {
            (backbone.right_overhang, backbone.left_overhang)
        }


def test_the_fragment_counts_sit_at_or_above_lund_s_measured_84_6_per_cent_point(plan):
    counted = dict(plan.pool.pool.fragment_counts())
    assert max(counted) == 5
    assert sum(blocks for pieces, blocks in counted.items() if pieces <= 2) == 56
    assert [row[2] for row in plan.pool.against_lund() if row[0] == 5] == [0.846]


def test_the_pool_reports_that_350_nt_has_no_slack_above_it(plan):
    assert "no slack above it at all" in "; ".join(str(one) for one in plan.pool.item.headroom)


def test_the_bill_carries_the_pool_row_with_its_band_and_the_demos_own_price(plan, protocol):
    """The largest line item is priced by the record the demo ships, and holed without one."""
    row = next(one for one in protocol.bill.rows if one.item.endswith("oligo pool"))
    assert (row.quantity, row.unit) == (131, "oligos")
    assert "131 count, 369 below the next band" in row.headroom
    assert "350 length, no slack above it at all" in row.headroom
    assert row.charge == "2575.00"
    assert row.citation == Citation("prices", "oligo-pool count 101-500; length 301-350")
    assert row.hole is None

    bare = next(
        one
        for one in replace(plan, prices=None).chain().bill.rows
        if one.item.endswith("oligo pool")
    )
    assert bare.charge == ""
    assert bare.hole is not None
    assert bare.hole.kind == "price"


def test_the_demos_price_record_prices_every_line_of_the_bill(plan, protocol):
    """Nine keys, nine rows: the demo ships the tariff a user would, so no money cell is a hole."""
    bill = protocol.bill
    assert [row.key for row in bill.rows] == [
        "oligo-pool",
        "pool-primers",
        "R3539",
        "R3733",
        "R0629",
        "R0560",
        "60242-2",
        "cuvettes",
        "plasmid prep",
    ]
    assert not [row.key for row in bill.rows if row.hole is not None]
    assert (bill.currency, bill.total) == ("USD", "3130.5875")
    record = plan.chain().sources["prices"]
    assert (record.document, record.edition) == ("AP-1 demo price record", "2026-10-08")
    assert record.read_as == "read from prices.csv"


def test_the_cargo_the_reads_run_across_is_the_cargo_the_release_digest_frees(plan):
    """Two derivations of one span, and they have to agree.

    The reads place the cargo by arithmetic over the rounds, because a destination is not
    obliged to free its own; the release step digests for it, because this destination does. A
    drift would hand the bench a pair that does not span what the digest put in the tube.
    """
    freed = released_cargo(plan.rounds[-1].product, plan.scheme)

    assert freed is not None
    over = plan.reads.linkage.amplicon_length - (freed.end - freed.start)
    assert 2 * (FLANK - ALLOWANCE) <= over <= 2 * (FLANK + ALLOWANCE)


def rerouted(plan, **changes):
    """Return the same plan with the project's validation keys replaced."""
    return replace(plan, project=replace(plan.project, **changes))


def read_back(plan, route="index PCR"):
    """Return the run's read-back on one route, of the two this demo offers."""
    return next(one for one in plan.validations if one.route.name == route)


def test_the_demo_reads_all_72_designs_back_as_288_wells(plan):
    """Both routes read the same designs into the same wells: they differ in how a well is marked."""
    assert [one.route.name for one in plan.validations] == ["barcode ligation", "index PCR"]
    for one in plan.validations:
        assert (one.floor, len(one.designs)) == (0, 72)
        assert one.wells == 288 == 72 * 4
        assert [len(picked.labels) for picked in one.picked] == [288]
    assert [index.name for index in read_back(plan).index] == ["index 1", "index 2", "index 3"]


def test_a_floor_above_a_design_shrinks_the_plates_by_exactly_what_it_leaves_out(plan):
    """The bench is sized from the designs read, not from the part list."""
    whole = read_back(plan)
    fewer = read_back(rerouted(plan, validate_from=2))
    left_out = len(whole.designs) - len(fewer.designs)
    assert left_out == 40
    assert fewer.wells == whole.wells - left_out * whole.colonies == 128
    assert sum(len(one.labels) for one in fewer.picked) == 128
    fewest = read_back(rerouted(plan, validate_from=3))
    assert (len(fewest.designs), fewest.wells) == (16, 64)
    assert len(fewest.index) == 1 < len(whole.index)
    assert rerouted(plan, validate_from=6).validations == ()


def test_a_design_is_read_in_the_pieces_the_pool_was_split_into(plan):
    """The count is the split's own, because arithmetic on the oligo length only bounds it."""
    counted = Counter(one.fragments for one in read_back(plan).designs)
    assert dict(plan.pool.pool.fragment_counts()) == counted
    assert max(counted) == 5


def test_a_project_with_no_floor_writes_a_protocol_with_no_validation(plan):
    """Cargo validation is optional, and a project that asks for none gets none."""
    polyclonal = rerouted(plan, validate_from=None, routes=())
    assert polyclonal.validations == ()
    titles = [step.title for step in whole(polyclonal.chain()).steps]
    assert "Pick 4 colonies of each design" not in titles
    assert "Order the oligo pool" in titles
    assert titles[titles.index("Order the oligo pool") + 5] == "Pool each part list"


def test_the_demo_carries_both_read_back_routes_as_the_two_ways_of_one_job(plan):
    """The demo offers both so a reader can read each; a run does one of them, never both."""
    chain = plan.chain()
    ways = [one for one in chain.protocols if one.choice]
    checks = {check.name: check for check in chain.audit()}

    assert [one.title for one in ways] == [
        "Cargo validation: barcode ligation",
        "Cargo validation: index PCR",
    ]
    assert {one.choice for one in ways} == {"read every well back"}
    left = {"clonal picked plate", "well calls"}
    assert [{item.name for item in one.produces} for one in ways] == [left, left]
    assert checks["choices"].status == "pass"
    assert (
        checks["choices"].detail
        == "2 ways to read every well back, each leaving the same 2 things."
    )
    assert checks["handoffs"].status == "pass"
    said = next(one for one in chain.background if one.title == "Read every well back")
    assert "Do one of them, never both." in said.body[0]
    assert [one.what for one in chain.inputs if "way only" in one.what] == [
        "the lab's own barcoding plasmids, one group a picked plate, for the barcode ligation "
        "way only",
        "the lab's own index primers, prepared once and called for by a run, for the index PCR "
        "way only",
    ]


def test_the_demo_emits_a_protocol_on_each_route(plan, protocol):
    """One set of parts, both routes: the ways of one job, never a branch the package picks."""
    index_pcr = protocol
    alone = rerouted(plan, routes=("barcode ligation",)).chain()
    # A build naming one route writes one page, offers no choice and is told nothing.
    assert [one.title for one in alone.protocols][3] == "Cargo validation: barcode ligation"
    assert not [one for one in alone.protocols if one.choice]
    assert [one.title for one in alone.background] == ["How this library is designed"]
    assert [check.name for check in alone.audit()] == ["handoffs", "sources"]
    ligation = whole(alone)
    assert "Amplify each well with its own pair" in [one.title for one in index_pcr.steps]
    assert "Barcode each well in lysate" in [one.title for one in ligation.steps]
    # Both cycle counts stay open: this build states neither, because nobody ran the pilot.
    pcrs = ["H29", "H30"]
    # Block assembly holds nothing open: every position has a destination presenting its own
    # entry overhang, and NEB's kit table sizes the reaction.
    blocks: list[str] = []
    # The linkage read and the final assembly are both answered by what this build states: its
    # own pass mark, its own working vector, and its own mass and ratio for the one-pot tube.
    linkage: list[str] = []
    final: list[str] = []
    assert [hole.id for step in index_pcr.steps for hole in step.holes] == [
        *pcrs,
        *blocks,
        "IDX1",
        *linkage,
        *final,
    ]
    assert [hole.id for step in ligation.steps for hole in step.holes] == [
        *pcrs,
        *blocks,
        *linkage,
        *final,
    ]
    for one in (ligation, index_pcr):
        assert [check.status for check in one.audit()] == ["pass", "pass", "pass", None]


def test_the_demo_is_left_with_the_four_numbers_no_input_of_its_own_can_give(plan):
    """The demo's floor, counted as the run page counts it: the bill's holes and each
    protocol's, one entry an id. Anything else here is a build field that stopped being read.
    """
    chain = plan.chain()
    found = {row.hole.id: row.hole for row in (chain.bill.rows if chain.bill else ()) if row.hole}
    for protocol in chain.protocols:
        for hole in protocol.all_holes:
            found.setdefault(hole.id, hole)

    assert {one: hole.kind for one, hole in found.items()} == {
        "H29": "unpublished",
        "H30": "unpublished",
        "IDX1": "lab",
        "H23": "unpublished",
    }


def test_a_chance_too_small_to_print_fixed_prints_as_a_power_of_ten(protocol):
    """A step builder formats it as the page does, so no line reads ``7.27e-07``."""
    grow = next(one for one in protocol.steps if one.title.startswith("Round 3: recover"))
    assert "At that count the chance a named product is missing is 7.27 × 10⁻⁷." in grow.expected
    assert not [line for one in protocol.steps for line in one.expected if re.search(r"\de-", line)]


def test_a_pools_picomoles_print_three_figures_and_not_six(protocol):
    """Every pmol a round states goes through the page's own formatter, not ``{:g}``."""
    said = {one.title: one.expected for one in protocol.steps}
    assert said["Round 1: clean both digests up"] == (
        "DMX-iGGA, opened: 0.00605 pmol is 20 ng at 5363 bp.",
        "N part list, released: 0.00605 pmol is 2.08 ng at 559 bp.",
    )
    assert said["Round 1: ligate the N part list into the library"] == (
        "Nothing visible. 0.00605 pmol of the released part list against 0.00605 pmol of the "
        "opened library is the 1:1 molar ratio.",
    )
    assert said["Round 1: clean the ligation up"] == (
        "Enough for one electroporation: 100 ng is 0.0264 pmol at 6161 bp.",
    )


def test_the_protocol_builds_the_blocks_it_has_a_pool_for_rather_than_ordering_them(protocol):
    """With a pool designed, nothing is ordered as a block: the pool is, resuspended, then used."""
    titles = [step.title for step in protocol.steps]
    start = titles.index("Order the oligo pool")
    assert titles[start : start + 5] == [
        "Order the oligo pool",
        "Resuspend the oligo pool",
        "PCR1: pull 1 batch out of the pool",
        "PCR2: pull each of the 72 blocks out of its batch",
        "Assemble each cargo into its position's destination, from its 1 to 5 pieces",
    ]
    note = next(one.note for one in protocol.materials if one.name == "N part list")
    assert note == "assembled from the oligo pool; pool.tsv says which oligos"


def test_the_pool_is_in_buffer_before_anything_amplifies_it(protocol):
    """PCR1 takes 20 ng/µL of template, so the step before it says how the pool got there."""
    made = next(one for one in protocol.steps if one.title == "Resuspend the oligo pool")
    assert made.instructions[0] == (
        "Divide the total yield in ng printed on the shipping tube label by 20, rounding down, "
        "to get the resuspension volume in µL."
    )
    assert "10 mM Tris buffer, pH 8.0" in made.instructions[1]
    assert made.expected == (
        "One tube of pool in solution at 20 ng/µL or above, with nothing left undissolved on the "
        "wall of the tube.",
    )
    assert not made.holes


def test_the_same_dna_is_billed_once(plan, protocol):
    """A pool buys oligos and primers; a project without one buys blocks. Never both."""
    pooled = [row.item for row in protocol.bill.rows]
    assert "Synthesised blocks" not in pooled
    assert pooled[:2] == ["AP-1 DESynR oligo pool", "Pool amplification primers"]
    primers = next(row for row in protocol.bill.rows if row.item.endswith("primers"))
    assert (primers.quantity, primers.unit) == (74, "primers")
    unpooled = [row.item for row in replace(plan, pool=None).chain().bill.rows]
    assert unpooled[0] == "Synthesised blocks"
    assert "Pool amplification primers" not in unpooled


def test_pcr1_cites_its_cycles_from_the_pool_length_and_pcr2_leaves_them_blank(protocol):
    """A 350 nt pool sits in Twist's top band; nothing sources the count for PCR2's template."""
    steps = {one.title.split(":")[0]: one for one in protocol.steps}
    first, second = steps["PCR1"], steps["PCR2"]

    cycling = first.programs[0].stages[1]
    assert cycling.cycles == 12
    assert cycling.citation is not None
    assert cycling.citation.source in protocol.sources
    assert second.programs[0].stages[1].cycles is None
    assert [hole.id for hole in second.holes] == ["H30"]


def test_the_pulse_prints_on_the_row_that_names_the_cells_manual(protocol):
    """The program belongs to the cells, so the settings and their source sit on one row."""
    cells = next(one for one in protocol.materials if one.catalog == "60242-2")
    assert cells.citation == Citation("MA133", "p. 4-5")
    assert "1800 V, 600 Ω and 10 µF" in cells.note
    step = next(one for one in protocol.steps if "electroporate" in one.title)
    assert any("1800 V, 600 Ω and 10 µF" in line for line in step.instructions)
    assert {"MA133", "Qian SI"} <= set(protocol.sources)


def test_the_assembly_step_names_a_destination_a_position_and_sizes_itself_from_nebs_table(
    plan, protocol
):
    """Each cargo closes into its own position's vector, in NEB's own kit reaction."""
    step = next(one for one in protocol.steps if one.title.startswith("Assemble"))
    opened = len(plan.rounds[0].destination) - plan.rounds[0].excised.length

    table = step.tables[0]
    assert "E1602" in table.title
    rows = {one.name: one for one in table.components}
    weighed = dna_amount("backbone", opened, pmol=0.05).nanograms
    assert rows["DMX-iGGA, opened"].final == f"0.05 pmol ({weighed:g} ng)"
    assert rows["T4 DNA Ligase Buffer"].stock == "10X"
    assert [one.name for one in table.components if one.name.startswith("PCR2 piece")] == [
        f"PCR2 piece {number}" for number in range(1, 6)
    ]
    assert step.programs[0].stages[-1].incubations[0].temperature_c == 60.0
    assert not step.holes
    assert step.instructions[0].startswith(
        "Open DMX-iGGA N (block-vector-1.dna), DMX-iGGA DBD (block-vector-2.dna) and "
        "DMX-iGGA C (block-vector-3.dna) with BbsI"
    )


def test_the_final_assembly_is_written_out_from_what_the_build_states(protocol):
    """Five steps, every one of them there, and no hole left where the build names its own."""
    steps = protocol.steps[-5:]

    assert [one.title for one in steps] == [
        "Pick the working vector and confirm PaqCI opens it",
        "Release the cargo with BsaI and PmeI",
        "Assemble the cargo into pLVX-TetOne-dom with PaqCI",
        "Clean the assembly up and electroporate into Endura ElectroCompetent Cells",
        "Read representation in the final vector",
    ]
    for step in steps:
        assert step.instructions, step.title
        assert step.expected, step.title
    assert "At least 195,386 net colonies: the floor for the 0.99 chance" in " ".join(
        steps[3].expected
    )
    # The paper runs its own last transfer by the methods it gives for a round, so the release
    # and the clean-up around the assembly are sized rather than held open.
    assert steps[1].tables[0].components[0].final.endswith("(1000 ng)")
    assert "at 1x the volume" in steps[3].instructions[0]
    assert steps[0].instructions[1].startswith("Have it made to working-vector-ccdb.dna")
    assert steps[2].instructions[0].startswith("Add 75 ng of working vector")
    # The cargo is never pipetted, so its row is the whole release and its weight follows from
    # the vector's and the ratio the build states.
    cargo, vector = steps[2].tables[0].components[:2]
    assert (cargo.final, cargo.volume_ul) == ("0.0246 pmol (16.31 ng)", 50.0)
    assert vector.final == "0.0123 pmol (75 ng)"
    assert steps[2].notes[-1] == (
        "75 ng of vector, at 2:1 cargo to vector, is what this run measured, not a published "
        "figure."
    )
    assert not [hole.id for step in steps for hole in step.holes]


def test_the_demo_puts_its_cassette_at_the_first_of_two_egfp_annotations(plan):
    """Addgene annotates EGFP twice on this deposit, so which one answers `--working-site` is set.

    The two spans share an end and differ by three bases at the start, and the cassette goes at
    the start of the first the record lists. Both annotations survive, shifted by what went in.
    """
    backbone = read_record(DEMO / "working-vector.gb")
    before = [one.segments[0].start for one in backbone.features if one.name == "EGFP"]
    assert before == [2513, 2516]

    after = [one.segments[0].start for one in plan.working.record.features if one.name == "EGFP"]
    cassette = len(plan.working.record) - len(backbone)
    assert after == [start + cassette for start in before]
    assert plan.working.enzyme.name == "PaqCI"


def test_a_build_naming_no_working_vector_carries_the_hole_one_fills(plan):
    """H31 is the demo's own input doing the work: take the vector away and the hole comes back."""
    steps = whole(replace(plan, working=None).chain()).steps[-5:]

    assert steps[0].title == "Pick the working vector"
    assert steps[2].title == "Assemble the cargo into the working vector"
    assert not steps[2].programs
    assert [hole.id for step in steps for hole in step.holes] == ["H31", "H31"]


def test_the_read_backs_plates_are_declared_and_every_well_resolves(protocol):
    """Index PCR pours picked and index plates, and each now belongs to a plate the page draws."""
    one = protocol
    named = {
        well.plate
        for step in one.steps
        for transfer in step.transfers
        for move in transfer.moves
        for well in (move.source, move.destination)
    }
    drawn = {plate.name: plate for plate in one.plates}

    assert named <= set(drawn)
    assert {"picked 1", "index 1"} <= set(drawn)
    assert drawn["picked 1"].labels["A1"] == "quarter 1"
    assert not drawn["picked 1"].seating
    assert [c.status for c in one.audit() if c.name == "wells"] == ["pass"]


def test_the_demo_names_kanamycin_wherever_the_paper_named_carbenicillin(protocol, tmp_path):
    """The destination is KanR, so no plate, well or broth of a round carries the drug D11 rules out.

    Asserted on the file the pipeline ships, so a vessel's or a plate's wording counts too. The
    working vector is AmpR and its own transfer does plate on carbenicillin, which is the one
    line the ban does not cover.
    """
    written = write_protocol(protocol, tmp_path / "protocol.json").read_text()
    rounds = [
        line
        for line in written.splitlines()
        if "carbenicillin" in line or "ampicillin" in line
        if "working vector's own marker" not in line
    ]

    assert "kanamycin" in written
    assert rounds == []


def test_the_run_is_one_protocol_a_sitting_and_every_handover_resolves(plan):
    """What one page hands the next is the point of the split, so the chain is checked."""
    chain = plan.chain()

    assert [one.title for one in chain.protocols] == [
        "Primer plates",
        "Cargo ordering and pool preparation",
        "Cargo creation",
        "Cargo validation: barcode ligation",
        "Cargo validation: index PCR",
        "Library assembly in rounds",
        "Final cargo ligation",
    ]
    assert [len(one.steps) for one in chain.protocols] == [5, 2, 3, 6, 6, 27, 5]
    assert [one.audit()[0].status for one in (chain,)] == ["pass"]
    handed = {item.name for item in chain.inputs}
    for one in chain.protocols:
        assert {item.name for item in one.consumes} <= handed, one.title
        handed |= {item.name for item in one.produces}
    assert "the library in its working vector" in handed


def test_a_repeated_caution_rides_its_material_and_no_step_of_the_run_stores_one(plan):
    """A caution the bench reads on four pages is one sentence, carried by the tube it is about."""
    chain = plan.chain()
    shown = {
        (one.title, step.key): one.cautions_for(step)
        for one in chain.protocols
        for step in one.steps
        if one.cautions_for(step)
    }

    assert shown == {
        ("Primer plates", "resuspend-primers"): (
            "Spin the plate down before taking the seal off.",
        ),
        ("Cargo ordering and pool preparation", "resuspend-pool"): (
            "Spin the tube down before taking the cap off.",
        ),
        ("Cargo creation", "pcr1"): (POLYMERASE_ON_ICE,),
        ("Cargo creation", "pcr2"): (POLYMERASE_ON_ICE,),
        ("Cargo validation: index PCR", "index-pcr"): (POLYMERASE_ON_ICE,),
        ("Library assembly in rounds", "round-1-electroporate"): (CUVETTE_ON_ICE,),
        ("Library assembly in rounds", "round-2-electroporate"): (CUVETTE_ON_ICE,),
        ("Library assembly in rounds", "round-3-electroporate"): (CUVETTE_ON_ICE,),
        ("Final cargo ligation", "electroporate-and-grow"): (CUVETTE_ON_ICE,),
    }
    written = {step.key for one in chain.protocols for step in one.steps if step.cautions}
    assert written == {"resuspend-primers", "resuspend-pool"}


def test_every_step_sits_under_a_stage_of_its_own_protocol(plan):
    """A page of 27 steps reads as rounds, so each step says which stage it belongs to."""
    chain = plan.chain()
    rounds = next(one for one in chain.protocols if one.title == "Library assembly in rounds")

    written = [one for one in chain.protocols if one.title != "Primer plates"]
    assert all(step.section for one in written for step in one.steps)
    assert list(dict.fromkeys(step.section for step in rounds.steps)) == [
        "Pool the part lists",
        "Round 1",
        "Round 2",
        "Round 3",
        "Read the library back",
    ]


def test_a_route_the_package_does_not_ship_is_refused(plan):
    with pytest.raises(ValueError, match="route is 'both'"):
        rerouted(plan, routes=("both",))


def test_the_primer_plates_are_written_only_where_the_project_says_how(plan):
    """The amounts are nobody's to guess, so a project stating none gets no such sitting."""
    plates = plan.chain().protocols[0]
    bare = rerouted(plan, primer_plates=None).chain()

    assert plates.title == "Primer plates"
    assert len(plates.oligos) == len(plan.pool.pool.primers) == 74
    assert [one.name for one in plates.produces] == [
        "primer stock plate",
        "primer working plate 1",
    ]
    assert [one.name for one in plates.plates] == [
        "primer stock plate",
        "primer working plate 1",
    ]
    assert bare.protocols[0].title == "Cargo ordering and pool preparation"
    assert bare.audit()[0].status == "pass"


def test_the_primers_are_ordered_once_however_the_run_is_split(plan):
    """A plated run orders its primers in that sitting, so the pool step stops ordering them."""
    plated = [step.title for step in whole(plan.chain()).steps]
    bare = [step.title for step in whole(rerouted(plan, primer_plates=None).chain()).steps]

    assert "Order the oligo pool" in plated
    assert "Order the primers" in plated
    assert "Order the oligo pool and the primers that amplify it" in bare
    assert "Order the primers" not in bare


def test_no_rendered_field_prints_a_parenthesised_plural(plan, tmp_path):
    """A reader follows a count and its noun, never ``tube(s)``: `counted` says both."""
    written = tmp_path / "project.json"
    write_project(plan.chain(), written)

    assert "(s)" not in written.read_text(encoding="utf-8")


def test_the_linkage_read_prints_the_block_1_based_and_inclusive(plan, protocol):
    """A span the bench reads counts from 1 and includes its last base, as the map prints one."""
    last = plan.rounds[-1]
    step = next(one for one in protocol.steps if one.title == "Read linkage")

    said = " ".join(step.expected)

    assert span_text(last.block.start, last.block.end, len(last.product)) in said
    assert f"{last.block.start}-{last.block.end}" not in said


def test_a_page_lists_only_the_reagents_its_own_steps_reach(plan):
    """A materials table sends nobody to a freezer for a reagent no step of that page uses."""
    listed = {page.title: {one.name for one in page.materials} for page in plan.chain().protocols}
    working = plan.working
    assert working is not None
    destination = f"{plan.vector.name or 'destination'} vector"
    carried = f"{working.record.name or 'Working'} vector"
    creation, rounds, final = listed[CREATION], listed[ASSEMBLY], listed[FINAL]

    # Protocol 03 opens a destination, cleans up and transforms; it preps nothing and never
    # reaches the vector the library ends in.
    assert {destination, CUTSMART, SPRI_BEADS, STRAIN} <= creation
    assert any(one.startswith(IGGA.internal.name) for one in creation)
    assert not {carried, PREP_KIT, CUVETTES, FINAL_SELECTIVE} & creation
    assert destination in rounds
    assert carried not in rounds
    assert carried in final
    assert destination not in final


def test_every_reagent_the_rounds_share_lands_on_a_page(plan):
    """Deriving a page's list from its steps may narrow a table; it may never lose an order."""
    working = plan.working
    assert working is not None
    bought = {one.name for page in plan.chain().protocols for one in page.materials}

    assert {
        f"{plan.vector.name or 'destination'} vector",
        f"{working.record.name or 'Working'} vector",
        CUTSMART,
        CUVETTES,
        FINAL_SELECTIVE,
        PREP_KIT,
        SELECTIVE,
        SPRI_BEADS,
        STRAIN,
    } <= bought
