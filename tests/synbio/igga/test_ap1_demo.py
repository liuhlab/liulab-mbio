"""The AP-1 demo planned end to end, from the amino acid sequences and nothing else.

The inputs are `docs/examples/ap1-library`, specified by `docs/research/ap1-demo-project.md`:
72 proteins, a project file and a destination. Every assertion here is on the planned product and
the blocks that make it, not on how the design reached them.
"""

from collections import Counter
from dataclasses import replace
from pathlib import Path

import pytest

from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.protocol.model import Citation
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import find_sites
from liulab_mbio.translate import translate
from liulab_synbio.igga import plan_igga

DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library"

#: The enzymes the method reserves: the three a round uses, and BsmBI, which seats a part in its
#: carrier and so must not cut a designed block either.
RESERVED = ("BsaI", "BsmBI", "BbsI", "SrfI", "PmeI")

#: What each block's own stuffers put there, and nothing else: BsaI releases the part, BbsI
#: opens the destination it becomes, SrfI and PmeI shred what a round throws away.
EXPECTED_SITES = {"BsaI": 2, "BbsI": 2, "SrfI": 1, "PmeI": 2}


@pytest.fixture(scope="module")
def plan():
    return plan_igga(DEMO / "project.json")


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


def test_the_pool_is_one_oligo_a_fragment_every_one_at_the_project_length(plan):
    pool = plan.pool.pool
    assert {len(one) for one in pool.oligos} == {plan.project.oligo_length}
    assert pool.spread == 0.0
    assert pool.count == len(pool.oligos)


def test_the_arithmetic_floor_is_152_and_the_design_spends_one_more_on_one_block(plan):
    """152 is what the length budget allows for; 153 is what a legal overhang set costs.

    N_ATF7 is 1,092 bp, exactly four spans of 276 less the three shared overhangs, so its four
    cut positions are forced and the first spells GCCG, which is refused as one base kind.
    """
    assert plan.pool.floor == 152
    assert plan.pool.pool.count == 153
    assert plan.pool.over_floor == ("N_ATF7",)


def test_every_oligo_traces_to_its_protein_and_its_fragment(plan):
    names = {one.name for one in plan.parts}
    for oligo in plan.pool.pool.oligos:
        assert oligo.source in names
        assert 1 <= oligo.fragment <= oligo.fragments


def test_every_block_reassembles_from_its_own_fragments(plan):
    for split, part in zip(plan.pool.splits, plan.parts, strict=True):
        assert split.reassembled() == part.sequence


def test_the_fragment_counts_sit_at_or_above_lund_s_measured_84_6_per_cent_point(plan):
    counted = dict(plan.pool.pool.fragment_counts())
    assert max(counted) == 5
    assert sum(blocks for pieces, blocks in counted.items() if pieces <= 2) == 55
    assert [row[2] for row in plan.pool.against_lund() if row[0] == 5] == [0.846]


def test_the_pool_reports_that_350_nt_has_no_slack_above_it(plan):
    assert "no slack above it at all" in "; ".join(str(one) for one in plan.pool.item.headroom)


def test_the_bill_carries_the_pool_row_with_its_band_and_a_money_hole(plan):
    """The largest line item is on the bill; with no tariff loaded its money cell is a hole."""
    row = next(one for one in plan.protocol().bill.rows if one.item.endswith("oligo pool"))
    assert (row.quantity, row.unit) == (153, "oligos")
    assert "153 count, 347 below the next band" in row.headroom
    assert "350 length, no slack above it at all" in row.headroom
    assert row.charge == ""
    assert row.hole is not None
    assert row.hole.kind == "price"


def rerouted(plan, **changes):
    """Return the same plan with the project's validation keys replaced."""
    return replace(plan, project=replace(plan.project, **changes))


def test_the_demo_reads_all_72_designs_back_as_288_wells(plan):
    one = plan.validation
    assert (one.route.name, one.floor, len(one.designs)) == ("B", 0, 72)
    assert one.wells == 288 == 72 * 4
    assert [len(picked.seating) for picked in one.picked] == [288]
    assert [index.name for index in one.index] == ["index 1", "index 2", "index 3"]


def test_a_floor_above_a_design_shrinks_the_plates_by_exactly_what_it_leaves_out(plan):
    """The bench is sized from the designs read, not from the part list."""
    whole = plan.validation
    fewer = rerouted(plan, validate_from=2).validation
    left_out = len(whole.designs) - len(fewer.designs)
    assert left_out == 21
    assert fewer.wells == whole.wells - left_out * whole.colonies == 204
    assert sum(len(one.seating) for one in fewer.picked) == 204
    fewest = rerouted(plan, validate_from=3).validation
    assert (len(fewest.designs), fewest.wells) == (17, 68)
    assert len(fewest.index) == 1 < len(whole.index)
    assert rerouted(plan, validate_from=6).validation is None


def test_a_design_is_read_in_the_pieces_the_pool_was_split_into(plan):
    """The count is the split's own, because arithmetic on the oligo length only bounds it."""
    counted = Counter(one.fragments for one in plan.validation.designs)
    assert dict(plan.pool.pool.fragment_counts()) == counted
    assert max(counted) == 5


def test_a_project_with_no_floor_writes_a_protocol_with_no_validation(plan):
    """Cargo validation is optional, and a project that asks for none gets none."""
    polyclonal = rerouted(plan, validate_from=None, route=None)
    assert polyclonal.validation is None
    titles = [step.title for step in polyclonal.protocol().steps]
    assert "Pick 4 colonies of each design" not in titles
    assert titles[0] == "Order the oligo pool and the primers that amplify it"
    assert titles[4] == "Pool each part list"


def test_the_demo_emits_a_protocol_on_each_route(plan):
    """One set of parts, two project files: a second project is never a second branch."""
    route_b = plan.protocol()
    route_a = rerouted(plan, route="A").protocol()
    assert "Amplify each well with its own pair" in [one.title for one in route_b.steps]
    assert "Barcode each well in lysate" in [one.title for one in route_a.steps]
    pcrs = ["H29", "H30"]
    reads = ["H27", "H28"]
    assert [hole.id for step in route_b.steps for hole in step.holes] == [
        *pcrs,
        "H25",
        "H26",
        "B1",
        "B2",
        *reads,
    ]
    assert [hole.id for step in route_a.steps for hole in step.holes] == [
        *pcrs,
        "H25",
        "H26",
        *reads,
    ]
    for one in (route_a, route_b):
        assert [check.status for check in one.audit()] == ["pass", "pass", "pass", None]


def test_the_protocol_builds_the_blocks_it_has_a_pool_for_rather_than_ordering_them(plan):
    """With a pool designed, nothing is ordered as a block: the pool is, and four steps follow."""
    titles = [step.title for step in plan.protocol().steps]
    assert titles[:4] == [
        "Order the oligo pool and the primers that amplify it",
        "PCR1: pull 1 batch out of the pool",
        "PCR2: pull each of the 72 blocks out of its batch",
        "Assemble each block from its 1 to 5 pieces",
    ]
    note = next(one.note for one in plan.protocol().materials if one.name == "N part list")
    assert note == "assembled from the oligo pool; pool.tsv says which oligos"


def test_the_same_dna_is_billed_once(plan):
    """A pool buys oligos and primers; a project without one buys blocks. Never both."""
    pooled = [row.item for row in plan.protocol().bill.rows]
    assert "Synthesised blocks" not in pooled
    assert pooled[:2] == ["AP-1 DESynR oligo pool", "Pool amplification primers"]
    primers = next(row for row in plan.protocol().bill.rows if row.item.endswith("primers"))
    assert (primers.quantity, primers.unit) == (74, "primers")
    unpooled = [row.item for row in replace(plan, pool=None).protocol().bill.rows]
    assert unpooled[0] == "Synthesised blocks"
    assert "Pool amplification primers" not in unpooled


def test_pcr1_cites_its_cycles_from_the_pool_length_and_pcr2_leaves_them_blank(plan):
    """A 350 nt pool sits in Twist's top band; nothing sources the count for PCR2's template."""
    protocol = plan.protocol()
    steps = {one.title.split(":")[0]: one for one in protocol.steps}
    first, second = steps["PCR1"], steps["PCR2"]

    cycling = first.programs[0].stages[1]
    assert cycling.cycles == 12
    assert cycling.citation is not None
    assert cycling.citation.source in protocol.sources
    assert second.programs[0].stages[1].cycles is None
    assert [hole.id for hole in second.holes] == ["H30"]


def test_the_pulse_prints_on_the_row_that_names_the_cells_manual(plan):
    """The program belongs to the cells, so the settings and their source sit on one row."""
    protocol = plan.protocol()
    cells = next(one for one in protocol.materials if one.catalog == "60242-2")
    assert cells.citation == Citation("MA133", "p. 4-5")
    assert "1800 V, 600 Ω and 10 µF" in cells.note
    step = next(one for one in protocol.steps if "electroporate" in one.title)
    assert any("1800 V, 600 Ω and 10 µF" in line for line in step.instructions)
    assert {"MA133", "Qian SI"} <= set(protocol.sources)


def test_the_assembly_step_names_what_nobody_decided_rather_than_a_number(plan):
    """A block is split for synthesis and nothing says what puts it back together."""
    step = next(one for one in plan.protocol().steps if one.title.startswith("Assemble"))
    assert [hole.id for hole in step.holes] == ["H25", "H26"]
    assert not step.tables
    assert not step.programs
    assert all(hole.kind == "undecided" for hole in step.holes)
