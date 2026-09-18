"""The whole pipeline, run on the fixtures with no agent.

`tests/data/pUC19.dna` and `tests/data/GFP.dna`: put GFP into the pUC19 multiple cloning site by
Gibson assembly. Nothing about the fixtures is hard-coded in the package; the numbers below are
read off the records and pinned here.
"""

import pytest

from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.bench.steps import COLONY_PCR_TITLE, SEQUENCING_TITLE, quantify_step
from liulab_mbio.bench.validation import SANGER_ALLOWANCE, SANGER_FLANK, ColonyCheck, SangerRead
from liulab_mbio.cloning.gibson import Plan, plan_gibson
from liulab_mbio.cloning.gibson.bench import IN_FUSION, NEBUILDER_HIFI
from liulab_mbio.cloning.gibson.design import (
    BRIDGE_HOMOLOGY_BP,
    STITCH_OLIGO_BASES,
    STITCH_OVERLAP_BP,
    wallace_tm,
)
from liulab_mbio.cloning.gibson.steps import (
    CLEANUP_FRAGMENTS,
    CORRECT_AT_FIVE,
    MOLECULES_PER_ERROR,
    SCREENED_COLONIES,
)
from liulab_mbio.protocol import OVERVIEW_CHARS, read_protocol, render_html
from liulab_mbio.sequence import SequenceRecord, Strand, reverse_complement
from liulab_mbio.snapgene import read_dna

#: Where the fixture's own MCS feature sits.
MCS = (395, 452)

#: Where GFP is cut in two, to be amplified and assembled as two inserts. Off any repeat and
#: far enough from either end for a primer.
SPLIT = 360


@pytest.fixture(scope="module")
def made(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP into the pUC19 multiple cloning site, every option left at its default."""
    return plan_gibson(puc19, gfp)


def test_each_junction_takes_its_overlap_from_the_vector_at_the_length_the_note_states(made, puc19):
    rule = NEBUILDER_HIFI.tiers[0].overlap
    floor = rule.tm_floor
    assert floor is not None
    assert floor == 48.0
    first, last = made.assembly.junctions
    assert [len(one) for one in made.overlaps] == [16, 16]
    for one in (first, last):
        assert rule.shortest <= one.length <= rule.longest
        assert wallace_tm(one.overlap) >= floor
        assert one.taken_from == made.linearised_vector.name
    # The bases are the vector's own, on either side of the span the insert replaces.
    assert first.overlap == puc19.sequence[MCS[0] - first.length : MCS[0]]
    assert last.overlap == puc19.sequence[MCS[1] : MCS[1] + last.length]


def test_the_insert_primers_carry_the_overlaps_as_tails_and_the_vector_primers_carry_none(made):
    first, last = made.assembly.junctions
    insert = made.insert_parts[0]
    assert (insert.left_tail, insert.right_tail) == (first.overlap, last.overlap)
    assert insert.forward.sequence.startswith(first.overlap)
    assert (made.linearised_vector.left_tail, made.linearised_vector.right_tail) == ("", "")
    # The amplicon is longer than what the part puts into the product, by the two tails.
    assert insert.length == insert.fragment_length + first.length + last.length


def test_the_tm_of_a_tailed_primer_is_read_from_its_annealing_region_alone(made):
    report = next(one for one in made.reports if one.primer.name == "GFP forward")
    annealing = report.primer.binding_sites[0]
    assert len(report.primer.sequence) == 42
    assert annealing.end - annealing.start == 26
    assert report["length"].value == 26
    # The whole oligo's melting temperature is reported beside it, and carries no verdict.
    assert report["tm_full"].value > report["tm"].value
    assert report["tm_full"].status is None
    # The structures are judged on the whole oligo, tail included.
    assert report["hairpin"].status is not None


def test_the_product_carries_every_feature_at_its_new_coordinates(made, puc19, gfp):
    product = made.plasmid
    assert len(product) == len(puc19) - (MCS[1] - MCS[0]) + len(gfp) == 3346
    features = {
        one.name: [(bit.start, bit.end) for bit in one.segments] for one in product.features
    }
    # The vector's reporter meets the replaced span twice, so it keeps a segment either side.
    reporter = next(name for name in features if name.startswith("lacZ"))
    assert features[reporter] == [(145, 395), (1112, 1129)]
    assert features["GFP"] == [(395, 1112)]
    # A feature lying wholly inside the span the insert replaced is gone.
    assert "MCS" not in features
    # One further round the plasmid, shifted by what the insert added and the span removed.
    assert features["M13 rev"] == [(1124, 1141)]


def test_the_product_keeps_the_vectors_own_origin_so_no_junction_sits_at_base_zero(made, puc19):
    assert made.plasmid.sequence[: MCS[0]] == puc19.sequence[: MCS[0]]
    assert made.plasmid.topology == "circular"
    assert made.assembly.junction_positions == (379, 1112)
    assert all(one.start > 0 for one in made.assembly.junctions)


def test_each_junction_is_marked_and_says_which_part_the_bases_came_from(made):
    marks = {one.name: one for one in made.plasmid.features if one.name.endswith("overlap")}
    assert set(marks) == {"pUC19 backbone-GFP overlap", "GFP-pUC19 backbone overlap"}
    for junction in made.assembly.junctions:
        mark = marks[f"{junction.before}-{junction.after} overlap"]
        assert (mark.segments[0].start, mark.segments[0].end) == (junction.start, junction.end)
        assert junction.taken_from in mark.qualifiers["note"][0]


def test_every_designed_primer_is_annotated_where_it_binds_on_the_product(made):
    placed = {one.name: one.binding_sites[0] for one in made.plasmid.primers}
    designed = [
        one.report.primer.name for one in made.designed_oligos if one.role == "amplification"
    ]
    assert designed == [
        "pUC19 backbone forward",
        "pUC19 backbone reverse",
        "GFP forward",
        "GFP reverse",
    ]
    assert set(designed) <= set(placed)
    for name in designed:
        assert 0 <= placed[name].start < placed[name].end <= len(made.plasmid)
    # The insert's forward primer anneals to the insert's own first bases, past its tail.
    assert (placed["GFP forward"].start, placed["GFP forward"].end) == (395, 421)


def test_a_part_gives_way_where_the_shared_bases_end_and_not_where_they_begin(made, gfp):
    # Both overlaps are the vector's own bases: the first sits at the end of the backbone's
    # share and the last at the start of it. Reading the junction starts as boundaries instead
    # would put a screening primer 16 bases out and turn vector bases round with the insert.
    first, last = made.assembly.junctions
    assert made.assembly.junction_positions == (first.start, last.start)
    assert made.assembly.boundaries == (first.end, last.start)
    assert made.assembly.boundaries == made.assembly.insert_span == (MCS[0], MCS[0] + len(gfp))


def test_the_colony_pcr_sizes_are_the_ones_the_simulated_clones_give(made, gfp):
    bands = {clone.name: clone.bands_bp for clone in made.colony.clones}
    start, end = made.assembly.boundaries
    removed = made.span[1] - made.span[0]
    forward, reverse, junction = (one.binding_sites[0] for one in made.colony.primers)
    # Where the three primers landed inside their placements, which is what the bands count off.
    ahead, behind, into = start - forward.start, reverse.end - end, junction.end - start
    assert forward.end <= start < junction.end <= end <= reverse.start
    # The junction primer reaches the near vector primer and the two vector primers span the
    # insert. A turned insert shares no overlap with the vector, so there is no lane for one.
    assert bands == {
        "Correct clone": (ahead + into, ahead + len(gfp) + behind),
        "Empty vector": (ahead + removed + behind,),
    }
    assert {lane.label: lane.bands_bp for lane in made.colony.gel.lanes} == bands


def test_the_sequencing_primers_read_in_from_outside_the_first_and_the_last_junction(made):
    first, last = made.assembly.junctions
    start, end = made.assembly.boundaries
    forward, reverse = made.reads
    ahead, behind = (read.primer.binding_sites[0] for read in made.reads)
    assert ahead.end <= first.start
    assert behind.start >= last.end
    for read in made.reads:
        assert SANGER_FLANK <= read.distance_bp <= SANGER_FLANK + SANGER_ALLOWANCE
    # Each has to carry as far as the other's junction for one read to confirm both.
    assert (forward.read_bp, reverse.read_bp) == (end - ahead.end, behind.start - start)


def test_the_page_says_what_the_plate_should_look_like_from_the_products_own_features(made):
    protocol = made.protocol()
    phenotype = made.phenotype
    assert phenotype.reporter is not None
    assert protocol.overview["Selection"] == phenotype.antibiotic
    assert phenotype.reporter.name in " ".join(protocol.highlights)
    plating = next(step for step in protocol.steps if step.title == "Transform and plate")
    told = " ".join((*plating.expected, *plating.notes))
    assert "white" in told
    assert "blue" in told
    # The lac promoter reads the other way and no ribosome binding site is annotated.
    assert not phenotype.expressed
    assert "not expected to make GFP" in told


def test_the_screening_steps_print_the_notes_numbers_and_cite_where_each_came_from(made):
    protocol = made.protocol()
    steps = {step.title: step for step in protocol.steps}
    purify = " ".join(steps["Purify every amplicon"].notes)
    assert f"below {CLEANUP_FRAGMENTS} PCR fragments" in purify
    assert steps[COLONY_PCR_TITLE].gels == (made.colony.gel,)
    screen = " ".join(steps[COLONY_PCR_TITLE].notes)
    assert f"{SCREENED_COLONIES} of {SCREENED_COLONIES} correct at two fragments" in screen
    assert f"{CORRECT_AT_FIVE} of {SCREENED_COLONIES} at five" in screen
    confirm = " ".join(steps[SEQUENCING_TITLE].notes)
    assert f"one error per {MOLECULES_PER_ERROR} molecules" in confirm
    citations = " ".join(one.text for one in protocol.references)
    assert "In-Fusion Cloning FAQs" in citations
    assert "Gibson, D.G." in citations
    for title in (COLONY_PCR_TITLE, SEQUENCING_TITLE):
        assert steps[title].troubleshooting


def test_the_screening_steps_are_the_shared_builders_and_not_a_second_copy(made):
    # `liulab_mbio.bench` is the only place these are built; this method is its second consumer.
    assert isinstance(made.colony, ColonyCheck)
    assert all(isinstance(read, SangerRead) for read in made.reads)
    steps = {step.title: step for step in made.protocol().steps}
    assert steps["Measure every concentration"] == quantify_step(made.amounts)


def test_the_plans_status_is_the_worst_of_its_checks(made):
    names = [one.name for one in made.checks]
    assert names == [
        "pUC19 backbone",
        "GFP",
        "junctions",
        "overlap length",
        "overlap tm",
        "overlap gc",
        "overlap hairpin",
        "overlap repeat",
        "overlap similarity",
        "fragment count",
        "assembly DNA",
        "primers",
    ]
    assert made.assembly.status == "pass"
    # Both AT-rich ends of GFP warn on GC wherever their placement allows, and nothing fails.
    assert made.status == "warn"
    assert {one.status for one in made.checks} == {"pass", "warn", None}
    # The two nobody has measured carry no verdict at all rather than a pass.
    assert [one.name for one in made.checks if one.status is None] == [
        "overlap gc",
        "overlap similarity",
    ]


def test_the_reaction_and_the_incubation_come_from_the_products_own_tier(made):
    assert made.product is NEBUILDER_HIFI
    assert [one.name for one in made.amounts] == [part.name for part in made.parts]
    program = next(
        one
        for step in made.protocol().steps
        for one in step.programs
        if one.title.startswith("Assembly")
    )
    assert program.stages[0].incubations[0].temperature_c == 50.0
    assert program.stages[0].incubations[0].seconds == 900


def test_the_four_outputs_land_in_the_directory_the_caller_names(made, tmp_path):
    outputs = made.write(tmp_path / "run")
    assert [(path.parent, path.name) for path in outputs.paths] == [
        (tmp_path / "run", "product.dna"),
        (tmp_path / "run", "primers.tsv"),
        (tmp_path / "run", "protocol.json"),
        (tmp_path / "run", "protocol.html"),
    ]
    assert all(path.stat().st_size > 0 for path in outputs.paths)
    assert read_dna(outputs.product) == made.plasmid
    assert read_protocol(outputs.protocol_data) == made.protocol()
    assert outputs.protocol.read_text(encoding="utf-8") == render_html(
        read_protocol(outputs.protocol_data)
    )


def test_the_same_inputs_write_the_same_bytes(made, puc19, gfp, tmp_path):
    first = made.write(tmp_path / "one")
    second = plan_gibson(puc19, gfp).write(tmp_path / "two")
    for one, other in zip(first.paths, second.paths, strict=True):
        assert one.read_bytes() == other.read_bytes()


def test_the_primer_sheet_carries_every_oligo_and_the_page_lists_the_same_rows(made):
    rows = primer_sheet(made.reports).splitlines()
    assert rows[0].split("\t") == ["name", "sequence", "length", "tm_c"]
    # Four amplification primers, three colony PCR primers and two sequencing primers.
    assert len(rows) == 1 + len(made.reports) == 10
    for row, oligo in zip(rows[1:], made.protocol().oligos, strict=True):
        name, sequence, length, tm = row.split("\t")
        assert (oligo.name, oligo.sequence) == (name, sequence)
        assert len(oligo.sequence) == int(length)
        assert oligo.tm_c == pytest.approx(float(tm), abs=0.05)
        assert oligo.purpose


def test_the_protocol_is_enough_to_run_the_experiment(made):
    protocol = made.protocol()
    assert [step.title for step in protocol.steps] == [
        "Amplify pUC19 backbone",
        "Amplify GFP",
        "Check the PCRs on a gel",
        "Digest the plasmid template with DpnI",
        "Purify every amplicon",
        "Measure every concentration",
        "Set up the NEBuilder HiFi DNA Assembly Master Mix reaction",
        "Incubate the assembly",
        "Transform and plate",
        "Screen colonies by PCR",
        "Confirm the clone by sequencing",
    ]
    for step in protocol.steps:
        assert step.expected, step.title
    assert all(len(value) <= OVERVIEW_CHARS for value in protocol.overview.values())
    # Every check reaches the page, the ones nothing judges included, with their detail.
    assert [(one.name, one.status, one.detail) for one in protocol.checks] == [
        (one.name, one.status, one.detail) for one in made.checks
    ]
    assert [one.name for one in protocol.checks if one.status is None] == [
        "overlap gc",
        "overlap similarity",
    ]


def test_the_dpni_digest_prints_nebs_own_dose_and_its_heat_inactivation(made):
    step = next(one for one in made.protocol().steps if one.title.endswith("DpnI"))
    instructions = " ".join(step.instructions)
    assert "20 units" in instructions
    assert "37 °C for 30 minutes" in instructions
    assert "80 °C for 20 minutes" in instructions
    assert [timer.label for timer in step.timers] == ["DpnI digest", "Heat-inactivate DpnI"]


def test_the_protocol_cites_the_documents_behind_its_numbers(made):
    citations = " ".join(reference.text for reference in made.protocol().references)
    assert "NEBuilder" in citations
    assert "protocols.io" in citations
    assert "REBASE" in citations


def test_a_plan_needs_at_least_one_insert(puc19):
    with pytest.raises(ValueError, match="at least one insert"):
        plan_gibson(puc19)


def test_there_is_one_orientation_for_each_insert(puc19, gfp):
    with pytest.raises(ValueError, match="2 values for 1 insert"):
        plan_gibson(puc19, gfp, orientation=["forward", "reverse"])


@pytest.fixture(scope="module")
def handed_in(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP the other way round, into a backbone that was cut last week and handed in as it is.

    Two things that do not touch each other, designed once: a plan whose product is new to the
    suite pays a second of primer design, and this is the only one here that has to.
    """
    cut = SequenceRecord(puc19.sequence[MCS[1] :] + puc19.sequence[: MCS[0]], name="pUC19 cut")
    return plan_gibson(cut, gfp, orientation="reverse")


def test_a_vector_already_linear_is_taken_as_the_opened_part(handed_in, puc19, gfp):
    backbone = handed_in.linearised_vector
    # No PCR opens it, so it has no primers, nothing to run on the gel and no template to cut.
    assert not backbone.amplified
    assert (backbone.forward, backbone.reverse, backbone.report) == (None, None, None)
    assert not backbone.dpni
    designed = [one.primer.name for one in handed_in.reports]
    assert designed[:2] == ["GFP forward", "GFP reverse"]
    assert not any(name.startswith("pUC19 cut backbone") for name in designed)
    titles = [step.title for step in handed_in.protocol().steps]
    assert titles[:2] == ["Amplify GFP", "Check the PCRs on a gel"]
    assert not any(title.endswith("DpnI") for title in titles)
    assert "DpnI" not in [one.name for one in handed_in.protocol().materials]
    # It still assembles into the same plasmid, with its two ends spelling the two overlaps.
    assert (
        len(handed_in.plasmid)
        == len(handed_in.vector) + len(gfp)
        == len(puc19) - (MCS[1] - MCS[0]) + len(gfp)
    )
    assert handed_in.span == (len(handed_in.vector), len(handed_in.vector))
    assert {one.taken_from for one in handed_in.assembly.junctions} == {backbone.name}
    # The inserts run from the end of the vector's own bases, across the origin.
    assert handed_in.assembly.insert_span == (len(handed_in.vector), len(handed_in.plasmid))


def test_an_insert_goes_in_on_the_strand_it_is_given_for(handed_in, gfp):
    assert handed_in.inserts[0].sequence == reverse_complement(gfp.sequence)
    assert handed_in.inserts[0].sequence in handed_in.plasmid.sequence
    assert gfp.sequence not in handed_in.plasmid.sequence


def test_a_site_is_refused_for_a_vector_that_is_already_linear(gfp):
    with pytest.raises(ValueError, match="already linear"):
        plan_gibson(gfp, gfp, site="MCS")


def test_a_vector_with_no_multiple_cloning_site_asks_where_to_put_the_insert(gfp):
    bare = SequenceRecord("ACGT" * 400, topology="circular", name="bare")
    with pytest.raises(ValueError, match="insertion site"):
        plan_gibson(bare, gfp)


def test_the_insertion_site_is_read_from_a_feature_name_or_given_as_coordinates(puc19, gfp):
    named = plan_gibson(puc19, gfp, site="MCS")
    given = plan_gibson(puc19, gfp, site=MCS)
    assert named.plasmid.sequence == given.plasmid.sequence
    assert named.span == given.span == MCS


def test_the_files_are_read_from_disk_when_a_path_is_given(puc19_file, gfp_file):
    assert plan_gibson(puc19_file, gfp_file).plasmid.name == "pUC19-GFP"


@pytest.fixture(scope="module")
def halves(gfp: SequenceRecord) -> tuple[SequenceRecord, SequenceRecord]:
    """GFP cut in two at `SPLIT`, to go in as two inserts."""
    return (
        SequenceRecord(gfp.sequence[:SPLIT], name="GFP 5'"),
        SequenceRecord(gfp.sequence[SPLIT:], name="GFP 3'"),
    )


@pytest.fixture(scope="module")
def several(puc19: SequenceRecord, halves: tuple[SequenceRecord, ...]) -> Plan:
    """The same GFP in two pieces, which is the commonest reason to join more than one insert."""
    return plan_gibson(puc19, *halves)


def test_the_inserts_go_round_the_product_in_the_order_they_are_given(several, made, puc19, gfp):
    assert [part.name for part in several.parts] == ["pUC19 backbone", "GFP 5'", "GFP 3'"]
    # Two inserts in the order given spell what the one they were cut from spells.
    assert several.plasmid.sequence == made.plasmid.sequence
    joined = [(one.before, one.after, one.taken_from) for one in several.assembly.junctions]
    assert joined == [
        ("pUC19 backbone", "GFP 5'", "pUC19 backbone"),
        ("GFP 5'", "GFP 3'", "GFP 5'"),
        ("GFP 3'", "pUC19 backbone", "pUC19 backbone"),
    ]
    # The two outer junctions are the vector's own bases, as NEB asks for a reusable backbone;
    # the inner one is the bases of the insert before it, which the next insert's primer tails.
    first, inner, last = several.assembly.junctions
    assert first.overlap == puc19.sequence[MCS[0] - first.length : MCS[0]]
    assert last.overlap == puc19.sequence[MCS[1] : MCS[1] + last.length]
    assert inner.overlap == gfp.sequence[SPLIT - inner.length : SPLIT]
    assert several.insert_parts[1].forward.sequence.startswith(inner.overlap)
    assert several.linearised_vector.left_tail == several.linearised_vector.right_tail == ""


def test_the_product_names_every_junction_and_the_reaction_counts_every_fragment(several):
    marks = sorted(one.name for one in several.plasmid.features if one.name.endswith("overlap"))
    assert marks == [
        "GFP 3'-pUC19 backbone overlap",
        "GFP 5'-GFP 3' overlap",
        "pUC19 backbone-GFP 5' overlap",
    ]
    assert [one.name for one in several.amounts] == [part.name for part in several.parts]
    judged = {one.name: one for one in several.checks}
    assert judged["fragment count"].value == 2
    assert judged["overlap similarity"].status is None


def test_the_assembly_product_the_caller_names_sets_the_design_and_the_reaction(puc19, gfp):
    made = plan_gibson(puc19, gfp, product=IN_FUSION)
    assert made.product is IN_FUSION
    # Takara's own overlap: 15 bp for one insert, and no melting temperature judges it.
    assert [len(one) for one in made.overlaps] == [15, 15]
    judged = {one.name: one for one in made.checks}
    assert judged["overlap tm"].status is None
    assert judged["assembly DNA"].status is None
    table = next(
        one for step in made.protocol().steps for one in step.tables if "reaction" in one.title
    )
    volumes = {one.name: one.volume_ul for one in table.components}
    assert sum(volumes.values()) == pytest.approx(10.0)
    assert volumes[IN_FUSION.name] == 2.0
    # Its picomoles and its weight are both printed, at the ratio Takara asks for.
    insert = next(one for one in table.components if one.name == "GFP")
    assert "pmol" in insert.final
    assert "ng" in insert.final
    assert made.amounts[1].pmol == pytest.approx(2.0 * made.amounts[0].pmol, rel=1e-3)


@pytest.fixture(scope="module")
def routed(puc19: SequenceRecord, halves: tuple[SequenceRecord, ...]) -> Plan:
    """The two inserts of `several`, the second stitched and the junction between them bridged.

    Every sequence is one that plan has already primed, so the two differ in route and in
    nothing else.
    """
    return plan_gibson(puc19, *halves, route=["amplify", "stitch"], bridge=[("GFP 5'", "GFP 3'")])


def test_a_stitched_part_is_oligos_tiling_both_strands_rather_than_a_pcr(routed):
    stitched = routed.insert_parts[1]
    assert stitched.stitched
    assert not stitched.amplified
    assert (stitched.forward, stitched.reverse, stitched.report, stitched.dpni) == (
        None,
        None,
        None,
        False,
    )
    # As few oligos as tile the whole molecule the reaction assembles, none longer than the
    # note's 60 bases, neighbours overlapping by exactly its 20 bp, and the strands alternating.
    tiled = stitched.length - STITCH_OVERLAP_BP
    step = STITCH_OLIGO_BASES - STITCH_OVERLAP_BP
    assert len(stitched.oligos) == -(-tiled // step)
    assert max(len(one.sequence) for one in stitched.oligos) <= STITCH_OLIGO_BASES
    assert [one.strand for one in stitched.oligos] == [Strand.FORWARD, Strand.REVERSE] * (
        len(stitched.oligos) // 2
    ) + [Strand.FORWARD] * (len(stitched.oligos) % 2)
    assert (stitched.oligos[0].start, stitched.oligos[-1].end) == (0, stitched.length)
    for one in stitched.oligos:
        written = one.sequence if one.strand == Strand.FORWARD else reverse_complement(one.sequence)
        assert written == stitched.amplicon.sequence[one.start : one.end]
    # Nothing amplifies it, so no PCR step and no gel lane names it.
    steps = {step.title: step for step in routed.protocol().steps}
    assert f"Amplify {stitched.name}" not in steps
    assert stitched.name not in [
        lane.label for gel in steps["Check the PCRs on a gel"].gels for lane in gel.lanes
    ]


def test_a_bridging_oligo_joins_two_fragments_and_neither_carries_a_tail(routed):
    first, second = routed.insert_parts
    bridged = next(one for one in routed.assembly.junctions if one.bridge)
    assert (bridged.before, bridged.after) == (first.name, second.name)
    assert (first.right_tail, second.left_tail) == ("", "")
    # It spells the end of one fragment and the start of the next, and nothing goes between them.
    oligo = next(one for one in routed.ordered_oligos if one.name == bridged.bridge)
    assert oligo.sequence == first.bases[-BRIDGE_HOMOLOGY_BP:] + second.bases[:BRIDGE_HOMOLOGY_BP]
    assert routed.plasmid.sequence.count(oligo.sequence) == 1
    # A map marks it as a bridge rather than as an overlap the two share.
    mark = next(
        one for one in routed.plasmid.features if one.name == f"{first.name}-{second.name} bridge"
    )
    assert (mark.segments[0].start, mark.segments[0].end) == (
        bridged.start,
        bridged.end + BRIDGE_HOMOLOGY_BP,
    )


def test_the_product_and_the_validation_are_the_same_whichever_route_made_a_part(routed, several):
    assert routed.plasmid.sequence == several.plasmid.sequence
    assert [(one.before, one.after) for one in routed.assembly.junctions] == [
        (one.before, one.after) for one in several.assembly.junctions
    ]
    assert routed.assembly.boundaries == several.assembly.boundaries
    assert routed.assembly.insert_span == several.assembly.insert_span
    assert [one.bands_bp for one in routed.colony.clones] == [
        one.bands_bp for one in several.colony.clones
    ]
    assert [one.read_bp for one in routed.reads] == [one.read_bp for one in several.reads]


def test_both_routes_put_their_oligos_on_the_one_sheet_with_no_verdict_and_why(routed):
    rows = primer_sheet(routed.reports, oligos=routed.ordered_oligos).splitlines()
    assert len(rows) == 1 + len(routed.reports) + len(routed.ordered_oligos)
    page = routed.protocol().oligos
    assert [one.name for one in page] == [row.split("\t")[0] for row in rows[1:]]
    ordered = page[len(routed.reports) :]
    assert [one.purpose for one in ordered] == ["Stitch GFP 3'"] * (len(ordered) - 1) + [
        "Bridge GFP 5' to GFP 3'"
    ]
    for oligo, row in zip(ordered, rows[1 + len(routed.reports) :], strict=True):
        assert row.split("\t") == [oligo.name, oligo.sequence, str(len(oligo.sequence)), ""]
        # Nothing measures an oligo that primes nothing, so the row says so and carries no verdict.
        assert (oligo.status, oligo.tm_c, oligo.checks) == (None, None, ())
        assert "carries no verdict" in oligo.note
        assert "desalted" in oligo.note
        assert oligo.stock


def test_an_oligo_with_no_verdict_survives_the_protocol_being_written_and_read_again(
    routed, tmp_path
):
    written = routed.write(tmp_path / "routed")
    assert read_protocol(written.protocol_data) == routed.protocol()
    rows = written.primers.read_text(encoding="utf-8").splitlines()
    assert rows[-1].split("\t")[0].endswith("bridge")


def test_the_assembly_step_doses_each_route_the_way_its_own_source_does(routed):
    step = next(one for one in routed.protocol().steps if one.title.startswith("Set up"))
    said = " ".join(step.notes)
    assert "45 nM of each" in said
    assert "no separate annealing step" in said.lower()
    assert "1 pmol" in said
    # The stitched part is still one of the fragments the table gives picomoles for.
    assert [one.name for one in routed.amounts] == [part.name for part in routed.parts]


def test_a_part_above_addgenes_window_is_told_and_still_stitched(routed):
    judged = {one.name: one for one in routed.checks}
    # Refused only above the oligo count Gibson states; above Addgene's window it is only told.
    assert judged["stitched size"].status == "warn"
    assert f"GFP 3' {routed.insert_parts[1].fragment_length} bp" in judged["stitched size"].detail
    assert "Addgene uses this route between 60 and 150 bp" in judged["stitched size"].detail


def test_a_product_supporting_neither_oligo_route_refuses_both_and_says_which_do(puc19, gfp):
    for asked in (
        lambda: plan_gibson(puc19, gfp, product=IN_FUSION, route="stitch"),
        lambda: plan_gibson(puc19, gfp, product=IN_FUSION, bridge=[("GFP", "pUC19 backbone")]),
    ):
        with pytest.raises(ValueError, match="no single-stranded oligo") as raised:
            asked()
        assert NEBUILDER_HIFI.name in str(raised.value)


def test_a_bridge_naming_a_junction_this_assembly_has_not_got_is_refused(puc19, gfp):
    with pytest.raises(ValueError, match="no junction joins"):
        plan_gibson(puc19, gfp, bridge=[("GFP", "linker")])


def test_there_is_one_route_for_each_insert(puc19, gfp):
    with pytest.raises(ValueError, match="2 values for 1 insert"):
        plan_gibson(puc19, gfp, route=["amplify", "stitch"])
