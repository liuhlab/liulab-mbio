"""The whole pipeline, run on the fixtures with no agent.

`tests/data/pUC19.dna` and `tests/data/GFP.dna`: put GFP into the pUC19 multiple cloning site by
Gibson assembly. Nothing about the fixtures is hard-coded in the package; the numbers below are
read off the records and pinned here.

One plan carries the common path, and one more each route the method takes differently: a
vector handed in linear, an insert the other way round, several inserts, and the two oligo
routes. What a function below the plan decides is checked there instead, in `test_assembly.py`
and `test_design.py`.
"""

import pytest

from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.bench.steps import COLONY_PCR_TITLE, SEQUENCING_TITLE, quantify_step
from liulab_mbio.bench.validation import ColonyCheck, SangerRead
from liulab_mbio.cloning.gibson import Plan, plan_gibson
from liulab_mbio.cloning.gibson.bench import IN_FUSION, NEBUILDER_HIFI
from liulab_mbio.cloning.gibson.design import BRIDGE_HOMOLOGY_BP
from liulab_mbio.cloning.gibson.steps import (
    CLEANUP_FRAGMENTS,
    CORRECT_AT_FIVE,
    MOLECULES_PER_ERROR,
    SCREENED_COLONIES,
)
from liulab_mbio.edits import flipped
from liulab_mbio.protocol import OVERVIEW_CHARS, read_protocol, render_html
from liulab_mbio.sequence import SequenceRecord
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


def test_the_junctions_are_the_vectors_own_bases_and_the_insert_carries_them_as_tails(made, puc19):
    first, last = made.assembly.junctions
    assert [len(one) for one in made.overlaps] == [16, 16]
    # The bases are the vector's own, on either side of the span the insert replaces, so the
    # opened backbone needs no tail and can be reused.
    assert first.overlap == puc19.sequence[MCS[0] - first.length : MCS[0]]
    assert last.overlap == puc19.sequence[MCS[1] : MCS[1] + last.length]
    assert {one.taken_from for one in (first, last)} == {made.linearised_vector.name}
    insert = made.insert_parts[0]
    assert (insert.left_tail, insert.right_tail) == (first.overlap, last.overlap)
    assert (made.linearised_vector.left_tail, made.linearised_vector.right_tail) == ("", "")


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
    for cited in ("In-Fusion Cloning FAQs", "Gibson, D.G.", "NEBuilder", "protocols.io", "REBASE"):
        assert cited in citations
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
    # The reaction and the incubation are the default product's, and every fragment is weighed.
    assert made.product is NEBUILDER_HIFI
    assert [one.name for one in made.amounts] == [part.name for part in made.parts]
    # The template the backbone was amplified off is digested, and the digest is timed.
    digest = next(one for one in protocol.steps if one.title.endswith("DpnI"))
    assert "37 °C for 30 minutes" in " ".join(digest.instructions)
    assert [timer.label for timer in digest.timers] == ["DpnI digest", "Heat-inactivate DpnI"]
    # Every check reaches the page, the ones nothing judges included, with their detail.
    assert [(one.name, one.status, one.detail) for one in protocol.checks] == [
        (one.name, one.status, one.detail) for one in made.checks
    ]


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


def test_a_plan_needs_at_least_one_insert(puc19):
    with pytest.raises(ValueError, match="at least one insert"):
        plan_gibson(puc19)


def test_there_is_one_orientation_for_each_insert(puc19, gfp):
    with pytest.raises(ValueError, match="2 values for 1 insert"):
        plan_gibson(puc19, gfp, orientation=["forward", "reverse"])


def test_an_insert_goes_in_on_the_strand_it_is_given_for(made, puc19, gfp):
    # The other strand of GFP, turned back by the orientation, is the plan the forward record
    # makes: every sequence here is one the run has already designed on.
    turned = plan_gibson(puc19, flipped(gfp), orientation="reverse")
    assert turned.inserts[0].sequence == gfp.sequence
    assert turned.plasmid.sequence == made.plasmid.sequence
    assert flipped(gfp).sequence not in turned.plasmid.sequence


@pytest.fixture(scope="module")
def handed_in(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP into a backbone that was cut last week and handed in as it is.

    The only plan here whose product is new to the run, so the only one paying a second of
    primer design: a linear vector puts a junction at base zero, which no other product has.
    """
    cut = SequenceRecord(puc19.sequence[MCS[1] :] + puc19.sequence[: MCS[0]], name="pUC19 cut")
    return plan_gibson(cut, gfp)


def test_a_vector_already_linear_is_opened_by_nothing_and_the_inserts_cross_the_origin(
    handed_in, puc19, gfp
):
    backbone = handed_in.linearised_vector
    # No PCR opens it, so nothing runs on the gel for it and no template has to be cut.
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
    # The reaction counts the inserts, and every fragment is weighed for it.
    assert {one.name: one.value for one in several.checks}["fragment count"] == 2
    assert [one.name for one in several.amounts] == [part.name for part in several.parts]


@pytest.fixture(scope="module")
def routed(puc19: SequenceRecord, halves: tuple[SequenceRecord, ...]) -> Plan:
    """The two inserts of `several`, the second stitched and the junction between them bridged.

    Every sequence is one that plan has already primed, so the two differ in route and in
    nothing else.
    """
    return plan_gibson(puc19, *halves, route=["amplify", "stitch"], bridge=[("GFP 5'", "GFP 3'")])


def test_a_stitched_part_is_oligos_rather_than_a_pcr_and_is_told_where_it_is_unusual(routed):
    stitched = routed.insert_parts[1]
    assert stitched.stitched
    assert not stitched.amplified
    assert (stitched.forward, stitched.reverse, stitched.report, stitched.dpni) == (
        None,
        None,
        None,
        False,
    )
    # Nothing amplifies it, so no PCR step and no gel lane names it.
    steps = {step.title: step for step in routed.protocol().steps}
    assert f"Amplify {stitched.name}" not in steps
    assert stitched.name not in [
        lane.label for gel in steps["Check the PCRs on a gel"].gels for lane in gel.lanes
    ]
    # Refused only above the oligo count Gibson states; above Addgene's window it is only told.
    judged = {one.name: one for one in routed.checks}
    assert judged["stitched size"].status == "warn"
    assert f"GFP 3' {stitched.fragment_length} bp" in judged["stitched size"].detail


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
