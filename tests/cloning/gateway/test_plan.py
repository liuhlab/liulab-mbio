"""`plan_gateway`, the entry point: the files, the status, the product and the protocol."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from liulab_mbio.bench.gels import agarose_percent
from liulab_mbio.bench.steps import COLONY_PCR_TITLE, SEQUENCING_TITLE
from liulab_mbio.bench.validation import CORRECT_CLONE, EMPTY_CLONE, SANGER_FLANK
from liulab_mbio.cloning.gateway import plan_gateway
from liulab_mbio.cloning.gateway.att import CROSSOVER, REGIONS, find_att_sites
from liulab_mbio.cloning.gateway.bench import (
    BP_CELSIUS,
    BP_VOLUME_UL,
    ENTRY_NG,
    LR_CELSIUS,
    LR_VOLUME_UL,
)
from liulab_mbio.cloning.gateway.design import attb_tails
from liulab_mbio.io import read_record
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
    from liulab_mbio.cloning.gateway import Plan


def test_the_plan_writes_the_product_the_sheet_and_the_protocol_pair(
    gateway_plan: Plan, tmp_path: Path
) -> None:
    written = gateway_plan.write(tmp_path / "run")

    assert [path.name for path in written.paths] == [
        "product.dna",
        "primers.tsv",
        "protocol.json",
        "protocol.html",
    ]
    assert all(path.exists() for path in written.paths)
    assert written.entry is None
    assert read_record(written.product).sequence == gateway_plan.product.sequence


def test_the_same_inputs_write_the_same_bytes(
    entry: SequenceRecord, destination: SequenceRecord, tmp_path: Path
) -> None:
    first = plan_gateway(entry, destination).write(tmp_path / "one")
    second = plan_gateway(entry, destination).write(tmp_path / "two")

    for one, other in zip(first.paths, second.paths, strict=True):
        assert one.read_bytes() == other.read_bytes()


def test_the_status_is_the_worst_of_the_checks_it_carries(gateway_plan: Plan) -> None:
    assert {check.name for check in gateway_plan.checks} == {
        "primers",
        "LR junctions",
        "LR att sites",
        "insert att sites",
        "ccdB host",
        "ccdB vector host",
        "LR markers",
    }
    assert gateway_plan.status == "pass"
    assert all(check.status in ("pass", None) for check in gateway_plan.checks)


def test_each_junction_spells_the_att_site_the_arithmetic_names(gateway_plan: Plan) -> None:
    product = gateway_plan.product
    first, second = gateway_plan.junctions

    assert (first.name, second.name) == ("attB1", "attB2")
    assert first.bases == REGIONS["attB1"]
    assert second.bases == REGIONS["attB2"]
    assert product.extract(first.span) == REGIONS["attB1"]
    assert reverse_complement(product.extract(second.span)) == REGIONS["attB2"]


def test_the_product_carries_the_insert_between_the_two_att_sites(
    gateway_plan: Plan, gfp: SequenceRecord
) -> None:
    product = gateway_plan.product
    first, second = gateway_plan.junctions

    assert product.topology == "circular"
    assert product.extract(Segment(first.end, second.start)) == gfp.sequence
    assert [one.name for one in find_att_sites(product)] == ["attB1", "attB2"]


def test_the_expression_clone_is_the_backbone_and_the_insert_and_nothing_else(
    gateway_plan: Plan, destination: SequenceRecord
) -> None:
    made = gateway_plan.lr.recombination
    dropped = destination.extract(Segment(*made.cassette))

    assert len(gateway_plan.product) == made.backbone.length + made.moved.length
    assert made.backbone.length + len(dropped) == len(destination)
    assert dropped not in gateway_plan.product.sequence
    assert "ccdB" not in {feature.name for feature in gateway_plan.product.features}


def test_every_feature_of_both_records_is_carried_to_its_new_coordinates(
    gateway_plan: Plan,
) -> None:
    product = gateway_plan.product
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
    carried = plan_gateway(
        SequenceRecord(
            entry.sequence,
            topology="circular",
            name=entry.name,
            features=entry.features,
            primers=(kept, dropped),
        ),
        destination,
    ).product

    annotated = {primer.name: primer for primer in carried.primers}
    assert set(annotated) == {"reads the insert"}
    site = annotated["reads the insert"].binding_sites[0]
    assert carried.sequence[site.start : site.end] == kept.sequence


def test_the_protocol_carries_the_reaction_its_incubation_its_stop_and_the_plating(
    gateway_plan: Plan,
) -> None:
    protocol = gateway_plan.protocol()
    titles = [step.title for step in protocol.steps]

    assert titles == [
        "Set up the LR reaction",
        "Run the LR reaction",
        "Stop the LR reaction with proteinase K",
        "Transform and plate the LR reaction",
        COLONY_PCR_TITLE,
        SEQUENCING_TITLE,
    ]
    table = protocol.steps[0].tables[0]
    assert sum(component.volume_ul for component in table.components) == LR_VOLUME_UL
    assert [component.name for component in table.components][-2:] == [
        "TE buffer, pH 8.0",
        "Gateway LR Clonase II enzyme mix",
    ]
    assert f"{LR_CELSIUS:g} °C for 1 hour" in " ".join(protocol.steps[1].instructions)
    assert "37 °C for 10 minutes" in " ".join(protocol.steps[2].instructions)
    plating = " ".join(protocol.steps[3].instructions)
    assert "Add 1 µL of the LR reaction" in plating
    assert "Spread 100 µL of a 1:10 dilution" in plating
    assert any("ampicillin" in material.name for material in protocol.materials)


def test_the_protocol_says_what_each_junction_spells_and_that_it_is_not_scarless(
    gateway_plan: Plan,
) -> None:
    protocol = gateway_plan.protocol()

    assert any("not scarless" in line for line in protocol.highlights)
    expected = " ".join(protocol.steps[1].expected)
    assert REGIONS["attB1"] in expected
    assert REGIONS["attB2"] in expected


def test_the_protocol_is_read_back_from_the_json_it_writes(
    gateway_plan: Plan, tmp_path: Path
) -> None:
    written = gateway_plan.write(tmp_path / "run")

    data = json.loads(written.protocol_data.read_text(encoding="utf-8"))
    page = written.protocol.read_text(encoding="utf-8")
    assert all(step["title"] in page for step in data["steps"])
    assert f"{ENTRY_NG:g} ng per reaction" in page


def test_a_record_carrying_no_att_site_is_refused_with_what_was_looked_for(
    destination: SequenceRecord, puc19: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="looked for attL1 and attL2 on either strand"):
        plan_gateway(puc19, destination)


def test_a_destination_vector_carrying_no_att_site_is_refused(
    entry: SequenceRecord, puc19: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="looked for attR1 and attR2"):
        plan_gateway(entry, puc19)


def test_an_entry_clone_read_on_the_other_strand_plans_the_same_product(
    entry: SequenceRecord, destination: SequenceRecord, gateway_plan: Plan
) -> None:
    turned = SequenceRecord(
        reverse_complement(entry.sequence), topology="circular", name=entry.name
    )

    assert plan_gateway(turned, destination).product.sequence == gateway_plan.product.sequence


def test_a_linear_destination_vector_is_refused(
    entry: SequenceRecord, destination: SequenceRecord
) -> None:
    opened = SequenceRecord(destination.sequence, name=destination.name)

    with pytest.raises(ValueError, match="must be circular"):
        plan_gateway(entry, opened)


def test_an_entry_clone_carrying_the_same_site_twice_is_refused(
    entry: SequenceRecord, destination: SequenceRecord
) -> None:
    repeated = SequenceRecord(
        entry.sequence + att_site("attL1"), topology="circular", name=entry.name
    )

    with pytest.raises(ValueError, match="carries 2 attL1 sites, needing one"):
        plan_gateway(repeated, destination)


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

    made = plan_gateway(entry, changed)

    assert made.junctions[0].bases == "T" + REGIONS["attB1"][1:]
    assert made.status == "pass"
    assert gfp.sequence in made.product.sequence


def test_the_bp_route_writes_the_entry_clone_as_a_file_of_its_own(
    staged_plan: Plan, tmp_path: Path
) -> None:
    written = staged_plan.write(tmp_path / "run")

    assert [path.name for path in written.paths] == [
        "entry-clone.dna",
        "product.dna",
        "primers.tsv",
        "protocol.json",
        "protocol.html",
    ]
    assert written.entry is not None
    assert read_record(written.entry).sequence == staged_plan.entry.sequence


def test_the_route_names_which_reactions_were_planned(
    gateway_plan: Plan, staged_plan: Plan
) -> None:
    assert (gateway_plan.route, staged_plan.route) == ("LR", "BP then LR")


def test_the_staged_inputs_write_the_same_bytes(
    insert: SequenceRecord, destination: SequenceRecord, donor: SequenceRecord, tmp_path: Path
) -> None:
    first = plan_gateway(insert, destination, donor=donor).write(tmp_path / "one")
    second = plan_gateway(insert, destination, donor=donor).write(tmp_path / "two")

    for one, other in zip(first.paths, second.paths, strict=True):
        assert one.read_bytes() == other.read_bytes()


def test_bp_writes_an_entry_clone_carrying_the_insert_between_its_attl_sites(
    staged_plan: Plan, gfp: SequenceRecord
) -> None:
    bp = staged_plan.bp
    assert bp is not None
    first, second = bp.junctions

    assert (first.name, second.name) == ("attL1", "attL2")
    assert first.bases == REGIONS["attL1"]
    assert second.bases == REGIONS["attL2"]
    assert bp.product.topology == "circular"
    assert bp.product.extract(Segment(first.end, second.start)) == gfp.sequence
    assert [one.name for one in find_att_sites(bp.product)] == ["attL1", "attL2"]
    assert "GGGG" + REGIONS["attB1"] not in bp.product.sequence


def test_lr_gives_back_the_attb_sites_the_insert_carried_in(
    staged_plan: Plan, gfp: SequenceRecord
) -> None:
    product = staged_plan.product
    first, second = staged_plan.junctions

    assert (first.name, second.name) == ("attB1", "attB2")
    assert product.extract(first.span) == REGIONS["attB1"]
    assert reverse_complement(product.extract(second.span)) == REGIONS["attB2"]
    assert product.extract(Segment(first.end, second.start)) == gfp.sequence
    assert "ccdB" not in {feature.name for feature in product.features}


def test_each_reaction_carries_its_own_verdicts_and_the_plan_its_own(staged_plan: Plan) -> None:
    assert [check.name for check in staged_plan.checks] == [
        "primers",
        "BP junctions",
        "BP att sites",
        "LR junctions",
        "LR att sites",
        "insert att sites",
        "ccdB host",
        "ccdB vector host",
        "LR markers",
    ]
    assert staged_plan.status == "pass"


def test_the_protocol_reads_as_two_staged_reactions_with_a_miniprep_between_them(
    staged_plan: Plan,
) -> None:
    protocol = staged_plan.protocol()
    steps = {step.title: step for step in protocol.steps}

    assert [step.title for step in protocol.steps] == [
        "Set up the BP reaction",
        "Run the BP reaction",
        "Stop the BP reaction with proteinase K",
        "Transform and plate the BP reaction",
        "Pick and miniprep the entry clone",
        "Set up the LR reaction",
        "Run the LR reaction",
        "Stop the LR reaction with proteinase K",
        "Transform and plate the LR reaction",
        COLONY_PCR_TITLE,
        SEQUENCING_TITLE,
    ]
    table = steps["Set up the BP reaction"].tables[0]
    assert table.title == "BP reaction"
    assert sum(component.volume_ul for component in table.components) == BP_VOLUME_UL
    assert [component.name for component in table.components][-2:] == [
        "TE buffer, pH 8.0",
        "Gateway BP Clonase II enzyme mix",
    ]
    assert f"{BP_CELSIUS:g} °C for 1 hour" in " ".join(steps["Run the BP reaction"].instructions)
    stop = " ".join(steps["Stop the BP reaction with proteinase K"].instructions)
    assert "37 °C for 10 minutes" in stop
    assert any(
        "not the stopped BP reaction" in note
        for note in steps["Pick and miniprep the entry clone"].notes
    )


def test_a_donor_vector_carrying_no_att_site_is_refused(
    insert: SequenceRecord, destination: SequenceRecord, puc19: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="looked for attP1 and attP2 on either strand"):
        plan_gateway(insert, destination, donor=puc19)


def test_an_insert_carrying_no_attb_end_is_refused(
    gfp: SequenceRecord, destination: SequenceRecord, donor: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="looked for attB1 and attB2 on either strand"):
        plan_gateway(gfp, destination, donor=donor)


def test_a_linear_donor_vector_is_refused(
    insert: SequenceRecord, destination: SequenceRecord, donor: SequenceRecord
) -> None:
    opened = SequenceRecord(donor.sequence, name=donor.name)

    with pytest.raises(ValueError, match="must be circular"):
        plan_gateway(insert, destination, donor=opened)


def test_a_plain_insert_is_amplified_onto_attb_ends_and_the_sheet_is_written(
    amplified_plan: Plan, gfp: SequenceRecord, tmp_path: Path
) -> None:
    made = amplified_plan.amplicon
    assert made is not None
    written = amplified_plan.write(tmp_path / "run")

    assert [path.name for path in written.paths] == [
        "entry-clone.dna",
        "product.dna",
        "primers.tsv",
        "protocol.json",
        "protocol.html",
    ]
    assert written.primers is not None
    assert made.record.sequence.startswith(attb_tails()[0] + gfp.sequence)
    assert [one.name for one in find_att_sites(made.record)] == ["attB1", "attB2"]


def test_the_entry_clone_is_the_same_whether_the_insert_arrived_attb_flanked_or_plain(
    amplified_plan: Plan, staged_plan: Plan
) -> None:
    assert amplified_plan.entry.sequence == staged_plan.entry.sequence
    assert amplified_plan.product.sequence == staged_plan.product.sequence
    assert "GGGG" + REGIONS["attB1"] not in amplified_plan.entry.sequence


def test_the_order_sheet_carries_the_name_sequence_length_and_tm_of_every_oligo(
    amplified_plan: Plan, tmp_path: Path
) -> None:
    written = amplified_plan.write(tmp_path / "run")
    assert written.primers is not None

    header, *rows = written.primers.read_text(encoding="utf-8").splitlines()
    assert header.split("\t") == ["name", "sequence", "length", "tm_c"]
    assert len(rows) == len(amplified_plan.reports)
    for row, report in zip(rows, amplified_plan.reports, strict=True):
        name, sequence, length, tm = row.split("\t")
        assert (name, sequence) == (report.primer.name, report.primer.sequence)
        assert int(length) == len(sequence)
        assert float(tm) == pytest.approx(report["tm"].value, abs=0.05)


def test_the_primers_are_judged_and_warn_only_where_nothing_allowed_does_better(
    amplified_plan: Plan,
) -> None:
    first, *_ = amplified_plan.checks
    verdict = next(check for check in amplified_plan.checks if check.name == "primers")

    assert first.name == "primers"
    assert amplified_plan.status == "warn"
    assert all(report.status != "fail" for report in amplified_plan.reports)
    assert verdict.value == 6


def test_the_protocol_carries_the_pcr_step_its_program_and_its_expected_band(
    amplified_plan: Plan,
) -> None:
    made = amplified_plan.amplicon
    assert made is not None
    protocol = amplified_plan.protocol()

    assert [step.title for step in protocol.steps][:2] == [
        f"Amplify {made.name}",
        "Purify every amplicon",
    ]
    step = protocol.steps[0]
    assert step.programs[0].title == f"{made.name} PCR"
    assert f"One band at {made.length} bp" in " ".join(step.expected)
    assert any(made.polymerase.name in component.name for component in step.tables[0].components)
    assert [oligo.name for oligo in protocol.oligos] == [
        report.primer.name for report in amplified_plan.reports
    ]


def test_an_insert_no_primer_can_be_placed_on_is_refused(
    destination: SequenceRecord, donor: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="takes no attB primer"):
        plan_gateway(
            SequenceRecord("ATGCATGCATGC", name="tiny"), destination, donor=donor, amplify=True
        )


def test_amplifying_an_insert_without_a_donor_vector_is_refused(
    gfp: SequenceRecord, destination: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="needs a donor vector as well"):
        plan_gateway(gfp, destination, amplify=True)


def test_the_colony_pcr_reads_from_the_part_boundary_not_the_att_region_edge(
    gateway_plan: Plan,
) -> None:
    made = gateway_plan.lr.recombination
    first, second = gateway_plan.junctions
    start, end = made.boundaries

    assert (start, end) == (first.start + CROSSOVER, second.end - CROSSOVER)
    assert gateway_plan.product.extract(Segment(start, end)) == made.moved.bases


def test_the_colony_pcr_crosses_both_attb_junctions_and_sits_outside_them(
    gateway_plan: Plan,
) -> None:
    colony = gateway_plan.colony
    first, second = gateway_plan.junctions
    forward, reverse = (primer.binding_sites[0] for primer in colony.primers)

    assert [primer.name for primer in colony.primers] == [
        "Colony PCR forward",
        "Colony PCR reverse",
    ]
    assert forward.end <= first.start
    assert reverse.start >= second.end
    assert all(report.status != "fail" for report in colony.reports)


def test_the_gel_carries_a_ladder_a_lane_per_candidate_and_no_reversed_insert(
    gateway_plan: Plan,
) -> None:
    colony = gateway_plan.colony
    correct, empty = colony.clones
    made = gateway_plan.lr.recombination
    cassette = made.cassette[1] - made.cassette[0]

    assert [lane.label for lane in colony.gel.lanes] == [CORRECT_CLONE, EMPTY_CLONE]
    assert colony.gel.ladder.bands_bp
    assert colony.agarose_percent == agarose_percent(correct.bands_bp + empty.bands_bp)
    assert correct.bands_bp[0] - empty.bands_bp[0] == made.moved.length - cassette


def test_the_sequencing_primers_read_in_from_outside_each_junction(gateway_plan: Plan) -> None:
    first, second = gateway_plan.junctions
    across = range(first.start, second.end)

    for read in gateway_plan.reads:
        site = read.primer.binding_sites[0]
        assert site.start not in across
        assert site.end not in across
        assert read.distance_bp >= SANGER_FLANK
        assert read.read_bp >= second.end - first.start


def test_the_protocol_carries_the_colony_pcr_its_program_the_gel_and_the_sequencing(
    gateway_plan: Plan,
) -> None:
    steps = {step.title: step for step in gateway_plan.protocol().steps}
    screen, confirm = steps[COLONY_PCR_TITLE], steps[SEQUENCING_TITLE]

    assert screen.programs
    assert screen.tables
    assert [lane.label for lane in screen.gels[0].lanes] == [CORRECT_CLONE, EMPTY_CLONE]
    assert f"{gateway_plan.colony.agarose_percent:g}% gel" in " ".join(screen.instructions)
    assert not any("reversed" in line.lower() for line in screen.expected)
    assert any("pCR8/GW/TOPO" in note for note in confirm.notes)
    assert REGIONS["attB1"] in " ".join(confirm.expected)


def test_both_sets_of_oligos_reach_the_order_sheet_with_what_each_is_for(
    amplified_plan: Plan,
) -> None:
    purposes = {oligo.name: oligo.purpose for oligo in amplified_plan.protocol().oligos}

    assert [one.role for one in amplified_plan.designed_oligos] == [
        "amplification",
        "amplification",
        "colony PCR",
        "colony PCR",
        "sequencing",
        "sequencing",
    ]
    assert purposes["Colony PCR forward"] == COLONY_PCR_TITLE
    assert purposes["Sequencing reverse"] == SEQUENCING_TITLE
    assert purposes[amplified_plan.designed_oligos[0].report.primer.name].startswith("Amplify")
