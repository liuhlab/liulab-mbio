"""`plan_gateway`, the entry point: the files, the status, the product and the protocol."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from liulab_mbio.cloning.gateway import plan_gateway
from liulab_mbio.cloning.gateway.att import REGIONS, find_att_sites
from liulab_mbio.cloning.gateway.bench import ENTRY_NG, LR_CELSIUS, LR_VOLUME_UL
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


def test_the_plan_writes_the_product_and_the_protocol_pair(
    gateway_plan: Plan, tmp_path: Path
) -> None:
    written = gateway_plan.write(tmp_path / "run")

    assert [path.name for path in written.paths] == [
        "product.dna",
        "protocol.json",
        "protocol.html",
    ]
    assert all(path.exists() for path in written.paths)
    assert read_record(written.product).sequence == gateway_plan.product.sequence


def test_the_same_inputs_write_the_same_bytes(
    entry: SequenceRecord, destination: SequenceRecord, tmp_path: Path
) -> None:
    first = plan_gateway(entry, destination).write(tmp_path / "one")
    second = plan_gateway(entry, destination).write(tmp_path / "two")

    for one, other in zip(first.paths, second.paths, strict=True):
        assert one.read_bytes() == other.read_bytes()


def test_the_status_is_the_worst_of_the_checks_it_carries(gateway_plan: Plan) -> None:
    assert {check.name for check in gateway_plan.checks} == {"junctions", "att sites"}
    assert gateway_plan.status == "pass"
    assert all(check.status == "pass" for check in gateway_plan.checks)


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
    made = gateway_plan.recombination
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
        "Stop the reaction with proteinase K",
        "Transform and plate",
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
