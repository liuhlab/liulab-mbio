"""The whole pipeline, run on the fixtures with no agent.

`tests/data/pUC19.dna` and `tests/data/GFP.dna`: put GFP into the pUC19 multiple cloning site by
Gibson assembly. Nothing about the fixtures is hard-coded in the package; the numbers below are
read off the records and pinned here.
"""

import pytest

from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.cloning.gibson import Plan, plan_gibson
from liulab_mbio.cloning.gibson.bench import NEBUILDER_HIFI
from liulab_mbio.cloning.gibson.design import wallace_tm
from liulab_mbio.protocol import OVERVIEW_CHARS, read_protocol, render_html
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.snapgene import read_dna

#: Where the fixture's own MCS feature sits.
MCS = (395, 452)


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
    designed = [one.primer.name for one in made.reports]
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


def test_the_plans_status_is_the_worst_of_its_checks(made):
    names = [one.name for one in made.checks]
    assert names == ["pUC19 backbone", "GFP", "junctions", "primers"]
    assert made.assembly.status == "pass"
    # Both AT-rich ends of GFP warn on GC wherever their placement allows, and nothing fails.
    assert made.status == "warn"
    assert {one.status for one in made.checks} == {"pass", "warn"}


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
    assert len(rows) == 1 + len(made.reports) == 5
    for row, oligo in zip(rows[1:], made.protocol().oligos, strict=True):
        name, sequence, length, tm = row.split("\t")
        assert (oligo.name, oligo.sequence) == (name, sequence)
        assert len(oligo.sequence) == int(length)
        assert oligo.tm_c == pytest.approx(float(tm), abs=0.05)
        assert oligo.purpose.startswith("Amplify")


def test_the_protocol_is_enough_to_run_the_experiment(made):
    protocol = made.protocol()
    assert [step.title for step in protocol.steps] == [
        "Amplify pUC19 backbone",
        "Amplify GFP",
        "Check the PCRs on a gel",
        "Digest the plasmid template with DpnI",
        "Set up the NEBuilder HiFi DNA Assembly Master Mix reaction",
        "Incubate the assembly",
        "Transform and plate",
    ]
    for step in protocol.steps:
        assert step.expected, step.title
    assert all(len(value) <= OVERVIEW_CHARS for value in protocol.overview.values())
    assert [(one.name, one.status) for one in protocol.checks] == [
        (one.name, one.status) for one in made.checks
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


@pytest.mark.parametrize("count", [0, 2])
def test_a_plan_takes_exactly_one_insert(puc19, gfp, count):
    with pytest.raises(ValueError, match="one insert"):
        plan_gibson(puc19, *((gfp,) * count))


def test_a_linear_vector_is_refused_and_says_why(gfp):
    with pytest.raises(ValueError, match="linear"):
        plan_gibson(gfp, gfp)


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
