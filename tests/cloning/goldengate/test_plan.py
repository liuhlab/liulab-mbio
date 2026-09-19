"""The whole pipeline, run on the fixtures with no agent.

`tests/data/pUC19.dna` and `tests/data/GFP.dna`: put GFP into the pUC19 multiple cloning site
and validate the insertion by colony PCR. Nothing about the fixtures is hard-coded in the
package; the numbers below are read off the records and pinned here.

Two plans cover the pipeline: the shared `plan`, which is the common path, and `four`, the
many-part case with one insert the other way round. Every other test here reads one of them, or
plans inputs one of them has already designed on.
"""

import dataclasses
from collections import Counter
from itertools import pairwise

import pytest

from liulab_mbio.bench import (
    COLONY_ALLOWANCE,
    COLONY_FLANK,
    JUNCTION_OFFSET,
    SANGER_ALLOWANCE,
    SANGER_FLANK,
)
from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.cloning.goldengate import plan_assembly
from liulab_mbio.cloning.goldengate.oligos import DesignedOligo
from liulab_mbio.edits import rotate
from liulab_mbio.protocol import OVERVIEW_CHARS, read_protocol, render_html
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand, reverse_complement
from liulab_mbio.snapgene import read_dna

#: Where the fixture's own MCS feature sits, and the vector bases past it.
MCS = (395, 452)


def test_the_vector_junction_moves_one_base_off_an_all_gc_overhang(plan, puc19):
    # The vector spells GGCG where the MCS ends, and an all-GC junction truncates.
    assert puc19.sequence[MCS[1] : MCS[1] + 4] == "GGCG"
    assert plan.overhangs.overhangs == ("ATGA", "TGGC")
    assert plan.span == (MCS[0], MCS[1] - 1)
    assert [rejection.rule for rejection in plan.overhangs.choices[1].rejected] == ["uniform"]


def test_an_enzyme_with_a_site_in_the_parts_is_refused(puc19, gfp):
    with pytest.raises(ValueError, match="BsaI"):
        plan_assembly(puc19, gfp, enzyme="BsaI")


def test_the_codon_table_chooses_the_codon_a_proposed_domestication_moves_to(puc19, gfp):
    # Met-Glu-Asp-His-Leu-Leu: the Glu and Asp codons spell a BbsI site and the three after
    # them a PaqCI one, so no candidate is free and each carries what it would change.
    peptide = "ATGGAAGACCACCTGCTTTAA"
    coding = SequenceRecord(
        peptide,
        name="Peptide",
        features=(Feature("Peptide", "CDS", (Segment(0, len(peptide)),), strand=Strand.FORWARD),),
    )

    # E. coli spells aspartate GAT more often than glutamate GAG, and human the other way round.
    for table, moved in (("e-coli-k12", "GAC -> GAT"), ("human", "GAA -> GAG")):
        with pytest.raises(ValueError, match=moved):
            plan_assembly(puc19, gfp, coding, codon_table=table)


def test_the_insertion_site_is_read_from_a_feature_name_or_given_as_coordinates(puc19, gfp):
    named = plan_assembly(puc19, gfp, site="MCS")
    given = plan_assembly(puc19, gfp, site=MCS)
    assert named.product.sequence == given.product.sequence
    assert named.span == given.span


def test_a_vector_with_no_multiple_cloning_site_asks_where_to_put_the_insert(gfp):
    bare = SequenceRecord("ACGT" * 400, topology="circular", name="bare")
    with pytest.raises(ValueError, match="insertion site"):
        plan_assembly(bare, gfp)


def test_the_files_are_read_from_disk_when_a_path_is_given(puc19_file, gfp_file):
    made = plan_assembly(puc19_file, gfp_file)
    assert made.product.name == "pUC19-GFP"
    assert len(made.product) == 3347


def test_every_designed_primer_passes_evaluation(plan):
    # Two PCRs, three colony PCR primers and two sequencing primers.
    assert len(plan.reports) == 9
    left = [
        (report.primer.name, check.name)
        for report in plan.reports
        for check in report.checks
        if check.status not in (None, "pass")
    ]
    # What is left is what no primer within its placement could avoid: the backbone reverse
    # primer's Tm, which only moving the vector cut would clear, and the GFP reverse primer,
    # anchored at its junction and warning at every length there.
    assert left == [("pUC19 backbone reverse", "tm"), ("GFP reverse", "gc_percent")]
    assert plan.status != "fail"


def test_the_colony_pcr_sizes_are_the_ones_the_simulated_product_gives(plan, gfp):
    bands = {clone.name: clone.bands_bp for clone in plan.colony.clones}
    insert_bp = plan.phenotype.insert[1] - plan.phenotype.insert[0]
    removed = plan.span[1] - plan.span[0]
    start, end = plan.assembly.junction_positions
    forward, reverse, junction = (primer.binding_sites[0] for primer in plan.colony.primers)
    # Where the three primers landed inside their placements, which is what the bands count off.
    ahead, behind, into = start - forward.start, reverse.end - end, junction.end - start
    assert insert_bp == len(gfp)
    # The junction primer reaches the near vector primer and the two vector primers span the
    # insert. No overhang is another's reverse complement, so the insert cannot turn round and
    # there is no lane for it.
    assert bands == {
        "Correct clone": (ahead + into, ahead + insert_bp + behind),
        "Empty vector": (ahead + removed + behind,),
    }


def test_the_four_files_land_where_they_are_named_and_hold_what_the_plan_holds(plan, tmp_path):
    outputs = plan.write(tmp_path / "run")
    paths = (outputs.product, outputs.primers, outputs.protocol_data, outputs.protocol)
    assert [(path.parent, path.name) for path in paths] == [
        (tmp_path / "run", "product.dna"),
        (tmp_path / "run", "primers.tsv"),
        (tmp_path / "run", "protocol.json"),
        (tmp_path / "run", "protocol.html"),
    ]
    assert all(path.stat().st_size > 0 for path in paths)
    assert read_dna(outputs.product) == plan.product
    assert read_protocol(outputs.protocol_data) == plan.protocol()
    page = outputs.protocol.read_text(encoding="utf-8")
    assert page == render_html(read_protocol(outputs.protocol_data))


def test_the_same_inputs_write_the_same_bytes(plan, puc19, gfp, tmp_path):
    first = plan.write(tmp_path / "one")
    second = plan_assembly(puc19, gfp).write(tmp_path / "two")
    for one, other in (
        (first.product, second.product),
        (first.primers, second.primers),
        (first.protocol_data, second.protocol_data),
        (first.protocol, second.protocol),
    ):
        assert one.read_bytes() == other.read_bytes()


def test_the_primer_sheet_and_the_oligo_table_are_the_same_sheet(plan):
    rows = primer_sheet(plan.reports).splitlines()
    assert rows[0].split("\t") == ["name", "sequence", "length", "tm_c"]
    assert len(rows) == 1 + len(plan.reports)
    oligos = plan.protocol().oligos
    for row, report, oligo in zip(rows[1:], plan.reports, oligos, strict=True):
        name, sequence, length, tm = row.split("\t")
        assert (name, sequence) == (report.primer.name, report.primer.sequence)
        assert int(length) == len(sequence)
        assert float(tm) == pytest.approx(report["tm"].value, abs=0.05)
        assert (oligo.name, oligo.sequence) == (name, sequence)
        assert len(oligo.sequence) == int(length)
        assert oligo.tm_c == pytest.approx(float(tm), abs=0.05)


def test_every_oligo_row_carries_its_verdict_and_a_warned_one_says_why(plan):
    rows = {oligo.name: oligo for oligo in plan.protocol().oligos}
    assert [row.status for row in rows.values()] == [report.status for report in plan.reports]
    warned = rows["pUC19 backbone reverse"]
    assert warned.status == "warn"
    # The check, the value and the band it missed, in a few words each. The Tm is of the part
    # that anneals, so the tail is none of the 36 bases it was read from.
    assert (len(warned.sequence), [(check.name, check.detail) for check in warned.checks]) == (
        36,
        [("Tm", "64.1 °C (band 60-64)")],
    )
    assert (rows["GFP forward"].status, rows["GFP forward"].checks) == ("pass", ())


def test_the_oligo_summary_names_the_kinds_and_counts_the_rows_it_names(plan):
    detail = next(check for check in plan.checks if check.name == "primers").detail
    oligos = plan.protocol().oligos
    warned = sum(1 for oligo in oligos if oligo.status == "warn")
    assert detail.startswith(f"{len(oligos)} designed, {warned} with a warning, 0 failing")
    fired = Counter(check.name for oligo in oligos for check in oligo.checks)
    assert fired
    for label, rows in fired.items():
        assert f"{label} on {rows}" in detail
    # A check no sourced threshold judges is named, not counted as a pass.
    assert detail.endswith("not judged: full-primer Tm, 3' end stability")


def test_every_oligo_names_the_step_that_uses_it_and_is_no_material(plan):
    protocol = plan.protocol()
    titles = {step.title for step in protocol.steps}
    listed = " ".join(material.name + material.note for material in protocol.materials)
    for oligo in protocol.oligos:
        assert oligo.purpose in titles
        assert oligo.stock
        assert oligo.name not in listed
        assert oligo.sequence not in listed


def test_a_material_carries_a_catalogue_number_only_where_the_package_knows_one(plan):
    materials = {material.name: material for material in plan.protocol().materials}
    enzyme = materials[plan.enzyme.commercial_name]
    assert (enzyme.supplier, enzyme.catalog) == (plan.enzyme.supplier, plan.enzyme.catalog_number)
    assert materials["NEBridge Ligase Master Mix"].catalog == "M1100"
    # Nothing is invented for a reagent no product name names.
    assert materials["Agarose and 1X TAE or TBE"].catalog == ""
    assert materials["PCR and gel cleanup spin columns"].supplier == ""


def test_the_equipment_is_named_apart_from_the_reagents(plan):
    protocol = plan.protocol()
    assert "Thermocycler with a heated lid" in protocol.equipment
    assert not any(item in {m.name for m in protocol.materials} for item in protocol.equipment)


def test_every_step_of_the_protocol_says_what_a_good_result_looks_like(plan):
    steps = plan.protocol().steps
    assert len(steps) >= 10
    for step in steps:
        assert step.expected, step.title


def test_the_protocol_states_the_colony_colour_and_the_host_it_needs(plan):
    text = _sentences(plan.protocol())
    assert "white" in text
    assert "blue" in text
    assert "lacZ" in text
    assert plan.host in text
    assert "alpha-complementing" in text


def test_the_protocol_states_the_orientation_and_that_no_protein_is_expected(plan):
    prose = " ".join(plan.protocol().highlights)
    assert "opposite strand" in prose
    assert "lac promoter" in prose
    assert "no ribosome binding site" in prose
    assert "not expected to make" in prose
    assert not plan.phenotype.expressed


def test_the_header_splits_short_facts_from_sentences_and_verdicts(plan):
    protocol = plan.protocol()
    assert list(protocol.overview) == [
        "Vector",
        "Insert",
        "Enzyme",
        "Fragments",
        "Overhangs",
        "Fidelity",
        "Product",
        "Selection",
    ]
    # A card is a fact of a few words; a sentence is prose and never a card.
    assert all(len(value) <= OVERVIEW_CHARS for value in protocol.overview.values())
    assert not any(value.endswith(".") for value in protocol.overview.values())
    assert all(sentence.endswith(".") for sentence in protocol.highlights)
    assert [(one.name, one.status) for one in protocol.checks] == [
        (one.name, one.status) for one in plan.checks
    ]


def test_the_protocol_carries_the_numbers_the_package_computed(plan):
    protocol = plan.protocol()
    text = _sentences(protocol)
    for part in plan.parts:
        assert f"{part.length} bp" in text
    assert plan.enzyme.name in _sentences(protocol)
    gels = [gel for step in protocol.steps for gel in step.gels]
    lanes = {lane.label: lane.bands_bp for gel in gels for lane in gel.lanes}
    for clone in plan.colony.clones:
        assert lanes[clone.name] == clone.bands_bp
    programs = [program for step in protocol.steps for program in step.programs]
    assert any(program.title == "Golden Gate assembly" for program in programs)
    assert any(
        incubation.temperature_c == 60.0
        for program in programs
        for stage in program.stages
        for incubation in stage.incubations
    )


def test_the_protocol_cites_the_data_its_fidelity_came_from(plan):
    citations = " ".join(reference.text for reference in plan.protocol().references)
    assert "Pryor" in citations
    assert "NEBridge" in citations


def test_one_insert_plans_exactly_what_it_did_before(plan):
    # The two-fragment pin: taking any number of inserts changes none of these.
    assert plan.enzyme.name == "BbsI"
    assert (len(plan.product), plan.assembly.junction_positions) == (3347, (395, 1112))
    assert {clone.name: clone.bands_bp for clone in plan.colony.clones} == {
        "Correct clone": (160, 838),
        "Empty vector": (177,),
    }
    # Each part primer stays anchored where its junction put it; only validation primers move.
    assert [
        (site.start, site.end)
        for oligo in plan.designed_oligos
        if oligo.role == "amplification"
        for site in oligo.report.primer.binding_sites
    ] == [(455, 479), (377, 395), (4, 29), (691, 717)]
    assert len(plan.parts) == 2
    assert len(plan.reports) == 9


def test_a_site_at_the_vectors_base_zero_is_read_across_the_insert_and_not_the_backbone(
    plan, puc19, gfp
):
    # Turned so its multiple cloning site begins at base 0, the insert lands at the product's
    # end and the junction set runs across the origin instead of starting at base zero.
    turned = plan_assembly(rotate(puc19, MCS[0]), gfp)
    length = len(turned.product)
    assert turned.assembly.junction_positions == (length - len(gfp), length)
    # Every band, read length and phenotype is what the plan on the unturned vector gives.
    assert {clone.name: clone.bands_bp for clone in turned.colony.clones} == {
        clone.name: clone.bands_bp for clone in plan.colony.clones
    }
    assert [(read.distance_bp, read.read_bp) for read in turned.reads] == [
        (read.distance_bp, read.read_bp) for read in plan.reads
    ]
    assert turned.phenotype.insert == (length - len(gfp), length)
    assert turned.phenotype.coding is not None
    assert turned.phenotype.coding.name == plan.phenotype.coding.name


def test_a_vector_with_nothing_to_put_in_it_is_refused(puc19):
    with pytest.raises(ValueError, match="insert"):
        plan_assembly(puc19)


def test_every_insert_reaches_the_product_in_the_order_it_was_given(four, puc19, gfp, linker, tag):
    going = (gfp.sequence, reverse_complement(linker.sequence), tag.sequence)
    assert [part.name for part in four.parts] == ["pUC19 backbone", "GFP", "Linker", "Tag"]
    assert four.span == (395, 453)
    assert four.product.sequence == puc19.sequence[:395] + "".join(going) + puc19.sequence[453:]
    assert len(four.product) == 3645
    # The linker was given the other way round, so only its other strand is in the product, and
    # its feature turned over with it.
    assert linker.sequence not in four.product.sequence
    turned = next(one for one in four.product.features if one.name == "Linker")
    assert turned.strand == Strand.REVERSE
    assert len(four.assembly.junctions) == len(four.overhangs.overhangs) == 4
    assert four.assembly.status == "pass"


def test_the_plan_names_its_linearised_vector_and_each_insert_in_insert_order(four):
    assert four.linearised_vector.name == "pUC19 backbone"
    assert four.linearised_vector.template is four.vector
    assert [part.name for part in four.insert_parts] == ["GFP", "Linker", "Tag"]
    assert four.parts == (four.linearised_vector, *four.insert_parts)


def test_every_oligo_says_what_it_is_for(plan):
    assert [(oligo.report.primer.name, oligo.role) for oligo in plan.designed_oligos] == [
        ("pUC19 backbone forward", "amplification"),
        ("pUC19 backbone reverse", "amplification"),
        ("GFP forward", "amplification"),
        ("GFP reverse", "amplification"),
        ("Colony PCR forward", "colony PCR"),
        ("Colony PCR reverse", "colony PCR"),
        ("Junction reverse", "colony PCR"),
        ("Sequencing forward", "sequencing"),
        ("Sequencing reverse", "sequencing"),
    ]
    for oligo in plan.designed_oligos:
        if oligo.role == "amplification":
            assert oligo.part is not None
            assert oligo.report in (oligo.part.report.forward, oligo.part.report.reverse)
        else:
            assert oligo.part is None


def test_an_oligo_keeps_its_purpose_wherever_it_is_listed(plan):
    turned = dataclasses.replace(plan, designed_oligos=plan.designed_oligos[::-1])
    written = {oligo.name: oligo.purpose for oligo in plan.protocol().oligos}
    assert {oligo.name: oligo.purpose for oligo in turned.protocol().oligos} == written
    assert written["GFP reverse"] == "Amplify GFP"
    assert written["Junction reverse"] == "Screen colonies by PCR"
    assert written["Sequencing reverse"] == "Confirm the clone by sequencing"


def test_an_oligo_amplifies_a_part_exactly_when_that_is_its_role(plan):
    colony, amplifying = plan.designed_oligos[4], plan.designed_oligos[0]
    with pytest.raises(ValueError, match="part"):
        DesignedOligo(colony.report, "colony PCR", amplifying.part)
    with pytest.raises(ValueError, match="part"):
        DesignedOligo(amplifying.report, "amplification")


def test_a_plan_whose_parts_are_not_vector_first_writes_the_same_protocol(four):
    turned = dataclasses.replace(
        four, assembly=dataclasses.replace(four.assembly, parts=four.assembly.parts[::-1])
    )
    # The checks judge the parts in the assembly's own order, so only their order may differ.
    written, again = four.protocol(), turned.protocol()
    assert sorted(again.checks, key=repr) == sorted(written.checks, key=repr)
    assert dataclasses.replace(again, checks=()) == dataclasses.replace(written, checks=())


def test_the_fidelity_of_the_set_falls_as_junctions_are_added(plan, four):
    assert four.overhangs.fidelity.measured
    assert four.overhangs.fidelity.enzyme == plan.overhangs.fidelity.enzyme
    assert four.overhangs.fidelity.value < plan.overhangs.fidelity.value


def test_the_reaction_takes_one_amount_for_every_part(four):
    assert [amount.name for amount in four.amounts] == [part.name for part in four.parts]
    assert len(four.amounts) == 4


def test_the_cycling_tier_follows_the_fragment_count(plan, four):
    two, many = _assembly_program(plan).stages[0], _assembly_program(four).stages[0]
    # NEB holds two fragments at 37 °C and cycles three or more of them.
    assert (two.cycles, two.incubations[0].seconds) == (1, 900)
    assert (many.cycles, many.incubations[0].seconds) == (30, 60)


def test_every_validation_primer_lies_where_its_placement_allows(plan, four):
    for made in (plan, four):
        junctions = made.assembly.junction_positions
        start, end = junctions[0], junctions[-1]
        forward, reverse, *inserts = (primer.binding_sites[0] for primer in made.colony.primers)
        ahead, behind = start - forward.start, reverse.end - end
        assert abs(ahead - COLONY_FLANK) <= COLONY_ALLOWANCE
        assert abs(behind - COLONY_FLANK) <= COLONY_ALLOWANCE
        for site, (first, last) in zip(inserts, pairwise(junctions), strict=True):
            assert abs((site.end - first) - JUNCTION_OFFSET) <= COLONY_ALLOWANCE
            assert first <= site.start < site.end <= last
        for read in made.reads:
            assert SANGER_FLANK <= read.distance_bp <= SANGER_FLANK + SANGER_ALLOWANCE


def test_the_colony_pcr_reads_every_junction(four):
    check = four.colony
    assert len(check.primers) == 5
    assert [clone.name for clone in check.clones] == ["Correct clone", "Empty vector"]
    correct = next(clone.bands_bp for clone in check.clones if clone.name == "Correct clone")
    assert len(correct) == 4


def test_the_protocol_names_every_part_and_every_junction(four):
    protocol = four.protocol()
    prose = " ".join(protocol.highlights)
    for part in four.parts:
        assert part.name in prose
    results = " ".join(line for step in protocol.steps for line in step.expected)
    for junction in four.assembly.junctions:
        assert f"{junction.start}" in results
        assert junction.overhang in results
    assert f"{four.overhangs.fidelity.value:.0%}" in protocol.overview["Fidelity"]


def test_a_step_is_per_experiment_and_not_per_pair_of_fragments(plan, four):
    # One more PCR step for each part, and the rest of the protocol the same length.
    assert len(four.protocol().steps) - len(plan.protocol().steps) == 2


def _assembly_program(plan):
    """The Golden Gate program the protocol prints."""
    return next(
        program
        for step in plan.protocol().steps
        for program in step.programs
        if program.title == "Golden Gate assembly"
    )


def _sentences(protocol) -> str:
    """Every sentence the protocol says, for a test to read."""
    parts = [
        protocol.title,
        protocol.summary,
        *protocol.overview.values(),
        *protocol.highlights,
    ]
    for material in protocol.materials:
        parts += [material.name, material.note]
    for oligo in protocol.oligos:
        parts += [oligo.name, oligo.purpose]
    for step in protocol.steps:
        parts += [step.title, *step.instructions, *step.cautions, *step.notes, *step.expected]
        for entry in step.troubleshooting:
            parts += [entry.problem, entry.solution]
    return " ".join(parts)
