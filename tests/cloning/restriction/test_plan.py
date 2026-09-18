"""The whole pipeline, run on the fixtures with no agent.

`tests/data/pUC19.dna` is the vector. The insert reaches it by both routes: cut out of the same
record carrying `tests/data/GFP.dna` between its own EcoRI and BamHI sites, which is what a lab
member subcloning out of one plasmid into another holds, and amplified from the GFP record
itself, which carries neither site. Nothing about the fixtures is hard-coded in the package; the
numbers below are read off the records and pinned here.
"""

import dataclasses

import pytest

from liulab_mbio.cloning.restriction import Plan, plan_restriction
from liulab_mbio.cloning.restriction.amplify import amplified
from liulab_mbio.cloning.restriction.bench import (
    BLUNT_SECONDS,
    COHESIVE_SECONDS,
    HIGH_LIGASE_UNITS_UL,
    LIGASE_UNITS_UL,
    PHOSPHATASE_KILL_SECONDS,
    PHOSPHATASE_SECONDS,
    ROOM_CELSIUS,
    phosphatase_units,
    shared_buffer,
)
from liulab_mbio.cloning.restriction.design import refusal
from liulab_mbio.cloning.restriction.digest import opened, resolve
from liulab_mbio.edits import carried, flipped, replace
from liulab_mbio.edits import insert as added
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.primers.evaluation import evaluate_primer
from liulab_mbio.protocol import OVERVIEW_CHARS, read_protocol, render_html
from liulab_mbio.sequence import Primer, SequenceRecord, reverse_complement
from liulab_mbio.sites import SPACER_LENGTH, find_sites
from liulab_mbio.snapgene import read_dna

#: Where pUC19's own EcoRI and BamHI sites begin, and the bases the two cuts leave between them.
ECORI, BAMHI, STUFFER = 395, 416, 21

#: Bases that spell no site of any enzyme the refusal tests name, for padding a made-up plasmid.
FILLER = "ACGT" * 5


def plasmid(sequence: str, name: str) -> SequenceRecord:
    """A circular record written in code, for a refusal that needs no real plasmid."""
    return SequenceRecord(sequence, topology="circular", name=name)


def carrying(vector: SequenceRecord, insert: SequenceRecord, name: str) -> SequenceRecord:
    """The vector with `insert` put between its own EcoRI and BamHI sites, both sites kept.

    What the insert annotates comes with it, as it would in a plasmid someone built.
    """
    first, second = sorted(find_sites(vector, ["EcoRI", "BamHI"]), key=lambda site: site.start)
    built, _ = replace(vector, first.end, second.start, insert.sequence)
    features, primers = carried(insert, 0, len(insert), offset=first.end)
    return dataclasses.replace(
        built,
        name=name,
        features=(*built.features, *features),
        primers=(*built.primers, *primers),
    )


def padded(record: SequenceRecord, bases: str) -> SequenceRecord:
    """The record with `bases` added at each end, what it annotates carried along."""
    one, _ = added(record, len(record), bases)
    one, _ = added(one, 0, bases)
    return one


@pytest.fixture(scope="module")
def source(puc19: SequenceRecord, gfp: SequenceRecord) -> SequenceRecord:
    """The plasmid GFP is cut out of."""
    return carrying(puc19, gfp, "pTrc-GFP")


@pytest.fixture(scope="module")
def blunt(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP cut out and put back with SmaI alone: one enzyme, blunt ends, both ways round.

    pUC19 reads one SmaI site, so the whole plasmid is the backbone; the plasmid GFP comes out
    of carries a SmaI site either side of it, which is the two cuts that release it.
    """
    return plan_restriction(
        puc19, carrying(puc19, padded(gfp, "CCCGGG"), "pSmaI-GFP"), enzymes=["SmaI"]
    )


@pytest.fixture(scope="module")
def made(puc19: SequenceRecord, source: SequenceRecord) -> Plan:
    """GFP cut out of its own plasmid and ligated into pUC19, every option left at its default."""
    return plan_restriction(puc19, source, enzymes=["EcoRI", "BamHI"])


def test_the_digest_keeps_the_backbone_and_the_insert_and_names_what_it_drops(made, gfp, puc19):
    assert [(one.name, one.length) for one in made.vector_pieces] == [
        ("pUC19 backbone", len(puc19) - STUFFER),
        ("pUC19 offcut", STUFFER),
    ]
    assert made.insert.name == "pTrc-GFP insert"
    # The insert is GFP plus the bases the two cuts leave either side of it.
    assert gfp.sequence in made.insert.bases
    assert (made.insert.left_enzyme.name, made.insert.right_enzyme.name) == ("EcoRI", "BamHI")
    assert made.span == (ECORI + 1, BAMHI + 1)


def test_the_product_carries_the_insert_once_with_its_features_at_their_new_coordinates(made, gfp):
    product = made.product
    assert len(product) == len(made.backbone.bases) + len(made.insert.bases)
    assert product.sequence.count(made.insert.bases) == 1
    coding = next(one for one in product.features if one.name == gfp.name)
    assert product.extract(coding.segments[0]) == gfp.sequence
    assert made.ligation.status == "pass"


def test_every_designed_oligo_is_annotated_where_it_binds_on_the_product(made):
    placed = {primer.name: primer for primer in made.product.primers}
    for report in made.reports:
        primer = report.primer
        assert primer.name in placed, primer.name
        site = placed[primer.name].binding_sites[0]
        read = made.product.extract(site)
        assert read == (primer.sequence if site.strand > 0 else reverse_complement(primer.sequence))
    # What the vector annotates came with it: the fixture draws its M13 sites as features.
    assert {"M13 fwd", "M13 rev"} <= {one.name for one in made.product.features}


def test_each_junction_is_marked_and_spells_the_site_its_two_ends_came_from(made):
    assert [(one.start, one.spells, one.enzyme) for one in made.junctions] == [
        (396, "GAATTC", "EcoRI"),
        (1119, "GGATCC", "BamHI"),
    ]
    drawn = {one.name: one for one in made.product.features if one.name.endswith("junction")}
    assert set(drawn) == {"GAATTC junction", "GGATCC junction"}
    for one in made.junctions:
        feature = drawn[f"{one.spells} junction"]
        assert made.product.extract(feature.segments[0]) == one.spells
        assert one.span is not None
        assert made.product.extract(one.span) == one.overhang


def test_the_insert_goes_in_the_same_way_round_whichever_strand_its_own_plasmid_wrote_it_on(
    puc19, source
):
    turned = plan_restriction(puc19, flipped(source), enzymes=["EcoRI", "BamHI"])
    made = plan_restriction(puc19, source, enzymes=["EcoRI", "BamHI"])
    assert turned.product.sequence == made.product.sequence


def test_a_pair_whose_ends_do_not_anneal_is_refused_naming_the_two_ends():
    # BsaI cuts outside its own site, so the bases it leaves are whatever the plasmid holds
    # there: AAAA in one of these and GGGG in the other, which cannot pair.
    vector = plasmid("GGTCTCAAAAA" + FILLER * 8 + "GAATTC" + FILLER * 2, "pBsaI-A")
    holder = plasmid("GGTCTCAGGGG" + FILLER * 2 + "GAATTC" + FILLER * 8, "pBsaI-G")
    with pytest.raises(ValueError, match=r"do not anneal.*(AAAA|GGGG)"):
        plan_restriction(vector, holder, enzymes=["BsaI", "EcoRI"])


def test_an_enzyme_with_a_site_inside_the_insert_is_refused_naming_the_site(puc19, gfp, source):
    inside = carrying(puc19, gfp, "pTrc-GFP")
    sites = find_sites(inside, "BsaI")
    assert len(sites) == 2
    with pytest.raises(
        ValueError, match=rf"2 BsaI sites \(at {sites[0].start} and {sites[1].start}\)"
    ):
        plan_restriction(puc19, inside, enzymes=["BsaI", "EcoRI"])


def test_with_no_enzymes_named_the_plan_chooses_the_pair_and_the_protocol_names_it(puc19, source):
    picked = plan_restriction(puc19, source)
    assert [one.name for one in picked.enzymes] == ["BamHI", "EcoRI"]
    assert (
        picked.product.sequence
        == plan_restriction(puc19, source, enzymes=["EcoRI", "BamHI"]).product.sequence
    )
    # Every pair refused is on the plan with the rule that refused it.
    assert {one.rule for one in picked.refusals} == {"vector site", "insert site", "backbone"}
    assert all(len(one.enzymes) == 2 for one in picked.refusals)
    said = " ".join(picked.protocol().highlights)
    assert "the pair was chosen: BamHI and EcoRI" in said
    assert "on the vector site rule" in said
    # A site a synonymous codon change could take away is reported as that, not as a flat no.
    domesticable = [one for one in picked.refusals if one.domesticable]
    assert domesticable
    assert all(one.rule == "insert site" for one in domesticable)


def test_naming_an_unusable_pair_refuses_in_the_choosers_own_words(puc19, source):
    with pytest.raises(ValueError, match="XhoI cuts pUC19 0 time") as named:
        plan_restriction(puc19, source, enzymes=["XhoI", "EcoRI"])
    chosen = refusal(puc19, source, resolve(["XhoI", "EcoRI"]))
    assert chosen is not None
    assert chosen.detail == str(named.value)


def test_one_enzyme_opens_the_vector_and_cuts_the_insert_out_on_its_own(blunt, puc19, gfp):
    assert [one.name for one in blunt.enzymes] == ["SmaI"]
    # One site opens the vector, so the whole plasmid is the backbone and the gel drops nothing.
    assert [(one.name, one.length) for one in blunt.vector_pieces] == [
        ("pUC19 backbone", len(puc19))
    ]
    assert gfp.sequence in blunt.insert.bases
    assert [(one.spells, one.enzyme) for one in blunt.junctions] == [
        ("CCCGGG", "SmaI"),
        ("CCCGGG", "SmaI"),
    ]
    assert blunt.status == "pass"


def test_a_vector_that_closes_on_itself_is_dephosphorylated_rather_than_refused(blunt):
    assert blunt.dephosphorylates
    verdict = next(one for one in blunt.checks if one.name == "self-ligation")
    assert verdict.status == "pass"
    assert "two blunt ends" in verdict.detail
    assert "anneal to each other" in verdict.detail
    protocol = blunt.protocol()
    assert "Shrimp Alkaline Phosphatase (rSAP)" in [one.name for one in protocol.materials]
    step = next(one for one in protocol.steps if one.title == "Dephosphorylate the cut pUC19")
    assert f"{phosphatase_units(blunt.digests[0]):g} units of rSAP" in step.instructions[0]
    assert [timer.seconds for timer in step.timers] == [
        PHOSPHATASE_SECONDS,
        PHOSPHATASE_KILL_SECONDS,
    ]
    assert any("rSAP" in one.text for one in protocol.references)
    # The page says the risk plainly, not only in the badge.
    assert "closes on itself with no insert" in " ".join(protocol.highlights)


def test_a_blunt_ligation_is_held_longer_and_the_page_says_what_that_costs(blunt):
    assert blunt.ligation.blunt
    step = next(one for one in blunt.protocol().steps if one.title.startswith("Ligate"))
    hold = step.programs[0].stages[0].incubations[0]
    assert (hold.temperature_c, hold.seconds) == (ROOM_CELSIUS, BLUNT_SECONDS)
    said = " ".join(step.notes)
    # The cost against a cohesive ligation is the incubation, which is all NEB states: the
    # longer hold, or the same short one with five times the ligase.
    assert (
        f"{BLUNT_SECONDS // 60} minutes at {ROOM_CELSIUS:g} °C where cohesive ends take "
        f"{COHESIVE_SECONDS // 60}" in said
    )
    assert f"{HIGH_LIGASE_UNITS_UL:g} U/µL ligase in place of the {LIGASE_UNITS_UL:g} U/µL" in said


def test_the_colony_pcr_tells_a_reversed_insert_from_a_correct_one_where_one_can_exist(blunt):
    lanes = {clone.name: clone.bands_bp for clone in blunt.colony.clones}
    assert set(lanes) == {"Correct clone", "Empty vector", "Reversed insert"}
    # A blunt end's two strands stop at the same base, so turning the product's own top strand
    # in place gives these bands too.
    assert lanes["Reversed insert"] == (213, 907)
    assert lanes["Reversed insert"] != lanes["Correct clone"]
    assert blunt.colony.tells_orientation
    # Two flanking vector primers cannot tell them apart; one reading out of the insert can.
    assert "Junction reverse" in [one.name for one in blunt.colony.primers]
    step = next(one for one in blunt.protocol().steps if one.title == "Screen colonies by PCR")
    said = " ".join(step.expected)
    for name, bands in lanes.items():
        assert f"{name}: {', '.join(f'{bp} bp' for bp in bands)}." in said


def test_a_reversed_insert_is_the_plasmid_its_own_ligation_makes_turned_over(puc19, gfp):
    cohesive = plan_restriction(puc19, gfp, enzymes=["EcoRI"])
    lanes = {clone.name: clone.bands_bp for clone in cohesive.colony.clones}
    # GFP turned over and ligated back into the same EcoRI ends, both sites put back. Turning
    # the product's own top strand over in place instead moves the four-base overhang to the
    # far side of each junction, restores neither site, and reads 231 bp here.
    assert lanes["Reversed insert"] == (227, 905)


def test_no_reversed_lane_is_invented_where_the_insert_cannot_go_in_backwards(made):
    assert not made.colony.reversed_clones
    assert [clone.name for clone in made.colony.clones] == ["Correct clone", "Empty vector"]
    step = next(one for one in made.protocol().steps if one.title == "Screen colonies by PCR")
    said = " ".join((*step.expected, *step.notes))
    assert "cannot go in the other way round" in said
    assert "Reversed insert" not in said


def test_the_four_outputs_land_in_the_directory_the_caller_names(made, tmp_path):
    outputs = made.write(tmp_path / "run")
    assert [(path.parent, path.name) for path in outputs.paths] == [
        (tmp_path / "run", "product.dna"),
        (tmp_path / "run", "primers.tsv"),
        (tmp_path / "run", "protocol.json"),
        (tmp_path / "run", "protocol.html"),
    ]
    assert all(path.stat().st_size > 0 for path in outputs.paths)
    assert read_dna(outputs.product) == made.product
    assert read_protocol(outputs.protocol_data) == made.protocol()
    assert outputs.protocol.read_text(encoding="utf-8") == render_html(
        read_protocol(outputs.protocol_data)
    )


def test_the_same_inputs_write_the_same_bytes(made, puc19, source, tmp_path):
    first = made.write(tmp_path / "one")
    again = plan_restriction(puc19, source, enzymes=["EcoRI", "BamHI"]).write(tmp_path / "two")
    for one, other in zip(first.paths, again.paths, strict=True):
        assert one.read_bytes() == other.read_bytes()


def test_the_protocol_runs_the_bench_from_the_digests_to_the_sequencing(made):
    protocol = made.protocol()
    assert [step.title for step in protocol.steps] == [
        "Digest pUC19 with EcoRI and BamHI",
        "Digest pTrc-GFP with EcoRI and BamHI",
        "Separate the digests on a gel and recover the two fragments",
        "Measure every concentration",
        "Ligate the insert into the backbone",
        "Transform and plate",
        "Screen colonies by PCR",
        "Check a miniprep by digesting it with EcoRI and BamHI",
        "Confirm the clone by sequencing",
    ]
    for step in protocol.steps:
        assert step.expected, step.title
    tables = {table.title for step in protocol.steps for table in step.tables}
    assert tables == {
        "pUC19 digest",
        "pTrc-GFP digest",
        "Ligation, T4 DNA Ligase (NEB #M0202)",
        "Colony PCR",
        "Diagnostic digest",
    }
    lanes = {
        lane.label: lane.bands_bp
        for step in protocol.steps
        for gel in step.gels
        for lane in gel.lanes
    }
    assert lanes["pUC19"] == (made.backbone.length, STUFFER)
    assert lanes["pTrc-GFP"] == (made.insert.length, made.source_pieces[1].length)


def test_the_protocol_carries_the_numbers_the_note_states_and_none_it_does_not(made):
    protocol = made.protocol()
    steps = {step.title: step for step in protocol.steps}
    digest = steps["Digest pUC19 with EcoRI and BamHI"]
    table = {component.name: component for component in digest.tables[0].components}
    # NEB's typical digest: 10 units of each enzyme, 5 µL of 10X buffer, 1 µg of DNA, 50 µL.
    assert table["EcoRI-HF (R3101)"].final == "10 units"
    assert table["10X rCutSmart Buffer"].volume_ul == 5.0
    assert table["pUC19"].final.endswith("(1000 ng)")
    assert sum(one.volume_ul for one in digest.tables[0].components) == 50.0
    ligation = steps["Ligate the insert into the backbone"]
    ligate, kill = (one for stage in ligation.programs[0].stages for one in stage.incubations)
    assert (ligate.temperature_c, ligate.seconds) == (25.0, 600)
    assert (kill.temperature_c, kill.seconds) == (65.0, 600)
    # The 1:3 in NEB's table is an example for a 4 kb vector and a 1 kb insert, so the ratio is
    # taken in picomoles and the masses fall out of the fragments' own lengths.
    backbone, insert = made.amounts
    assert (backbone.pmol, insert.pmol) == (0.02, 0.06)
    assert insert.nanograms < backbone.nanograms
    plate = " ".join(steps["Transform and plate"].expected)
    assert "No supplier states a colony count" in plate
    assert "under 1% of the colonies the uncut vector gave" in plate


def test_the_buffer_is_named_only_where_both_records_are_supplied_in_one(made):
    assert shared_buffer([get_enzyme("EcoRI"), get_enzyme("BamHI")]) == "rCutSmart Buffer"
    assert shared_buffer([get_enzyme("EcoRI"), get_enzyme("BsmBI")]) is None
    said = " ".join(note for step in made.protocol().steps for note in step.notes)
    assert "Both enzymes are supplied in rCutSmart Buffer" in said


def test_the_page_says_what_each_junction_now_spells(made):
    protocol = made.protocol()
    assert protocol.overview["Junctions"] == "GAATTC and GGATCC"
    assert all(len(value) <= OVERVIEW_CHARS for value in protocol.overview.values())
    prose = " ".join(protocol.highlights)
    assert "not scarless" in prose
    assert "GAATTC at 396 (EcoRI)" in prose
    assert "GGATCC at 1119 (BamHI)" in prose


def test_the_plan_status_is_the_worst_of_its_checks(made):
    assert [check.name for check in made.checks] == [
        "buffer",
        "digest temperature",
        "digest clean-up",
        "methylation",
        "self-ligation",
        "pUC19 backbone",
        "pTrc-GFP insert",
        "junctions",
        "reading frame",
        "ligation ratio",
        "diagnostic digest",
        "primers",
    ]
    assert made.status == "pass"
    # Every check reaches the page, the ones carrying no verdict with it.
    assert [one.status for one in made.protocol().checks] == [one.status for one in made.checks]


def test_the_diagnostic_digest_tells_a_correct_clone_from_the_vector_it_went_into(made, puc19):
    assert made.diagnostic.clone == (made.backbone.length, made.insert.length)
    assert made.diagnostic.empty == (len(puc19) - STUFFER, STUFFER)
    step = next(one for one in made.protocol().steps if "miniprep" in one.title)
    lanes = {lane.label: lane.bands_bp for gel in step.gels for lane in gel.lanes}
    assert lanes == {
        made.product.name: made.diagnostic.clone,
        "pUC19": made.diagnostic.empty,
    }


def test_the_protocol_cites_the_note_the_bench_numbers_came_from(made):
    citations = " ".join(reference.text for reference in made.protocol().references)
    assert "Optimizing Restriction Endonuclease Reactions" in citations
    assert "T4 DNA Ligase" in citations
    assert "Monarch Spin DNA Gel Extraction Kit" in citations
    assert "Troubleshooting Guide for Cloning" in citations


@pytest.fixture(scope="module")
def tailed(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP amplified with a site on each tail, no plasmid holding it between two of them."""
    return plan_restriction(puc19, gfp, enzymes=["EcoRI", "BamHI"])


def test_an_insert_no_plasmid_holds_is_amplified_with_a_tailed_primer_at_each_end(tailed, gfp):
    amplicon = tailed.amplicon
    assert amplicon is not None
    forward, reverse = amplicon.primers
    assert forward.sequence.startswith(amplicon.left_tail)
    assert reverse.sequence.startswith(reverse_complement(amplicon.right_tail))
    assert amplicon.record.sequence == amplicon.left_tail + gfp.sequence + amplicon.right_tail
    # Each tail is the spacer the note asks for and then the site, and the plan says which
    # enzyme each end is for.
    for _, tail, enzyme in amplicon.ends:
        assert enzyme.site in tail
        assert len(tail) - len(enzyme.site) == SPACER_LENGTH
    ends = (amplicon.left_enzyme.name, amplicon.right_enzyme.name)
    assert ends == ("EcoRI", "BamHI")
    assert (tailed.insert.left_enzyme.name, tailed.insert.right_enzyme.name) == ends
    # The two tails spell one site of each enzyme and no more, so the amplicon cuts cleanly.
    found = [site for site in find_sites(amplicon.record, ["EcoRI", "BamHI"]) if site.cuts]
    assert sorted(site.enzyme.name for site in found) == ["BamHI", "EcoRI"]


def test_both_routes_to_the_insert_reach_the_same_product(tailed, made):
    assert tailed.product.sequence == made.product.sequence


def test_a_fragment_already_carrying_the_sites_is_taken_as_it_is(tailed, puc19):
    given = plan_restriction(puc19, tailed.amplicon.record, enzymes=["EcoRI", "BamHI"])
    assert given.amplicon is None
    assert given.product.sequence == tailed.product.sequence


def test_the_tm_is_the_annealing_region_and_the_structures_are_the_whole_oligo(tailed):
    designed = [one for one in tailed.designed_oligos if one.role == "amplification"]
    assert [one.report.primer.name for one in designed] == ["GFP forward", "GFP reverse"]
    for report in (one.report for one in designed):
        primer = report.primer
        alone = evaluate_primer(Primer(primer.name, primer.sequence))
        assert report["tm"].value != alone["tm"].value
        assert report["tm_full"].value == alone["tm"].value
        assert report["hairpin"].value == alone["hairpin"].value
        assert report["self_dimer"].value == alone["self_dimer"].value


def test_the_protocol_gains_the_pcr_its_program_the_amplicon_gel_and_the_template_removal(tailed):
    protocol = tailed.protocol()
    pcr, gel, cleanup, *rest = protocol.steps
    assert [pcr.title, gel.title, cleanup.title] == [
        "Amplify GFP",
        "Check the PCRs on a gel",
        "Purify every amplicon",
    ]
    assert rest[0].title == "Digest pUC19 with EcoRI and BamHI"
    assert [stage.incubations[0].label for stage in pcr.programs[0].stages]
    assert [lane.bands_bp for lane in gel.gels[0].lanes] == [(tailed.amplicon.length,)]
    # A linear template neither transforms nor ligates, so the column is what takes it away.
    assert "take GFP away" in " ".join(cleanup.notes)
    assert [row.purpose for row in protocol.oligos][:2] == ["Amplify GFP"] * 2
    said = " ".join(protocol.highlights)
    assert "forward primer, 6 spacer bases and the EcoRI site" in said
    assert "reverse primer, 6 spacer bases and the BamHI site" in said


def test_a_plasmid_template_is_taken_away_with_dpni(puc19, gfp):
    made = plan_restriction(
        puc19,
        SequenceRecord(gfp.sequence, topology="circular", name="pGFP"),
        enzymes=["EcoRI", "BamHI"],
    )
    amplicon = made.amplicon
    assert amplicon is not None
    assert amplicon.dpni
    titles = [step.title for step in made.protocol().steps]
    assert "Digest the plasmid template with DpnI" in titles


def test_one_enzyme_puts_its_own_site_on_both_tails(puc19, gfp):
    # A tail spells its own enzyme's site by design, so that enzyme is not its own avoid.
    amplicon = amplified(gfp, into=opened(puc19, resolve(["SmaI"]))[0])
    assert amplicon.left_enzyme.name == amplicon.right_enzyme.name == "SmaI"
    assert len([site for site in find_sites(amplicon.record, "SmaI") if site.cuts]) == 2


def test_a_tail_spelling_a_second_site_against_the_insert_is_refused_naming_the_enzyme(puc19, gfp):
    # An EcoRI tail ends GAATTC, and an insert beginning TAGA completes TCTAGA across the join:
    # an XbaI site no spacer can be moved off, because the spacer sits the other side of it.
    edge = SequenceRecord("TAGA" + gfp.sequence, name="GFP")
    with pytest.raises(ValueError, match=r"XbaI cuts GFP amplicon 2 time"):
        plan_restriction(puc19, edge, enzymes=["EcoRI", "XbaI"])
