"""The whole pipeline, run on the fixtures with no agent.

`tests/data/pUC19.dna` is the vector. The insert reaches it by both routes: cut out of the
`source` record carrying `tests/data/GFP.dna` between its own EcoRI and BamHI sites, which is
what a lab member subcloning out of one plasmid into another holds, and amplified from the GFP
record itself, which carries neither site. What each function below the plan decides is tested
on that function; what is here is the plan whole, the files it writes, the page it renders, and
the route it takes for each of its critical edge cases.
"""

import pytest

from liulab_mbio.cloning.restriction import Plan, plan_restriction
from liulab_mbio.cloning.restriction.bench import (
    BLUNT_SECONDS,
    COHESIVE_SECONDS,
    HIGH_LIGASE_UNITS_UL,
    LIGASE_UNITS_UL,
    PHOSPHATASE_KILL_SECONDS,
    PHOSPHATASE_SECONDS,
    ROOM_CELSIUS,
    phosphatase_units,
)
from liulab_mbio.cloning.restriction.design import refusal
from liulab_mbio.cloning.restriction.digest import resolve
from liulab_mbio.edits import flipped
from liulab_mbio.protocol import OVERVIEW_CHARS, read_protocol, render_html
from liulab_mbio.sequence import SequenceRecord, reverse_complement
from liulab_mbio.snapgene import read_dna

from .records import BAMHI, ECORI, STUFFER, carrying, padded

# --------------------------------------------------------------------------------------
# The plans: the common path, and one for each route the method takes differently
# --------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def made(puc19: SequenceRecord, source: SequenceRecord) -> Plan:
    """GFP cut out of its own plasmid and ligated into pUC19, every option left at its default."""
    return plan_restriction(puc19, source, enzymes=["EcoRI", "BamHI"])


@pytest.fixture(scope="module")
def tailed(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP amplified with a site on each tail, no plasmid holding it between two of them."""
    return plan_restriction(puc19, gfp, enzymes=["EcoRI", "BamHI"])


@pytest.fixture(scope="module")
def blunt(puc19: SequenceRecord, gfp: SequenceRecord) -> Plan:
    """GFP cut out and put back with SmaI alone: one enzyme, blunt ends, both ways round.

    pUC19 reads one SmaI site, so the whole plasmid is the backbone; the plasmid GFP comes out
    of carries a SmaI site either side of it, which is the two cuts that release it.
    """
    return plan_restriction(
        puc19, carrying(puc19, padded(gfp, "CCCGGG"), "pSmaI-GFP"), enzymes=["SmaI"]
    )


# --------------------------------------------------------------------------------------
# The common path
# --------------------------------------------------------------------------------------


def test_the_product_replaces_the_cloning_site_and_carries_every_designed_oligo_where_it_binds(
    made,
):
    # The vector bases the insert replaced, which is what the phenotype is read against.
    assert made.span == (ECORI + 1, BAMHI + 1)
    placed = {primer.name: primer for primer in made.product.primers}
    for report in made.reports:
        primer = report.primer
        assert primer.name in placed, primer.name
        site = placed[primer.name].binding_sites[0]
        read = made.product.extract(site)
        assert read == (primer.sequence if site.strand > 0 else reverse_complement(primer.sequence))


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
    steps = {step.title: step for step in protocol.steps}
    assert list(steps) == [
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
    # No supplier states a colony count for this method, so the plate promises a ratio.
    plate = " ".join(steps["Transform and plate"].expected)
    assert "No supplier states a colony count" in plate
    assert "under 1% of the colonies the uncut vector gave" in plate


def test_the_page_says_what_each_junction_now_spells_and_which_buffer_both_enzymes_share(made):
    protocol = made.protocol()
    assert protocol.overview["Junctions"] == "GAATTC and GGATCC"
    assert all(len(value) <= OVERVIEW_CHARS for value in protocol.overview.values())
    prose = " ".join(protocol.highlights)
    assert "not scarless" in prose
    assert "GAATTC at 396 (EcoRI)" in prose
    assert "GGATCC at 1119 (BamHI)" in prose
    said = " ".join(note for step in protocol.steps for note in step.notes)
    assert "Both enzymes are supplied in rCutSmart Buffer" in said


def test_the_diagnostic_digest_tells_a_correct_clone_from_the_vector_it_went_into(made):
    step = next(one for one in made.protocol().steps if "miniprep" in one.title)
    lanes = {lane.label: lane.bands_bp for gel in step.gels for lane in gel.lanes}
    assert lanes == {
        made.product.name: made.diagnostic.clone,
        "pUC19": made.diagnostic.empty,
    }


def test_no_reversed_lane_is_invented_where_the_insert_cannot_go_in_backwards(made):
    assert not made.colony.reversed_clones
    assert [clone.name for clone in made.colony.clones] == ["Correct clone", "Empty vector"]
    step = next(one for one in made.protocol().steps if one.title == "Screen colonies by PCR")
    said = " ".join((*step.expected, *step.notes))
    assert "cannot go in the other way round" in said
    assert "Reversed insert" not in said


def test_the_protocol_cites_the_note_the_bench_numbers_came_from(made):
    citations = " ".join(reference.text for reference in made.protocol().references)
    assert "Optimizing Restriction Endonuclease Reactions" in citations
    assert "T4 DNA Ligase" in citations
    assert "Monarch Spin DNA Gel Extraction Kit" in citations
    assert "Troubleshooting Guide for Cloning" in citations


# --------------------------------------------------------------------------------------
# The insert written on the other strand, and the pair chosen rather than named
# --------------------------------------------------------------------------------------


def test_the_insert_goes_in_the_same_way_round_whichever_strand_its_own_plasmid_wrote_it_on(
    puc19, source, made
):
    turned = plan_restriction(puc19, flipped(source), enzymes=["EcoRI", "BamHI"])
    assert turned.product.sequence == made.product.sequence


def test_with_no_enzymes_named_the_plan_chooses_the_pair_and_the_protocol_names_it(
    puc19, source, made
):
    picked = plan_restriction(puc19, source)
    assert [one.name for one in picked.enzymes] == ["BamHI", "EcoRI"]
    assert picked.product.sequence == made.product.sequence
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


# --------------------------------------------------------------------------------------
# One enzyme, blunt ends: a backbone that closes on itself and an insert that turns round
# --------------------------------------------------------------------------------------


def test_one_enzyme_leaves_a_backbone_that_closes_on_itself_and_the_plan_dephosphorylates_it(
    blunt, puc19, gfp
):
    assert [one.name for one in blunt.enzymes] == ["SmaI"]
    # One site opens the vector, so the whole plasmid is the backbone and the gel drops nothing.
    assert [(one.name, one.length) for one in blunt.vector_pieces] == [
        ("pUC19 backbone", len(puc19))
    ]
    assert gfp.sequence in blunt.insert.bases
    assert blunt.dephosphorylates
    verdict = next(one for one in blunt.checks if one.name == "self-ligation")
    assert verdict.status == "pass"
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
    assert blunt.status == "pass"


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


# --------------------------------------------------------------------------------------
# The insert amplified rather than cut out
# --------------------------------------------------------------------------------------


def test_both_routes_to_the_insert_reach_the_same_product(tailed, made):
    assert tailed.amplicon is not None
    assert tailed.product.sequence == made.product.sequence


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


def test_a_tail_spelling_a_second_site_against_the_insert_is_refused_naming_the_enzyme(puc19, gfp):
    # An EcoRI tail ends GAATTC, and an insert beginning TAGA completes TCTAGA across the join:
    # an XbaI site no spacer can be moved off, because the spacer sits the other side of it.
    edge = SequenceRecord("TAGA" + gfp.sequence, name="GFP")
    with pytest.raises(ValueError, match=r"XbaI cuts GFP amplicon 2 time"):
        plan_restriction(puc19, edge, enzymes=["EcoRI", "XbaI"])
