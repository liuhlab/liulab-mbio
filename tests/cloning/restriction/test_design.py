"""Choosing the pair, over records small enough to read.

Every record here is written in code, so what refuses a pair and what ranks one above another is
on the page beside the assertion. The plan seam is in `test_plan.py`.
"""

import dataclasses

import pytest

from liulab_mbio.cloning.restriction.design import candidates, choose_pair, refusal
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord

#: Bases spelling no site of any enzyme named below, for padding a record written in code.
FILLER = "ACGT" * 6

#: An enzyme wanting its own buffer and its own temperature, so that a pair holding it warns on
#: one verdict and carries no verdict on the other. Nothing shipped is Type II and warm at once.
WARM = Enzyme(
    "WarmI",
    "TTCGAA",
    top_cut=1,
    bottom_cut=5,
    incubation_celsius=55,
    supplied_buffer="NEBuffer r3.1",
)

#: The same enzyme in the buffer and at the temperature everything else here wants.
COOL = dataclasses.replace(
    WARM, name="CoolI", incubation_celsius=37, supplied_buffer="rCutSmart Buffer"
)


def plasmid(sequence: str, name: str, features: tuple[Feature, ...] = ()) -> SequenceRecord:
    """A circular record written in code."""
    return SequenceRecord(sequence, topology="circular", name=name, features=features)


@pytest.fixture(scope="module")
def vector() -> SequenceRecord:
    """A vector whose cloning site reads EcoRI, then BamHI, then SalI close behind it."""
    return plasmid(
        FILLER * 6 + "GAATTC" + FILLER + "GGATCC" + "ACGTACGT" + "GTCGAC" + FILLER * 6, "pVec"
    )


@pytest.fixture(scope="module")
def source() -> SequenceRecord:
    """The same three sites, with something annotated between the first two."""
    bases = FILLER * 4 + "GAATTC" + FILLER * 2 + "GGATCC" + "ACGTACGT" + "GTCGAC" + FILLER * 4
    cargo = Feature("cargo", "misc_feature", (Segment(len(FILLER * 4) + 6, len(FILLER * 6) + 6),))
    return plasmid(bases, "pIns", (cargo,))


def test_the_chosen_pair_moves_what_the_source_annotates_rather_than_the_least_vector(
    vector, source
):
    # BamHI and SalI lie eight bases apart, so that pair gives up less of the vector than any
    # other -- and carries none of what the source annotates, which is what the cloning is for.
    assert [one.name for one in choose_pair(vector, source).enzymes] == ["BamHI", "EcoRI"]
    bare = dataclasses.replace(source, features=())
    assert [one.name for one in choose_pair(vector, bare).enzymes] == ["BamHI", "SalI"]


def test_a_pair_the_verdicts_rank_down_loses_to_one_they_do_not():
    # One layout, two enzymes differing in nothing but their buffer and their temperature: the
    # pair holding the warm one gives up the least vector and is passed over all the same.
    site = "GAATTC" + "ACGTACGT" + "TTCGAA" + FILLER + "GGATCC"
    into = plasmid(FILLER * 6 + site + FILLER * 6, "pWarmVec")
    holder = plasmid(FILLER * 4 + site + FILLER * 4, "pWarmIns")
    warm = choose_pair(into, holder, enzymes=["EcoRI", WARM, "BamHI"])
    assert [one.name for one in warm.enzymes] == ["EcoRI", "BamHI"]
    cool = choose_pair(into, holder, enzymes=["EcoRI", COOL, "BamHI"])
    assert [one.name for one in cool.enzymes] == ["EcoRI", "CoolI"]


def test_an_enzyme_reading_no_site_in_the_vector_refuses_every_pair_it_is_in(vector, source):
    holder = plasmid(source.sequence + "AAGCTT" + FILLER, "pIns2")
    chosen = choose_pair(vector, holder, enzymes=["EcoRI", "BamHI", "HindIII"])
    assert [(one.rule, one.names) for one in chosen.refusals] == [
        ("vector site", "BamHI and HindIII"),
        ("vector site", "EcoRI and HindIII"),
    ]
    assert "HindIII cuts pVec 0 time(s)" in chosen.refusals[0].detail


def test_an_insert_neither_route_reaches_says_what_a_codon_change_would_cost(vector):
    # One EcoRI site and no BamHI site: too few sites to cut the insert out, and one too many to
    # amplify it with the sites on its tails.
    coding = "ATG" + "GAATTC" + "AAACCTTTAGGCATT" + "TAA"
    gene = SequenceRecord(
        "ACGTAC" + coding + "GTACGT",
        name="cargo",
        features=(Feature("cargo", "CDS", (Segment(6, 6 + len(coding)),)),),
    )
    pair = [get_enzyme("EcoRI"), get_enzyme("BamHI")]
    inside = refusal(vector, gene, pair)
    assert inside is not None
    assert (inside.rule, inside.domesticable) == ("insert site", True)
    assert "1 EcoRI site (at 9) and no BamHI site" in inside.detail
    assert "EcoRI at 9 by TTC to TTT in cargo" in inside.detail
    assert [(one.old_codon, one.new_codon) for one in inside.changes] == [("TTC", "TTT")]

    # The same bases annotating no coding sequence: a flat refusal, because taking the site out
    # would change what the record spells.
    outside = refusal(vector, dataclasses.replace(gene, features=()), pair)
    assert outside is not None
    assert (outside.rule, outside.domesticable) == ("insert site", False)
    assert "no synonymous codon change reaches EcoRI at 9" in outside.detail


def test_a_backbone_that_closes_on_itself_and_ends_that_do_not_anneal_name_their_rules():
    # SalI and XhoI both leave a 5' TCGA, so the backbone they open closes again with nothing in.
    into = plasmid(FILLER * 6 + "GTCGAC" + FILLER + "CTCGAG" + FILLER * 6, "pTcga")
    holder = plasmid(FILLER * 4 + "GTCGAC" + FILLER * 2 + "CTCGAG" + FILLER * 4, "pTcgaIns")
    pair = [get_enzyme("SalI"), get_enzyme("XhoI")]
    closed = refusal(into, holder, pair, choosing=True)
    assert closed is not None
    assert closed.rule == "backbone"
    assert "anneal to each other" in closed.detail
    # The rule is the chooser's alone: naming that pair plans it, with a dephosphorylation.
    assert refusal(into, holder, pair) is None

    # BsaI cuts outside its own site, so what it leaves is whatever the record spells there.
    apart = refusal(
        plasmid("GGTCTCAAAAA" + FILLER * 8 + "GAATTC" + FILLER * 2, "pBsaI-A"),
        plasmid("GGTCTCAGGGG" + FILLER * 2 + "GAATTC" + FILLER * 8, "pBsaI-G"),
        [get_enzyme("BsaI"), get_enzyme("EcoRI")],
    )
    assert apart is not None
    assert apart.rule == "ends"
    assert "do not anneal" in apart.detail


def test_where_no_pair_can_be_found_the_reasons_are_named(vector, source):
    with pytest.raises(ValueError, match="every pair was refused, 1 on the vector site rule"):
        choose_pair(vector, source, enzymes=["HindIII", "NotI"])


def test_a_usable_pair_is_refused_by_nothing_and_every_candidate_puts_its_own_site_back(
    vector, source
):
    assert refusal(vector, source, [get_enzyme("EcoRI"), get_enzyme("BamHI")]) is None
    assert all(one.type == "II" for one in candidates())
    assert {"EcoRI", "BamHI"} <= {one.name for one in candidates()}
