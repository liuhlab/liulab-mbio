"""Accepting a destination vector, and making one that is not compatible.

Every method here is laid out from the shipped enzyme definitions: the internal stuffer puts that
enzyme's two sites where its own cut offsets ask for them, so the overhangs a test asserts are read
back off the DNA rather than copied from whatever wrote it.

The working vector is the one real record here. It is pLVX with its ccdB cassette where EGFP was,
and the gate judges the tube that opens it.
"""

from pathlib import Path
from typing import Any

import pytest

from liulab_mbio.edits import replace, rotate
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.io import read_record
from liulab_mbio.reaction import Pool, Reaction
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)
from liulab_mbio.sites import find_sites
from liulab_mbio.translate import translate
from liulab_synbio.igga.gate import check_reaction
from liulab_synbio.igga.method import IGGA, INTERFACE_OVERHANGS, Scheme
from liulab_synbio.igga.project import read_project
from liulab_synbio.igga.vector import (
    B0034,
    BIOBRICK_SCAR,
    CCDB,
    CCDB_PAYLOAD,
    J23119,
    Destination,
    cargo_candidates,
    cargo_enzyme,
    ccdb_cassette,
    destination_vector,
    domesticate_vector,
    donor_cassette,
    round_cassette,
    working_vector,
)

#: The jobs the test schemes give their enzymes. All four are free of sites in pUC19.
INTERNAL = "BbsI"
EXTERNAL = "PaqCI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"

#: The overhang a part enters on, and the cloning scar every part's 3' end leaves.
ENTRY = ("CTCC",)
SCAR = "AGCG"

#: A stretch of pUC19 annotated by no feature, between the lac promoter and the origin.
GAP = (600, 700)

#: The enzyme that admits cargo to pLVX, which carries no site of it: #256's choice, which
#: `test_the_cargo_enzyme_search_leaves_plvx_only_paqci` measures again.
CARGO = "PaqCI"

#: Where the project the gate reads lives. A digest is judged on the acting enzymes alone, so any
#: project answers, and this is the one the suite already has.
DEMO = Path(__file__).parents[3] / "docs" / "examples" / "ap1-library" / "project.json"


def pad(length: int) -> str:
    """`length` bases spelling no site any of the methods' enzymes reads."""
    return ("TA" * length)[:length]


def core(cutter: Enzyme, chopper: Enzyme, *, scar: str = SCAR) -> str:
    """A shared internal stuffer core, carrying both of the cuts that open a vector.

    The reverse cut reaches back into the prefix and the forward one lands on `scar`, so the piece
    excised from a vector is bounded by the two overhangs a first-round part enters and leaves on.
    Its length is chosen to leave the retained region a whole number of codons.
    """
    reach = pad(cutter.top_cut - len(cutter.site))
    head = (
        pad(cutter.bottom_cut - len(cutter.site) - len(scar))
        + reverse_complement(cutter.site)
        + pad(1)
        + chopper.site
    )
    end = cutter.site + reach + scar
    fill = (len(scar) - len(ENTRY[0]) - len(head) - len(end)) % 3
    return head + pad(fill) + end


def external_5(overhang: str, cutter: Enzyme, chopper: Enzyme) -> str:
    """A 5' external stuffer whose `cutter` cut leaves `overhang` as its last bases."""
    reach = pad(cutter.top_cut - len(cutter.site))
    return pad(3) + chopper.site + pad(3) + cutter.site + reach + overhang


def external_3(cutter: Enzyme, chopper: Enzyme) -> str:
    """A 3' external stuffer whose `cutter` cut leaves the cloning scar as its first bases."""
    reach = pad(cutter.top_cut - len(cutter.site))
    return SCAR + reach + reverse_complement(cutter.site) + pad(3) + chopper.site + pad(3)


def scheme(*, internal: str = INTERNAL, **changes: Any) -> Scheme:
    """A valid method, with any field replaced."""
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    fields: dict[str, Any] = {
        "internal_enzyme": internal,
        "external_enzyme": EXTERNAL,
        "blunt_enzymes": (CORE_CHOPPER, FLANK_CHOPPER),
        "internal_stuffer_prefix": ENTRY[0],
        "internal_stuffer_core": core(get_enzyme(internal), get_enzyme(CORE_CHOPPER)),
        "external_stuffer_5": external_5(ENTRY[0], cutter, chopper),
        "external_stuffer_3": external_3(cutter, chopper),
        "cloning_scar": SCAR,
    }
    fields.update(changes)
    return Scheme("test", **fields)


def carrier(made: Scheme, *, flank: int = 80) -> SequenceRecord:
    """A circular vector already carrying the method's internal stuffer."""
    return SequenceRecord(
        pad(flank) + made.internal_stuffer + pad(flank),
        topology="circular",
        name="carrier",
    )


def test_a_vector_carrying_a_stuffer_is_accepted_unchanged():
    made = scheme()
    vector = carrier(made)

    taken = destination_vector(vector, made)

    assert taken.record is vector
    assert taken.edit is None
    excised = vector.extract(taken.stuffer)
    assert excised.startswith(made.entry_overhang)
    assert vector.sequence[taken.stuffer.end : taken.stuffer.end + len(SCAR)] == made.scar_overhang


def test_a_stuffer_across_the_origin_is_found_where_it_lies():
    made = scheme()
    vector = rotate(carrier(made), 90)

    taken = destination_vector(vector, made)

    assert taken.stuffer.end > len(vector)
    assert taken.record is vector
    assert vector.extract(taken.stuffer).startswith(made.entry_overhang)


def test_puc19_is_made_compatible_at_a_named_span(puc19):
    made = scheme()
    stuffer = made.internal_stuffer

    taken = destination_vector(puc19, made, site=GAP)

    assert len(taken.record) == len(puc19) + len(stuffer)
    assert taken.record.sequence[GAP[0] : GAP[0] + len(stuffer)] == stuffer
    assert taken.edit is not None
    # The retrofit is held to the rule that accepts a vector already carrying one.
    again = destination_vector(taken.record, made)
    assert again.record is taken.record
    assert again.stuffer == taken.stuffer


def test_an_edit_names_the_features_and_binding_sites_it_changed():
    made = scheme()
    at = 40
    spanning = Feature("spanning", "misc_feature", (Segment(at - 10, at + 10),))
    primer = Primer("reader", pad(20), binding_sites=(BindingSite(at - 5, at + 5, Strand.FORWARD),))
    vector = SequenceRecord(pad(200), topology="circular", features=(spanning,), primers=(primer,))

    taken = destination_vector(vector, made, site=(at, at + 10))

    assert taken.edit is not None
    assert taken.edit.changed == (spanning,)
    assert taken.edit.dropped_sites == ((primer, primer.binding_sites[0]),)


def test_an_edit_across_the_origin_leaves_what_it_annotates_reading_the_same():
    made = scheme()
    length = 200
    wrapping = Feature("wrapping", "misc_feature", (Segment(length - 8, length + 8),))
    primer = Primer(
        "reader", pad(16), binding_sites=(BindingSite(length - 8, length + 8, Strand.FORWARD),)
    )
    vector = SequenceRecord(
        pad(length), topology="circular", features=(wrapping,), primers=(primer,)
    )
    before = vector.extract(wrapping)

    taken = destination_vector(vector, made, site=(100, 110))

    moved = taken.record.features[0].segments[0]
    assert moved.end > len(taken.record)
    assert taken.record.extract(taken.record.features[0]) == before
    site = taken.record.primers[0].binding_sites[0]
    assert (site.start, site.end) == (moved.start, moved.end)


def test_a_stuffer_that_is_not_whole_codons_inside_a_coding_sequence_is_refused(puc19):
    made = scheme()

    with pytest.raises(ValueError, match="reading frame"):
        destination_vector(puc19, made, site="MCS")


def test_an_enzyme_site_outside_the_stuffer_is_refused(puc19):
    made = scheme(internal="BsaI")

    # pUC19's own BsaI site is at 1765, and the stuffer put at 600 moved it along.
    with pytest.raises(ValueError, match="BsaI reads a site at 1796 on the reverse strand"):
        destination_vector(puc19, made, site=GAP)


def test_a_vector_with_no_stuffer_and_no_site_named_is_refused():
    made = scheme()

    with pytest.raises(ValueError, match="name the site to put one at"):
        destination_vector(SequenceRecord(pad(200), topology="circular"), made)


def test_a_site_naming_no_feature_is_refused():
    made = scheme()

    with pytest.raises(ValueError, match="annotates no feature called 'nope'"):
        destination_vector(SequenceRecord(pad(200), topology="circular"), made, site="nope")


def test_the_cargo_enzyme_search_leaves_plvx_only_paqci(plvx: SequenceRecord) -> None:
    """The measurement `docs/research/working-vector-plvx-tetone.md` made, as the rule reads it."""
    chosen = cargo_enzyme([plvx])
    assert chosen.enzyme is not None
    assert chosen.enzyme.name == "PaqCI"
    assert [one.name for one in chosen.search.free] == ["PaqCI"]
    assert chosen.check.status == "pass"
    assert "BsmBI" not in cargo_candidates()
    assert "Esp3I" not in cargo_candidates()


def test_a_three_base_overhang_is_blocked_before_its_sites_are_counted(
    plvx: SequenceRecord,
) -> None:
    blocked = {one.enzyme.name: one for one in cargo_enzyme([plvx]).search.blocked}
    assert not blocked["SapI"].sites
    assert "3 bases" in blocked["SapI"].reason


def test_bsai_is_rejected_naming_the_two_ltr_sites(plvx: SequenceRecord) -> None:
    """#256's collision: BsaI shared the final-assembly pot with the LTRs at 454 and 7112."""
    blocked = {one.enzyme.name: one for one in cargo_enzyme([plvx]).search.blocked}
    assert [site.start for site in blocked["BsaI"].sites] == [
        454,
        3724,
        5730,
        6472,
        7112,
        8720,
    ]


def test_the_two_ltr_bsai_sites_are_refused_as_past_every_oligos_reach(
    plvx: SequenceRecord,
) -> None:
    held = domesticate_vector(plvx, (get_enzyme("BsaI"),))
    assert [site.start for site in held.unreachable] == [454, 7112]
    assert held.check.status == "fail"
    assert "No oligo reaches 454, 7112" in held.check.detail
    # #289: the two coding sites go, so what is left is the four in no coding sequence.
    assert [site.start for site in held.remaining] == [454, 3724, 6472, 7112]


def test_a_vector_free_of_the_enzymes_passes_domestication_unchanged(
    puc19: SequenceRecord,
) -> None:
    held = domesticate_vector(puc19, (get_enzyme("PaqCI"),))
    assert held.check.status == "pass"
    assert held.record.sequence == puc19.sequence
    assert held.unreachable == ()
    assert held.remaining == ()


def test_an_empty_search_is_a_finding_and_not_an_error(plvx: SequenceRecord) -> None:
    chosen = cargo_enzyme([plvx], candidates=("BsaI", "BbsI"))
    assert chosen.enzyme is None
    assert chosen.check.status == "fail"
    assert "cannot be one pot" in chosen.check.detail
    assert {site.start for site in chosen.sites} >= {454, 7112, 24}


def test_the_ccdb_payload_is_the_359_bases_287_decided():
    """#287: promoter, ribosome binding site, scar and the toxin, clean but for one SrfI."""
    assert CCDB_PAYLOAD == J23119 + B0034 + BIOBRICK_SCAR + CCDB
    assert len(CCDB_PAYLOAD) == 359
    assert translate(CCDB) == (
        "MQFKVYTYKRESRYRLFVDVQSDIIDTPGRRMVIPLASARLLSDKVSRELYPVVHIGDESWRMMTTDMASVP"
        "VSVIGEEVADLSHRENDIKNAINLMFWGI*"
    )
    read = SequenceRecord(CCDB_PAYLOAD)
    assert find_sites(read, ("BsaI", "BsmBI", "BbsI", "PaqCI", "SapI", "PmeI")) == ()
    assert [site.start for site in find_sites(read, CORE_CHOPPER)] == [133]


def test_a_skipped_side_of_the_working_cassette_carries_the_end_its_part_would_have():
    start, end = INTERFACE_OVERHANGS["working cassette"]
    both = ccdb_cassette(get_enzyme(CARGO))
    neither = ccdb_cassette(get_enzyme(CARGO), n_part=False, c_part=False)

    assert both.bases.startswith(IGGA.entry_overhang)
    assert both.bases.endswith(IGGA.scar_overhang)
    assert neither.bases.startswith(start + "GG")
    assert neither.bases.endswith("GG" + end)
    assert neither.bases[len(start) + 2 : -len(end) - 2] == both.bases
    assert both.free_of == (get_enzyme(CARGO), *IGGA.reserved_enzymes)


def test_an_enzyme_reading_a_site_inside_the_payload_is_refused():
    """SrfI's site is the one the payload carries, so it cannot be the enzyme that releases it."""
    with pytest.raises(ValueError, match="reads 3 site\\(s\\) in this ccdB cassette"):
        ccdb_cassette(get_enzyme(CORE_CHOPPER))


@pytest.fixture(scope="module")
def working(plvx: SequenceRecord) -> Destination:
    """pLVX as a working vector: no BsmBI left, and the ccdB cassette where EGFP was.

    `docs/research/working-vector-plvx-tetone.md` counts two BsmBI sites. The one in PuroR is
    coding, so `domesticate_vector` takes it out; #256 left the other, in the hPGK promoter, for
    a person to edit, and this stands in for that edit.
    """
    held = domesticate_vector(plvx, (get_enzyme("BsmBI"),))
    (left,) = held.remaining
    cleared, _ = replace(
        held.record, left.start, left.start + len(left.enzyme.site), pad(len(left.enzyme.site))
    )
    return destination_vector(cleared, IGGA, site="EGFP", cassette=ccdb_cassette(get_enzyme(CARGO)))


def test_a_working_vector_gives_up_its_ccdb_cassette_to_the_cargo_enzyme(working: Destination):
    excised = working.record.extract(working.stuffer)
    after = working.stuffer.end

    assert excised.startswith(IGGA.entry_overhang)
    assert CCDB_PAYLOAD in excised
    assert working.record.sequence[after : after + len(IGGA.scar_overhang)] == IGGA.scar_overhang
    assert working.edit is not None
    # The working vector is held to the rule that accepts a vector already carrying a cassette.
    again = destination_vector(working.record, IGGA, cassette=ccdb_cassette(get_enzyme(CARGO)))
    assert again.record is working.record
    assert again.stuffer == working.stuffer


def test_the_gate_passes_the_tube_the_cargo_enzyme_opens_a_working_vector_in(working: Destination):
    reaction = Reaction(
        "the final assembly",
        "digest",
        pools=[Pool("destination", [working.record])],
        enzymes=[CARGO],
    )

    judged = check_reaction(reaction, project=read_project(DEMO))

    assert [one.status for one in judged] == ["pass"]
    assert "in 2 places" in judged[0].check.detail


def test_a_working_vector_keeps_the_sites_no_tube_it_meets_bars(working: Destination):
    """BsaI never shares a tube with this vector, so its six sites are counted and left."""
    bsai = get_enzyme("BsaI")
    assert len(find_sites(working.record, bsai)) == 6
    assert bsai not in ccdb_cassette(get_enzyme(CARGO)).free_of
    assert bsai not in round_cassette(IGGA).free_of
    assert donor_cassette(IGGA).enzyme == bsai


def test_a_bsmbi_site_the_working_vector_still_reads_is_refused(plvx: SequenceRecord):
    """#256 left one site in the hPGK promoter; without it the backbone is cut twice."""
    held = domesticate_vector(plvx, (get_enzyme("BsmBI"),))

    with pytest.raises(ValueError, match="BsmBI reads a site at 4265 on the forward strand"):
        destination_vector(
            held.record, IGGA, site="EGFP", cassette=ccdb_cassette(get_enzyme(CARGO))
        )


def test_a_working_vector_chooses_its_enzyme_before_it_builds_the_cassette(plvx):
    """`working_vector` is the chain in one call: choose, write the cassette, put it in."""
    held = domesticate_vector(plvx, (get_enzyme("BsmBI"),))
    (left,) = held.remaining
    cleared, _ = replace(
        held.record, left.start, left.start + len(left.enzyme.site), pad(len(left.enzyme.site))
    )

    made = working_vector(cleared, [], scheme=IGGA, site="EGFP")

    assert made.enzyme == get_enzyme(CARGO)
    assert made.cassette.bases == ccdb_cassette(get_enzyme(CARGO)).bases
    assert made.record.extract(made.destination.stuffer).startswith(IGGA.entry_overhang)


def test_a_pot_no_candidate_is_free_of_refuses_rather_than_choosing_one(plvx):
    """The AP-1 library spells PaqCI twice, so nothing is left to admit it to pLVX."""
    product = read_record(DEMO.parent / "product.dna")

    with pytest.raises(ValueError, match="no candidate is free to admit cargo"):
        working_vector(plvx, [product], scheme=IGGA, site="EGFP")


def donor_carrier(made: Scheme, *, flank: int = 80) -> SequenceRecord:
    """A circular donor backbone: the method's stuffer between the two external stuffers.

    What the external enzyme frees from it is the same piece the internal enzyme frees from a
    round's destination, which is the point: one stuffer, two tubes.
    """
    entry, scar = made.entry_overhang, made.scar_overhang
    between = made.internal_stuffer[len(entry) : len(made.internal_stuffer) - len(scar)]
    return SequenceRecord(
        pad(flank) + made.external_stuffer_5 + between + made.external_stuffer_3 + pad(flank),
        topology="circular",
        name="donor",
    )


def test_each_cassette_bars_the_enzymes_of_its_own_tube_and_no_others():
    made = scheme()
    internal, external = made.internal, made.external
    core_chopper, flank_chopper = get_enzyme(CORE_CHOPPER), get_enzyme(FLANK_CHOPPER)

    assert round_cassette(made).free_of == (internal, core_chopper)
    assert made.blunt_for_the_destination == (core_chopper,)
    # The donor throws its backbone away, so nothing it reads outside the cassette is a defect.
    assert donor_cassette(made).free_of == ()
    assert donor_cassette(made).enzyme == external
    assert flank_chopper not in round_cassette(made).free_of


def test_a_destination_carrying_the_sites_that_release_its_cargo_is_accepted():
    """The outboard external and blunt sites are what frees the cargo once the rounds are done."""
    made = scheme()
    vector = donor_carrier(made)

    taken = destination_vector(vector, made)

    assert taken.record is vector
    assert vector.extract(taken.stuffer).startswith(made.entry_overhang)
    assert find_sites(vector, made.external)


def test_the_same_record_is_a_donor_to_the_other_tube():
    made = scheme()
    vector = donor_carrier(made)

    taken = destination_vector(vector, made, cassette=donor_cassette(made))

    end = taken.stuffer.end
    assert vector.extract(taken.stuffer).startswith(made.entry_overhang)
    assert vector.sequence[end : end + len(SCAR)] == made.scar_overhang


def test_a_refusal_names_the_tube_the_site_would_be_cut_in():
    made = scheme()
    at = 20
    vector = donor_carrier(made)
    strayed, _ = replace(vector, at, at + len(made.internal.site), made.internal.site)

    with pytest.raises(ValueError, match="acts in the digest that opens a round's destination"):
        destination_vector(strayed, made)


def test_a_donor_backbone_keeps_the_releasing_and_blunt_sites_it_is_built_on():
    """Both lie outboard of the cargo, which is the whole point of holding a part in a backbone."""
    made = scheme()
    vector = donor_carrier(made)
    held = donor_cassette(made)
    taken = destination_vector(vector, made, cassette=held)

    outboard = [
        site
        for site in find_sites(vector, (made.external, *made.blunt))
        if not vector.covers(taken.stuffer, site.span)
    ]

    assert {site.enzyme for site in outboard} == {made.external, get_enzyme(FLANK_CHOPPER)}
