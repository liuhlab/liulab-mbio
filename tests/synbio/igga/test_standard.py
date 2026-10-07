"""What the overhang standard chooses, what it charges the proteins, and what it refuses.

The method here is laid out from the shipped enzyme definitions, as `test_method.py` lays the
shipped one out, so nothing states an overhang that the DNA does not spell.
"""

from itertools import product

import pytest

from liulab_mbio.codons import amino_acid, codon_usage
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.overhangs import refusal
from liulab_mbio.sequence import reverse_complement
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.standard import (
    _attempt,
    _Reading,
    _sites,
    design_standard,
    junction_residues,
)

INTERNAL = "BsaI"
EXTERNAL = "BbsI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"

#: The overhang a part enters on, the cloning scar, and the positions a build fills. Every
#: position carries the method's own stuffers, so one entry overhang is stated.
ENTRY = "CTCC"
SCAR = "AGCG"
TWO = ("p1", "p2")

#: Two part lists whose members disagree about their junction: one ends in D and the other in K,
#: and no single codon spells both, so no overhang joins them without changing a residue.
FIRST = {"a": "MKTD", "b": "MKTK"}
SECOND = {"x": "MKTQ", "y": "WKTQ"}

#: Three N-domain endings the paper reports, and a following list far enough from E that the
#: junction's codon is cheapest charged to the N domain.
ENDINGS = {"ATF2": "MKTQED", "ATF4": "MKTQEK", "ATF6": "MKTQIA"}
FOLLOWING = {f"b{index}": "WQAAA" for index in range(5)}
LAST = {"c1": "MKTQ"}


def pad(length: int) -> str:
    """`length` bases spelling no site any of the method's enzymes reads."""
    return ("TA" * length)[:length]


def external_5(overhang: str, cutter: Enzyme, chopper: Enzyme) -> str:
    """A 5' external stuffer whose `cutter` cut leaves `overhang` as its last bases."""
    reach = pad(cutter.top_cut - len(cutter.site))
    return pad(3) + chopper.site + pad(3) + cutter.site + reach + overhang


def external_3(scar: str, cutter: Enzyme, chopper: Enzyme) -> str:
    """A 3' external stuffer whose `cutter` cut leaves `scar` as its first bases."""
    reach = pad(cutter.top_cut - len(cutter.site))
    return scar + reach + reverse_complement(cutter.site) + pad(3) + chopper.site + pad(3)


def core(cutter: Enzyme, chopper: Enzyme, scar: str = SCAR) -> str:
    """A shared internal stuffer core, carrying both of the cuts that open a stuffer.

    Each cut is placed from `cutter`'s own offsets: one cuts back into whatever prefix precedes
    the core, and one leaves `scar` as the stuffer's last bases. The filler length is what keeps
    the retained region a whole number of codons.
    """
    lead = cutter.bottom_cut - len(cutter.site) - cutter.overhang_length
    reach = cutter.top_cut - len(cutter.site)
    return (
        pad(lead)
        + reverse_complement(cutter.site)
        + pad(3)
        + chopper.site
        + pad(5)
        + cutter.site
        + pad(reach)
        + scar
    )


def scheme() -> Scheme:
    """A valid method, every position carrying its stuffers."""
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    return Scheme(
        "test",
        internal_enzyme=INTERNAL,
        external_enzyme=EXTERNAL,
        blunt_enzymes=(CORE_CHOPPER, FLANK_CHOPPER),
        internal_stuffer_prefix=pad(2) + ENTRY,
        internal_stuffer_core=core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER)),
        external_stuffer_5=external_5(ENTRY, cutter, chopper),
        external_stuffer_3=external_3(SCAR, cutter, chopper),
        cloning_scar=SCAR,
    )


def every_overhang(length: int) -> tuple[str, ...]:
    """Every overhang of this length, in no particular order."""
    return tuple("".join(bases) for bases in product("ACGT", repeat=length))


def retained_codons(made: Scheme, candidate: str) -> tuple[str, ...]:
    """The codons a first entry overhang's own bases fall in, in the stuffer the product retains.

    The terminal stuffer is never excised, so the product reads it in frame from its first base,
    and the prefix ends with the overhang admitting position one.
    """
    prefix = made.internal_stuffer_prefix
    head = prefix[: len(prefix) - len(candidate)]
    stuffer = head + candidate + made.internal_stuffer_core
    return tuple(
        stuffer[at : at + 3] for at in range((len(head) // 3) * 3, len(head) + len(candidate), 3)
    )


@pytest.fixture(scope="module")
def made():
    """A method laid out from the shipped enzyme definitions."""
    return scheme()


@pytest.fixture(scope="module")
def standard(made):
    """What the designer returns over the two part lists, with nothing pinned."""
    return design_standard(made, TWO, (FIRST, SECOND))


def test_the_standard_returns_one_overhang_a_position_and_the_method_s_scar(made, standard):
    assert len(standard.entry_overhangs) == len(TWO)
    assert standard.scar_overhang == made.cloning_scar
    assert len(standard.choices) == len(TWO) + 1


def test_every_chosen_overhang_holds_the_reused_rules(made, standard):
    chosen = (*standard.entry_overhangs, standard.scar_overhang)

    for index, overhang in enumerate(chosen):
        others = (*chosen[:index], *chosen[index + 1 :])
        assert (
            refusal(
                overhang,
                made.internal,
                taken=others,
                avoid=(made.external, *made.blunt),
            )
            is None
        )


def test_no_overhang_spells_a_stop_where_the_product_reads_through_it(standard):
    for overhang in standard.entry_overhangs:
        donated, _ = junction_residues(len(overhang))
        start = (3 - donated) % 3
        codons = [overhang[at : at + 3] for at in range(start, len(overhang), 3)]
        assert all(codon not in ("TAA", "TAG", "TGA") for codon in codons)


def test_the_published_terminal_changes_are_reproduced():
    made = scheme()

    standard = design_standard(
        made, ("N", "bZIP", "C"), (ENDINGS, FOLLOWING, LAST), pinned={"bZIP": "GGAG"}
    )

    charged = {
        one.part: (one.wild_type, one.synthesised)
        for one in standard.termini
        if one.position == "N" and one.end == "3'"
    }
    assert charged == {
        "ATF2": ("ED", "EE"),
        "ATF4": ("EK", "EE"),
        "ATF6": ("IA", "ME"),
    }


def test_a_junction_no_overhang_spells_unchanged_is_reported(standard):
    assert "p2" in standard.forced
    assert standard.cost > 0
    assert {one.part for one in standard.changes}


def test_the_first_entry_overhang_spells_no_stop_in_the_retained_stuffer(made, standard):
    spelled = retained_codons(made, standard.entry_overhangs[0])

    assert spelled
    assert all(amino_acid(one) != "*" for one in spelled)


def test_a_pinned_overhang_spelling_a_retained_stop_is_refused(made):
    with pytest.raises(ValueError, match="stop"):
        design_standard(made, TWO, (FIRST, SECOND), pinned={"p1": "AGGA"})


def test_a_pinned_overhang_is_kept_and_charges_no_part_list_its_5_terminus(made):
    starts_with_s = {"a": "SKTD", "b": "SKTK"}

    pinned = design_standard(made, TWO, (starts_with_s, SECOND), pinned={"p1": "CTCC"})

    assert pinned.entry_overhangs[0] == "CTCC"
    assert not [one for one in pinned.termini if one.position == "p1" and one.end == "5'"]


def test_no_cheaper_standard_is_available(made, standard):
    """Price every set the rules allow, and show none beats what the designer returned.

    A junction's cost never falls, so a set whose first junctions already cost more than the
    answer cannot win and is not priced further. The cost model itself is pinned by
    `test_the_published_terminal_changes_are_reproduced`.
    """
    lists = (FIRST, SECOND)
    enzyme, avoid = made.internal, (made.external, *made.blunt)
    usage = codon_usage()
    sites = _sites(
        made, TWO, lists, {}, junction_residues(enzyme.overhang_length)[1], enzyme.overhang_length
    )
    pool = [
        one
        for one in every_overhang(enzyme.overhang_length)
        if refusal(one, enzyme, avoid=avoid) is None
    ]
    priced = [
        {one: read.cost for one in pool if isinstance(read := _attempt(site, one, usage), _Reading)}
        for site in sites[:-1]
    ]

    scar = made.cloning_scar
    cheapest = None
    for first, first_cost in priced[0].items():
        if first_cost > standard.cost:
            continue
        if refusal(first, enzyme, taken=(scar,), avoid=avoid) is not None:
            continue
        for second, second_cost in priced[1].items():
            total = first_cost + second_cost
            if total > standard.cost:
                continue
            if refusal(second, enzyme, taken=(scar, first), avoid=avoid) is not None:
                continue
            cheapest = total if cheapest is None else min(cheapest, total)

    assert cheapest == standard.cost


def test_a_pinned_overhang_is_held_to_every_rule(made):
    with pytest.raises(ValueError, match="palindrome"):
        design_standard(made, TWO, (FIRST, SECOND), pinned={"p1": "GGCC"})


def test_a_pinned_name_that_is_not_a_position_is_refused(made):
    with pytest.raises(ValueError, match="not positions of this build"):
        design_standard(made, TWO, (FIRST, SECOND), pinned={"nowhere": "AGGT"})


def test_part_lists_must_match_the_positions(made):
    with pytest.raises(ValueError, match="part list"):
        design_standard(made, TWO, (FIRST,))


def test_a_protein_shorter_than_a_junction_spells_is_refused(made):
    with pytest.raises(ValueError, match="amino acid"):
        design_standard(made, TWO, ({"tiny": "M"}, SECOND))


def test_the_rejection_trail_names_the_rule_that_refused_a_candidate(made, standard):
    wanted = standard.entry_overhangs[0]

    # A pinned junction settles first, so pinning the second position to what the first would
    # otherwise take leaves the first to refuse it. The method fixed that overhang, so the
    # trail says it is reserved rather than that it duplicates a junction free to move.
    again = design_standard(made, TWO, (FIRST, SECOND), pinned={"p2": wanted})
    trail = [one for choice in again.choices for one in choice.rejected]

    assert ("reserved", wanted) in {(one.rule, one.overhang) for one in trail}
    assert all(one.detail for one in trail)


def test_an_overhang_the_caller_reserves_is_held_out_and_refused_as_reserved(made, standard):
    wanted = standard.entry_overhangs[0]

    held = design_standard(made, TWO, (FIRST, SECOND), reserved=[wanted])
    refused = refusal(wanted, made.internal, reserved=[wanted])

    assert wanted not in (*held.entry_overhangs, held.scar_overhang)
    assert refused is not None
    assert refused.rule == "reserved"
