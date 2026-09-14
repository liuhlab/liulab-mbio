"""What the overhang standard chooses, what it charges the proteins, and what it refuses.

The schemes here are laid out from the shipped enzyme definitions, as `test_scheme.py` lays its
own out, so nothing states an overhang that the DNA does not spell.
"""

from collections.abc import Sequence
from itertools import product
from pathlib import Path

import pytest

from liulab_mbio.codons import amino_acid, codon_usage
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.goldengate.design import refusal
from liulab_mbio.library.scheme import Position, Scheme, read_scheme
from liulab_mbio.library.standard import (
    _attempt,
    _Reading,
    _sites,
    design_standard,
    junction_residues,
)
from liulab_mbio.sequence import reverse_complement

#: The paper's scheme, as a user supplies one.
EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "protein-library" / "scheme.json"

INTERNAL = "BsaI"
EXTERNAL = "BbsI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"

#: Entry overhangs and a cloning scar, to lay out a scheme of one to three positions.
MORE = ("CTCC", "GGAG", "CCGA")
SCAR = "AGCG"

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
    """`length` bases spelling no site any of the scheme's enzymes reads."""
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


def scheme(overhangs: Sequence[str] = MORE) -> Scheme:
    """A valid scheme of one position per overhang given."""
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    made = [
        Position(
            f"p{index + 1}",
            internal_stuffer_prefix=(pad(2) if index == len(overhangs) - 1 else "")
            + overhangs[(index + 1) % len(overhangs)],
            external_stuffer_5=external_5(overhang, cutter, chopper),
            external_stuffer_3=external_3(SCAR, cutter, chopper),
        )
        for index, overhang in enumerate(overhangs)
    ]
    return Scheme(
        "test",
        positions=tuple(made),
        internal_enzyme=INTERNAL,
        external_enzyme=EXTERNAL,
        blunt_enzymes=(CORE_CHOPPER, FLANK_CHOPPER),
        internal_stuffer_core=core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER)),
        cloning_scar=SCAR,
        barcode_length=11,
    )


def every_overhang(length: int) -> tuple[str, ...]:
    """Every overhang of this length, in no particular order."""
    return tuple("".join(bases) for bases in product("ACGT", repeat=length))


def retained_codons(made: Scheme, candidate: str) -> tuple[str, ...]:
    """The codons a first entry overhang's own bases fall in, in the stuffer the product retains.

    The terminal stuffer is never excised, so the product reads it in frame from its first base,
    and the prefix ends with the overhang admitting position one.
    """
    prefix = made.positions[-1].internal_stuffer_prefix
    head = prefix[: len(prefix) - len(candidate)]
    stuffer = head + candidate + made.internal_stuffer_core
    return tuple(
        stuffer[at : at + 3] for at in range((len(head) // 3) * 3, len(head) + len(candidate), 3)
    )


def test_the_standard_returns_one_overhang_a_position_and_the_scheme_s_scar():
    made = scheme(MORE[:2])

    standard = design_standard(made, (FIRST, SECOND))

    assert len(standard.entry_overhangs) == made.position_count
    assert standard.scar_overhang == made.cloning_scar
    assert len(standard.choices) == made.position_count + 1


def test_every_chosen_overhang_holds_the_reused_rules():
    made = scheme(MORE[:2])
    standard = design_standard(made, (FIRST, SECOND))
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


def test_no_overhang_spells_a_stop_where_the_product_reads_through_it():
    made = scheme(MORE[:2])
    standard = design_standard(made, (FIRST, SECOND))

    for overhang in standard.entry_overhangs:
        donated, _ = junction_residues(len(overhang))
        start = (3 - donated) % 3
        codons = [overhang[at : at + 3] for at in range(start, len(overhang), 3)]
        assert all(codon not in ("TAA", "TAG", "TGA") for codon in codons)


def test_the_published_terminal_changes_are_reproduced():
    made = read_scheme(EXAMPLE)

    standard = design_standard(made, (ENDINGS, FOLLOWING, LAST), pinned={"bZIP": "GGAG"})

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


def test_a_junction_no_overhang_spells_unchanged_is_reported():
    made = scheme(MORE[:2])

    standard = design_standard(made, (FIRST, SECOND))

    assert "p2" in standard.forced
    assert standard.cost > 0
    assert {one.part for one in standard.changes}


def test_the_first_entry_overhang_spells_no_stop_in_the_retained_stuffer():
    made = scheme(MORE[:2])

    standard = design_standard(made, (FIRST, SECOND))

    spelled = retained_codons(made, standard.entry_overhangs[0])
    assert spelled
    assert all(amino_acid(one) != "*" for one in spelled)


def test_a_retained_stop_names_its_rule_in_the_trail():
    made = scheme(MORE[:2])

    standard = design_standard(made, (FIRST, SECOND))
    refused = {one.overhang for one in standard.choices[0].rejected if one.rule == "stop"}

    assert refused
    for overhang in refused:
        assert any(amino_acid(one) == "*" for one in retained_codons(made, overhang))


def test_a_pinned_overhang_spelling_a_retained_stop_is_refused():
    made = scheme(MORE[:2])

    with pytest.raises(ValueError, match="stop"):
        design_standard(made, (FIRST, SECOND), pinned={"p1": "AGGA"})


def test_the_vector_junction_charges_no_part_list_its_5_terminus():
    made = scheme(MORE[:2])
    starts_with_s = {"a": "SKTD", "b": "SKTK"}

    standard = design_standard(made, (starts_with_s, SECOND), pinned={"p1": "CTCC"})

    assert standard.entry_overhangs[0] == "CTCC"
    assert not [one for one in standard.termini if one.position == "p1" and one.end == "5'"]


def test_no_cheaper_standard_is_available():
    """Price every set the rules allow, and show none beats what the designer returned.

    A junction's cost never falls, so a set whose first junctions already cost more than the
    answer cannot win and is not priced further. The cost model itself is pinned by
    `test_the_published_terminal_changes_are_reproduced`.
    """
    made = scheme(MORE[:2])
    lists = (FIRST, SECOND)
    standard = design_standard(made, lists)
    enzyme, avoid = made.internal, (made.external, *made.blunt)
    usage = codon_usage()
    sites = _sites(
        made, lists, {}, junction_residues(enzyme.overhang_length)[1], enzyme.overhang_length
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


def test_a_pinned_overhang_is_kept():
    made = scheme(MORE[:2])

    standard = design_standard(made, (FIRST, SECOND), pinned={"p1": "CTCC"})

    assert standard.entry_overhangs[0] == "CTCC"


def test_a_pinned_overhang_is_held_to_every_rule():
    made = scheme(MORE[:2])

    with pytest.raises(ValueError, match="palindrome"):
        design_standard(made, (FIRST, SECOND), pinned={"p1": "GGCC"})


def test_a_pinned_name_that_is_not_a_position_is_refused():
    made = scheme(MORE[:2])

    with pytest.raises(ValueError, match="not positions of this scheme"):
        design_standard(made, (FIRST, SECOND), pinned={"nowhere": "AGGT"})


def test_part_lists_must_match_the_positions():
    made = scheme(MORE[:2])

    with pytest.raises(ValueError, match="part list"):
        design_standard(made, (FIRST,))


def test_a_protein_shorter_than_a_junction_spells_is_refused():
    made = scheme(MORE[:2])

    with pytest.raises(ValueError, match="amino acid"):
        design_standard(made, ({"tiny": "M"}, SECOND))


def test_the_rejection_trail_names_the_rule_that_refused_a_candidate():
    made = scheme(MORE[:2])
    wanted = design_standard(made, (FIRST, SECOND)).entry_overhangs[0]

    # A pinned junction settles first, so pinning the second position to what the first would
    # otherwise take leaves the first to refuse it as a repeat and say so.
    standard = design_standard(made, (FIRST, SECOND), pinned={"p2": wanted})
    trail = [one for choice in standard.choices for one in choice.rejected]

    assert ("repeat", wanted) in {(one.rule, one.overhang) for one in trail}
    assert all(one.detail for one in trail)
