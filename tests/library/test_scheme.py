"""What a scheme derives from its enzymes and stuffers, and what it refuses.

Every scheme built here is laid out from the shipped enzyme definitions: a stuffer puts the
enzyme's site at the offset its own cut offsets ask for, so an overhang a test asserts is read
back off the DNA rather than copied from whatever wrote it.
"""

from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.library.scheme import Position, Scheme, read_scheme
from liulab_mbio.sequence import reverse_complement
from liulab_mbio.translate import translate

#: The paper's scheme, as a user supplies one.
EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "protein-library" / "scheme.json"

#: The four enzymes the paper's scheme names, by their jobs.
INTERNAL = "BsaI"
EXTERNAL = "BbsI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"

#: Four entry overhangs and a cloning scar, to build test schemes of one to four positions.
MORE = ("CTCC", "GGAG", "CCGA", "TTAC")
SCAR = "AGCG"
BARCODE = 11

#: The spacer before the reverse cut. It falls on a codon boundary of the retained stuffer, where
#: `pad`'s T would spell a stop with the site that follows it.
LEAD = "C"


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
        LEAD * lead
        + reverse_complement(cutter.site)
        + pad(3)
        + chopper.site
        + pad(5)
        + cutter.site
        + pad(reach)
        + scar
    )


def positions(overhangs: Sequence[str] = MORE[:3], scar: str = SCAR) -> tuple[Position, ...]:
    """One position an entry overhang, each admitting the next and the last admitting the first."""
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    made: list[Position] = []
    for index, overhang in enumerate(overhangs):
        following = overhangs[(index + 1) % len(overhangs)]
        terminal = index == len(overhangs) - 1
        made.append(
            Position(
                f"p{index + 1}",
                internal_stuffer_prefix=pad(2) + following if terminal else following,
                external_stuffer_5=external_5(overhang, cutter, chopper),
                external_stuffer_3=external_3(scar, cutter, chopper),
            )
        )
    return tuple(made)


def scheme(**changes: Any) -> Scheme:
    """A valid scheme, with any field replaced."""
    fields: dict[str, Any] = {
        "positions": positions(),
        "internal_enzyme": INTERNAL,
        "external_enzyme": EXTERNAL,
        "blunt_enzymes": (CORE_CHOPPER, FLANK_CHOPPER),
        "internal_stuffer_core": core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER)),
        "cloning_scar": SCAR,
        "barcode_length": BARCODE,
    }
    fields.update(changes)
    return Scheme("test", **fields)


def test_the_worked_example_derives_what_a_build_reads_off_it():
    made = read_scheme(EXAMPLE)
    terminal = made.positions[-1]
    # Only the terminal stuffer is read in the product. Every other one is excised by the round
    # that opens it, which leaves nothing of it behind but the overhang at either end.
    retained = made.internal_stuffer(-1)

    assert made.entry_overhangs == ("CTCC", "GGAG", "CCGA")
    assert made.scar_overhang == "AGCG"
    assert tuple(position.name for position in made.positions) == ("N", "bZIP", "C")
    assert made.internal.site == "GGTCTC"
    assert made.external.site == "GAAGAC"
    assert tuple(enzyme.end for enzyme in made.blunt) == ("blunt", "blunt")
    assert made.retained_length % 3 == 0
    assert len(terminal.internal_stuffer_prefix) > len(made.entry_overhang(0))
    assert "*" not in translate(retained[: len(retained) // 3 * 3])
    assert made.barcode_length == 11
    assert made.source


@pytest.mark.parametrize("count", [1, 2, 4])
def test_any_number_of_positions_validates(count: int):
    made = scheme(positions=positions(MORE[:count]))

    assert made.position_count == count
    assert made.entry_overhangs == MORE[:count]
    assert made.scar_overhang == SCAR


def test_a_prefix_not_ending_with_the_next_entry_overhang_is_refused():
    made = positions()
    broken = (replace(made[0], internal_stuffer_prefix="TTTT"), *made[1:])

    with pytest.raises(ValueError, match="internal-stuffer-prefix"):
        scheme(positions=broken)


def test_a_5_external_stuffer_the_external_enzyme_does_not_cut_is_refused():
    made = positions()
    broken = (replace(made[0], external_stuffer_5=pad(20)), *made[1:])

    with pytest.raises(ValueError, match="external-stuffer-5"):
        scheme(positions=broken)


def test_a_3_external_stuffer_yielding_another_overhang_is_refused():
    made = positions()
    odd = external_3("TTTT", get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER))
    broken = (replace(made[0], external_stuffer_3=odd), *made[1:])

    with pytest.raises(ValueError, match="external-stuffer-3"):
        scheme(positions=broken)


def test_a_barcode_not_completing_a_codon_with_the_scar_is_refused():
    with pytest.raises(ValueError, match="barcode-frame"):
        scheme(barcode_length=BARCODE - 1)


def test_a_terminal_prefix_leaving_the_product_out_of_frame_is_refused():
    made = positions()
    broken = (*made[:-1], replace(made[-1], internal_stuffer_prefix=MORE[0]))

    with pytest.raises(ValueError, match="terminal-frame"):
        scheme(positions=broken)


def test_the_internal_enzyme_in_an_external_stuffer_is_refused():
    made = positions()
    trespass = get_enzyme(INTERNAL).site + made[0].external_stuffer_5
    broken = (replace(made[0], external_stuffer_5=trespass), *made[1:])

    with pytest.raises(ValueError, match="enzyme-regions"):
        scheme(positions=broken)


def test_a_blunt_enzyme_chopping_nothing_is_refused():
    with pytest.raises(ValueError, match="enzyme-regions"):
        scheme(blunt_enzymes=(CORE_CHOPPER, FLANK_CHOPPER, "EcoRI"))


def test_a_scheme_missing_a_key_is_refused():
    data = {"name": "thin"}

    with pytest.raises(ValueError, match="missing barcode_length"):
        Scheme.from_dict(data)


def test_an_internal_stuffer_leaving_no_scar_overhang_is_refused():
    whole = core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER))
    without = whole[: -(len(SCAR) + 1)]

    with pytest.raises(ValueError, match="internal-stuffer-cuts"):
        scheme(internal_stuffer_core=without)
