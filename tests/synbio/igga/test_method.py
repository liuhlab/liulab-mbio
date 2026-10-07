"""What the method derives from its enzymes and stuffers, and what it refuses.

`IGGA` is checked when it is imported, so a test of it is a test of what the import already ran.
Every other method built here is laid out from the shipped enzyme definitions: a stuffer puts the
enzyme's site at the offset its own cut offsets ask for, so an overhang a test asserts is read
back off the DNA rather than copied from whatever wrote it.
"""

from typing import Any

import pytest

from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import reverse_complement
from liulab_mbio.translate import translate
from liulab_synbio.igga.method import IGGA, INTERFACE_OVERHANGS, Scheme

#: The four enzymes a built method gives each job, and the overhangs its stuffers spell.
INTERNAL = "BsaI"
EXTERNAL = "BbsI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"
ENTRY = "CTCC"
SCAR = "AGCG"

#: The spacer before the reverse cut. It falls on a codon boundary of the retained stuffer, where
#: `pad`'s T would spell a stop with the site that follows it.
LEAD = "C"


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

    Each cut is placed from `cutter`'s own offsets: one cuts back into the prefix, and one leaves
    `scar` as the stuffer's last bases. The filler length is what keeps the retained region a
    whole number of codons.
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


def scheme(**changes: Any) -> Scheme:
    """A valid method, with any field replaced."""
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    fields: dict[str, Any] = {
        "internal_enzyme": INTERNAL,
        "external_enzyme": EXTERNAL,
        "blunt_enzymes": (CORE_CHOPPER, FLANK_CHOPPER),
        "internal_stuffer_prefix": pad(2) + ENTRY,
        "internal_stuffer_core": core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER)),
        "external_stuffer_5": external_5(ENTRY, cutter, chopper),
        "external_stuffer_3": external_3(SCAR, cutter, chopper),
        "cloning_scar": SCAR,
    }
    fields.update(changes)
    return Scheme("test", **fields)


def test_igga_derives_what_a_build_reads_off_it():
    retained = IGGA.internal_stuffer

    assert IGGA.entry_overhang == "AGGA"
    assert IGGA.scar_overhang == IGGA.cloning_scar == "TTCC"
    assert IGGA.internal.site == "GAAGAC"
    assert IGGA.external.site == "GGTCTC"
    assert tuple(enzyme.end for enzyme in IGGA.blunt) == ("blunt", "blunt")
    assert [one.name for one in IGGA.reserved_enzymes] == ["BsmBI"]
    # The 34-mer the method page states, and the frame it leaves everything downstream in.
    assert len(retained) == 34
    assert IGGA.retained_length(11, 3) % 3 == 0
    assert len(IGGA.internal_stuffer_prefix) == len(IGGA.entry_overhang)
    assert "*" not in translate(retained[: len(retained) // 3 * 3])
    assert IGGA.source


def test_igga_reads_the_cargo_pair_off_its_own_dna():
    assert (IGGA.entry_overhang, IGGA.scar_overhang) == INTERFACE_OVERHANGS["cargo"]
    # Four overhangs over four pairs, each one base closing a codon plus a codon.
    assert len(INTERFACE_OVERHANGS) == 4
    assert {one for pair in INTERFACE_OVERHANGS.values() for one in pair} == {
        "TATG",
        "AGGA",
        "TTCC",
        "CTAA",
    }


def test_the_barcode_block_is_whole_codons_however_many_positions_there_are():
    assert [IGGA.retained_length(11, count) % 3 for count in (1, 2, 3, 8)] == [0, 0, 0, 0]


def test_a_prefix_not_ending_with_the_entry_overhang_is_refused():
    with pytest.raises(ValueError, match="internal-stuffer-prefix"):
        scheme(internal_stuffer_prefix=pad(2) + "TTTT")


def test_a_5_external_stuffer_the_external_enzyme_does_not_cut_is_refused():
    with pytest.raises(ValueError, match="external-stuffer-5"):
        scheme(external_stuffer_5=pad(20))


def test_a_3_external_stuffer_yielding_another_overhang_is_refused():
    odd = external_3("TTTT", get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER))

    with pytest.raises(ValueError, match="external-stuffer-3"):
        scheme(external_stuffer_3=odd)


def test_a_stuffer_leaving_the_product_out_of_frame_is_refused():
    with pytest.raises(ValueError, match="terminal-frame"):
        scheme(internal_stuffer_prefix=ENTRY)


def test_the_internal_enzyme_in_an_external_stuffer_is_refused():
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    trespass = get_enzyme(INTERNAL).site + external_5(ENTRY, cutter, chopper)

    with pytest.raises(ValueError, match="enzyme-regions"):
        scheme(external_stuffer_5=trespass)


def test_a_blunt_enzyme_chopping_nothing_is_refused():
    with pytest.raises(ValueError, match="enzyme-regions"):
        scheme(blunt_enzymes=(CORE_CHOPPER, FLANK_CHOPPER, "EcoRI"))


def test_a_reserved_enzyme_reading_a_stuffer_is_refused():
    with pytest.raises(ValueError, match="enzyme-regions"):
        scheme(reserved=(FLANK_CHOPPER,))


def test_an_internal_stuffer_leaving_no_scar_overhang_is_refused():
    whole = core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER))

    with pytest.raises(ValueError, match="internal-stuffer-cuts"):
        scheme(internal_stuffer_core=whole[: -(len(SCAR) + 1)])


def test_an_entry_overhang_that_is_the_cloning_scar_is_refused():
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)

    with pytest.raises(ValueError, match="internal-stuffer-cuts") as refused:
        scheme(
            internal_stuffer_prefix=pad(2) + SCAR,
            external_stuffer_5=external_5(SCAR, cutter, chopper),
        )

    assert "either way round" in str(refused.value)
