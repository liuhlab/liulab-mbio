"""Accepting a destination vector, and making one that is not compatible.

Every scheme here is laid out from the shipped enzyme definitions: the internal stuffer puts that
enzyme's two sites where its own cut offsets ask for them, so the overhangs a test asserts are read
back off the DNA rather than copied from whatever wrote it.
"""

from typing import Any

import pytest

from liulab_mbio.edits import rotate
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.library.scheme import Position, Scheme
from liulab_mbio.library.vector import destination_vector
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)

#: The jobs the test schemes give their enzymes. All four are free of sites in pUC19.
INTERNAL = "BbsI"
EXTERNAL = "PaqCI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"

#: Three entry overhangs and the cloning scar every part's 3' end leaves.
ENTRY = ("CTCC", "GGAG", "CCGA")
SCAR = "AGCG"

#: A stretch of pUC19 annotated by no feature, between the lac promoter and the origin.
GAP = (600, 700)


def pad(length: int) -> str:
    """`length` bases spelling no site any of the schemes' enzymes reads."""
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
    """A valid scheme of three positions, with any field replaced."""
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    positions = tuple(
        Position(
            f"p{index + 1}",
            internal_stuffer_prefix=ENTRY[(index + 1) % len(ENTRY)],
            external_stuffer_5=external_5(overhang, cutter, chopper),
            external_stuffer_3=external_3(cutter, chopper),
        )
        for index, overhang in enumerate(ENTRY)
    )
    fields: dict[str, Any] = {
        "positions": positions,
        "internal_enzyme": internal,
        "external_enzyme": EXTERNAL,
        "blunt_enzymes": (CORE_CHOPPER, FLANK_CHOPPER),
        "internal_stuffer_core": core(get_enzyme(internal), get_enzyme(CORE_CHOPPER)),
        "cloning_scar": SCAR,
        "barcode_length": 11,
    }
    fields.update(changes)
    return Scheme("test", **fields)


def carrier(made: Scheme, *, flank: int = 80) -> SequenceRecord:
    """A circular vector already carrying the scheme's internal stuffer."""
    return SequenceRecord(
        pad(flank) + made.internal_stuffer(-1) + pad(flank),
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
    assert excised.startswith(made.entry_overhang(0))
    assert vector.sequence[taken.stuffer.end : taken.stuffer.end + len(SCAR)] == made.scar_overhang


def test_a_stuffer_across_the_origin_is_found_where_it_lies():
    made = scheme()
    vector = rotate(carrier(made), 90)

    taken = destination_vector(vector, made)

    assert taken.stuffer.end > len(vector)
    assert taken.record is vector
    assert vector.extract(taken.stuffer).startswith(made.entry_overhang(0))


def test_puc19_is_made_compatible_at_a_named_span(puc19):
    made = scheme()
    stuffer = made.internal_stuffer(-1)

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
