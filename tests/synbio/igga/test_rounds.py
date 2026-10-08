"""What one round joins, what the records annotate, and what a mismatch is refused with.

The method here is laid out from the shipped enzyme definitions: each stuffer puts an enzyme's
site at the offset that enzyme's own cut offsets ask for, so every overhang a test asserts is read
back off the DNA rather than copied from whatever wrote it. Position one's entry overhang is
pinned, which is what keeps the vector's own stuffer and the parts speaking the same standard.
"""

from pathlib import Path
from typing import Any

import pytest

from liulab_mbio.edits import rotate
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.io import read_record
from liulab_mbio.sequence import Segment, SequenceRecord, reverse_complement
from liulab_mbio.sites import digest, find_sites
from liulab_mbio.translate import translate
from liulab_synbio.igga.gate import check_product
from liulab_synbio.igga.method import Scheme
from liulab_synbio.igga.parts import barcode_rules, design_parts
from liulab_synbio.igga.project import Barcode, Build
from liulab_synbio.igga.rounds import (
    PRODUCT_FILE,
    ROUND_FILE,
    assemble_round,
    assemble_rounds,
    representative,
    write_records,
)
from liulab_synbio.igga.standard import design_standard

#: The jobs the test scheme gives its enzymes.
INTERNAL = "BsaI"
EXTERNAL = "BbsI"
CORE_CHOPPER = "SrfI"
FLANK_CHOPPER = "PmeI"

#: The positions one build fills, the overhang a part enters on, and the cloning scar its 3' end
#: leaves. Every position carries the method's own stuffers, so one entry overhang is stated.
POSITIONS = ("p1", "p2", "p3")
ENTRY = ("CTCC",)
SCAR = "AGCG"

HOST = "e-coli-k12"

#: How many bases name one part, as a build states.
BARCODE = 11

#: Three part lists, two members each, so a whole build stays small.
LISTS = (
    {"n1": "MKTAEK", "n2": "MKTCEK"},
    {"z1": "WQAFAK", "z2": "WQAYAK"},
    {"y1": "MKTHGK", "y2": "MKTWGK"},
)


def pad(length: int) -> str:
    """`length` bases spelling no site any of the method's enzymes reads."""
    return ("TA" * length)[:length]


def core(cutter: Enzyme, chopper: Enzyme) -> str:
    """The shared internal stuffer core, both of the cuts that open it laid from the offsets.

    One cut reaches back into whatever prefix precedes the core and the other leaves the cloning
    scar as the stuffer's last bases, so a stuffer carries both of the cuts a round needs.
    """
    lead = cutter.bottom_cut - len(cutter.site) - cutter.overhang_length
    reach = cutter.top_cut - len(cutter.site)
    return (
        pad(lead) + reverse_complement(cutter.site) + chopper.site + cutter.site + pad(reach) + SCAR
    )


def external_5(overhang: str, cutter: Enzyme, chopper: Enzyme) -> str:
    """A 5' external stuffer whose `cutter` cut leaves `overhang` as its last bases."""
    reach = pad(cutter.top_cut - len(cutter.site))
    return pad(3) + chopper.site + pad(3) + cutter.site + reach + overhang


def external_3(cutter: Enzyme, chopper: Enzyme) -> str:
    """A 3' external stuffer whose `cutter` cut leaves the cloning scar as its first bases."""
    reach = pad(cutter.top_cut - len(cutter.site))
    return SCAR + reach + reverse_complement(cutter.site) + pad(3) + chopper.site + pad(3)


def scheme(**changes: Any) -> Scheme:
    """A valid method, with any field replaced.

    The stuffer prefix is longer than an overhang, which is what buys the frame back over the
    stuffer the product keeps whole.
    """
    cutter, chopper = get_enzyme(EXTERNAL), get_enzyme(FLANK_CHOPPER)
    fields: dict[str, Any] = {
        "internal_enzyme": INTERNAL,
        "external_enzyme": EXTERNAL,
        "blunt_enzymes": (CORE_CHOPPER, FLANK_CHOPPER),
        "internal_stuffer_prefix": pad(1) + ENTRY[0],
        "internal_stuffer_core": core(get_enzyme(INTERNAL), get_enzyme(CORE_CHOPPER)),
        "external_stuffer_5": external_5(ENTRY[0], cutter, chopper),
        "external_stuffer_3": external_3(cutter, chopper),
        "cloning_scar": SCAR,
    }
    fields.update(changes)
    return Scheme("test", **fields)


def carrier(made: Scheme, *, flank: int = 80) -> SequenceRecord:
    """A circular vector already carrying the method's internal stuffer."""
    return SequenceRecord(
        pad(flank) + made.internal_stuffer + pad(flank), topology="circular", name="carrier"
    )


def build_of(made: Scheme) -> Build:
    """A build made under `made`, for the gate to read its own choices off."""
    return Build(
        "test",
        positions=POSITIONS,
        parts=Path("parts.fasta"),
        vector=Path("vector.gb"),
        host=HOST,
        oligo_length=350,
        batch_size=96,
        completeness=0.99,
        barcode=Barcode(length=BARCODE, min_distance=3),
        scheme=made,
    )


def built(made: Scheme):
    """The standard and the parts one build works from, position one's overhang pinned."""
    one = design_standard(made, POSITIONS, LISTS, pinned={"p1": ENTRY[0]})
    rules = barcode_rules(made, BARCODE)
    return one, design_parts(made, POSITIONS, LISTS, one, host=HOST, rules=rules)


def same_circle(one: SequenceRecord, other: SequenceRecord) -> bool:
    """Whether two circular records are the same molecule, read from wherever each begins."""
    return len(one) == len(other) and other.sequence in one.sequence + one.sequence


@pytest.fixture(scope="module")
def made():
    return scheme()


@pytest.fixture(scope="module")
def build(made):
    return built(made)


@pytest.fixture(scope="module")
def standard(build):
    return build[0]


@pytest.fixture(scope="module")
def parts(build, made):
    return representative(build[1], POSITIONS)


@pytest.fixture(scope="module")
def vector(made):
    return carrier(made)


@pytest.fixture(scope="module")
def rounds(vector, parts, made):
    return assemble_rounds(vector, parts, made, POSITIONS, name="library")


def test_a_round_replaces_the_destination_s_stuffer_with_the_part(vector, parts, made):
    one = assemble_round(vector, parts[0], made)

    assert one.product.topology == "circular"
    assert len(one.product) == len(vector) - one.excised.length + one.released.length
    assert one.product.extract(one.coding) == parts[0].coding_sequence
    assert one.product.extract(one.barcode) == parts[0].barcode
    assert one.product.extract(one.entry) == one.entry_overhang == ENTRY[0]
    assert one.product.extract(one.scar) == one.scar_overhang == SCAR


def test_the_product_is_the_two_digest_fragments_joined(vector, parts, made):
    one = assemble_round(vector, parts[0], made)
    keep = next(
        piece
        for piece in digest(vector, made.internal)
        if (piece.left_overhang, piece.right_overhang) == (one.scar_overhang, one.entry_overhang)
    )
    released = SequenceRecord(parts[0].sequence).extract(
        Segment(one.released.start, one.released.end)
    )

    joined = released + vector.extract(Segment(keep.start, keep.end))

    assert same_circle(one.product, SequenceRecord(joined, topology="circular"))


def test_overhangs_that_do_not_match_are_refused_naming_both(vector, parts, made, standard):
    # The vector opens on position one's overhangs; the second round's part enters on another.
    with pytest.raises(ValueError, match="ligates only where both overhangs match") as caught:
        assemble_round(vector, parts[1], made, number=2)

    message = str(caught.value)
    assert standard.entry_overhangs[1] in message
    assert ENTRY[0] in message
    assert parts[1].name in message


def test_each_round_appends_its_barcode_to_the_front_of_the_block(rounds, parts, made):
    for number, one in enumerate(rounds, start=1):
        block = one.product.extract(one.block)

        assert block.startswith(one.product.extract(one.barcode))
        assert block == SCAR.join(part.barcode for part in reversed(parts[:number]))
        assert one.block.end - one.block.start == len(block)
    # Three barcodes and two scars: the vector's own barcode is not one of them.
    block = made.barcode_block_length(BARCODE, len(POSITIONS))
    assert len(rounds[-1].product.extract(rounds[-1].block)) == block == 41


def test_the_retained_region_is_frame_correct_and_free_of_in_phase_stops(rounds, made):
    final = rounds[-1]

    bases = final.product.extract(final.retained)

    assert len(bases) == made.retained_length(BARCODE, len(POSITIONS))
    assert len(bases) % 3 == 0
    assert "*" not in translate(bases)


def test_a_stop_in_the_retained_region_is_the_gate_s_to_report_not_a_round_s():
    # A terminal prefix spelling TGA where the product reads through it, which nothing excises.
    stopping = scheme(internal_stuffer_prefix="TGAT" + ENTRY[0])
    parts = built(stopping)[1]

    made = assemble_rounds(carrier(stopping), representative(parts, POSITIONS), stopping, POSITIONS)

    final = made[-1]
    assert translate(final.product.extract(final.retained)).startswith("*")
    judged = check_product(
        final.product,
        build=build_of(stopping),
        barcodes={one.position: [one.barcode] for one in representative(parts, POSITIONS)},
    )
    stops = next(one for one in judged if one.name == "terminal stop")
    assert stops.status == "fail"
    assert stops.check.value >= 1


def test_every_product_keeps_the_internal_sites_and_loses_the_external_ones(rounds, made):
    for one in rounds:
        assert find_sites(one.product, made.internal)
        assert find_sites(one.product, made.external) == ()


def test_the_last_product_opens_on_position_one_s_overhang_again(rounds, made, standard):
    final = rounds[-1]

    opened = [
        piece
        for piece in digest(final.product, made.internal)
        if (piece.start, piece.end) == (final.stuffer.start, final.stuffer.end)
    ]

    assert len(opened) == 1
    assert opened[0].left_overhang == standard.entry_overhangs[0]
    assert opened[0].right_overhang == SCAR


def test_the_representative_annotates_every_part_stuffer_barcode_and_junction(rounds, parts):
    final = rounds[-1]

    named = [feature.name for feature in final.product.features]

    for part in parts:
        assert part.name in named
        assert f"{part.name} barcode" in named
        assert f"{part.position} entry junction" in named
    assert named.count("cloning scar") == len(parts)
    # Only the stuffer the product keeps is left: each round excised the one before it.
    assert [one for one in named if one.endswith("internal stuffer")] == [
        f"{parts[-1].position} internal stuffer"
    ]
    for feature in final.product.features:
        assert final.product.extract(feature)


def test_the_records_are_written_and_read_back_with_their_features(rounds, tmp_path):
    written = write_records(rounds, tmp_path / "records")

    assert [path.name for path in written] == [
        ROUND_FILE.format(number=1),
        ROUND_FILE.format(number=2),
        PRODUCT_FILE,
    ]
    back = read_record(written[-1])
    final = rounds[-1].product
    assert back.sequence == final.sequence
    assert back.topology == "circular"
    assert [feature.name for feature in back.features] == [
        feature.name for feature in final.features
    ]
    for feature, before in zip(back.features, final.features, strict=True):
        assert back.extract(feature) == final.extract(before)


def test_a_stuffer_across_the_origin_is_opened_where_it_lies(vector, parts, made, rounds):
    turned = rotate(vector, 90)
    assert turned.sequence != vector.sequence

    made_rounds = assemble_rounds(turned, parts, made, POSITIONS)

    assert same_circle(made_rounds[-1].product, rounds[-1].product)
    final = made_rounds[-1]
    assert final.product.extract(final.block) == rounds[-1].product.extract(rounds[-1].block)


def test_parts_out_of_the_build_s_order_are_refused(vector, parts, made):
    with pytest.raises(ValueError, match="the rounds run in the build's own order"):
        assemble_rounds(vector, (parts[1], parts[0], parts[2]), made, POSITIONS)


def test_one_part_a_position_is_required(vector, parts, made):
    with pytest.raises(ValueError, match="one a round"):
        assemble_rounds(vector, parts[:2], made, POSITIONS)


def test_a_linear_destination_is_refused(vector, parts, made):
    linear = SequenceRecord(vector.sequence, name="linear")

    with pytest.raises(ValueError, match="circular destination"):
        assemble_round(linear, parts[0], made)


def test_a_representative_needs_a_part_for_every_position(build, made):
    # Both of these fill position one, so nothing is left to fill position two.
    with pytest.raises(ValueError, match="no part fills position 'p2'"):
        representative(build[1][:2], POSITIONS)


def test_no_round_is_nothing_to_write(tmp_path):
    with pytest.raises(ValueError, match="at least one round"):
        write_records((), tmp_path)
