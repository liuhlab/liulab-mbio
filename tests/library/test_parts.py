"""What a part's synthesis sequence holds, what a digest leaves, and what the sheets say.

The scheme is the worked example a user supplies, so every overhang a test asserts is read back
off the DNA the parts are built from rather than stated here.
"""

from pathlib import Path

import pytest

from liulab_mbio.barcodes import BarcodeRules, check_barcodes
from liulab_mbio.library.parts import (
    BARCODE_COLUMNS,
    CHANGE_COLUMNS,
    SHEET_COLUMNS,
    barcode_phase,
    barcode_rules,
    barcode_table,
    change_table,
    design_parts,
    synthesis_sheet,
)
from liulab_mbio.library.scheme import read_scheme
from liulab_mbio.library.standard import design_standard, junction_residues
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import digest, find_sites
from liulab_mbio.translate import reverse_translate, translate

#: The paper's scheme, as a user supplies one.
EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "protein-library" / "scheme.json"

HOST = "e-coli-k12"

#: Three part lists whose junction residues disagree, so the standard has to charge somebody.
LISTS = (
    {"ATF2": "MKTQED", "ATF4": "MKTQEK"},
    {"b1": "WQAAAK", "b2": "WQAAAR"},
    {"c1": "MKTQGH", "c2": "MKTQGY"},
)

#: Three part lists agreeing at every junction, so the standard need charge nobody.
AGREED = (
    {"n1": "MKTAEK", "n2": "MKTCEK"},
    {"z1": "WQAFAK", "z2": "WQAYAK"},
    {"y1": "MKTHGK", "y2": "MKTWGK"},
)


@pytest.fixture(scope="module")
def scheme():
    return read_scheme(EXAMPLE)


@pytest.fixture(scope="module")
def standard(scheme):
    return design_standard(scheme, LISTS)


@pytest.fixture(scope="module")
def parts(scheme, standard):
    return design_parts(scheme, LISTS, standard, host=HOST)


@pytest.fixture(scope="module")
def agreed(scheme):
    """The standard over the part lists that agree at every junction, which charges nobody."""
    return design_standard(scheme, AGREED)


def enzymes(scheme):
    """Every enzyme the scheme names."""
    return (scheme.internal, scheme.external, *scheme.blunt)


def charged(standard, part, end):
    """What the standard charges one end of one part, or None where it charges nothing."""
    return next(
        (
            one
            for one in standard.termini
            if one.part == part.name and one.position == part.position and one.end == end
        ),
        None,
    )


def stuffer_span(scheme, part):
    """Where the internal stuffer lies in a part's block, read off the block's own shape."""
    ext3 = len(scheme.positions[part.index].external_stuffer_3)
    return part.coding.end, part.length - ext3 - len(part.barcode)


def test_a_part_is_its_stuffers_its_coding_bases_and_its_barcode_in_order(scheme, standard, parts):
    assert len(parts) == sum(len(one) for one in LISTS)
    for part in parts:
        position = scheme.positions[part.index]
        assert part.sequence.startswith(
            position.external_stuffer_5[: -len(standard.entry_overhangs[part.index])]
        )
        assert part.sequence.endswith(position.external_stuffer_3[len(scheme.cloning_scar) :])
        start, end = stuffer_span(scheme, part)
        assert part.sequence[start:end].endswith(scheme.internal_stuffer_core)
        assert part.sequence[end : end + len(part.barcode)] == part.barcode
        assert part.length == len(part.sequence)


def test_a_part_enters_on_the_standard_s_overhang_and_not_the_scheme_s(scheme, standard, parts):
    # The scheme's own stuffers spell CTCC, GGAG and CCGA; the standard chooses for the proteins.
    assert standard.entry_overhangs != scheme.entry_overhangs
    for part in parts:
        entry = standard.entry_overhangs[part.index]
        assert part.sequence[part.coding.start - len(entry) : part.coding.start] == entry


def test_every_part_spells_only_the_sites_its_stuffers_do(scheme, parts):
    for part in parts:
        found = find_sites(SequenceRecord(part.sequence), enzymes(scheme))
        assert found
        for site in found:
            inside_coding = site.start < part.coding.end and site.end > part.coding.start
            assert not inside_coding, (part.name, site.enzyme.name, site.start)
            _, barcode_at = stuffer_span(scheme, part)
            touches_barcode = site.start < barcode_at + len(part.barcode) and site.end > barcode_at
            assert not touches_barcode, (part.name, site.enzyme.name, site.start)


def test_the_external_digest_releases_the_part_on_the_standard_s_two_overhangs(
    scheme, standard, parts
):
    for part in parts:
        pieces = digest(SequenceRecord(part.sequence), scheme.external)
        released = [
            piece
            for piece in pieces
            if piece.left_overhang == standard.entry_overhangs[part.index]
            and piece.right_overhang == standard.scar_overhang
        ]
        assert len(released) == 1, (
            part.name,
            [(p.left_overhang, p.right_overhang) for p in pieces],
        )


def test_the_internal_digest_excises_the_stuffer_on_the_next_round_s_overhangs(
    scheme, standard, parts
):
    for part in parts:
        following = standard.entry_overhangs[(part.index + 1) % scheme.position_count]
        pieces = digest(SequenceRecord(part.sequence), scheme.internal)
        excised = [
            piece
            for piece in pieces
            if piece.left_overhang == following and piece.right_overhang == standard.scar_overhang
        ]
        assert len(excised) == 1, (part.name, [(p.left_overhang, p.right_overhang) for p in pieces])
        start, end = stuffer_span(scheme, part)
        assert start <= excised[0].start < excised[0].end <= end


def test_the_product_reads_every_part_s_protein_in_one_frame(scheme, standard, parts):
    """Join one part a position the way the rounds do, and read the whole insert as codons.

    This is the frame the scheme charges the proteins for: each part's own bases plus the
    overhangs either side spell its amino acids exactly, and the retained stuffer carries the
    frame to a barcode block that begins `barcode_phase` bases into a codon and holds no stop.
    """
    chosen = [next(one for one in parts if one.index == at) for at in range(scheme.position_count)]
    donated, codons = junction_residues(len(standard.entry_overhangs[0]))
    reading = "".join(standard.entry_overhangs[one.index] + one.coding_sequence for one in chosen)
    last = chosen[-1]
    reading += last.sequence[slice(*stuffer_span(scheme, last))] + last.barcode
    for one in reversed(chosen[:-1]):
        reading += scheme.cloning_scar + one.barcode

    trimmed = reading[(3 - donated) % 3 :]

    assert len(trimmed) % 3 == 0
    joined = "".join(one.protein for one in chosen)
    assert translate(trimmed)[codons : codons + len(joined)] == joined
    # Scoped to the block: what the retained stuffer itself spells is the scheme's own data, and
    # this worked example's synthetic stuffer does spell a stop.
    block = len(trimmed) - scheme.barcode_block_length - barcode_phase(scheme)
    assert block % 3 == 0
    assert "*" not in translate(trimmed[block:])


def test_the_barcode_block_begins_one_base_into_a_codon(scheme, parts):
    assert barcode_phase(scheme) == 1
    # The scheme's own two frame invariants make this the same number as the scar's length.
    assert barcode_phase(scheme) == len(scheme.cloning_scar) % 3
    assert barcode_rules(scheme).phase == 1
    assert check_barcodes([one.barcode for one in parts], barcode_rules(scheme)) == ()


def test_the_sheet_holds_one_row_a_part(parts):
    sheet = synthesis_sheet(parts)

    rows = sheet.splitlines()
    assert rows[0].split("\t") == list(SHEET_COLUMNS)
    assert len(rows) == len(parts) + 1
    for part, row in zip(parts, rows[1:], strict=True):
        name, sequence, length, position, barcode = row.split("\t")
        assert (name, sequence, position, barcode) == (
            part.name,
            part.sequence,
            part.position,
            part.barcode,
        )
        assert int(length) == part.length
    assert sheet.endswith("\n")


def test_the_barcode_table_names_every_part_and_where_it_reads(scheme, parts):
    rows = barcode_table(parts, scheme).splitlines()

    assert rows[0].split("\t") == list(BARCODE_COLUMNS)
    assert len(rows) == len(parts) + 1
    for part, row in zip(parts, rows[1:], strict=True):
        name, position, number, slot, barcode = row.split("\t")
        assert (name, position, barcode) == (part.name, part.position, part.barcode)
        # The block reads newest first, so the last round's part is slot one.
        assert int(number) == part.index + 1
        assert int(slot) == scheme.position_count - part.index


def test_coding_bases_given_are_kept_but_for_a_codon_a_site_moves(scheme, standard, parts):
    # Coded for another host, so every codon differs from the one this host would have written.
    supplied = {
        name: reverse_translate(protein, host="human")
        for one in LISTS
        for name, protein in one.items()
    }
    coding = tuple({name: supplied[name] for name in one} for one in LISTS)

    made = design_parts(scheme, LISTS, standard, host=HOST, coding=coding)

    assert [one.coding_sequence for one in made] != [one.coding_sequence for one in parts]
    for part in made:
        head = charged(standard, part, "5'")
        tail = charged(standard, part, "3'")
        first = len(head.wild_type) if head is not None else 0
        whole = len(part.protein) - first - (len(tail.wild_type) if tail is not None else 0)
        body = supplied[part.name][3 * first : 3 * (first + whole)]
        assert body
        # Every codon is the one that was handed over, except the ones the part reports moving
        # to take a forbidden site out. Nothing else was written again.
        expected = [body[at : at + 3] for at in range(0, len(body), 3)]
        for change in part.changes:
            assert expected[change.codon_index] == change.old_codon, part.name
            expected[change.codon_index] = change.new_codon
        assert part.coding_sequence.startswith("".join(expected)), part.name


def test_the_change_table_names_every_part_whose_residues_moved(scheme, standard):
    table = change_table(standard)

    rows = table.splitlines()
    assert rows[0].split("\t") == list(CHANGE_COLUMNS)
    assert len(rows) == len(standard.changes) + 1
    assert standard.changes
    for one, row in zip(standard.changes, rows[1:], strict=True):
        assert row.split("\t") == [
            one.part,
            one.position,
            one.end,
            one.wild_type,
            one.synthesised,
        ]


def test_the_change_table_is_empty_where_the_standard_moved_nothing(agreed):
    assert agreed.cost == 0
    assert change_table(agreed) == "\t".join(CHANGE_COLUMNS) + "\n"


def test_the_same_inputs_write_the_same_bytes(scheme, standard):
    once = design_parts(scheme, LISTS, standard, host=HOST)
    twice = design_parts(scheme, LISTS, standard, host=HOST)

    assert synthesis_sheet(once) == synthesis_sheet(twice)
    assert [one.sequence for one in once] == [one.sequence for one in twice]


def test_another_seed_draws_other_barcodes(scheme, standard, parts):
    other = design_parts(scheme, LISTS, standard, host=HOST, seed=7)

    assert [one.barcode for one in other] != [one.barcode for one in parts]


def test_part_lists_must_match_the_scheme(scheme, standard):
    with pytest.raises(ValueError, match="part list"):
        design_parts(scheme, LISTS[:2], standard, host=HOST)


def test_a_standard_charging_other_part_lists_is_refused(scheme, agreed):
    with pytest.raises(ValueError, match="design it over these same part lists"):
        design_parts(scheme, LISTS, agreed, host=HOST)


def test_a_part_spelling_a_site_the_scheme_does_not_expect_is_refused(scheme, standard):
    """The rules keep a scheme's sites out of a barcode; drop that rule and this seed lets one in.

    The refusal names the part and the site, which is what someone has to act on.
    """
    loose = BarcodeRules(
        scheme.barcode_length,
        scar=scheme.cloning_scar,
        phase=barcode_phase(scheme),
        forbidden=(),
    )

    with pytest.raises(ValueError, match=r"part 'ATF4'") as caught:
        design_parts(scheme, LISTS, standard, host=HOST, rules=loose, seed=62)

    assert "BsaI" in str(caught.value)
    assert "where the scheme expects none" in str(caught.value)


def test_an_unshipped_host_is_refused(scheme, standard):
    with pytest.raises(KeyError, match="codon usage table"):
        design_parts(scheme, LISTS, standard, host="nowhere")
