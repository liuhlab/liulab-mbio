"""Golden Gate assembly simulated on the fixtures.

pUC19 (circular) and the GFP CDS (linear) are read through `liulab_mbio.io`. BsaI reads a site
in both, so the enzyme here is BbsI or PaqCI. Nothing about the fixtures is hard-coded into the
package: the spans and overhangs below are read off the records and pinned here.
"""

import dataclasses

import pytest

from liulab_mbio import bench, edits
from liulab_mbio.goldengate.assembly import (
    JUNCTION_COLOR,
    amplify,
    assemble,
    dam_sites,
    open_vector,
)
from liulab_mbio.sequence import (
    BindingSite,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)
from liulab_mbio.sites import find_sites, has_site, insert_site
from liulab_mbio.snapgene import read_dna, write_dna

#: How many bases `primer_tail` puts 5' of the recognition site.
SPACER = 6

#: Addgene's 23-mer M13/pUC pair, which #13 pins its own expected bands against.
M13_FORWARD = Primer("M13/pUC Forward", "CCCAGTCACGACGTTGTAAAACG")
M13_REVERSE = Primer("M13/pUC Reverse", "AGCGGATAACAATTTCACACAGG")


@pytest.fixture(scope="module")
def mcs(puc19: SequenceRecord) -> tuple[int, int]:
    """The span the assembly replaces, read off the fixture's own MCS feature."""
    segment = next(one for one in puc19.features if one.name == "MCS").segments[0]
    return segment.start, segment.end


@pytest.fixture(scope="module")
def overhangs(puc19: SequenceRecord, gfp: SequenceRecord, mcs: tuple[int, int]) -> tuple[str, str]:
    """A scarless pair: the insert's own first bases, and the vector's own first bases past the MCS."""
    return gfp.sequence[:4], puc19.sequence[mcs[1] : mcs[1] + 4]


@pytest.fixture(scope="module")
def backbone(puc19: SequenceRecord, mcs: tuple[int, int], overhangs: tuple[str, str]):
    return open_vector(puc19, "BbsI", *mcs, overhangs=overhangs, name="pUC19 backbone")


@pytest.fixture(scope="module")
def insert(gfp: SequenceRecord, overhangs: tuple[str, str]):
    return amplify(
        gfp,
        "BbsI",
        0,
        len(gfp),
        left_overhang=overhangs[0],
        right_overhang=overhangs[1],
        name="GFP",
    )


def test_the_fixtures_are_what_the_design_assumes(puc19, gfp, mcs, overhangs):
    assert (len(puc19), puc19.topology) == (2686, "circular")
    assert (len(gfp), gfp.topology) == (717, "linear")
    assert mcs == (395, 452)
    assert overhangs == ("ATGA", "GGCG")
    assert not has_site(puc19, "BbsI")
    assert not has_site(gfp, "BbsI")


def test_vector_primers_face_outward_from_the_span(backbone, puc19, mcs, overhangs):
    forward = backbone.forward.binding_sites[0]
    reverse = backbone.reverse.binding_sites[0]
    # The forward primer reads away from the span, starting past the overhang the backbone keeps.
    assert (forward.start, forward.strand) == (mcs[1] + len(overhangs[1]), Strand.FORWARD)
    # The reverse primer reads the other way, ending where the span begins.
    assert (reverse.end, reverse.strand) == (mcs[0], Strand.REVERSE)
    assert backbone.span == (mcs[1], mcs[0] + len(puc19))


def test_each_vector_primer_carries_one_site_and_its_overhang(backbone, overhangs):
    # The forward primer copies the top strand and the reverse the bottom, so each tail ends
    # with the overhang as its own strand spells it.
    for primer, overhang in (
        (backbone.forward, overhangs[1]),
        (backbone.reverse, reverse_complement(overhangs[0])),
    ):
        tail = primer.sequence[: len(primer.sequence) - _annealed(primer)]
        sites = find_sites(SequenceRecord(tail), "BbsI")
        assert len(sites) == 1
        assert len(tail) == SPACER + sites[0].enzyme.bottom_cut
        assert tail.endswith(overhang)


def test_the_backbone_amplicon_is_the_vector_without_the_span(backbone, puc19, mcs, overhangs):
    kept = len(puc19) - (mcs[1] - mcs[0]) - len(overhangs[1])
    assert len(backbone.amplicon) == 2 * (SPACER + 12) + kept
    assert len(backbone.amplicon) == 2661
    assert backbone.amplicon.topology == "linear"
    # What #8 says the pair makes is what this simulation built.
    assert backbone.report.amplicon_length == len(backbone.amplicon)
    assert backbone.report["products"].value == 1


def test_the_insert_amplicon_carries_the_whole_coding_sequence(insert, gfp):
    assert len(insert.amplicon) == 2 * (SPACER + 12) + len(gfp) - 4
    assert insert.amplicon.sequence.count(gfp.sequence[4:]) == 1
    assert insert.report.amplicon_length == len(insert.amplicon)


def test_dpni_removes_a_plasmid_template_and_nothing_else(backbone, insert, puc19, gfp):
    assert backbone.dpni is True
    # The linear fixture carries Dam sites too; being no plasmid is what spares it.
    assert insert.dpni is False
    assert dam_sites(puc19) == 15
    assert dam_sites(gfp) == 2


def test_dam_sites_are_counted_across_the_origin():
    assert dam_sites(SequenceRecord("TCAAAAAAGA", topology="circular")) == 1
    assert dam_sites(SequenceRecord("TCAAAAAAGA")) == 0


def test_open_vector_refuses_a_linear_vector(gfp, overhangs):
    with pytest.raises(ValueError, match="circular"):
        open_vector(gfp, "BbsI", 10, 20, overhangs=overhangs)


def test_amplify_refuses_a_span_the_overhang_swallows(gfp, overhangs):
    with pytest.raises(ValueError, match="shorter"):
        amplify(gfp, "BbsI", 0, 3, left_overhang=overhangs[0], right_overhang=overhangs[1])


def test_features_are_carried_into_the_amplicon(backbone, insert, puc19, gfp):
    carried = {one.name: one for one in backbone.amplicon.features}
    # The MCS was removed with the span, and everything outside it survives.
    assert "MCS" not in carried
    assert _spans(carried["AmpR"]) == [(1187, 1979), (1979, 2048)]
    assert carried["AmpR"].strand == Strand.REVERSE
    # The disrupted CDS meets the amplified span twice, once either side of the removed MCS.
    assert _spans(carried["lacZα"]) == [(14, 31), (2393, 2643)]  # noqa: RUF001
    assert [one.name for one in insert.amplicon.features] == ["GFP"]


def test_designed_primers_are_annotated_on_the_amplicon(backbone):
    annotated = {one.name: one for one in backbone.amplicon.primers}
    assert set(annotated) == {"pUC19 backbone forward", "pUC19 backbone reverse"}
    for primer in annotated.values():
        site = primer.binding_sites[0]
        annealed = backbone.amplicon.extract(Segment(site.start, site.end))
        assert primer.sequence.endswith(
            annealed if site.strand == Strand.FORWARD else reverse_complement(annealed)
        )


def test_a_template_primer_is_kept_only_where_its_whole_site_survives():
    kept = Primer("kept", "CCCCGGGGAAAA", binding_sites=(BindingSite(4, 16, Strand.FORWARD),))
    cut = Primer("cut", "TTTTTTTT", binding_sites=(BindingSite(0, 8, Strand.FORWARD),))
    record = SequenceRecord(
        "TTTTAAAACCCCGGGGAAAATTTTCCCCGGGG", topology="circular", primers=(kept, cut)
    )
    part = amplify(record, "BbsI", 4, 32, left_overhang="AAAA", right_overhang="CCCC")
    # "cut" annealed across the edge of the span, so it has nowhere left to sit.
    assert [one.name for one in part.amplicon.primers] == ["kept", "forward", "reverse"]


@pytest.fixture(scope="module")
def assembly(backbone, insert):
    return assemble((backbone, insert), "BbsI", name="pUC19-GFP")


def test_the_product_replaces_the_span_with_the_insert(assembly, puc19, gfp, mcs):
    product = assembly.product
    assert product.topology == "circular"
    assert product.sequence == puc19.sequence[: mcs[0]] + gfp.sequence + puc19.sequence[mcs[1] :]
    assert len(product) == 3346
    assert product.name == "pUC19-GFP"


def test_the_product_keeps_the_vector_origin(assembly, puc19, mcs):
    # Rotating back to the backbone's own first base keeps the vector's coordinates readable
    # and keeps a junction off base zero, which is what #13's validation reads across.
    assert assembly.product.sequence[: mcs[0]] == puc19.sequence[: mcs[0]]


def test_junctions_are_the_overhangs_standing_in_the_product(assembly, overhangs):
    assert assembly.junction_positions == (395, 1112)
    at = {one.start: one for one in assembly.junctions}
    assert (at[395].overhang, at[1112].overhang) == overhangs
    for junction in assembly.junctions:
        assert assembly.product.extract(junction.span) == junction.overhang
    assert (at[395].before, at[395].after) == ("pUC19 backbone", "GFP")
    assert (at[1112].before, at[1112].after) == ("GFP", "pUC19 backbone")


def test_features_carry_over_to_their_new_coordinates(assembly):
    carried = {one.name: one for one in assembly.product.features}
    assert _spans(carried["GFP"]) == [(395, 1112)]
    assert _spans(carried["lacZα"]) == [(145, 395), (1112, 1129)]  # noqa: RUF001
    assert _spans(carried["AmpR"]) == [(2285, 3077), (3077, 3146)]
    assert _spans(carried["ori"]) == [(1526, 2115)]
    # The MCS went with the span the insert replaced.
    assert "MCS" not in carried


def test_every_feature_of_the_product_carries_a_colour(assembly):
    # SnapGene writes its own default grey for a feature that has none.
    assert all(one.color for one in assembly.product.features)


def test_each_junction_is_annotated(assembly, overhangs):
    drawn = [one for one in assembly.product.features if one.color == JUNCTION_COLOR]
    assert [one.name for one in drawn] == [f"{one} junction" for one in overhangs]
    assert _spans(drawn[0]) == [(395, 399)]


def test_the_designed_primers_are_annotated_where_they_anneal(assembly):
    annotated = {one.name: one for one in assembly.product.primers}
    assert {
        "pUC19 backbone forward",
        "pUC19 backbone reverse",
        "GFP forward",
        "GFP reverse",
    } <= set(annotated)
    for primer in annotated.values():
        site = primer.binding_sites[0]
        annealed = assembly.product.extract(Segment(site.start, site.end))
        assert primer.sequence.endswith(
            annealed if site.strand == Strand.FORWARD else reverse_complement(annealed)
        )


def test_the_product_writes_and_reads_back_unchanged(assembly, tmp_path):
    path = tmp_path / "product.dna"
    write_dna(assembly.product, path)
    assert read_dna(path) == assembly.product


def test_another_enzyme_reaching_further_leaves_the_same_product(
    assembly, puc19, gfp, mcs, overhangs
):
    # PaqCI reads a longer site and cuts further from it, so the tails differ and the product
    # must not.
    parts = (
        open_vector(puc19, "PaqCI", *mcs, overhangs=overhangs, name="pUC19 backbone"),
        amplify(
            gfp,
            "PaqCI",
            0,
            len(gfp),
            left_overhang=overhangs[0],
            right_overhang=overhangs[1],
            name="GFP",
        ),
    )
    assert len(parts[0].amplicon) != len(assembly.parts[0].amplicon)
    other = assemble(parts, "PaqCI", name="pUC19-GFP")
    assert other.product.sequence == assembly.product.sequence
    assert other.junction_positions == assembly.junction_positions


def test_assemble_refuses_parts_that_do_not_close_the_circle(backbone, gfp, overhangs):
    stranger = amplify(
        gfp, "BbsI", 0, len(gfp), left_overhang=overhangs[0], right_overhang="TTAG", name="GFP"
    )
    with pytest.raises(ValueError, match="TTAG"):
        assemble((backbone, stranger), "BbsI")


def test_assemble_refuses_one_part(backbone):
    with pytest.raises(ValueError, match="two"):
        assemble((backbone,), "BbsI")


def test_the_bands_bench_expects_are_the_ones_this_product_gives(assembly, puc19):
    check = bench.colony_pcr_check(
        assembly.product,
        assembly.junction_positions,
        vector=puc19,
        primers=(M13_FORWARD, M13_REVERSE),
    )
    bands = {one.name: one.bands_bp for one in check.clones}
    assert bands["Correct clone"] == (797,)
    assert bands["Empty vector"] == (137,)


def test_a_part_puts_exactly_its_own_span_into_the_product(backbone, insert, gfp):
    # The insert is scarless, so its overhang is its own first bases and nothing is added.
    assert insert.bases == gfp.sequence
    assert len(backbone.bases) == backbone.fragment_length


def test_the_product_passes_every_check(assembly):
    assert assembly.status == "pass"
    assert [one.name for one in assembly.checks] == [
        "sites",
        "pUC19 backbone",
        "GFP",
        "junctions",
    ]
    assert assembly["sites"].value == 0
    assert assembly["GFP"].value == 1
    assert assembly["pUC19 backbone"].value == 1
    assert assembly["junctions"].value == 2


def test_a_site_left_for_the_enzyme_fails_the_check(assembly):
    spoiled, _ = insert_site(assembly.product, "BbsI", 2000)
    checked = dataclasses.replace(assembly, product=spoiled)
    assert checked["sites"].value == 1
    assert checked["sites"].status == "fail"
    assert checked.status == "fail"


def test_an_insert_found_twice_fails_its_check(assembly, gfp):
    doubled, _ = edits.insert(assembly.product, 2000, gfp.sequence)
    checked = dataclasses.replace(assembly, product=doubled)
    assert checked["GFP"].value == 2
    assert checked["GFP"].status == "fail"


def test_a_junction_that_lost_its_bases_fails_the_check(assembly):
    changed, _ = edits.replace(assembly.product, 395, 399, "TTTT")
    checked = dataclasses.replace(assembly, product=changed)
    assert checked["junctions"].value == 1
    assert checked["junctions"].status == "fail"


def test_an_unknown_check_is_a_key_error(assembly):
    with pytest.raises(KeyError):
        assembly["nonsense"]


def _annealed(primer: Primer) -> int:
    return primer.binding_sites[0].end - primer.binding_sites[0].start


def _spans(feature) -> list[tuple[int, int]]:
    return [(one.start, one.end) for one in feature.segments]
