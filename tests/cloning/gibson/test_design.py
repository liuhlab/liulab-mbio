"""The overlap rule and the two oligo routes, as units, against the numbers the note states.

These are the fast tests, and the ones carrying the values of `docs/research/gibson-assembly.md`:
the Wallace rule of §4, the NEBuilder bands of §3, the tiers of §7 and §9, the amounts and ratios
of §8, the stitching shape of §14 and the bridging homology of §15.
"""

from itertools import pairwise

import pytest

from liulab_mbio.cloning.gibson.assembly import Part, stitch
from liulab_mbio.cloning.gibson.bench import (
    ASSEMBLY_PRODUCTS,
    GIBSON_MASTER_MIX,
    IN_FUSION,
    NEBUILDER_HIFI,
    AssemblyProduct,
    OverlapRule,
    Tier,
    assembly_amounts,
    assembly_dna_check,
    assembly_product,
    assembly_program,
    assembly_reaction,
    fragment_check,
)
from liulab_mbio.cloning.gibson.design import (
    BRIDGE_HOMOLOGY_BP,
    STITCH_OLIGO_BASES,
    STITCH_OLIGOS,
    STITCH_OVERLAP_BP,
    STITCH_WINDOW_BP,
    bridging_oligo,
    overlap_after,
    overlap_before,
    overlap_checks,
    requires_oligos,
    stitch_checks,
    stitch_limit_bp,
    stitch_oligos,
    wallace_tm,
)
from liulab_mbio.sequence import SequenceRecord, Strand, reverse_complement

#: The band and floor NEBuilder HiFi documents for two or three fragments.
NEB_RULE = NEBUILDER_HIFI.tiers[0].overlap


def test_the_wallace_rule_is_two_degrees_a_weak_pair_and_four_a_strong_one():
    assert wallace_tm("AT") == 4.0
    assert wallace_tm("GC") == 8.0
    # NEB's own worked figure: a 15 bp overlap reaches 48 °C only at nine G or C bases.
    assert wallace_tm("GGGGGGGGGAAAAAA") == 48.0
    assert wallace_tm("GGGGGGGGAAAAAAA") == 46.0


def test_the_overlap_is_the_shortest_in_the_band_that_reaches_the_tm_floor():
    # Sixteen bases of this reach 48 °C and fifteen do not, so sixteen is what NEB's rule takes.
    record = SequenceRecord("A" * 40 + "TAAAACGACGGCCAGT" + "C" * 40)
    floor = NEB_RULE.tm_floor
    assert floor is not None
    assert floor == 48.0
    taken = overlap_before(record, 56, NEB_RULE)
    assert taken == "TAAAACGACGGCCAGT"
    assert len(taken) == 16
    assert wallace_tm(taken) >= floor
    assert wallace_tm(taken[1:]) < floor


def test_an_overlap_that_cannot_reach_the_floor_takes_the_longest_the_band_allows():
    weak = SequenceRecord("AT" * 40)
    taken = overlap_after(weak, 10, NEB_RULE)
    assert len(taken) == NEB_RULE.longest
    assert wallace_tm(taken) < (NEB_RULE.tm_floor or 0.0)


def test_a_product_stating_no_melting_temperature_takes_the_shortest_its_band_allows():
    silent = OverlapRule(15, 20, None)
    assert len(overlap_after(SequenceRecord("AT" * 40), 10, silent)) == 15


def test_an_overlap_is_read_across_the_origin_of_a_circular_vector():
    circle = SequenceRecord("GCGCGCGC" + "A" * 40 + "CGCGCGCG", topology="circular")
    rule = OverlapRule(4, 6, 24.0)
    assert overlap_before(circle, 2, rule) == "CGCGGC"
    assert overlap_after(circle, len(circle) - 3, rule) == "GCGGCG"


def test_a_vector_too_short_for_the_band_is_refused_with_the_band_named():
    with pytest.raises(ValueError, match="15-20 bp"):
        overlap_before(SequenceRecord("ACGT" * 3), 4, NEB_RULE)


def test_the_tier_steps_with_the_fragment_count():
    two, five = NEBUILDER_HIFI.tier(3), NEBUILDER_HIFI.tier(5)
    # Note §3, §7 and §8: the overlap, the incubation and the ratio all move together.
    assert (two.overlap.shortest, two.overlap.longest) == (15, 20)
    assert (five.overlap.shortest, five.overlap.longest) == (20, 30)
    assert (two.incubation_seconds, five.incubation_seconds) == (900, 3600)
    assert (two.insert_ratio, five.insert_ratio) == (2.0, 1.0)
    # More fragments than any tier covers still take the most the manual says.
    assert NEBUILDER_HIFI.tier(20) is five
    with pytest.raises(ValueError, match="at least two"):
        NEBUILDER_HIFI.tier(1)


def test_the_vector_goes_in_at_nebs_fifty_nanograms_and_the_insert_at_the_tiers_ratio():
    vector, insert = assembly_amounts(("backbone", 2629), [("insert", 717)], product=NEBUILDER_HIFI)
    assert vector.nanograms == pytest.approx(50.0, abs=0.1)
    assert insert.pmol == pytest.approx(2.0 * vector.pmol, rel=1e-3)
    # The picomoles stay inside the total NEB's table takes for two or three fragments.
    band = NEBUILDER_HIFI.tiers[0].total_pmol
    assert band is not None
    low, high = band
    assert low <= vector.pmol + insert.pmol <= high


def test_an_insert_under_two_hundred_bases_goes_in_at_a_fivefold_excess():
    vector, short = assembly_amounts(("backbone", 2629), [("linker", 120)], product=NEBUILDER_HIFI)
    assert short.pmol == pytest.approx(5.0 * vector.pmol, rel=1e-3)


def test_the_reaction_is_twenty_microlitres_of_which_half_is_master_mix():
    amounts = assembly_amounts(("backbone", 2629), [("insert", 717)], product=NEBUILDER_HIFI)
    table = assembly_reaction(NEBUILDER_HIFI, amounts)
    volumes = {one.name: one.volume_ul for one in table.components}
    assert sum(volumes.values()) == pytest.approx(NEBUILDER_HIFI.reaction_ul)
    assert volumes[NEBUILDER_HIFI.name] == NEBUILDER_HIFI.master_mix_ul == 10.0
    assert NEBUILDER_HIFI.unpurified_ul == 4.0
    with pytest.raises(ValueError, match="at least two"):
        assembly_reaction(NEBUILDER_HIFI, amounts[:1])


def test_in_fusions_table_prints_both_the_picomoles_and_the_weight_of_every_fragment():
    amounts = assembly_amounts(("backbone", 2629), [("GFP", 717)], product=IN_FUSION)
    table = assembly_reaction(IN_FUSION, amounts)
    volumes = {one.name: one.volume_ul for one in table.components}
    assert sum(volumes.values()) == pytest.approx(IN_FUSION.reaction_ul)
    assert volumes[IN_FUSION.name] == IN_FUSION.master_mix_ul
    # Takara asks for its inserts by weight and by picomoles, so the table prints both.
    insert = next(one for one in table.components if one.name == "GFP")
    assert "pmol" in insert.final
    assert "ng" in insert.final
    assert amounts[1].pmol == pytest.approx(2.0 * amounts[0].pmol, rel=1e-3)


def test_the_incubation_is_one_temperature_and_no_cycling():
    program = assembly_program(NEBUILDER_HIFI, fragments=2)
    assembly, hold = program.stages
    assert [stage.cycles for stage in program.stages] == [1, 1]
    assert (assembly.incubations[0].temperature_c, assembly.incubations[0].seconds) == (50.0, 900)
    assert hold.incubations[0].seconds is None
    assert assembly_program(NEBUILDER_HIFI, fragments=5).stages[0].incubations[0].seconds == 3600


def test_the_default_product_cites_the_documents_its_numbers_came_from():
    assert isinstance(NEBUILDER_HIFI, AssemblyProduct)
    assert NEBUILDER_HIFI.inserts_limit == 5
    assert NEBUILDER_HIFI.shortest_overlap_bp == 12
    citations = " ".join(reference.text for reference in NEBUILDER_HIFI.references)
    assert "E2621" in citations
    assert "CC BY" in citations
    assert all(isinstance(tier, Tier) for tier in NEBUILDER_HIFI.tiers)


def test_each_product_carries_its_own_rules_and_in_fusion_inherits_none_of_nebs():
    # Note §3, §6, §7 and §8: the bands, the reaction and the totals differ by product.
    neb, gibson, takara = (product.tiers for product in ASSEMBLY_PRODUCTS)
    assert [(one.overlap.shortest, one.overlap.longest) for one in neb] == [(15, 20), (20, 30)]
    assert [(one.overlap.shortest, one.overlap.longest) for one in gibson] == [(15, 25), (20, 40)]
    assert [one.total_pmol for one in gibson] == [(0.02, 0.5), (0.2, 1.0)]
    # In-Fusion states one length for one insert and another above two fragments, no melting
    # temperature for either, and the same incubation whatever the fragment count.
    assert [(one.overlap.shortest, one.overlap.longest) for one in takara] == [(15, 15), (20, 20)]
    assert {one.overlap.tm_floor for one in takara} == {None}
    assert {one.incubation_seconds for one in takara} == {900}
    assert {one.total_pmol for one in takara} == {None}
    assert IN_FUSION.tier(2) is takara[0]
    assert IN_FUSION.tier(3) is takara[1]
    # Its reaction is 10 µL of which 2 µL is a 5X mix, and it documents no unpurified allowance.
    assert (IN_FUSION.reaction_ul, IN_FUSION.master_mix_ul, IN_FUSION.master_mix_fold) == (
        10.0,
        2.0,
        "5X",
    )
    assert IN_FUSION.unpurified_fraction is None
    assert IN_FUSION.unpurified_ul is None
    assert IN_FUSION.supplier == "Takara Bio"
    assert all(product.inserts_limit == 5 for product in ASSEMBLY_PRODUCTS)
    for product in ASSEMBLY_PRODUCTS:
        assert product.references
        assert all(one.url for one in product.references)


def test_a_product_is_named_by_the_beginning_of_its_name():
    assert assembly_product("nebuilder") is NEBUILDER_HIFI
    assert assembly_product("Gibson") is GIBSON_MASTER_MIX
    assert assembly_product("in-fusion") is IN_FUSION
    with pytest.raises(ValueError, match="no assembly product"):
        assembly_product("KLD")


def test_more_inserts_than_the_product_documents_warns_and_names_the_limit():
    assert fragment_check(NEBUILDER_HIFI, 5).status == "pass"
    over = fragment_check(NEBUILDER_HIFI, 6)
    assert (over.status, over.value) == ("warn", 6)
    assert "5 or fewer" in over.detail
    assert NEBUILDER_HIFI.supplier in over.detail


def test_the_dna_is_judged_only_where_the_supplier_gives_a_picomole_band():
    amounts = assembly_amounts(("backbone", 2629), [("insert", 717)], product=NEBUILDER_HIFI)
    judged = assembly_dna_check(NEBUILDER_HIFI, amounts)
    assert (judged.status, judged.value) == ("pass", 0.0927)
    assert "0.03-0.2 pmol" in judged.detail
    # Takara states a mass and no band, so nothing judges the same DNA.
    unjudged = assembly_dna_check(IN_FUSION, amounts)
    assert unjudged.status is None
    assert "no picomole band" in unjudged.detail


def test_every_overlap_is_judged_on_the_length_and_melting_temperature_the_product_states():
    named = (("one", "TAAAACGACGGCCAGT"), ("two", "GGCGTAATCATGGTCA"))
    judged = {one.name: one for one in overlap_checks(named, NEBUILDER_HIFI)}
    assert judged["overlap length"].status == "pass"
    assert "15-20 bp" in judged["overlap length"].detail
    assert judged["overlap tm"].status == "pass"
    # Below the band warns; below the shortest overlap the product documents at all fails.
    short = overlap_checks((("one", "ACGTACGTACGTAC"), ("two", "GGCGTAATCATGGTCA")), NEBUILDER_HIFI)
    assert next(one for one in short if one.name == "overlap length").status == "warn"
    tiny = overlap_checks((("one", "ACGTACGTACG"), ("two", "GGCGTAATCATGGTCA")), NEBUILDER_HIFI)
    assert next(one for one in tiny if one.name == "overlap length").status == "fail"
    # A junction too cool for NEB's floor warns and is named.
    cool = overlap_checks(
        (("one", "ATATATATATATATAT"), ("two", "GGCGTAATCATGGTCA")), NEBUILDER_HIFI
    )
    assert next(one for one in cool if one.name == "overlap tm").status == "warn"
    # In-Fusion states no floor, so nothing judges the same overlaps.
    assert (
        next(one for one in overlap_checks(named, IN_FUSION) if one.name == "overlap tm").status
        is None
    )


def test_a_repeat_or_a_palindrome_at_a_junction_is_reported_and_names_it():
    named = (("his tag", "CACCACCACCACGTTA"), ("clean", "GGCGTAATCATGGTCA"))
    judged = {one.name: one for one in overlap_checks(named, NEBUILDER_HIFI)}
    assert judged["overlap repeat"].status == "warn"
    assert "his tag 12 bp" in judged["overlap repeat"].detail
    assert "clean 0 bp" in judged["overlap repeat"].detail
    folded = (("hairpin", "TTGCGGCCGCAAGGTA"), ("clean", "GGCGTAATCATGGTCA"))
    assert {one.name: one.status for one in overlap_checks(folded, NEBUILDER_HIFI)}[
        "overlap hairpin"
    ] == "warn"
    assert judged["overlap hairpin"].status == "pass"


def test_two_overlaps_alike_are_reported_and_judged_by_nothing():
    same = (("one", "GGCGTAATCATGGTCA"), ("two", "GGCGTAATCATGGTCA"))
    judged = {one.name: one for one in overlap_checks(same, NEBUILDER_HIFI)}
    alike = judged["overlap similarity"]
    assert alike.status is None
    assert alike.value == 16
    assert "one and two" in alike.detail
    assert "nobody quantifies" in alike.detail
    # GC is measured for every junction and judged by nothing either.
    assert judged["overlap gc"].status is None
    assert "no source quantifies a GC band" in judged["overlap gc"].detail


def test_stitching_oligos_tile_both_strands_at_the_shape_the_note_states():
    part = "".join("ACGTTGCA"[index % 8] for index in range(150))
    oligos = stitch_oligos(part, name="linker", product=NEBUILDER_HIFI)
    # n oligos of length L overlapping by v tile n(L-v) + v bases, so 150 takes four 60-mers.
    assert len(oligos) == 4
    assert [one.name for one in oligos] == [f"linker oligo {number}" for number in (1, 2, 3, 4)]
    assert [one.strand for one in oligos] == [
        Strand.FORWARD,
        Strand.REVERSE,
        Strand.FORWARD,
        Strand.REVERSE,
    ]
    assert max(len(one.sequence) for one in oligos) <= STITCH_OLIGO_BASES
    # Every neighbour overlaps by exactly the note's 20 bp, and the set tiles the whole part.
    assert [one.end - other.start for one, other in pairwise(oligos)] == [STITCH_OVERLAP_BP] * 3
    assert (oligos[0].start, oligos[-1].end) == (0, len(part))
    for one in oligos:
        written = one.sequence if one.strand == Strand.FORWARD else reverse_complement(one.sequence)
        assert written == part[one.start : one.end]


def test_a_short_part_takes_two_oligos_so_the_set_still_tiles_both_strands():
    oligos = stitch_oligos("AT" * 30, product=NEBUILDER_HIFI)
    assert [one.strand for one in oligos] == [Strand.FORWARD, Strand.REVERSE]


def test_a_part_too_long_to_stitch_is_refused_with_the_ceiling_and_its_source():
    assert stitch_limit_bp() == STITCH_OLIGOS * (STITCH_OLIGO_BASES - STITCH_OVERLAP_BP) + (
        STITCH_OVERLAP_BP
    )
    assert len(stitch_oligos("A" * stitch_limit_bp(), product=NEBUILDER_HIFI)) == STITCH_OLIGOS
    with pytest.raises(ValueError, match="too long to stitch") as raised:
        stitch_oligos("A" * (stitch_limit_bp() + 1), name="promoter", product=NEBUILDER_HIFI)
    said = str(raised.value)
    assert f"at most {STITCH_OLIGOS} {STITCH_OLIGO_BASES}-base oligos" in said
    assert f"tile {stitch_limit_bp()} bases" in said
    assert "Gibson 2011" in said


def test_a_part_no_longer_than_one_overlap_cannot_be_tiled():
    with pytest.raises(ValueError, match="no longer than"):
        stitch_oligos("A" * STITCH_OVERLAP_BP, product=NEBUILDER_HIFI)


def test_a_bridging_oligo_takes_its_homology_from_the_end_of_each_fragment():
    before, after = "AAAACCCCGGGGTTTTACGT" * 2, "TTTTGGGGCCCCAAAATGCA" * 2
    oligo = bridging_oligo(before, after, name="bridge", product=NEBUILDER_HIFI)
    assert oligo.homology_bp == BRIDGE_HOMOLOGY_BP == 20
    assert oligo.sequence == before[-20:] + after[:20]
    assert len(oligo.sequence) == 2 * BRIDGE_HOMOLOGY_BP


def test_a_fragment_too_short_for_the_homology_is_refused():
    with pytest.raises(ValueError, match="too short for the 20 bp of homology"):
        bridging_oligo("ACGT", "T" * 40, product=NEBUILDER_HIFI)


@pytest.mark.parametrize("product", [NEBUILDER_HIFI, GIBSON_MASTER_MIX])
def test_both_neb_products_document_single_stranded_oligos_in_the_reaction(product):
    assert product.single_stranded_oligos
    requires_oligos(product, "stitch")


def test_in_fusion_supports_neither_route_and_says_which_products_do():
    assert not IN_FUSION.single_stranded_oligos
    for call in (
        lambda: stitch_oligos("AT" * 40, product=IN_FUSION),
        lambda: bridging_oligo("A" * 40, "C" * 40, product=IN_FUSION),
    ):
        with pytest.raises(ValueError, match="no single-stranded oligo") as raised:
            call()
        said = str(raised.value)
        assert IN_FUSION.name in said
        assert NEBUILDER_HIFI.name in said
        assert GIBSON_MASTER_MIX.name in said


def test_a_stitched_part_outside_addgenes_window_is_told_and_not_refused():
    def sized(bases: int) -> Part:
        record = SequenceRecord("AC" * bases, name="part")
        return stitch(record, stitch_oligos(record.sequence, product=NEBUILDER_HIFI))

    low, high = STITCH_WINDOW_BP
    assert stitch_checks([sized(low // 2), sized(high // 2)])[0].status == "pass"
    outside = stitch_checks([sized((high + 2) // 2)])[0]
    assert outside.status == "warn"
    assert f"between {low} and {high} bp" in outside.detail
    # Nothing stitched, nothing to judge.
    assert stitch_checks([]) == ()
