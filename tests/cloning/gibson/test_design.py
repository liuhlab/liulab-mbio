"""The overlap rule, as units, against the numbers `docs/research/gibson-assembly.md` states.

These are the fast tests, and the ones carrying the note's values: the Wallace rule of §4, the
NEBuilder bands of §3, the tiers of §7 and §9, and the amounts and ratios of §8.
"""

import pytest

from liulab_mbio.cloning.gibson.bench import (
    NEBUILDER_HIFI,
    AssemblyProduct,
    OverlapRule,
    Tier,
    assembly_amounts,
    assembly_program,
    assembly_reaction,
)
from liulab_mbio.cloning.gibson.design import overlap_after, overlap_before, wallace_tm
from liulab_mbio.sequence import SequenceRecord

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
    vector, insert = assembly_amounts(
        ("backbone", 2629), [("insert", 717)], tier=NEBUILDER_HIFI.tiers[0]
    )
    assert vector.nanograms == pytest.approx(50.0, abs=0.1)
    assert insert.pmol == pytest.approx(2.0 * vector.pmol, rel=1e-3)
    # The picomoles stay inside the total NEB's table takes for two or three fragments.
    low, high = NEBUILDER_HIFI.tiers[0].total_pmol
    assert low <= vector.pmol + insert.pmol <= high


def test_an_insert_under_two_hundred_bases_goes_in_at_a_fivefold_excess():
    vector, short = assembly_amounts(
        ("backbone", 2629), [("linker", 120)], tier=NEBUILDER_HIFI.tiers[0]
    )
    assert short.pmol == pytest.approx(5.0 * vector.pmol, rel=1e-3)


def test_the_reaction_is_twenty_microlitres_of_which_half_is_master_mix():
    amounts = assembly_amounts(("backbone", 2629), [("insert", 717)], tier=NEBUILDER_HIFI.tiers[0])
    table = assembly_reaction(NEBUILDER_HIFI, amounts)
    volumes = {one.name: one.volume_ul for one in table.components}
    assert sum(volumes.values()) == pytest.approx(NEBUILDER_HIFI.reaction_ul)
    assert volumes[NEBUILDER_HIFI.name] == NEBUILDER_HIFI.master_mix_ul == 10.0
    assert NEBUILDER_HIFI.unpurified_ul == 4.0
    with pytest.raises(ValueError, match="at least two"):
        assembly_reaction(NEBUILDER_HIFI, amounts[:1])


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
