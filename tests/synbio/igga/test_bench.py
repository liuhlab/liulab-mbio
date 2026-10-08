"""Library bench amounts and reactions, computed from the lengths and the sourced ratios.

The defaults are the method's own, through `docs/research/protein-library-assembly.md`.
"""

import pytest

from liulab_mbio.bench.amounts import dna_amount
from liulab_mbio.enzymes import get_enzyme
from liulab_synbio.igga.bench import (
    DIGEST_CELSIUS,
    DIGEST_SECONDS,
    DIGEST_SOURCE_KEY,
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    GROWTH_CELSIUS,
    LIGATION_VOLUME_UL,
    OUTGROWTH_SECONDS,
    RECOVERY_SECONDS,
    SPRI_AFTER_DIGEST,
    SPRI_AFTER_LIGATION,
    choppers,
    digest_amount,
    digest_program,
    digest_reaction,
    growth_program,
    ligation_amounts,
    ligation_reaction,
    transformation_amount,
)
from liulab_synbio.igga.method import IGGA


def test_a_digest_takes_one_microgram_as_picomoles_of_its_own_length() -> None:
    amount = digest_amount(("N part list", 5000))

    assert amount.nanograms == pytest.approx(1000.0)
    assert amount.pmol == pytest.approx(0.325, abs=0.001)


def test_a_digest_refuses_dna_too_dilute_for_its_reaction() -> None:
    with pytest.raises(ValueError, match=r"concentrate .* or scale the reaction up"):
        digest_amount(("N part list", 5000), concentration_ng_ul=10.0)


def test_a_digest_fits_when_the_dna_is_concentrated_enough() -> None:
    assert digest_amount(("N part list", 5000), concentration_ng_ul=100.0).volume_ul == 10.0


def test_the_ligation_matches_molecules_and_not_masses() -> None:
    opened, released = ligation_amounts(("library", 5000), ("C part list", 1200))

    assert opened.pmol == pytest.approx(released.pmol)
    assert released.nanograms < opened.nanograms
    assert opened.nanograms == pytest.approx(20.0)


def test_a_wider_molar_ratio_gives_the_donor_more_molecules() -> None:
    _, one = ligation_amounts(("library", 5000), ("C part list", 1200))
    _, three = ligation_amounts(("library", 5000), ("C part list", 1200), ratio=3.0)

    assert three.pmol == pytest.approx(3.0 * one.pmol)


def test_the_ligation_refuses_dna_that_does_not_fit_its_volume() -> None:
    with pytest.raises(ValueError, match="exceeds the 200 µL reaction"):
        ligation_amounts(
            ("library", 5000), ("C part list", 1200), destination_ng_ul=0.05, donor_ng_ul=0.05
        )


def test_the_ligation_refuses_a_ratio_that_is_not_positive() -> None:
    with pytest.raises(ValueError, match="molar ratio must be positive"):
        ligation_amounts(("library", 5000), ("C part list", 1200), ratio=0.0)


def test_an_electroporation_takes_at_most_a_hundred_nanograms() -> None:
    assert transformation_amount(("round 1 library", 6200)).nanograms == pytest.approx(100.0)


def test_the_two_spri_ratios_differ() -> None:
    assert (SPRI_AFTER_DIGEST, SPRI_AFTER_LIGATION) == (2.0, 1.0)


def test_a_digest_fills_its_volume_and_gives_each_enzyme_its_own_line() -> None:
    dna = dna_amount("library", 5000, pmol=0.3)
    enzymes = (get_enzyme("BsaI"), get_enzyme("SrfI"))

    table = digest_reaction(dna, enzymes)

    assert round(sum(one.volume_ul for one in table.components), 2) == DIGEST_VOLUME_UL
    assert [one.volume_ul for one in table.components[1:3]] == [ENZYME_UL, ENZYME_UL]
    # The DNA goes in each tube, the rest into the mix.
    assert table.components[0].master_mix is False


def test_a_digest_line_is_worth_a_unit_count_its_catalogue_number_determines() -> None:
    """The paper gives volumes; NEB's specification for the product it names gives the rest."""
    enzymes = (get_enzyme("BsaI"), get_enzyme("PmeI"))

    lines = digest_reaction(dna_amount("library", 5000, pmol=0.3), enzymes).components[1:3]

    assert [one.name for one in lines] == ["BsaI-HFv2 (R3733L)", "PmeI (R0560L)"]
    assert [one.stock for one in lines] == ["20 U/µL", "10 U/µL"]
    # 2.5 µL of each, so 50 units of the 20 U/µL product and 25 of the 10 U/µL one.
    assert [one.final for one in lines] == ["50 units", "25 units"]
    assert [one.citation and one.citation.source for one in lines] == [DIGEST_SOURCE_KEY] * 2


def test_a_digest_runs_the_first_enzyme_alone_and_then_the_second() -> None:
    program = digest_program((get_enzyme("BsaI"), get_enzyme("SrfI")))

    assert len(program.stages) == 2
    assert program.duration_seconds == 2 * DIGEST_SECONDS
    for stage in program.stages:
        for incubation in stage.incubations:
            assert incubation.temperature_c == DIGEST_CELSIUS
    assert program.stages[0].incubations[0].label == "BsaI"
    assert "SrfI" in program.stages[1].incubations[0].label


def test_a_ligation_fills_its_volume_and_leaves_the_ligase_to_the_supplier() -> None:
    amounts = (
        dna_amount("library, opened", 5000, pmol=0.006),
        dna_amount("part, released", 400, pmol=0.006),
    )

    table = ligation_reaction(amounts)

    assert round(sum(one.volume_ul for one in table.components), 2) == LIGATION_VOLUME_UL
    assert [one.master_mix for one in table.components] == [False, False, True]
    assert "T7 DNA Ligase" in table.components[-1].name


def test_both_growth_steps_run_at_thirty_degrees() -> None:
    program = growth_program()

    temperatures = {
        incubation.temperature_c for stage in program.stages for incubation in stage.incubations
    }
    assert temperatures == {GROWTH_CELSIUS}
    assert GROWTH_CELSIUS == 30.0
    assert program.duration_seconds == RECOVERY_SECONDS + OUTGROWTH_SECONDS[0]


def test_each_blunt_enzyme_belongs_to_the_digest_whose_piece_it_cuts() -> None:
    inside, outside = choppers(IGGA)

    assert [one.name for one in inside] == ["SrfI"]
    assert [one.name for one in outside] == ["PmeI"]
