"""Library bench amounts, computed from the lengths, with the two sourced traps pinned.

The defaults are the method's own, through `docs/research/protein-library-assembly.md`.
"""

import pytest

from liulab_mbio.library.bench import (
    GROWTH_CELSIUS,
    SPRI_AFTER_DIGEST,
    SPRI_AFTER_LIGATION,
    digest_amount,
    ligation_amounts,
    transformation_amount,
)


def test_a_digest_takes_one_microgram_as_picomoles_of_its_own_length() -> None:
    amount = digest_amount(("N part list", 5000))

    assert amount.nanograms == pytest.approx(1000.0)
    assert amount.pmol == pytest.approx(0.325, abs=0.001)


def test_a_longer_plasmid_is_fewer_picomoles_for_the_same_microgram() -> None:
    assert digest_amount(("long", 10000)).pmol < digest_amount(("short", 5000)).pmol


def test_a_digest_refuses_dna_too_dilute_for_its_reaction() -> None:
    with pytest.raises(ValueError, match=r"concentrate .* or scale the reaction up"):
        digest_amount(("N part list", 5000), concentration_ng_ul=10.0)


def test_a_digest_fits_when_the_dna_is_concentrated_enough() -> None:
    assert digest_amount(("N part list", 5000), concentration_ng_ul=100.0).volume_ul == 10.0


def test_the_ligation_matches_molecules_and_not_masses() -> None:
    opened, released = ligation_amounts(("library", 5000), ("C part list", 1200))

    assert opened.pmol == pytest.approx(released.pmol)
    assert released.nanograms < opened.nanograms


def test_the_ligation_gives_the_destination_twenty_nanograms() -> None:
    opened, _ = ligation_amounts(("library", 5000), ("C part list", 1200))

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


def test_both_growth_steps_run_at_thirty_degrees() -> None:
    assert GROWTH_CELSIUS == 30.0


def test_the_two_spri_ratios_differ() -> None:
    assert (SPRI_AFTER_DIGEST, SPRI_AFTER_LIGATION) == (2.0, 1.0)
