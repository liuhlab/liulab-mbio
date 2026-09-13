"""Bench numbers pinned against NEB's own tables.

Sources are `docs/research/golden-gate-assembly.md` for the reactions and cycling, and
`docs/research/primer-design-and-pcr.md` for ng to pmol, the PCR programs and the ladders.
"""

import pytest

from liulab_mbio.goldengate.bench import (
    DNA_VOLUME_UL,
    Amount,
    Fragment,
    assembly_amounts,
    molecular_weight,
    to_nanograms,
    to_pmol,
)


def test_a_kilobase_pair_weighs_what_nebiocalculator_says() -> None:
    assert molecular_weight(1000) == pytest.approx(36.04 + 1000 * 615.94)


def test_one_microgram_of_puc19_is_the_picomoles_neb_quotes() -> None:
    # NEB Nucleic Acid Data gives 0.57 pmol by the 650 Da rule; NEBioCalculator gives 0.604.
    assert to_pmol(1000, 2686) == pytest.approx(0.604, abs=0.001)


def test_fifty_nanograms_of_five_kilobases_is_about_the_pmol_neb_quotes() -> None:
    assert to_pmol(50, 5000) == pytest.approx(0.016, abs=0.001)


def test_mass_and_moles_are_inverse() -> None:
    assert to_nanograms(to_pmol(750, 2686), 2686) == pytest.approx(750)


def test_every_fragment_takes_the_same_picomoles() -> None:
    amounts = assembly_amounts(Fragment("pUC19", 2686), (Fragment("GFP", 717),))
    assert [amount.pmol for amount in amounts] == [0.05, 0.05]
    assert [amount.name for amount in amounts] == ["pUC19", "GFP"]
    assert amounts[0].nanograms == pytest.approx(82.72, abs=0.01)


def test_a_molar_ratio_scales_the_inserts_and_not_the_vector() -> None:
    amounts = assembly_amounts(Fragment("pUC19", 2686), (Fragment("GFP", 717),), insert_ratio=2.0)
    assert amounts[0].pmol == pytest.approx(0.05)
    assert amounts[1].pmol == pytest.approx(0.10)


def test_a_known_concentration_gives_a_volume_to_pipette() -> None:
    amounts = assembly_amounts(Fragment("pUC19", 2686, concentration_ng_ul=100.0), ())
    assert amounts[0].volume_ul == pytest.approx(0.83, abs=0.01)


def test_an_unknown_concentration_falls_back_to_one_microlitre() -> None:
    amounts = assembly_amounts(Fragment("pUC19", 2686), ())
    assert amounts[0].volume_ul == DNA_VOLUME_UL


def test_a_fragment_refuses_a_length_that_is_not_positive() -> None:
    with pytest.raises(ValueError, match="length_bp"):
        Fragment("empty", 0)


def test_an_amount_refuses_picomoles_that_are_not_positive() -> None:
    with pytest.raises(ValueError, match="pmol"):
        Amount("x", 100, pmol=0.0, nanograms=1.0, volume_ul=1.0)
