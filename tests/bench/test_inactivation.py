"""Heat inactivation, read from the enzyme record."""

from liulab_mbio.bench.inactivation import heat_inactivation, heat_inactivations
from liulab_mbio.enzymes import get_enzyme


def test_heat_inactivation_comes_from_the_enzyme_record() -> None:
    program = heat_inactivation(get_enzyme("BsaI"))
    assert program is not None
    step = program.stages[0].incubations[0]
    assert (step.temperature_c, step.seconds) == (80.0, 1200)


def test_an_enzyme_with_no_recorded_heat_inactivation_gets_no_program() -> None:
    assert heat_inactivation(get_enzyme("AarI")) is None


def test_one_tube_kills_enzymes_that_agree_together() -> None:
    programs = heat_inactivations([get_enzyme(name) for name in ("BsaI", "EcoRI", "BsmBI")])
    assert [program.title for program in programs] == [
        "BsaI and BsmBI heat inactivation",
        "EcoRI heat inactivation",
    ]


def test_a_tube_of_enzymes_nobody_states_an_inactivation_for_gets_no_program() -> None:
    assert heat_inactivations([get_enzyme("AarI")]) == ()


def test_a_hotter_hold_does_not_stand_in_for_a_cooler_one() -> None:
    """BsaI-HFv2 and PmeI in one tube are held twice, 80 °C then 65 °C.

    `docs/research/heat-inactivation-order.md`: NEB assays each enzyme at one temperature, and
    nothing published says the 80 °C hold inactivates an enzyme listed at 65 °C. Collapsing the
    two needs a source, so the tube keeps both holds.
    """
    programs = heat_inactivations([get_enzyme(name) for name in ("BsaI-HFv2", "PmeI")])
    held = [
        (program.stages[0].incubations[0].temperature_c, program.stages[0].incubations[0].seconds)
        for program in programs
    ]
    assert held == [(80.0, 1200), (65.0, 1200)]
