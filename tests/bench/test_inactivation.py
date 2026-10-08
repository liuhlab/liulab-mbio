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
