"""Heat inactivation, read from the enzyme record."""

from liulab_mbio.bench.inactivation import heat_inactivation
from liulab_mbio.enzymes import get_enzyme


def test_heat_inactivation_comes_from_the_enzyme_record() -> None:
    program = heat_inactivation(get_enzyme("BsaI"))
    assert program is not None
    step = program.stages[0].incubations[0]
    assert (step.temperature_c, step.seconds) == (80.0, 1200)


def test_an_enzyme_with_no_recorded_heat_inactivation_gets_no_program() -> None:
    assert heat_inactivation(get_enzyme("AarI")) is None
