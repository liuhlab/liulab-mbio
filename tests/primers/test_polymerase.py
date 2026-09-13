"""Pinned Tm and Ta values come from NEB's own calculator and API, through
``docs/research/primer-design-and-pcr.md``.
"""

import pytest

from liulab_mbio.primers import ONETAQ, PHUSION, Q5, TAQ, melting_temperature

from .sequences import M13_FWD, M13_REV, PUC_FWD, PUC_REV


@pytest.mark.parametrize(
    ("polymerase", "forward_tm", "reverse_tm"),
    [
        (Q5, 62.27, 56.31),
        (PHUSION, 56.57, 50.89),
        (TAQ, 53.96, 47.99),
        (ONETAQ, 53.82, 47.85),
    ],
    ids=lambda value: getattr(value, "name", value),
)
def test_tm_matches_nebs_calculator_for_the_polymerase(polymerase, forward_tm, reverse_tm) -> None:
    assert melting_temperature(M13_FWD, polymerase) == pytest.approx(forward_tm, abs=0.01)
    assert melting_temperature(M13_REV, polymerase) == pytest.approx(reverse_tm, abs=0.01)


def test_tm_defaults_to_q5() -> None:
    assert melting_temperature(M13_REV) == melting_temperature(M13_REV, Q5)


@pytest.mark.parametrize(
    ("polymerase", "annealing_temperature"),
    [(Q5, 57.3), (PHUSION, 54.8), (TAQ, 43.0), (ONETAQ, 42.8)],
    ids=lambda value: getattr(value, "name", value),
)
def test_annealing_temperature_follows_nebs_rule_for_the_polymerase(
    polymerase, annealing_temperature
) -> None:
    tms = (melting_temperature(M13_FWD, polymerase), melting_temperature(M13_REV, polymerase))
    assert polymerase.annealing_temperature(*tms) == pytest.approx(annealing_temperature, abs=0.1)


def test_the_annealing_temperature_is_capped() -> None:
    long_pair = ("CGCCAGGGTTTTCCCAGTCACGACGTTG", "GCGGATAACAATTTCACACAGGAAACAGCTATGAC")
    tms = tuple(melting_temperature(primer, Q5) for primer in long_pair)
    assert Q5.annealing_temperature(*tms) == 72.0
    assert TAQ.annealing_temperature(80.0, 75.0) == 68.0


def test_the_twenty_three_mer_pair_anneals_above_nebs_floor_in_onetaq() -> None:
    tms = (melting_temperature(PUC_FWD, ONETAQ), melting_temperature(PUC_REV, ONETAQ))
    assert ONETAQ.annealing_temperature(*tms) == pytest.approx(51.5, abs=0.1)
    assert ONETAQ.annealing_temperature(*tms) > ONETAQ.annealing_min


def test_extension_time_rounds_the_amplicon_up_to_whole_kilobases() -> None:
    assert Q5.extension_seconds(2686) == 60
    assert Q5.extension_seconds(1000) == 20
    assert PHUSION.extension_seconds(1000) == 15
    assert TAQ.extension_seconds(103) == 60
    assert ONETAQ.extension_seconds(1001) == 120
