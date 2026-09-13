import pytest

from liulab_mbio.primers import ONETAQ, Q5, TAQ, melting_temperature

M13_FWD = "GTAAAACGACGGCCAGT"
M13_REV = "CAGGAAACAGCTATGAC"


@pytest.mark.parametrize(
    ("polymerase", "monovalent", "magnesium", "primer_nm"),
    [(Q5, 50.0, 2.0, 500.0), (TAQ, 50.0, 1.5, 200.0), (ONETAQ, 44.0, 1.8, 200.0)],
    ids=["Q5", "Taq", "OneTaq"],
)
def test_tm_is_santalucia_nearest_neighbour_in_the_polymerase_buffer(
    polymerase, monovalent, magnesium, primer_nm
) -> None:
    from Bio.SeqUtils import MeltingTemp

    # An independent implementation; primer3 takes a quarter of the primer concentration.
    expected = MeltingTemp.Tm_NN(
        M13_FWD,
        nn_table=MeltingTemp.DNA_NN3,
        Na=monovalent,
        Mg=magnesium,
        dNTPs=0.8,
        dnac1=primer_nm / 4,
        dnac2=0,
        saltcorr=5,
    )
    assert melting_temperature(M13_FWD, polymerase) == pytest.approx(expected, abs=0.05)


def test_tm_defaults_to_q5() -> None:
    assert melting_temperature(M13_REV) == melting_temperature(M13_REV, Q5)


def test_q5_anneals_3_degrees_above_the_lower_tm_and_never_above_72() -> None:
    assert Q5.annealing_temperature(61.0, 58.0) == 61.0
    assert Q5.annealing_temperature(71.0, 74.0) == 72.0


@pytest.mark.parametrize("polymerase", [TAQ, ONETAQ], ids=["Taq", "OneTaq"])
def test_taq_polymerases_anneal_5_degrees_below_the_lower_tm_and_never_above_68(
    polymerase,
) -> None:
    assert polymerase.annealing_temperature(61.0, 58.0) == 53.0
    assert polymerase.annealing_temperature(80.0, 75.0) == 68.0


def test_extension_time_rounds_the_amplicon_up_to_whole_kilobases() -> None:
    assert Q5.extension_seconds(2686) == 90
    assert Q5.extension_seconds(1000) == 30
    assert TAQ.extension_seconds(103) == 60
    assert ONETAQ.extension_seconds(1001) == 120
