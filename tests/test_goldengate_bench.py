"""Bench numbers pinned against NEB's own tables.

Sources are `docs/research/golden-gate-assembly.md` for the reactions and cycling, and
`docs/research/primer-design-and-pcr.md` for ng to pmol, the PCR programs and the ladders.
"""

import pytest

from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.goldengate.bench import (
    COLONY_PCR_MASTER_MIX,
    DNA_VOLUME_UL,
    KIT,
    Amount,
    Fragment,
    ReactionTable,
    assembly_amounts,
    assembly_program,
    assembly_reaction,
    colony_pcr_program,
    colony_pcr_reaction,
    golden_gate_temperature,
    heat_inactivation,
    molecular_weight,
    pcr_program,
    pcr_reaction,
    to_nanograms,
    to_pmol,
)
from liulab_mbio.primers import Q5, TAQ


def volumes(table: ReactionTable) -> dict[str, float]:
    return {component.name: component.volume_ul for component in table.components}


def total(table: ReactionTable) -> float:
    return sum(component.volume_ul for component in table.components)


def two_fragments() -> tuple[Amount, ...]:
    return assembly_amounts(Fragment("pUC19", 2686), (Fragment("GFP", 717),))


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
    amounts = two_fragments()
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


def test_the_ligase_master_mix_reaction_holds_fifteen_microlitres_for_two_fragments() -> None:
    table = assembly_reaction(get_enzyme("BbsI"), two_fragments())
    assert total(table) == pytest.approx(15.0)
    assert volumes(table)["NEBridge Ligase Master Mix"] == 5.0
    assert volumes(table)["BbsI-HF (R3539)"] == 1.0


def test_the_dna_goes_in_each_tube_and_the_rest_into_the_mix() -> None:
    table = assembly_reaction(get_enzyme("BbsI"), two_fragments())
    per_tube = [c.name for c in table.components if not c.master_mix]
    assert per_tube == ["pUC19", "GFP"]


def test_paqci_carries_an_activator_line_and_bsai_does_not() -> None:
    paqci = assembly_reaction(get_enzyme("PaqCI"), two_fragments())
    bsai = assembly_reaction(get_enzyme("BsaI"), two_fragments())
    assert volumes(paqci)["PaqCI Activator"] == 0.5
    assert not [c for c in bsai.components if "Activator" in c.name]


def test_seven_fragments_double_the_reaction_and_the_bsmbi_enzyme() -> None:
    amounts = assembly_amounts(
        Fragment("pUC19", 2686), tuple(Fragment(f"part {n}", 700) for n in range(6))
    )
    table = assembly_reaction(get_enzyme("BsmBI"), amounts)
    assert total(table) == pytest.approx(30.0)
    assert volumes(table)["NEBridge Ligase Master Mix"] == 10.0
    assert volumes(table)["BsmBI-v2 (R0739)"] == 6.0


def test_the_kit_reaction_is_twenty_microlitres_with_its_own_enzyme_mix() -> None:
    table = assembly_reaction(get_enzyme("BsaI"), two_fragments(), system=KIT)
    assert total(table) == pytest.approx(20.0)
    assert volumes(table)["T4 DNA Ligase Buffer"] == 2.0
    assert volumes(table)["NEBridge Golden Gate Enzyme Mix"] == 1.0
    assert "E1601" in table.title


def test_no_kit_carries_bbsi() -> None:
    with pytest.raises(ValueError, match="no NEBridge kit"):
        assembly_reaction(get_enzyme("BbsI"), two_fragments(), system=KIT)


def test_an_enzyme_neb_gives_no_golden_gate_protocol_for_is_refused() -> None:
    with pytest.raises(ValueError, match="no Golden Gate protocol"):
        assembly_reaction(get_enzyme("BtgZI"), two_fragments())


def test_dna_that_does_not_fit_the_reaction_is_refused() -> None:
    dilute = assembly_amounts(
        Fragment("pUC19", 2686, concentration_ng_ul=1.0),
        (Fragment("GFP", 717, concentration_ng_ul=1.0),),
    )
    with pytest.raises(ValueError, match="µL reaction"):
        assembly_reaction(get_enzyme("BbsI"), dilute)


def test_golden_gate_temperature_splits_two_isoschizomers() -> None:
    assert golden_gate_temperature(get_enzyme("Esp3I")) == 37.0
    assert golden_gate_temperature(get_enzyme("BsmBI")) == 42.0


def test_two_fragments_of_a_thirty_seven_degree_enzyme_are_one_incubation() -> None:
    program = assembly_program(get_enzyme("BbsI"), fragments=2)
    held = program.stages[0].incubations[0]
    assert (held.temperature_c, held.seconds) == (37.0, 900)
    assert program.stages[0].cycles == 1


def test_two_fragments_of_a_forty_two_degree_enzyme_cycle_fifteen_times() -> None:
    program = assembly_program(get_enzyme("BsmBI"), fragments=2)
    assert program.stages[0].cycles == 15
    assert [i.temperature_c for i in program.stages[0].incubations] == [42.0, 16.0]


def test_a_library_assembly_incubates_for_an_hour() -> None:
    program = assembly_program(get_enzyme("BbsI"), fragments=2, library=True)
    assert program.stages[0].incubations[0].seconds == 3600


def test_more_fragments_lengthen_and_then_multiply_the_cycles() -> None:
    tiers = {
        fragments: assembly_program(get_enzyme("BsaI"), fragments=fragments).stages[0]
        for fragments in (6, 13, 14)
    }
    assert (tiers[6].cycles, tiers[6].incubations[0].seconds) == (30, 60)
    assert (tiers[13].cycles, tiers[13].incubations[0].seconds) == (30, 300)
    assert (tiers[14].cycles, tiers[14].incubations[0].seconds) == (60, 300)


def test_every_assembly_program_ends_with_the_sixty_degree_soak() -> None:
    for fragments in (2, 4, 8):
        program = assembly_program(get_enzyme("PaqCI"), fragments=fragments)
        soak = program.stages[-1].incubations[-1]
        assert (soak.label, soak.temperature_c, soak.seconds) == ("End soak", 60.0, 300)


def test_the_kit_program_counts_inserts_and_not_fragments() -> None:
    single = assembly_program(get_enzyme("BsaI"), fragments=2, system=KIT).stages[0]
    eleven = assembly_program(get_enzyme("BsaI"), fragments=12, system=KIT).stages[0]
    assert (single.cycles, single.incubations[0].seconds) == (1, 300)
    assert (eleven.cycles, eleven.incubations[0].seconds) == (30, 300)


def test_heat_inactivation_comes_from_the_enzyme_record() -> None:
    program = heat_inactivation(get_enzyme("BsaI"))
    assert program is not None
    step = program.stages[0].incubations[0]
    assert (step.temperature_c, step.seconds) == (80.0, 1200)


def test_an_enzyme_with_no_recorded_heat_inactivation_gets_no_program() -> None:
    assert heat_inactivation(get_enzyme("AarI")) is None


def test_the_q5_reaction_fills_fifty_microlitres_as_nebs_table_does() -> None:
    table = pcr_reaction(Q5)
    assert total(table) == pytest.approx(50.0)
    assert volumes(table) == {
        "Q5 Reaction Buffer": 10.0,
        "dNTP mix": 1.0,
        "Forward primer": 2.5,
        "Reverse primer": 2.5,
        "Template DNA": 1.0,
        "Q5 DNA Polymerase": 0.5,
        "Nuclease-free water": 32.5,
    }


def test_taq_takes_a_tenfold_buffer_and_less_enzyme() -> None:
    table = pcr_reaction(TAQ)
    assert total(table) == pytest.approx(50.0)
    assert volumes(table)["Standard Taq Reaction Buffer"] == 5.0
    assert volumes(table)["Forward primer"] == 1.0
    assert volumes(table)["Taq DNA Polymerase"] == 0.25


def test_a_pcr_program_anneals_at_the_pair_temperature_and_extends_by_length() -> None:
    program = pcr_program(Q5, annealing_temperature=57.3, amplicon_length=800)
    cycled = program.stages[1]
    assert cycled.cycles == 30
    assert [i.label for i in cycled.incubations] == ["Denature", "Anneal", "Extend"]
    assert cycled.incubations[1].temperature_c == 57.3
    assert cycled.incubations[2].seconds == 20
    assert program.stages[-1].incubations[0].seconds is None


def test_a_high_annealing_temperature_combines_annealing_and_extension() -> None:
    program = pcr_program(Q5, annealing_temperature=72.0, amplicon_length=800)
    cycled = program.stages[1]
    assert [i.label for i in cycled.incubations] == ["Denature", "Anneal and extend"]
    assert cycled.incubations[1].temperature_c == 72.0


def test_taq_goes_two_step_above_sixty_five_degrees_and_not_at_it() -> None:
    at = pcr_program(TAQ, annealing_temperature=65.0, amplicon_length=500)
    above = pcr_program(TAQ, annealing_temperature=65.1, amplicon_length=500)
    assert len(at.stages[1].incubations) == 3
    assert len(above.stages[1].incubations) == 2


def test_the_colony_pcr_reaction_is_half_master_mix() -> None:
    table = colony_pcr_reaction()
    assert total(table) == pytest.approx(25.0)
    assert volumes(table)[COLONY_PCR_MASTER_MIX] == 12.5
    assert volumes(table)["Forward primer"] == 0.5
    assert volumes(table)["Nuclease-free water"] == 11.5


def test_the_colony_pcr_program_lyses_the_cells_and_holds_at_ten() -> None:
    program = colony_pcr_program(annealing_temperature=51.5, amplicon_length=137)
    lysis = program.stages[0].incubations[0]
    assert (lysis.label, lysis.temperature_c, lysis.seconds) == ("Lysis", 94.0, 300)
    assert program.stages[1].incubations[2].seconds == 60
    assert program.stages[-1].incubations[0].temperature_c == 10.0
