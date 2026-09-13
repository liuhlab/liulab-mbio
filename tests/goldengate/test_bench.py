"""Golden Gate bench numbers pinned against NEB's own tables.

The source is `docs/research/golden-gate-assembly.md` for the reactions and cycling.
"""

import pytest

from liulab_mbio.bench.amounts import Amount, dna_amount
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.goldengate.bench import (
    FRAGMENT_PMOL,
    KIT,
    assembly_amounts,
    assembly_program,
    assembly_reaction,
    enzyme_component,
    golden_gate_temperature,
    ligase_master_mix_component,
)

from ..reactions import total, volumes


def two_fragments() -> tuple[Amount, ...]:
    return assembly_amounts(("pUC19", 2686), (("GFP", 717),))


def test_every_fragment_takes_the_same_picomoles() -> None:
    amounts = two_fragments()
    assert [amount.pmol for amount in amounts] == [0.05, 0.05]
    assert [amount.name for amount in amounts] == ["pUC19", "GFP"]
    assert amounts[0].nanograms == pytest.approx(82.72, abs=0.01)


def test_a_molar_ratio_scales_the_inserts_and_not_the_vector() -> None:
    amounts = assembly_amounts(("pUC19", 2686), (("GFP", 717),), insert_ratio=2.0)
    assert amounts[0].pmol == pytest.approx(0.05)
    assert amounts[1].pmol == pytest.approx(0.10)


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
    amounts = assembly_amounts(("pUC19", 2686), tuple((f"part {n}", 700) for n in range(6)))
    table = assembly_reaction(get_enzyme("BsmBI"), amounts)
    assert total(table) == pytest.approx(30.0)
    assert volumes(table)["NEBridge Ligase Master Mix"] == 10.0
    assert volumes(table)["BsmBI-v2 (R0739)"] == 6.0


@pytest.mark.parametrize(("name", "fragments"), [("BbsI", 2), ("PaqCI", 4), ("BsmBI", 7)])
def test_the_master_mix_and_enzyme_components_are_the_ones_the_reaction_prints(
    name: str, fragments: int
) -> None:
    enzyme = get_enzyme(name)
    amounts = assembly_amounts(
        ("pUC19", 2686), tuple((f"part {n}", 700) for n in range(fragments - 1))
    )
    components = assembly_reaction(enzyme, amounts).components
    assert ligase_master_mix_component(fragments) in components
    assert enzyme_component(enzyme, fragments) in components


def test_an_enzyme_with_no_row_in_nebs_table_has_no_component() -> None:
    with pytest.raises(ValueError, match="no row"):
        enzyme_component(get_enzyme("BtgZI"), 2)


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
    dilute = (
        dna_amount("pUC19", 2686, pmol=FRAGMENT_PMOL, concentration_ng_ul=1.0),
        dna_amount("GFP", 717, pmol=FRAGMENT_PMOL, concentration_ng_ul=1.0),
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
