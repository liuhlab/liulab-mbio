"""Bench numbers pinned against NEB's own tables.

Sources are `docs/research/golden-gate-assembly.md` for the reactions and cycling, and
`docs/research/primer-design-and-pcr.md` for ng to pmol, the PCR programs and the ladders.
"""

from pathlib import Path

import pytest

from liulab_mbio import edits
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.goldengate.bench import (
    COLONY_PCR_MASTER_MIX,
    DNA_VOLUME_UL,
    KIT,
    LADDER_1_KB_PLUS,
    LADDER_100_BP,
    Amount,
    ColonyCheck,
    Fragment,
    ReactionTable,
    agarose_percent,
    assembly_amounts,
    assembly_program,
    assembly_reaction,
    choose_ladder,
    colony_pcr_check,
    colony_pcr_program,
    colony_pcr_reaction,
    golden_gate_temperature,
    heat_inactivation,
    molecular_weight,
    pcr_program,
    pcr_reaction,
    sanger_primers,
    to_nanograms,
    to_pmol,
)
from liulab_mbio.io import read_record
from liulab_mbio.primers import ONETAQ, PHUSION, Q5, TAQ, Polymerase, amplicon_sizes
from liulab_mbio.sequence import Primer, SequenceRecord

DATA = Path(__file__).parent / "data"

#: The multiple cloning site of the fixture, 0-based and half-open, and the GFP that replaces it.
MCS = (395, 452)
#: Addgene's 23-mer M13/pUC pair, which the note recommends over the 17-mers for colony PCR.
M13_FORWARD = Primer("M13/pUC Forward", "CCCAGTCACGACGTTGTAAAACG")
M13_REVERSE = Primer("M13/pUC Reverse", "AGCGGATAACAATTTCACACAGG")


@pytest.fixture(scope="module")
def puc19() -> SequenceRecord:
    return read_record(DATA / "pUC19.dna")


@pytest.fixture(scope="module")
def gfp() -> SequenceRecord:
    return read_record(DATA / "GFP.dna")


@pytest.fixture(scope="module")
def product(puc19: SequenceRecord, gfp: SequenceRecord) -> SequenceRecord:
    """pUC19 with GFP in place of the whole multiple cloning site."""
    edited, _ = edits.replace(puc19, *MCS, gfp.sequence)
    return edited


@pytest.fixture(scope="module")
def junctions(gfp: SequenceRecord) -> tuple[int, int]:
    return (MCS[0], MCS[0] + len(gfp))


@pytest.fixture(scope="module")
def three_junctions(gfp: SequenceRecord) -> tuple[int, int, int]:
    """The same product read as two inserts, split inside GFP."""
    return (MCS[0], MCS[0] + 350, MCS[0] + len(gfp))


def volumes(table: ReactionTable) -> dict[str, float]:
    return {component.name: component.volume_ul for component in table.components}


def total(table: ReactionTable) -> float:
    return sum(component.volume_ul for component in table.components)


def bands(check: ColonyCheck, name: str) -> tuple[int, ...]:
    return next(clone.bands_bp for clone in check.clones if clone.name == name)


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


#: NEB's 50 µL PCR for each shipped polymerase, line by line as (name, µL, stock, final); its
#: program for 1.5 kb annealing at 55 °C as (label, °C, seconds); and the lowest two-step Ta.
SHIPPED_PCR = {
    "Q5": (
        [
            ("Q5 Reaction Buffer", 10.0, "5X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 2.5, "10 µM", "500 nM"),
            ("Reverse primer", 2.5, "10 µM", "500 nM"),
            ("Template DNA", 1.0, "", ""),
            ("Q5 DNA Polymerase", 0.5, "2 U/µL", "1 units"),
            ("Nuclease-free water", 32.5, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 98.0, 30),
            ("Denature", 98.0, 10),
            ("Anneal", 55.0, 20),
            ("Extend", 72.0, 40),
            ("Final extension", 72.0, 120),
            ("Hold", 4.0, None),
        ],
        72.0,
    ),
    "Phusion": (
        [
            ("Phusion HF Buffer", 10.0, "5X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 2.5, "10 µM", "500 nM"),
            ("Reverse primer", 2.5, "10 µM", "500 nM"),
            ("Template DNA", 1.0, "", ""),
            ("Phusion DNA Polymerase", 0.5, "2 U/µL", "1 units"),
            ("Nuclease-free water", 32.5, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 98.0, 30),
            ("Denature", 98.0, 10),
            ("Anneal", 55.0, 20),
            ("Extend", 72.0, 30),
            ("Final extension", 72.0, 300),
            ("Hold", 4.0, None),
        ],
        72.0,
    ),
    "Taq": (
        [
            ("Standard Taq Reaction Buffer", 5.0, "10X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 1.0, "10 µM", "200 nM"),
            ("Reverse primer", 1.0, "10 µM", "200 nM"),
            ("Template DNA", 1.0, "", ""),
            ("Taq DNA Polymerase", 0.25, "5 U/µL", "1.25 units"),
            ("Nuclease-free water", 40.75, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 95.0, 30),
            ("Denature", 95.0, 30),
            ("Anneal", 55.0, 30),
            ("Extend", 68.0, 120),
            ("Final extension", 68.0, 300),
            ("Hold", 4.0, None),
        ],
        65.1,
    ),
    "OneTaq": (
        [
            ("OneTaq Standard Reaction Buffer", 10.0, "5X", "1X"),
            ("dNTP mix", 1.0, "10 mM each", "200 µM each"),
            ("Forward primer", 1.0, "10 µM", "200 nM"),
            ("Reverse primer", 1.0, "10 µM", "200 nM"),
            ("Template DNA", 1.0, "", ""),
            ("OneTaq DNA Polymerase", 0.25, "5 U/µL", "1.25 units"),
            ("Nuclease-free water", 35.75, "", "to 50 µL"),
        ],
        [
            ("Initial denaturation", 94.0, 30),
            ("Denature", 94.0, 30),
            ("Anneal", 55.0, 30),
            ("Extend", 68.0, 120),
            ("Final extension", 68.0, 300),
            ("Hold", 4.0, None),
        ],
        68.0,
    ),
}


@pytest.mark.parametrize("polymerase", [Q5, PHUSION, TAQ, ONETAQ], ids=lambda one: one.name)
def test_each_shipped_polymerase_gets_nebs_reaction_and_program(polymerase: Polymerase) -> None:
    lines, program, two_step = SHIPPED_PCR[polymerase.name]
    table = pcr_reaction(polymerase)
    assert [(one.name, one.volume_ul, one.stock, one.final) for one in table.components] == lines
    cycled = pcr_program(polymerase, annealing_temperature=55.0, amplicon_length=1500)
    assert [stage.cycles for stage in cycled.stages] == [1, 30, 1, 1]
    assert [
        (one.label, one.temperature_c, one.seconds)
        for stage in cycled.stages
        for one in stage.incubations
    ] == program
    for annealing, steps in ((round(two_step - 0.1, 1), 3), (two_step, 2)):
        split = pcr_program(polymerase, annealing_temperature=annealing, amplicon_length=1500)
        assert len(split.stages[1].incubations) == steps


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


def test_the_m13_pair_gives_nebs_empty_vector_band(puc19: SequenceRecord) -> None:
    assert amplicon_sizes(M13_FORWARD, M13_REVERSE, puc19) == (137,)


def test_the_product_band_is_the_empty_band_less_what_the_insert_replaced(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    # 137 bp empty, less the 57 bp multiple cloning site, plus 717 bp of GFP.
    assert bands(check, "Correct clone") == (797,)
    assert bands(check, "Empty vector") == (137,)


def test_two_flanking_primers_cannot_tell_the_insert_round_the_other_way(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    assert bands(check, "Reversed insert") == bands(check, "Correct clone")
    assert not check.tells_orientation


def test_a_junction_primer_tells_orientation_when_the_flanks_differ(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(
        product,
        junctions,
        vector=puc19,
        primers=(M13_FORWARD, M13_REVERSE),
        insert_primer=True,
    )
    assert len(check.primers) == 3
    assert bands(check, "Reversed insert") != bands(check, "Correct clone")
    assert check.tells_orientation
    assert bands(check, "Empty vector") == (137,)


def test_a_designed_pair_flanks_both_junctions(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, flank=60)
    assert bands(check, "Correct clone") == (837,)
    assert bands(check, "Empty vector") == (177,)
    assert all(report.status != "fail" for report in check.reports)


def test_the_gel_carries_one_lane_per_clone_and_a_ladder_for_the_range(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    gel = check.gel
    assert [lane.label for lane in gel.lanes] == [
        "Correct clone",
        "Empty vector",
        "Reversed insert",
    ]
    assert gel.ladder is LADDER_100_BP
    assert check.agarose_percent == 2.0


def test_a_band_over_a_kilobase_takes_the_wider_ladder() -> None:
    assert choose_ladder((137, 797)) is LADDER_100_BP
    assert choose_ladder((137, 3000)) is LADDER_1_KB_PLUS
    assert agarose_percent((137, 797)) == 2.0
    assert agarose_percent((137, 3000)) == 1.0


def test_the_annealing_temperature_and_extension_come_from_the_polymerase(
    product: SequenceRecord, puc19: SequenceRecord, junctions: tuple[int, int]
) -> None:
    check = colony_pcr_check(product, junctions, vector=puc19, primers=(M13_FORWARD, M13_REVERSE))
    # OneTaq: the 23-mer M13/pUC pair anneals at 51.5 C, and 797 bp extends for a minute.
    assert check.annealing_temperature == pytest.approx(51.5, abs=0.1)
    assert check.extension_seconds == 60


def test_sanger_primers_read_from_outside_each_junction(
    product: SequenceRecord, junctions: tuple[int, int], gfp: SequenceRecord
) -> None:
    forward, reverse = sanger_primers(product, junctions)
    assert forward.distance_bp >= 100
    assert reverse.distance_bp >= 100
    assert forward.read_bp == forward.distance_bp + len(gfp)
    assert reverse.read_bp == reverse.distance_bp + len(gfp)


def test_colony_pcr_check_refuses_a_junction_pair_it_cannot_place(
    product: SequenceRecord, puc19: SequenceRecord
) -> None:
    with pytest.raises(ValueError, match="junction"):
        colony_pcr_check(product, (395,), vector=puc19)


def test_a_junction_primer_is_designed_for_every_insert(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    check = colony_pcr_check(product, three_junctions, vector=puc19, flank=60, insert_primer=True)
    assert len(check.primers) == 4
    assert [clone.name for clone in check.clones] == [
        "Correct clone",
        "Empty vector",
        "Reversed insert 1",
        "Reversed insert 2",
    ]


def test_the_correct_clone_shows_one_band_reading_each_junction(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    check = colony_pcr_check(product, three_junctions, vector=puc19, flank=60, insert_primer=True)
    # 60 bases of vector then 100 into each insert, and the flanking pair across the whole span.
    assert bands(check, "Correct clone") == (160, 510, 837)
    assert bands(check, "Empty vector") == (177,)


def test_each_insert_gets_its_own_reversed_lane(
    product: SequenceRecord, puc19: SequenceRecord, three_junctions: tuple[int, ...]
) -> None:
    check = colony_pcr_check(product, three_junctions, vector=puc19, flank=60, insert_primer=True)
    correct = bands(check, "Correct clone")
    assert bands(check, "Reversed insert 1") != correct
    assert bands(check, "Reversed insert 2") != correct
    assert check.tells_orientation


def test_sanger_primers_read_from_outside_the_whole_inserted_span(
    product: SequenceRecord, three_junctions: tuple[int, ...], gfp: SequenceRecord
) -> None:
    forward, reverse = sanger_primers(product, three_junctions)
    assert forward.read_bp == forward.distance_bp + len(gfp)
    assert reverse.read_bp == reverse.distance_bp + len(gfp)
