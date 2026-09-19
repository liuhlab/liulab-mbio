"""The reactions this method runs, and the numbers the note states for them.

Every number is NEB's, through ``docs/research/restriction-ligation.md``. `test_plan.py` is
where the protocol these reach is checked.
"""

from liulab_mbio.cloning.restriction.bench import (
    BLUNT_SECONDS,
    COHESIVE_SECONDS,
    LIGASE_KILL_CELSIUS,
    LIGASE_KILL_SECONDS,
    ROOM_CELSIUS,
    digest_amount,
    digest_reaction,
    ligation_amounts,
    ligation_program,
    shared_buffer,
)
from liulab_mbio.enzymes import get_enzyme


def test_a_double_digest_is_nebs_typical_reaction_in_the_buffer_both_enzymes_share():
    pair = [get_enzyme("EcoRI"), get_enzyme("BamHI")]
    reaction = digest_reaction(pair, digest_amount(("pUC19", 2686)))
    table = {component.name: component for component in reaction.components}
    # NEB's typical digest: 10 units of each enzyme, 5 µL of 10X buffer, 1 µg of DNA, 50 µL.
    assert table["EcoRI-HF (R3101)"].final == "10 units"
    assert table["10X rCutSmart Buffer"].volume_ul == 5.0
    assert table["pUC19"].final.endswith("(1000 ng)")
    assert sum(one.volume_ul for one in reaction.components) == 50.0


def test_the_buffer_is_named_only_where_both_enzymes_are_supplied_in_one():
    assert shared_buffer([get_enzyme("EcoRI"), get_enzyme("BamHI")]) == "rCutSmart Buffer"
    assert shared_buffer([get_enzyme("EcoRI"), get_enzyme("BsmBI")]) is None


def test_the_ligation_takes_its_ratio_in_picomoles_so_the_masses_fall_out_of_the_lengths():
    # The 1:3 in NEB's table is an example for a 4 kb vector and a 1 kb insert, so the ratio is
    # taken in picomoles and the masses fall out of the fragments' own lengths.
    backbone, insert = ligation_amounts(("pUC19 backbone", 2665), ("pTrc-GFP insert", 723))
    assert (backbone.pmol, insert.pmol) == (0.02, 0.06)
    assert insert.nanograms < backbone.nanograms


def test_a_blunt_ligation_is_held_longer_than_a_cohesive_one_and_stopped_the_same_way():
    cohesive, blunt = (ligation_program(blunt=one) for one in (False, True))
    for program, seconds in ((cohesive, COHESIVE_SECONDS), (blunt, BLUNT_SECONDS)):
        hold, kill = (one for stage in program.stages for one in stage.incubations)
        assert (hold.temperature_c, hold.seconds) == (ROOM_CELSIUS, seconds)
        assert (kill.temperature_c, kill.seconds) == (LIGASE_KILL_CELSIUS, LIGASE_KILL_SECONDS)
