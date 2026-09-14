"""The reactions and the programs one round runs, and which digest each blunt enzyme belongs to."""

from pathlib import Path

import pytest

from liulab_mbio.bench.amounts import dna_amount
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.library.bench import (
    DIGEST_CELSIUS,
    DIGEST_SECONDS,
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    GROWTH_CELSIUS,
    LIGATION_VOLUME_UL,
    OUTGROWTH_SECONDS,
    RECOVERY_SECONDS,
)
from liulab_mbio.library.scheme import read_scheme
from liulab_mbio.library.steps import (
    choppers,
    digest_program,
    digest_reaction,
    growth_program,
    ligation_reaction,
)

EXAMPLE = Path(__file__).parents[2] / "docs" / "examples" / "protein-library" / "scheme.json"


@pytest.fixture(scope="module")
def scheme():
    return read_scheme(EXAMPLE)


def test_a_digest_fills_its_volume_and_gives_each_enzyme_its_own_line():
    dna = dna_amount("library", 5000, pmol=0.3)
    enzymes = (get_enzyme("BsaI"), get_enzyme("SrfI"))

    table = digest_reaction(dna, enzymes)

    assert round(sum(one.volume_ul for one in table.components), 2) == DIGEST_VOLUME_UL
    assert [one.volume_ul for one in table.components[1:3]] == [ENZYME_UL, ENZYME_UL]
    # The DNA goes in each tube, the rest into the mix.
    assert table.components[0].master_mix is False


def test_a_digest_runs_the_first_enzyme_alone_and_then_the_second():
    program = digest_program((get_enzyme("BsaI"), get_enzyme("SrfI")))

    assert len(program.stages) == 2
    assert program.duration_seconds == 2 * DIGEST_SECONDS
    for stage in program.stages:
        for incubation in stage.incubations:
            assert incubation.temperature_c == DIGEST_CELSIUS
    assert program.stages[0].incubations[0].label == "BsaI"
    assert "SrfI" in program.stages[1].incubations[0].label


def test_a_ligation_fills_its_volume_and_leaves_the_ligase_to_the_supplier():
    amounts = (
        dna_amount("library, opened", 5000, pmol=0.006),
        dna_amount("part, released", 400, pmol=0.006),
    )

    table = ligation_reaction(amounts)

    assert round(sum(one.volume_ul for one in table.components), 2) == LIGATION_VOLUME_UL
    assert [one.master_mix for one in table.components] == [False, False, True]
    assert "T7 DNA Ligase" in table.components[-1].name


def test_both_growth_steps_run_at_thirty_degrees():
    program = growth_program()

    temperatures = {
        incubation.temperature_c for stage in program.stages for incubation in stage.incubations
    }
    assert temperatures == {GROWTH_CELSIUS}
    assert GROWTH_CELSIUS == 30.0
    assert program.duration_seconds == RECOVERY_SECONDS + OUTGROWTH_SECONDS[0]


def test_each_blunt_enzyme_belongs_to_the_digest_whose_piece_it_cuts(scheme):
    inside, outside = choppers(scheme)

    assert [one.name for one in inside] == ["SrfI"]
    assert [one.name for one in outside] == ["PmeI"]
