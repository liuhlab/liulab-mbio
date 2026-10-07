"""The reactions and the programs one round runs, and which digest each blunt enzyme belongs to."""

from liulab_mbio.bench.amounts import dna_amount
from liulab_mbio.enzymes import get_enzyme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord
from liulab_synbio.igga import stages
from liulab_synbio.igga.bench import (
    DIGEST_CELSIUS,
    DIGEST_SECONDS,
    DIGEST_VOLUME_UL,
    ENZYME_UL,
    GROWTH_CELSIUS,
    LIGATION_VOLUME_UL,
    OUTGROWTH_SECONDS,
    RECOVERY_SECONDS,
)
from liulab_synbio.igga.method import IGGA
from liulab_synbio.igga.steps import (
    choppers,
    digest_program,
    digest_reaction,
    growth_program,
    ligation_reaction,
)


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


def test_each_blunt_enzyme_belongs_to_the_digest_whose_piece_it_cuts():
    inside, outside = choppers(IGGA)

    assert [one.name for one in inside] == ["SrfI"]
    assert [one.name for one in outside] == ["PmeI"]


def marked(marker: str | None) -> SequenceRecord:
    """A vector annotating `marker` as its selection marker, or nothing where it is `None`."""
    features = (Feature(marker, "CDS", (Segment(10, 100),)),) if marker else ()
    return SequenceRecord("TA" * 200, topology="circular", name="vector", features=features)


def test_a_round_is_plated_on_the_vectors_own_marker_and_not_the_papers():
    """The method rebuilt its DMX vector KanR, so the published carbenicillin does not carry."""
    assert stages.selection_for(marked("KanR")) == "50 µg/mL kanamycin"
    assert stages.selection_for(marked("AmpR")) == "100 µg/mL carbenicillin"


def test_a_vector_naming_no_marker_names_no_drug_and_keeps_the_hole():
    """H22 stands only where nothing in the record answers it."""
    assert stages.selection_for(marked(None)) == ""
    assert stages.ROUND_SELECTION in stages.holes_for(marked(None))
    assert stages.ROUND_SELECTION not in stages.holes_for(marked("KanR"))
    assert next(one.id for one in stages.holes_for(marked(None))) == "H22"
