"""What a protocol carries beyond its steps: a number's source, a rule, a well and a hole."""

from pathlib import Path

import pytest

from liulab_mbio.bench import materials, plates
from liulab_mbio.protocol import model
from liulab_mbio.protocol.model import (
    Check,
    Citation,
    Component,
    Hole,
    Incubation,
    Material,
    Move,
    Oligo,
    Protocol,
    ReactionTable,
    Source,
    Stage,
    Stamp,
    Step,
    ThermocyclerProgram,
    Transfer,
    Troubleshooting,
    Vessel,
    Well,
    citing,
    read_protocol,
    write_protocol,
)
from liulab_mbio.protocol.render import NO_NUMBER, render_html

from ..html import parse

LIGASE = materials.material(
    "T7 DNA Ligase", catalog="#M0318L", citation=Citation("M0318", "reaction conditions")
)
BUFFER = materials.material("StickTogether DNA Ligase Buffer", catalog="#B0535S")
FRESH_BUFFER = Troubleshooting(
    "Few colonies",
    "Repeat the ligation with fresh buffer.",
    citation=Citation("NEB cloning", "ligation"),
)


def ligation(*, extra: tuple[Component, ...] = (), title: str = "Ligate") -> Step:
    return Step(
        title,
        instructions=("Mix on ice, then hold 30 minutes at room temperature.",),
        tables=(
            ReactionTable(
                (
                    Component("T7 DNA Ligase", 1.0, citation=Citation("M0318", "PS v2.0")),
                    Component("StickTogether DNA Ligase Buffer", 100.0),
                    *extra,
                )
            ),
        ),
    )


def protocol(*steps: Step, **rest: object) -> Protocol:
    return Protocol(
        "iGGA round",
        materials=(LIGASE, BUFFER),
        steps=steps,
        sources=dict(materials.SOURCES),
        **rest,  # pyright: ignore[reportArgumentType]
    )


def test_both_t7_ligase_rules_reach_the_rendered_page() -> None:
    page = render_html(protocol(ligation()))
    assert "Never add PEG to a T7 ligase reaction" in page
    assert "Do not heat-inactivate T7 DNA ligase in PEG" in page


def test_the_rule_travels_with_the_material_and_no_step_stores_it() -> None:
    step = ligation()
    assert step.cautions == ()
    assert [rule.subject for _, rule in protocol(step).rules_for(step)] == [
        "PEG",
        "heat inactivation",
    ]


def test_the_caution_travels_with_the_polymerase_and_no_step_stores_it() -> None:
    tube = materials.material("Q5 DNA Polymerase", catalog="M0491")
    amplifying = Step("Amplify", tables=(ReactionTable((Component("Q5 DNA Polymerase", 0.5),)),))
    one = Protocol("PCR", materials=(tube,), steps=(amplifying,))
    assert amplifying.cautions == ()
    assert one.cautions_for(amplifying) == (materials.POLYMERASE_ON_ICE,)
    assert materials.POLYMERASE_ON_ICE in render_html(one)


def test_a_step_writing_out_a_caution_its_material_carries_shows_it_once() -> None:
    tube = materials.material("Q5 DNA Polymerase", catalog="M0491")
    written = Step(
        "Amplify",
        cautions=(materials.POLYMERASE_ON_ICE, "Spin the plate down."),
        tables=(ReactionTable((Component("Q5 DNA Polymerase", 0.5),)),),
    )
    one = Protocol("PCR", materials=(tube,), steps=(written,))
    assert one.cautions_for(written) == (materials.POLYMERASE_ON_ICE, "Spin the plate down.")


def test_adding_peg_to_the_t7_reaction_fails_the_protocols_own_check() -> None:
    broken = protocol(ligation(extra=(Component("PEG 6000", 5.0),)))
    (check,) = [one for one in broken.audit() if one.name == "rules"]
    assert check.status == "fail"
    assert "forbids PEG" in check.detail


def test_heat_inactivating_in_peg_fails_and_the_same_step_without_peg_does_not() -> None:
    program = ThermocyclerProgram((Stage((Incubation("Heat inactivation", 65.0, 600),)),))
    inactivating = Step(
        "Stop the reaction",
        tables=ligation().tables,
        programs=(program,),
    )
    assert any(
        one.name == "rules" and one.status == "fail" for one in protocol(inactivating).audit()
    )
    # The same program on a reaction with no PEG in it is nobody's business: `when` gates it.
    dry = Protocol(
        "dry",
        materials=(materials.material("T7 DNA Ligase", catalog="#M0318L"),),
        steps=(
            Step(
                "Stop the reaction",
                tables=(ReactionTable((Component("T7 DNA Ligase", 1.0),)),),
                programs=(program,),
            ),
        ),
    )
    (check,) = [one for one in dry.audit() if one.name == "rules"]
    assert check.status == "pass"


def test_a_step_references_a_well_rather_than_describing_one() -> None:
    picked = plates.plate("picked", 384, seating=plates.seat(["UMI-1"], 384))
    compressed = plates.plate("lysate", 1536)
    moved = plates.compact(
        [Well("picked", "A1")], compressed, 2.0, title="Compress", instrument="ECHO 525"
    )
    one = protocol(
        Step("Barcode in lysate", transfers=(moved,)),
        plates=(picked, compressed),
        oligos=(Oligo("UMI-1", "ACGTACGTAC"),),
    )
    page = render_html(one)
    assert "picked A1" in page
    assert "lysate A1" in page
    assert "ECHO 525" in page
    assert [check.status for check in one.audit() if check.name == "wells"] == ["pass"]


def test_a_well_naming_nothing_the_protocol_holds_is_reported() -> None:
    one = protocol(ligation(), plates=(plates.plate("picked", 96, seating={"A1": "ghost"}),))
    (check,) = [c for c in one.audit() if c.name == "wells"]
    assert check.status == "fail"
    assert "ghost" in check.detail


def test_a_hole_renders_and_never_reads_as_a_value() -> None:
    hole = Hole(
        "H23",
        "no document gives the ligase units",
        "unpublished",
        where="ligase units per reaction",
        issue="liuhlab/liulab-mbio#264",
    )
    page = render_html(protocol(Step("Ligate", holes=(hole,))))
    assert NO_NUMBER in page
    assert "H23" in page
    assert "nobody published it" in page
    assert "liuhlab/liulab-mbio#264" not in page
    assert (
        "1 number has no source"
        in [check.detail for check in protocol(Step("Ligate", holes=(hole,))).audit()][3]
    )


def test_a_cycle_count_is_cited_like_any_number_and_a_blank_one_prints_as_none() -> None:
    cited = Stage((Incubation("Denature", 98, 10),), cycles=12, citation=Citation("nowhere"))
    blank = Stage((Incubation("Denature", 98, 10),), cycles=None)
    one = protocol(Step("PCR", programs=(ThermocyclerProgram((cited, blank)),)))
    (check,) = [c for c in one.audit() if c.name == "sources"]
    assert check.status == "fail"
    assert "nowhere" in check.detail
    assert NO_NUMBER in render_html(one)


def test_a_protocol_keeps_the_sources_it_cites_and_drops_the_rest() -> None:
    """A reference list a reader cannot follow back to a row on the page is what this stops."""
    catalogue = {"M0318": Source("NEB #M0318"), "C3020": Source("NEB #C3020")}
    one = citing(Protocol("x", materials=(LIGASE,), sources=catalogue))
    assert list(one.sources) == ["M0318"]
    assert [c.status for c in one.audit() if c.name == "sources"] == ["pass"]


def test_a_hole_the_run_and_a_step_both_carry_is_collected_once() -> None:
    hole = Hole("H24", "no mass for the one-pot assembly", "undecided")
    one = protocol(Step("Assemble", holes=(hole,)), holes=(hole,))
    assert [each.id for each in one.all_holes] == ["H24"]
    detail = next(c.detail for c in one.audit() if c.name == "holes")
    assert "1 number has no source" in detail
    assert "1 number this protocol would otherwise have to invent" in render_html(one)


def test_the_banner_counts_bench_numbers_apart_from_prices() -> None:
    """A missing price is a missing input, so it never reads as an unsourced bench number."""
    one = protocol(
        Step(
            "Assemble",
            holes=(
                Hole("H24", "no mass is sourced", "undecided"),
                Hole("H31", "no working vector is named", "lab"),
            ),
        ),
        holes=(Hole("P1", "nothing prices the kit", "price"),),
    )
    page = render_html(one)
    assert "<strong>3</strong> numbers in this protocol have no source" in page
    assert "2 bench numbers and 1 price." in page
    assert "ready to run" not in page


def holes_verdict(*holes: Hole) -> Check:
    (check,) = [c for c in protocol(Step("Assemble", holes=holes)).audit() if c.name == "holes"]
    return check


def test_a_protocol_missing_no_number_passes_its_own_hole_check() -> None:
    assert holes_verdict().status == "pass"


def test_a_hole_no_source_closes_leaves_the_verdict_open_rather_than_failing() -> None:
    """ADR 0017: a plan is finished with these still standing, so they are no failure."""
    standing = tuple(
        Hole(f"H{n}", "nothing sources it", kind)
        for n, kind in enumerate(("undecided", "unpublished", "lab", "price"), 1)
    )
    assert holes_verdict(*standing).status is None


def test_a_hole_waiting_on_a_source_nobody_read_fails_and_is_named() -> None:
    """ADR 0017: a hole naming a source nobody read is work left undone, not a finished plan."""
    check = holes_verdict(
        Hole("H24", "nothing sources it", "undecided"),
        Hole("IDX2", "the units one reaction takes", "unread"),
    )
    assert check.status == "fail"
    assert "IDX2" in check.detail
    assert "H24" not in check.detail


def test_a_price_hole_names_no_issue() -> None:
    with pytest.raises(ValueError, match="missing price is a missing input"):
        Hole("P1", "nothing prices it", "price", issue="#264")


def test_a_citation_naming_no_source_is_reported() -> None:
    one = Protocol("x", materials=(Material("Water", citation=Citation("nowhere")),))
    (check,) = [c for c in one.audit() if c.name == "sources"]
    assert check.status == "fail"
    assert "nowhere" in check.detail


def test_a_troubleshooting_answer_cites_its_source_beside_it() -> None:
    one = citing(
        Protocol(
            "x",
            steps=(Step("Plate", troubleshooting=(FRESH_BUFFER,)),),
            sources={"NEB cloning": Source("NEB cloning troubleshooting guide")},
        )
    )
    assert list(one.sources) == ["NEB cloning"]
    (answer,) = parse(render_html(one)).find_all(cls="trouble")[0].find_all("dd")
    assert [a.attrs["href"] for a in answer.find_all("a", cls="cite")] == ["#source-neb-cloning"]


def test_a_troubleshooting_citation_naming_no_source_is_reported() -> None:
    one = Protocol("x", steps=(Step("Plate", troubleshooting=(FRESH_BUFFER,)),))
    (check,) = [c for c in one.audit() if c.name == "sources"]
    assert check.status == "fail"
    assert "NEB cloning" in check.detail


def test_everything_new_round_trips_through_json(tmp_path: Path) -> None:
    one = protocol(
        ligation(),
        Step("Plate", troubleshooting=(FRESH_BUFFER,)),
        plates=(plates.plate("picked", 96, seating={"A1": "T7 DNA Ligase"}),),
        vessels=(Vessel("reservoir", kind="trough"),),
        holes=(Hole("H1", "no polymerase is named", "undecided", issue="liuhlab/liulab-mbio#264"),),
    )
    path = write_protocol(one, tmp_path / "protocol.json")
    assert read_protocol(path) == one


def test_the_row_label_follows_the_plate_past_the_alphabet() -> None:
    assert [model.row_label(n) for n in (0, 25, 26, 31)] == ["A", "Z", "AA", "AF"]


def test_a_labelled_well_describes_itself_and_resolves_against_nothing() -> None:
    """A label is what `holds` is for a plate, said a well at a time, so it never dangles."""
    one = protocol(
        ligation(),
        plates=(plates.plate("picked", 96, labels={"A1": "quarter 1", "A2": "quarter 2"}),),
    )

    (check,) = [c for c in one.audit() if c.name == "wells"]

    assert check.status == "pass"
    assert "quarter 1" in render_html(one)


def test_a_label_does_not_excuse_a_seated_well_naming_nothing_declared() -> None:
    """The one defect the check catches survives the new channel."""
    one = protocol(
        ligation(),
        plates=(plates.plate("picked", 96, seating={"A1": "ghost"}, labels={"A2": "quarter 1"}),),
    )

    (check,) = [c for c in one.audit() if c.name == "wells"]

    assert check.status == "fail"
    assert "ghost" in check.detail


def test_a_label_off_the_array_is_refused_as_a_seating_is() -> None:
    with pytest.raises(ValueError, match="has no well"):
        plates.plate("picked", 96, labels={"Z1": "quarter 1"})


def test_a_transfer_that_repeats_one_pattern_answers_with_the_pattern() -> None:
    """96 moves of a quarter-sampling are one stride and one start, so a page can draw them."""
    moved = Transfer(
        "Sample a quarter",
        tuple(
            Move(Well("picked", at), Well("index", to), 1.0)
            for at, to in (("A2", "A1"), ("A4", "A2"), ("C2", "B1"))
        ),
    )

    stamp = moved.stamp

    assert stamp is not None
    assert stamp == Stamp(2, 0, 1)
    assert stamp.words == "every other row and column, starting A2"


def test_moves_that_repeat_no_one_pattern_stay_a_list() -> None:
    """A compaction leaves out the wells that failed, so no stride describes it."""
    pooled = plates.plate("pooled", 96)
    moved = plates.compact(
        [Well("picked", "A1"), Well("picked", "B4")], pooled, 2.0, title="Compact"
    )

    assert moved.stamp is None
