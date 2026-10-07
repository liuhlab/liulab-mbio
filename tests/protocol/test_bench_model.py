"""What a protocol carries beyond its steps: a number's source, a rule, a well and a hole."""

from pathlib import Path

import pytest

from liulab_mbio.bench import materials, plates
from liulab_mbio.protocol import model
from liulab_mbio.protocol.model import (
    Citation,
    Component,
    Hole,
    Incubation,
    Material,
    Oligo,
    Protocol,
    ReactionTable,
    Stage,
    Step,
    ThermocyclerProgram,
    Vessel,
    Well,
    read_protocol,
    write_protocol,
)
from liulab_mbio.protocol.render import NO_NUMBER, render_html

LIGASE = materials.material(
    "T7 DNA Ligase", catalog="#M0318L", citation=Citation("M0318", "reaction conditions")
)
BUFFER = materials.material("StickTogether DNA Ligase Buffer", catalog="#B0535S")


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
    assert (
        "1 number has no source"
        in [check.detail for check in protocol(Step("Ligate", holes=(hole,))).audit()][3]
    )


def test_a_price_hole_names_no_issue() -> None:
    with pytest.raises(ValueError, match="missing price is a missing input"):
        Hole("P1", "nothing prices it", "price", issue="#264")


def test_a_citation_naming_no_source_is_reported() -> None:
    one = Protocol("x", materials=(Material("Water", citation=Citation("nowhere")),))
    (check,) = [c for c in one.audit() if c.name == "sources"]
    assert check.status == "fail"
    assert "nowhere" in check.detail


def test_everything_new_round_trips_through_json(tmp_path: Path) -> None:
    one = protocol(
        ligation(),
        plates=(plates.plate("picked", 96, seating={"A1": "T7 DNA Ligase"}),),
        vessels=(Vessel("reservoir", kind="trough"),),
        holes=(Hole("H1", "no polymerase is named", "undecided"),),
    )
    path = write_protocol(one, tmp_path / "protocol.json")
    assert read_protocol(path) == one


def test_the_row_label_follows_the_plate_past_the_alphabet() -> None:
    assert [model.row_label(n) for n in (0, 25, 26, 31)] == ["A", "Z", "AA", "AF"]
