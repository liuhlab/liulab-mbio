"""What a protocol declares, and what a project chains by name."""

from pathlib import Path

import pytest

from liulab_mbio.cloning.plan import as_project
from liulab_mbio.protocol.model import (
    Bill,
    BillRow,
    Check,
    Item,
    Project,
    Protocol,
    Step,
    read_project,
    write_project,
)

ENTRY = Item("entry clone", "pENTR carrying the insert", spec=("≥100 ng/µL",), storage="-20 °C")


def bp(**kwargs: object) -> Protocol:
    """The protocol that makes the entry clone, with fields changed."""
    return Protocol("BP reaction", produces=(ENTRY,), **kwargs)  # pyright: ignore[reportArgumentType]


def lr(**kwargs: object) -> Protocol:
    """The protocol that spends the entry clone, with fields changed."""
    return Protocol(
        "LR reaction",
        consumes=(ENTRY,),
        produces=(Item("expression clone", "the finished plasmid"),),
        **kwargs,  # pyright: ignore[reportArgumentType]
    )


def test_an_item_carries_the_bench_prose_beside_the_name_that_is_the_contract() -> None:
    assert ENTRY.name == "entry clone"
    assert ENTRY.spec == ("≥100 ng/µL",)
    assert ENTRY.storage == "-20 °C"


@pytest.mark.parametrize("build", [lambda: Item("", "x"), lambda: Item("x", " ")])
def test_an_item_that_names_nothing_or_says_nothing_is_refused(build) -> None:
    with pytest.raises(ValueError):  # noqa: PT011
        build()


def test_a_project_chains_one_protocols_output_into_the_next() -> None:
    one = Project("Gateway cloning", protocols=(bp(), lr()))
    assert one.audit() == (Check("handoffs", "pass", "1 consumed item resolves"),)


def test_what_the_project_was_handed_satisfies_a_consumer() -> None:
    one = Project("Gateway cloning", inputs=(ENTRY,), protocols=(lr(),))
    assert one.audit()[0].status == "pass"


def test_a_consumed_name_nothing_produces_is_a_badge_and_not_an_exception() -> None:
    (check,) = Project("Gateway cloning", protocols=(lr(),)).audit()
    assert check.status == "fail"
    assert check.detail == "LR reaction consumes 'entry clone', which nothing hands it"


def test_a_later_protocol_does_not_feed_an_earlier_one() -> None:
    (check,) = Project("Gateway cloning", protocols=(lr(), bp())).audit()
    assert check.status == "fail"


def test_a_step_carries_a_section_label_and_the_steps_stay_one_flat_list() -> None:
    one = Protocol("Demo", steps=(Step("Thaw the cells", section="Day 1"), Step("Plate them")))
    assert [step.section for step in one.steps] == ["Day 1", ""]


def test_a_project_reads_back_equal_from_its_json(tmp_path: Path) -> None:
    one = Project(
        "Gateway cloning",
        summary="Two reactions, two days apart.",
        inputs=(Item("insert", "the amplicon, attB-tailed"),),
        protocols=(bp(), lr(steps=(Step("Set up the LR reaction", section="Day 2"),))),
        checks=(Check("junctions", "pass", "4 of 4 read in frame"),),
        bill=Bill((BillRow("LR Clonase II", 1, unit="reaction"),)),
    )
    assert read_project(write_project(one, tmp_path / "project.json")) == one


def test_one_protocol_wraps_as_a_project_of_one() -> None:
    one = as_project(lr(summary="Move the insert into the destination vector."))
    assert one.title == "LR reaction"
    assert one.summary == "Move the insert into the destination vector."
    assert one.inputs == (ENTRY,)
    assert one.audit()[0].status == "pass"
