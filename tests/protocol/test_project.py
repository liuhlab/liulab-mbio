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
    Topic,
    read_project,
    write_project,
)

ENTRY = Item("entry clone", "pENTR carrying the insert", spec=("≥100 ng/µL",), storage="-20 °C")
PLATE = Item("archive plate", "the picked colonies")
CALLS = Item("well calls", "one design name per well")
POOL = Item("indexed pool", "every well's amplicon, pooled")
JOB = "read every well back"


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


def way(title: str, *produces: Item) -> Protocol:
    """One way of reading every well back, leaving whatever it is given to leave."""
    return Protocol(title, choice=JOB, consumes=(PLATE,), produces=produces or (CALLS,))


def readback(*protocols: Protocol, explained: bool = True) -> Project:
    """A run whose middle place is a choice, with the guidance the reader needs or without it."""
    return Project(
        "DMX",
        background=(Topic(JOB, ("Ligation scales; index PCR is quicker to set up.",)),)
        if explained
        else (),
        inputs=(PLATE,),
        protocols=protocols,
    )


def named(project: Project, name: str) -> Check:
    """The one check of `project` called `name`."""
    (check,) = [one for one in project.audit() if one.name == name]
    return check


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
    assert one.audit()[0] == Check("handoffs", "pass", "1 consumed item resolves")


def test_what_the_project_was_handed_satisfies_a_consumer() -> None:
    one = Project("Gateway cloning", inputs=(ENTRY,), protocols=(lr(),))
    assert one.audit()[0].status == "pass"


def test_a_consumed_name_nothing_produces_is_a_badge_and_not_an_exception() -> None:
    check = Project("Gateway cloning", protocols=(lr(),)).audit()[0]
    assert check.status == "fail"
    assert check.detail == "LR reaction consumes 'entry clone', which nothing hands it"


def test_a_later_protocol_does_not_feed_an_earlier_one() -> None:
    check = Project("Gateway cloning", protocols=(lr(), bp())).audit()[0]
    assert check.status == "fail"


def test_ways_of_one_job_with_another_protocol_between_them_are_refused() -> None:
    """The ways take one place in the run, and a place is contiguous."""
    with pytest.raises(ValueError, match="stand apart"):
        readback(way("Barcode ligation"), Protocol("Pick"), way("Index PCR"))


def test_a_way_is_judged_against_what_was_handed_over_before_the_choice() -> None:
    """A sibling hands a sibling nothing: the bench did one of them, not both."""
    sibling = Protocol("Index PCR", choice=JOB, consumes=(CALLS,), produces=(CALLS,))
    check = named(readback(way("Barcode ligation"), sibling), "handoffs")
    assert check.status == "fail"
    assert check.detail == "Index PCR consumes 'well calls', which nothing hands it"


def test_only_a_name_every_way_leaves_passes_out_of_the_choice() -> None:
    ways = (way("Barcode ligation"), way("Index PCR", CALLS, POOL))
    assert (
        named(readback(*ways, Protocol("Report", consumes=(CALLS,))), "handoffs").status == "pass"
    )
    check = named(readback(*ways, Protocol("Report", consumes=(POOL,))), "handoffs")
    assert check.status == "fail"
    assert check.detail == "Report consumes 'indexed pool', which nothing hands it"


def test_a_run_offering_a_choice_says_how_many_ways_it_has_and_what_they_leave() -> None:
    check = named(readback(way("Barcode ligation"), way("Index PCR")), "choices")
    assert check.status == "pass"
    assert check.detail == "2 ways to read every well back, each leaving the same thing."


def test_a_job_only_one_protocol_names_fails_the_choices_check() -> None:
    check = named(readback(way("Barcode ligation")), "choices")
    assert check.status == "fail"
    assert (
        check.detail
        == "Only one way to read every well back is written, so there is nothing to choose."
    )


def test_ways_leaving_the_bench_holding_different_things_fail_the_choices_check() -> None:
    check = named(readback(way("Barcode ligation"), way("Index PCR", CALLS, POOL)), "choices")
    assert check.status == "fail"
    assert check.detail == (
        'The ways to read every well back do not leave the same things: "indexed pool" comes '
        "from only one of them."
    )


def test_a_choice_the_overview_says_nothing_about_warns_rather_than_fails() -> None:
    """The run is still followable; a reader met with two ways and no help picking is not."""
    ways = (way("Barcode ligation"), way("Index PCR"))
    check = named(readback(*ways, explained=False), "choices")
    assert check.status == "warn"
    assert check.detail == "Nothing on the overview says how to pick a way to read every well back."


def test_a_step_carries_a_section_label_and_the_steps_stay_one_flat_list() -> None:
    one = Protocol("Demo", steps=(Step("Thaw the cells", section="Day 1"), Step("Plate them")))
    assert [step.section for step in one.steps] == ["Day 1", ""]


def test_a_project_carries_the_background_no_one_protocol_explains() -> None:
    one = Project(
        "Gateway cloning",
        background=(Topic("Why two reactions", ("The entry clone is reused.",)),),
    )
    assert one.background[0].body == ("The entry clone is reused.",)
    with pytest.raises(ValueError, match="title"):
        Topic("", ("x",))


def test_a_project_reads_back_equal_from_its_json(tmp_path: Path) -> None:
    one = Project(
        "Gateway cloning",
        summary="Two reactions, two days apart.",
        background=(Topic("Why two reactions", ("The entry clone is reused.",)),),
        inputs=(Item("insert", "the amplicon, attB-tailed"),),
        protocols=(bp(), lr(steps=(Step("Set up the LR reaction", section="Day 2"),))),
        checks=(Check("junctions", "pass", "4 of 4 read in frame"),),
        bill=Bill((BillRow("LR Clonase II", 1, unit="reaction", charge="129.00"),)),
    )
    assert read_project(write_project(one, tmp_path / "project.json")) == one


def test_one_protocol_wraps_as_a_project_of_one() -> None:
    one = as_project(lr(summary="Move the insert into the destination vector."))
    assert one.title == "LR reaction"
    assert one.summary == "Move the insert into the destination vector."
    assert one.inputs == (ENTRY,)
    assert one.audit()[0].status == "pass"
