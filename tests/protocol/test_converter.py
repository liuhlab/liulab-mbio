"""A protocol's figures, gels and plates handed to molviz, held to its own checks of each format."""

import json
from pathlib import Path

import molviz

from mbio.protocol.converter import molviz_figure, molviz_gel, molviz_plate, molviz_stamp
from mbio.protocol.model import (
    Figure,
    Gel,
    Ladder,
    Lane,
    Move,
    Plate,
    Protocol,
    Transfer,
    Well,
    read_project,
    read_protocol,
)

EXAMPLES = Path(__file__).parents[2] / "docs" / "examples"


def _protocols() -> list[tuple[Path, Protocol]]:
    """Every protocol the docs show, each with the folder its figures' records are named from."""
    found = []
    for path in sorted(EXAMPLES.rglob("*.json")):
        data = json.loads(path.read_text(encoding="utf-8"))
        if "protocols" in data:
            found += [(path.parent, one) for one in read_project(path).protocols]
        elif "steps" in data:
            found.append((path.parent, read_protocol(path)))
    return found


def test_every_figure_gel_and_plate_the_docs_show_converts_to_one_molviz_accepts() -> None:
    protocols = _protocols()
    assert len(protocols) > 5
    figures = gels = plates = 0
    for base, protocol in protocols:
        for step in protocol.steps:
            for figure in step.figures:
                converted = molviz_figure(figure, base)
                for record in converted["records"]:
                    molviz.check(record)
                molviz.draw_figure(converted["records"], **converted["options"])
                figures += 1
            for gel in step.gels:
                molviz.check_gel(molviz_gel(gel))
                gels += 1
            for transfer in step.transfers:
                for plate in molviz_stamp(transfer, protocol.plates):
                    molviz.check_plate(plate)
                    plates += 1
        for plate in protocol.plates:
            molviz.check_plate(molviz_plate(plate))
            plates += 1
    assert figures
    assert gels
    assert plates


def test_a_figure_names_its_records_and_says_what_differs_from_molvizs_defaults(
    data_dir: Path,
) -> None:
    figure = Figure(
        ("pUC19.dna", "GFP.dna"),
        "pUC19, then GFP",
        span=(390, 460),
        linear=True,
        sequence_view=True,
        enzymes=("EcoRI",),
        highlight=("AmpR",),
    )
    converted = molviz_figure(figure, data_dir)
    assert [record["name"] for record in converted["records"]] == ["pUC19", "GFP"]
    assert converted["options"] == {
        "span": {"start": 390, "end": 460},
        "map": "line",
        "sequence_view": True,
        "enzymes": ["EcoRI"],
        "highlight": ["AmpR"],
    }
    assert molviz_figure(Figure(("pUC19.dna",), "pUC19"), data_dir)["options"] == {}


def test_a_gel_names_its_ladder_and_lanes_with_their_band_sizes() -> None:
    gel = Gel(Ladder("100 bp", (100, 200)), (Lane("GFP", (741,)), Lane("water")), title="PCR")
    assert molviz_gel(gel) == {
        "title": "PCR",
        "ladder": {"name": "100 bp", "bands": [100, 200]},
        "lanes": [{"name": "GFP", "bands": [741]}, {"name": "water", "bands": []}],
    }


def test_a_plate_holds_its_labels_and_seating_in_reading_order_the_seating_winning() -> None:
    plate = Plate(
        "picked", 96, seating={"B1": "clone 2", "A2": "clone 1"}, labels={"A1": "q1", "A2": "q2"}
    )
    converted = molviz_plate(plate)
    assert converted == {
        "title": "picked",
        "rows": 8,
        "columns": 12,
        "wells": {"A1": "q1", "A2": "clone 1", "B1": "clone 2"},
    }
    assert list(converted["wells"]) == ["A1", "A2", "B1"]


def test_a_stamp_is_its_two_plates_every_touched_well_holding_the_transfers_title() -> None:
    plates = (Plate("picked", 96), Plate("index", 384))
    # Every other well of the first row, into the first wells of the next plate's.
    moves = tuple(
        Move(Well("picked", f"A{2 * n - 1}"), Well("index", f"A{n}"), 1.0) for n in (1, 2)
    )
    from_, into = molviz_stamp(Transfer("Index the clones", moves), plates)
    assert from_ == {
        "title": "From picked",
        "rows": 8,
        "columns": 12,
        "wells": {"A1": "Index the clones", "A3": "Index the clones"},
    }
    assert (into["title"], into["rows"], list(into["wells"])) == ("Into index", 16, ["A1", "A2"])
    # Moves at two volumes repeat no one pattern, so they are no stamp and draw no plate.
    uneven = (moves[0], Move(moves[1].source, moves[1].destination, 2.0))
    assert molviz_stamp(Transfer("Index the clones", uneven), plates) == []
