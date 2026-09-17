"""A record drawn through `draw_map` and written as a page: what the page carries, and refusals."""

import base64
import dataclasses
import json
import math
import re
import struct
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pypdf
import pytest
import vl_convert
from pypdf.generic import DictionaryObject

from liulab_mbio.plot import (
    Drawing,
    circular,
    convert,
    draw_map,
    layers,
    linear,
    sequence_view,
    svg,
)
from liulab_mbio.plot.fonts import BOLD, MONO, SANS
from liulab_mbio.plot.labels import Box
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand

from ..html import Node, parse
from . import crowds


@pytest.fixture(scope="module")
def colour_test(data_dir: Path) -> Drawing:
    """SnapGene's GenBank export of pUC19 with joined, white, grey and origin-crossing features."""
    return draw_map(data_dir / "colour-test-snapgene.gbk")


def _page(drawing: Drawing, path: Path) -> tuple[str, Node]:
    html = drawing.write(path).read_text(encoding="utf-8")
    return html, parse(html)


@pytest.fixture(scope="module")
def colour_test_page(colour_test: Drawing, tmp_path_factory: pytest.TempPathFactory) -> Node:
    """The colour test's page, written once."""
    return _page(colour_test, tmp_path_factory.mktemp("colour-test") / "map.html")[1]


@pytest.fixture(scope="module")
def puc19_page(puc19_file: Path, tmp_path_factory: pytest.TempPathFactory) -> tuple[str, Node]:
    """pUC19's page with nothing asked for, written once."""
    return _page(draw_map(puc19_file), tmp_path_factory.mktemp("puc19") / "map.html")


def _shapes(page: Node) -> dict[str, Node]:
    """The map in each shape the page carries, by shape, the one shown first first."""
    shapes = page.find_all("div", cls="shape")
    shapes.sort(key=lambda shape: "hidden" in shape.attrs)
    return {shape.attrs["data-shape"]: shape for shape in shapes}


def _map(page: Node) -> Node:
    """The map the page shows first."""
    return next(iter(_shapes(page).values()))


def _items(page: Node, kind: str = "feature") -> dict[str, list[Node]]:
    """Each item's groups on the page, the arrows' first and any boxed label's after, by name."""
    groups: dict[str, list[Node]] = {}
    for group in page.find_all("g", data_kind=kind):
        groups.setdefault(group.attrs["data-name"], []).append(group)
    return groups


def _labels(page: Node, kind: str) -> list[list[tuple[str, str]]]:
    """Each label of a kind, as its stretches of text, each with its font weight."""
    return [
        [(str(text.children[0]), text.attrs["font-weight"]) for text in group.find_all("text")]
        for group in page.find_all("g", cls="label", data_kind=kind)
    ]


def _fills(group: Node) -> list[str]:
    return [path.attrs["fill"] for path in group.find_all("path") if path.attrs["fill"] != "none"]


def test_a_record_or_any_file_the_pipelines_read_is_drawn(
    puc19: SequenceRecord, puc19_file: Path, tmp_path: Path
) -> None:
    fasta = tmp_path / "insert.fasta"
    fasta.write_text(">insert\nACGTACGTAC\n")
    for source in (puc19, puc19_file, str(fasta)):
        written = draw_map(source).write(tmp_path / "map.html")
        assert written == tmp_path / "map.html"
        assert parse(written.read_text(encoding="utf-8")).find_all("svg")


@pytest.mark.parametrize("name", ["map.svg", "map.jpg", "map"])
def test_a_suffix_it_does_not_write_is_refused(
    puc19: SequenceRecord, tmp_path: Path, name: str
) -> None:
    with pytest.raises(ValueError, match=r"the suffix must be \.html, \.png, \.pdf"):
        draw_map(puc19).write(tmp_path / name)
    assert not (tmp_path / name).exists()


@pytest.fixture(scope="module")
def small() -> Drawing:
    """The smallest map to convert: a short circular record with no features."""
    return draw_map(SequenceRecord("ACGT" * 10, topology="circular", name="small"))


def _png_size_and_dpi(path: Path) -> tuple[int, int, float]:
    data = path.read_bytes()
    assert data.startswith(b"\x89PNG\r\n\x1a\n")
    width, height = struct.unpack(">II", data[16:24])
    per_metre, _, unit = struct.unpack(">IIB", data[data.index(b"pHYs") + 4 :][:9])
    assert unit == 1
    return width, height, per_metre * 0.0254


def test_a_png_is_drawn_at_300_dpi_unless_told_otherwise(small: Drawing, tmp_path: Path) -> None:
    extent = small.layout.extent
    written = {300: small.write(tmp_path / "300.png"), 36: small.write(tmp_path / "36.png", dpi=36)}
    for dpi, path in written.items():
        width, height, recorded = _png_size_and_dpi(path)
        assert width == pytest.approx(extent.width * dpi / 72, abs=1)
        assert height == pytest.approx(extent.height * dpi / 72, abs=1)
        assert recorded == pytest.approx(dpi, abs=0.1)


def _lines(shapes: Iterable[svg.Shape]) -> list[svg.Text | svg.Letters]:
    """Every line of text among `shapes`, in groups or not."""
    found: list[svg.Text | svg.Letters] = []
    for shape in shapes:
        if isinstance(shape, svg.Group):
            found += _lines(shape.shapes)
        elif isinstance(shape, svg.Text | svg.Letters):
            found.append(shape)
    return found


@pytest.mark.parametrize("opened", [False, True], ids=["circle", "line"])
def test_a_pdf_is_one_page_the_size_laid_out_every_letter_handed_over_as_its_outline(
    colour_test: Drawing, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, opened: bool
) -> None:
    if opened:
        colour_test = draw_map(colour_test.record, linear=True)
    handed: list[str] = []
    # vl-convert's linux-64 conda build ships no type stub, so pyright there sees no functions.
    to_pdf = vl_convert.svg_to_pdf  # pyright: ignore[reportAttributeAccessIssue]

    def recorded(image: str) -> bytes:
        handed.append(image)
        return to_pdf(image)

    monkeypatch.setattr(vl_convert, "svg_to_pdf", recorded)
    # A dpi counts for a PNG alone.
    [page] = pypdf.PdfReader(colour_test.write(tmp_path / "map.pdf", dpi=1)).pages
    extent = colour_test.layout.extent
    assert float(page.mediabox.width) == pytest.approx(extent.width, abs=0.01)
    assert float(page.mediabox.height) == pytest.approx(extent.height, abs=0.01)
    lines = _lines(colour_test.layout.shapes)
    assert any(isinstance(line, svg.Letters) for line in lines), "no name on an arrow"
    letters = [letter for line in lines for letter in line.font.letters(line.text, line.size)]
    [image] = handed
    assert not parse(image).find_all("text")
    assert len(parse(image).find_all("use")) == len([one for one in letters if one.outline])
    # No font to look up, and so no text to find.
    resources = page["/Resources"].get_object()
    assert isinstance(resources, DictionaryObject)
    assert "/Font" not in resources
    assert page.extract_text() == ""


@pytest.fixture(scope="module")
def rows() -> Drawing:
    """A short record's map and sequence view at ten bases a row, one row too tall for a page."""
    stacked = tuple(Feature(f"f{n}", "misc_feature", (Segment(20, 30),)) for n in range(12))
    record = SequenceRecord("ACGT" * 10, name="rows", features=stacked)
    return draw_map(record, sequence_view=True, bases_per_row=10, cut_sites=False)


def test_a_png_is_one_image_the_sequence_view_under_the_map(rows: Drawing, tmp_path: Path) -> None:
    assert rows.sequence_view is not None
    top, under = rows.layout.extent, rows.sequence_view.extent
    width, height, _ = _png_size_and_dpi(rows.write(tmp_path / "rows.png", dpi=36))
    assert width == pytest.approx(max(top.width, under.width) / 2, abs=1)
    assert height == pytest.approx((top.height + under.height) / 2, abs=1)


def test_a_pdf_has_the_map_then_as_many_whole_rows_to_a_page_as_fit(
    rows: Drawing, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    view = rows.sequence_view
    assert view is not None
    handed: list[str] = []
    to_pdf = convert.pdf

    def recorded(pages: Iterable[str]) -> bytes:
        handed.extend(pages)
        return to_pdf(handed)

    monkeypatch.setattr(convert, "pdf", recorded)
    pages = pypdf.PdfReader(rows.write(tmp_path / "rows.pdf")).pages
    boxes = [
        [float(one) for one in parse(page).find_all("svg")[0].attrs["viewbox"].split()]
        for page in handed
    ]
    assert len(pages) == len(boxes)
    for page, (_, _, width, height) in zip(pages, boxes, strict=True):
        assert float(page.mediabox.width) == pytest.approx(width, abs=0.01)
        assert float(page.mediabox.height) == pytest.approx(height, abs=0.01)
    extent = rows.layout.extent
    assert boxes[0] == pytest.approx([extent.x, extent.y, extent.width, extent.height], abs=0.01)
    # Each row lies whole on one page after the map's, in order, on pages the view's width across.
    tall = view.extent.width * math.sqrt(2)
    on = [
        [index for index, box in enumerate(boxes[1:]) if _inside(row.extent, box)]
        for row in view.rows
    ]
    assert all(len(found) == 1 for found in on)
    first = [found[0] for found in on]
    assert first == sorted(first)
    assert set(first) == set(range(len(boxes) - 1))
    assert all(box[2] == pytest.approx(view.extent.width, abs=0.01) for box in boxes[1:])
    assert all(box[3] >= tall - 0.01 for box in boxes[1:])
    # A page ends only where the next row would not fit it, and a row that fits no page has one of
    # its own, as tall as it needs.
    margin = view.rows[0].extent.y - view.extent.y
    for index, row in enumerate(view.rows[1:], start=1):
        if first[index] != first[index - 1]:
            top = boxes[1 + first[index - 1]][1]
            assert row.extent.y + row.extent.height + margin > top + tall
    assert all(
        boxes[1 + page][3] == pytest.approx(tall, abs=0.01)
        for page in first
        if first.count(page) > 1
    )
    assert max(first.count(page) for page in first) > 1
    assert any(box[3] > tall + 1 for box in boxes[1:])


def _inside(inner: Box, outer: list[float]) -> bool:
    x, y, width, height = outer
    return (
        x - 0.01 <= inner.x
        and inner.x + inner.width <= x + width + 0.01
        and y - 0.01 <= inner.y
        and inner.y + inner.height <= y + height + 0.01
    )


@pytest.mark.parametrize(
    ("drawn", "dpi", "advice"),
    [
        ("small", 10**7, "and a PDF at any size"),
        ("small", 0, "and a PDF at any size"),
        ("small", None, "and a PDF at any size"),
        ("rows", 10**7, "name a region to draw fewer bases, or write a PDF"),
    ],
    ids=["a side too long", "no pixels", "too many pixels", "with its sequence view"],
)
def test_a_png_too_large_or_small_to_draw_is_refused_before_it_is_drawn(
    request: pytest.FixtureRequest,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    drawn: str,
    dpi: float | None,
    advice: str,
) -> None:
    drawing: Drawing = request.getfixturevalue(drawn)
    if dpi is None:
        # More pixels than memory allows, though no side is too long.
        extent = drawing.layout.extent
        dpi = 1.01 * 72 * math.sqrt(convert.PNG_PIXELS / (extent.width * extent.height))
        assert max(extent.width, extent.height) * dpi / 72 < convert.PNG_SIDE
    monkeypatch.delattr(convert, "png")
    with pytest.raises(ValueError, match=rf"can be drawn at .* dpi.*{advice}"):
        drawing.write(tmp_path / "map.png", dpi=dpi)
    assert not (tmp_path / "map.png").exists()


def test_a_record_with_no_bases_is_refused() -> None:
    with pytest.raises(ValueError, match="no bases"):
        draw_map(SequenceRecord("", name="empty"))


def test_a_file_no_reader_reads_is_refused(tmp_path: Path) -> None:
    notes = tmp_path / "notes.txt"
    notes.write_text("ACGT")
    with pytest.raises(ValueError, match="no reader"):
        draw_map(notes)


def test_the_page_loads_nothing_over_the_network(puc19_page: tuple[str, Node]) -> None:
    html, page = puc19_page
    assert not page.find_all("link")
    assert not [node for node in page.iter() if {"src", "href", "xlink:href"} & set(node.attrs)]
    assert "@import" not in html
    assert set(re.findall(r"url\(([a-z]+):", html)) == {"data"}
    # The one address is the SVG namespace, which names the markup and is never fetched.
    assert set(re.findall(r"https?://[^\s\"'<>]+", html)) == {"http://www.w3.org/2000/svg"}


def test_the_page_embeds_the_font_subsets_it_measures_with(puc19_page: tuple[str, Node]) -> None:
    html, _ = puc19_page
    embedded = re.findall(r"url\(data:font/woff2;base64,([A-Za-z0-9+/=]+)\)", html)
    assert [base64.b64decode(data) for data in embedded] == [
        font.woff2() for font in (SANS, BOLD, MONO)
    ]


def test_the_page_stays_white_in_a_dark_colour_scheme(puc19_page: tuple[str, Node]) -> None:
    _, page = puc19_page
    [scheme] = [meta for meta in page.find_all("meta") if meta.attrs.get("name") == "color-scheme"]
    assert scheme.attrs["content"] == "light"


def test_the_centre_names_the_record_in_bold_over_its_length_inside_a_bp_scale(
    puc19_page: tuple[str, Node],
) -> None:
    shown = _map(puc19_page[1])
    weights = {text.text: text.attrs["font-weight"] for text in shown.find_all("text")}
    assert weights["pUC19"] == "700"
    assert weights["2686 bp"] == "400"
    [scale] = shown.find_all("g", cls="scale")
    assert [text.text for text in scale.find_all("text")] == ["500", "1000", "1500", "2000", "2500"]


def test_hovering_over_a_feature_shows_its_name_type_span_and_length_one_based(
    colour_test_page: Node,
) -> None:
    page = colour_test_page
    hover = {
        name: {key: groups[0].attrs[f"data-{key}"] for key in ("type", "span", "length")}
        for name, groups in _items(page).items()
    }
    assert hover["AmpR"] == {"type": "CDS", "span": "1626 .. 2486", "length": "861 bp"}
    assert hover["across"] == {"type": "misc_feature", "span": "2601 .. 64", "length": "150 bp"}
    assert hover["split"]["span"] == "101 .. 160, 201 .. 260, 301 .. 360"


def test_every_feature_is_drawn_with_its_name_on_its_arrow_or_boxed_and_source_switched_off(
    colour_test: Drawing, colour_test_page: Node
) -> None:
    page = colour_test_page
    shown = _map(page)
    [source] = [f.name for f in colour_test.record.features if f.type == "source"]
    items = _items(shown)
    assert sorted(items) == sorted(f.name for f in colour_test.record.features)
    [arrow] = items.pop(source)
    assert arrow.attrs["class"] == "feature off"
    assert not [g for g in shown.find_all("g", cls="off") if g is not arrow]
    on_arrows = set()
    for name, (arrows, *boxed) in items.items():
        assert _fills(arrows)
        if texts := arrows.find_all("text"):
            # One line of text, to search and copy, with each letter placed as it was laid out.
            [text] = texts
            assert (text.text, boxed) == (name, [])
            assert {len(text.attrs[key].split()) for key in ("x", "y", "rotate")} == {len(name)}
            on_arrows.add(name)
            continue
        [label] = boxed
        assert label.attrs["class"] == "feature label"
        assert [text.text for text in label.find_all("text")] == [name]
        assert len(label.find_all("rect")) == 1
        assert len(label.find_all("line")) == 1
    # The names SnapGene Viewer draws on arrows in its map of the same record.
    assert on_arrows == {"AmpR", "across", "lacZ-alpha", "ori"}


def test_features_draw_in_their_files_colours_segment_by_segment(
    colour_test_page: Node,
) -> None:
    page = colour_test_page
    items = _items(page)
    assert _fills(items["split"][0]) == ["#ff0000", "#00ff00", "#0000ff"]
    assert _fills(items["split-same-first"][0]) == ["#ffcc00", "#123abc"]
    assert _fills(items["white"][0]) == ["#ffffff"]
    # SnapGene's grey is a colour the file gives, not a missing one.
    assert _fills(items["uncoloured"][0]) == ["#a6acb3"]


def test_every_arrow_is_outlined_and_every_name_contrasts_with_its_fill(
    colour_test_page: Node,
) -> None:
    page = colour_test_page
    items = _items(page)
    for arrows, *_ in items.values():
        assert all(
            path.attrs["stroke"] not in ("none", path.attrs["fill"])
            for path in arrows.find_all("path")
            if path.attrs["fill"] != "none"
        )
    # The last group holds the name, on its arrow or boxed.
    fills = {
        name: [text.attrs["fill"] for text in groups[-1].find_all("text")]
        for name, groups in items.items()
    }
    assert fills["white"] == fills["ori"] == ["#000000"]
    assert fills["M13 rev"] == fills["lacZ-alpha"] == ["#ffffff"]


def test_a_name_on_its_arrow_contrasts_with_the_segment_it_sits_on(tmp_path: Path) -> None:
    white, dark = "#ffffff", "#1f3a93"
    record = SequenceRecord(
        "A" * 1000,
        topology="circular",
        features=(
            Feature("white promoter", "promoter", (Segment(0, 300),), color=white),
            Feature(
                "dark where it fits",
                "CDS",
                (Segment(400, 420), Segment(500, 950, color=dark)),
                strand=Strand.FORWARD,
                color=white,
            ),
        ),
    )
    items = _items(_page(draw_map(record), tmp_path / "map.html")[1])
    fills = {
        name: [text.attrs["fill"] for text in groups[0].find_all("text")]
        for name, groups in items.items()
    }
    assert fills == {"white promoter": ["#000000"], "dark where it fits": ["#ffffff"]}


def test_a_feature_with_no_colour_takes_its_groups_colour(tmp_path: Path) -> None:
    types = {
        "CDS": "#77aadd",
        "promoter": "#ee8866",
        "terminator": "#ffaabb",
        "rep_origin": "#eedd88",
        "protein_bind": "#99ddff",
        "ncRNA": "#44bb99",
        "LTR": "#bbcc33",
        "misc_feature": "#dddddd",
    }
    features = [
        Feature(kind, kind, (Segment(100 * i, 100 * i + 50),), color=None if i % 2 else "noColor")
        for i, kind in enumerate(types)
    ]
    features.append(
        Feature(
            "class terminator",
            "regulatory",
            (Segment(900, 950),),
            strand=Strand.FORWARD,
            qualifiers={"regulatory_class": ("terminator",)},
        )
    )
    record = SequenceRecord("A" * 1000, topology="circular", features=tuple(features))
    items = _items(_page(draw_map(record), tmp_path / "map.html")[1])
    assert {kind: _fills(items[kind][0])[0] for kind in types} == types
    assert _fills(items["class terminator"][0]) == ["#ffaabb"]


def test_a_name_is_written_as_the_face_draws_it_and_hovers_as_written(tmp_path: Path) -> None:
    name = '<b>"lac" & ﻿Z</b>'
    # One long enough to carry its name on its arrow, and one whose name is boxed.
    features = (Feature(name, "CDS", (Segment(0, 30),)), Feature(name, "CDS", (Segment(50, 51),)))
    record = SequenceRecord("A" * 100, topology="circular", features=features)
    _, page = _page(draw_map(record), tmp_path / "map.html")
    groups = _items(_map(page))[name]
    assert {group.attrs["data-name"] for group in groups} == {name}
    assert not page.find_all("b")
    texts = [text.text for group in groups for text in group.find_all("text")]
    assert texts == ['<b>"lac" & Z</b>'] * 2


def test_the_shipped_unique_cutters_are_labelled_where_snapgene_numbers_their_cuts(
    puc19_page: tuple[str, Node],
) -> None:
    _, page = puc19_page
    labels = _labels(_map(page), "cut_site")
    # SnapGene Viewer's own map of this file numbers each of these cuts the same, BsaI's on the
    # bottom strand and SapI's outside its site included.
    assert ["".join(text for text, _ in label) for label in labels] == [
        "NdeI (184)",
        "EcoRI (396)",
        "SacI (406)",
        "KpnI - XmaI (412)",
        "SmaI (414)",
        "BamHI (417)",
        "XbaI (423)",
        "SalI (429)",
        "PstI - SbfI (439)",
        "SphI (445)",
        "HindIII (447)",
        "BspQI - SapI (690)",
        "BsaI (1760)",
    ]
    assert labels[3] == [("KpnI", "700"), (" - ", "400"), ("XmaI", "700"), (" (412)", "400")]


def test_named_enzymes_draw_every_cut_site_bold_only_for_one_that_cuts_once(
    puc19: SequenceRecord, tmp_path: Path
) -> None:
    drawing = draw_map(puc19, enzymes=["BsmBI", "EcoRI-HF", "EcoRI"])
    _, page = _page(drawing, tmp_path / "map.html")
    assert _labels(_map(page), "cut_site") == [
        [("BsmBI", "400"), (" (3)", "400")],
        [("BsmBI", "400"), (" (45)", "400")],
        [("EcoRI", "700"), (" (396)", "400")],
    ]


def test_an_enzyme_no_shipped_enzyme_answers_to_is_refused(puc19: SequenceRecord) -> None:
    with pytest.raises(KeyError, match="EcoRJ"):
        draw_map(puc19, enzymes=["EcoRI", "EcoRJ"])


def test_each_primer_is_drawn_in_purple_at_every_binding_site_labelled_with_its_span(
    tmp_path: Path,
) -> None:
    sites = (BindingSite(990, 1010, Strand.FORWARD), BindingSite(100, 120, Strand.REVERSE))
    record = SequenceRecord(
        "A" * 1000,
        topology="circular",
        features=(
            Feature("uncoloured", "primer_bind", (Segment(300, 320),)),
            Feature("coloured", "primer_bind", (Segment(400, 420),), color="#a020f0"),
        ),
        primers=(Primer("M13 fwd", "GTAAAACGACGGCCAGT", binding_sites=sites),),
    )
    page = _map(_page(draw_map(record), tmp_path / "map.html")[1])
    [(across, reverse, *labels)] = _items(page, "primer").values()
    assert [group.attrs["data-span"] for group in (across, reverse)] == ["991 .. 10", "101 .. 120"]
    assert across.attrs["data-type"] == "primer"
    assert across.attrs["data-length"] == "20 bp"
    assert _fills(across) == _fills(reverse) == ["#aa3377"]
    assert [[text.attrs["fill"] for text in label.find_all("text")] for label in labels] == [
        ["#aa3377"],
        ["#aa3377"],
    ]
    assert sorted(label.find_all("text")[0].text for label in labels) == [
        "M13 fwd (101 .. 120)",
        "M13 fwd (991 .. 10)",
    ]
    features = _items(page)
    assert _fills(features["uncoloured"][0]) == ["#aa3377"]
    assert _fills(features["coloured"][0]) == ["#a020f0"]


@pytest.mark.parametrize(
    ("switches", "off"),
    [
        ({}, {"feature": {"source"}}),
        ({"features": False}, {"feature": {"CDS", "rep_origin", "source"}}),
        ({"primers": False}, {"feature": {"source"}, "primer": {"primer"}}),
        ({"cut_sites": False}, {"feature": {"source"}, "cut_site": {"cut site"}}),
        ({"hide_types": ["CDS"], "source": True}, {"feature": {"CDS"}}),
        ({"hide_types": ["CDS", "source"], "source": True}, {"feature": {"CDS", "source"}}),
    ],
)
def test_the_page_carries_every_layer_and_type_and_switches_off_only_what_was_asked(
    tmp_path: Path, switches: dict, off: dict[str, set[str]]
) -> None:
    drawing = draw_map(_layered(), **switches)
    _, page = _page(drawing, tmp_path / "map.html")
    everything = {
        "feature": {"CDS", "rep_origin", "source"},
        "primer": {"primer"},
        "cut_site": {"cut site"},
    }
    for shown in _shapes(page).values():
        found: dict[str, set[str]] = {}
        switched: dict[str, set[str]] = {}
        for group in shown.find_all("g"):
            if "data-kind" in group.attrs:
                kind, kind_type = group.attrs["data-kind"], group.attrs["data-type"]
                found.setdefault(kind, set()).add(kind_type)
                if "off" in group.attrs["class"].split():
                    switched.setdefault(kind, set()).add(kind_type)
        assert (found, switched) == (everything, off)
    kinds = {"feature": "features", "primer": "primers", "cut_site": "cut_sites"}
    types_off = {*switches.get("hide_types", ()), *(() if switches.get("source") else ["source"])}
    assert {
        (box.attrs["name"], box.attrs["value"]): "checked" in box.attrs
        for box in page.find_all("input")
        if box.attrs["type"] == "checkbox"
    } == {
        **{("kind", kind): switches.get(word, True) for kind, word in kinds.items()},
        **{("type", one): one not in types_off for one in everything["feature"]},
    }
    # A PNG or a PDF draws only what shows.
    layout = drawing.layout
    drawn = {(one.item.kind, one.item.type) for one in (*layout.arrows, *layout.labels)}
    assert drawn == {
        (kind, kind_type)
        for kind, types in everything.items()
        for kind_type in types
        if kind_type not in off.get(kind, ())
    }


@pytest.mark.parametrize(
    ("switches", "counts"),
    [
        ({}, {"circular": 1, "linear": 1, "sequence_view": 1}),
        ({"primers": False}, {"circular": 2, "linear": 1, "sequence_view": 2}),
        ({"linear": True, "primers": False}, {"circular": 1, "linear": 2, "sequence_view": 2}),
        ({"region": "mcs"}, {"linear": 1, "sequence_view": 1}),
        ({"region": "mcs", "cut_sites": False}, {"linear": 2, "sequence_view": 2}),
    ],
)
def test_each_distinct_layout_runs_once_when_first_needed_and_is_kept(
    puc19: SequenceRecord,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    switches: dict,
    counts: dict[str, int],
) -> None:
    primer = Primer("M13 fwd", "GTAAAACG", binding_sites=(BindingSite(378, 395, Strand.FORWARD),))
    record = dataclasses.replace(puc19, primers=(primer,))
    counted: dict[str, int] = {}
    for module in (circular, linear, sequence_view):

        def lay_out(
            *args: object, module: Any = module, real: Any = module.layout, **kwargs: object
        ):
            key = module.__name__.rsplit(".", 1)[-1]
            counted[key] = counted.get(key, 0) + 1
            return real(*args, **kwargs)

        monkeypatch.setattr(module, "layout", lay_out)
    drawing = draw_map(record, sequence_view=True, bases_per_row=200, **switches)
    assert counted == {}
    drawing.write(tmp_path / "map.html")
    for _ in range(2):
        # What a PNG and a PDF draw, and the notice the command prints.
        assert drawing.layout.extent
        assert drawing.sequence_view is not None
        assert drawing.hidden == ()
    assert counted == counts


def _texts(page: Node, cls: str = "") -> list[str]:
    """Each line of text on the page, or in the groups of a class, in the order drawn."""
    groups = page.find_all("g", cls=cls) if cls else [page]
    return [text.text for group in groups for text in group.find_all("text")]


def _dots(page: Node) -> int:
    """How many dots mark the ends of a line the molecule carries on past."""
    return sum(len(group.find_all("circle")) for group in page.find_all("g", cls="ends"))


@pytest.mark.parametrize(
    ("record", "opened", "line", "dots"),
    [
        ("puc19", False, False, 0),
        ("puc19", True, True, 6),
        ("gfp", False, True, 0),
        ("gfp", True, True, 0),
    ],
)
def test_a_linear_record_is_always_a_line_and_a_circular_one_opens_with_dots_at_each_end(
    request: pytest.FixtureRequest, tmp_path: Path, record: str, opened: bool, line: bool, dots: int
) -> None:
    drawing = draw_map(request.getfixturevalue(record), linear=opened, cut_sites=False)
    kind = linear.LinearMap if line else circular.CircularMap
    assert isinstance(drawing.layout, kind)
    assert _dots(_map(_page(drawing, tmp_path / "map.html")[1])) == dots


@pytest.mark.parametrize(
    ("record", "switches", "shapes"),
    [
        ("puc19", {}, ["circle", "line"]),
        ("puc19", {"linear": True}, ["line", "circle"]),
        ("puc19", {"region": "mcs"}, ["line"]),
        ("gfp", {}, ["line"]),
    ],
)
def test_a_circular_record_drawn_whole_flips_between_circle_and_line_at_top_right(
    request: pytest.FixtureRequest, tmp_path: Path, record: str, switches: dict, shapes: list[str]
) -> None:
    _, page = _page(draw_map(request.getfixturevalue(record), **switches), tmp_path / "map.html")
    carried = _shapes(page)
    assert list(carried) == shapes
    assert [len(shape.find_all("svg")) for shape in carried.values()] == [1] * len(shapes)
    radios = {
        box.attrs["value"]: "checked" in box.attrs for box in page.find_all("input", name="shape")
    }
    if len(shapes) == 1:
        assert not radios
    else:
        assert radios == {"circle": shapes[0] == "circle", "line": shapes[0] == "line"}
        [bar] = page.find_all("form", cls="switches")
        assert [one.attrs["class"] for one in bar.find_all("fieldset")][-1] == "shapes"


def test_a_region_named_by_a_feature_is_a_line_keeping_the_records_numbering(
    puc19: SequenceRecord, tmp_path: Path
) -> None:
    drawing = draw_map(puc19, region="mcs")
    assert isinstance(drawing.layout, linear.LinearMap)
    _, page = _page(drawing, tmp_path / "map.html")
    assert "396 .. 452 (57 bp)" in _texts(page)
    assert _texts(page, "scale") == ["400", "410", "420", "430", "440", "450"]
    assert ["".join(text for text, _ in label) for label in _labels(page, "cut_site")] == [
        "EcoRI (396)",
        "SacI (406)",
        "KpnI - XmaI (412)",
        "SmaI (414)",
        "BamHI (417)",
        "XbaI (423)",
        "SalI (429)",
        "PstI - SbfI (439)",
        "SphI (445)",
        "HindIII (447)",
    ]
    features = {name: groups[0].attrs["data-span"] for name, groups in _items(page).items()}
    assert features == {"lacZ\N{GREEK SMALL LETTER ALPHA}": "146 .. 469", "MCS": "396 .. 452"}
    # A region at the start of a linear record carries on past its end only.
    line = draw_map(dataclasses.replace(puc19, topology="linear"), region=(0, 100))
    assert _dots(_page(line, tmp_path / "line.html")[1]) == 3


@pytest.mark.parametrize(
    ("region", "title", "scale"),
    [
        ((2679, 2696), "2680 .. 10 (17 bp)", ["2680", "2682", "2684", "2", "4", "6", "8", "10"]),
        ("across", "2601 .. 64 (150 bp)", ["2620", "2640", "2660", "2680", "20", "40", "60"]),
    ],
)
def test_a_region_runs_across_the_origin_of_a_circular_record(
    colour_test: Drawing,
    tmp_path: Path,
    region: tuple[int, int] | str,
    title: str,
    scale: list[str],
) -> None:
    _, page = _page(draw_map(colour_test.record, region=region), tmp_path / "map.html")
    assert title in _texts(page)
    assert _texts(page, "scale") == scale
    assert _dots(page) == 6


def test_a_region_a_feature_name_several_features_share_is_the_first_one(tmp_path: Path) -> None:
    features = tuple(
        Feature(name, "protein_bind", (Segment(start, start + 34),))
        for name, start in (("loxP", 100), ("LoxP", 500))
    )
    record = SequenceRecord("A" * 1000, topology="circular", features=features)
    _, page = _page(draw_map(record, region="LOXP"), tmp_path / "map.html")
    assert "101 .. 134 (34 bp)" in _texts(page)


def test_a_cutter_unique_in_a_region_but_not_in_the_record_is_not_bold_or_shown_by_default(
    puc19: SequenceRecord, tmp_path: Path
) -> None:
    # BsmBI cuts pUC19 after bases 3 and 45.
    _, named = _page(draw_map(puc19, region=(0, 20), enzymes=["BsmBI"]), tmp_path / "named.html")
    assert _labels(named, "cut_site") == [[("BsmBI", "400"), (" (3)", "400")]]
    _, shipped = _page(draw_map(puc19, region=(0, 20)), tmp_path / "shipped.html")
    assert not _labels(shipped, "cut_site")


@pytest.mark.parametrize(
    ("record", "region", "message"),
    [
        ("gfp", "nothing", "record 'GFP' has no feature called 'nothing'"),
        ("gfp", (0, 718), "does not lie on the linear record 'GFP' of 717 bases"),
        ("gfp", (700, 10), "does not lie on the linear record"),
        ("gfp", (-1, 10), "does not lie on the linear record"),
        ("puc19", (2680, 10), "does not lie on the circular record 'pUC19' of 2686 bases"),
        ("puc19", (2686, 2690), "does not lie on the circular record"),
        ("puc19", (100, 2787), "does not lie on the circular record"),
    ],
)
def test_a_region_naming_no_feature_or_lying_off_the_record_is_refused(
    request: pytest.FixtureRequest, record: str, region: str | tuple[int, int], message: str
) -> None:
    with pytest.raises(ValueError, match=re.escape(message)):
        draw_map(request.getfixturevalue(record), region=region)


def _rows(page: Node) -> list[Node]:
    """The sequence view's rows, top to bottom."""
    [view] = page.find_all("figure", cls="sequence-view")
    return view.find_all("g", cls="row")


def test_the_sequence_view_goes_beside_the_map_only_when_asked_for(
    gfp: SequenceRecord, tmp_path: Path
) -> None:
    alone = draw_map(gfp)
    assert alone.sequence_view is None
    assert not _page(alone, tmp_path / "alone.html")[1].find_all("figure", cls="sequence-view")
    _, page = _page(draw_map(gfp, sequence_view=True), tmp_path / "view.html")
    [main] = page.find_all("main")
    assert [figure.attrs["class"] for figure in main.find_all("figure")] == ["map", "sequence-view"]


def test_each_row_shows_a_ruler_both_strands_and_its_last_position_as_many_bases_as_asked(
    gfp: SequenceRecord, tmp_path: Path
) -> None:
    rows = _rows(
        _page(draw_map(gfp, sequence_view=True, bases_per_row=100), tmp_path / "a.html")[1]
    )
    assert len(rows) == 8
    top = gfp.sequence[:100]
    assert _texts(rows[0], "top") == [top]
    assert _texts(rows[0], "bottom") == [top.translate(str.maketrans("ACGT", "TGCA"))]
    assert _texts(rows[0], "ruler") == ["10 20 30 40 50 60 70 80 90 100"]
    assert _texts(rows[0], "position") == ["100"]
    assert _texts(rows[-1], "top") == [gfp.sequence[700:]]
    assert _texts(rows[-1], "ruler") == ["710"]
    assert _texts(rows[-1], "position") == ["717"]
    [bar] = rows[0].find_all("g", data_kind="feature")
    assert (bar.attrs["data-name"], bar.attrs["data-span"]) == ("GFP", "1 .. 717")
    one = _rows(
        _page(draw_map(gfp, sequence_view=True, both_strands=False), tmp_path / "b.html")[1]
    )
    assert len(one) == 12
    assert _texts(one[0], "top") == [gfp.sequence[:60]]
    assert not [row for row in one if row.find_all("g", cls="bottom")]


def _amino_acids(drawing: Drawing, name: str) -> list[str]:
    """A feature's amino acids in the order the sequence view reads them along its strand."""
    assert drawing.sequence_view is not None
    found = []
    for index, row in enumerate(drawing.sequence_view.rows):
        lines = [
            line
            for translation in row.translations
            if translation.item.name == name
            for line in (translation.letters, translation.stops)
            if line
        ]
        for line in lines:
            at = 0
            for word in line.text.split(" "):
                found.append((index, line.places[at].x, word, line.fill))
                at += len(word) + 1
    [feature] = [one for one in drawing.record.features if one.name == name]
    found.sort(reverse=feature.strand == Strand.REVERSE)
    assert {fill for *_, word, fill in found if word == "Ter"} <= {sequence_view.STOP}
    assert {fill for *_, word, fill in found if word != "Ter"} <= {"#252525"}
    return [word for *_, word, _ in found]


def _three_letters(protein: str) -> list[str]:
    from Bio.SeqUtils import seq3

    three = seq3(protein)
    return [three[at : at + 3] for at in range(0, len(three), 3)]


def test_every_cds_shows_its_three_letter_translation_as_snapgene_translates_it(
    puc19: SequenceRecord, colour_test: Drawing, tmp_path: Path
) -> None:
    drawing = draw_map(puc19, sequence_view=True)
    for feature in (one for one in puc19.features if one.type == "CDS"):
        # SnapGene marks where AmpR's segments join with a comma.
        protein = str(feature.qualifiers["translation"][0]).replace(",", "")
        assert _amino_acids(drawing, feature.name) == _three_letters(protein)
    # Read across three joined segments, stops and all.
    joined = draw_map(colour_test.record, sequence_view=True)
    [split] = [one for one in colour_test.record.features if one.name == "split"]
    assert _amino_acids(joined, "split") == _three_letters(str(split.qualifiers["translation"][0]))
    [view] = _page(drawing, tmp_path / "map.html")[1].find_all("figure", cls="sequence-view")
    fills = {
        text.attrs["fill"]
        for group in view.find_all("g", cls="translation")
        for text in group.find_all("text")
    }
    assert fills == {"#252525", sequence_view.STOP}


def test_a_regions_sequence_view_keeps_the_records_numbering_across_the_origin(
    puc19: SequenceRecord, tmp_path: Path
) -> None:
    drawing = draw_map(puc19, region=(2679, 2696), sequence_view=True, bases_per_row=10)
    rows = _rows(_page(drawing, tmp_path / "map.html")[1])
    assert [_texts(row, "top") for row in rows] == [
        [puc19.sequence[2679:] + puc19.sequence[:3]],
        [puc19.sequence[3:10]],
    ]
    assert [_texts(row, "ruler") for row in rows] == [["2680"], ["10"]]
    assert [_texts(row, "position") for row in rows] == [["3"], ["10"]]


def _layered() -> SequenceRecord:
    """A small circular record with a feature of two types, its source, a primer, and a site each
    for two shipped enzymes that cut it once."""
    bases = list("ACGT" * 75)
    bases[100:106] = "GAATTC"
    bases[200:206] = "AAGCTT"
    return SequenceRecord(
        "".join(bases),
        topology="circular",
        name="layered",
        features=(
            Feature("layered", "source", (Segment(0, 300),)),
            Feature("coding", "CDS", (Segment(10, 70),), strand=Strand.FORWARD),
            Feature("origin", "rep_origin", (Segment(150, 190),)),
        ),
        primers=(
            Primer(
                "forward",
                "ACGTACGTACGTACGTAC",
                binding_sites=(BindingSite(120, 138, Strand.FORWARD),),
            ),
        ),
    )


@pytest.mark.parametrize(
    ("switches", "off"),
    [
        ({}, {"feature": {"layered"}}),
        ({"features": False}, {"feature": {"coding", "origin", "layered"}}),
        (
            {"primers": False, "enzymes": ["HindIII"]},
            {"feature": {"layered"}, "primer": {"forward"}},
        ),
        ({"cut_sites": False}, {"feature": {"layered"}, "cut_site": {"EcoRI", "HindIII"}}),
        ({"hide_types": ["CDS"], "source": True}, {"feature": {"coding"}}),
    ],
)
def test_the_sequence_view_takes_the_maps_enzymes_layers_and_feature_types(
    tmp_path: Path, switches: dict, off: dict[str, set[str]]
) -> None:
    drawing = draw_map(_layered(), sequence_view=True, **switches)
    everything = {
        "feature": {"coding", "origin", "layered"},
        "primer": {"forward"},
        "cut_site": set(switches.get("enzymes", ["EcoRI", "HindIII"])),
    }
    page = _page(drawing, tmp_path / "map.html")[1]
    [view] = page.find_all("figure", cls="sequence-view")
    found: dict[str, set[str]] = {}
    switched: dict[str, set[str]] = {}
    for group in view.find_all("g"):
        if "data-kind" in group.attrs:
            kind, name = group.attrs["data-kind"], group.attrs["data-name"]
            found.setdefault(kind, set()).add(name)
            if "off" in group.attrs["class"].split():
                switched.setdefault(kind, set()).add(name)
    assert (found, switched) == (everything, off)
    # Each item is known by the same details in both views, for a switch or a click to reach.
    details = ["data-kind", "data-name", "data-type", "data-span", "data-length"]

    def known(figure: Node) -> set[tuple[str | None, ...]]:
        return {
            tuple(group.attrs.get(one) for one in details)
            for group in figure.find_all("g")
            if "data-kind" in group.attrs
        }

    assert known(view) == known(_map(page))
    # A PNG or a PDF draws only what shows.
    assert drawing.sequence_view is not None
    drawn = {
        (one.item.kind, one.item.name)
        for row in drawing.sequence_view.rows
        for one in (*row.bars, *row.arrows, *row.cuts, *row.labels)
    }
    assert drawn == {
        (kind, name)
        for kind, names in everything.items()
        for name in names
        if name not in off.get(kind, ())
    }


def test_the_sequence_view_stacks_the_names_at_one_cut_bold_where_one_cuts_once(
    puc19: SequenceRecord, tmp_path: Path
) -> None:
    drawing = draw_map(
        puc19, region=(660, 720), sequence_view=True, enzymes=["SapI", "BspQI", "BsaI"]
    )
    [view] = _page(drawing, tmp_path / "map.html")[1].find_all("figure", cls="sequence-view")
    assert _labels(view, "cut_site") == [[("BspQI", "700"), ("SapI", "700")]]


def test_the_sequence_view_draws_each_primers_tail_and_mismatches_where_snapgene_binds_it(
    data_dir: Path,
) -> None:
    drawing = draw_map(data_dir / "primer-test.dna", sequence_view=True, cut_sites=False)
    assert drawing.sequence_view is not None
    cell = sequence_view.CELL
    tails: dict[str, tuple[float, float]] = {}
    marks = set()
    for row in drawing.sequence_view.rows:
        for tail in row.tails:
            ends = [row.start + point.x / cell for point in tail.points]
            tails[tail.item.name] = (min(ends), max(ends))
        marks |= {(mark.item.name, row.start + int(mark.box.x // cell)) for mark in row.mismatches}
    # Each tail is 9 bases, 5' of the binding site on the primer's own strand.
    assert tails == {"fwd-BsaI-tail": (21, 30), "rev-mismatch": (272, 281)}
    assert marks == {("rev-mismatch", 263)}


def test_a_sequence_view_past_100_kb_is_refused_asking_for_a_region_and_a_map_is_not() -> None:
    record = SequenceRecord("ACGT" * 25_001, name="long")
    with pytest.raises(ValueError, match=r"at most 100,000 bases.*100,004: name a region"):
        draw_map(record, sequence_view=True)
    assert isinstance(draw_map(record, cut_sites=False).layout, linear.LinearMap)
    region = draw_map(record, region=(99_000, 99_060), sequence_view=True, cut_sites=False)
    assert region.sequence_view is not None
    assert [(row.start, row.end) for row in region.sequence_view.rows] == [(99_000, 99_060)]


@pytest.mark.parametrize("bases_per_row", [0, -60])
def test_a_row_of_no_bases_is_refused(gfp: SequenceRecord, bases_per_row: int) -> None:
    with pytest.raises(ValueError, match="holds at least 1 base"):
        draw_map(gfp, sequence_view=True, bases_per_row=bases_per_row)


@pytest.mark.parametrize("opened", [False, True], ids=["circle", "line"])
def test_a_crowded_map_hides_cut_sites_before_the_primer_among_them_and_says_what_it_hid(
    tmp_path: Path, opened: bool
) -> None:
    drawing = draw_map(crowds.ecori_crowd(), enzymes=["EcoRI", "HindIII"], linear=opened)
    hidden = drawing.hidden
    assert hidden
    assert {(item.kind, item.name) for item in hidden} == {("cut_site", "EcoRI")}
    # Labels of one length hide in the record's order.
    cuts = [item.spans[0].start for item in hidden]
    assert cuts == sorted(cuts)
    assert set(cuts) < set(range(101, 500, 8))
    page = _map(_page(drawing, tmp_path / "map.html")[1])
    [notice] = page.find_all("g", cls="notice")
    assert notice.text == f"{len(hidden)} enzyme sites are hidden"
    assert json.loads(notice.attrs["data-hidden"]) == [
        {"label": f"EcoRI ({cut})", "kind": "cut_site", "type": "cut site"} for cut in cuts
    ]
    shown = {
        kind: ["".join(text for text, _ in label) for label in _labels(page, kind)]
        for kind in ("cut_site", "primer")
    }
    # The HindIII sites lie apart from the crowd, so they stay however often HindIII cuts.
    assert {"HindIII (1701)", "HindIII (2001)"} <= set(shown["cut_site"])
    assert len(shown["cut_site"]) + len(hidden) == 52
    assert shown["primer"] == ["among them (201 .. 220)"]


def _placed(page: Node) -> list[tuple[str, str, str]]:
    """Where each line of text the map draws starts, but for the notice's."""
    notices = [id(text) for notice in page.find_all("g", cls="notice") for text in notice.iter()]
    return [
        (text.attrs["x"], text.attrs["y"], text.text)
        for text in page.find_all("text")
        if id(text) not in notices
    ]


@pytest.fixture(scope="module")
def primed_crowd(tmp_path_factory: pytest.TempPathFactory) -> tuple[SequenceRecord, Node]:
    """The EcoRI crowd with twenty primers among its sites, and the map its page shows whole."""
    primers = tuple(
        Primer(f"p{i}", "ACGT", binding_sites=(BindingSite(at, at + 20, Strand.FORWARD),))
        for i, at in enumerate(range(150, 310, 8))
    )
    record = dataclasses.replace(crowds.ecori_crowd(), primers=primers)
    path = tmp_path_factory.mktemp("crowd") / "map.html"
    return record, _map(_page(draw_map(record, enzymes=["EcoRI", "HindIII"]), path)[1])


@pytest.mark.parametrize(
    ("switches", "on_the_page", "in_a_pdf"),
    [
        (
            {},
            "50 enzyme sites and 3 primers are hidden",
            "50 enzyme sites and 3 primers are hidden",
        ),
        ({"primers": False}, "50 enzyme sites are hidden", "24 enzyme sites are hidden"),
        ({"cut_sites": False}, "3 primers are hidden", "3 primers are hidden"),
        ({"cut_sites": False, "primers": False}, "", ""),
    ],
)
def test_the_page_says_what_its_own_map_hid_that_shows_and_the_drawing_what_a_pdf_hid(
    primed_crowd: tuple[SequenceRecord, Node],
    tmp_path: Path,
    switches: dict,
    on_the_page: str,
    in_a_pdf: str,
) -> None:
    record, whole = primed_crowd
    drawing = draw_map(record, enzymes=["EcoRI", "HindIII"], **switches)
    assert layers.notice(drawing.hidden) == in_a_pdf
    page = _map(_page(drawing, tmp_path / "map.html")[1])
    [notice] = page.find_all("g", cls="notice")
    assert ("off" in notice.attrs["class"].split()) == (not on_the_page)
    if on_the_page:
        assert notice.text == on_the_page
    # The page lays out every label whatever shows, so none moves when a switch flips.
    [listed] = whole.find_all("g", cls="notice")
    assert notice.attrs["data-hidden"] == listed.attrs["data-hidden"]
    assert _placed(page) == _placed(whole)


def test_a_map_with_room_for_every_label_hides_none_and_says_nothing(
    puc19_file: Path, puc19_page: tuple[str, Node]
) -> None:
    assert draw_map(puc19_file).hidden == ()
    _, page = puc19_page
    assert not page.find_all("g", cls="notice")
