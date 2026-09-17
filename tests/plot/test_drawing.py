"""A record drawn through `draw_map` and written as a page: what the page carries, and refusals."""

import base64
import re
from pathlib import Path

import pytest

from liulab_mbio.plot import Drawing, draw_map
from liulab_mbio.plot.fonts import BOLD, MONO, SANS
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand

from ..html import Node, parse


@pytest.fixture(scope="module")
def colour_test(data_dir: Path) -> Drawing:
    """SnapGene's GenBank export of pUC19 with joined, white, grey and origin-crossing features."""
    return draw_map(data_dir / "colour-test-snapgene.gbk")


def _page(drawing: Drawing, path: Path) -> tuple[str, Node]:
    html = drawing.write(path).read_text(encoding="utf-8")
    return html, parse(html)


def _items(page: Node) -> dict[str, list[Node]]:
    """Each item's groups on the page, the arrows' first and the label's after, by name."""
    groups: dict[str, list[Node]] = {}
    for group in page.find_all("g"):
        if "data-kind" in group.attrs:
            groups.setdefault(group.attrs["data-name"], []).append(group)
    return groups


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


@pytest.mark.parametrize("name", ["map.png", "map.pdf", "map.svg", "map"])
def test_a_suffix_it_does_not_write_is_refused(
    puc19: SequenceRecord, tmp_path: Path, name: str
) -> None:
    with pytest.raises(ValueError, match=r"the suffix must be \.html"):
        draw_map(puc19).write(tmp_path / name)
    assert not (tmp_path / name).exists()


def test_a_record_with_no_bases_is_refused() -> None:
    with pytest.raises(ValueError, match="no bases"):
        draw_map(SequenceRecord("", name="empty"))


def test_a_file_no_reader_reads_is_refused(tmp_path: Path) -> None:
    notes = tmp_path / "notes.txt"
    notes.write_text("ACGT")
    with pytest.raises(ValueError, match="no reader"):
        draw_map(notes)


def test_the_page_loads_nothing_over_the_network(puc19_file: Path, tmp_path: Path) -> None:
    html, page = _page(draw_map(puc19_file), tmp_path / "map.html")
    assert not page.find_all("link")
    assert not [node for node in page.iter() if {"src", "href", "xlink:href"} & set(node.attrs)]
    assert "@import" not in html
    assert set(re.findall(r"url\(([a-z]+):", html)) == {"data"}
    # The one address is the SVG namespace, which names the markup and is never fetched.
    assert set(re.findall(r"https?://[^\s\"'<>]+", html)) == {"http://www.w3.org/2000/svg"}


def test_the_page_embeds_the_font_subsets_it_measures_with(
    puc19_file: Path, tmp_path: Path
) -> None:
    html, _ = _page(draw_map(puc19_file), tmp_path / "map.html")
    embedded = re.findall(r"url\(data:font/woff2;base64,([A-Za-z0-9+/=]+)\)", html)
    assert [base64.b64decode(data) for data in embedded] == [
        font.woff2() for font in (SANS, BOLD, MONO)
    ]


def test_the_page_stays_white_in_a_dark_colour_scheme(puc19_file: Path, tmp_path: Path) -> None:
    _, page = _page(draw_map(puc19_file), tmp_path / "map.html")
    [scheme] = [meta for meta in page.find_all("meta") if meta.attrs.get("name") == "color-scheme"]
    assert scheme.attrs["content"] == "light"


def test_the_centre_names_the_record_in_bold_over_its_length_inside_a_bp_scale(
    puc19_file: Path, tmp_path: Path
) -> None:
    _, page = _page(draw_map(puc19_file), tmp_path / "map.html")
    weights = {text.text: text.attrs["font-weight"] for text in page.find_all("text")}
    assert weights["pUC19"] == "700"
    assert weights["2686 bp"] == "400"
    [scale] = page.find_all("g", cls="scale")
    assert [text.text for text in scale.find_all("text")] == ["500", "1000", "1500", "2000", "2500"]


def test_hovering_over_a_feature_shows_its_name_type_span_and_length_one_based(
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
    hover = {
        name: {key: groups[0].attrs[f"data-{key}"] for key in ("type", "span", "length")}
        for name, groups in _items(page).items()
    }
    assert hover["AmpR"] == {"type": "CDS", "span": "1626 .. 2486", "length": "861 bp"}
    assert hover["across"] == {"type": "misc_feature", "span": "2601 .. 64", "length": "150 bp"}
    assert hover["split"]["span"] == "101 .. 160, 201 .. 260, 301 .. 360"


def test_every_feature_but_source_is_drawn_with_its_name_boxed(
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
    items = _items(page)
    names = [f.name for f in colour_test.record.features if f.type != "source"]
    assert sorted(items) == sorted(names)
    assert not [g for g in page.find_all("g") if g.attrs.get("data-type") == "source"]
    for name, (arrows, label) in items.items():
        assert label.attrs["class"] == "feature label"
        assert [text.text for text in label.find_all("text")] == [name]
        assert len(label.find_all("rect")) == 1
        assert len(label.find_all("line")) == 1
        assert _fills(arrows)


def test_features_draw_in_their_files_colours_segment_by_segment(
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
    items = _items(page)
    assert _fills(items["split"][0]) == ["#ff0000", "#00ff00", "#0000ff"]
    assert _fills(items["split-same-first"][0]) == ["#ffcc00", "#123abc"]
    assert _fills(items["white"][0]) == ["#ffffff"]
    # SnapGene's grey is a colour the file gives, not a missing one.
    assert _fills(items["uncoloured"][0]) == ["#a6acb3"]


def test_every_arrow_is_outlined_and_its_label_text_contrasts_with_its_box(
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
    items = _items(page)
    for arrows, _ in items.values():
        assert all(
            path.attrs["stroke"] not in ("none", path.attrs["fill"])
            for path in arrows.find_all("path")
            if path.attrs["fill"] != "none"
        )
    assert items["white"][1].find_all("text")[0].attrs["fill"] == "#000000"
    assert items["lacZ-alpha"][1].find_all("text")[0].attrs["fill"] == "#ffffff"


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
    record = SequenceRecord(
        "A" * 100, topology="circular", features=(Feature(name, "CDS", (Segment(0, 10),)),)
    )
    _, page = _page(draw_map(record), tmp_path / "map.html")
    [arrows, label] = _items(page)[name]
    assert arrows.attrs["data-name"] == name
    assert not page.find_all("b")
    assert label.find_all("text")[0].text == '<b>"lac" & Z</b>'
