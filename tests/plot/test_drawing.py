"""A record drawn through `draw_map` and written as a page: what the page carries, and refusals."""

import base64
import dataclasses
import re
from pathlib import Path

import pytest

from liulab_mbio.plot import Drawing, draw_map
from liulab_mbio.plot.fonts import BOLD, MONO, SANS
from liulab_mbio.sequence import BindingSite, Feature, Primer, Segment, SequenceRecord, Strand

from ..html import Node, parse


@pytest.fixture(scope="module")
def colour_test(data_dir: Path) -> Drawing:
    """SnapGene's GenBank export of pUC19 with joined, white, grey and origin-crossing features."""
    return draw_map(data_dir / "colour-test-snapgene.gbk")


def _page(drawing: Drawing, path: Path) -> tuple[str, Node]:
    html = drawing.write(path).read_text(encoding="utf-8")
    return html, parse(html)


def _items(page: Node, kind: str = "feature") -> dict[str, list[Node]]:
    """Each item's groups on the page, the arrows' first and the label's after, by name."""
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


def test_the_shipped_unique_cutters_are_labelled_where_snapgene_numbers_their_cuts(
    puc19_file: Path, tmp_path: Path
) -> None:
    _, page = _page(draw_map(puc19_file), tmp_path / "map.html")
    labels = _labels(page, "cut_site")
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
    assert _labels(page, "cut_site") == [
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
    _, page = _page(draw_map(record), tmp_path / "map.html")
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
    ("switches", "drawn"),
    [
        ({}, {"feature": {"CDS", "rep_origin"}, "primer": {"primer"}, "cut_site": {"cut site"}}),
        ({"features": False}, {"primer": {"primer"}, "cut_site": {"cut site"}}),
        ({"primers": False}, {"feature": {"CDS", "rep_origin"}, "cut_site": {"cut site"}}),
        ({"cut_sites": False}, {"feature": {"CDS", "rep_origin"}, "primer": {"primer"}}),
        (
            {"hide_types": ["CDS"], "source": True},
            {"feature": {"rep_origin", "source"}, "primer": {"primer"}, "cut_site": {"cut site"}},
        ),
    ],
)
def test_each_layer_and_feature_type_is_drawn_only_when_asked_for(
    puc19: SequenceRecord, tmp_path: Path, switches: dict, drawn: dict[str, set[str]]
) -> None:
    record = dataclasses.replace(
        puc19,
        features=(
            Feature("pUC19", "source", (Segment(0, len(puc19)),)),
            *(feature for feature in puc19.features if feature.type in ("CDS", "rep_origin")),
        ),
        primers=(
            Primer(
                "M13 fwd",
                "GTAAAACGACGGCCAGT",
                binding_sites=(BindingSite(378, 395, Strand.FORWARD),),
            ),
        ),
    )
    _, page = _page(draw_map(record, **switches), tmp_path / "map.html")
    found: dict[str, set[str]] = {}
    for group in page.find_all("g"):
        if "data-kind" in group.attrs:
            found.setdefault(group.attrs["data-kind"], set()).add(group.attrs["data-type"])
    assert found == drawn
