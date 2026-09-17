"""A record drawn through `draw_map` and written as a page: what the page carries, and refusals."""

import base64
import dataclasses
import re
import struct
from collections.abc import Iterable
from pathlib import Path

import pypdf
import pytest
import vl_convert
from pypdf.generic import DictionaryObject

from liulab_mbio.plot import Drawing, circular, draw_map, linear, svg
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
    to_pdf = vl_convert.svg_to_pdf

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


@pytest.mark.parametrize("dpi", [10**7, 0])
def test_a_png_too_large_or_small_to_draw_is_refused(
    small: Drawing, tmp_path: Path, dpi: float
) -> None:
    with pytest.raises(ValueError, match=r"can be drawn at .* dpi, and a PDF at any size"):
        small.write(tmp_path / "map.png", dpi=dpi)
    assert not (tmp_path / "map.png").exists()


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


def test_every_feature_but_source_is_drawn_with_its_name_on_its_arrow_or_boxed(
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
    items = _items(page)
    names = [f.name for f in colour_test.record.features if f.type != "source"]
    assert sorted(items) == sorted(names)
    assert not [g for g in page.find_all("g") if g.attrs.get("data-type") == "source"]
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
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
    items = _items(page)
    assert _fills(items["split"][0]) == ["#ff0000", "#00ff00", "#0000ff"]
    assert _fills(items["split-same-first"][0]) == ["#ffcc00", "#123abc"]
    assert _fills(items["white"][0]) == ["#ffffff"]
    # SnapGene's grey is a colour the file gives, not a missing one.
    assert _fills(items["uncoloured"][0]) == ["#a6acb3"]


def test_every_arrow_is_outlined_and_every_name_contrasts_with_its_fill(
    colour_test: Drawing, tmp_path: Path
) -> None:
    _, page = _page(colour_test, tmp_path / "map.html")
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
    groups = _items(page)[name]
    assert {group.attrs["data-name"] for group in groups} == {name}
    assert not page.find_all("b")
    texts = [text.text for group in groups for text in group.find_all("text")]
    assert texts == ['<b>"lac" & Z</b>'] * 2


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
    assert _dots(_page(drawing, tmp_path / "map.html")[1]) == dots


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
