"""A project written as a folder of pages: the nav bar, the two columns and the links between.

What each page holds beyond its frame is tested where that content is built; this file tests
the frame and the folder.
"""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from liulab_mbio.cli import app
from liulab_mbio.cloning.plan import write_protocol_files
from liulab_mbio.protocol import (
    INDEX_FILE,
    PROJECT_DATA_FILE,
    REAGENTS_FILE,
    REFERENCES_FILE,
    Figure,
    Item,
    Project,
    Protocol,
    Step,
    page_key,
    write_project_files,
)

from ..html import Node, parse

NAV = ("Overview", "Protocols", "Reagents and equipment", "References")

PLASMID = Item("entry clone", "pENTR carrying the insert")
CELLS = Item("colonies", "transformed cells on an agar plate")


def chain() -> Project:
    """A run of three protocols, each handed what the one before it made."""
    return Project(
        "AP-1 library",
        summary="Build the library in three sittings.",
        inputs=(Item("template", "the plasmid the insert is read from"),),
        protocols=(
            Protocol(
                "BP reaction",
                consumes=(Item("template", "the plasmid the insert is read from"),),
                produces=(PLASMID,),
                steps=(Step("Mix the reaction"), Step("Incubate overnight")),
            ),
            Protocol(
                "Transformation",
                consumes=(PLASMID,),
                produces=(CELLS,),
                steps=(Step("Thaw the cells"),),
            ),
            Protocol("Colony check", consumes=(CELLS,), steps=(Step("Pick ten colonies"),)),
        ),
    )


@pytest.fixture(scope="module")
def folder(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The run of three, written once as a folder."""
    out = tmp_path_factory.mktemp("project")
    write_project_files(chain(), out)
    return out


def pages(folder: Path) -> dict[str, Node]:
    """Every HTML file in the folder, parsed."""
    return {p.name: parse(p.read_text(encoding="utf-8")) for p in sorted(folder.glob("*.html"))}


def test_a_project_writes_its_data_an_index_a_page_per_protocol_and_the_two_shared_pages(
    folder: Path,
) -> None:
    names = sorted(p.name for p in folder.iterdir())
    assert names == [
        "01-bp-reaction.html",
        "02-transformation.html",
        "03-colony-check.html",
        INDEX_FILE,
        PROJECT_DATA_FILE,
        REAGENTS_FILE,
        REFERENCES_FILE,
    ]


def test_the_written_files_are_reported_in_the_order_they_were_written(folder: Path) -> None:
    written = write_project_files(chain(), folder)
    assert [path.name for path in written.paths] == [
        PROJECT_DATA_FILE,
        "01-bp-reaction.html",
        "02-transformation.html",
        "03-colony-check.html",
        INDEX_FILE,
        REAGENTS_FILE,
        REFERENCES_FILE,
    ]


def test_every_link_between_pages_resolves_to_a_file_written_in_the_same_run(
    folder: Path,
) -> None:
    written = {path.name for path in folder.iterdir()}
    for name, page in pages(folder).items():
        for anchor in page.find_all("a"):
            href = anchor.attrs["href"]
            if href.startswith(("#", "http:", "https:", "mailto:")):
                continue
            assert "/" not in href, f"{name} links out of the folder"
            assert href.split("#")[0] in written, f"{name} links to {href}, which nothing wrote"


def test_no_page_in_the_folder_loads_anything_over_the_network(folder: Path) -> None:
    for name, page in pages(folder).items():
        assert not page.find_all("link"), name
        assert not [n for n in page.iter() if "src" in n.attrs], name


def test_the_nav_bar_on_every_page_holds_the_four_items_and_marks_the_page_it_is_on(
    folder: Path,
) -> None:
    for name, page in pages(folder).items():
        [bar] = page.find_all("nav", cls="site")
        links = bar.find_all("a")
        assert tuple(link.text for link in links) == NAV, name
        [current] = [link for link in links if "aria-current" in link.attrs]
        assert current.text == ("Overview" if name == INDEX_FILE else _item(name)), name


def _item(name: str) -> str:
    return {REAGENTS_FILE: "Reagents and equipment", REFERENCES_FILE: "References"}.get(
        name, "Protocols"
    )


def test_every_page_stands_between_two_columns(folder: Path) -> None:
    for name, page in pages(folder).items():
        [frame] = page.find_all("div", cls="frame")
        kinds = [n.attrs.get("class", "") for n in frame.children if isinstance(n, Node)]
        assert len(kinds) == 3, name
        assert "chain" in kinds[0], name
        assert kinds[1] == "page", name
        assert "within" in kinds[2], name


def test_the_left_column_switches_protocol_and_the_right_column_jumps_within_one(
    folder: Path,
) -> None:
    page = pages(folder)["02-transformation.html"]
    [left] = page.find_all("nav", cls="chain")
    assert [a.text for a in left.find_all("a")] == ["BP reaction", "Transformation", "Colony check"]
    assert left.find_all("li", cls="is-here")[0].text == "Transformation"
    [right] = page.find_all("nav", cls="within")
    assert [a.attrs["href"] for a in right.find_all("a")] == ["#step-1"]


def test_a_page_names_its_neighbours_in_words_where_the_page_itself_prints(folder: Path) -> None:
    middle = pages(folder)["02-transformation.html"]
    [main] = middle.find_all("main", cls="page")
    [line] = main.find_all("p", cls="neighbours")
    assert line.text == "Protocol 2 of 3. Comes after BP reaction. Next is Colony check."
    [first] = pages(folder)["01-bp-reaction.html"].find_all("p", cls="neighbours")
    assert first.text == "Protocol 1 of 3. The run starts here. Next is Transformation."
    [last] = pages(folder)["03-colony-check.html"].find_all("p", cls="neighbours")
    assert last.text == "Protocol 3 of 3. Comes after Transformation. The run ends here."


def test_the_index_carries_the_summary_and_every_protocol_with_its_step_count(
    folder: Path,
) -> None:
    index = pages(folder)[INDEX_FILE]
    [main] = index.find_all("main", cls="page")
    assert main.find_all("h1")[0].text == "AP-1 library"
    assert main.find_all("p", cls="summary")[0].text == "Build the library in three sittings."
    [listed] = main.find_all("ol", cls="chain-pages")
    rows = [item.text for item in listed.find_all("li")]
    assert rows == ["BP reaction 2 steps", "Transformation 1 step", "Colony check 1 step"]


def test_each_page_remembers_its_own_checks_under_the_key_the_index_can_read(
    folder: Path,
) -> None:
    keys = {
        name: page.find_all("body")[0].attrs["data-protocol"]
        for name, page in pages(folder).items()
    }
    assert len(set(keys.values())) == len(keys)
    for protocol, name in zip(
        chain().protocols,
        ["01-bp-reaction.html", "02-transformation.html", "03-colony-check.html"],
        strict=True,
    ):
        assert keys[name] == page_key(protocol)


def test_the_two_shared_pages_say_so_where_no_protocol_of_the_run_lists_anything(
    folder: Path,
) -> None:
    for name, heading, empty in (
        (REAGENTS_FILE, "Reagents and equipment", "reagent"),
        (REFERENCES_FILE, "References", "cites a document"),
    ):
        [main] = pages(folder)[name].find_all("main", cls="page")
        assert [n.tag for n in main.children if isinstance(n, Node)] == ["h1", "p"]
        assert main.find_all("h1")[0].text == heading
        assert empty in main.find_all("p")[0].text


def test_a_protocol_written_on_its_own_carries_no_frame(tmp_path: Path) -> None:
    written = write_protocol_files(chain().protocols[0], tmp_path)
    page = parse(written.page.read_text(encoding="utf-8"))
    assert not page.find_all("nav", cls="site")
    assert not page.find_all("div", cls="frame")
    assert page.find_all("nav", cls="toc")


def test_the_cli_renders_a_whole_folder_again_after_an_edit(tmp_path: Path) -> None:
    write_project_files(chain(), tmp_path)
    data = tmp_path / PROJECT_DATA_FILE
    data.write_text(
        data.read_text(encoding="utf-8").replace("Thaw the cells", "Thaw the cells on ice"),
        encoding="utf-8",
    )
    result = CliRunner().invoke(app, ["protocol", "render", str(tmp_path)])
    assert result.exit_code == 0
    assert "Thaw the cells on ice" in (tmp_path / "02-transformation.html").read_text(
        encoding="utf-8"
    )
    assert INDEX_FILE in result.output


def test_the_cli_renders_the_folder_a_project_data_file_stands_in(tmp_path: Path) -> None:
    write_project_files(chain(), tmp_path)
    result = CliRunner().invoke(app, ["protocol", "render", str(tmp_path / PROJECT_DATA_FILE)])
    assert result.exit_code == 0
    assert (tmp_path / INDEX_FILE).exists()


def test_the_cli_refuses_one_output_file_for_a_folder_of_pages(tmp_path: Path) -> None:
    write_project_files(chain(), tmp_path)
    result = CliRunner().invoke(
        app, ["protocol", "render", str(tmp_path), "-o", str(tmp_path / "one.html")]
    )
    assert result.exit_code == 1


def test_a_figure_resolves_beside_the_project_data_file(data_dir: Path, tmp_path: Path) -> None:
    """Every page of a folder is written beside the data, so a figure's record is there too."""
    (tmp_path / "pUC19.dna").write_bytes((data_dir / "pUC19.dna").read_bytes())
    figure = Figure(("pUC19.dna",), "The vector", span=(400, 700))
    one = Protocol("Cut", steps=(Step("Cut the vector", figures=(figure,)),))

    files = write_project_files(Project("Run", protocols=(one,)), tmp_path)

    [page] = files.protocols
    assert "<svg" in page.read_text(encoding="utf-8")
