import json
import re
from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app


def test_the_cli_renders_a_protocol_to_the_file_it_is_given(data_dir: Path, tmp_path: Path) -> None:
    out = tmp_path / "bench.html"
    example = data_dir / "pcr-protocol.json"
    result = CliRunner().invoke(app, ["protocol", "render", str(example), "-o", str(out)])
    assert result.exit_code == 0
    assert "Colony check by PCR" in out.read_text(encoding="utf-8")
    assert out.name in result.output


def test_the_cli_writes_beside_the_protocol_by_default(data_dir: Path, tmp_path: Path) -> None:
    source = tmp_path / "run.json"
    source.write_text(
        (data_dir / "pcr-protocol.json").read_text(encoding="utf-8"), encoding="utf-8"
    )
    result = CliRunner().invoke(app, ["protocol", "render", str(source)])
    assert result.exit_code == 0
    assert (tmp_path / "run.html").exists()


def test_the_cli_says_which_key_it_could_not_read(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text('{"title": "t", "steps": [{"title": "s", "name": "x"}]}', encoding="utf-8")
    result = CliRunner().invoke(app, ["protocol", "render", str(bad)])
    assert result.exit_code == 1
    assert "name" in result.output


def test_the_cli_says_in_one_line_where_a_value_has_the_wrong_type(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(
        '{"title": "t", "steps": [{"title": "s", "instructions": "Mix well."}]}', encoding="utf-8"
    )
    result = CliRunner().invoke(app, ["protocol", "render", str(bad)])
    assert result.exit_code == 1
    [line] = re.sub(r"\x1b\[[0-9;]*m", "", result.output).splitlines()
    assert "protocol.steps[0].instructions" in line


def test_the_cli_refuses_a_protocol_file_that_is_not_there(tmp_path: Path) -> None:
    result = CliRunner().invoke(app, ["protocol", "render", str(tmp_path / "absent.json")])
    assert result.exit_code != 0


def test_the_cli_draws_a_figure_from_beside_the_protocol_it_read(
    data_dir: Path, tmp_path: Path
) -> None:
    """A page written elsewhere still finds the records the protocol was written with."""
    source = tmp_path / "run.json"
    figure = {"records": ["pUC19.dna"], "caption": "The vector", "span": [400, 700]}
    source.write_text(
        json.dumps({"title": "Clone", "steps": [{"title": "Cut", "figures": [figure]}]}),
        encoding="utf-8",
    )
    (tmp_path / "pUC19.dna").write_bytes((data_dir / "pUC19.dna").read_bytes())
    out = tmp_path / "pages" / "bench.html"
    out.parent.mkdir()

    result = CliRunner().invoke(app, ["protocol", "render", str(source), "-o", str(out)])

    assert result.exit_code == 0
    assert "<svg" in out.read_text(encoding="utf-8")
