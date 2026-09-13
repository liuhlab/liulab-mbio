from pathlib import Path

from typer.testing import CliRunner

from liulab_mbio.cli import app

EXAMPLE = Path(__file__).parents[1] / "data" / "pcr-protocol.json"


def test_the_cli_renders_a_protocol_to_the_file_it_is_given(tmp_path: Path) -> None:
    out = tmp_path / "bench.html"
    result = CliRunner().invoke(app, ["protocol", "render", str(EXAMPLE), "-o", str(out)])
    assert result.exit_code == 0
    assert "Colony check by PCR" in out.read_text(encoding="utf-8")
    assert out.name in result.output


def test_the_cli_writes_beside_the_protocol_by_default(tmp_path: Path) -> None:
    source = tmp_path / "run.json"
    source.write_text(EXAMPLE.read_text(encoding="utf-8"), encoding="utf-8")
    result = CliRunner().invoke(app, ["protocol", "render", str(source)])
    assert result.exit_code == 0
    assert (tmp_path / "run.html").exists()


def test_the_cli_says_which_key_it_could_not_read(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text('{"title": "t", "steps": [{"name": "no title here"}]}', encoding="utf-8")
    result = CliRunner().invoke(app, ["protocol", "render", str(bad)])
    assert result.exit_code == 1
    assert "name" in result.output


def test_the_cli_refuses_a_protocol_file_that_is_not_there(tmp_path: Path) -> None:
    result = CliRunner().invoke(app, ["protocol", "render", str(tmp_path / "absent.json")])
    assert result.exit_code != 0
