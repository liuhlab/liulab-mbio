"""The `protocol` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.protocol.model import read_project, read_protocol
from liulab_mbio.protocol.render import PROJECT_DATA_FILE, write_html, write_project_files

app = typer.Typer(help="Render bench protocols.", no_args_is_help=True)


@app.command()
def render(
    source: Annotated[
        Path,
        typer.Argument(exists=True, readable=True, help="Protocol JSON file, or a project folder."),
    ],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="HTML file to write; default: SOURCE with .html."),
    ] = None,
) -> None:
    """Render a protocol JSON file to one page, or a whole project folder to its pages again."""
    folder = _folder(source)
    if folder is not None:
        _render_folder(folder, output)
        return
    try:
        protocol = read_protocol(source)
    except ValueError as error:
        typer.echo(f"error: {source}: {error}", err=True)
        raise typer.Exit(1) from error
    out = output or source.with_suffix(".html")
    # A figure's record sits beside the protocol it was read from, wherever the page is written.
    typer.echo(str(write_html(protocol, out, base=source.parent)))


def _folder(source: Path) -> Path | None:
    """Return the project folder `source` names, or `None` where it names one protocol."""
    if source.is_dir():
        return source
    return source.parent if source.name == PROJECT_DATA_FILE else None


def _render_folder(folder: Path, output: Path | None) -> None:
    """Write every page of a project folder again from the data standing in it."""
    if output is not None:
        typer.echo("error: a project folder is rendered in place, so --output is not one", err=True)
        raise typer.Exit(1)
    data = folder / PROJECT_DATA_FILE
    try:
        project = read_project(data)
    except (OSError, ValueError) as error:
        typer.echo(f"error: {data}: {error}", err=True)
        raise typer.Exit(1) from error
    for path in write_project_files(project, folder).paths:
        typer.echo(str(path))
