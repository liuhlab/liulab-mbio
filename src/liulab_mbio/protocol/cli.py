"""The `protocol` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.protocol.model import read_protocol
from liulab_mbio.protocol.render import write_html

app = typer.Typer(help="Render bench protocols.", no_args_is_help=True)


@app.command()
def render(
    source: Annotated[
        Path,
        typer.Argument(exists=True, dir_okay=False, readable=True, help="Protocol JSON file."),
    ],
    output: Annotated[
        Path | None,
        typer.Option("--output", "-o", help="HTML file to write; default: SOURCE with .html."),
    ] = None,
) -> None:
    """Render a protocol JSON file to one self-contained HTML file."""
    try:
        protocol = read_protocol(source)
    except ValueError as error:
        typer.echo(f"error: {source}: {error}", err=True)
        raise typer.Exit(1) from error
    typer.echo(str(write_html(protocol, output or source.with_suffix(".html"))))
