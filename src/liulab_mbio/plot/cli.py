"""The `plot` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.plot.drawing import draw_map

app = typer.Typer(help="Draw sequence records.", no_args_is_help=True)


@app.command("map")
def map_(
    record: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Sequence file: .dna, GenBank or FASTA.",
        ),
    ],
    output: Annotated[
        Path,
        typer.Option("--output", "-o", dir_okay=False, help="File to write: .html."),
    ],
) -> None:
    """Draw RECORD as a map and write it to the file named, printing the file written."""
    try:
        written = draw_map(record).write(output)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    typer.echo(str(written))
