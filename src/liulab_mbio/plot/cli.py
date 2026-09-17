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
        list[Path],
        typer.Option(
            "--output",
            "-o",
            dir_okay=False,
            help="File to write: .html, .png or .pdf. Repeat it to write several.",
        ),
    ],
    no_features: Annotated[
        bool, typer.Option("--no-features", help="Leave the features off.")
    ] = False,
    no_primers: Annotated[
        bool, typer.Option("--no-primers", help="Leave the primers off.")
    ] = False,
    no_cut_sites: Annotated[
        bool, typer.Option("--no-cut-sites", help="Leave the cut sites off.")
    ] = False,
    enzyme: Annotated[
        list[str] | None,
        typer.Option(
            "--enzyme",
            help="An enzyme whose every cut site to draw; once per enzyme. The shipped enzymes "
            "that cut once when not given.",
        ),
    ] = None,
    hide_type: Annotated[
        list[str] | None,
        typer.Option("--hide-type", help="A feature type to leave off; once per type."),
    ] = None,
    source: Annotated[bool, typer.Option("--source", help="Draw the source feature.")] = False,
    dpi: Annotated[float, typer.Option(help="Resolution of a PNG, in dots per inch.")] = 300,
) -> None:
    """Draw RECORD as a map and write it to each file named, printing each file written."""
    try:
        drawing = draw_map(
            record,
            features=not no_features,
            primers=not no_primers,
            cut_sites=not no_cut_sites,
            enzymes=enzyme,
            hide_types=hide_type or (),
            source=source,
        )
        for out in output:
            typer.echo(str(drawing.write(out, dpi=dpi)))
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
