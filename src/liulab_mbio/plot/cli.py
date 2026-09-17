"""The `plot` verbs, mounted on the package command line."""

import re
from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.io import read_record
from liulab_mbio.plot.drawing import Region, draw_map
from liulab_mbio.plot.layers import notice
from liulab_mbio.sequence import SequenceRecord

#: A span as a person types it: `START..END`, 1-based and inclusive, as the map prints one.
_SPAN = re.compile(r"\s*(\d+)\s*\.\.\s*(\d+)\s*")

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
    region: Annotated[
        str | None,
        typer.Option(
            help="A feature's name, or START..END counted from 1 with both ends included, to draw "
            "as a line numbered as the record is. END before START runs across the origin.",
        ),
    ] = None,
    linear: Annotated[
        bool, typer.Option("--linear", help="Draw a circular record opened as a line.")
    ] = False,
    sequence_view: Annotated[
        bool,
        typer.Option(
            "--sequence-view", help="Draw the bases beside the map, at most 100,000 of them."
        ),
    ] = False,
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
    bases_per_row: Annotated[
        int, typer.Option(help="How many bases a row of the sequence view holds.")
    ] = 60,
    one_strand: Annotated[
        bool, typer.Option("--one-strand", help="Draw only the top strand in the sequence view.")
    ] = False,
    dpi: Annotated[float, typer.Option(help="Resolution of a PNG, in dots per inch.")] = 300,
) -> None:
    """Draw RECORD as a map and write it to each file named, printing each file written.

    Then prints how many labels the map hid for want of room, if any did.
    """
    try:
        read = read_record(record)
        drawing = draw_map(
            read,
            region=_region(region, read),
            linear=linear,
            sequence_view=sequence_view,
            features=not no_features,
            primers=not no_primers,
            cut_sites=not no_cut_sites,
            enzymes=enzyme,
            hide_types=hide_type or (),
            source=source,
            bases_per_row=bases_per_row,
            both_strands=not one_strand,
        )
        for out in output:
            typer.echo(str(drawing.write(out, dpi=dpi)))
        if said := notice(drawing.hidden):
            typer.echo(said)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error


def _region(text: str | None, record: SequenceRecord) -> Region | None:
    """Read `--region`: a feature's name, or `START..END` as the span `draw_map` takes.

    Text of that form is always a span, even where a feature is named so.

    Raises
    ------
    ValueError
        If a position lies off the record, or a span runs across the origin of a linear record.
    """
    match = None if text is None else _SPAN.fullmatch(text)
    if match is None:
        return text
    first, last, length = int(match[1]), int(match[2]), len(record)
    if not (1 <= first <= length and 1 <= last <= length):
        raise ValueError(f"--region {text}: a position runs from 1 to {length}")
    if last >= first:
        return first - 1, last
    if record.topology == "linear":
        raise ValueError(
            f"--region {text} runs across the origin, which the linear record {record.name!r} "
            "does not have"
        )
    return first - 1, last + length
