"""The `goldengate` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.goldengate.plan import DEFAULT_HOST, Orientation, Site, plan_assembly
from liulab_mbio.primers import ONETAQ, PHUSION, Q5, TAQ, Polymerase

app = typer.Typer(help="Plan Golden Gate assemblies.", no_args_is_help=True)

#: The polymerases a caller can name on the command line.
POLYMERASES: dict[str, Polymerase] = {one.name.lower(): one for one in (Q5, PHUSION, TAQ, ONETAQ)}


@app.command()
def plan(
    vector: Annotated[
        Path,
        typer.Argument(exists=True, dir_okay=False, readable=True, help="Vector sequence file."),
    ],
    inserts: Annotated[
        list[Path],
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Insert sequence files, in the order they go round the product.",
        ),
    ],
    out: Annotated[
        Path,
        typer.Option(
            "--out", "-o", file_okay=False, help="Directory to write the three outputs into."
        ),
    ],
    site: Annotated[
        str,
        typer.Option(help="Feature name, or START-END, that the inserts replace."),
    ] = "",
    orientation: Annotated[
        list[str] | None,
        typer.Option(help="Which way round an insert goes: forward or reverse; once per insert."),
    ] = None,
    in_frame: Annotated[
        bool,
        typer.Option("--in-frame", help="Hold every insert's junction on a codon boundary."),
    ] = False,
    enzyme: Annotated[
        str, typer.Option(help="Type IIS enzyme to use; the best free one when not given.")
    ] = "",
    polymerase: Annotated[str, typer.Option(help="Polymerase for the two PCRs.")] = Q5.name,
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Plan an assembly and write the product, the primer sheet and the protocol into OUT."""
    try:
        made = plan_assembly(
            vector,
            *inserts,
            site=_site(site),
            orientation=_orientations(orientation, len(inserts)),
            in_frame=in_frame,
            enzyme=enzyme or None,
            polymerase=_polymerase(polymerase),
            host=host,
            name=name,
        )
        outputs = made.write(out)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    typer.echo(
        f"{made.product.name}: {len(made.product)} bp, {made.enzyme.name}, "
        f"{len(made.parts)} fragments, overhangs {', '.join(made.overhangs.overhangs)}, "
        f"checks {made.status}"
    )
    for path in (outputs.product, outputs.primers, outputs.protocol):
        typer.echo(str(path))


def _site(text: str) -> Site:
    """Read the insertion site: nothing, a feature name, or a START-END span."""
    if not text:
        return None
    start, sep, end = text.partition("-")
    if sep and start.strip().isdigit() and end.strip().isdigit():
        return int(start), int(end)
    return text


def _orientations(given: list[str] | None, count: int) -> tuple[Orientation, ...]:
    """Read one orientation per insert, spreading a single value over them all.

    Raises
    ------
    ValueError
        If a value is neither ``forward`` nor ``reverse``, or there is more than one and not
        one per insert.
    """
    values = given or ["forward"]
    if len(values) == 1:
        values = values * count
    if len(values) != count:
        raise ValueError(f"--orientation given {len(values)} times for {count} insert(s)")
    return tuple(_orientation(text) for text in values)


def _orientation(text: str) -> Orientation:
    """Read the orientation, refusing anything else.

    Raises
    ------
    ValueError
        If it is neither ``forward`` nor ``reverse``.
    """
    if text == "forward":
        return "forward"
    if text == "reverse":
        return "reverse"
    raise ValueError(f"orientation is 'forward' or 'reverse', got {text!r}")


def _polymerase(name: str) -> Polymerase:
    """Read one of the polymerases this package ships.

    Raises
    ------
    ValueError
        If it ships no polymerase of that name.
    """
    found = POLYMERASES.get(name.lower())
    if found is None:
        shipped = ", ".join(one.name for one in POLYMERASES.values())
        raise ValueError(f"no polymerase called {name!r}; this package ships {shipped}")
    return found
