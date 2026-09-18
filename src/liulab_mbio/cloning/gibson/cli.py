"""The `gibson` verbs, mounted under the `cloning` group on the package command line."""

import typing
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.cloning.cli import plan_command, read_orientations
from liulab_mbio.cloning.gibson.bench import NEBUILDER_HIFI, assembly_product
from liulab_mbio.cloning.gibson.plan import DEFAULT_HOST, Plan, Route, plan_gibson
from liulab_mbio.cloning.plan import Site
from liulab_mbio.primers.polymerase import POLYMERASES, Q5, Polymerase

app = typer.Typer(help="Plan Gibson assemblies.", no_args_is_help=True)


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
            "--out", "-o", file_okay=False, help="Directory to write the four outputs into."
        ),
    ],
    site: Annotated[
        str, typer.Option(help="Feature name, or START-END, that the inserts replace.")
    ] = "",
    orientation: Annotated[
        list[str] | None,
        typer.Option(help="Which way round an insert goes: forward or reverse; once per insert."),
    ] = None,
    route: Annotated[
        list[str] | None,
        typer.Option(help="How an insert is made: amplify or stitch; once per insert."),
    ] = None,
    bridge: Annotated[
        list[str] | None,
        typer.Option(help="A junction joined by one oligo, as BEFORE:AFTER; once per junction."),
    ] = None,
    product: Annotated[
        str, typer.Option(help="Assembly product on the bench.")
    ] = NEBUILDER_HIFI.name,
    polymerase: Annotated[str, typer.Option(help="Polymerase for the PCRs.")] = Q5.name,
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Plan an assembly and write the product, the oligo sheet and the protocol into OUT."""
    plan_command(
        lambda: plan_gibson(
            vector,
            *inserts,
            site=_site(site),
            orientation=read_orientations(orientation, len(inserts)),
            route=_routes(route, len(inserts)),
            bridge=_bridges(bridge),
            product=assembly_product(product),
            polymerase=_polymerase(polymerase),
            host=host,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: Plan) -> str:
    """Report the assembly in one line: what it makes, and what joins it."""
    overlaps = ", ".join(f"{len(one)} bp" for one in made.overlaps)
    return (
        f"{made.plasmid.name}: {len(made.plasmid)} bp, {len(made.parts)} fragments, "
        f"overlaps {overlaps}, {made.product.name}, checks {made.status}"
    )


def _site(text: str) -> Site:
    """Read the insertion site: nothing, a feature name, or a START-END span."""
    if not text:
        return None
    start, sep, end = text.partition("-")
    if sep and start.strip().isdigit() and end.strip().isdigit():
        return int(start), int(end)
    return text


def _routes(given: Sequence[str] | None, count: int) -> tuple[Route, ...]:
    """Read one route per insert, spreading a single value over them all.

    Raises
    ------
    ValueError
        If there is more than one and not one per insert. `plan_gibson` refuses an unknown one.

    Examples
    --------
    >>> _routes(["stitch"], 2)
    ('stitch', 'stitch')
    """
    values = list(given) if given else ["amplify"]
    if len(values) == 1:
        values = values * count
    if len(values) != count:
        raise ValueError(f"--route given {len(values)} times for {count} insert(s)")
    return tuple(typing.cast("Route", text.strip().lower()) for text in values)


def _bridges(given: Sequence[str] | None) -> tuple[tuple[str, str], ...]:
    """Read each bridged junction, written as the two parts it joins with a colon between.

    Raises
    ------
    ValueError
        If one does not name two parts. `plan_gibson` refuses a junction it has not got.

    Examples
    --------
    >>> _bridges(["GFP:pUC19 backbone"])
    (('GFP', 'pUC19 backbone'),)
    """
    pairs = []
    for text in given or ():
        before, sep, after = text.partition(":")
        if not (sep and before.strip() and after.strip()):
            raise ValueError(
                f"--bridge names the two parts a junction joins, as BEFORE:AFTER, got {text!r}"
            )
        pairs.append((before.strip(), after.strip()))
    return tuple(pairs)


def _polymerase(name: str) -> Polymerase:
    """Read one of the polymerases this package ships, whatever the case of its name.

    Raises
    ------
    ValueError
        If it ships no polymerase of that name.
    """
    found = next((one for one in POLYMERASES if one.name.lower() == name.lower()), None)
    if found is None:
        shipped = ", ".join(one.name for one in POLYMERASES)
        raise ValueError(f"no polymerase called {name!r}; this package ships {shipped}")
    return found
