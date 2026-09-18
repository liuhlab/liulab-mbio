"""The `gibson` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.cloning.cli import plan_command
from liulab_mbio.cloning.gibson.plan import DEFAULT_HOST, Plan, plan_gibson
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
        typer.Argument(exists=True, dir_okay=False, readable=True, help="Insert sequence file."),
    ],
    out: Annotated[
        Path,
        typer.Option(
            "--out", "-o", file_okay=False, help="Directory to write the four outputs into."
        ),
    ],
    site: Annotated[
        str, typer.Option(help="Feature name, or START-END, that the insert replaces.")
    ] = "",
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
