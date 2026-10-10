"""The `gateway` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from mbio.cloning.cli import plan_command
from mbio.cloning.gateway.design import FUSIONS, Fusion
from mbio.cloning.gateway.plan import DEFAULT_HOST, Plan, plan_gateway
from mbio.primers.polymerase import Q5, get_polymerase
from mbio.sequence import position_text

app = typer.Typer(help="Plan Gateway cloning.", no_args_is_help=True)


@app.command()
def plan(
    carrier: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Entry clone, or the insert when --donor is given.",
        ),
    ],
    destination: Annotated[
        Path,
        typer.Argument(
            exists=True, dir_okay=False, readable=True, help="Destination vector sequence file."
        ),
    ],
    out: Annotated[
        Path,
        typer.Option("--out", "-o", file_okay=False, help="Directory to write the outputs into."),
    ],
    donor: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Donor vector, to plan the BP reaction that makes the entry clone.",
        ),
    ] = None,
    amplify: Annotated[
        bool,
        typer.Option(
            "--amplify/--no-amplify",
            help="Design attB primers for an insert carrying no att site, and amplify it.",
        ),
    ] = False,
    fusion: Annotated[
        str, typer.Option(help=f"Tag the insert is read into: {', '.join(FUSIONS)}.")
    ] = "none",
    polymerase: Annotated[str, typer.Option(help="Polymerase for the attB PCR.")] = Q5.name,
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    cleanup_kit: Annotated[
        str,
        typer.Option(
            "--cleanup-kit",
            help="Spin-column kit the protocol names, by catalogue number or by name.",
        ),
    ] = "",
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Plan the Gateway reactions and write the clones and the protocol into OUT."""
    plan_command(
        lambda: plan_gateway(
            carrier,
            destination,
            donor=donor,
            amplify=amplify,
            fusion=_fusion(fusion),
            polymerase=get_polymerase(polymerase),
            host=host,
            cleanup_kit=cleanup_kit,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: Plan) -> str:
    """Report the plan in one line: the route it took, and what stands at each junction."""
    junctions = ", ".join(
        f"{one.name} at {position_text(one.start, len(made.product))}" for one in made.junctions
    )
    return (
        f"{made.product.name}: {len(made.product)} bp, {made.route}, "
        f"{made.lr.recombination.moved.length} bp insert, junctions {junctions}, "
        f"checks {made.status}"
    )


def _fusion(text: str) -> Fusion:
    """Read which tag the insert is read into, whatever the case of the name.

    Raises
    ------
    ValueError
        If it names no fusion this method plans.
    """
    for one in FUSIONS:
        if one.lower() == text.lower():
            return one
    raise ValueError(f"no fusion called {text!r}; name one of {', '.join(FUSIONS)}")
