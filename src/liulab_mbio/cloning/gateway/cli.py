"""The `gateway` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.cloning.cli import plan_command
from liulab_mbio.cloning.gateway.plan import DEFAULT_HOST, Plan, plan_gateway

app = typer.Typer(help="Plan Gateway cloning.", no_args_is_help=True)


@app.command()
def plan(
    entry: Annotated[
        Path,
        typer.Argument(
            exists=True, dir_okay=False, readable=True, help="Entry clone sequence file."
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
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Plan an LR reaction and write the expression clone and the protocol into OUT."""
    plan_command(lambda: plan_gateway(entry, destination, host=host, name=name), out, _summary)


def _summary(made: Plan) -> str:
    """Report the reaction in one line: what it makes, and what stands at each junction."""
    junctions = ", ".join(f"{one.name} at {one.start}" for one in made.junctions)
    return (
        f"{made.product.name}: {len(made.product)} bp, LR reaction, "
        f"{made.recombination.moved.length} bp insert, junctions {junctions}, "
        f"checks {made.status}"
    )
