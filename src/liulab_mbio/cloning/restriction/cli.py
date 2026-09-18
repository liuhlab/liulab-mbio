"""The `restriction` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.cloning.cli import plan_command
from liulab_mbio.cloning.restriction.plan import DEFAULT_HOST, Plan, plan_restriction

app = typer.Typer(help="Plan restriction and ligation cloning.", no_args_is_help=True)


@app.command()
def plan(
    vector: Annotated[
        Path,
        typer.Argument(exists=True, dir_okay=False, readable=True, help="Vector sequence file."),
    ],
    insert: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="The plasmid the insert is cut out of.",
        ),
    ],
    out: Annotated[
        Path,
        typer.Option(
            "--out", "-o", file_okay=False, help="Directory to write the four outputs into."
        ),
    ],
    enzyme: Annotated[
        list[str] | None,
        typer.Option(help="An enzyme both digests use; once per enzyme, at most twice."),
    ] = None,
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Cut an insert out of one plasmid, ligate it into another, and write the plan into OUT."""
    plan_command(
        lambda: plan_restriction(
            vector,
            insert,
            enzymes=_enzymes(enzyme),
            host=host,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: Plan) -> str:
    """Report the cloning in one line: what it cuts, what it makes, and what each junction spells."""
    junctions = ", ".join(f"{one.label} at {one.start}" for one in made.junctions)
    return (
        f"{made.product.name}: {len(made.product)} bp, "
        f"{', '.join(one.name for one in made.enzymes)}, "
        f"{made.insert.length} bp insert into a {made.backbone.length} bp backbone, "
        f"junctions {junctions}, checks {made.status}"
    )


def _enzymes(given: list[str] | None) -> tuple[str, ...]:
    """Read the enzymes named, refusing a run with none.

    Raises
    ------
    ValueError
        If no enzyme was named.
    """
    if not given:
        raise ValueError("name the enzymes to cut with, once each: --enzyme EcoRI --enzyme BamHI")
    return tuple(given)
