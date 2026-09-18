"""The `restriction` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.cloning.cli import plan_command
from liulab_mbio.cloning.restriction.plan import DEFAULT_HOST, Plan, plan_restriction
from liulab_mbio.primers.polymerase import POLYMERASES, Q5, Polymerase

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
            help="The insert, or the plasmid it is cut out of.",
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
        typer.Option(
            help="An enzyme both digests use; once per enzyme, at most twice. "
            "Chosen for you when you name none."
        ),
    ] = None,
    polymerase: Annotated[
        str, typer.Option(help="Polymerase for the insert's PCR, where one is run.")
    ] = Q5.name,
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Ligate an insert into a vector between two sites, and write the plan into OUT."""
    plan_command(
        lambda: plan_restriction(
            vector,
            insert,
            enzymes=tuple(enzyme or ()),
            polymerase=_polymerase(polymerase),
            host=host,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: Plan) -> str:
    """Report the cloning in one line: what it cuts, what it makes, and what each junction spells."""
    junctions = ", ".join(f"{one.label} at {one.start}" for one in made.junctions)
    named = ", ".join(one.name for one in made.enzymes)
    weighed = f" chosen over {len(made.refusals)} refused pairs" if made.refusals else ""
    return (
        f"{made.product.name}: {len(made.product)} bp, "
        f"{named}{weighed}, "
        f"{made.insert.length} bp insert into a {made.backbone.length} bp backbone, "
        f"junctions {junctions}, checks {made.status}"
    )


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
