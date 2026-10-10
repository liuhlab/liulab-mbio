"""The `restriction` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from mbio.bench.materials import DEFAULT_CLEANUP_KIT
from mbio.cloning.cli import CleanupKit, plan_command
from mbio.cloning.restriction.plan import DEFAULT_HOST, Plan, plan_restriction
from mbio.primers.polymerase import Q5, get_polymerase
from mbio.sequence import position_text

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
        typer.Option("--out", "-o", file_okay=False, help="Directory to write the outputs into."),
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
    cleanup_kit: CleanupKit = DEFAULT_CLEANUP_KIT,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Ligate an insert into a vector between two sites, and write the plan into OUT."""
    plan_command(
        lambda: plan_restriction(
            vector,
            insert,
            enzymes=tuple(enzyme or ()),
            polymerase=get_polymerase(polymerase),
            host=host,
            cleanup_kit=cleanup_kit,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: Plan) -> str:
    """Report the cloning in one line: what it cuts, what it makes, and what each junction spells."""
    junctions = ", ".join(
        f"{one.label} at {position_text(one.start, len(made.product))}" for one in made.junctions
    )
    named = ", ".join(one.name for one in made.enzymes)
    weighed = f" chosen over {len(made.refusals)} refused pairs" if made.refusals else ""
    return (
        f"{made.product.name}: {len(made.product)} bp, "
        f"{named}{weighed}, "
        f"{made.insert.length} bp insert into a {made.backbone.length} bp backbone, "
        f"junctions {junctions}, checks {made.status}"
    )
