"""The `dmx` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from mbio.checks import counted
from mbio.cloning.cli import plan_command
from synbio.dmx.plan import ReadBackPlan, plan_dmx

app = typer.Typer(help="Read designs back one well at a time.", no_args_is_help=True)


@app.command()
def plan(
    build: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Build JSON: the designs sheet, the archive plate, the route and the floor.",
        ),
    ],
    out: Annotated[
        Path,
        typer.Option("--out", "-o", file_okay=False, help="Directory to write the outputs into."),
    ],
) -> None:
    """Plan the read-back BUILD asks for, and write its protocol into OUT."""
    plan_command(lambda: plan_dmx(build), out, _summary)


def _summary(made: ReadBackPlan) -> str:
    """Report the read-back in one line: what is read, how, and what it takes at the bench."""
    one = made.validation
    return (
        f"{made.build.name}: {len(one.designs)} of "
        f"{counted(len(made.designs), 'design')} read back on {one.route.name}, "
        f"{counted(one.wells, 'well')} over {counted(len(one.picked), 'picked plate')}"
    )
