"""The `cloning` group, and the spine every method's plan verb follows.

One cloning method is one sub-app under this group, mounted by `liulab_mbio.cli`, so
`cloning --help` lists the methods this package supports and nothing else. `plan_command` is
what each of those verbs repeats: call the planner, write what it planned, turn a refusal into
an error line and exit 1, then print one summary line and one line per file written. A method
supplies its own options and its own summary.

The library pipeline is not a cloning method -- `docs/adr/0004-library-rounds.md` says why --
but its verb is that same spine, so it uses this module too.
"""

from collections.abc import Callable
from pathlib import Path

import typer

from liulab_mbio.cloning.plan import Planned

app = typer.Typer(help="Plan a cloning experiment, one command per method.", no_args_is_help=True)


def plan_command[Plan: Planned](
    planner: Callable[[], Plan], out: Path, summary: Callable[[Plan], str]
) -> None:
    """Run `planner`, write what it planned into `out`, then report the plan and its files.

    The options a verb read are `planner` arguments, so reading one, planning and writing all
    refuse the same way: a `KeyError` or `ValueError` from any of them becomes one error line
    on stderr and exit 1.
    """
    try:
        made = planner()
        written = made.write(out)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    typer.echo(summary(made))
    for path in written.paths:
        typer.echo(str(path))
