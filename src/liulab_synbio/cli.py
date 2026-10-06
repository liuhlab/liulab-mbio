"""The command line: a version, and a sub-app for each pipeline this package holds.

Shaped like `liulab_mbio.cli`, because the two are read side by side: one `typer.Typer` named
`app`, `no_args_is_help=True` so a bare invocation prints help instead of nothing, a `version`
command, and a pipeline mounted as a sub-app with `app.add_typer`.
"""

import typer

from liulab_synbio import __version__ as _package_version

#: What `[project.scripts]` registers. Typer builds the parser from the signatures below, so
#: a verb is a function and its help is the docstring.
app = typer.Typer(
    help="Method-specific design pipelines, built on liulab_mbio.",
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Keep this a group of verbs.

    Typer collapses an app holding one command and no callback into that command alone, so
    without this `version` stops being a verb and a bare invocation runs it. The callback goes
    when a second pipeline mounts here, or stays if it grows an option every verb shares.
    """


@app.command()
def version() -> None:
    """Print the installed package version."""
    typer.echo(_package_version)
