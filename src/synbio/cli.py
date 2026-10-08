"""The command line: a version, and a sub-app for each pipeline this package holds.

Shaped like `mbio.cli`, because the two are read side by side: one `typer.Typer` named
`app`, `no_args_is_help=True` so a bare invocation prints help instead of nothing, a `version`
command, and a pipeline mounted as a sub-app with `app.add_typer`.
"""

import typer

from synbio import __version__ as _package_version
from synbio.dmx.cli import app as _dmx_app
from synbio.igga.cli import app as _igga_app

#: What `[project.scripts]` registers. Typer builds the parser from the signatures below, so
#: a verb is a function and its help is the docstring.
app = typer.Typer(
    help="Method-specific design pipelines, built on mbio.",
    no_args_is_help=True,
)

app.add_typer(_dmx_app, name="dmx")
app.add_typer(_igga_app, name="igga")


@app.command()
def version() -> None:
    """Print the installed package version."""
    typer.echo(_package_version)
