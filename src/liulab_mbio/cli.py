"""The command line: a version, and a sub-app for each thing the package does.

Typer, because every lab repo that ships a command line uses it: one `typer.Typer` named
`app`, `no_args_is_help=True` so a bare invocation prints help instead of nothing, and a
`version` command. A verb over a pipeline is mounted as a sub-app with `app.add_typer`.
"""

import typer

from liulab_mbio import __version__ as _package_version
from liulab_mbio.goldengate.cli import app as _goldengate_app
from liulab_mbio.protocol.cli import app as _protocol_app

#: What `[project.scripts]` registers. Typer builds the parser from the signatures below, so
#: a verb is a function and its help is the docstring.
app = typer.Typer(
    help="Molecular biology design tools for sequences, enzymes, primers and cloning.",
    no_args_is_help=True,
)

app.add_typer(_goldengate_app, name="goldengate")
app.add_typer(_protocol_app, name="protocol")


@app.command()
def version() -> None:
    """Print the installed package version."""
    typer.echo(_package_version)
