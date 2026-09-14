"""The command line: a version, a codon-optimise verb, and a sub-app for each pipeline.

Typer, because every lab repo that ships a command line uses it: one `typer.Typer` named
`app`, `no_args_is_help=True` so a bare invocation prints help instead of nothing, and a
`version` command. A verb over a pipeline is mounted as a sub-app with `app.add_typer`.
"""

from typing import Annotated

import typer

from liulab_mbio import __version__ as _package_version
from liulab_mbio.goldengate.cli import app as _goldengate_app
from liulab_mbio.protocol.cli import app as _protocol_app
from liulab_mbio.translate import (
    DEFAULT_NAME,
    CodingSequence,
    optimize_coding_sequence,
    optimize_protein,
)

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


@app.command("codon-optimize")
def codon_optimize(
    sequence: Annotated[
        str, typer.Argument(help="The protein or the coding sequence, as letters.")
    ],
    kind: Annotated[str, typer.Option("--kind", help="What SEQUENCE is: protein or dna.")],
    host: Annotated[str, typer.Option("--host", help="Codon usage table to write for.")],
    forbid: Annotated[
        list[str] | None,
        typer.Option("--forbid", help="An enzyme whose site must not appear; once per enzyme."),
    ] = None,
    name: Annotated[str, typer.Option("--name", help="What to call the sequence.")] = DEFAULT_NAME,
) -> None:
    """Write a protein as DNA for a host, or check a coded one, free of forbidden sites."""
    try:
        gene = _coded(sequence, kind, host=host, forbidden=forbid or [], name=name)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    typer.echo(
        f"{gene.name}: {len(gene.dna)} bp, {len(gene.protein)} aa, host {gene.host}, "
        f"{len(gene.changes)} codon change(s), free of "
        f"{', '.join(gene.forbidden) or 'nothing named'}"
    )
    for change in gene.changes:
        typer.echo(
            f"codon {change.codon_index} {change.old_codon} -> {change.new_codon} "
            f"({change.amino_acid}), {change.site.enzyme.name} site removed"
        )
    typer.echo(gene.dna)


def _coded(
    sequence: str, kind: str, *, host: str, forbidden: list[str], name: str
) -> CodingSequence:
    """Write or check one sequence, whichever `kind` names.

    Raises
    ------
    ValueError
        If `kind` is neither ``protein`` nor ``dna``.
    """
    if kind == "protein":
        return optimize_protein(sequence, host=host, forbidden=forbidden, name=name)
    if kind == "dna":
        return optimize_coding_sequence(sequence, host=host, forbidden=forbidden, name=name)
    raise ValueError(f"--kind is 'protein' or 'dna', got {kind!r}")
