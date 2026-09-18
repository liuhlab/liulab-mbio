"""The `library` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.barcodes import SEED
from liulab_mbio.cloning.cli import plan_command
from liulab_mbio.library.plan import NAME_PATTERN, Kind, LibraryPlan, plan_library
from liulab_mbio.library.vector import Site

app = typer.Typer(help="Plan combinatorial protein libraries.", no_args_is_help=True)


@app.command()
def plan(
    parts: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="FASTA of every part list, each record named for the position it fills.",
        ),
    ],
    scheme: Annotated[
        Path,
        typer.Option(
            "--scheme",
            exists=True,
            dir_okay=False,
            readable=True,
            help="Scheme JSON the build is given.",
        ),
    ],
    vector: Annotated[
        Path,
        typer.Option(
            "--vector",
            exists=True,
            dir_okay=False,
            readable=True,
            help="Destination vector sequence file, circular.",
        ),
    ],
    host: Annotated[
        str, typer.Option("--host", help="Codon usage table the coding bases are written for.")
    ],
    coverage: Annotated[
        float,
        typer.Option("--coverage", help="Colonies over products each round is sized for."),
    ],
    out: Annotated[
        Path,
        typer.Option("--out", "-o", file_okay=False, help="Directory to write the outputs into."),
    ],
    kind: Annotated[
        str, typer.Option("--kind", help="What the FASTA holds: protein or dna.")
    ] = "protein",
    site: Annotated[
        str,
        typer.Option(
            "--site",
            help="Feature name, or START-END, to put an internal stuffer at where the vector "
            "carries none.",
        ),
    ] = "",
    pattern: Annotated[
        str,
        typer.Option("--pattern", help="Regex matching a record name, {position} for a position."),
    ] = NAME_PATTERN,
    seed: Annotated[int, typer.Option("--seed", help="Seed the barcodes are drawn with.")] = SEED,
    name: Annotated[str, typer.Option("--name", help="What to call each round's product.")] = "",
) -> None:
    """Plan a library and write the sheets, the records and the protocol into OUT."""
    plan_command(
        lambda: plan_library(
            parts,
            scheme,
            vector,
            host=host,
            coverage=coverage,
            kind=_kind(kind),
            site=_site(site),
            pattern=pattern,
            seed=seed,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: LibraryPlan) -> str:
    """Report the library in one line: what it makes, what it costs, and how it is judged."""
    return (
        f"{made.product.name}: {len(made.product)} bp, {len(made.parts)} part(s) in "
        f"{len(made.part_lists)} list(s), {made.constructs} construct(s), "
        f"{len(made.rounds)} round(s), entry overhangs "
        f"{', '.join(made.standard.entry_overhangs)}, scar {made.standard.scar_overhang}, "
        f"{made.standard.cost} amino acid change(s), checks {made.status}"
    )


def _kind(text: str) -> Kind:
    """Read what the FASTA holds, refusing anything else.

    Raises
    ------
    ValueError
        If it is neither ``protein`` nor ``dna``.
    """
    if text == "protein":
        return "protein"
    if text == "dna":
        return "dna"
    raise ValueError(f"--kind is 'protein' or 'dna', got {text!r}")


def _site(text: str) -> Site | None:
    """Read where a stuffer goes: nothing, a feature name, or a START-END span."""
    if not text:
        return None
    start, sep, end = text.partition("-")
    if sep and start.strip().isdigit() and end.strip().isdigit():
        return int(start), int(end)
    return text
