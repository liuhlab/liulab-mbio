"""The `igga` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.bench.prices import PRICES_ENV
from liulab_mbio.cloning.cli import plan_command
from liulab_mbio.ligase import LIGASE_MATRIX_ENV
from liulab_synbio.igga.plan import NAME_PATTERN, Kind, LibraryPlan, Site, plan_igga

app = typer.Typer(help="Plan combinatorial protein libraries.", no_args_is_help=True)


@app.command()
def plan(
    build: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="The build, as JSON: the positions, the parts and vector files, and the dials.",
        ),
    ],
    out: Annotated[
        Path,
        typer.Option("--out", "-o", file_okay=False, help="Directory to write the outputs into."),
    ],
    kind: Annotated[
        str, typer.Option("--kind", help="What the parts FASTA holds: protein or dna.")
    ] = "protein",
    site: Annotated[
        str,
        typer.Option(
            "--site",
            help="Feature name, or START-END, to put an internal stuffer at where the vector "
            "carries none.",
        ),
    ] = "",
    working_site: Annotated[
        str,
        typer.Option(
            "--working-site",
            help="Feature name, or START-END, to put the ccdB cassette at where the working "
            "vector carries none.",
        ),
    ] = "",
    pattern: Annotated[
        str,
        typer.Option("--pattern", help="Regex matching a record name, {position} for a position."),
    ] = NAME_PATTERN,
    prices: Annotated[
        Path | None,
        typer.Option(
            "--prices",
            envvar=PRICES_ENV,
            exists=True,
            dir_okay=False,
            readable=True,
            help="A banded price record you hold (.csv), to cost the protocol's bill. The "
            "quantities compute without one; the money is then a hole.",
        ),
    ] = None,
    ligase_matrix: Annotated[
        Path | None,
        typer.Option(
            "--ligase-matrix",
            envvar=LIGASE_MATRIX_ENV,
            exists=True,
            dir_okay=False,
            readable=True,
            help="A ligase fidelity matrix you hold (.xlsx or .csv), to report how often the "
            "ligase it measured joins each round's overhangs. Nothing is designed on it.",
        ),
    ] = None,
    ligase_sheet: Annotated[
        str,
        typer.Option(
            "--ligase-sheet",
            help="Which sheet of that matrix to read, for a workbook holding several.",
        ),
    ] = "",
) -> None:
    """Plan the library BUILD asks for, and write its sheets, records and protocol into OUT."""
    plan_command(
        lambda: plan_igga(
            build,
            kind=_kind(kind),
            site=_site(site),
            working_site=_site(working_site),
            pattern=pattern,
            prices=prices,
            profile=ligase_matrix,
            profile_sheet=ligase_sheet or None,
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
    """Read where a cassette goes: nothing, a feature name, or a START-END span."""
    if not text:
        return None
    start, sep, end = text.partition("-")
    if sep and start.strip().isdigit() and end.strip().isdigit():
        return int(start), int(end)
    return text
