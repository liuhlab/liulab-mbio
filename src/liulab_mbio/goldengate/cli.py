"""The `goldengate` verbs, mounted on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.goldengate.plan import DEFAULT_HOST, Orientation, Site, plan_assembly
from liulab_mbio.primers import POLYMERASES, Q5, Polymerase

app = typer.Typer(help="Plan Golden Gate assemblies.", no_args_is_help=True)

#: Where a ligase fidelity matrix is read from when `--ligase-matrix` names none. The package
#: ships no such matrix; this points at a copy the user holds.
LIGASE_MATRIX_ENV = "LIULAB_MBIO_LIGASE_MATRIX"


@app.command()
def plan(
    vector: Annotated[
        Path,
        typer.Argument(exists=True, dir_okay=False, readable=True, help="Vector sequence file."),
    ],
    inserts: Annotated[
        list[Path],
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Insert sequence files, in the order they go round the product.",
        ),
    ],
    out: Annotated[
        Path,
        typer.Option(
            "--out", "-o", file_okay=False, help="Directory to write the four outputs into."
        ),
    ],
    site: Annotated[
        str,
        typer.Option(help="Feature name, or START-END, that the inserts replace."),
    ] = "",
    orientation: Annotated[
        list[str] | None,
        typer.Option(help="Which way round an insert goes: forward or reverse; once per insert."),
    ] = None,
    in_frame: Annotated[
        bool,
        typer.Option("--in-frame", help="Hold every insert's junction on a codon boundary."),
    ] = False,
    enzyme: Annotated[
        str, typer.Option(help="Type IIS enzyme to use; the best free one when not given.")
    ] = "",
    ligase_matrix: Annotated[
        Path | None,
        typer.Option(
            "--ligase-matrix",
            envvar=LIGASE_MATRIX_ENV,
            exists=True,
            dir_okay=False,
            readable=True,
            help="A ligase fidelity matrix you hold (.xlsx or .csv), to score the overhangs of "
            "an enzyme nobody has measured.",
        ),
    ] = None,
    prefer_ligase_matrix: Annotated[
        bool,
        typer.Option(
            "--prefer-ligase-matrix",
            help="Score with that matrix even where the enzyme has a measured one.",
        ),
    ] = False,
    polymerase: Annotated[str, typer.Option(help="Polymerase for the two PCRs.")] = Q5.name,
    host: Annotated[str, typer.Option(help="Strain the protocol names.")] = DEFAULT_HOST,
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Plan an assembly and write the product, the primer sheet and the protocol into OUT."""
    try:
        made = plan_assembly(
            vector,
            *inserts,
            site=_site(site),
            orientation=_orientations(orientation, len(inserts)),
            in_frame=in_frame,
            enzyme=enzyme or None,
            profile=ligase_matrix,
            prefer_profile=prefer_ligase_matrix,
            polymerase=_polymerase(polymerase),
            host=host,
            name=name,
        )
        outputs = made.write(out)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    scored = made.overhangs.fidelity
    typer.echo(
        f"{made.product.name}: {len(made.product)} bp, {made.enzyme.name}, "
        f"{len(made.parts)} fragments, overhangs {', '.join(made.overhangs.overhangs)}, "
        f"fidelity {scored.value:.0%} ({scored.label}), checks {made.status}"
    )
    for path in (outputs.product, outputs.primers, outputs.protocol_data, outputs.protocol):
        typer.echo(str(path))


def _site(text: str) -> Site:
    """Read the insertion site: nothing, a feature name, or a START-END span."""
    if not text:
        return None
    start, sep, end = text.partition("-")
    if sep and start.strip().isdigit() and end.strip().isdigit():
        return int(start), int(end)
    return text


def _orientations(given: list[str] | None, count: int) -> tuple[Orientation, ...]:
    """Read one orientation per insert, spreading a single value over them all.

    Raises
    ------
    ValueError
        If a value is neither ``forward`` nor ``reverse``, or there is more than one and not
        one per insert.
    """
    values = given or ["forward"]
    if len(values) == 1:
        values = values * count
    if len(values) != count:
        raise ValueError(f"--orientation given {len(values)} times for {count} insert(s)")
    return tuple(_orientation(text) for text in values)


def _orientation(text: str) -> Orientation:
    """Read the orientation, refusing anything else.

    Raises
    ------
    ValueError
        If it is neither ``forward`` nor ``reverse``.
    """
    if text == "forward":
        return "forward"
    if text == "reverse":
        return "reverse"
    raise ValueError(f"orientation is 'forward' or 'reverse', got {text!r}")


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
