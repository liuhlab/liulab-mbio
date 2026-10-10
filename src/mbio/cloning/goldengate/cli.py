"""The `goldengate` verbs, mounted under the `cloning` group on the package command line."""

from pathlib import Path
from typing import Annotated

import typer

from mbio.cloning.cli import plan_command, read_orientations, read_site
from mbio.cloning.goldengate.plan import DEFAULT_HOST, Plan, plan_assembly
from mbio.codons import DEFAULT_TABLE
from mbio.ligase import LIGASE_MATRIX_ENV
from mbio.primers.polymerase import Q5, get_polymerase

app = typer.Typer(help="Plan Golden Gate assemblies.", no_args_is_help=True)


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
    codon_table: Annotated[
        str,
        typer.Option(
            "--codon-table",
            help="Whose codon usage a proposed domestication picks from; not the strain.",
        ),
    ] = DEFAULT_TABLE,
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
    cleanup_kit: Annotated[
        str,
        typer.Option(
            "--cleanup-kit",
            help="Spin-column kit the protocol names, by catalogue number or by name.",
        ),
    ] = "",
    name: Annotated[str, typer.Option(help="What to call the product.")] = "",
) -> None:
    """Plan an assembly and write the product, the primer sheet and the protocol into OUT."""
    plan_command(
        lambda: plan_assembly(
            vector,
            *inserts,
            site=read_site(site),
            orientation=read_orientations(orientation, len(inserts)),
            in_frame=in_frame,
            enzyme=enzyme or None,
            codon_table=codon_table,
            profile=ligase_matrix,
            prefer_profile=prefer_ligase_matrix,
            polymerase=get_polymerase(polymerase),
            host=host,
            cleanup_kit=cleanup_kit,
            name=name,
        ),
        out,
        _summary,
    )


def _summary(made: Plan) -> str:
    """Report the assembly in one line: what it makes, and how well its set should ligate."""
    scored = made.overhangs.fidelity
    return (
        f"{made.product.name}: {len(made.product)} bp, {made.enzyme.name}, "
        f"{len(made.parts)} fragments, overhangs {', '.join(made.overhangs.overhangs)}, "
        f"fidelity {scored.value:.0%} ({scored.label}), checks {made.status}"
    )
