"""The `primers` verbs, mounted on the package command line."""

import re
from pathlib import Path
from typing import Annotated

import typer

from liulab_mbio.bench.oligos import primer_sheet
from liulab_mbio.io import read_record
from liulab_mbio.primers.design import design_pair
from liulab_mbio.primers.evaluation import PairReport, evaluate_pair
from liulab_mbio.primers.genome import GenomeReport, Locus, design_pair_on_genome
from liulab_mbio.primers.polymerase import Q5, Polymerase, get_polymerase
from liulab_mbio.primers.thresholds import TARGET_TM

#: What the verb calls the sheet it orders its oligos from, and the genome report beside it.
SHEET_FILE = "primers.tsv"
GENOME_FILE = "genome.tsv"

#: The columns of the genome report: every amplicon the pair makes on the genome.
GENOME_COLUMNS = ("location", "length", "made_by", "mismatches", "intended")

#: A region as a person types it: `START..END`, 1-based and inclusive, with the sequence it
#: lies on in front of it where the template is a genome.
_REGION = re.compile(r"(?:(?P<name>[^\s:]+):)?\s*(?P<start>\d+)\s*\.\.\s*(?P<end>\d+)\s*")

app = typer.Typer(help="Design PCR primers.", no_args_is_help=True)


@app.command()
def design(
    template: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            readable=True,
            help="Template: a .dna, GenBank or FASTA record, or an indexed genome FASTA "
            "with --assembly.",
        ),
    ],
    out: Annotated[
        Path,
        typer.Option("--out", "-o", file_okay=False, help="Directory to write the outputs into."),
    ],
    region: Annotated[
        str,
        typer.Option(
            "--region",
            help="The amplicon, START..END counted from 1 with both ends included, and "
            "NAME:START..END on a genome. The whole record when not given.",
        ),
    ] = "",
    assembly: Annotated[
        str,
        typer.Option(
            "--assembly",
            help="The genome's name, such as hg38. Given, TEMPLATE is that genome's FASTA and "
            "the pair is designed to be specific on it.",
        ),
    ] = "",
    flank: Annotated[
        int,
        typer.Option("--flank", help="How far either side of a genome region a primer may lie."),
    ] = 200,
    forward_tail: Annotated[
        str, typer.Option("--forward-tail", help="Bases joined to the forward primer's 5' end.")
    ] = "",
    reverse_tail: Annotated[
        str, typer.Option("--reverse-tail", help="Bases joined to the reverse primer's 5' end.")
    ] = "",
    forward_name: Annotated[
        str, typer.Option("--forward-name", help="What to order the forward primer under.")
    ] = "forward",
    reverse_name: Annotated[
        str, typer.Option("--reverse-name", help="What to order the reverse primer under.")
    ] = "reverse",
    target_tm: Annotated[
        float, typer.Option("--target-tm", help="The melting temperature the design aims for, °C.")
    ] = TARGET_TM,
    polymerase: Annotated[
        str, typer.Option("--polymerase", help="Polymerase whose rules judge the pair.")
    ] = Q5.name,
) -> None:
    """Design a primer pair over REGION and write the sheet, and a genome report, into OUT.

    On a genome the pair is searched against that genome and moved off any off-target
    amplicon, so the amplicons it does make are written beside the sheet.
    """
    try:
        chosen = get_polymerase(polymerase)
        tails, names = (forward_tail, reverse_tail), (forward_name, reverse_name)
        written = (
            _on_genome(
                template,
                assembly,
                _locus(region),
                out,
                flank=flank,
                tails=tails,
                names=names,
                target_tm=target_tm,
                polymerase=chosen,
            )
            if assembly
            else _on_template(
                template,
                region,
                out,
                tails=tails,
                names=names,
                target_tm=target_tm,
                polymerase=chosen,
            )
        )
    except (KeyError, ValueError, OSError, RuntimeError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    for path in written:
        typer.echo(str(path))


def _on_template(
    template: Path,
    region: str,
    out: Path,
    *,
    tails: tuple[str, str],
    names: tuple[str, str],
    target_tm: float,
    polymerase: Polymerase,
) -> tuple[Path, ...]:
    """Design on one record, report the pair, and write the sheet."""
    record = read_record(template)
    start, end = _span(region, len(record))
    forward, reverse = design_pair(
        record,
        start,
        end,
        forward_tail=tails[0],
        reverse_tail=tails[1],
        forward_name=names[0],
        reverse_name=names[1],
        target_tm=target_tm,
        polymerase=polymerase,
    )
    report = evaluate_pair(forward, reverse, record, polymerase=polymerase)
    typer.echo(f"{record.name or template.name}: {_said(report)}")
    out.mkdir(parents=True, exist_ok=True)
    sheet = out / SHEET_FILE
    sheet.write_text(primer_sheet([report.forward, report.reverse]), encoding="utf-8")
    return (sheet,)


def _on_genome(
    fasta: Path,
    assembly: str,
    region: Locus,
    out: Path,
    *,
    flank: int,
    tails: tuple[str, str],
    names: tuple[str, str],
    target_tm: float,
    polymerase: Polymerase,
) -> tuple[Path, ...]:
    """Design a pair specific on a genome, report it, and write the sheet and the amplicons."""
    design = design_pair_on_genome(
        fasta,
        assembly,
        region,
        flank=flank,
        forward_tail=tails[0],
        reverse_tail=tails[1],
        forward_name=names[0],
        reverse_name=names[1],
        target_tm=target_tm,
        polymerase=polymerase,
    )
    found = (
        "specific" if design.specific else f"{len(design.genome.off_target)} off-target amplicon(s)"
    )
    near = "" if design.genome.near_matches_checked else ", near matches not checked"
    typer.echo(
        f"{assembly} {_where(region.sequence_name, region.start, region.end)}: {_said(design.pair)}"
    )
    typer.echo(
        f"on {assembly}: {found} after {design.rounds} search(es){near}, "
        f"checks {design.genome.status}"
    )
    out.mkdir(parents=True, exist_ok=True)
    sheet, amplicons = out / SHEET_FILE, out / GENOME_FILE
    sheet.write_text(primer_sheet([design.pair.forward, design.pair.reverse]), encoding="utf-8")
    amplicons.write_text(_genome_table(design.genome), encoding="utf-8")
    return sheet, amplicons


def _said(report: PairReport) -> str:
    """Report the pair in one line: what it makes, how it runs, and how it was judged."""
    length = (
        "no single amplicon" if report.amplicon_length is None else f"{report.amplicon_length} bp"
    )
    return (
        f"{report.forward.primer.sequence} / {report.reverse.primer.sequence}, {length}, "
        f"Ta {report.annealing_temperature:.0f} °C, checks {report.status}"
    )


def _genome_table(report: GenomeReport) -> str:
    """Return every amplicon on the genome, the intended one first, the header alone for none."""
    rows = ["\t".join(GENOME_COLUMNS)]
    rows.extend(
        "\t".join(
            (
                _where(one.sequence_name, one.start, one.end),
                str(one.length),
                one.made_by,
                "/".join(str(count) for count in one.mismatches),
                "yes" if one.intended else "no",
            )
        )
        for one in sorted(report.amplicons, key=lambda one: not one.intended)
    )
    return "\n".join(rows) + "\n"


def _where(name: str, start: int, end: int) -> str:
    """Spell a span as a person reads one: `NAME:START..END`, 1-based and both ends included."""
    return f"{name}:{start + 1}..{end}"


def _span(region: str, length: int) -> tuple[int, int]:
    """Read `--region` on one record, the whole of it where none was given.

    Raises
    ------
    ValueError
        If it is not a span, names a sequence a single record does not have, or lies off the
        record.
    """
    if not region:
        return 0, length
    match = _REGION.fullmatch(region)
    if match is None or match["name"]:
        raise ValueError(f"--region is START..END on a record, got {region!r}")
    first, last = int(match["start"]), int(match["end"])
    if not 1 <= first <= last <= length:
        raise ValueError(f"--region {region}: a span runs from 1 to {length}")
    return first - 1, last


def _locus(region: str) -> Locus:
    """Read `--region` on a genome: the sequence it lies on, and the span.

    Raises
    ------
    ValueError
        If it does not name a sequence and a span running forwards.
    """
    match = _REGION.fullmatch(region)
    if match is None or not match["name"]:
        raise ValueError(f"--region is NAME:START..END on a genome, got {region!r}")
    first, last = int(match["start"]), int(match["end"])
    if not 1 <= first <= last:
        raise ValueError(f"--region {region}: a span runs from 1 up, forwards")
    return Locus(match["name"], first - 1, last)
