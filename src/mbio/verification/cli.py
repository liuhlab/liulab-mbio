"""`sequence-verify`: one clone's sequencing results held against the product it should be.

`mbio.cli` mounts `sequence_verify` on the root app, as it writes `codon-optimize` there. The
regions judged are the ones `regions` reads off the product, or the features `--feature` names.
`--out` also writes the clone's page, `draw_page`'s, each trace drawn under its bases.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Annotated

import typer

from mbio.checks import Check
from mbio.io import read_record
from mbio.verification.judge import CONSENSUS_CAVEAT, Verification, regions, verify
from mbio.verification.page import draw_page
from mbio.verification.result import SequencingResult
from mbio.verification.trace import Channels, read_channels, read_trace

#: What `--out` names the clone's page in the directory it is given.
PAGE = "verification.html"

#: The extensions a Sanger trace is written under: both name the one ABIF format.
TRACE_SUFFIXES = (".ab1", ".abi")


def read_result(path: Path) -> SequencingResult:
    """Read one sequencing result by its file's extension, named for the file.

    An ``.ab1`` or ``.abi`` is a Sanger trace. Any file `mbio.io.read_record` reads, a ``.dna``, GenBank or
    FASTA file, is a whole-plasmid consensus, trusted whole.

    Raises
    ------
    NotImplementedError
        For a whole-plasmid folder or a ``.tsv``, whose per-base table this package does not read.
    ValueError
        If the file reads as neither.
    """
    suffix = path.suffix.lower()
    if path.is_dir() or suffix == ".tsv":
        raise NotImplementedError(
            f"{path.name}: the per-base table is not read yet, so give the consensus FASTA or "
            "GenBank file"
        )
    if suffix in TRACE_SUFFIXES:
        return read_trace(path)
    return SequencingResult(path.name, read_record(path).sequence)


def sequence_verify(
    product: Annotated[
        Path,
        typer.Argument(
            exists=True, dir_okay=False, help="The record the clone should be: .dna or GenBank."
        ),
    ],
    results: Annotated[
        list[Path],
        typer.Argument(
            exists=True,
            metavar="RESULT...",
            help="Its sequencing results: .ab1 or .abi traces, or a whole-plasmid consensus.",
        ),
    ],
    feature: Annotated[
        list[str] | None,
        typer.Option(
            "--feature",
            help="A feature to judge in place of the junctions and inserts a plan tagged; once "
            "per feature.",
        ),
    ] = None,
    out: Annotated[
        Path | None,
        typer.Option(
            "--out",
            file_okay=False,
            help=f"A directory to write the clone's page into, as {PAGE}; made if it is not there.",
        ),
    ] = None,
) -> None:
    """Give each junction and insert of PRODUCT a verdict from RESULT, then the clone one."""
    try:
        record = read_record(product)
        judged = regions(record, feature or ())
        given = [read_result(one) for one in results]
        made = verify(record, given, judged)
        channels = _channels(results) if out is not None else []
    except (KeyError, ValueError, NotImplementedError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    for line in report(made, given, len(record)):
        typer.echo(line)
    if out is not None:
        typer.echo(draw_page(record, given, made, channels=channels).write(out / PAGE))
    if not made.verified:
        raise typer.Exit(1)


def _channels(paths: Sequence[Path]) -> list[Channels | None]:
    """Return each result's channels in the order given, ``None`` for one that is no trace."""
    return [
        read_channels(path) if path.suffix.lower() in TRACE_SUFFIXES else None for path in paths
    ]


def report(made: Verification, results: Sequence[SequencingResult], length: int) -> list[str]:
    """Return the lines the verb prints: the regions, then the results, then the clone.

    Each disagreement no region holds is listed after the regions. A result that does not read
    as the product prints its own line and nothing else; where every result is one, so does the
    whole report, but for the clone's line.

    Examples
    --------
    >>> from mbio.sequence import SequenceRecord
    >>> one = SequencingResult("one.fasta", "GATTACAGGCATTCGACCTA")
    >>> made = verify(SequenceRecord(one.bases), [one], [])
    >>> report(made, [one], len(one.bases))[-1]
    'not verified'
    """
    # A check and a placement per result, in the order the results were given, so two results
    # of one name are still told apart.
    gated = [
        check.status == "fail" and placement.span is None
        for check, placement in zip(made.result_checks, made.placements, strict=True)
    ]
    width = max((len(one.name) for one in (*made.checks, *made.result_checks)), default=0)
    lines = []
    if not (gated and all(gated)):
        lines += [_line(one, width) for one in made.checks]
        lines += [
            f"{one.said(length)}, in {', '.join(one.features) or 'no feature'}: outside every "
            "region, so no verdict"
            for one in made.disagreements
            if not one.regions
        ]
    for check, withheld, one in zip(made.result_checks, gated, results, strict=True):
        lines.append(_line(check, width))
        if not withheld and one.quality is None and one.depth is None:
            lines.append(f"{one.name}: {CONSENSUS_CAVEAT}")
    lines.append("verified" if made.verified else "not verified")
    return lines


def _line(check: Check, width: int) -> str:
    """Return one check as its name, its verdict and its detail, in columns."""
    return f"{check.name:<{width}}  {check.status or 'no verdict':<10}  {check.detail}"
