"""`sequence-verify`: one clone's sequencing results held against the product it should be.

`mbio.cli` mounts `sequence_verify` on the root app, as it writes `codon-optimize` there. The
regions judged are the ones `regions` reads off the product, or the features `--feature` names.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Annotated

import typer

from mbio.checks import Check
from mbio.io import read_record
from mbio.verification.judge import Verification, regions, verify
from mbio.verification.result import SequencingResult

#: The extensions a whole-plasmid consensus is read from, as `mbio.io.read_record` reads them.
CONSENSUS_SUFFIXES = frozenset({".dna", ".fa", ".fasta", ".fna", ".gb", ".gbk", ".genbank"})

#: What a consensus read with no per-base support cannot show, printed beside each one: the
#: evidence note's section 3.5.
CONSENSUS_CAVEAT = (
    "a consensus alone cannot show a mixed sample, because the commonest plasmid becomes the "
    "consensus"
)


def read_result(path: Path) -> SequencingResult:
    """Read one sequencing result by its file's extension, named for the file.

    An ``.ab1`` is a Sanger trace. A ``.dna``, GenBank or FASTA file is a whole-plasmid
    consensus, trusted whole.

    Raises
    ------
    NotImplementedError
        For a whole-plasmid folder or a ``.tsv``: the per-base table is not read yet.
    ValueError
        If the extension names neither, or the file does not read as one.
    """
    suffix = path.suffix.lower()
    if path.is_dir() or suffix == ".tsv":
        raise NotImplementedError(
            f"{path.name}: the per-base table is not read yet, so give the consensus FASTA or "
            "GenBank file"
        )
    if suffix == ".ab1":
        from mbio.verification.trace import read_trace

        return read_trace(path)
    if suffix in CONSENSUS_SUFFIXES:
        return SequencingResult(path.name, read_record(path).sequence)
    raise ValueError(
        f"{path.name}: a result is an .ab1 trace, or a consensus as .dna, GenBank or FASTA"
    )


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
            exists=True, help="Its sequencing results: .ab1 traces, or a whole-plasmid consensus."
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
) -> None:
    """Give each junction and insert of PRODUCT a verdict from RESULT, then the clone one."""
    try:
        record = read_record(product)
        judged = regions(record, feature or ())
        read = [read_result(one) for one in results]
        made = verify(record, read, judged)
    except (KeyError, ValueError, NotImplementedError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    for line in report(made, read, len(record)):
        typer.echo(line)
    if not made.verified:
        raise typer.Exit(1)


def report(made: Verification, read: Sequence[SequencingResult], length: int) -> list[str]:
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
    placed = {one.result: one.span is not None for one in made.placements}
    gated = {
        one.name for one in made.result_checks if one.status == "fail" and not placed[one.name]
    }
    width = max((len(one.name) for one in (*made.checks, *made.result_checks)), default=0)
    lines = []
    if len(gated) < len(read):
        lines += [_line(one, width) for one in made.checks]
        lines += [
            f"{one.said(length)}, in {', '.join(one.features) or 'no feature'}: outside every "
            "region, so no verdict"
            for one in made.disagreements
            if not one.regions
        ]
    bare = {one.name for one in read if one.quality is None and one.depth is None}
    for one in made.result_checks:
        lines.append(_line(one, width))
        if one.name in bare - gated:
            lines.append(f"{one.name}: {CONSENSUS_CAVEAT}")
    lines.append("verified" if made.verified else "not verified")
    return lines


def _line(check: Check, width: int) -> str:
    """Return one check as its name, its verdict and its detail, in columns."""
    return f"{check.name:<{width}}  {check.status or 'no verdict':<10}  {check.detail}"
