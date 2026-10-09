"""The command line: a version, the verbs over a flat module, and a sub-app for each feature.

Typer, because every lab repo that ships a command line uses it: one `typer.Typer` named
`app`, `no_args_is_help=True` so a bare invocation prints help instead of nothing, and a
`version` command. A feature with a package of its own mounts its sub-app here with
`app.add_typer`, and one cloning method is mounted under the `cloning` group, so
`cloning --help` lists the methods that exist. A feature that is one flat module -- `translate`
and `barcodes` -- has its verb here instead, because the command line is the layer that may
import typer.

A verb writing more than one file takes `--out DIR` and prints the summary and then each file
written, as every pipeline's plan verb does.
"""

import re
from pathlib import Path
from typing import Annotated

import typer

from mbio import __version__ as _package_version
from mbio.barcodes import (
    GC_BAND,
    MAX_HOMOPOLYMER,
    METRIC,
    MIN_DISTANCE,
    SEED,
    BarcodeRules,
    Metric,
    check_barcodes,
    deletion_ambiguity,
    design_barcodes,
    separation,
)
from mbio.checks import Check, worst_of
from mbio.cloning.cli import app as _cloning_app
from mbio.cloning.gateway.cli import app as _gateway_app
from mbio.cloning.gibson.cli import app as _gibson_app
from mbio.cloning.goldengate.cli import app as _goldengate_app
from mbio.cloning.restriction.cli import app as _restriction_app
from mbio.plot.cli import app as _plot_app
from mbio.primers.cli import app as _primers_app
from mbio.protocol.cli import app as _protocol_app
from mbio.snapgene import write_dna
from mbio.translate import (
    DEFAULT_NAME,
    CodingSequence,
    coding_record,
    optimize_coding_sequence,
    optimize_protein,
)

#: What `[project.scripts]` registers. Typer builds the parser from the signatures below, so
#: a verb is a function and its help is the docstring.
app = typer.Typer(
    help="Molecular biology design tools for sequences, enzymes, primers and cloning.",
    no_args_is_help=True,
)

#: What `codon-optimize --out` calls the coding sequence it wrote, and the codons it moved.
CODING_FILE = "coding-sequence.dna"
CHANGE_FILE = "codon-changes.tsv"

#: The columns of the codon change table.
CHANGE_COLUMNS = ("codon_index", "old_codon", "new_codon", "amino_acid", "enzyme")

#: What `barcode-design` calls the set it designed, and the checks on it.
BARCODE_FILE = "barcodes.tsv"
CHECK_FILE = "checks.tsv"

#: The columns of the barcode set, and of a check table.
BARCODE_COLUMNS = ("index", "barcode")
CHECK_COLUMNS = ("name", "status", "value", "detail")

#: A GC band as a person types it: `LOW-HIGH`, each a share of 1.
_BAND = re.compile(r"\s*(\d*\.?\d+)\s*-\s*(\d*\.?\d+)\s*")

_cloning_app.add_typer(_gateway_app, name="gateway")
_cloning_app.add_typer(_gibson_app, name="gibson")
_cloning_app.add_typer(_goldengate_app, name="goldengate")
_cloning_app.add_typer(_restriction_app, name="restriction")

app.add_typer(_cloning_app, name="cloning")
app.add_typer(_protocol_app, name="protocol")
app.add_typer(_plot_app, name="plot")
app.add_typer(_primers_app, name="primers")


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
    out: Annotated[
        Path | None,
        typer.Option(
            "--out",
            "-o",
            file_okay=False,
            help="Directory to write the coding sequence and the codon changes into; "
            "printed instead when it is not given.",
        ),
    ] = None,
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
    if out is None:
        for change in gene.changes:
            typer.echo(
                f"codon {change.codon_index} {change.old_codon} -> {change.new_codon} "
                f"({change.amino_acid}), {change.site.enzyme.name} site removed"
            )
        typer.echo(gene.dna)
        return
    for path in _write_gene(gene, out):
        typer.echo(str(path))


@app.command("barcode-design")
def barcode_design(
    count: Annotated[int, typer.Argument(help="How many barcodes the part list needs.")],
    length: Annotated[int, typer.Option("--length", help="How many bases name one part.")],
    out: Annotated[
        Path,
        typer.Option(
            "--out", "-o", file_okay=False, help="Directory to write the set and its checks into."
        ),
    ],
    scar: Annotated[
        str, typer.Option("--scar", help="The cloning scar joining one barcode to the next.")
    ] = "",
    phase: Annotated[
        int,
        typer.Option(
            "--phase", help="Bases of the barcode's first codon read before it: 0, 1 or 2."
        ),
    ] = 0,
    untranslated: Annotated[
        bool,
        typer.Option(
            "--untranslated", help="The construct never translates the barcode: no frame rules."
        ),
    ] = False,
    forbid: Annotated[
        list[str] | None,
        typer.Option("--forbid", help="An enzyme whose site must not appear; once per enzyme."),
    ] = None,
    distance: Annotated[
        int, typer.Option("--distance", help="The fewest units two barcodes stand apart.")
    ] = MIN_DISTANCE,
    metric: Annotated[str, typer.Option("--metric", help="How a distance is counted.")] = METRIC,
    max_homopolymer: Annotated[
        int,
        typer.Option(
            "--max-homopolymer", help="The longest run of one base allowed; 0 for no cap."
        ),
    ] = MAX_HOMOPOLYMER,
    gc_band: Annotated[
        str,
        typer.Option(
            "--gc-band", help="The share of G and C allowed, as LOW-HIGH; none by default."
        ),
    ] = "",
    seed: Annotated[int, typer.Option("--seed", help="The seed the draw is made with.")] = SEED,
) -> None:
    """Design one part list's barcodes and write the set and its checks into OUT."""
    try:
        rules = BarcodeRules(
            length,
            scar=scar,
            phase=None if untranslated else phase,
            forbidden=forbid or (),
            distance=distance,
            max_homopolymer=max_homopolymer or None,
            gc_band=_band(gc_band),
            metric=_metric(metric),
        )
        barcodes = design_barcodes(count, rules, seed=seed)
    except (KeyError, ValueError) as error:
        typer.echo(f"error: {error}", err=True)
        raise typer.Exit(1) from error
    checks = _barcode_checks(barcodes, rules)
    typer.echo(
        f"{len(barcodes)} barcode(s) of {rules.length} bases, scar {rules.scar or 'none'}, "
        f"{rules.distance} apart by {rules.metric}, free of "
        f"{', '.join(one.name for one in rules.enzymes) or 'nothing named'}, seed {seed}, "
        f"checks {worst_of(checks)}"
    )
    out.mkdir(parents=True, exist_ok=True)
    set_file, check_file = out / BARCODE_FILE, out / CHECK_FILE
    set_file.write_text(_barcode_table(barcodes), encoding="utf-8")
    check_file.write_text(_check_table(checks), encoding="utf-8")
    for path in (set_file, check_file):
        typer.echo(str(path))


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


def _write_gene(gene: CodingSequence, directory: Path) -> tuple[Path, ...]:
    """Write the coding sequence as a record, and every codon it moved as a table."""
    directory.mkdir(parents=True, exist_ok=True)
    record = directory / CODING_FILE
    write_dna(coding_record(gene.dna, gene.name), record)
    rows = [
        (
            str(change.codon_index),
            change.old_codon,
            change.new_codon,
            change.amino_acid,
            change.site.enzyme.name,
        )
        for change in gene.changes
    ]
    changes = directory / CHANGE_FILE
    changes.write_text(_tsv(CHANGE_COLUMNS, rows), encoding="utf-8")
    return record, changes


def _barcode_checks(barcodes: tuple[str, ...], rules: BarcodeRules) -> tuple[Check, ...]:
    """Judge the set: the rules each barcode holds to, how far apart it stands, and its cost.

    Deletion ambiguity carries no verdict: no sourced threshold judges it, and it is the
    number to read beside a set designed on mismatches.
    """
    broken = check_barcodes(barcodes, rules)
    checks = [
        Check(
            "rules",
            "pass" if not broken else "fail",
            len(broken),
            "; ".join(broken) or "every barcode holds",
        )
    ]
    if len(barcodes) > 1:
        apart = separation(barcodes, metric=rules.metric)
        checks.append(
            Check(
                "separation",
                "pass" if apart >= rules.distance else "fail",
                apart,
                f"{rules.metric}, against a rule of {rules.distance}",
            )
        )
        checks.append(
            Check(
                "deletion_ambiguity",
                None,
                deletion_ambiguity(barcodes),
                "share of one-base deletions another barcode could leave too",
            )
        )
    return tuple(checks)


def _barcode_table(barcodes: tuple[str, ...]) -> str:
    """Return the set as a table, each barcode numbered in the order it was accepted."""
    return _tsv(
        BARCODE_COLUMNS, [(str(at), barcode) for at, barcode in enumerate(barcodes, start=1)]
    )


def _check_table(checks: tuple[Check, ...]) -> str:
    """Return the checks as a table, a check nothing judged carrying no verdict."""
    return _tsv(
        CHECK_COLUMNS,
        [(check.name, check.status or "", f"{check.value:g}", check.detail) for check in checks],
    )


def _tsv(columns: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
    """Return a tab-separated table, the header alone where there are no rows."""
    return "\n".join(["\t".join(columns), *("\t".join(row) for row in rows)]) + "\n"


def _band(text: str) -> tuple[float, float] | None:
    """Read `--gc-band` as the share low and high, or nothing where none was given.

    Raises
    ------
    ValueError
        If it is not two numbers separated by a dash.
    """
    if not text:
        return GC_BAND
    match = _BAND.fullmatch(text)
    if match is None:
        raise ValueError(f"--gc-band is LOW-HIGH, each 0 to 1, got {text!r}")
    return float(match[1]), float(match[2])


def _metric(text: str) -> Metric:
    """Read the distance metric, refusing anything else.

    Raises
    ------
    ValueError
        If it is neither ``hamming`` nor ``sequence-levenshtein``.
    """
    if text == "hamming":
        return "hamming"
    if text == "sequence-levenshtein":
        return "sequence-levenshtein"
    raise ValueError(f"--metric is 'hamming' or 'sequence-levenshtein', got {text!r}")
