"""Check primer pairs against a genome FASTA, and design a pair specific on one.

A check reports every amplicon a pair makes there and which is intended; a design searches a
region's flanks for the best-ranked pair that makes only the intended one. The in-silico PCR
tool `ipcr` does the search; ``docs/adr/0003-genome-check.md`` says why.
"""

import json
import shutil
import subprocess
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Any, Literal

from liulab_mbio.checks import Check, Status, worst
from liulab_mbio.io import read_region
from liulab_mbio.primers.design import ranked_pairs
from liulab_mbio.primers.evaluation import PairReport
from liulab_mbio.primers.placement import Placement
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import TARGET_TM, THRESHOLDS, Thresholds
from liulab_mbio.sequence import BindingSite, Primer, Segment, Strand, reverse_complement

#: Which primers make an amplicon: the two together, or one alone.
type MadeBy = Literal["pair", "forward", "reverse"]

#: How many off-target amplicons a check's detail lists by place before only counting the rest.
_LISTED = 10


@dataclass(frozen=True, slots=True)
class Locus:
    """A span of one sequence in a genome FASTA.

    Parameters
    ----------
    sequence_name
        As the FASTA spells it: its header up to the first space.
    start, end
        0-based, half-open. `end` passes the sequence's length across the origin of a circular
        one.
    """

    sequence_name: str
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class Amplicon:
    """One amplicon on a genome.

    Parameters
    ----------
    sequence_name, start, end
        Where it lies, from one binding site's 5' end to the other's, as in `Locus`.
    length
        Tails included, as in `PairReport.amplicon_length`.
    made_by
        ``"pair"``, or the one primer that makes it alone.
    mismatches
        At each binding site: the forward primer's then the reverse primer's for a pair, and
        the site at `start` then the one at `end` for a primer alone.
    intended
        Whether it is the amplicon named as intended.
    """

    sequence_name: str
    start: int
    end: int
    length: int
    made_by: MadeBy
    mismatches: tuple[int, int]
    intended: bool


@dataclass(frozen=True, slots=True)
class GenomeReport:
    """What one primer pair makes on a genome, and the checks on it.

    Parameters
    ----------
    assembly
        The genome's name, such as ``"hg38"``.
    forward, reverse
        The pair checked.
    amplicons
        Every amplicon, by sequence name and then position.
    mismatches
        The most a binding site could carry: zero is a perfect-match screen.
    terminal_window
        The 3'-terminal bases held to a perfect match.
    checks
        Off-target amplicons, and the intended amplicon where one was named.
    """

    assembly: str
    forward: Primer
    reverse: Primer
    amplicons: tuple[Amplicon, ...]
    mismatches: int
    terminal_window: int
    checks: tuple[Check, ...]

    @property
    def near_matches_checked(self) -> bool:
        """Return whether the search found binding sites that do not match perfectly."""
        return self.mismatches > 0

    @property
    def off_target(self) -> tuple[Amplicon, ...]:
        """Return every amplicon but the intended one."""
        return tuple(amplicon for amplicon in self.amplicons if not amplicon.intended)

    @property
    def status(self) -> Status:
        """Return the worst status of any check."""
        return worst(check.status for check in self.checks)

    def __getitem__(self, name: str) -> Check:
        """Return the check of that name.

        Raises
        ------
        KeyError
            If no check has it.
        """
        for check in self.checks:
            if check.name == name:
                return check
        raise KeyError(name)


def evaluate_pair_on_genome(
    forward: Primer,
    reverse: Primer,
    fasta: str | Path,
    assembly: str,
    *,
    intended: Locus | None = None,
    circular: bool = False,
    thresholds: Thresholds = THRESHOLDS,
    timeout: float = 300.0,
) -> GenomeReport:
    """Check one primer pair against a genome FASTA, as `evaluate_on_genome` checks several.

    Examples
    --------
    >>> forward = Primer("fwd", "ACCTAGGAGAAGTGGCCAGC")
    >>> reverse = Primer("rev", "GCCTTGGAGAGGAGTTCTGC")
    >>> target = Locus("chr1", 156710301, 156710801)
    >>> report = evaluate_pair_on_genome(forward, reverse, "hg38.fa", "hg38", intended=target)  # doctest: +SKIP
    >>> report.status, report.near_matches_checked  # doctest: +SKIP
    ('pass', False)
    """
    return evaluate_on_genome(
        [(forward, reverse)],
        fasta,
        assembly,
        intended=[intended],
        circular=circular,
        thresholds=thresholds,
        timeout=timeout,
    )[0]


def evaluate_on_genome(
    pairs: Sequence[tuple[Primer, Primer]],
    fasta: str | Path,
    assembly: str,
    *,
    intended: Sequence[Locus | None] = (),
    circular: bool = False,
    thresholds: Thresholds = THRESHOLDS,
    timeout: float = 300.0,
) -> tuple[GenomeReport, ...]:
    """Check primer pairs against a genome FASTA in one search, with a report for each pair.

    Every pair is its own reaction: primers of different pairs are never combined. Each primer
    is searched by its annealing region, so a tail only lengthens an amplicon. On a FASTA up to
    `Thresholds.genome_size_limit` bytes on disk a binding site may carry
    `Thresholds.genome_mismatches`, none in its `Thresholds.genome_terminal_window`; on a larger
    one it must match perfectly, and each report says near-match sites were not checked. A site
    must also melt within `Thresholds.off_target_margin` of a perfect match, or its amplicon is
    dropped. Amplicons run up to `Thresholds.genome_max_amplicon` bases, or to the intended
    one's length where that is longer.

    Every amplicon other than the intended one is an off-target amplicon, so where no amplicon
    is named as intended, any amplicon warns. `ipcr` keeps 10,000 sites per primer on one
    sequence, so a primer in a repeat that common lists only some of its amplicons.

    Parameters
    ----------
    intended
        One entry for each pair, ``None`` for a pair with none; empty names none at all.
    circular
        Whether every record in the FASTA is circular, as a plasmid is.
    timeout
        Seconds the search may run before it is killed.

    Raises
    ------
    RuntimeError
        If `ipcr` is not installed, or refuses the search.
    TimeoutError
        If the search runs past `timeout`.
    ValueError
        If `intended` is neither empty nor one entry for each pair.
    """
    loci = tuple(intended) or (None,) * len(pairs)
    if len(loci) != len(pairs):
        raise ValueError(f"{len(loci)} intended amplicons for {len(pairs)} pairs")
    fasta = Path(fasta)
    below = fasta.stat().st_size <= thresholds.genome_size_limit
    mismatches = thresholds.genome_mismatches if below else thresholds.genome_mismatches_above_limit
    longest = [
        max(thresholds.genome_max_amplicon, 0 if locus is None else locus.end - locus.start)
        for locus in loci
    ]
    annealing = [(_annealing(forward), _annealing(reverse)) for forward, reverse in pairs]
    hits = _search(
        annealing,
        fasta,
        mismatches=mismatches,
        window=thresholds.genome_terminal_window,
        max_length=max(longest, default=0),
        circular=circular,
        timeout=timeout,
    )
    found: list[list[Amplicon]] = [[] for _ in pairs]
    for hit in hits:
        index = hit.pair
        if hit.end - hit.start > longest[index]:
            continue
        forward, reverse = pairs[index]
        left_primer, right_primer = {
            "pair": (forward, reverse) if hit.forward_first else (reverse, forward),
            "forward": (forward, forward),
            "reverse": (reverse, reverse),
        }[hit.made_by]
        left, right = _annealing(left_primer), _annealing(right_primer)
        if hit.bases and not (
            _melts(left, reverse_complement(hit.bases[: len(left)]), thresholds)
            and _melts(right, hit.bases[-len(right) :], thresholds)
        ):
            continue
        counts = (hit.left_mismatches, hit.right_mismatches)
        if not hit.forward_first:
            counts = (hit.right_mismatches, hit.left_mismatches)
        locus = loci[index]
        found[index].append(
            Amplicon(
                hit.sequence_name,
                hit.start,
                hit.end,
                hit.end - hit.start + _tail(left_primer) + _tail(right_primer),
                hit.made_by,
                counts,
                hit.made_by == "pair" and locus == Locus(hit.sequence_name, hit.start, hit.end),
            )
        )
    return tuple(
        GenomeReport(
            assembly,
            forward,
            reverse,
            tuple(amplicons),
            mismatches,
            thresholds.genome_terminal_window,
            _checks(amplicons, locus, mismatches, thresholds),
        )
        for (forward, reverse), amplicons, locus in zip(pairs, found, loci, strict=True)
    )


@dataclass(frozen=True, slots=True)
class GenomeDesign:
    """A primer pair designed for a genome region, and what it makes on that genome.

    Parameters
    ----------
    template
        Where the region and its flanks lie, so a binding site's position counts from its start.
    pair
        Every check on the pair, on that template.
    genome
        Every amplicon it makes on the genome, and the checks on them.
    rounds
        Genome searches run.
    specific
        Whether it makes no off-target amplicon.
    exhausted
        Whether the pairs ran out, so nothing was left to check.
    """

    template: Locus
    pair: PairReport
    genome: GenomeReport
    rounds: int
    specific: bool
    exhausted: bool

    @property
    def forward(self) -> Primer:
        """Return the forward primer."""
        return self.pair.forward.primer

    @property
    def reverse(self) -> Primer:
        """Return the reverse primer."""
        return self.pair.reverse.primer


def design_pair_on_genome(
    fasta: str | Path,
    assembly: str,
    region: Locus,
    *,
    flank: int,
    forward_tail: str = "",
    reverse_tail: str = "",
    forward_name: str = "",
    reverse_name: str = "",
    target_tm: float = TARGET_TM,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
    timeout: float = 300.0,
) -> GenomeDesign:
    """Design the best-ranked pair amplifying a genome region that is specific on that genome.

    The template is the region and a flank each side, read through the FASTA's index and
    clipped to the sequence. The forward primer is placed in the left flank and the reverse in
    the right, so the amplicon covers the region whole and ties go to the pair nearest it, as
    `design_pair` ranks one. `Thresholds.genome_pairs_per_search` pairs are then checked in one
    genome search, best-ranked first, and the first making no off-target amplicon wins.

    Where none does, the round leaves the search narrower: a primer that makes an off-target
    amplicon on its own is kept out of every later pair, as is every binding site sharing its
    3' end, which primes the same place; a pair that makes one only together is passed over,
    each of its primers still free to pair with another. After `Thresholds.genome_rounds`
    searches, or once the pairs run out, the best-ranked pair checked comes back with its
    off-target amplicons, marked not specific.

    Raises
    ------
    FileNotFoundError
        If no index lies beside the FASTA.
    RuntimeError
        If `ipcr` is not installed, or refuses a search.
    TimeoutError
        If a search runs past `timeout`.
    ValueError
        If the sequence does not hold the region, or a flank holds no annealing region.

    Examples
    --------
    >>> region = Locus("chr1", 156710401, 156710701)
    >>> design = design_pair_on_genome("hg38.fa", "hg38", region, flank=100)  # doctest: +SKIP
    >>> design.specific, design.rounds, design.forward.sequence  # doctest: +SKIP
    (True, 1, 'ACCTAGGAGAAGTGGCCAGC')
    """
    start = max(region.start - flank, 0)
    template = read_region(fasta, region.sequence_name, start, region.end + flank)
    where = Locus(region.sequence_name, start, start + len(template))
    left, right = region.start - where.start, region.end - where.start
    if left < 1 or right >= len(template):
        raise ValueError(
            f"{region.sequence_name} holds no flank either side of {region.start}-{region.end}"
        )
    excluded = _ThreePrimeEnds()
    pairs = ranked_pairs(
        template,
        left,
        right,
        forward_placement=Placement(five_prime=Segment(0, left), three_prime=Segment(0, left + 1)),
        reverse_placement=Placement(
            five_prime=Segment(right + 1, len(template) + 1),
            three_prime=Segment(right, len(template)),
        ),
        forward_tail=forward_tail,
        reverse_tail=reverse_tail,
        forward_name=forward_name,
        reverse_name=reverse_name,
        target_tm=target_tm,
        polymerase=polymerase,
        thresholds=thresholds,
        exclude=excluded,
    )
    best: tuple[PairReport, GenomeReport] | None = None
    rounds = 0
    exhausted = False
    while rounds < thresholds.genome_rounds and not exhausted:
        batch = list(islice(pairs, thresholds.genome_pairs_per_search))
        exhausted = len(batch) < thresholds.genome_pairs_per_search
        if not batch:
            break
        rounds += 1
        reports = evaluate_on_genome(
            [(pair.forward.primer, pair.reverse.primer) for pair in batch],
            fasta,
            assembly,
            intended=[_intended(pair, where) for pair in batch],
            thresholds=thresholds,
            timeout=timeout,
        )
        for pair, report in zip(batch, reports, strict=True):
            best = best or (pair, report)
            if not report.off_target:
                return GenomeDesign(where, pair, report, rounds, True, exhausted)
            for amplicon in report.off_target:
                if amplicon.made_by != "pair":
                    alone = pair.forward if amplicon.made_by == "forward" else pair.reverse
                    excluded.add(alone.primer.binding_sites[0])
    if best is None:
        raise ValueError(f"no pair amplifies {region.sequence_name}:{region.start}-{region.end}")
    return GenomeDesign(where, *best, rounds, False, exhausted)


def _intended(pair: PairReport, template: Locus) -> Locus:
    """Return where a designed pair's amplicon lies on the genome.

    The one place a template coordinate becomes a genome one.
    """
    forward = pair.forward.primer.binding_sites[0]
    reverse = pair.reverse.primer.binding_sites[0]
    return Locus(
        template.sequence_name, template.start + forward.start, template.start + reverse.end
    )


class _ThreePrimeEnds:
    """The 3' ends a design has ruled out, matching a binding site of any length.

    A primer primes an off-target site by its 3' end, so a site sharing that end primes there
    too: a shorter one matches whatever the longer one did, and a longer one only adds bases
    away from the end. Ruling out the end rather than the site keeps every one of them out.
    """

    def __init__(self) -> None:
        self._ends: set[tuple[Strand, int]] = set()

    def add(self, site: BindingSite) -> None:
        """Rule out every binding site with this one's 3' end."""
        self._ends.add(_three_prime(site))

    def __contains__(self, site: object) -> bool:
        """Return whether a binding site's 3' end is ruled out."""
        return isinstance(site, BindingSite) and _three_prime(site) in self._ends


def _three_prime(site: BindingSite) -> tuple[Strand, int]:
    """Return the strand a binding site lies on and where its 3' end is."""
    return site.strand, site.end if site.strand is Strand.FORWARD else site.start


@dataclass(frozen=True, slots=True)
class _Hit:
    """One amplicon as the search reports it, before the Tm rule and the intended one are judged.

    `left_mismatches` is the binding site's at `start`, `right_mismatches` the one's at `end`.
    `forward_first` says the forward primer's site is the one at `start`, or would be for a
    primer alone. `bases` are the top strand from `start` to `end`, or empty where the search
    allowed no mismatch and so read none.
    """

    pair: int
    made_by: MadeBy
    sequence_name: str
    start: int
    end: int
    left_mismatches: int
    right_mismatches: int
    forward_first: bool
    bases: str


def _search(
    pairs: Sequence[tuple[str, str]],
    fasta: Path,
    *,
    mismatches: int,
    window: int,
    max_length: int,
    circular: bool,
    timeout: float,
) -> list[_Hit]:
    """Run `ipcr` once over every pair and return each amplicon once, by name and position.

    The one place the tool runs. `ipcr` misses a binding site across the origin of a circular
    record, so each circular record is searched written twice over, and only an amplicon
    starting in its first copy is kept, its end passing the record's length as ADR 0001 has it.
    """
    tool = shutil.which("ipcr")
    if tool is None:
        raise RuntimeError(
            "ipcr is not installed: it is a pixi dependency of liulab-mbio, from bioconda, "
            "so run `pixi install`"
        )
    if not pairs:
        return []
    with tempfile.TemporaryDirectory() as scratch:
        table = Path(scratch) / "pairs.tsv"
        table.write_text("".join(f"{i}\t{f}\t{r}\n" for i, (f, r) in enumerate(pairs)))
        lengths: dict[str, int] = {}
        source = fasta
        if circular:
            source = Path(scratch) / "circular.fa"
            lengths = _write_twice(fasta, source)
        command = [
            tool,
            *("--primers", str(table), "--sequences", str(source), "--output", "jsonl"),
            *("--mismatches", str(mismatches), "--terminal-window", str(window)),
            *("--max-length", str(max_length), "--quiet"),
            *(("--products",) if mismatches else ()),
        ]
        try:
            run = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            # `ipcr` ignores SIGTERM; `subprocess.run` sends SIGKILL, which it cannot.
            raise TimeoutError(f"ipcr ran past {timeout} s and was killed") from None
    if run.returncode != 0:
        raise RuntimeError(f"ipcr refused the search: {run.stderr.strip()}")
    hits: dict[tuple[int, str, str, int, int], _Hit] = {}
    for line in run.stdout.splitlines():
        if line:
            hit = _hit(json.loads(line))
            size = lengths.get(hit.sequence_name)
            if size is None or (hit.start < size and hit.end - hit.start <= size):
                hits.setdefault((hit.pair, hit.made_by, hit.sequence_name, hit.start, hit.end), hit)
    return sorted(hits.values(), key=lambda hit: (hit.sequence_name, hit.start, hit.end))


def _hit(row: dict[str, Any]) -> _Hit:
    """Read one line of `ipcr`'s JSON output.

    An amplicon one primer makes alone carries the pair's id and ``+A:self`` for the forward
    primer or ``+B:self`` for the reverse, and comes once for each strand. Mismatch counts
    belong to the sites at the start and the end, and are absent where zero.
    """
    pair, _, alone = str(row["experiment_id"]).partition("+")
    made_by: MadeBy = "pair" if not alone else "forward" if alone.startswith("A") else "reverse"
    return _Hit(
        int(pair),
        made_by,
        str(row["sequence_id"]),
        int(row["start"]),
        int(row["end"]),
        int(row.get("fwd_mm", 0)),
        int(row.get("rev_mm", 0)),
        made_by != "pair" or row["type"] == "forward",
        str(row.get("seq", "")).upper(),
    )


def _write_twice(fasta: Path, target: Path) -> dict[str, int]:
    """Write each record of `fasta` to `target` twice over, and return each one's length."""
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    lengths = {}
    with fasta.open() as source, target.open("w") as out:
        for title, bases in SimpleFastaParser(source):
            lengths[title.split(maxsplit=1)[0] if title.strip() else ""] = len(bases)
            out.write(f">{title}\n{bases}{bases}\n")
    return lengths


def _melts(primer: str, target: str, thresholds: Thresholds) -> bool:
    """Return whether `primer` melts off `target` within the off-target margin of a match.

    `target` is the strand it pairs with, 5' to 3'. The template check's rule, from
    `find_priming_sites`.
    """
    import primer3

    perfect = primer3.calc_end_stability(primer, reverse_complement(primer)).tm
    return primer3.calc_end_stability(primer, target).tm >= perfect - thresholds.off_target_margin


def _annealing(primer: Primer) -> str:
    """Return the part of a primer its binding sites cover, or all of it where it carries none."""
    size = max((site.end - site.start for site in primer.binding_sites), default=0)
    return primer.sequence[-size:] if size else primer.sequence


def _tail(primer: Primer) -> int:
    return len(primer.sequence) - len(_annealing(primer))


def _checks(
    amplicons: list[Amplicon], locus: Locus | None, mismatches: int, thresholds: Thresholds
) -> tuple[Check, ...]:
    off_target = [amplicon for amplicon in amplicons if not amplicon.intended]
    notes = [_described(amplicon) for amplicon in off_target[:_LISTED]]
    if len(off_target) > _LISTED:
        notes.append(f"and {len(off_target) - _LISTED} more")
    if not mismatches:
        notes.append("perfect matches only: near-match sites were not checked")
    checks = [
        Check(
            "off_target_amplicons",
            thresholds.off_target_amplicons.grade(len(off_target)),
            len(off_target),
            "; ".join(notes),
        )
    ]
    if locus is not None:
        present = sum(amplicon.intended for amplicon in amplicons)
        checks.append(
            Check(
                "intended_amplicon",
                thresholds.intended_amplicon.grade(present),
                present,
                "" if present else f"none at {locus.sequence_name}:{locus.start}-{locus.end}",
            )
        )
    return tuple(checks)


def _described(amplicon: Amplicon) -> str:
    who = "the pair" if amplicon.made_by == "pair" else f"{amplicon.made_by} primer alone"
    where = f"{amplicon.sequence_name}:{amplicon.start}-{amplicon.end}"
    return f"{where}, {amplicon.length} bp, by {who}"
