"""Choosing the enzyme pair this method cuts with, and saying why every other pair was refused.

A pair has to do four things, and each is a rule the plan is already held to: read a site in the
vector, get the insert to two cut ends, leave a backbone that cannot close on itself, and leave
ends that anneal. So a candidate is weighed by running the plan's own `opened` and `excised` over
it and keeping what they refuse, and a pair the chooser refuses reads the same as the pair a user
names. `Refusal` is one pair and the rule that refused it, as `liulab_mbio.overhangs.refusal`
already reads for one overhang.

Where a site in the insert is what refused the pair and a synonymous codon change could take it
away, the refusal carries the change: `liulab_mbio.sites.domesticate` reports it, and clearing
the insert of both enzymes' sites is what leaves it to be amplified with the sites on its tails.
A site lying in no coding sequence is a flat refusal, because taking it out would change what the
record spells.

The chooser stays inside this method. Which pair opens a vector where its insert goes is this
method's judgement and nothing below it asks the question -- `docs/adr/0007-cloning-methods.md`
says why.
"""

import dataclasses
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from itertools import combinations
from typing import Literal

from liulab_mbio.bench.steps import listed
from liulab_mbio.cloning.restriction.digest import Piece, cut, excised, opened, resolve
from liulab_mbio.cloning.restriction.verdicts import (
    buffer_check,
    cleanup_check,
    methylation_check,
    temperature_check,
)
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.enzymes import enzymes as shipped
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_mbio.sites import CutSite, Domestication, DomesticationReport, EnzymeLike, domesticate
from liulab_mbio.sites import find_sites as sites_in

#: Why a candidate pair will not do. The first two are one enzyme's doing and the other two the
#: pair's, which is the order they are weighed in and the reverse of how far a pair got.
type PairRule = Literal["vector site", "insert site", "backbone", "ends"]

#: The order refusals are read in: the pairs that got furthest first.
READING_ORDER: tuple[PairRule, ...] = ("ends", "backbone", "insert site", "vector site")


@dataclass(frozen=True, slots=True)
class Refusal:
    """One candidate pair a rule refused, and what keeping its enzymes would cost.

    Parameters
    ----------
    enzymes
        The pair that was weighed.
    rule
        Which rule refused it.
    detail
        Why that rule refused this pair, in the words the plan itself refuses with.
    changes
        The synonymous codon changes that would clear the insert of both enzymes' sites and leave
        it to be amplified. Empty where no site is in the way, and empty where one lies somewhere
        no change reaches, which is the flat refusal.
    """

    enzymes: tuple[Enzyme, ...]
    rule: PairRule
    detail: str
    changes: tuple[Domestication, ...] = ()

    @property
    def names(self) -> str:
        """The pair, as a sentence names it."""
        return listed([one.name for one in self.enzymes])

    @property
    def domesticable(self) -> bool:
        """Whether a synonymous codon change could take away what refused this pair."""
        return bool(self.changes)


@dataclass(frozen=True, slots=True)
class Chosen:
    """The pair a plan cuts with, and every pair it refused on the way.

    Parameters
    ----------
    enzymes
        The pair to cut both records with.
    refusals
        Every other pair weighed, with the rule that refused it, in `READING_ORDER`.
    """

    enzymes: tuple[Enzyme, ...]
    refusals: tuple[Refusal, ...] = ()


def candidates() -> tuple[Enzyme, ...]:
    """Every shipped enzyme this method can choose from, ordered by name.

    Type II only. This method's junction is the recognition site the two ends came from, and an
    enzyme cutting outside its own site leaves an overhang the record happens to spell there
    rather than one the junction puts back.
    """
    return tuple(one for one in shipped() if one.type == "II")


def choose_pair(
    vector: SequenceRecord,
    source: SequenceRecord,
    *,
    enzymes: Iterable[EnzymeLike] | None = None,
) -> Chosen:
    """Choose the pair to cut `vector` and `source` with, and report every pair refused.

    Every pair of `enzymes` is weighed by the rules the plan itself is held to. What passes is
    ranked on what a digest of the pair would carry: the verdicts
    `liulab_mbio.cloning.restriction.verdicts` gives it, worst counted first and a verdict nothing
    sourced can give counted as its own thing rather than as a pass; then the source's features
    the insert carries whole, which is what the cloning moves; then the vector bases the backbone
    gives up; then the longer insert. A source annotating nothing leaves the third rule deciding,
    so the plan keeps as much of the vector as it can.

    Parameters
    ----------
    vector
        The circular plasmid the insert goes into.
    source
        The insert, or the plasmid it is cut out of.
    enzymes
        The candidates, each an `liulab_mbio.enzymes.Enzyme` or a name the package ships.
        `candidates` by default.

    Raises
    ------
    KeyError
        If a name is not one the package ships.
    ValueError
        If no pair can be found, naming how many were refused, by which rules, and the detail of
        the pair that got furthest.
    """
    pool = candidates() if enzymes is None else resolve(enzymes)
    weighing = _Weighing(vector, source)
    refused: list[Refusal] = []
    ranked: list[tuple[tuple[int, ...], str, tuple[Enzyme, ...]]] = []
    for pair in combinations(pool, 2):
        weighed = weighing.weigh(pair)
        if isinstance(weighed, Refusal):
            refused.append(weighed)
        else:
            ranked.append((*_rank(vector, source, pair, *weighed), pair))
    reported = tuple(sorted(refused, key=_read_first))
    if not ranked:
        raise ValueError(_stuck(pool, reported))
    return Chosen(min(ranked, key=lambda entry: entry[:2])[2], reported)


def refusal(
    vector: SequenceRecord, source: SequenceRecord, enzymes: Sequence[Enzyme]
) -> Refusal | None:
    """Why this pair will not cut `vector` and `source`, or ``None`` when it will.

    The rules are the plan's own, run in the order it runs them, so a pair a user names is
    refused in the same words the chooser refuses it in.

    Examples
    --------
    >>> from liulab_mbio.enzymes import get_enzyme
    >>> plasmid = SequenceRecord("GAATTC" + "ACGT" * 8, topology="circular", name="p")
    >>> refusal(plasmid, plasmid, [get_enzyme("BamHI")]).rule
    'vector site'
    """
    weighed = _Weighing(vector, source).weigh(tuple(enzymes))
    return weighed if isinstance(weighed, Refusal) else None


@dataclass(frozen=True, slots=True)
class _Weighing:
    """One vector and one insert, and what has already been read of each candidate enzyme.

    An enzyme is read once however many pairs it appears in: the sites it reads in the vector,
    and what a synonymous codon change could do about the sites it reads in the insert.
    """

    vector: SequenceRecord
    source: SequenceRecord
    blocked: dict[str, Refusal | None] = field(default_factory=dict, hash=False)
    reports: dict[str, DomesticationReport] = field(default_factory=dict, hash=False)

    def weigh(self, enzymes: tuple[Enzyme, ...]) -> tuple[Piece, Piece | None] | Refusal:
        """Return the backbone and the insert this pair leaves, or why it will not do.

        The insert is ``None`` where the record reads no site of either enzyme, because it is
        amplified with the sites on its primer tails and no digest of it is designed here.
        """
        for one in enzymes:
            if (found := self._vector_refusal(one)) is not None:
                return dataclasses.replace(found, enzymes=enzymes)
        if (found := self._insert_refusal(enzymes)) is not None:
            return found
        try:
            backbone = opened(self.vector, enzymes)[0]
        except ValueError as error:
            return Refusal(enzymes, "backbone", str(error))
        if not _reads_a_site(self.source, enzymes):
            return backbone, None
        try:
            return backbone, excised(self.source, enzymes, into=backbone)[0]
        except ValueError as error:
            return Refusal(enzymes, "ends", str(error))

    def _vector_refusal(self, enzyme: Enzyme) -> Refusal | None:
        """Refuse an enzyme that reads no site in the vector, in `cut`'s own words."""
        if enzyme.name not in self.blocked:
            self.blocked[enzyme.name] = None
            try:
                cut(self.vector, (enzyme,), unique=False)
            except ValueError as error:
                self.blocked[enzyme.name] = Refusal((enzyme,), "vector site", str(error))
        return self.blocked[enzyme.name]

    def _insert_refusal(self, enzymes: tuple[Enzyme, ...]) -> Refusal | None:
        """Refuse an insert record neither route can get to two cut ends.

        One site of each enzyme is cut out and goes in as it is; none of either is amplified with
        the sites on its primer tails. Anything between is neither, and the refusal says what a
        synonymous codon change could do about the sites in the way.
        """
        found = {one.name: _cut_sites(self.source, one) for one in enzymes}
        counts = {len(sites) for sites in found.values()}
        if counts <= {0} or counts == {1}:
            return None
        said = listed([_reads(name, sites) for name, sites in found.items()])
        what = self.source.name or "the insert"
        return _with_domestication(
            Refusal(
                enzymes,
                "insert site",
                f"{what} reads {said}; this method cuts the insert out between one site of each "
                "enzyme, and a record reading none of them is amplified with the sites on its "
                "primer tails instead, so neither route reaches two cut ends",
            ),
            self._report(enzymes),
        )

    def _report(self, enzymes: tuple[Enzyme, ...]) -> DomesticationReport:
        """Report what a synonymous codon change could do about their sites in the insert.

        Each enzyme is weighed on its own, so its report is read once however many pairs it is
        in, and the pair's is the two put together.
        """
        for one in enzymes:
            if one.name not in self.reports:
                _, self.reports[one.name] = domesticate(self.source, one)
        return DomesticationReport(
            tuple(change for one in enzymes for change in self.reports[one.name].changes),
            tuple(site for one in enzymes for site in self.reports[one.name].outside_cds),
            tuple(site for one in enzymes for site in self.reports[one.name].unchanged),
        )


def _with_domestication(refused: Refusal, report: DomesticationReport) -> Refusal:
    """Add to a refusal what keeping its enzymes would cost the insert's coding sequences."""
    stuck = (*report.outside_cds, *report.unchanged)
    if stuck:
        where = listed([f"{site.enzyme.name} at {site.start}" for site in stuck])
        return dataclasses.replace(
            refused,
            detail=f"{refused.detail}; no synonymous codon change reaches {where}, and taking "
            "those bases out would change what the record spells",
        )
    changed = listed(
        [
            f"{one.site.enzyme.name} at {one.site.start} by {one.old_codon} to {one.new_codon} "
            f"in {one.feature.name}"
            for one in report.changes
        ]
    )
    return dataclasses.replace(
        refused,
        detail=f"{refused.detail}; a synonymous codon change reaches every one of them, leaving "
        f"the insert to be amplified instead: {changed}",
        changes=report.changes,
    )


def _rank(
    vector: SequenceRecord,
    source: SequenceRecord,
    enzymes: tuple[Enzyme, ...],
    backbone: Piece,
    insert: Piece | None,
) -> tuple[tuple[int, ...], str]:
    """Rank one pair that passed, least first. `choose_pair` documents what each place means."""
    verdicts = [
        check.status
        for check in (
            buffer_check(enzymes),
            temperature_check(enzymes),
            cleanup_check(enzymes),
            methylation_check(enzymes, (vector, source)),
        )
    ]
    return (
        verdicts.count("fail"),
        verdicts.count("warn"),
        verdicts.count(None),
        -_carries(insert),
        len(vector) - backbone.length,
        -(0 if insert is None else insert.length),
    ), " ".join(one.name for one in enzymes)


def _carries(insert: Piece | None) -> int:
    """How many of its own record's features the insert carries whole, which is what it moves."""
    if insert is None:
        return 0
    length = len(insert.source)
    return sum(
        1
        for feature in insert.source.features
        if all(_inside(segment, insert, length) for segment in feature.segments)
    )


def _inside(segment: Segment, piece: Piece, length: int) -> bool:
    """Whether a span of the piece's own record lies inside the piece, across the origin or not."""
    return (segment.start - piece.start) % length + segment.end - segment.start <= piece.length


def _reads_a_site(record: SequenceRecord, enzymes: Sequence[Enzyme]) -> bool:
    """Whether any of these enzymes cuts `record`, which is what decides the insert's route."""
    return any(site.cuts for site in sites_in(record, enzymes))


def _cut_sites(record: SequenceRecord, enzyme: Enzyme) -> tuple[CutSite, ...]:
    """Every site of one enzyme the record is really cut at."""
    return tuple(site for site in sites_in(record, enzyme) if site.cuts)


def _reads(name: str, sites: Sequence[CutSite]) -> str:
    """Say how many sites of one enzyme a record reads, and where."""
    if not sites:
        return f"no {name} site"
    where = listed([str(site.start) for site in sites])
    return f"{len(sites)} {name} site{'' if len(sites) == 1 else 's'} (at {where})"


def _read_first(refused: Refusal) -> tuple[int, str]:
    """Order refusals by how far the pair got, then by name, so a report reads the same twice."""
    return READING_ORDER.index(refused.rule), " ".join(one.name for one in refused.enzymes)


def _stuck(pool: Sequence[Enzyme], refused: Sequence[Refusal]) -> ValueError:
    """Return the error for a vector and an insert no pair of these enzymes can join."""
    if not refused:
        return ValueError(
            f"a pair takes two enzymes and {len(pool)} were offered: "
            f"{listed([one.name for one in pool]) or 'none'}"
        )
    counted = Counter(one.rule for one in refused)
    rules = listed([f"{count} on the {rule} rule" for rule, count in counted.most_common()])
    nearest = refused[0]
    return ValueError(
        f"no pair of these {len(pool)} enzymes cuts both records; every pair was refused, "
        f"{rules}. The pair that got furthest was {nearest.names}, and {nearest.detail}"
    )
