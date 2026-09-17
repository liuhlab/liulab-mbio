"""What a map draws from a record: its items, with their names, colours and hover details.

A view takes these items and never the record. A map draws a record's features, each primer at
its binding sites, and the cut sites of the shipped unique cutters or of the enzymes named. A
primer carries its 5' tail and the bases it does not pair with, found by comparing it with the
record, and each enzyme at a cut site how its cuts in the two strands stagger.

A feature draws in its file's colour, and a segment in its own where that differs. One the file
gives no colour takes Paul Tol's light scheme by the group its type falls in, and pale grey for a
type in none. Primers are purple and enzyme names black.

Positions a person reads, in labels and hover details, are 1-based and inclusive. A cut site is
numbered as SnapGene numbers one: by the base after which its enzymes cut the top strand.

Asked for, a CDS carries its translation, read with the standard genetic code across its joined
segments from its `/codon_start`, the whole feature through, stops included.

Where labels crowd past what a map grows to, they hide in the order `hiding` sorts them, and
`notice` says how many hid.
"""

import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import cache
from typing import Literal

from liulab_mbio.codons import amino_acid
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.enzymes import enzymes as shipped
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)
from liulab_mbio.sites import find_sites

#: What an item is: a feature, a primer at one binding site, or the enzymes cutting at one position.
type Kind = Literal["feature", "primer", "cut_site"]

_HEX = re.compile(r"#[0-9a-fA-F]{6}")

# The colours below are from Paul Tol's `tol_colors.py`, shipped under its licence:
#
# Copyright (c) 2022, Paul Tol
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
#    conditions and the following disclaimer in the documentation and/or other materials
#    provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
#    endorse or promote products derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
# AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

#: Each group of feature types, its colour in Paul Tol's light scheme, and the keys it holds:
#: INSDC's current and retired keys, and SnapGene's `LTR`. A type in no group takes `OTHER`.
GROUPS: Mapping[str, tuple[str, frozenset[str]]] = {
    "coding": (
        "#77AADD",
        frozenset({
            "CDS", "gene", "mat_peptide", "sig_peptide", "transit_peptide", "propeptide",
            "proprotein", "exon", "C_region", "D_segment", "J_segment", "N_region", "S_region",
            "V_region", "V_segment",
        }),
    ),
    "promoters and enhancers": (
        "#EE8866",
        frozenset({
            "promoter", "enhancer", "RBS", "-10_signal", "-35_signal", "CAAT_signal",
            "TATA_signal", "GC_signal", "attenuator", "misc_signal", "5'UTR", "regulatory",
        }),
    ),
    "terminators and polyA": (
        "#FFAABB",
        frozenset({"terminator", "polyA_signal", "polyA_site", "3'UTR"}),
    ),
    "replication origins": ("#EEDD88", frozenset({"rep_origin", "oriT"})),
    "binding sites": ("#99DDFF", frozenset({"protein_bind", "misc_binding", "primer_bind"})),
    "RNA": (
        "#44BB99",
        frozenset({
            "ncRNA", "misc_RNA", "mRNA", "tRNA", "rRNA", "tmRNA", "precursor_RNA",
            "prim_transcript",
        }),
    ),
    "repeats and recombination": (
        "#BBCC33",
        frozenset({
            "LTR", "repeat_region", "mobile_element", "misc_recomb", "stem_loop", "D-loop",
        }),
    ),
}  # fmt: skip

#: Pale grey, for a type in no group.
OTHER = "#DDDDDD"

#: Tol bright purple, for primers and for a `primer_bind` feature the file gives no colour.
PRIMER = "#AA3377"

#: The colour of an enzyme's name.
ENZYME = "#000000"

#: What joins the names of the enzymes cutting at one position.
SEPARATOR = " - "

#: The `/regulatory_class` values that put a `regulatory` feature among the terminators.
TERMINATING_CLASSES = frozenset({"terminator", "polyA_signal_sequence"})

# The kinds in the order their labels hide.
_HIDING: tuple[Kind, ...] = ("cut_site", "primer", "feature")


@dataclass(frozen=True, slots=True)
class Span:
    """One stretch of an item, 0-based and half-open as `Segment` is, and the colour it fills."""

    start: int
    end: int
    color: str


@dataclass(frozen=True, slots=True)
class Cutter:
    """An enzyme a cut site names.

    Parameters
    ----------
    name
        The enzyme's.
    cuts
        How many times it cuts the whole record.
    stagger
        How many bases past its cut in the top strand it cuts the bottom strand: more than 0 where
        it leaves a 5' overhang, less than 0 where it leaves a 3' one, and 0 where it cuts blunt.
    """

    name: str
    cuts: int
    stagger: int = 0

    @property
    def unique(self) -> bool:
        """Whether it cuts the record once, which is what makes its name bold."""
        return self.cuts == 1


@dataclass(frozen=True, slots=True)
class Codon:
    """One codon of a CDS's translation: where its middle base lies, and what it spells.

    Parameters
    ----------
    middle
        The 0-based position of its middle base, less than the record's length.
    amino_acid
        One letter, ``"*"`` for a stop and ``"X"`` for a codon holding a base other than ACGT.
    """

    middle: int
    amino_acid: str

    @property
    def name(self) -> str:
        """The amino acid's three-letter code, ``Ter`` for a stop.

        Examples
        --------
        >>> Codon(1, "M").name, Codon(4, "*").name
        ('Met', 'Ter')
        """
        return _three_letters()[self.amino_acid]

    @property
    def stop(self) -> bool:
        """Whether it is a stop codon."""
        return self.amino_acid == "*"


@dataclass(frozen=True, slots=True)
class Item:
    """One thing a map draws.

    Parameters
    ----------
    kind
        What it is.
    name
        As the file writes it; for a cut site, its enzymes' names joined as its label joins them.
    type
        The feature type; ``primer`` or ``cut site`` for the other kinds.
    strand
        The strand it reads along.
    spans
        In top-strand order, as the feature's segments are. A primer's is its binding site. A cut
        site's is empty, at the boundary where its enzymes cut the top strand.
    label
        The text a map labels it with.
    hover
        What hovering over it shows: its name, type, span and length, as a person reads them.
    cutters
        A cut site's enzymes, in the order its label names them.
    translation
        A CDS's codons in the order they are read, when translations were asked for.
    tail
        How many of a primer's bases lie 5' of its binding site, pairing with nothing there.
    mismatches
        Where a primer's base does not pair with the record's in its binding site, found by
        comparing the two: 0-based positions less than the record's length, in order along the
        site.
    """

    kind: Kind
    name: str
    type: str
    strand: Strand
    spans: tuple[Span, ...]
    label: str
    hover: Mapping[str, str] = field(default_factory=dict, hash=False)
    cutters: tuple[Cutter, ...] = ()
    translation: tuple[Codon, ...] = field(default=(), hash=False)
    tail: int = 0
    mismatches: tuple[int, ...] = ()

    @property
    def color(self) -> str:
        """The colour of its label: its first span's."""
        return self.spans[0].color

    @property
    def runs(self) -> tuple[tuple[str, bool], ...]:
        """Its label in stretches of text, each marked bold or not: a unique cutter's name is."""
        if not self.cutters:
            return ((self.label, False),)
        runs: list[tuple[str, bool]] = []
        for cutter in self.cutters:
            if runs:
                runs.append((SEPARATOR, False))
            runs.append((cutter.name, cutter.unique))
        runs.append((self.label.removeprefix(self.name), False))
        return tuple(runs)


@dataclass(frozen=True, slots=True)
class Piece:
    """What of an item lies in a stretch of a record: part of one of its spans, or of a gap.

    Parameters
    ----------
    start, end
        Counted as the stretch counts them, running past the record's length where the stretch
        does. A cut site's piece is empty.
    span
        The span it is part of, or ``None`` for part of a gap a joined item leaves between two.
    starts, ends
        Whether the item itself starts at `start`, and ends at `end`, rather than the stretch
        cutting it there.
    """

    start: int
    end: int
    span: Span | None
    starts: bool
    ends: bool


def items(
    record: SequenceRecord,
    *,
    features: bool = True,
    primers: bool = True,
    cut_sites: bool = True,
    enzymes: str | Iterable[str] | None = None,
    hide_types: Iterable[str] = (),
    source: bool = False,
    translations: bool = False,
) -> tuple[Item, ...]:
    """Return what a map of `record` draws: its features, its primers, then its cut sites.

    Parameters
    ----------
    record
        What is drawn.
    features, primers, cut_sites
        Whether each is drawn.
    enzymes
        The names of the enzymes whose every cut site is drawn; the shipped unique cutters when
        ``None``. How many times an enzyme cuts is counted over the whole record.
    hide_types
        The feature types left off.
    source
        Whether a `source` feature is drawn.
    translations
        Whether each CDS carries its translation, as a sequence view draws it.

    Raises
    ------
    KeyError
        If no shipped enzyme answers to a name in `enzymes`, as `sites.find_sites` raises.
    """
    chosen = None if enzymes is None else _chosen(enzymes)
    left_off = set(hide_types) if source else {*hide_types, "source"}
    drawn = []
    if features:
        drawn.extend(
            _feature(feature, record, translations)
            for feature in record.features
            if feature.type not in left_off
        )
    if primers:
        drawn.extend(
            _primer(primer, site, record)
            for primer in record.primers
            for site in primer.binding_sites
        )
    if cut_sites:
        drawn.extend(_cut_sites(record, chosen))
    return tuple(drawn)


def merge_cuts(
    cuts: Iterable[tuple[str, int]], length: int, *, staggers: Mapping[str, int] | None = None
) -> tuple[Item, ...]:
    """Return a cut site for each position cut, naming every enzyme that cuts there.

    Parameters
    ----------
    cuts
        Each cut as an enzyme's name and the 0-based boundary where it cuts the top strand, one
        for every site, so a name given twice cuts twice.
    length
        The record's.
    staggers
        Each enzyme's `Cutter.stagger`, by name; 0 for a name it does not give.

    Returns
    -------
    tuple[Item, ...]
        In order of position, each naming its enzymes alphabetically, as SnapGene merges them.

    Examples
    --------
    >>> [site.label for site in merge_cuts([("SacI", 406), ("BanII", 406), ("EcoRI", 396)], 2686)]
    ['EcoRI (396)', 'BanII - SacI (406)']
    """
    listed = list(cuts)
    counts = Counter(name for name, _ in listed)
    at: dict[int, set[str]] = {}
    for name, position in listed:
        at.setdefault(position % length, set()).add(name)
    merged = []
    for position, names in sorted(at.items()):
        cutters = tuple(
            Cutter(one, counts[one], (staggers or {}).get(one, 0))
            for one in sorted(names, key=lambda one: (one.casefold(), one))
        )
        name = SEPARATOR.join(cutter.name for cutter in cutters)
        shown = str((position - 1) % length + 1)
        merged.append(
            Item(
                "cut_site",
                name,
                "cut site",
                Strand.NONE,
                (Span(position, position, ENZYME),),
                f"{name} ({shown})",
                {"name": name, "type": "cut site", "span": shown},
                cutters,
            )
        )
    return tuple(merged)


def hiding(item: Item) -> tuple[int, int, int]:
    """Return where an item's label comes in the order labels hide: sort by it, first to hide first.

    Cut sites hide first, those whose enzymes cut most often before the rest, then primers, then
    features, and within each the longest label first. A cut site naming several enzymes hides as
    late as the one among them that cuts least often.

    Examples
    --------
    >>> sites = merge_cuts([("EcoRI", 396), ("BsmBI", 3), ("BsmBI", 45)], 2686)
    >>> [site.label for site in sorted(sites, key=hiding)]
    ['BsmBI (45)', 'BsmBI (3)', 'EcoRI (396)']
    """
    kind = _HIDING.index(item.kind)
    cuts = min((cutter.cuts for cutter in item.cutters), default=0)
    return kind, -cuts, -len(item.label)


def notice(hidden: Iterable[Item]) -> str:
    """Return what a map says of the labels it hid, counting each kind, or nothing if none hid.

    Examples
    --------
    >>> sites = merge_cuts([("EcoRI", 396), ("BamHI", 417), ("SacI", 406)], 2686)
    >>> primer = Item("primer", "M13 fwd", "primer", Strand.FORWARD, (), "M13 fwd (378 .. 394)")
    >>> notice([*sites, primer])
    '3 enzyme sites and 1 primer are hidden'
    >>> notice(sites[:1]), notice([])
    ('1 enzyme site is hidden', '')
    """
    counts = Counter(item.kind for item in hidden)
    parts = [
        f"{counts[kind]} {word if counts[kind] == 1 else word + 's'}"
        for kind, word in zip(_HIDING, ("enzyme site", "primer", "feature"), strict=True)
        if counts[kind]
    ]
    if not parts:
        return ""
    listed = parts[0] if len(parts) == 1 else f"{', '.join(parts[:-1])} and {parts[-1]}"
    verb = "is" if sum(counts.values()) == 1 else "are"
    return f"{listed} {verb} hidden"


def default_color(feature: Feature) -> str:
    """Return the colour a feature with none of its own takes, by the group its type falls in.

    A `primer_bind` feature takes the primers' purple instead. A `regulatory` feature is grouped
    by its `/regulatory_class`: a terminator's or a polyA signal's among the terminators, any other
    among the promoters and enhancers.

    Examples
    --------
    >>> from liulab_mbio.sequence import Segment
    >>> default_color(Feature("ori", "rep_origin", (Segment(0, 10),)))
    '#eedd88'
    """
    if feature.type == "primer_bind":
        return PRIMER.lower()
    if feature.type == "regulatory" and TERMINATING_CLASSES & set(
        map(str, feature.qualifiers.get("regulatory_class", ()))
    ):
        return GROUPS["terminators and polyA"][0].lower()
    for color, types in GROUPS.values():
        if feature.type in types:
            return color.lower()
    return OTHER.lower()


def unwrapped(spans: Sequence[Span | Segment], length: int) -> tuple[tuple[int, int], ...]:
    """Return each span's start and end, counted round a circular record from the first span.

    A span that starts before the one before it ends lies past the origin, so it moves on by the
    record's length, as SnapGene reads a feature's segments.
    """
    counted: list[tuple[int, int]] = []
    for span in spans:
        start = span.start
        if counted and start < counted[-1][1]:
            start += length
        counted.append((start, start + span.end - span.start))
    return tuple(counted)


def pieces(item: Item, start: int, end: int, length: int, *, circular: bool) -> tuple[Piece, ...]:
    """Return what of `item` lies in a stretch of a record `length` bases long, in order along it.

    The stretch is 0-based and half-open, ending past `length` where it runs across the origin of
    a circular record, and at most one turn long. An item it cuts, or one lying across where the
    stretch opens a circular record, comes back in pieces. A cut site at either end of the stretch
    lies in it. The two ends of a whole turn are one boundary, and a cut there lies at the end.

    Examples
    --------
    A feature across the origin of a circular record opened there lies at both ends:

    >>> ori = Item("feature", "ori", "rep_origin", Strand.FORWARD, (Span(90, 110, "#eedd88"),), "")
    >>> opened = pieces(ori, 0, 100, 100, circular=True)
    >>> [(one.start, one.end, one.starts, one.ends) for one in opened]
    [(0, 10, False, True), (90, 100, True, False)]
    """
    counted = unwrapped(item.spans, length)
    runs: list[tuple[int, int, Span | None]] = []
    for index, (span, (low, high)) in enumerate(zip(item.spans, counted, strict=True)):
        if index and low > counted[index - 1][1]:
            runs.append((counted[index - 1][1], low, None))
        runs.append((low, high, span))
    first, last = counted[0][0], counted[-1][1]
    whole = circular and end - start == length
    found = []
    for low, high, span in runs:
        for shift in (-length, 0, length) if circular else (0,):
            if low == high:
                at = low + shift
                if start <= at <= end and not (whole and at == start):
                    found.append(Piece(at, at, span, starts=True, ends=True))
                continue
            clipped = max(low + shift, start), min(high + shift, end)
            if clipped[0] < clipped[1]:
                found.append(
                    Piece(*clipped, span, clipped[0] == first + shift, clipped[1] == last + shift)
                )
    return tuple(sorted(found, key=lambda piece: (piece.start, piece.end)))


def span_text(start: int, end: int, length: int) -> str:
    """Return a span as a person reads it, 1-based and inclusive, across the origin as needed.

    Examples
    --------
    >>> span_text(2683, 2689, 2686)
    '2684 .. 3'
    """
    return f"{start % length + 1} .. {(end - 1) % length + 1}"


def text_color(fill: str) -> str:
    """Return black or white, whichever contrasts more with `fill`, black when they tie.

    Examples
    --------
    >>> text_color("#ffffff"), text_color("#993366")
    ('#000000', '#ffffff')
    """
    luminance = _luminance(fill)
    return "#000000" if (luminance + 0.05) / 0.05 >= 1.05 / (luminance + 0.05) else "#ffffff"


def outline_color(fill: str) -> str:
    """Return the colour an arrow filled with `fill` is outlined in: the fill, darkened.

    Darkened this far, even a white arrow's outline contrasts with a white page.
    """
    red, green, blue = _channels(fill)
    return "#{:02x}{:02x}{:02x}".format(*(round(0.45 * channel) for channel in (red, green, blue)))


def translation(feature: Feature, record: SequenceRecord) -> tuple[Codon, ...]:
    """Return the codons a feature of `record` reads as a CDS, in the order they are read.

    The feature's segments are joined, on its strand, and read from its `/codon_start`, the first
    base when it gives none; bases left over at the end make no codon.

    Examples
    --------
    Read from its second base, across the origin:

    >>> record = SequenceRecord("TAAATGC", topology="circular")
    >>> cds = Feature("orf", "CDS", (Segment(3, 10),), qualifiers={"codon_start": (2,)})
    >>> [(codon.middle, codon.name) for codon in translation(cds, record)]
    [(5, 'Cys'), (1, 'Ter')]
    """
    length = len(record)
    positions = [
        at % length
        for start, end in unwrapped(feature.segments, length)
        for at in range(start, end)
    ]
    bases = "".join(record.sequence[at] for at in positions)
    if feature.strand == Strand.REVERSE:
        positions.reverse()
        bases = reverse_complement(bases)
    first = _codon_start(feature) - 1
    return tuple(
        Codon(positions[at + 1], _amino_acid(bases[at : at + 3]))
        for at in range(first, len(bases) - 2, 3)
    )


def _codon_start(feature: Feature) -> int:
    """Return the base, 1 to 3, a feature's `/codon_start` reads from: 1 when it gives none."""
    try:
        start = int(feature.qualifiers.get("codon_start", (1,))[0])
    except (IndexError, ValueError):
        return 1
    return start if start in (1, 2, 3) else 1


def _amino_acid(codon: str) -> str:
    try:
        return amino_acid(codon)
    except KeyError:
        return "X"


@cache
def _three_letters() -> Mapping[str, str]:
    """Each amino acid's one letter and three, Biopython's, with ``Ter`` for a stop."""
    from Bio.Data.IUPACData import protein_letters_1to3_extended

    return {**protein_letters_1to3_extended, "*": "Ter"}


def _feature(feature: Feature, record: SequenceRecord, translations: bool) -> Item:
    length = len(record)
    own = _given(feature.color)
    default = default_color(feature)
    spans = tuple(
        Span(segment.start, segment.end, _given(segment.color) or own or default)
        for segment in feature.segments
    )
    runs: list[list[int]] = []
    for start, end in unwrapped(spans, length):
        if runs and start == runs[-1][1]:
            runs[-1][1] = end
        else:
            runs.append([start, end])
    bases = sum(segment.end - segment.start for segment in feature.segments)
    hover = {
        "name": feature.name,
        "type": feature.type,
        "span": ", ".join(span_text(start, end, length) for start, end in runs),
        "length": f"{bases} bp",
    }
    codons = translation(feature, record) if translations and feature.type == "CDS" else ()
    return Item(
        "feature",
        feature.name,
        feature.type,
        feature.strand,
        spans,
        feature.name,
        hover,
        translation=codons,
    )


def _primer(primer: Primer, site: BindingSite, record: SequenceRecord) -> Item:
    length = len(record)
    span = span_text(site.start, site.end, length)
    hover = {
        "name": primer.name,
        "type": "primer",
        "span": span,
        "length": f"{site.end - site.start} bp",
    }
    return Item(
        "primer",
        primer.name,
        "primer",
        site.strand,
        (Span(site.start, site.end, PRIMER.lower()),),
        f"{primer.name} ({span})",
        hover,
        tail=max(0, len(primer.sequence) - (site.end - site.start)),
        mismatches=_mismatches(primer, site, record),
    )


def _mismatches(primer: Primer, site: BindingSite, record: SequenceRecord) -> tuple[int, ...]:
    """Return where a primer's bases differ from the record's they pair with at `site`.

    The primer's 3' end pairs with the site's 3' end on the primer's strand, and each base 5' of it
    with the site's base as far from that end, as far as the site or the primer runs.
    """
    length = len(record)
    pairing = primer.sequence[-(site.end - site.start) :]
    if site.strand == Strand.REVERSE:
        reads, first = reverse_complement(pairing), site.start
    else:
        reads, first = pairing, site.end - len(pairing)
    return tuple(
        (first + index) % length
        for index, base in enumerate(reads)
        if base != record.sequence[(first + index) % length]
    )


def _chosen(names: str | Iterable[str]) -> tuple[Enzyme, ...]:
    """Return the shipped enzymes `names` answer to, each once however many names it has."""
    named = (names,) if isinstance(names, str) else names
    return tuple({enzyme.name: enzyme for enzyme in map(get_enzyme, named)}.values())


def _cut_sites(record: SequenceRecord, chosen: tuple[Enzyme, ...] | None) -> tuple[Item, ...]:
    """Return the cut sites of `chosen`, or of the shipped unique cutters when that is ``None``."""
    # A site whose cut falls off a linear record cuts nothing there.
    found = [
        site for site in find_sites(record, shipped() if chosen is None else chosen) if site.cuts
    ]
    counts = Counter(site.enzyme.name for site in found)
    return merge_cuts(
        (
            (site.enzyme.name, site.top_cut)
            for site in found
            if chosen is not None or counts[site.enzyme.name] == 1
        ),
        len(record),
        staggers={site.enzyme.name: site.enzyme.bottom_cut - site.enzyme.top_cut for site in found},
    )


def _given(color: str | None) -> str | None:
    """Return a file's colour lower-cased, or `None` for none, SnapGene's `noColor` included."""
    return color.lower() if color is not None and _HEX.fullmatch(color) else None


def _channels(color: str) -> tuple[int, int, int]:
    return int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)


def _luminance(color: str) -> float:
    """WCAG 2's relative luminance."""
    linear = [
        value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4
        for value in (channel / 255 for channel in _channels(color))
    ]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]
