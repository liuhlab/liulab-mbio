"""An oligo pool as a vendor takes it: one oligo a fragment, padded to one length, and its bill.

A fragment of a split cargo is not yet an oligo. It needs the primer sites one PCR pulls it out
by, the recognition sites that cut it back out of the amplicon, and filler, because a vendor
synthesises a pool at one length and holds every member to a narrow band around it. This builds
that oligo, keeps what it came from on it, and writes the sheet an order is placed from.

The filler is screened rather than trusted. It sits outboard of the 3' cut, so it never enters
the product, but a recognition site inside it is still cut in the tube and an oligo spelling an
amplification primer still mis-primes. So the padding is drawn base by base against both, and
the finished oligo is read once more as a whole, which is what catches a motif the join spells.

**Price steers nothing here.** The pool reports its count, the band each quantity falls in and
the headroom to that band's top edge, and `liulab_mbio.bench.prices` prices it where the user
holds a tariff. Headroom is reported, never optimised.
"""

import random
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import KW_ONLY, dataclass
from decimal import Decimal

from liulab_mbio.bench.prices import Band, Headroom, Item
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
    span_text,
)
from liulab_mbio.sites import EnzymeLike, find_sites
from liulab_mbio.split import Budget, CargoSplit, Fragment

#: The columns of the pool order sheet. `source` and `fragment` are the provenance: which cargo
#: an oligo carries and which piece of it, so an oligo on a plate is traceable to its protein.
POOL_COLUMNS = (
    "name",
    "sequence",
    "length",
    "source",
    "fragment",
    "fragments",
    "span",
    "overhangs",
    "primers",
    "pad",
)

#: The columns of the primer inventory beside a pool: what each primer does and what to order.
PRIMER_COLUMNS = ("name", "role", "sequence", "length", "oligos")

#: How far apart a pool's longest and shortest oligo may stand, as a share of the longest, before
#: a vendor bins it. Padding to one length leaves this at zero, which is the point of padding.
LENGTH_TOLERANCE = 0.15


@dataclass(frozen=True, slots=True)
class PrimerSite:
    """One amplification primer written onto every oligo it serves.

    Parameters
    ----------
    role
        What the primer does, in the method's own words.
    name
        What it is called, so an inventory and an oligo agree.
    sequence
        The primer itself, 5' to 3' as it is ordered.
    """

    role: str
    name: str
    sequence: str

    def __post_init__(self) -> None:
        """Hold the sequence to the one case the rest of this package writes."""
        object.__setattr__(self, "sequence", self.sequence.upper())

    def __len__(self) -> int:
        """How many bases the oligo spends on it."""
        return len(self.sequence)


#: What stands between a Type IIS site and the span it cuts, where no caller names another. How
#: long it is the enzyme fixes; which base it spells nothing does, so one base serves everywhere.
SPACER = "A"


@dataclass(frozen=True, slots=True)
class OligoLayout:
    """Where everything sits on one oligo, and what that leaves for the cargo.

    An oligo reads, 5' to 3':
    ``[forward primers][site][spacer][span][spacer][site reversed][filler][reverse primers]``.
    Both recognition sites face inward, so the cuts fall at the span's two ends and the filler
    is outboard of the 3' one.

    Parameters
    ----------
    length
        Every base the vendor synthesises, the one length the pool is padded to.
    enzyme
        The Type IIS enzyme cutting the span out.
    primers
        How many primer sites the oligo carries.
    primer_length
        How long each is.
    spacer
        Between a recognition site and the span it cuts.

    Raises
    ------
    ValueError
        If the overhead leaves no span, or the spacer does not reach the enzyme's cut.
    """

    length: int
    _: KW_ONLY
    enzyme: EnzymeLike
    primers: int = 3
    primer_length: int = 20
    spacer: str = SPACER

    def __post_init__(self) -> None:
        """Check the spacer puts the cut where the span begins."""
        reach = self.cutter.top_cut - len(self.cutter.site)
        if len(self.spacer) != reach:
            raise ValueError(
                f"{self.cutter.name} cuts {reach} base(s) past its site, so the spacer is "
                f"{reach} base(s) long and not {len(self.spacer)}"
            )

    @property
    def cutter(self) -> Enzyme:
        """The enzyme itself."""
        return get_enzyme(self.enzyme) if isinstance(self.enzyme, str) else self.enzyme

    @property
    def overhead(self) -> int:
        """Every base of the oligo that is not cargo: the primer sites, the sites and spacers."""
        return self.primers * self.primer_length + 2 * (len(self.cutter.site) + len(self.spacer))

    @property
    def budget(self) -> Budget:
        """What one oligo of this layout carries, which is what `split_cargo` is given."""
        return Budget(self.length, self.overhead)


@dataclass(frozen=True, slots=True)
class Oligo:
    """One member of the pool, and what it came from.

    Parameters
    ----------
    name
        What the vendor's plate calls it.
    sequence
        Every base ordered, 5' to 3'.
    source
        The cargo this is a piece of.
    fragment, fragments
        Which piece, counting from one, and of how many.
    span
        Where that piece lies in the cargo, its overhangs included.
    cargo_length
        How long that cargo is, which is what the span is read against.
    overhangs
        What its two cuts leave, 5' then 3'.
    primers
        The primers written onto it, in the order the oligo spells them.
    pad
        How many bases of filler it carries.
    """

    name: str
    sequence: str
    _: KW_ONLY
    source: str
    fragment: int
    fragments: int
    span: Segment
    cargo_length: int
    overhangs: tuple[str, str]
    primers: tuple[str, ...]
    pad: int

    def __len__(self) -> int:
        """How many bases are ordered."""
        return len(self.sequence)


@dataclass(frozen=True, slots=True)
class Pool:
    """Every oligo one order holds.

    Parameters
    ----------
    name
        What the order is called.
    layout
        What every oligo of it was built to.
    oligos
        In the order they are listed.
    primers
        Every primer the pool is amplified by, once each.
    """

    name: str
    _: KW_ONLY
    layout: OligoLayout
    oligos: tuple[Oligo, ...]
    primers: tuple[PrimerSite, ...] = ()

    @property
    def count(self) -> int:
        """How many oligos are ordered."""
        return len(self.oligos)

    @property
    def spread(self) -> float:
        """How far the longest and shortest oligo stand apart, over the longest.

        Zero where every oligo is padded to one length, which is what a vendor's uniformity
        band asks for.
        """
        if not self.oligos:
            return 0.0
        lengths = [len(one) for one in self.oligos]
        return (max(lengths) - min(lengths)) / max(lengths)

    @property
    def sources(self) -> tuple[str, ...]:
        """Each cargo the pool carries, once, in the order the oligos list them."""
        return tuple(dict.fromkeys(one.source for one in self.oligos))

    def fragment_counts(self) -> Mapping[int, int]:
        """How many cargoes took each number of fragments, which is what success is read on."""
        counted: dict[int, int] = {}
        for source in self.sources:
            pieces = next(one.fragments for one in self.oligos if one.source == source)
            counted[pieces] = counted.get(pieces, 0) + 1
        return dict(sorted(counted.items()))


def pad_bases(
    before: str,
    after: str,
    count: int,
    *,
    avoid: Iterable[EnzymeLike] = (),
    avoid_bases: Iterable[str] = (),
    seed: int = 0,
) -> str:
    """Draw `count` filler bases that spell nothing forbidden, in the context they sit in.

    The bases either side are given, so a motif the join creates is refused as surely as one
    inside the filler. The draw is seeded, so the same inputs give the same padding.

    Raises
    ------
    ValueError
        If no base is left at some position, which names how far it got.

    Examples
    --------
    >>> len(pad_bases("AAAA", "TTTT", 20, avoid=["BsaI"], seed=1))
    20
    """
    if count <= 0:
        return ""
    motifs = _motifs(avoid, avoid_bases)
    reach = max((len(one) for one in motifs), default=1)
    draw = random.Random(seed)
    filler = ""
    for place in range(count):
        for base in draw.sample("ACGT", 4):
            if not _spells(before + filler + base, after, motifs, reach):
                filler += base
                break
        else:
            raise ValueError(
                f"no base spells a clean filler at position {place} of {count}: every one of "
                "ACGT finishes a forbidden motif"
            )
    return filler


def build_oligo(
    split: CargoSplit,
    fragment: Fragment,
    *,
    name: str,
    source: str,
    layout: OligoLayout,
    forward: Sequence[PrimerSite],
    reverse: Sequence[PrimerSite],
    avoid: Iterable[EnzymeLike] = (),
    seed: int = 0,
) -> Oligo:
    """Build one oligo of the pool from one fragment, padded to the layout's one length.

    The reverse primers are written onto the top strand reverse-complemented, which is how an
    oligo spells a primer that reads the other way.

    Raises
    ------
    ValueError
        If the primer sites are not the layout's, if the fragment does not fit, or if the
        finished oligo spells a forbidden site anywhere but its two designed ones.
    """
    site = layout.cutter.site
    head = "".join(one.sequence for one in forward)
    tail = "".join(reverse_complement(one.sequence) for one in reverse)
    if len(forward) + len(reverse) != layout.primers:
        raise ValueError(
            f"this layout carries {layout.primers} primer site(s) and {len(forward)} forward "
            f"with {len(reverse)} reverse were given"
        )
    span = split.sequence(fragment)
    opening = head + site + layout.spacer
    closing = reverse_complement(layout.spacer) + reverse_complement(site)
    pad = layout.length - len(opening) - len(span) - len(closing) - len(tail)
    if pad < 0:
        raise ValueError(
            f"fragment {fragment.index + 1} of {source} is {len(span)} bases, {-pad} more than "
            f"an oligo of {layout.length} nt holds at this layout"
        )
    filler = pad_bases(
        opening + span + closing,
        tail,
        pad,
        avoid=avoid,
        avoid_bases=[one.sequence for one in (*forward, *reverse)],
        seed=seed,
    )
    sequence = opening + span + closing + filler + tail
    _check(sequence, span, name, layout, avoid)
    return Oligo(
        name,
        sequence,
        source=source,
        fragment=fragment.index + 1,
        fragments=split.pieces,
        span=fragment.span,
        cargo_length=len(split.cargo),
        overhangs=fragment.overhangs,
        primers=tuple(one.name for one in (*forward, *reverse)),
        pad=pad,
    )


def oligo_record(
    oligo: Oligo, *, layout: OligoLayout, primers: Sequence[PrimerSite]
) -> SequenceRecord:
    """One oligo as a record: its cargo, its filler, and every primer at the site it reads from.

    What a map of an oligo is drawn over. Each primer is found where the oligo spells it, a
    forward one along the top strand and a reverse one reverse-complemented at the far end, so
    the record cannot disagree with the bases ordered. The two recognition sites are left to a
    map's cut sites, which draw where the enzyme cuts and not only where it binds.

    Parameters
    ----------
    oligo
        The member of the pool to draw.
    layout
        What it was built to, which says how long a primer site is.
    primers
        Every primer the pool is amplified by: `Pool.primers` serves.

    Raises
    ------
    ValueError
        If a primer the oligo names is not in `primers`, or is not spelt where the layout puts
        it, which says the oligo and the layout disagree.

    Examples
    --------
    >>> from liulab_mbio.split import split_cargo
    >>> cargo = SequenceRecord("ATG" + "ACGTTGCA" * 6, name="block")
    >>> layout = OligoLayout(120, enzyme="BsaI", primers=2, primer_length=20)
    >>> sites = [PrimerSite("forward", "F1", "A" * 20), PrimerSite("inner", "R1", "C" * 20)]
    >>> split = split_cargo(cargo, layout.cutter, budget=layout.budget)
    >>> one = build_oligo(
    ...     split, split.fragments[0], name="block_f1", source="block", layout=layout,
    ...     forward=sites[:1], reverse=sites[1:],
    ... )
    >>> [p.name for p in oligo_record(one, layout=layout, primers=sites).primers]
    ['F1', 'R1']
    """
    known = {one.name: one for one in primers}
    missing = [name for name in oligo.primers if name not in known]
    if missing:
        listed = ", ".join(repr(name) for name in missing)
        raise ValueError(f"oligo {oligo.name!r} names {listed}, which the pool does not carry")
    sites = [known[name] for name in oligo.primers]
    leading = 0
    while leading < len(sites) and oligo.sequence.startswith(
        "".join(one.sequence for one in sites[: leading + 1])
    ):
        leading += 1
    drawn, at = [], 0
    for one in sites[:leading]:
        drawn.append(_bound(one, at, Strand.FORWARD))
        at += len(one)
    tail = sites[leading:]
    at = len(oligo.sequence) - sum(len(one) for one in tail)
    for one in tail:
        spelt = reverse_complement(one.sequence)
        if oligo.sequence[at : at + len(one)] != spelt:
            raise ValueError(
                f"oligo {oligo.name!r} does not spell {one.name!r} at base {at + 1}, where this "
                "layout puts it"
            )
        drawn.append(_bound(one, at, Strand.REVERSE))
        at += len(one)
    return SequenceRecord(
        oligo.sequence,
        name=oligo.name,
        features=_oligo_features(oligo, layout, leading),
        primers=tuple(drawn),
        notes={"description": f"{oligo.source}, fragment {oligo.fragment} of {oligo.fragments}"},
    )


def _bound(site: PrimerSite, at: int, strand: Strand) -> Primer:
    """Return the primer bound where the oligo spells it, described by the role it serves."""
    return Primer(
        site.name,
        site.sequence,
        binding_sites=(BindingSite(at, at + len(site), strand),),
        description=site.role,
    )


def _oligo_features(oligo: Oligo, layout: OligoLayout, leading: int) -> tuple[Feature, ...]:
    """Return the cargo this oligo carries and the filler padding it, where the layout puts them."""
    flank = len(layout.cutter.site) + len(layout.spacer)
    start = leading * layout.primer_length + flank
    trailing = (layout.primers - leading) * layout.primer_length
    end = len(oligo.sequence) - trailing - oligo.pad - flank
    found = [
        Feature(
            f"{oligo.source} fragment {oligo.fragment}",
            "misc_feature",
            (Segment(start, end),),
            strand=Strand.FORWARD,
        )
    ]
    if oligo.pad:
        found.append(
            Feature("filler", "misc_feature", (Segment(end + flank, end + flank + oligo.pad),))
        )
    return tuple(found)


def pool_sheet(pool: Pool) -> str:
    """Return the pool as a tab-separated order sheet, one row an oligo."""
    lines = ["\t".join(POOL_COLUMNS)]
    for one in pool.oligos:
        lines.append(
            "\t".join(
                (
                    one.name,
                    one.sequence,
                    str(len(one)),
                    one.source,
                    str(one.fragment),
                    str(one.fragments),
                    span_text(one.span.start, one.span.end, one.cargo_length),
                    ",".join(one.overhangs),
                    ",".join(one.primers),
                    str(one.pad),
                )
            )
        )
    return "\n".join(lines) + "\n"


def primer_inventory(pool: Pool) -> str:
    """Return every primer the pool is amplified by, with its role, its bases and its reach."""
    lines = ["\t".join(PRIMER_COLUMNS)]
    for one in pool.primers:
        serves = sum(1 for oligo in pool.oligos if one.name in oligo.primers)
        lines.append("\t".join((one.name, one.role, one.sequence, str(len(one)), str(serves))))
    return "\n".join(lines) + "\n"


def headroom(quantity: str, value: float, bands: Iterable[Band]) -> Headroom | None:
    """How far `value` sits from the top of the band holding it, or ``None`` for no band.

    Examples
    --------
    >>> str(headroom("length", 350, [Band("length", Decimal(301), Decimal(350))]))
    '350 length, no slack above it at all'
    """
    asked = Decimal(str(value))
    found = next((band for band in bands if band.holds(asked)), None)
    return None if found is None else Headroom(quantity, asked, found)


def pool_item(
    pool: Pool,
    *,
    key: str,
    bands: Mapping[str, Sequence[Band]] | None = None,
    item: str = "",
) -> Item:
    """Return the pool as one line of a bill, with the band each of its quantities falls in.

    The quantities are the count and the oligo length, which are the two a vendor bands a pool
    by. They compute whether or not anyone holds a tariff, so the band and the headroom are
    reported either way and the money cell is left to `liulab_mbio.bench.prices.bill`.
    """
    quantities = {"count": float(pool.count), "length": float(pool.layout.length)}
    stated = bands or {}
    found = tuple(
        gap
        for gap in (
            headroom(name, value, stated.get(name, ())) for name, value in quantities.items()
        )
        if gap is not None
    )
    return Item(
        item or f"{pool.name} oligo pool",
        pool.count,
        unit="oligos",
        key=key,
        quantities=quantities,
        units=pool.count,
        headroom=found,
    )


def _motifs(avoid: Iterable[EnzymeLike], avoid_bases: Iterable[str]) -> tuple[str, ...]:
    """Every run of bases filler must not spell, both orientations of each, once each."""
    sites = [(get_enzyme(one) if isinstance(one, str) else one).site.upper() for one in avoid]
    both = [*sites, *(one.upper() for one in avoid_bases)]
    return tuple(dict.fromkeys(one for each in both for one in (each, reverse_complement(each))))


def _spells(prefix: str, suffix: str, motifs: Sequence[str], reach: int) -> bool:
    """Whether the last base of `prefix` completes a motif, with `suffix` following it.

    Only a motif covering that base is new, so the bases already settled either side are context
    and never a refusal of their own -- which is what lets the filler sit against a recognition
    site the layout put there on purpose.
    """
    place = len(prefix) - 1
    start = max(0, place - reach + 1)
    window = prefix[start:] + suffix[: reach - 1]
    at = place - start
    for motif in motifs:
        size = len(motif)
        for index in range(max(0, at - size + 1), min(at, len(window) - size) + 1):
            if window[index : index + size] == motif:
                return True
    return False


def _check(
    sequence: str, span: str, name: str, layout: OligoLayout, avoid: Iterable[EnzymeLike]
) -> None:
    """Refuse an oligo cut anywhere but its two designed sites, or dressed into a new one.

    What the span already carries is the cargo's own and was settled before the split, so this
    counts sites rather than forbidding them: a site the layout added is one the span did not
    have, wherever the padding, a primer or a join spelled it.

    Raises
    ------
    ValueError
        Naming the oligo and what reads it.
    """
    record = SequenceRecord(sequence)
    cutter = layout.cutter
    found = find_sites(record, cutter)
    if len(found) != 2:
        raise ValueError(
            f"oligo {name!r} carries {len(found)} {cutter.name} site(s) and a layout puts two "
            "there, one facing each end of the span"
        )
    others = tuple(
        one
        for one in (get_enzyme(each) if isinstance(each, str) else each for each in avoid)
        if one.name != cutter.name
    )
    if not others:
        return
    gained = len(find_sites(record, others)) - len(find_sites(SequenceRecord(span), others))
    if gained:
        raise ValueError(
            f"dressing fragment {name!r} as an oligo spells {gained} site(s) its span did not "
            "carry, which would cut the pool"
        )
