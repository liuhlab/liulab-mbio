"""Design the barcodes of one part list, and check a set someone already holds.

A barcode names one part, so a set has to hold three things at once: no two barcodes of one part
list may read as one another, none may spell a site a cloning enzyme reads, and none may put a
stop in the frame the construct is translated in. Distance is checked **within a part list**,
which is where two barcodes are told apart and what the published set holds.

The composition dials are the other half, and `docs/research/barcode-design.md` is the evidence
for each: a GC band is convention that no measurement supports at this length, while a
homopolymer cap covers an indel no distance rule can see. Candidates are drawn in a seeded
random order rather than walked in order, which is what keeps a set from opening on a run of one
base, and the seed is a parameter, so one request returns one set.
"""

import random
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from dataclasses import KW_ONLY, dataclass
from itertools import combinations, groupby

from liulab_mbio.codons import amino_acid
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import SequenceRecord
from liulab_mbio.sites import EnzymeLike, find_sites

#: The fewest mismatches two barcodes of one part list stand apart. Measured on the published
#: set in `docs/research/barcode-design.md`.
MIN_DISTANCE = 3

#: The longest run of one base a barcode may carry, and the one composition dial on by default:
#: a homopolymer indel is the dominant error of the long reads a barcode block is read with, and
#: a distance rule cannot see an indel whatever its distance. `docs/research/barcode-design.md`
#: weighs both. The number itself is a judgement and not a measurement — one base longer would
#: reject nothing that has been published and is equally defensible.
MAX_HOMOPOLYMER = 5

#: No GC band, because no measurement supports one at this length: every band in the barcode
#: literature descends from a single uncited sentence, and the one study that imposed a band and
#: then measured its consequences found GC was not what separated its probes.
#: `docs/research/barcode-design.md` traces both, and the vendors document no band at all.
GC_BAND: tuple[float, float] | None = None

#: The seed a design draws with when the caller names none.
SEED = 0

#: The bases a barcode is spelled from, and the order an index spells them in.
BASES = "ACGT"

#: How many candidates in a row may be rejected before a design gives up on the space.
_TRIES = 4096


class SpaceExhaustedError(ValueError):
    """No further barcode could be found under the rules asked for.

    A `ValueError`, so a caller catching that catches this. The message says how many were
    found and which rule rejected the most candidates.
    """


@dataclass(frozen=True, slots=True)
class BarcodeRules:
    """What every barcode of one part list holds to, checked as it is read.

    Parameters
    ----------
    length
        How many bases name one part.
    scar
        The cloning scar joining one barcode to the next, ``""`` where a barcode stands alone.
        A barcode is read with the scar on either side of it, so a forbidden site or a stop
        spanning that junction is found rather than missed.
    phase
        The reading frame at a barcode's first base: how many bases of that codon were read
        before it, 0, 1 or 2. ``None`` where the construct never translates the barcode, which
        takes away the stop-codon rule and the whole-codon rule together.
    forbidden
        Enzymes whose recognition sites no barcode may spell, on either strand. Each an
        `Enzyme` or a name `liulab_mbio.enzymes.get_enzyme` answers to.
    distance
        The fewest mismatches two barcodes of one part list may stand apart.
    max_homopolymer
        The longest run of one base a barcode may carry, or ``None`` for no cap.
    gc_band
        The share of G and C a barcode must hold, low and high inclusive, each 0 to 1. ``None``
        is no band at all, which is the default.

    Raises
    ------
    ValueError
        Naming what cannot be designed under: a length or a distance that is not positive, a
        scar that is not definite bases, a phase outside 0 to 2, a barcode and scar that are not
        a whole number of codons together, or a GC band outside 0 to 1.
    KeyError
        If a forbidden enzyme is not one the package ships.
    """

    length: int
    _: KW_ONLY
    scar: str = ""
    phase: int | None = 0
    forbidden: Iterable[EnzymeLike] = ()
    distance: int = MIN_DISTANCE
    max_homopolymer: int | None = MAX_HOMOPOLYMER
    gc_band: tuple[float, float] | None = GC_BAND

    def __post_init__(self) -> None:
        """Normalise the scar, then refuse rules no set could be designed under."""
        object.__setattr__(self, "scar", self.scar.upper())
        object.__setattr__(self, "forbidden", tuple(self.forbidden))
        if self.length <= 0:
            raise ValueError(f"a barcode needs a positive length, got {self.length}")
        if self.distance < 1:
            raise ValueError(f"two barcodes stand at least one mismatch apart, got {self.distance}")
        if bad := sorted(set(self.scar) - set(BASES)):
            raise ValueError(
                f"a cloning scar needs definite bases, and {''.join(bad)} is not one of ACGT"
            )
        if self.max_homopolymer is not None and self.max_homopolymer < 1:
            raise ValueError(
                f"a homopolymer cap is at least one base, got {self.max_homopolymer}; "
                "pass None for no cap"
            )
        if self.gc_band is not None:
            low, high = self.gc_band
            if not 0 <= low <= high <= 1:
                raise ValueError(
                    f"a GC band runs low to high inside 0 to 1, got {low} to {high}; "
                    "pass None for no band"
                )
        self._check_frame()
        _resolve(self.forbidden)

    @property
    def enzymes(self) -> tuple[Enzyme, ...]:
        """The forbidden enzymes, each read into a record."""
        return _resolve(self.forbidden)

    def _check_frame(self) -> None:
        """Refuse a barcode the construct could not read as codons."""
        if self.phase is None:
            return
        if self.phase not in (0, 1, 2):
            raise ValueError(
                f"a phase is 0, 1 or 2 bases into a codon, got {self.phase}; pass None where "
                "the construct does not translate the barcode"
            )
        unit = self.length + len(self.scar)
        if unit % 3:
            raise ValueError(
                f"a barcode of {self.length} bases and a cloning scar of {len(self.scar)} make "
                f"{unit}, which is not a whole number of codons"
            )
        if self.phase > len(self.scar):
            raise ValueError(
                f"a phase of {self.phase} puts {self.phase} base(s) before the barcode in its "
                f"first codon, and the {len(self.scar)}-base cloning scar does not hold them"
            )


def design_barcodes(count: int, rules: BarcodeRules, *, seed: int = SEED) -> tuple[str, ...]:
    """Design `count` barcodes for one part list, holding every rule of `rules`.

    Candidates are drawn from the whole space in a seeded random order and kept when they pass,
    so one seed returns one set. Walking the space in order instead opens the set on a run of
    one base and leaves it poor in G and C, which no filter here repairs — the reason the draw
    is shuffled is measured in `docs/research/barcode-design.md`.

    Parameters
    ----------
    count
        How many barcodes the part list needs.
    rules
        What each barcode, and each pair of them, must hold to.
    seed
        The seed the draw is made with.

    Returns
    -------
    tuple of str
        The barcodes, 5' to 3', in the order they were accepted.

    Raises
    ------
    SpaceExhaustedError
        If no further barcode can be found, naming which rule rejected the most candidates.
    ValueError
        If `count` is negative.
    KeyError
        If a forbidden enzyme is not one the package ships.

    Examples
    --------
    >>> design_barcodes(3, BarcodeRules(6))
    ('TACCAT', 'TCCTCC', 'ACCAGT')
    >>> design_barcodes(3, BarcodeRules(6)) == design_barcodes(3, BarcodeRules(6))
    True
    """
    if count < 0:
        raise ValueError(f"a part list needs no fewer than no barcodes, got {count}")
    enzymes = rules.enzymes
    space = 4**rules.length
    generator = random.Random(seed)
    chosen: list[str] = []
    drawn: set[int] = set()
    rejected: Counter[str] = Counter()
    misses = 0
    while len(chosen) < count:
        if len(drawn) == space or misses >= _TRIES:
            raise SpaceExhaustedError(_exhausted(count, chosen, rejected, rules, drawn, space))
        index = generator.randrange(space)
        if index in drawn:
            continue
        drawn.add(index)
        barcode = _spell(index, rules.length)
        broken = _problem(barcode, rules, enzymes)
        if broken is None and (near := _nearest(barcode, chosen, rules.distance)) is not None:
            broken = near
        if broken is None:
            chosen.append(barcode)
            misses = 0
        else:
            rejected[broken[0]] += 1
            misses += 1
    return tuple(chosen)


def check_barcodes(barcodes: Iterable[str], rules: BarcodeRules) -> tuple[str, ...]:
    """Return one message per rule these barcodes break, empty where the set holds.

    Each barcode is judged on its own and each pair on its distance, so a set from anywhere —
    a published one, or one a user wrote by hand — is read by the rules a design works to.

    Examples
    --------
    >>> check_barcodes(["TACCAT", "TCCTCC"], BarcodeRules(6))
    ()
    >>> check_barcodes(["AAAAAA", "AAAAAC"], BarcodeRules(6))[0]
    'AAAAAA carries a run of 6 A, over the cap of 5'
    """
    held = tuple(barcode.upper() for barcode in barcodes)
    enzymes = rules.enzymes
    problems = [
        f"{barcode} {broken[1]}"
        for barcode in held
        if (broken := _problem(barcode, rules, enzymes)) is not None
    ]
    problems.extend(
        f"{one} and {other} {_apart(_distance(one, other), rules.distance)}"
        for one, other in combinations(held, 2)
        if len(one) == len(other) and _distance(one, other) < rules.distance
    )
    return tuple(problems)


def _problem(
    barcode: str, rules: BarcodeRules, enzymes: tuple[Enzyme, ...]
) -> tuple[str, str] | None:
    """Which rule this one barcode breaks, as a key and a phrase, or ``None`` where it holds.

    The first rule broken is the one reported, cheapest and most telling first, so that a
    refusal names a filter rejecting everything rather than the distance rule that follows it.
    """
    if len(barcode) != rules.length:
        return "length", f"is {len(barcode)} bases, not the {rules.length} asked for"
    if bad := sorted(set(barcode) - set(BASES)):
        return "bases", f"holds {''.join(bad)}, which is not one of ACGT"
    if rules.max_homopolymer is not None:
        base, run = _longest_run(barcode)
        if run > rules.max_homopolymer:
            return (
                "homopolymer",
                f"carries a run of {run} {base}, over the cap of {rules.max_homopolymer}",
            )
    if rules.gc_band is not None:
        low, high = rules.gc_band
        share = sum(base in "GC" for base in barcode) / len(barcode)
        if not low <= share <= high:
            return (
                "gc",
                f"is {share:.0%} GC, outside the band {low:.0%} to {high:.0%}",
            )
    for codon in _codons(barcode, rules):
        if amino_acid(codon) == "*":
            return "stop", f"spells {codon}, a stop, in the frame the construct reads"
    if enzymes and (site := _site(barcode, rules, enzymes)) is not None:
        return site
    return None


def _site(barcode: str, rules: BarcodeRules, enzymes: tuple[Enzyme, ...]) -> tuple[str, str] | None:
    """Whether a forbidden enzyme reads a site in this barcode, the cloning scar included.

    The barcode is searched with the scar on either side of it, so a site spanning the junction
    the ligation makes is found. A site needing bases from two different barcodes is not: those
    are neighbours only once the block is built.
    """
    found = find_sites(SequenceRecord(rules.scar + barcode + rules.scar), enzymes)
    if not found:
        return None
    site = found[0]
    crosses = bool(rules.scar) and not (
        site.start >= len(rules.scar) and site.end <= len(rules.scar) + rules.length
    )
    where = " across the junction with the cloning scar" if crosses else ""
    return "site", f"spells a {site.enzyme.name} site{where}"


def _codons(barcode: str, rules: BarcodeRules) -> list[str]:
    """Every codon the construct reads that this barcode puts a base in.

    A barcode and its scar are a whole number of codons, so the block repeats one frame: the
    bases before a barcode's first codon are the end of the scar that joins it to the last.
    """
    if rules.phase is None:
        return []
    lead = rules.scar[len(rules.scar) - rules.phase :] if rules.phase else ""
    unit = lead + barcode + rules.scar
    return [unit[at : at + 3] for at in range(0, len(unit) - rules.phase, 3)]


def _nearest(barcode: str, chosen: Iterable[str], distance: int) -> tuple[str, str] | None:
    """Whether this barcode stands too close to one already chosen."""
    for other in chosen:
        if (apart := _distance(barcode, other)) < distance:
            return "distance", _apart(apart, distance)
    return None


def _apart(apart: int, distance: int) -> str:
    """How two barcodes standing too close are described."""
    return f"stand {apart} mismatch(es) apart, under the {distance} one part list needs"


def _distance(one: str, other: str) -> int:
    """How many positions two barcodes of one length differ in."""
    return sum(a != b for a, b in zip(one, other, strict=True))


def _longest_run(barcode: str) -> tuple[str, int]:
    """Return the base with the longest run in this barcode, and how long that run is."""
    runs = [(base, len(tuple(same))) for base, same in groupby(barcode)]
    base, run = max(runs, key=lambda one: one[1])
    return base, run


def _spell(index: int, length: int) -> str:
    """Spell one index of the ``4 ** length`` space as bases."""
    spelled = []
    for _ in range(length):
        index, remainder = divmod(index, 4)
        spelled.append(BASES[remainder])
    return "".join(reversed(spelled))


def _resolve(forbidden: Iterable[EnzymeLike]) -> tuple[Enzyme, ...]:
    """Read each forbidden enzyme, by name or by record, refusing a name nothing answers to."""
    return tuple(one if isinstance(one, Enzyme) else get_enzyme(one) for one in forbidden)


#: How each rule is named in a refusal, so that whoever reads one knows which dial to turn.
_WORDING: Mapping[str, Callable[[BarcodeRules], str]] = {
    "distance": lambda rules: (
        f"the distance rule (at least {rules.distance} mismatches within a part list)"
    ),
    "homopolymer": lambda rules: f"the homopolymer cap (no run over {rules.max_homopolymer})",
    "gc": lambda rules: (
        f"the GC band ({rules.gc_band[0]:.0%} to {rules.gc_band[1]:.0%})"
        if rules.gc_band is not None
        else "the GC band"
    ),
    "site": lambda rules: f"the forbidden site(s) ({', '.join(one.name for one in rules.enzymes)})",
    "stop": lambda _: "the stop-codon rule",
    "length": lambda rules: f"the length rule ({rules.length} bases)",
    "bases": lambda _: "the definite-bases rule",
}


def _exhausted(
    count: int,
    chosen: list[str],
    rejected: Counter[str],
    rules: BarcodeRules,
    drawn: set[int],
    space: int,
) -> str:
    """Say how far the design got and which rule took the space away from it."""
    ranked = sorted(rejected.items(), key=lambda entry: (-entry[1], entry[0]))
    named = ", ".join(f"{_WORDING[key](rules)} rejected {number}" for key, number in ranked)
    where = (
        f"every one of the {space} barcode(s) the space holds was drawn"
        if len(drawn) == space
        else f"no further barcode passed in {_TRIES} candidates"
    )
    return (
        f"{count} barcodes of {rules.length} bases were asked for and {len(chosen)} found: "
        f"{where}, and of {sum(rejected.values())} candidate(s) rejected, "
        f"{named or 'none was rejected by any rule'}"
    )
