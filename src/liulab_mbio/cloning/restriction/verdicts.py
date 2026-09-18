"""The verdicts a restriction and ligation plan carries, and the one nothing sourced can give.

`liulab_mbio.cloning.restriction.plan` hangs these off its `Plan` and takes the worst of them as
its status. `liulab_mbio.checks` states the rule every one of them follows: a check carries the
value it judged, and a check no sourced threshold judges carries **no verdict** rather than a
pass. Each threshold below is NEB's, through ``docs/research/restriction-ligation.md``, and says
which section it was read from.

The buffer is the one nothing here can judge. Two products supplied in one buffer digest
together in it, which is NEB's own rule; a mixed pair takes the activity percentages, which no
shippable source states for the enzymes this package ships. So that check carries no verdict,
names the supplier's own table, and the protocol page shows it as unjudged.
"""

from collections.abc import Mapping, Sequence

from liulab_mbio.bench.amounts import Amount
from liulab_mbio.bench.inactivation import heat_inactivation
from liulab_mbio.bench.steps import listed
from liulab_mbio.checks import Check
from liulab_mbio.cloning.restriction.bench import (
    BUFFER_FINDER,
    METHYLATION_FREE_HOST,
    RATIO_RANGE,
    shared_buffer,
)
from liulab_mbio.cloning.restriction.digest import Diagnostic
from liulab_mbio.cloning.restriction.ligation import Junction
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.sequence import Segment, SequenceRecord
from liulab_mbio.sites import DAM_SITE, DCM_SITE, CutSite, find_sites

#: The two methylases as something to search a record for. `liulab_mbio.sites.find_sites` is what
#: reads a motif through an IUPAC code and across a circular origin, and a methylase cuts
#: nothing. Keyed as `liulab_mbio.enzymes.Enzyme.methylation` keys its own answers.
METHYLASES: Mapping[str, Enzyme] = {
    "dam": Enzyme("Dam", DAM_SITE, top_cut=0, bottom_cut=0),
    "dcm": Enzyme("Dcm", DCM_SITE, top_cut=0, bottom_cut=0),
}

#: What a supplier says of an enzyme its methylation leaves alone.
INSENSITIVE = "not sensitive"

#: The word in a supplier's answer that makes it a property of the flanking bases rather than of
#: the enzyme. "Blocked by overlapping" is true only where a methylase site really does overlap,
#: so it is computed from the record and never raised from the enzyme's name. §11.
OVERLAPPING = "overlapping"

#: How many bases a codon is, which is the whole of the reading-frame arithmetic below.
CODON = 3


def buffer_check(enzymes: Sequence[Enzyme]) -> Check:
    """Whether these enzymes digest together in one buffer, or nothing sourced can say. §3, §4.

    Compared per product and not per enzyme name: BamHI and BamHI-HF are supplied in different
    buffers, and a supplier outside NEB's buffer system states none at all.
    """
    supplied = {enzyme.supplied_buffer for enzyme in enzymes}
    buffer = shared_buffer(enzymes)
    if buffer is not None:
        return Check(
            "buffer",
            "pass",
            len(supplied),
            f"every one is supplied in {buffer}, which is NEB's own rule for digesting two of "
            "them in one tube",
        )
    named = listed(
        [
            f"{enzyme.supplier_label} in {enzyme.supplied_buffer}"
            if enzyme.supplied_buffer
            else f"{enzyme.supplier_label}, no buffer stated"
            for enzyme in enzymes
        ]
    )
    return Check(
        "buffer",
        None,
        len(supplied),
        f"{named}. Judging a mixed pair takes the activity percentages, which no shippable "
        f"source states for these enzymes, so this carries no verdict rather than a pass: look "
        f"the pair up in {BUFFER_FINDER} before putting both in one tube",
    )


def temperature_check(enzymes: Sequence[Enzyme]) -> Check:
    """Whether the enzymes agree on an incubation temperature. §3."""
    wanted = sorted(
        {enzyme.incubation_celsius for enzyme in enzymes if enzyme.incubation_celsius is not None}
    )
    if not wanted:
        return Check(
            "digest temperature", None, 0, "no record states an incubation temperature for these"
        )
    if len(wanted) == 1:
        return Check("digest temperature", "pass", 1, f"every one incubates at {wanted[0]} °C")
    named = listed([f"{enzyme.name} at {enzyme.incubation_celsius} °C" for enzyme in enzymes])
    return Check(
        "digest temperature",
        "warn",
        len(wanted),
        f"{named}; NEB answers a pair wanting different temperatures with two digests, the "
        "first enzyme heat inactivated before the second goes in",
    )


def cleanup_check(enzymes: Sequence[Enzyme]) -> Check:
    """How each enzyme is stopped before the ligase goes in. §10.

    Heat where the supplier gives one, and the gel purification where it does not. Both are
    routes NEB names, so the verdict is a pass either way and the detail says which was taken.
    """
    unstoppable = [enzyme for enzyme in enzymes if heat_inactivation(enzyme) is None]
    if not unstoppable:
        named = listed(
            [
                f"{enzyme.supplier_label} at {enzyme.heat_inactivation_celsius} °C for "
                f"{enzyme.heat_inactivation_minutes} minutes"
                for enzyme in enzymes
            ]
        )
        return Check("digest clean-up", "pass", 0, f"heat stops every one: {named}")
    them = "them" if len(unstoppable) > 1 else "it"
    return Check(
        "digest clean-up",
        "pass",
        len(unstoppable),
        f"the supplier states no heat inactivation for "
        f"{listed([enzyme.supplier_label for enzyme in unstoppable])}, so the gel purification "
        f"is what takes {them} away; do not swap that step for a heat step",
    )


def methylation_check(enzymes: Sequence[Enzyme], records: Sequence[SequenceRecord]) -> Check:
    """Whether the host's Dam or Dcm methylation blocks a chosen enzyme on these plasmids. §11.

    A plasmid grown in an ordinary laboratory strain is Dam and Dcm methylated. An enzyme its
    supplier calls blocked or impaired is reported wherever it is used; one blocked only by an
    *overlapping* methylase site is reported only where the record really carries one, because
    that is a property of the flanking bases and not of the enzyme.
    """
    blocked = [
        said for enzyme in enzymes for record in records if (said := _blocked(enzyme, record))
    ]
    if not blocked:
        return Check(
            "methylation",
            "pass",
            0,
            "no chosen enzyme is blocked or impaired by Dam or Dcm methylation of these plasmids",
        )
    return Check(
        "methylation",
        "warn",
        len(blocked),
        f"{listed(blocked)}; prepare that DNA from {METHYLATION_FREE_HOST} before this digest, "
        "and transform the ligation itself into the ordinary host",
    )


def frame_check(product: SequenceRecord, junctions: Sequence[Junction]) -> Check:
    """Report the frame a junction's added bases leave, where one falls inside a coding sequence.

    A coding sequence the insert lands in comes through as one feature in two parts, and what
    stands between them is the insert and both junctions. A whole number of codons keeps the
    frame across it; anything else shifts what is read after it. The arithmetic is the genetic
    code's and no supplier's, so nothing here is cited.
    """
    interrupted = _interrupted(product, junctions)
    if not interrupted:
        return Check("reading frame", "pass", 0, "no junction falls inside a coding sequence")
    said = listed(
        [
            f"{name} gains {gap} bases between its two parts, "
            + (
                "a whole number of codons"
                if gap % CODON == 0
                else f"{gap % CODON} past a whole number of codons, so the frame shifts"
            )
            for name, gap in interrupted
        ]
    )
    shifted = sum(1 for _, gap in interrupted if gap % CODON)
    return Check("reading frame", "warn" if shifted else "pass", len(interrupted), said)


def ratio_check(backbone: Amount, insert: Amount) -> Check:
    """Judge the insert-to-vector molar ratio, in picomoles as well as nanograms. §7."""
    ratio = insert.pmol / backbone.pmol
    low, high = RATIO_RANGE
    said = (
        f"{ratio:g}:1 insert to vector -- {insert.pmol:g} pmol ({insert.nanograms:g} ng) of "
        f"{insert.name} against {backbone.pmol:g} pmol ({backbone.nanograms:g} ng) of "
        f"{backbone.name}"
    )
    return Check(
        "ligation ratio",
        "pass" if low <= ratio <= high else "warn",
        ratio,
        f"{said}; NEB calls a vector-to-insert ratio between 1:{low:g} and 1:{high:g} "
        "optimal for a single insertion",
    )


def diagnostic_check(diagnostic: Diagnostic) -> Check:
    """Whether the diagnostic digest tells a correct clone from the vector it went into."""
    clone, empty = diagnostic.names
    named = listed([enzyme.name for enzyme in diagnostic.enzymes])
    said = (
        f"{named} cut a correct {clone} to {_bp(diagnostic.clone)}; {empty} gives "
        f"{_bp(diagnostic.empty)}"
    )
    if diagnostic.tells_them_apart:
        return Check("diagnostic digest", "pass", len(diagnostic.clone), said)
    return Check(
        "diagnostic digest",
        "fail",
        len(diagnostic.clone),
        f"{said} -- the same bands, so this digest cannot tell one from the other",
    )


def _bp(bands: Sequence[int]) -> str:
    """Name a lane's bands, largest first as a gel reads them."""
    return listed([f"{band} bp" for band in bands])


def _interrupted(
    product: SequenceRecord, junctions: Sequence[Junction]
) -> tuple[tuple[str, int], ...]:
    """Each coding sequence a junction falls inside, and the bases standing between its parts."""
    found: list[tuple[str, int]] = []
    for feature in product.features:
        if feature.type != "CDS":
            continue
        for left, right in zip(feature.segments, feature.segments[1:], strict=False):
            gap = right.start - left.end
            if gap > 0 and any(left.end <= one.start < right.start for one in junctions):
                found.append((feature.name, gap))
    return tuple(found)


def _blocked(enzyme: Enzyme, record: SequenceRecord) -> str:
    """Say how this record's methylation blocks this enzyme, or nothing where it does not."""
    said: list[str] = []
    for key, methylase in METHYLASES.items():
        word = enzyme.methylation.get(key, INSENSITIVE)
        if word == INSENSITIVE:
            continue
        if OVERLAPPING not in word:
            said.append(
                f"{methylase.name} methylation of {record.name} leaves {enzyme.name} {word}"
            )
            continue
        at = [
            position
            for site in find_sites(record, enzyme)
            for position in _overlapping(record, site, methylase)
        ]
        if at:
            where = listed([str(position) for position in sorted(set(at))])
            said.append(
                f"{enzyme.name} is {word} and {record.name} carries a {methylase.name} site "
                f"over its own at {where}"
            )
    return listed(said)


def _overlapping(record: SequenceRecord, site: CutSite, methylase: Enzyme) -> tuple[int, ...]:
    """Where a methylase site overlaps this recognition site, in the record's own coordinates."""
    reach = len(methylase.site) - 1
    length = len(record)
    if record.topology == "circular":
        start, here = (site.start - reach) % length, reach
        width = min(reach + len(site.enzyme.site) + reach, length)
    else:
        start = max(site.start - reach, 0)
        here = site.start - start
        width = min(site.start + len(site.enzyme.site) + reach, length) - start
    window = record.extract(Segment(start, start + width))
    return tuple(
        (start + hit.start) % length
        for hit in find_sites(SequenceRecord(window), methylase)
        if hit.start < here + len(site.enzyme.site) and hit.start + len(methylase.site) > here
    )
