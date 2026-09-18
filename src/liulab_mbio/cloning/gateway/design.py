"""The attB tail a plain insert is amplified with, and the PCR that puts it there.

An insert carrying no att site gives the BP reaction nothing to recombine. The vendor's answer
is a PCR: each primer is an attB tail and an annealing region, and the amplicon is the DNA the
reaction takes.

The tail is the whole one. MAN0000470 page 13 asks for "Four guanine (G) residues at the
5' end", then "The 25-bp attB1 site", then "At least 18-25 bp of template- or gene-specific
sequences", and the troubleshooting row for few or no colonies is a tail short of that
(note §7, §16). A fusion adds the frame bases note §8 counts and nothing else.

The shortest figure the note reports is not a shorter tail: 11 bp is a *gene-specific overlap*
in a two-step adapter PCR at 2 pmol against 40, whose product still carries the whole 29-base
tail. Nothing here treats it as one.

Primers are designed and judged through `liulab_mbio.primers`, so the Tm is read from the
annealing region and the structures from the whole oligo. Both 5' ends are anchored at the
insert's own ends, because moving one would leave bases of the insert out of the product.
"""

from collections.abc import Mapping
from dataclasses import dataclass, replace
from itertools import groupby
from typing import Literal

from liulab_mbio.cloning.gateway.att import REGIONS
from liulab_mbio.edits import annealed, carried
from liulab_mbio.primers.design import design_pair
from liulab_mbio.primers.evaluation import PairReport, PrimerReport, evaluate_pair
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS, Thresholds
from liulab_mbio.sequence import Primer, SequenceRecord, Strand, reverse_complement

#: The guanines an attB primer carries 5' of the att site (note §7). The patent measured "four
#: or five" and both manuals say four. They never reach the clone: BP leaves them in the
#: by-product (note §5).
SPACER = "GGGG"

#: What a fusion adds to hold the reading frame: two bases after attB1 on the forward primer,
#: one before attB2 on the reverse (note §8). The note fixes how many and rules out only `AA`,
#: `AG` and `GA` for the pair, each of which makes a stop with attB1's final `T`; the bases
#: themselves are this package's choice inside that rule.
N_TERMINAL_FRAME = "CA"
C_TERMINAL_FRAME = "A"

#: Which tag the insert is read into, which is what decides the frame bases each tail carries.
type Fusion = Literal["none", "N-terminal", "C-terminal", "both"]

#: Every fusion there is, for a caller reading one off a name.
FUSIONS: tuple[Fusion, ...] = ("none", "N-terminal", "C-terminal", "both")

#: The frame bases each fusion adds to the forward tail and to the reverse one.
_FRAME: Mapping[Fusion, tuple[str, str]] = {
    "none": ("", ""),
    "N-terminal": (N_TERMINAL_FRAME, ""),
    "C-terminal": ("", C_TERMINAL_FRAME),
    "both": (N_TERMINAL_FRAME, C_TERMINAL_FRAME),
}


def attb_tails(fusion: Fusion = "none") -> tuple[str, str]:
    """Return the whole 5' tail of the forward primer and of the reverse one.

    Each is the spacer, the 25 bp att site, and the frame bases this fusion needs: 29 bases
    with no fusion, 31 on the forward primer for an N-terminal tag and 30 on the reverse for a
    C-terminal one.

    Raises
    ------
    KeyError
        If `fusion` is not one of `FUSIONS`.

    Examples
    --------
    >>> attb_tails()
    ('GGGGACAAGTTTGTACAAAAAAGCAGGCT', 'GGGGACCACTTTGTACAAGAAAGCTGGGT')
    >>> [len(tail) for tail in attb_tails("both")]
    [31, 30]
    """
    forward, reverse = _FRAME[fusion]
    return SPACER + REGIONS["attB1"] + forward, SPACER + REGIONS["attB2"] + reverse


#: The longest run of one base an attB tail spells, measured off the tails themselves: attB1
#: carries six A's. No design choice reaches it, so a band that fails a run that long fails the
#: vendor's own site rather than the oligo.
TAIL_RUN = max(len(list(run)) for tail in attb_tails("both") for _, run in groupby(tail))


def attb_thresholds(base: Thresholds = THRESHOLDS) -> Thresholds:
    """Return `base` with its run band widened to take the run an attB tail itself spells.

    Only the warning edge moves: the run is real and worth a warning, which the manual answers
    with an HPLC- or PAGE-purified oligo (note §16), and a run the annealing region adds beyond
    the tail's still fails. Every other band is the one the caller gave.

    Examples
    --------
    >>> attb_thresholds().mononucleotide_run.warn_high
    6
    """
    band = base.mononucleotide_run
    return replace(base, mononucleotide_run=replace(band, warn_high=max(band.warn_high, TAIL_RUN)))


@dataclass(frozen=True, slots=True)
class Amplicon:
    """The attB PCR that turns a plain insert into the DNA a BP reaction takes.

    Parameters
    ----------
    template
        The insert, as it was given.
    record
        What the PCR makes: both tails, the insert, its features carried to their new
        coordinates and each primer annotated where it anneals.
    tails
        The forward primer's 5' tail and the reverse primer's, as each is written.
    report
        What the pair scored on the insert, carrying the annealing temperature and the
        extension time the PCR needs.
    polymerase, thresholds
        What it was designed and judged with.
    """

    template: SequenceRecord
    record: SequenceRecord
    tails: tuple[str, str]
    report: PairReport
    polymerase: Polymerase
    thresholds: Thresholds

    @property
    def name(self) -> str:
        """What the amplicon and its two primers are called."""
        return self.record.name

    @property
    def length(self) -> int:
        """Bases of amplicon, which is what a gel measures and a PCR program times."""
        return len(self.record)

    @property
    def primers(self) -> tuple[Primer, ...]:
        """The two primers, forward first."""
        return tuple(report.primer for report in self.reports)

    @property
    def reports(self) -> tuple[PrimerReport, ...]:
        """What each primer scored on its own, forward first."""
        return (self.report.forward, self.report.reverse)


def amplify_attb(
    insert: SequenceRecord,
    *,
    fusion: Fusion = "none",
    name: str = "",
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Amplicon:
    """Design the PCR that puts an attB site on each end of `insert`.

    The whole record is amplified, so both 5' ends are anchored and only the annealing regions'
    lengths vary; each primer is then judged whole, structures and off-target sites included,
    and warns only where no length allowed does better.

    Parameters
    ----------
    insert
        The record to amplify. It needs no att site, and gains one at each end.
    fusion
        Which tag it is read into, which decides the frame bases each tail carries.
    name
        What to call the amplicon and its primers.
    polymerase, thresholds
        Passed to `liulab_mbio.primers`, the thresholds through `attb_thresholds`.

    Returns
    -------
    Amplicon
        The primers, what they scored, and the attB DNA they make.

    Raises
    ------
    KeyError
        If `fusion` is not one of `FUSIONS`.
    ValueError
        If no annealing region fits either end of the insert.
    """
    forward_tail, reverse_tail = attb_tails(fusion)
    judged = attb_thresholds(thresholds)
    label = name or (f"attB-{insert.name}" if insert.name else "attB product")
    try:
        forward, reverse = design_pair(
            insert,
            0,
            len(insert),
            forward_tail=forward_tail,
            reverse_tail=reverse_tail,
            forward_name=f"{label} forward",
            reverse_name=f"{label} reverse",
            polymerase=polymerase,
            thresholds=judged,
        )
    except ValueError as reason:
        raise ValueError(
            f"{insert.name or 'this insert'} takes no attB primer: {reason}. Each end needs an "
            f"annealing region of {judged.length.warn_low:g}-{judged.length.warn_high:g} bases "
            f"under its tail, and this record is {len(insert)} bases long"
        ) from reason
    bases = forward_tail + insert.sequence + reverse_complement(reverse_tail)
    features, kept = carried(insert, 0, len(insert), offset=len(forward_tail))
    placed = (
        annealed(forward, len(forward_tail), Strand.FORWARD),
        annealed(reverse, len(bases) - len(reverse_tail), Strand.REVERSE),
    )
    return Amplicon(
        insert,
        SequenceRecord(bases, name=label, features=features, primers=kept + placed),
        (forward_tail, reverse_tail),
        evaluate_pair(forward, reverse, insert, polymerase=polymerase, thresholds=judged),
        polymerase,
        judged,
    )
