"""Amplify an insert with a recognition site on each 5' tail, for an insert no plasmid holds.

Where the insert does not already sit between two usable sites, the sites are put on it: each
primer carries a spacer and the recognition site the end it makes has to be cut at, which
`liulab_mbio.sites.primer_tail` builds. The amplicon is then digested like any other record, so
the piece it releases reaches the ligation with the same two ends an excised insert has.

Two things the tails are held to:

- **The site has to be far enough from the end to be cut at all.** NEB measures per-enzyme
  cleavage at one to five spacer bases and answers six for the enzymes it does not list
  (``docs/research/restriction-ligation.md`` §2), which is what `liulab_mbio.sites.SPACER_LENGTH`
  already carries.
- **A tail spells one site and no other.** `primer_tail` moves the spacer off a second site of
  either enzyme; a site the tail completes against the insert's own first bases cannot be moved
  at all, and `liulab_mbio.cloning.restriction.digest.cut` refuses the amplicon, naming it.

Coordinates are the model's, 0-based and half-open.
"""

from dataclasses import dataclass

from liulab_mbio.bench.steps import dam_sites
from liulab_mbio.cloning.restriction.digest import Piece
from liulab_mbio.edits import annealed, carried
from liulab_mbio.enzymes import Enzyme
from liulab_mbio.overhangs import End
from liulab_mbio.primers.design import design_pair
from liulab_mbio.primers.evaluation import PairReport, evaluate_pair
from liulab_mbio.primers.polymerase import Q5, Polymerase
from liulab_mbio.primers.thresholds import THRESHOLDS, Thresholds
from liulab_mbio.sequence import Primer, SequenceRecord, Strand, reverse_complement
from liulab_mbio.sites import SPACER_LENGTH, primer_tail


@dataclass(frozen=True, slots=True)
class Amplicon:
    """One insert as a PCR product, with the site each end is cut at on its own tail.

    Parameters
    ----------
    name
        What the tube and its gel lane are labelled.
    template
        The record the insert is amplified from.
    record
        What the PCR makes: both tails, the insert, and the template's features and primers
        carried to their new coordinates.
    left_tail, right_tail
        The bases each tail adds, written on the top strand. The reverse primer carries
        `right_tail` reverse-complemented.
    left_enzyme, right_enzyme
        The enzyme each tail's site is for, in the same order the tails are written.
    report
        What the pair scored on the template, carrying the annealing temperature and the
        extension time the PCR needs. A tail moves neither annealing region, so each melting
        temperature is of the annealing region and the structure checks are of the whole oligo.
    """

    name: str
    template: SequenceRecord
    record: SequenceRecord
    left_tail: str
    right_tail: str
    left_enzyme: Enzyme
    right_enzyme: Enzyme
    report: PairReport

    @property
    def length(self) -> int:
        """Bases of amplicon, which is what a gel measures and a PCR program times."""
        return len(self.record)

    @property
    def primers(self) -> tuple[Primer, Primer]:
        """The two oligos to order, the forward one first."""
        return self.report.forward.primer, self.report.reverse.primer

    @property
    def dpni(self) -> bool:
        """Whether DpnI should take the template away afterwards.

        A circular template is a plasmid from a Dam-positive host, and `dam_sites` says whether
        DpnI can cut it. A linear fragment was never methylated.
        """
        return self.template.topology == "circular" and dam_sites(self.template) > 0

    @property
    def ends(self) -> tuple[tuple[str, str, Enzyme], ...]:
        """Each end: the primer that carries its tail, the tail, and the enzyme the site is for."""
        return (
            ("forward", self.left_tail, self.left_enzyme),
            ("reverse", self.right_tail, self.right_enzyme),
        )


def amplified(
    insert: SequenceRecord,
    *,
    into: Piece,
    name: str = "",
    spacer_length: int = SPACER_LENGTH,
    polymerase: Polymerase = Q5,
    thresholds: Thresholds = THRESHOLDS,
) -> Amplicon:
    """Design the PCR that puts the sites `into` needs on each end of `insert`.

    The backbone's right end meets the insert's left, so the forward primer carries the enzyme
    that cut the backbone's right end and the reverse primer the one that cut its left. Each
    tail leaves, when cut, exactly the overhang that end has to anneal to.

    Parameters
    ----------
    insert
        The record to amplify. It goes into the product whole.
    into
        The backbone the insert is going into, which names the enzyme for each end.
    name
        Names the amplicon and its two primers.
    spacer_length
        How many bases to put 5' of each recognition site, so the site is cut at all.
    polymerase, thresholds
        Passed to `liulab_mbio.primers`.

    Raises
    ------
    ValueError
        If `insert` is empty, if no spacer leaves a tail spelling one site, or if no annealing
        region fits either end.
    """
    if not len(insert):
        raise ValueError("an insert needs bases to amplify")
    left_enzyme, right_enzyme = into.right_enzyme, into.left_enzyme
    left = primer_tail(
        left_enzyme,
        _wanted(left_enzyme, into.right_end),
        spacer_length=spacer_length,
        avoid=_besides(left_enzyme, right_enzyme),
    )
    reverse_tail = primer_tail(
        right_enzyme,
        _wanted(right_enzyme, into.left_end, flip=True),
        spacer_length=spacer_length,
        avoid=_besides(right_enzyme, left_enzyme),
    )
    right = reverse_complement(reverse_tail)
    called = name or insert.name or "insert"
    forward_primer, reverse_primer = design_pair(
        insert,
        0,
        len(insert),
        forward_tail=left,
        reverse_tail=reverse_tail,
        forward_name=f"{called} forward",
        reverse_name=f"{called} reverse",
        polymerase=polymerase,
        thresholds=thresholds,
    )
    bases = left + insert.sequence + right
    features, kept = carried(insert, 0, len(insert), offset=len(left))
    placed = (
        annealed(forward_primer, len(left), Strand.FORWARD),
        annealed(reverse_primer, len(bases) - len(right), Strand.REVERSE),
    )
    return Amplicon(
        called,
        insert,
        SequenceRecord(
            bases, name=f"{called} amplicon", features=features, primers=(*kept, *placed)
        ),
        left,
        right,
        left_enzyme,
        right_enzyme,
        evaluate_pair(
            forward_primer, reverse_primer, insert, polymerase=polymerase, thresholds=thresholds
        ),
    )


def _besides(enzyme: Enzyme, other: Enzyme) -> tuple[Enzyme, ...]:
    """Return the other end's enzyme, and nothing where it is the same one.

    A tail spells its own enzyme's site by design, so that enzyme is never its own avoid.
    """
    return () if other == enzyme else (other,)


def _wanted(enzyme: Enzyme, end: End, *, flip: bool = False) -> str:
    """Return the overhang a tail is built to leave, for the end it has to anneal to.

    An enzyme cutting inside its own site leaves whatever that site spells, so it is not one to
    choose. `flip` is for the reverse primer, whose tail is written on the other strand.
    """
    if enzyme.type == "II":
        return ""
    return reverse_complement(end.overhang) if flip else end.overhang
