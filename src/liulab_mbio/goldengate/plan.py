"""One Golden Gate experiment, planned from a vector and the inserts that go round it.

`plan_assembly` joins as many inserts as the overhangs allow, given in the order they go round
the product, and runs the whole design: choose the enzyme, design the overhangs, simulate the
PCRs and the ligation, work out the bench quantities, and design the colony PCR and sequencing
that validate the clone. `Plan.write` puts the three things a bench needs in one directory --
the annotated product, a primer order sheet, and the interactive HTML protocol.

Every number the protocol prints is computed here or by the modules this one calls. What the
protocol says about the phenotype -- what drives the inserts, whether anything should be
translated, and how a plate reads -- is `Phenotype`, read off the product's own features.
"""

import dataclasses
import os
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from liulab_mbio.bench import (
    COLONY_FLANK,
    Amount,
    ColonyCheck,
    SangerRead,
    colony_pcr_check,
    sanger_primers,
)
from liulab_mbio.checks import Check, Status, worst
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.goldengate.assembly import Assembly, Part, amplify, assemble, open_vector
from liulab_mbio.goldengate.bench import assembly_amounts
from liulab_mbio.goldengate.design import (
    EnzymeChoice,
    Junction,
    OverhangSet,
    choose_enzyme,
    design_overhangs,
)
from liulab_mbio.goldengate.ligase import LigaseProfile, read_profile
from liulab_mbio.goldengate.steps import DEFAULT_HOST
from liulab_mbio.goldengate.steps import protocol as protocol_for
from liulab_mbio.io import read_record
from liulab_mbio.primers import (
    ONETAQ,
    Q5,
    THRESHOLDS,
    Polymerase,
    PrimerReport,
    Thresholds,
    design_pair,
    evaluate_primer,
    reading,
)
from liulab_mbio.protocol import Protocol, write_html
from liulab_mbio.sequence import (
    BindingSite,
    Feature,
    Primer,
    Segment,
    SequenceRecord,
    Strand,
    reverse_complement,
)
from liulab_mbio.sites import EnzymeLike
from liulab_mbio.snapgene import write_dna

#: Which way round an insert goes into the vector.
type Orientation = Literal["forward", "reverse"]

#: Where the inserts go: a feature name, a span, or `None` to look for `MCS_FEATURE`.
type Site = str | tuple[int, int] | None

#: The feature a vector names its cloning site with, looked for when the caller names none.
MCS_FEATURE = "MCS"

#: How far the vector junction may slide to get past an overhang rule. It moves where the vector
#: is cut inside the span the assembly replaces, so the product keeps a base or two more of it.
VECTOR_WINDOW = 6

#: Bases of vector between the downstream colony PCR primer and its own junction. Deliberately
#: not `COLONY_FLANK`: two primers the same distance from their junctions give a reversed insert
#: the same bands as a correct one, so the gel could not tell them apart.
REVERSE_FLANK = 2 * COLONY_FLANK

#: Selection markers this package can name an antibiotic for, keyed by the feature name lowered.
#: pUC19's is the one `docs/research/golden-gate-assembly.md` §8 states a plate recipe for; a
#: marker absent from here is named rather than translated.
SELECTION: Mapping[str, str] = {
    "ampr": "ampicillin or carbenicillin",
    "bla": "ampicillin or carbenicillin",
}

#: What `Plan.write` calls the three files it writes.
PRODUCT_FILE = "product.dna"
PRIMER_FILE = "primers.tsv"
PROTOCOL_FILE = "protocol.html"

#: The columns of the primer order sheet.
SHEET_COLUMNS = ("name", "sequence", "length", "tm_c")


@dataclass(frozen=True, slots=True)
class Files:
    """The three files a plan writes.

    Parameters
    ----------
    product
        The annotated product, as a SnapGene ``.dna`` file.
    primers
        Every designed oligo, as a tab-separated sheet to order from.
    protocol
        The interactive bench protocol, as one self-contained HTML page.
    """

    product: Path
    primers: Path
    protocol: Path


@dataclass(frozen=True, slots=True)
class Phenotype:
    """What the product says about itself, read off its own features.

    Parameters
    ----------
    insert
        The span the inserts occupy in the product, between the first junction and the last.
    coding
        The longest coding sequence in that span, or ``None`` when it annotates none.
    promoter
        The promoter nearest the insert on the promoter's own reading direction, or ``None``.
    gap_bp
        Bases between that promoter and the insert.
    driven
        Whether that promoter reads along the strand the insert is coded on.
    ribosome_binding_site
        Whether one is annotated between that promoter and the insert.
    reporter
        The vector coding sequence the insertion interrupts, or ``None``.
    marker
        The vector's selection marker, or ``None`` when it annotates none this package knows.
    """

    insert: tuple[int, int]
    coding: Feature | None
    promoter: Feature | None
    gap_bp: int
    driven: bool
    ribosome_binding_site: bool
    reporter: Feature | None
    marker: Feature | None

    @property
    def expressed(self) -> bool:
        """Whether the product should make the insert's protein."""
        return self.coding is not None and self.driven and self.ribosome_binding_site

    @property
    def blue_white(self) -> bool:
        """Whether X-gal and IPTG tell a correct clone from an empty vector.

        True when the insertion interrupts a lacZ fragment, which is then not there to
        complement the host's own.
        """
        return self.reporter is not None and self.reporter.name.lower().startswith("lacz")

    @property
    def antibiotic(self) -> str:
        """What to select transformants on, or an empty string when the marker is unknown."""
        if self.marker is None:
            return ""
        return SELECTION.get(self.marker.name.lower(), "")


#: What an oligo is for: amplifying a part, colony PCR, or sequencing the clone.
type OligoRole = Literal["amplification", "colony PCR", "sequencing"]


@dataclass(frozen=True, slots=True)
class DesignedOligo:
    """One oligo a plan orders, and what it is for.

    Parameters
    ----------
    report
        What it scored, the primer included.
    role
        What it is for.
    part
        The part it amplifies, given for an amplification primer and for nothing else.

    Raises
    ------
    ValueError
        If a part is given for any role but amplification, or not given for amplification.
    """

    report: PrimerReport
    role: OligoRole
    part: Part | None = None

    def __post_init__(self) -> None:
        """Refuse a part on anything but an amplification primer, and one missing from it."""
        if (self.role == "amplification") != (self.part is not None):
            raise ValueError(
                f"oligo {self.report.primer.name!r}: an amplification primer names its part, "
                "and no other role does"
            )


@dataclass(frozen=True, slots=True)
class Plan:
    """One planned Golden Gate experiment.

    Parameters
    ----------
    vector, inserts
        The records the plan was made from, the inserts in the order they go round the product
        and each on the strand that goes in.
    span
        The vector bases the assembly replaces, after any junction slide.
    choice
        The enzyme the plan uses, and what using it costs.
    ranking
        Every candidate enzyme, best first.
    overhangs
        The overhang every junction takes, with the fidelity of the whole set. An assembly of n
        inserts has n + 1 junctions.
    linearised_vector
        The part the vector is opened into.
    insert_parts
        The part each insert is amplified into, in insert order.
    assembly
        The simulated product, the parts and the junctions.
    colony
        The colony PCR that reads every junction and tells a correct clone from an empty vector
        or from one carrying an insert the other way round.
    reads
        A sequencing primer reading in from outside the first junction and the last.
    amounts
        What to put in the assembly reaction, vector first.
    phenotype
        What the product says about itself.
    designed_oligos
        Every designed oligo and what it is for, in order: each part's PCR, the colony PCR, the
        reads.
    host, polymerase
        The choices the protocol names.
    thresholds
        What those oligos were judged by, so a page prints the band beside the value.
    """

    vector: SequenceRecord
    inserts: tuple[SequenceRecord, ...]
    span: tuple[int, int]
    choice: EnzymeChoice
    ranking: tuple[EnzymeChoice, ...]
    overhangs: OverhangSet
    linearised_vector: Part
    insert_parts: tuple[Part, ...]
    assembly: Assembly
    colony: ColonyCheck
    reads: tuple[SangerRead, SangerRead]
    amounts: tuple[Amount, ...]
    phenotype: Phenotype
    designed_oligos: tuple[DesignedOligo, ...]
    host: str
    polymerase: Polymerase
    thresholds: Thresholds = THRESHOLDS

    @property
    def enzyme(self) -> Enzyme:
        """The Type IIS enzyme the assembly is cut with."""
        return self.assembly.enzyme

    @property
    def product(self) -> SequenceRecord:
        """The circular plasmid the assembly makes."""
        return self.assembly.product

    @property
    def parts(self) -> tuple[Part, ...]:
        """The parts that go into the reaction, the linearised vector first."""
        return (self.linearised_vector, *self.insert_parts)

    @property
    def oligos(self) -> tuple[Primer, ...]:
        """Every oligo the plan designs, in the order the sheet lists them."""
        return tuple(report.primer for report in self.reports)

    @property
    def reports(self) -> tuple[PrimerReport, ...]:
        """Every designed oligo's evaluation, in the order the sheet lists them."""
        return tuple(oligo.report for oligo in self.designed_oligos)

    @property
    def checks(self) -> tuple[Check, ...]:
        """The product's checks, with one more for the oligos."""
        return (
            *self.assembly.checks,
            Check(
                "primers",
                worst(report.status for report in self.reports),
                len(self.reports),
                _primer_detail(self.reports, self.thresholds),
            ),
        )

    @property
    def status(self) -> Status:
        """The worst status of any check."""
        return worst(check.status for check in self.checks)

    def protocol(self) -> Protocol:
        """Return the bench protocol for this plan."""
        return protocol_for(self)

    def write(self, directory: str | os.PathLike[str]) -> Files:
        """Write the product, the primer sheet and the protocol into `directory`.

        The directory is made when it is not there. The three files are named by
        `PRODUCT_FILE`, `PRIMER_FILE` and `PROTOCOL_FILE`, and a second run over the same
        inputs writes the same bytes.
        """
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        product = out / PRODUCT_FILE
        write_dna(self.product, product)
        sheet = out / PRIMER_FILE
        sheet.write_text(primer_sheet(self), encoding="utf-8")
        return Files(product, sheet, write_html(self.protocol(), out / PROTOCOL_FILE))


def plan_assembly(
    vector: SequenceRecord | str | os.PathLike[str],
    *inserts: SequenceRecord | str | os.PathLike[str],
    site: Site = None,
    orientation: Orientation | Sequence[Orientation] = "forward",
    in_frame: bool | Sequence[bool] = False,
    enzyme: EnzymeLike | None = None,
    profile: LigaseProfile | str | os.PathLike[str] | None = None,
    prefer_profile: bool = False,
    polymerase: Polymerase = Q5,
    host: str = DEFAULT_HOST,
    name: str = "",
    window: int = VECTOR_WINDOW,
    thresholds: Thresholds = THRESHOLDS,
) -> Plan:
    """Plan one Golden Gate experiment putting `inserts` into `vector`.

    One reaction joins as many inserts as the overhangs allow, given in the order they go round
    the product. The vector is opened by PCR across the span they replace, each insert is
    amplified with tails of its own, and every junction is scarless: it takes the bases the part
    already spells there. Only the vector junction may slide, by up to `window` bases, to get
    past an overhang rule, which moves where the vector is cut and not what the inserts spell.

    Parameters
    ----------
    vector, inserts
        A record, or a path to a ``.dna``, GenBank or FASTA file holding one.
    site
        Where the inserts go: a feature name, a ``(start, end)`` span of the vector, or
        ``None`` to use the vector's own `MCS_FEATURE` feature.
    orientation
        ``"reverse"`` puts the other strand of an insert into the product. One value covers
        every insert; a sequence gives one for each.
    in_frame
        Hold an insert's junction on a codon boundary of the coding sequence it lies in. One
        value covers every insert; a sequence gives one for each.
    enzyme
        The Type IIS enzyme to use. Chosen by `choose_enzyme` when not given, and refused
        either way if it reads a site in any part.
    profile
        A ligase fidelity matrix the caller holds, as a path or an already read `LigaseProfile`.
        It scores the overhangs where no shipped matrix covers the enzyme. The package ships
        none: see `liulab_mbio.goldengate.ligase`.
    prefer_profile
        Use it even where a shipped matrix covers the enzyme.
    polymerase
        For the PCRs. The colony PCR uses OneTaq, which is what NEB's protocol asks for.
    host, name
        The strain the protocol names, and what to call the product.
    window
        How far the vector junction may slide.
    thresholds
        Passed to `liulab_mbio.primers`.

    Returns
    -------
    Plan
        The design, the simulated product and the validation.

    Raises
    ------
    ValueError
        If no insert is given, if no insertion site is named and the vector annotates none, if
        the enzyme reads a site in a part, if no overhang passes every rule, if `profile` names
        a file that is not a count matrix, or if the parts do not assemble.
    """
    if not inserts:
        raise ValueError("an assembly needs a vector and at least one insert")
    one = _record(vector)
    ways = _orientations(orientation, len(inserts))
    frames = _frames(in_frame, len(inserts))
    going = [
        flipped(read) if way == "reverse" else read
        for read, way in ((_record(record), way) for record, way in zip(inserts, ways, strict=True))
    ]
    labels = [record.name or f"insert {number}" for number, record in enumerate(going, start=1)]
    start, end = _span(one, site)
    ranking = choose_enzyme([one, *going])
    choice = _chosen(ranking, enzyme, (one, *going))
    chosen = choice.enzyme
    designed = design_overhangs(
        (
            *(
                Junction(label, record=record, position=0, scarless=not frame, in_frame=frame)
                for label, record, frame in zip(labels, going, frames, strict=True)
            ),
            Junction(one.name or "vector", record=one, position=end, scarless=True, window=window),
        ),
        chosen,
        profile=_profile(profile),
        prefer_profile=prefer_profile,
    )
    overhangs = designed.overhangs
    span = (start, end + designed.choices[-1].offset)
    linearised_vector = open_vector(
        one,
        chosen,
        *span,
        overhangs=(overhangs[0], overhangs[-1]),
        name=f"{one.name} backbone".strip(),
        polymerase=polymerase,
        thresholds=thresholds,
    )
    insert_parts = tuple(
        amplify(
            record,
            chosen,
            0,
            len(record),
            left_overhang=overhangs[number],
            right_overhang=overhangs[number + 1],
            name=label,
            polymerase=polymerase,
            thresholds=thresholds,
        )
        for number, (label, record) in enumerate(zip(labels, going, strict=True))
    )
    parts = (linearised_vector, *insert_parts)
    built = assemble(parts, chosen, name=name or "-".join([one.name, *labels]).strip("-"))
    junctions = built.junction_positions
    first, last = junctions[0], junctions[-1]
    colony = colony_pcr_check(
        built.product,
        junctions,
        vector=one,
        primers=design_pair(
            built.product,
            first - COLONY_FLANK,
            last + REVERSE_FLANK,
            forward_name="Colony PCR forward",
            reverse_name="Colony PCR reverse",
            polymerase=ONETAQ,
            thresholds=thresholds,
        ),
        insert_primer=True,
        polymerase=ONETAQ,
        thresholds=thresholds,
    )
    reads = sanger_primers(built.product, junctions, thresholds=thresholds)
    return Plan(
        one,
        tuple(going),
        span,
        choice,
        ranking,
        designed,
        linearised_vector,
        insert_parts,
        built,
        colony,
        reads,
        assembly_amounts(
            (linearised_vector.name, linearised_vector.length),
            tuple((part.name, part.length) for part in insert_parts),
        ),
        _phenotype(one, built, span),
        (
            *(
                DesignedOligo(report, "amplification", part)
                for part in parts
                for report in (part.report.forward, part.report.reverse)
            ),
            *(DesignedOligo(report, "colony PCR") for report in colony.reports),
            *(
                DesignedOligo(
                    evaluate_primer(read.primer, built.product, thresholds=thresholds),
                    "sequencing",
                )
                for read in reads
            ),
        ),
        host,
        polymerase,
        thresholds,
    )


def _primer_detail(reports: tuple[PrimerReport, ...], thresholds: Thresholds) -> str:
    """Return what the oligos' verdicts say: the counts, the kinds, and what nothing judged.

    A count alone cannot be acted on, so each kind that fired is named with the rows it covers.
    """
    warned = sum(1 for report in reports if report.status == "warn")
    failed = sum(1 for report in reports if report.status == "fail")
    counted = Counter(
        reading(check, thresholds).label
        for report in reports
        for check in report.checks
        if check.status not in (None, "pass")
    )
    unjudged = dict.fromkeys(
        reading(check, thresholds).label
        for report in reports
        for check in report.checks
        if check.status is None
    )
    said = f"{len(reports)} designed, {warned} with a warning, {failed} failing"
    if counted:
        kinds = ", ".join(
            f"{label} on {rows}"
            for label, rows in sorted(counted.items(), key=lambda one: (-one[1], one[0]))
        )
        said += f": {kinds}"
    return f"{said}; not judged: {', '.join(unjudged)}" if unjudged else said


def _orientations(
    orientation: Orientation | Sequence[Orientation], count: int
) -> tuple[Orientation, ...]:
    """Spread one orientation over every insert, or take the one given for each.

    Raises
    ------
    ValueError
        If a value is neither ``forward`` nor ``reverse``, or there is not one per insert.
    """
    given = (orientation,) * count if isinstance(orientation, str) else tuple(orientation)
    if len(given) != count:
        raise ValueError(f"orientation has {len(given)} values for {count} insert(s)")
    for one in given:
        if one not in ("forward", "reverse"):
            raise ValueError(f"orientation is 'forward' or 'reverse', got {one!r}")
    return given


def _frames(in_frame: bool | Sequence[bool], count: int) -> tuple[bool, ...]:
    """Spread one in-frame choice over every insert, or take the one given for each.

    Raises
    ------
    ValueError
        If there is not one per insert.
    """
    given = (in_frame,) * count if isinstance(in_frame, bool) else tuple(in_frame)
    if len(given) != count:
        raise ValueError(f"in_frame has {len(given)} values for {count} insert(s)")
    return given


def primer_sheet(plan: Plan) -> str:
    """Return every designed oligo as a tab-separated sheet, one row each.

    The columns are `SHEET_COLUMNS`: the name to order it under, the sequence 5' to 3', its
    length, and the Tm of the part that anneals.
    """
    rows = ["\t".join(SHEET_COLUMNS)]
    for report in plan.reports:
        primer = report.primer
        rows.append(
            "\t".join(
                (
                    primer.name,
                    primer.sequence,
                    str(len(primer.sequence)),
                    f"{report['tm'].value:.1f}",
                )
            )
        )
    return "\n".join(rows) + "\n"


def flipped(record: SequenceRecord) -> SequenceRecord:
    """Return `record` read from the other strand, features and binding sites turned with it.

    Raises
    ------
    ValueError
        If a span runs across the origin, which has no place on the other strand of a record
        this turns end for end.

    Examples
    --------
    >>> flipped(SequenceRecord("AAAACCCG")).sequence
    'CGGGTTTT'
    """
    length = len(record)
    other = {Strand.FORWARD: Strand.REVERSE, Strand.REVERSE: Strand.FORWARD}
    spans = [
        (segment.start, segment.end) for feature in record.features for segment in feature.segments
    ]
    spans += [(site.start, site.end) for primer in record.primers for site in primer.binding_sites]
    if any(end > length for _, end in spans):
        raise ValueError("a record with a span across its origin cannot be turned end for end")
    features = tuple(
        dataclasses.replace(
            feature,
            segments=tuple(
                sorted(
                    (
                        Segment(
                            length - segment.end,
                            length - segment.start,
                            name=segment.name,
                            color=segment.color,
                        )
                        for segment in feature.segments
                    ),
                    key=lambda segment: (segment.start, segment.end),
                )
            ),
            strand=other.get(feature.strand, feature.strand),
        )
        for feature in record.features
    )
    primers = tuple(
        dataclasses.replace(
            primer,
            binding_sites=tuple(
                BindingSite(length - site.end, length - site.start, other[site.strand])
                for site in primer.binding_sites
            ),
        )
        for primer in record.primers
    )
    return dataclasses.replace(
        record,
        sequence=reverse_complement(record.sequence),
        features=features,
        primers=primers,
        extras={},
    )


def _record(value: SequenceRecord | str | os.PathLike[str]) -> SequenceRecord:
    """Read a record, or take one already read."""
    return value if isinstance(value, SequenceRecord) else read_record(value)


def _profile(value: LigaseProfile | str | os.PathLike[str] | None) -> LigaseProfile | None:
    """Read a ligase profile, or take one already read."""
    if value is None or isinstance(value, LigaseProfile):
        return value
    return read_profile(value)


def _span(vector: SequenceRecord, site: Site) -> tuple[int, int]:
    """Return the vector bases the assembly replaces."""
    if isinstance(site, tuple):
        start, end = site
    else:
        if site is None:
            found = _named(vector, MCS_FEATURE)
            if found is None:
                raise ValueError(
                    f"name the insertion site: this vector annotates no {MCS_FEATURE!r} "
                    "feature. Pass a feature name or a (start, end) span"
                )
        else:
            found = _named(vector, site)
            if found is None:
                raise ValueError(f"this vector annotates no feature called {site!r}")
        start, end = found.segments[0].start, found.segments[-1].end
    if not 0 <= start < end <= len(vector):
        raise ValueError(
            f"the insertion site {start}-{end} does not lie inside {len(vector)} bases"
        )
    return start, end


def _named(record: SequenceRecord, name: str) -> Feature | None:
    """Return the first feature of that name, whatever its case."""
    return next(
        (feature for feature in record.features if feature.name.lower() == name.lower()), None
    )


def _chosen(
    ranking: tuple[EnzymeChoice, ...],
    enzyme: EnzymeLike | None,
    parts: tuple[SequenceRecord, ...],
) -> EnzymeChoice:
    """Return the enzyme to use, refusing one that would cut the product open.

    Raises
    ------
    ValueError
        If the enzyme reads a site in either part.
    """
    if enzyme is None:
        choice = ranking[0]
    else:
        wanted = get_enzyme(enzyme) if isinstance(enzyme, str) else enzyme
        choice = next(
            (one for one in ranking if one.enzyme.name == wanted.name),
            choose_enzyme(parts, enzymes=[wanted])[0],
        )
    if choice.free:
        return choice
    free = [one.enzyme.name for one in ranking if one.free]
    if free:
        rest = f"free here: {', '.join(free)}"
    else:
        rest = (
            f"no candidate is free, and of this one's sites {len(choice.changes)} could go by a "
            f"synonymous codon change and {len(choice.outside_cds)} lie outside a coding sequence"
        )
    raise ValueError(
        f"{choice.enzyme.name} reads {choice.sites} site(s) in the parts, so it would cut the "
        f"product open again; {rest}"
    )


def _phenotype(vector: SequenceRecord, built: Assembly, span: tuple[int, int]) -> Phenotype:
    """Read what the product says about itself off its own features."""
    junctions = built.junction_positions
    first, last = junctions[0], junctions[-1]
    coding = _coding(built.product, first, last)
    promoter, gap = _promoter(built.product, first, last)
    return Phenotype(
        (first, last),
        coding,
        promoter,
        gap,
        promoter is not None and coding is not None and promoter.strand == coding.strand,
        _ribosome_binding_site(built.product, promoter, first, last),
        _interrupted(vector, span),
        _marker(vector),
    )


def _coding(product: SequenceRecord, first: int, last: int) -> Feature | None:
    """Return the longest coding sequence lying wholly between the outer two junctions."""
    inside = [
        feature
        for feature in product.features
        if feature.type == "CDS"
        and all(first <= segment.start and segment.end <= last for segment in feature.segments)
    ]
    return max(
        inside,
        key=lambda feature: sum(segment.end - segment.start for segment in feature.segments),
        default=None,
    )


def _promoter(product: SequenceRecord, first: int, last: int) -> tuple[Feature | None, int]:
    """Return the promoter nearest the insert along its own reading direction, and the gap."""
    length = len(product)
    found: Feature | None = None
    gap = length
    for feature in product.features:
        if feature.type != "promoter":
            continue
        low = min(segment.start for segment in feature.segments)
        high = max(segment.end for segment in feature.segments)
        distance = (
            (low - last) % length if feature.strand == Strand.REVERSE else (first - high) % length
        )
        if distance < gap:
            found, gap = feature, distance
    return found, gap if found is not None else 0


def _ribosome_binding_site(
    product: SequenceRecord, promoter: Feature | None, first: int, last: int
) -> bool:
    """Whether one is annotated between the promoter and the insert."""
    if promoter is None:
        return False
    length = len(product)
    if promoter.strand == Strand.REVERSE:
        low, high = last, min(segment.start for segment in promoter.segments)
    else:
        low, high = max(segment.end for segment in promoter.segments), first
    return any(
        feature.type == "RBS"
        and any(
            (segment.start - low) % length < (high - low) % length for segment in feature.segments
        )
        for feature in product.features
    )


def _interrupted(vector: SequenceRecord, span: tuple[int, int]) -> Feature | None:
    """Return the vector coding sequence the insertion breaks, or ``None``."""
    start, end = span
    return next(
        (
            feature
            for feature in vector.features
            if feature.type == "CDS"
            and any(segment.start < end and start < segment.end for segment in feature.segments)
        ),
        None,
    )


def _marker(vector: SequenceRecord) -> Feature | None:
    """Return the vector's selection marker, or ``None`` when it annotates none."""
    return next(
        (
            feature
            for feature in vector.features
            if feature.type == "CDS" and feature.name.lower() in SELECTION
        ),
        None,
    )
