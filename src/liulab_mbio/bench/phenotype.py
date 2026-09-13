"""What a clone is expected to show, read off the product's own features.

What drives the inserts, whether anything should be translated, and how a plate reads: a
protocol states these from the features rather than a person asserting them.
"""

from collections.abc import Mapping
from dataclasses import dataclass

from liulab_mbio.sequence import Feature, SequenceRecord, Strand

#: Selection markers this package can name an antibiotic for, keyed by the feature name lowered.
#: pUC19's is the one `docs/research/golden-gate-assembly.md` §8 states a plate recipe for; a
#: marker absent from here is named rather than translated.
SELECTION: Mapping[str, str] = {
    "ampr": "ampicillin or carbenicillin",
    "bla": "ampicillin or carbenicillin",
}


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


def read_phenotype(
    product: SequenceRecord,
    insert: tuple[int, int],
    *,
    vector: SequenceRecord,
    span: tuple[int, int],
) -> Phenotype:
    """Read what the product says about itself off its own features.

    Parameters
    ----------
    product
        The plasmid a clone carries.
    insert
        Where the inserts lie in it, from the first junction to the last.
    vector
        The plasmid they went into.
    span
        The vector bases they replace.
    """
    first, last = insert
    coding = _coding(product, first, last)
    promoter, gap = _promoter(product, first, last)
    return Phenotype(
        (first, last),
        coding,
        promoter,
        gap,
        promoter is not None and coding is not None and promoter.strand == coding.strand,
        _ribosome_binding_site(product, promoter, first, last),
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
