"""The verdicts a Gateway plan carries beyond one reaction's own.

`liulab_mbio.cloning.gateway.recombination` judges each reaction from its product's bases. Here
are the ones that span the experiment: an att site inside the DNA that moves, the strain each
vector is grown and selected in, the markers each plate has to tell apart, and the reading frame
where a fusion reads through a junction.

A check no sourced threshold judges carries no verdict rather than a pass. That is what carries
the ccdB-resistant host: no manual forbids plating a reaction on one, and that the selection
would be lost there is an inference `docs/research/gateway-cloning.md` marks as one, so the plan
states it as guidance and judges nothing.

`liulab_mbio.checks` is the model this builds on -- a `Check` is one verdict and the value it
judged. Coordinates are the model's, 0-based and half-open.
"""

import re

from liulab_mbio.bench.phenotype import SELECTION, selection_marker
from liulab_mbio.checks import Check
from liulab_mbio.cloning.gateway.att import REGION_BP, find_att_sites
from liulab_mbio.cloning.gateway.bench import HOSTS, PROPAGATION_HOST
from liulab_mbio.cloning.gateway.design import Fusion
from liulab_mbio.cloning.gateway.recombination import Junction, Piece, PlannedReaction
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.translate import translate

#: How much of a C-terminal fusion junction the destination vector supplies, base pairs: attB2's
#: 25 bases carry only eight of the nine codons, and the vector's next base finishes the last one
#: (Hartley 2000 figure 1D, MAN0000470 pages 14-15; note §8). It is why this end cannot be judged
#: from the insert and attB2 alone.
C_VECTOR_BP = 1


def plan_checks(
    *, lr: PlannedReaction, bp: PlannedReaction | None, host: str, fusion: Fusion
) -> tuple[Check, ...]:
    """Return the verdicts that span the plan rather than one reaction.

    Each argument is the `liulab_mbio.cloning.gateway.plan.Plan` field or property of that name.
    The frame is judged only at an end `fusion` names, because the vendor states both rules as
    intentions -- "If you wish to fuse your PCR product in frame with an N-terminal tag"
    (MAN0000470 page 13) -- and no record says which coding sequence is a tag.
    """
    made = [
        _internal_att_sites((bp or lr).recombination.moved),
        _selection_host(host),
        _vector_host(host),
        _markers(lr=lr, bp=bp),
    ]
    if fusion in ("N-terminal", "both"):
        made.append(_n_frame(lr))
    if fusion in ("C-terminal", "both"):
        made.append(_c_frame(lr))
    return tuple(made)


def _internal_att_sites(moved: Piece) -> Check:
    """Count the att sites inside the DNA that moves, read on either strand."""
    found = find_att_sites(SequenceRecord(moved.bases, name=moved.name))
    where = ", ".join(
        f"{one.name} at {one.start} on the {_strand(one.strand)} strand" for one in found
    )
    return Check(
        "insert att sites",
        "pass" if not found else "fail",
        len(found),
        f"no att site inside the {moved.length} bp that move"
        if not found
        else f"{where}, which would recombine where nobody meant it to",
    )


def _selection_host(host: str) -> Check:
    """Whether the strain the reactions are transformed into still counter-selects ccdB.

    The value is how many of the prohibitions the manuals state the named strain trips, which is
    one or none. The F' rule is the only one they state, so it is the only one that judges: a
    strain the manuals name neither way is not judged at all, because no genotype can be read off
    a name.
    """
    kind = _host_kind(host)
    if kind == "f-prime":
        return Check(
            "ccdB host",
            "fail",
            1,
            f"{host} carries the F' episome, whose ccdA cancels the ccdB selection, so unreacted "
            "vector grows",
        )
    if kind == "sensitive":
        return Check(
            "ccdB host",
            "pass",
            0,
            f"CcdB kills {host}, so unreacted vector and by-product do not grow",
        )
    why = (
        f"{host} resists CcdB, and no manual forbids plating a reaction on such a strain"
        if kind == "resistant"
        else f"the manuals say nothing about {host}"
    )
    return Check(
        "ccdB host",
        None,
        0,
        f"{why}; confirm it carries no F' episome, whose ccdA cancels the selection",
    )


def _vector_host(host: str) -> Check:
    """Name what grows a ccdB vector, and carry the inference the sources do not.

    Never a verdict: the manuals require the propagation strain and forbid F', and none of them
    says what plating a reaction on a ccdB-resistant strain does.
    """
    said = (
        f"grow the donor and destination vectors in {PROPAGATION_HOST}, the only kind of strain "
        "that grows one"
    )
    if _host_kind(host) == "resistant":
        said += (
            f"; {host} is one of those, and no manual says whether a reaction plated on it still "
            "counter-selects"
        )
    return Check("ccdB vector host", None, 0, said)


def _markers(*, lr: PlannedReaction, bp: PlannedReaction | None) -> Check:
    """Whether the LR plate tells the expression clone from the entry clone that went into it.

    Every marker is read off its own record's features. The vendor states the rule as a
    requirement on the entry clone: "you will need to perform the LR recombination reaction with
    an entry clone that carries a selection marker other than the kanamycin resistance gene"
    (MAN0000470 page 33; note §13). The value is how many drugs the two plate on.
    """
    product = selection_marker(lr.product)
    entry = selection_marker(lr.recombination.moved.record)
    if product is None or entry is None:
        missing = "expression clone" if product is None else "entry clone"
        return Check(
            "LR markers", None, 0, f"the {missing} annotates no selection marker to select on"
        )
    detail = f"{_drug(product)} on the expression clone, {_drug(entry)} on the entry clone"
    donor = None if bp is None else selection_marker(bp.recombination.backbone.record)
    if donor is not None:
        detail += f", from the donor vector's {donor.name}"
    apart = _drug(product) != _drug(entry)
    return Check(
        "LR markers",
        "pass" if apart else "warn",
        2 if apart else 1,
        detail if apart else f"{detail}: one plate, so unreacted entry clone grows on it too",
    )


def _n_frame(lr: PlannedReaction) -> Check:
    """Whether an N-terminal fusion reads in frame through attB1 without a stop.

    Nine codons: attB1's 25 bases and the two the primer adds, the ninth being attB1's final
    ``T`` and those two, which is why they may not be ``AA``, ``AG`` or ``GA``. The frame runs on
    from what the vector carries upstream, so the nearest coding sequence ahead of the junction
    is measured against too. The value is how many bases stand between attB1 and the insert.
    """
    coding, why = _coding(lr)
    if coding is None:
        return Check("N-terminal frame", None, 0, why)
    junction, product = lr.junctions[0], lr.product
    added = (coding.segments[0].start - junction.end) % len(product)
    window = product.extract(Segment(junction.start, junction.start + REGION_BP + added))
    spelled = _spelled(window)
    said = f"{added} base(s) added after attB1"
    if spelled is None:
        return Check(
            "N-terminal frame",
            "fail",
            added,
            f"{said}, so the {len(window)} bases to {coding.name} are not whole codons",
        )
    if "*" in spelled:
        return Check("N-terminal frame", "fail", added, f"{said}, spelling {spelled}, a stop")
    upstream = _upstream(product, junction, coding)
    if upstream is not None and (junction.start - upstream.segments[0].start) % 3:
        return Check(
            "N-terminal frame",
            "fail",
            added,
            f"{said}, spelling {spelled}, but out of frame with {upstream.name} upstream",
        )
    with_it = "" if upstream is None else f", in frame with {upstream.name} upstream"
    return Check("N-terminal frame", "pass", added, f"{said}, spelling {spelled}{with_it}")


def _c_frame(lr: PlannedReaction) -> Check:
    """Whether a C-terminal fusion reads in frame through attB2 without a stop.

    Nine codons again: the one base the primer adds, attB2's 25, and `C_VECTOR_BP` from the
    destination vector. The value is how many bases stand between the insert and attB2. "Any
    in-frame stop codons between the attB2 site and your gene of interest must be removed"
    (MAN0000470 page 14), and the gene's own last codon is one of those.
    """
    coding, why = _coding(lr)
    if coding is None:
        return Check("C-terminal frame", None, 0, why)
    junction, product = lr.junctions[1], lr.product
    end = coding.segments[-1].end
    added = (junction.start - end) % len(product)
    if _spelled(product.extract(Segment(end - 3, end))) == "*":
        return Check(
            "C-terminal frame",
            "fail",
            added,
            f"{coding.name} ends in a stop codon, which a C-terminal fusion has to lose",
        )
    window = product.extract(Segment(end, end + added + REGION_BP + C_VECTOR_BP))
    spelled = _spelled(window)
    said = f"{added} base(s) added before attB2, and one from the vector after it"
    if spelled is None:
        return Check(
            "C-terminal frame",
            "fail",
            added,
            f"{said}, so the {len(window)} bases past {coding.name} are not whole codons",
        )
    if "*" in spelled:
        return Check("C-terminal frame", "fail", added, f"{said}, spelling {spelled}, a stop")
    return Check("C-terminal frame", "pass", added, f"{said}, spelling {spelled}")


def _coding(lr: PlannedReaction) -> tuple[Feature | None, str]:
    """Return the insert's coding sequence, or why no frame reads across a junction."""
    coding = lr.phenotype.coding
    if coding is None:
        return None, "the insert annotates no coding sequence to read a frame from"
    if coding.strand == Strand.REVERSE:
        return None, f"{coding.name} reads on the other strand, so no fusion reads through here"
    return coding, ""


def _spelled(window: str) -> str | None:
    """Return what a junction window translates as, or ``None`` where it is not whole codons."""
    return translate(window) if window and len(window) % 3 == 0 else None


def _upstream(product: SequenceRecord, junction: Junction, coding: Feature) -> Feature | None:
    """Return the coding sequence standing nearest ahead of the junction, on the insert's strand."""
    length = len(product)
    found, gap = None, length
    for feature in product.features:
        if feature.type != "CDS" or feature is coding or feature.strand != Strand.FORWARD:
            continue
        distance = (junction.start - max(one.end for one in feature.segments)) % length
        if distance < gap:
            found, gap = feature, distance
    return found


def _drug(marker: Feature) -> str:
    """Return what a marker plates on, or its own name where `SELECTION` names no drug."""
    return SELECTION.get(marker.name.lower(), marker.name)


def _host_kind(host: str) -> str:
    """Return what the manuals say about this strain, or ``""`` where they name none.

    The longest name matches first, so ``TOP10F'`` is not read as ``TOP10``.
    """
    name = _plain(host)
    for known in sorted(HOSTS, key=len, reverse=True):
        if _plain(known) in name:
            return HOSTS[known]
    return ""


def _plain(name: str) -> str:
    """Return a strain's name with the punctuation and spacing two catalogues differ in gone."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def _strand(strand: Strand) -> str:
    """Name a strand the way a check's detail reads it."""
    return "reverse" if strand == Strand.REVERSE else "top"
