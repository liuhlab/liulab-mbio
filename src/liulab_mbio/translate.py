"""Write a protein as DNA for one host, and take the forbidden sites out of a coding sequence.

Reverse translation spells each amino acid with the codon a host's usage table counts most
often, so one protein and one host always give the same bases. A forbidden site is then taken
away by synonymous codon change, on either strand, and a change that would spell another
forbidden site is passed over. A site no synonymous change can reach is refused rather than
left in, because a part carrying it would cut itself.

The host is named and never assumed: `liulab_mbio.codons.codon_tables` lists the tables that
ship. A coding sequence the caller already holds is checked and left as it is, not written
again.
"""

from collections.abc import Iterable, Mapping
from dataclasses import KW_ONLY, dataclass

from liulab_mbio.codons import CodonUsage, amino_acid, codon_usage
from liulab_mbio.enzymes import Enzyme, get_enzyme
from liulab_mbio.sequence import Feature, Segment, SequenceRecord, Strand
from liulab_mbio.sites import CutSite, Domestication, EnzymeLike, domesticate, find_sites

#: What a sequence is called when the caller names none.
DEFAULT_NAME = "coding sequence"


class SiteNotRemovableError(ValueError):
    """A forbidden site no synonymous codon change could take away.

    A `ValueError`, so a caller catching that catches this. The message names the sequence,
    the site and the reason.
    """


@dataclass(frozen=True, slots=True)
class CodingSequence:
    """One synthesis-ready coding sequence, and what writing it changed.

    Parameters
    ----------
    name
        What this sequence is called, as the report and any refusal name it.
    dna
        The coding sequence, 5' to 3', a whole number of codons and free of every forbidden
        site on both strands.
    protein
        The amino acids `dna` spells, a stop written ``"*"``.
    host
        The codon usage table it was written for.
    forbidden
        The enzymes whose sites it is free of, by name.
    changes
        Each synonymous codon change that took a forbidden site away, in the order they were
        made. `liulab_mbio.sites.Domestication` carries the codon, the amino acid and the site.
    recoded
        ``True`` when the codons were written from a protein, ``False`` when a coded sequence
        was checked and its codons left alone.
    """

    name: str
    dna: str
    protein: str
    _: KW_ONLY
    host: str
    forbidden: tuple[str, ...] = ()
    changes: tuple[Domestication, ...] = ()
    recoded: bool = False


def translate(dna: str) -> str:
    """Return the amino acids `dna` spells, a stop written ``"*"``.

    Raises
    ------
    ValueError
        If `dna` is empty, is not a whole number of codons, or holds a base outside ACGT.

    Examples
    --------
    >>> translate("ATGTGGTAA")
    'MW*'
    """
    bases = _checked_dna(dna)
    return "".join(amino_acid(bases[at : at + 3]) for at in range(0, len(bases), 3))


def reverse_translate(protein: str, *, host: str) -> str:
    """Write `protein` as DNA, each amino acid as the codon `host` counts most often.

    Parameters
    ----------
    protein
        One-letter amino acids, a stop written ``"*"``, read whatever the case.
    host
        The name of a codon usage table, such as ``"e-coli-k12"``.

    Raises
    ------
    KeyError
        If no shipped table is called `host`.
    ValueError
        If `protein` is empty, holds a letter that is not an amino acid, or holds a stop
        anywhere but its last position.

    Examples
    --------
    >>> reverse_translate("MW*", host="e-coli-k12")
    'ATGTGGTAA'
    """
    favourites = _favourites(codon_usage(host))
    return "".join(favourites[one] for one in _checked_protein(protein, favourites))


def optimize_protein(
    protein: str,
    *,
    host: str,
    forbidden: Iterable[EnzymeLike] = (),
    name: str = DEFAULT_NAME,
) -> CodingSequence:
    """Write `protein` as a coding sequence for `host` with no forbidden site in it.

    The codons are the ones `host` counts most often, and a forbidden site they spell is taken
    away by synonymous change, so the protein that comes back is the protein that went in.

    Parameters
    ----------
    protein
        One-letter amino acids, a stop written ``"*"``.
    host
        The name of a codon usage table.
    forbidden
        Enzymes whose recognition sites the DNA must not spell, on either strand. Each an
        `Enzyme` or a name `liulab_mbio.enzymes.get_enzyme` answers to.
    name
        What to call this sequence in the report and in a refusal.

    Returns
    -------
    CodingSequence
        The DNA, the protein it spells, and every codon change that took a site away.

    Raises
    ------
    KeyError
        If no shipped table is called `host`, or no shipped enzyme answers to a name.
    ValueError
        If `protein` is not one. `SiteNotRemovableError` when a forbidden site is left that no
        synonymous change can reach.

    Examples
    --------
    >>> gene = optimize_protein("MSGRW*", host="e-coli-k12", forbidden=["NotI"])
    >>> gene.protein
    'MSGRW*'
    >>> (gene.changes[0].old_codon, gene.changes[0].new_codon)
    ('CGC', 'CGT')
    """
    return _cleaned(
        reverse_translate(protein, host=host),
        host=host,
        forbidden=forbidden,
        name=name,
        recoded=True,
    )


def optimize_coding_sequence(
    dna: str,
    *,
    host: str,
    forbidden: Iterable[EnzymeLike] = (),
    name: str = DEFAULT_NAME,
) -> CodingSequence:
    """Check a coding sequence someone already holds, and take any forbidden site out of it.

    A codon spelling no forbidden site is left exactly as it stands: a sequence already coded
    is checked, never written again. Only the codons a site sits on move, and they move to a
    synonym.

    Parameters
    ----------
    dna
        The coding sequence, 5' to 3', definite bases and a whole number of codons.
    host
        The name of a codon usage table, which decides which synonym a change reaches for.
    forbidden
        Enzymes whose recognition sites the DNA must not spell, on either strand.
    name
        What to call this sequence in the report and in a refusal.

    Returns
    -------
    CodingSequence
        The DNA, the protein it spells, and every codon change that took a site away.

    Raises
    ------
    KeyError
        If no shipped table is called `host`, or no shipped enzyme answers to a name.
    ValueError
        If `dna` is not a coding sequence: an odd length, a base outside ACGT, or a stop
        before the end. `SiteNotRemovableError` when a forbidden site cannot be taken away.

    Examples
    --------
    >>> gene = optimize_coding_sequence("ATGGGTCTCAAGACCTAA", host="e-coli-k12",
    ...                                forbidden=["BsaI", "BbsI"])
    >>> gene.dna
    'ATGGGCCTCAAGACCTAA'
    >>> gene.recoded
    False
    """
    return _cleaned(_checked_dna(dna), host=host, forbidden=forbidden, name=name, recoded=False)


def _cleaned(
    dna: str,
    *,
    host: str,
    forbidden: Iterable[EnzymeLike],
    name: str,
    recoded: bool,
) -> CodingSequence:
    """Take every forbidden site out of `dna` by synonymous change, or refuse to."""
    usage = codon_usage(host)
    enzymes = _resolved(forbidden)
    protein = translate(dna)
    _no_internal_stop(protein)
    record, report = domesticate(_coding_record(dna, name), enzymes, usage=usage)
    if left := find_sites(record, enzymes):
        raise SiteNotRemovableError(_why(left[0], record.sequence, usage, name))
    return CodingSequence(
        name,
        record.sequence,
        protein,
        host=usage.name,
        forbidden=tuple(one.name for one in enzymes),
        changes=report.changes,
        recoded=recoded,
    )


def _coding_record(dna: str, name: str) -> SequenceRecord:
    """Annotate the whole sequence as one coding sequence, which is what a synthesis part is."""
    return SequenceRecord(
        dna,
        name=name,
        features=(Feature(name, "CDS", (Segment(0, len(dna)),), strand=Strand.FORWARD),),
    )


def _why(site: CutSite, dna: str, usage: CodonUsage, name: str) -> str:
    """Name the sequence, the site and why no synonymous change took it away."""
    indices = sorted({at // 3 for at in range(site.start, site.end)})
    codons = [dna[3 * at : 3 * at + 3] for at in indices]
    covered = ", ".join(
        f"{at} {codon} ({usage.amino_acid(codon)})"
        for at, codon in zip(indices, codons, strict=True)
    )
    cause = (
        "each spells its amino acid the only way there is"
        if all(len(usage.synonymous(codon)) == 1 for codon in codons)
        else "every synonymous change there spells another forbidden site"
    )
    strand = "forward" if site.strand == Strand.FORWARD else "reverse"
    return (
        f"{name}: the {site.enzyme.name} site at {site.start}-{site.end} on the {strand} "
        f"strand cannot be removed: it covers codon(s) {covered}, and {cause}"
    )


def _resolved(forbidden: Iterable[EnzymeLike]) -> tuple[Enzyme, ...]:
    """Read each forbidden enzyme, by name or by record, refusing a name nothing answers to."""
    return tuple(one if isinstance(one, Enzyme) else get_enzyme(one) for one in forbidden)


def _favourites(usage: CodonUsage) -> Mapping[str, str]:
    """Each amino acid, and the codon this host counts most often for it."""
    favourites: dict[str, str] = {}
    for codon in sorted(usage.counts):
        favourites.setdefault(usage.amino_acid(codon), usage.synonymous(codon)[0])
    return favourites


def _checked_protein(protein: str, amino_acids: Iterable[str]) -> str:
    """Upper-case a protein, or refuse something that is not one."""
    upper = protein.upper()
    if not upper:
        raise ValueError("a protein of no amino acids cannot be written as DNA")
    if bad := sorted(set(upper) - set(amino_acids)):
        raise ValueError(f"not one-letter amino acids: {''.join(bad)}")
    _no_internal_stop(upper)
    return upper


def _checked_dna(dna: str) -> str:
    """Upper-case a coding sequence, or refuse something that cannot be read as codons."""
    upper = dna.upper()
    if bad := sorted(set(upper) - set("ACGT")):
        raise ValueError(
            f"a sequence to synthesise needs definite bases, and {''.join(bad)} is not one of ACGT"
        )
    if not upper or len(upper) % 3:
        raise ValueError(f"{len(upper)} bases is not a whole number of codons")
    return upper


def _no_internal_stop(protein: str) -> None:
    """Refuse a stop before the end, which would truncate the protein."""
    if "*" in protein[:-1]:
        raise ValueError(
            f"a stop at amino acid {protein.index('*')} would truncate the protein; "
            "a stop belongs at the end or nowhere"
        )
