"""A sequencing result: one sample's called bases, and what the route gave to trust them by.

Every route hands one over in this shape: a reader of a Sanger trace, of a whole-plasmid
consensus, or a DMX well's call. Each converts at its own boundary, so nothing past here knows
which file a result came from.
"""

from dataclasses import KW_ONLY, dataclass

from mbio.sequence import IUPAC_DNA


def _bases(sequence: str, owner: str) -> str:
    upper = sequence.upper()
    if bad := set(upper) - IUPAC_DNA:
        raise ValueError(f"{owner} holds letters that are not IUPAC DNA: {''.join(sorted(bad))}")
    return upper


@dataclass(frozen=True, slots=True)
class SequencingResult:
    """One sample's sequencing, read into the bases called for it.

    Parameters
    ----------
    name
        What the sample is called, such as the read's file name.
    bases
        The called bases, stored upper-case. An ``N`` is no call; any other ambiguity code is a
        mixed base.
    quality
        A Phred value per base, where the route gives one.
    depth
        How many reads stand behind each base, where the route gives it.
    reads
        How many reads stood behind the whole.
    trusted
        The stretch worth trusting, half-open in the result's own coordinates. It may be empty,
        when a reader trusts nothing; ``None`` trusts all of it.
    others
        Any further consensus the same sample gave, which makes it mixed. Stored upper-case.

    Raises
    ------
    ValueError
        If `bases` or one of `others` is not IUPAC DNA, `quality` or `depth` is not one value a
        base, or `trusted` falls outside the bases.

    Examples
    --------
    >>> SequencingResult("A1", "acgt", trusted=(1, 3)).bases
    'ACGT'
    """

    name: str
    bases: str
    _: KW_ONLY
    quality: tuple[int, ...] | None = None
    depth: tuple[int, ...] | None = None
    reads: int | None = None
    trusted: tuple[int, int] | None = None
    others: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Upper-case the bases and refuse a field that does not fit them."""
        object.__setattr__(self, "bases", _bases(self.bases, f"result {self.name!r}"))
        object.__setattr__(
            self,
            "others",
            tuple(_bases(one, f"a further consensus of {self.name!r}") for one in self.others),
        )
        length = len(self.bases)
        for field, values in (("quality", self.quality), ("depth", self.depth)):
            if values is not None and len(values) != length:
                raise ValueError(
                    f"result {self.name!r} carries {len(values)} {field} values for {length} bases"
                )
        if self.trusted is not None:
            start, end = self.trusted
            if not 0 <= start <= end <= length:
                raise ValueError(
                    f"result {self.name!r} trusts {start}-{end}, which falls outside its "
                    f"{length} bases"
                )

    @property
    def trusted_span(self) -> tuple[int, int]:
        """The trusted stretch, half-open: all of the bases when nothing narrowed it."""
        return self.trusted if self.trusted is not None else (0, len(self.bases))
