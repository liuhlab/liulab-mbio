"""How often a host spells each codon, and the genetic code that groups them.

A table here counts codons over the complete coding sequences of one genome, one transcript per
gene where a gene has several, so it is a measurement of that genome rather than a copy of a
published compilation. `scripts/build_codon_usage.py` rebuilds it and
`docs/research/codon-usage.md` says where the sequences came from.

A whole-genome table is the background the genome itself uses. It is not a highly expressed
reference set, which is what a codon adaptation index wants and is a different object.
"""

import json
from collections.abc import Mapping
from dataclasses import KW_ONLY, dataclass, field
from functools import cache
from importlib.resources import files

#: The table used when a caller names none.
DEFAULT_TABLE = "e-coli-k12"


@dataclass(frozen=True, slots=True)
class CodonUsage:
    """How often one genome spells each of the 64 codons.

    Parameters
    ----------
    name
        The short name this table is asked for by, such as ``"e-coli-k12"`` or ``"human"``.
    organism
        The organism as its genome record names it.
    taxid, accession
        The NCBI taxonomy identifier, and the sequence record or genome assembly counted.
    counts
        Codon to the number of times the coding sequences spell it. All 64 are present.
    cds_count, codon_count
        How many coding sequences were counted, and how many codons they held.
    note
        What the table is and is not, for a reader choosing between tables.
    """

    name: str
    organism: str
    _: KW_ONLY
    taxid: int
    accession: str
    counts: Mapping[str, int] = field(hash=False)
    cds_count: int
    codon_count: int
    note: str = ""

    def amino_acid(self, codon: str) -> str:
        """Return the amino acid this codon spells, or ``"*"`` for a stop.

        Raises
        ------
        KeyError
            If `codon` is not one of the 64.
        """
        return amino_acid(codon)

    def fraction(self, codon: str) -> float:
        """Return this codon's share of the codons spelling the same amino acid, 0 to 1."""
        family = _families()[self.amino_acid(codon)]
        return self.counts[_checked(codon)] / sum(self.counts[one] for one in family)

    def per_thousand(self, codon: str) -> float:
        """How often this codon occurs in a thousand codons of the genome."""
        return 1000 * self.counts[_checked(codon)] / self.codon_count

    def synonymous(self, codon: str) -> tuple[str, ...]:
        """Every codon for the same amino acid, the one this genome uses most often first.

        Examples
        --------
        >>> codon_usage().synonymous("GAC")
        ('GAT', 'GAC')
        >>> codon_usage("human").synonymous("GAC")
        ('GAC', 'GAT')
        """
        family = _families()[self.amino_acid(codon)]
        return tuple(sorted(family, key=lambda one: (-self.counts[one], one)))


def codon_tables() -> tuple[str, ...]:
    """Return the name of every codon usage table the package ships.

    Examples
    --------
    >>> codon_tables()
    ('e-coli-k12', 'human', 'mouse')
    """
    return tuple(_shipped())


def codon_usage(name: str = DEFAULT_TABLE) -> CodonUsage:
    """Return one host's codon usage.

    Raises
    ------
    KeyError
        If no shipped table is called `name`.

    Examples
    --------
    >>> codon_usage().organism
    'Escherichia coli str. K-12 substr. MG1655'
    """
    shipped = _shipped()
    if name not in shipped:
        raise KeyError(f"no shipped codon usage table is called {name!r}")
    return shipped[name]


def amino_acid(codon: str) -> str:
    """Return the amino acid this codon spells, or ``"*"`` for a stop.

    The genetic code groups the codons every table counts, so this asks for no host.

    Raises
    ------
    KeyError
        If `codon` is not one of the 64.

    Examples
    --------
    >>> amino_acid("atg")
    'M'
    """
    return _code()[_checked(codon)]


@cache
def _code() -> Mapping[str, str]:
    """Return the standard genetic code, stops written ``"*"``.

    It serves bacterial hosts too: the bacterial table assigns the same amino acids and differs
    only in which codons may start a gene, which is not what a codon usage table is asked about.
    """
    from Bio.Data.CodonTable import unambiguous_dna_by_id

    table = unambiguous_dna_by_id[1]
    return {**table.forward_table, **dict.fromkeys(table.stop_codons, "*")}


@cache
def _families() -> Mapping[str, tuple[str, ...]]:
    """Each amino acid, and the codons that spell it."""
    families: dict[str, list[str]] = {}
    for codon, amino in sorted(_code().items()):
        families.setdefault(amino, []).append(codon)
    return {amino: tuple(codons) for amino, codons in families.items()}


def _checked(codon: str) -> str:
    """Upper-case a codon, or refuse something that is not one."""
    upper = codon.upper()
    if upper not in _code():
        raise KeyError(f"{codon!r} is not one of the 64 codons")
    return upper


@cache
def _shipped() -> Mapping[str, CodonUsage]:
    data = json.loads((files("liulab_mbio") / "data" / "codon_usage.json").read_text())
    return {
        entry["name"]: CodonUsage(
            entry["name"],
            entry["organism"],
            taxid=entry["taxid"],
            accession=entry["accession"],
            counts=dict(entry["counts"]),
            cds_count=entry["cds_count"],
            codon_count=entry["codon_count"],
            note=entry["note"],
        )
        for entry in data["tables"]
    }
