"""The inserts the many-part Golden Gate case is built from, and which way each goes in.

They are written here rather than read from a file: neither carries a BbsI site, and each begins
on bases no other junction of the set can take. The first four bases of an insert are the
overhang its junction takes, and a reversed one's junction takes the four its other strand
spells.
"""

from typing import Literal

#: A flexible glycine-serine linker.
LINKER = (
    "AACGGTTCAGGTGGATCTGGCGGTTCTGGAGGCAGCGGTTCAGGAGGTTCTGGCGGATCA"
    "GGTGGTTCAGGAGGCTCAGGTTCTGGAGGATCTGGCGGTTCAGGAGGTTCTGGATCAGGT"
    "TCTGGAGGCAGCGGTTCAGGAGGATCTGGT"
)

#: A His-tagged spacer.
TAG = (
    "CTTGGTCACCATCACCATCACCATGGTTCTGGATCAGGTTCTGCTTGGAGCCATCCGCAA"
    "TTCGAAAAAGGTGGTTCTGGCGGATCAGGTTCTGGAGGCAGCTTCGGTTCAGGAGGATCT"
    "GGTGGTTCAGGAAGCTCAGGTTCTGGAGGT"
)

#: Which way each of GFP, the linker and the tag goes in. The linker turns round, so one plan
#: covers the many-part route and the reversed-insert route.
ORIENTATIONS: tuple[Literal["forward", "reverse"], ...] = ("forward", "reverse", "forward")


def fasta(name: str, sequence: str, width: int = 60) -> str:
    """Return the record as FASTA text, wrapped at `width` bases a line.

    Examples
    --------
    >>> fasta("Tiny", "ACGTACGT", width=4)
    '>Tiny\\nACGT\\nACGT\\n'
    """
    lines = (sequence[at : at + width] for at in range(0, len(sequence), width))
    return f">{name}\n" + "".join(f"{line}\n" for line in lines)
