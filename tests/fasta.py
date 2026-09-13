"""A FASTA and the `.fai` index samtools would write beside it, as the genome tests need one."""

from collections.abc import Mapping
from pathlib import Path


def write_fasta(path: Path, records: Mapping[str, str], width: int = 60) -> None:
    """Write the records wrapped at `width` bases a line, with their index beside them.

    The index columns are samtools': name, length, the byte the bases start at, bases per line
    and bytes per line.
    """
    text = ""
    index = []
    for name, bases in records.items():
        text += f">{name} planted\n"
        index.append(f"{name}\t{len(bases)}\t{len(text)}\t{width}\t{width + 1}")
        text += "".join(f"{bases[at : at + width]}\n" for at in range(0, len(bases), width))
    path.write_text(text)
    path.with_name(path.name + ".fai").write_text("\n".join(index) + "\n")
