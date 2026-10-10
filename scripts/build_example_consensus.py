#!/usr/bin/env python3
"""Write the sequencing result the sequence verification page checks.

    pixi run python scripts/build_example_consensus.py [--from PATH] [--out PATH]

It is the worked Golden Gate product's own bases with two changed: one inside GFP, the insert,
and one inside AmpR, in the backbone. Held against the product, that one consensus gives all
three outcomes the page explains: both junctions pass, GFP fails, and the AmpR change is listed
with its feature and no verdict.

No real clone of the worked example was ever sequenced, so the result is made rather than
downloaded. The file is a FASTA, as a whole-plasmid service sends its consensus.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mbio.io import read_record
from mbio.sequence import SequenceRecord

REPO = Path(__file__).resolve().parents[1]

#: The record the consensus is made from: the worked Golden Gate plan's product.
PRODUCT = REPO / "docs/examples/pUC19-GFP/product.dna"

#: Where the consensus goes, beside the product it is checked against.
OUT = REPO / "docs/examples/pUC19-GFP/verification/consensus.fasta"

#: Each change: the feature it falls in, and how far past that feature's first base it sits.
PLANTED = (("GFP", 300), ("AmpR", 400))

#: Each base and the one a transition turns it into.
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}

#: Bases per line, as a FASTA is usually wrapped.
WIDTH = 60


def planted(record: SequenceRecord) -> tuple[str, list[tuple[str, int, str, str]]]:
    """Return the record's bases with each change in `PLANTED` made, and what each one was."""
    bases = list(record.sequence.upper())
    made = []
    for name, offset in PLANTED:
        (feature,) = [one for one in record.features if one.name == name]
        position = min(segment.start for segment in feature.segments) + offset
        if not any(segment.start <= position < segment.end for segment in feature.segments):
            raise ValueError(f"{position} does not lie inside {name}")
        # A feature across the origin ends past the record's length (ADR 0001).
        position %= len(bases)
        before = bases[position]
        bases[position] = TRANSITION[before]
        made.append((name, position, before, bases[position]))
    return "".join(bases), made


def main(argv: list[str] | None = None) -> int:
    """Write the consensus, printing each change it carries."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from", dest="source", type=Path, default=PRODUCT, help="the product")
    parser.add_argument("--out", type=Path, default=OUT, help="where the consensus goes")
    args = parser.parse_args(argv)
    if not args.source.exists():
        print(f"{args.source}: no such file; run the Golden Gate plan first", file=sys.stderr)
        return 1
    record = read_record(args.source)
    bases, made = planted(record)
    lines = [bases[at : at + WIDTH] for at in range(0, len(bases), WIDTH)]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(f">{record.name} consensus\n" + "\n".join(lines) + "\n", encoding="utf-8")
    print(f"{args.out}: {record.name} consensus, {len(bases)} bp")
    for name, position, before, after in made:
        print(f"  {name} at {position}: {before} -> {after}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
