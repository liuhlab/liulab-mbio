#!/usr/bin/env python3
"""Write the page sequence verification writes for the worked example's clone.

    pixi run python scripts/build_example_verification.py [--result PATH] [--out PATH]

It holds the consensus `build_example_consensus.py` writes against the worked Golden Gate
product, as `mbio sequence-verify --out` does, and writes the page through `draw_page`. The
command exits with 1 on a clone it does not verify, which this one is not, so the page is
written here rather than by the command.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mbio.io import read_record
from mbio.verification import draw_page, regions, verify
from mbio.verification.cli import PAGE, read_result

REPO = Path(__file__).resolve().parents[1]

#: The record the clone should be: the worked Golden Gate plan's product.
PRODUCT = REPO / "docs/examples/pUC19-GFP/product.dna"

#: The consensus held against it, and where the page goes, beside it.
RESULT = REPO / "docs/examples/pUC19-GFP/verification/consensus.fasta"
OUT = RESULT.parent / PAGE


def main(argv: list[str] | None = None) -> int:
    """Write the page, printing where it went and the clone's verdict."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT, help="the consensus")
    parser.add_argument("--out", type=Path, default=OUT, help="where the page goes")
    args = parser.parse_args(argv)
    for needed in (PRODUCT, args.result):
        if not needed.exists():
            print(f"{needed}: no such file; run the plan and the consensus first", file=sys.stderr)
            return 1
    record = read_record(PRODUCT)
    result = read_result(args.result)
    made = verify(record, [result], regions(record))
    written = draw_page(record, [result], made).write(args.out)
    print(f"{written}: {'verified' if made.verified else 'not verified'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
