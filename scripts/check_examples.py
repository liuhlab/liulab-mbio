#!/usr/bin/env python3
"""Check every published example against the command that writes it.

    pixi run examples-check

A page under `docs/examples/` tells its reader that each file beside it was written by the
command printed above it. Nothing else holds that promise, so the code can move and the
committed bytes stay behind. This regenerates each example into a temporary directory and
compares the bytes, reporting every difference rather than the first.

A red run is fixed by running the example's own command and committing what it writes. The
check never writes into `docs/examples/`.

The generators are named below rather than discovered. One of them is a script rather than a
documented command, and another writes files that `tests/synbio/igga/test_gate.py` reads as its
known-good corpus; neither is reachable by reading a page. A fifth generator is a line added
here, and an example that grows a file fails here until `writes` names it. That is the check
working: knowing what a command writes without running it is what no discovery rule can do.

This is a step of the `docs` CI job and not of `pixi run check`, so the gate pays nothing for
the seconds the example commands take.
"""

from __future__ import annotations

import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[1]
PUC19 = "docs/examples/pUC19-GFP"
AP1 = "docs/examples/ap1-library"
READBACK = "docs/examples/ap1-readback"

#: Stands in a command for the directory that run writes into. Each command below is otherwise
#: the one its example's page prints, so the two can be read against each other.
OUT = "<out>"


@dataclass(frozen=True)
class Generator:
    """One command set that writes part of a published example, and the files it owns."""

    what: str
    directory: Path
    commands: tuple[str, ...]
    writes: tuple[str, ...]


GENERATORS: tuple[Generator, ...] = (
    Generator(
        what="the pUC19-GFP Golden Gate plan",
        directory=REPO / PUC19,
        commands=(
            f"liulab_mbio cloning goldengate plan "
            f"tests/data/pUC19.dna tests/data/GFP.dna --out {OUT}",
        ),
        writes=("primers.tsv", "product.dna", "protocol.html", "protocol.json"),
    ),
    Generator(
        what="the AP-1 library plan",
        directory=REPO / AP1,
        commands=(
            f"liulab_synbio igga plan {AP1}/project.json --out {OUT} "
            f"--working-site EGFP --prices {AP1}/prices.csv",
        ),
        writes=(
            "barcodes.tsv",
            "block-vector-1.dna",
            "block-vector-2.dna",
            "block-vector-3.dna",
            "changes.tsv",
            "library-read-primers.tsv",
            "oligo.dna",
            "parts.tsv",
            "pool-primers.tsv",
            "pool.tsv",
            "product.dna",
            "protocol/01-primer-plates.html",
            "protocol/02-cargo-ordering-and-pool-preparation.html",
            "protocol/03-cargo-creation.html",
            "protocol/04-cargo-validation-barcode-ligation.html",
            "protocol/04-cargo-validation-index-pcr.html",
            "protocol/05-library-assembly-in-rounds.html",
            "protocol/06-final-cargo-ligation.html",
            "protocol/index.html",
            "protocol/project.json",
            "protocol/reagents.html",
            "protocol/references.html",
            "round-1.dna",
            "round-2.dna",
        ),
    ),
    Generator(
        what="the AP-1 figures",
        directory=REPO / AP1,
        commands=(
            f"liulab_mbio plot map {AP1}/vector.gb -o {OUT}/vector-map.pdf",
            f"liulab_mbio plot map {AP1}/round-1.dna -o {OUT}/round-1-map.pdf",
            f"liulab_mbio plot map {AP1}/round-2.dna -o {OUT}/round-2-map.pdf",
            f"liulab_mbio plot map {AP1}/product.dna -o {OUT}/product-map.pdf",
            f"liulab_mbio plot map {AP1}/product.dna "
            f"--region 1368..1442 --sequence-view -o {OUT}/barcode-block.pdf",
        ),
        writes=(
            "barcode-block.pdf",
            "product-map.pdf",
            "round-1-map.pdf",
            "round-2-map.pdf",
            "vector-map.pdf",
        ),
    ),
    Generator(
        what="the iGGA destination",
        directory=REPO / AP1,
        commands=(f"python scripts/build_dmx_vector.py --out {OUT}/vector.gb",),
        writes=("vector.gb",),
    ),
    Generator(
        what="the domesticated working vector and its protocol",
        directory=REPO / AP1,
        commands=(f"python scripts/build_working_vector.py --out {OUT}",),
        writes=(
            "working-vector-domestication.html",
            "working-vector-domestication.json",
            "working-vector.gb",
        ),
    ),
    Generator(
        what="the AP-1 cargo read-back",
        directory=REPO / READBACK / "protocol",
        commands=(f"liulab_synbio dmx plan {READBACK}/build.json --out {OUT}",),
        writes=(
            "01-design-read-back-index-pcr.html",
            "index.html",
            "project.json",
            "reagents.html",
            "references.html",
        ),
    ),
)


def run(generator: Generator, into: Path) -> list[str]:
    """Run every command of `generator`, returning what each one that failed printed."""
    failures = []
    for command in generator.commands:
        spelled = [part.replace(OUT, str(into)) for part in shlex.split(command)]
        # `python` means this interpreter, so a script runs in the environment the check runs in.
        spelled[:1] = [sys.executable if spelled[0] == "python" else spelled[0]]
        done = subprocess.run(spelled, cwd=REPO, capture_output=True, text=True, check=False)
        if done.returncode != 0:
            printed = (done.stderr or done.stdout).strip()
            failures.append(f"`{command}` exited {done.returncode}: {printed}")
    return failures


def compare(generator: Generator, into: Path) -> list[str]:
    """Return one line for every file of `generator` the command no longer writes as committed."""
    differences = []
    for name in generator.writes:
        fresh, committed = into / name, generator.directory / name
        if not fresh.exists():
            differences.append(f"{_said(committed)}: the command wrote no such file")
        elif not committed.exists():
            differences.append(f"{_said(committed)}: written by the command, not committed")
        elif fresh.read_bytes() != committed.read_bytes():
            differences.append(f"{_said(committed)}: the command writes other bytes")
    for found in sorted(one for one in into.rglob("*") if one.is_file()):
        if str(found.relative_to(into)) not in generator.writes:
            differences.append(
                f"{_said(generator.directory / found.relative_to(into))}: written by the "
                "command, and this check does not name it"
            )
    return differences


def _said(path: Path) -> str:
    """Spell a path the way a reader standing in the repository would type it."""
    return str(path.relative_to(REPO))


def main() -> int:
    """Check every generator, print all of them, and fail if any example has drifted."""
    drifted = 0
    for generator in GENERATORS:
        with TemporaryDirectory() as made:
            into = Path(made)
            reported = run(generator, into) or compare(generator, into)
        if reported:
            drifted += 1
            print(f"DRIFTED    {generator.what}")
            for line in reported:
                print(f"           {line}")
        else:
            print(f"ok         {generator.what}, {len(generator.writes)} files")
    if drifted:
        print(
            f"\n{drifted} of {len(GENERATORS)} generators no longer write what is committed. "
            "Run the command on the example's own page and commit what it writes.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
