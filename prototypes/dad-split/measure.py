"""Run the two splitters over the AP-1 example and print the tables issue #261 asks for.

Usage, from the repository root::

    pixi run python prototypes/dad-split/measure.py

Every number it prints is produced here, on the inputs named in the heading above each table.
Nothing is copied from a source except Lund's clone rates, which are labelled as theirs.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from collections.abc import Callable, Sequence
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from split import (
    Budget,
    Split,
    assembles,
    cargo_lengths,
    clever_split,
    fewest_fragments,
    fragment_windows,
    simple_split,
)

from liulab_mbio.io import read_record
from liulab_mbio.ligase import LigaseProfile, read_profile

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "docs/examples/ap1-library"
PROFILE = ROOT / "reference_docs/ligation-fidelity/potapov2018/FileS03_T4_18h_25C.xlsx"
RESERVED = ("AGGA", "TTCC")

# Lund et al. 2024's measured share of clones with no error, by fragment count. Theirs, not ours.
LUND = {2: 100.0, 3: 93.8, 5: 84.6, 8: 66.7, 12: 40.0, 16: 0.0}


def parts() -> list[tuple[str, str]]:
    """Read the 72 AP-1 synthesis blocks, name and sequence."""
    with (EXAMPLE / "parts.tsv").open() as handle:
        return [(row["name"], row["sequence"]) for row in csv.DictReader(handle, delimiter="\t")]


def product() -> str:
    """Read the largest library member, one assembled plasmid's insert span."""
    return str(read_record(EXAMPLE / "product.dna").sequence)


def profile() -> LigaseProfile | None:
    """Read the user-held T4 matrix, or ``None`` where this machine does not hold it."""
    return read_profile(PROFILE) if PROFILE.exists() else None


def frame(sequence: str) -> int:
    """Find the frame of the longest open reading frame, for the synonymous-edit run."""
    best, at = 0, 0
    for offset in range(3):
        start = None
        for index in range(offset, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]
            if codon == "ATG" and start is None:
                start = index
            if codon in ("TAA", "TAG", "TGA") and start is not None:
                if index - start > best:
                    best, at = index - start, start
                start = None
        if start is not None and len(sequence) - start > best:
            best, at = len(sequence) - start, start
    return at


def table(heading: str, columns: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    """Print one Markdown table."""
    print(f"\n### {heading}\n")
    print("| " + " | ".join(columns) + " |")
    print("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        print("| " + " | ".join(str(cell) for cell in row) + " |")


def sweep_fragments(sequence: str, ligase: LigaseProfile | None) -> None:
    """Fidelity against fragment count against runtime, on one 2.3 kb cargo.

    The oligo is widened per row to whatever the count needs, so the count is the only thing
    that moves. A 300 nt oligo cannot reach the low counts at all, which is the first finding.
    """
    span = len(sequence) - 4
    rows = []
    for count in range(2, 17):
        budget = Budget(oligo=54 + 24 + -(-span // count))
        plain = simple_split(
            sequence, budget=budget, profile=ligase, count=count, reserved=RESERVED
        )
        clever = clever_split(
            sequence,
            budget=budget,
            profile=ligase,
            count=count,
            reserved=RESERVED,
            cap=12,
            seconds_limit=10.0,
        )
        widest = max(len(window) for window in fragment_windows(len(sequence), budget, count))
        rows.append(
            [
                count,
                budget.oligo,
                widest,
                show(plain),
                f"{plain.seconds * 1000:.0f}",
                show(clever),
                f"{clever.seconds * 1000:.0f}",
                lund(count),
            ]
        )
    table(
        "1. Fidelity, fragment count and runtime — 2,276 bp cargo, oligo widened per row so "
        "that the fragment count is the only thing that moves",
        [
            "fragments",
            "oligo needed (nt)",
            "widest window (nt)",
            "simple fidelity",
            "simple ms",
            "clever fidelity",
            "clever ms",
            "Lund error-free clones (%)",
        ],
        rows,
    )


def lund(count: int) -> str:
    """Lund's measured rate at this fragment count, or a blank where they measured none."""
    return f"{LUND[count]:.1f}" if count in LUND else ""


def show(answer: Split) -> str:
    """One cell: the fidelity, or why there is none."""
    return f"{answer.value:.4f}" if answer.feasible else "none"


def sweep_caps(sequence: str, ligase: LigaseProfile | None) -> None:
    """Measure what the one optimality knob buys, at counts the window makes hard."""
    rows = []
    for count in (10, 12):
        for cap in (2, 4, 8, 16, 32, None):
            answer = clever_split(
                sequence,
                profile=ligase,
                count=count,
                reserved=RESERVED,
                cap=cap,
                seconds_limit=10.0,
            )
            rows.append(
                [
                    count,
                    "all" if cap is None else cap,
                    show(answer),
                    f"{answer.seconds * 1000:.0f}",
                    answer.scored,
                    answer.note or "proved optimal over the kept candidates",
                ]
            )
    table(
        "2. What the candidate cap buys and costs — same 2,276 bp cargo, 300 nt oligo",
        ["fragments", "candidates per cut", "fidelity", "ms", "fidelity calls", "outcome"],
        rows,
    )


def sweep_parts(ligase: LigaseProfile | None) -> None:
    """Both searches over the real range: 72 AP-1 parts, 133 to 1,149 bp."""
    rows = []
    for name, run in (
        ("simple (greedy first-fit, as `design_overhangs` is today)", runner(simple_split)),
        ("clever (interval DP + branch and bound, cap 12)", runner(clever_split, cap=12)),
    ):
        began = time.perf_counter()
        answers = [(sequence, run(sequence, ligase)) for _, sequence in parts()]
        seconds = time.perf_counter() - began
        good = [answer for _, answer in answers if answer.feasible]
        values = [answer.value for answer in good]
        spreads = [spread(sequence, answer) for sequence, answer in answers if answer.feasible]
        rows.append(
            [
                name,
                f"{len(good)}/{len(answers)}",
                sum(answer.fragments for _, answer in answers),
                f"{min(values):.4f}",
                f"{statistics.median(values):.4f}",
                f"{seconds * 1000 / len(answers):.1f}",
                f"{max(spreads):.0f}%",
                sum(assembles(sequence, answer) for sequence, answer in answers),
            ]
        )
    table(
        "3. The real range — all 72 AP-1 parts, 133 to 1,149 bp, 300 nt oligo",
        [
            "search",
            "parts solved",
            "oligos in total",
            "worst fidelity",
            "median fidelity",
            "ms per part",
            "widest length spread",
            "parts that reassemble",
        ],
        rows,
    )


def spread(sequence: str, answer: Split) -> float:
    """How far the shortest oligo of one part falls below the longest, as a percentage."""
    lengths = cargo_lengths(sequence, answer)
    return 100 * (1 - min(lengths) / max(lengths)) if lengths else 0.0


def runner(
    search: Callable[..., Split], **extra: object
) -> Callable[[str, LigaseProfile | None], Split]:
    """Bind one search to the shared inputs."""

    def run(sequence: str, ligase: LigaseProfile | None) -> Split:
        return search(sequence, profile=ligase, reserved=RESERVED, **extra)

    return run


def sweep_oligo(ligase: LigaseProfile | None) -> None:
    """Where it falls over: shorten the oligo until the cargo no longer fits."""
    sequence = product()
    rows = []
    for oligo in (120, 150, 200, 250, 300, 350):
        budget = Budget(oligo=oligo)
        try:
            least = fewest_fragments(len(sequence), budget)
        except ValueError as error:
            rows.append([oligo, budget.ceiling, "-", "no partition", "-", str(error)[:40]])
            continue
        answer = clever_split(sequence, budget=budget, profile=ligase, reserved=RESERVED, cap=12)
        rows.append(
            [
                oligo,
                budget.ceiling,
                least,
                show(answer),
                f"{answer.seconds * 1000:.0f}",
                answer.note or "solved",
            ]
        )
    table(
        "4. Where it falls over — the oligo shortened under a fixed 2,276 bp cargo",
        ["oligo (nt)", "cargo per oligo (nt)", "fewest fragments", "fidelity", "ms", "outcome"],
        rows,
    )


def sweep_pressure(ligase: LigaseProfile | None) -> None:
    """Squeeze the oligo until no set exists, then try each of the four answers to that.

    No cap here. A cap ranks unedited candidates first, so it starves the edited ones and the
    edit column would measure the cap rather than the edit.
    """
    name, sequence = parts()[0]
    reading = frame(sequence)
    rows = []
    for oligo in (300, 150, 120, 110, 100):
        budget = Budget(oligo=oligo)
        shared = {
            "budget": budget,
            "profile": ligase,
            "reserved": RESERVED,
            "cap": None,
            "seconds_limit": 8.0,
        }
        least = fewest_fragments(len(sequence), budget)
        hard = clever_split(sequence, **shared)
        wider = clever_split(
            sequence,
            profile=ligase,
            reserved=RESERVED,
            cap=None,
            seconds_limit=8.0,
            budget=Budget(oligo=oligo + 60),
        )
        more = clever_split(sequence, count=least + 1, **shared)
        edited = clever_split(sequence, frame=reading, **shared)
        rows.append(
            [
                oligo,
                least,
                show(hard),
                show(wider),
                show(more),
                f"{show(edited)}, {edited.edits} edit(s)",
                "yes" if assembles(sequence, edited) else "-",
            ]
        )
    table(
        f"5. What it does when it cannot win — {name}, {len(sequence)} bp, the oligo squeezed",
        [
            "oligo (nt)",
            "fewest fragments",
            "fail (refuse, as asked)",
            "widen the oligo by 60 nt",
            "allow one more fragment",
            "allow synonymous edits",
            "edited cargo reassembles",
        ],
        rows,
    )


def sweep_reserved(ligase: LigaseProfile | None) -> None:
    """Hold more overhangs out by name until the search runs out of room."""
    sequence = product()
    held: list[str] = list(RESERVED)
    extra = ["AATG", "GCTT", "CCAT", "TGAC", "GGTA", "CATC", "ATCC", "AAGT", "CAGG", "TTAC"]
    rows = []
    for take in (0, 2, 4, 6, 8, 10):
        reserved = tuple(held + extra[:take])
        answer = clever_split(sequence, profile=ligase, reserved=reserved, cap=12)
        rows.append(
            [len(reserved), show(answer), f"{answer.seconds * 1000:.0f}", answer.note or "solved"]
        )
    table(
        "6. Overhangs held out by name — 2,276 bp cargo at its fewest fragments",
        ["reserved overhangs", "fidelity", "ms", "outcome"],
        rows,
    )


def sweep_scoring(sequence: str) -> None:
    """Compare the shipped BsaI matrix with the user-held T4 profile, same cargo."""
    rows = []
    for name, ligase in (("shipped BsaI matrix", None), ("FileS03_T4_18h_25C.xlsx", profile())):
        if ligase is None and name != "shipped BsaI matrix":
            rows.append([name, "not held on this machine", "-", "-"])
            continue
        answer = clever_split(sequence, profile=ligase, reserved=RESERVED, cap=12)
        rows.append([name, show(answer), "".join(answer.overhangs), f"{answer.seconds * 1000:.0f}"])
    table(
        "7. Which measurement scores the set — 2,276 bp cargo, same length bound",
        ["scored by", "fidelity", "overhangs chosen", "ms"],
        rows,
    )


def oligo_count() -> None:
    """Report the whole library's oligo count beside the pool bands. Reported, never optimised."""
    counts = [fewest_fragments(len(sequence), Budget()) for _, sequence in parts()]
    total = sum(counts)
    rows = [
        ["72 AP-1 parts at 300 nt", total, nearest(total)],
        [
            "72 AP-1 parts at 250 nt",
            (lesser := sum(fewest_fragments(len(s), Budget(oligo=250)) for _, s in parts())),
            nearest(lesser),
        ],
    ]
    table(
        "8. Oligos the split orders, beside the nearest pool band (reported, not optimised)",
        ["order", "oligos", "nearest band and headroom"],
        rows,
    )


def nearest(total: int) -> str:
    """Which pool band this order falls in, and how far it sits from the edge."""
    for band in (100, 500):
        if total <= band:
            return f"{band} ({band - total} spare)"
    return f"above 500 ({total - 500} over)"


def main() -> None:
    """Print every table."""
    ligase = profile()
    print("# Prototype measurements for issue #261\n")
    print(f"- Ligase profile: {'held' if ligase else 'ABSENT, shipped BsaI matrix used'}")
    print(f"- Reserved overhangs held out by name: {', '.join(RESERVED)}")
    print(f"- Budget: {Budget()}, cargo per oligo {Budget().ceiling} nt")
    large = product()
    print(f"- Largest cargo: product.dna, {len(large)} bp")
    sweep_fragments(large, ligase)
    sweep_caps(large, ligase)
    sweep_parts(ligase)
    sweep_oligo(ligase)
    sweep_pressure(ligase)
    sweep_reserved(ligase)
    sweep_scoring(large)
    oligo_count()


if __name__ == "__main__":
    main()
