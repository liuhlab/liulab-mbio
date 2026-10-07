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
    greedy_set_split,
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


def runner(search: Callable[..., Split], **extra: object) -> Callable[..., Split]:
    """Bind one search to the shared inputs, letting a caller override per run."""

    def run(sequence: str, ligase: LigaseProfile | None, **more: object) -> Split:
        return search(sequence, profile=ligase, reserved=RESERVED, **{**extra, **more})

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


SEARCHES = (
    ("proxy greedy", runner(simple_split)),
    ("set-fidelity greedy", runner(greedy_set_split)),
    ("branch and bound, cap 8", runner(clever_split, cap=8, seconds_limit=10.0)),
)


def sweep_three_parts(ligase: LigaseProfile | None, oligo: int, number: str) -> None:
    """All three searches over the 72 AP-1 parts, at one oligo length (M1a)."""
    budget = Budget(oligo=oligo)
    rows = []
    for name, run in SEARCHES:
        began = time.perf_counter()
        answers = [(sequence, run(sequence, ligase, budget=budget)) for _, sequence in parts()]
        seconds = time.perf_counter() - began
        values = [answer.value for _, answer in answers if answer.feasible]
        rows.append(
            [
                name,
                f"{sum(answer.feasible for _, answer in answers)}/{len(answers)}",
                sum(answer.fragments for _, answer in answers),
                f"{min(values):.4f}",
                f"{statistics.median(values):.4f}",
                f"{seconds * 1000 / len(answers):.1f}",
                f"{sum(assembles(sequence, answer) for sequence, answer in answers)}/{len(answers)}",
            ]
        )
    table(
        f"{number}. Three searches over all 72 AP-1 parts, 133 to 1,149 bp, {oligo} nt oligo",
        [
            "search",
            "parts solved",
            "oligos in total",
            "worst fidelity",
            "median fidelity",
            "ms per part",
            "parts that reassemble",
        ],
        rows,
    )


def sweep_three_counts(ligase: LigaseProfile | None, number: str) -> None:
    """All three searches across the 2,276 bp cargo's fragment-count range (M1b)."""
    sequence = product()
    span = len(sequence) - 4
    rows = []
    for count in range(2, 17):
        budget = Budget(oligo=54 + 24 + -(-span // count))
        row: list[object] = [count, budget.oligo]
        for _, run in SEARCHES:
            answer = run(sequence, ligase, budget=budget, count=count)
            row += [show(answer), f"{answer.seconds * 1000:.0f}", yes(assembles(sequence, answer))]
        row.append(lund(count))
        rows.append(row)
    table(
        f"{number}. Three searches against fragment count — 2,276 bp cargo, oligo widened per "
        "row so the fragment count is the only thing that moves",
        [
            "fragments",
            "oligo needed (nt)",
            "proxy fidelity",
            "proxy ms",
            "round trip",
            "set-greedy fidelity",
            "set-greedy ms",
            "round trip",
            "B&B fidelity",
            "B&B ms",
            "round trip",
            "Lund error-free clones (%)",
        ],
        rows,
    )


def yes(ok: bool) -> str:
    """Show whether a design rebuilt its own input."""
    return "yes" if ok else "no"


def outcome(answer: Split) -> str:
    """Say how a search ended: solved, out of time, or no legal set at all."""
    if not answer.feasible:
        return "refused"
    return "stopped early" if answer.note else "ok"


def sweep_budgets(ligase: LigaseProfile | None, number: str) -> None:
    """Measure what a bigger wall-clock budget buys the branch and bound (M2).

    The question is whether more time turns `stopped early` into `proved optimal`.

    One row per budget at the fragment count and cap the earlier run showed stopping early.
    Node count is reported for every cell, abandoned or not: its growth is the finding.
    """
    sequence = product()
    rows = []
    for count, cap in ((12, 8),):
        for limit in (10.0, 60.0, 300.0):
            answer = clever_split(
                sequence,
                profile=ligase,
                count=count,
                reserved=RESERVED,
                cap=cap,
                seconds_limit=limit,
            )
            rows.append(
                [
                    count,
                    cap,
                    f"{limit:.0f}",
                    show(answer),
                    "stopped early" if answer.note else "PROVED optimal over kept candidates",
                    answer.nodes,
                    answer.scored,
                    f"{answer.seconds:.1f}",
                ]
            )
            print(f"<progress> {count} frag cap {cap} limit {limit}: {rows[-1]}", flush=True)
    table(
        f"{number}. What a generous time budget buys — 2,276 bp cargo, 300 nt oligo",
        [
            "fragments",
            "candidates per cut",
            "budget (s)",
            "fidelity",
            "outcome",
            "nodes",
            "fidelity calls",
            "wall (s)",
        ],
        rows,
    )


def grown(target: int) -> str:
    """Realistic cargo of about `target` bases: AP-1 parts concatenated in frame, cycling."""
    blocks = [sequence[: len(sequence) // 3 * 3] for _, sequence in parts()]
    built: list[str] = []
    total = 0
    at = 0
    while total < target:
        block = blocks[at % len(blocks)]
        built.append(block)
        total += len(block)
        at += 1
    return "".join(built)


def sweep_scale(ligase: LigaseProfile | None, number: str) -> None:
    """How every stage scales with cargo length, at two oligo lengths (M3)."""
    import resource

    rows = []
    for target in (2300, 5000, 10000, 20000):
        cargo = grown(target)
        for oligo in (300, 1000):
            budget = Budget(oligo=oligo)
            began = time.perf_counter()
            try:
                least = fewest_fragments(len(cargo), budget)
            except ValueError:
                rows.append([len(cargo), oligo, "none", "-", "-", "-", "-", "-", "-", "-", "-"])
                continue
            dp = time.perf_counter() - began
            row: list[object] = [len(cargo), oligo, least, f"{dp * 1000:.0f}"]
            for name, run in SEARCHES:
                answer = run(
                    cargo,
                    ligase,
                    budget=budget,
                    **({"seconds_limit": 60.0} if "bound" in name else {}),
                )
                row += [
                    show(answer),
                    f"{answer.seconds * 1000:.0f}",
                    outcome(answer),
                ]
            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
            row += [f"{peak:.0f}", lund(least) or "off their scale"]
            rows.append(row)
            print(f"<progress> {len(cargo)} bp {oligo} nt: {row}", flush=True)
    table(
        f"{number}. Scaling — AP-1 parts concatenated in frame, B&B at cap 8 with a 60 s budget",
        [
            "cargo (bp)",
            "oligo (nt)",
            "fewest fragments",
            "DP ms",
            "proxy fidelity",
            "proxy ms",
            "outcome",
            "set-greedy fidelity",
            "set-greedy ms",
            "outcome",
            "B&B fidelity",
            "B&B ms",
            "outcome",
            "peak RSS (MB)",
            "Lund error-free clones (%)",
        ],
        rows,
    )


def sweep_four_parts(ligase: LigaseProfile | None, oligo: int, number: str) -> None:
    """All four searches over the 72 AP-1 parts, the fourth uncapped (M4).

    The uncapped branch and bound is the one that can prove optimality rather than prove it
    over kept candidates, so the proved and stopped columns are the point of the table.
    """
    budget = Budget(oligo=oligo)
    uncapped = ("branch and bound, uncapped", runner(clever_split, cap=None, seconds_limit=15.0))
    rows = []
    for name, run in (*SEARCHES, uncapped):
        began = time.perf_counter()
        answers = []
        for part, sequence in parts():
            answer = run(sequence, ligase, budget=budget)
            answers.append((part, sequence, answer))
            if "uncapped" in name:
                print(
                    f"<progress> {oligo} {part} {len(sequence)}bp {answer.fragments}f "
                    f"{answer.value:.4f} {answer.seconds:.1f}s {outcome(answer)}",
                    flush=True,
                )
        seconds = time.perf_counter() - began
        values = [answer.value for _, _, answer in answers if answer.feasible]
        stopped = sum(1 for _, _, answer in answers if answer.feasible and answer.note)
        searching = "bound" in name
        slow = max(answers, key=lambda one: one[2].seconds)
        rows.append(
            [
                name,
                f"{sum(answer.feasible for _, _, answer in answers)}/{len(answers)}",
                f"{len(answers) - stopped}" if searching else "n/a",
                f"{stopped}" if searching else "n/a",
                sum(answer.fragments for _, _, answer in answers),
                f"{min(values):.4f}",
                f"{statistics.median(values):.4f}",
                f"{seconds * 1000 / len(answers):.1f}",
                f"{seconds:.1f}",
                f"{slow[0]}, {len(slow[1])} bp, {slow[2].fragments} frags, {slow[2].seconds:.1f} s",
            ]
        )
    table(
        f"{number}. Four searches over all 72 AP-1 parts, {oligo} nt oligo, the uncapped "
        "branch and bound held to 15 s per part",
        [
            "search",
            "parts solved",
            "proved optimal",
            "stopped early",
            "oligos in total",
            "worst fidelity",
            "median fidelity",
            "ms per part",
            "sweep (s)",
            "slowest part",
        ],
        rows,
    )


def sweep_window(ligase: LigaseProfile | None, number: str) -> None:
    """Whether a wider window removes the set-fidelity greedy's refusals (M5)."""
    sequence = product()
    span = len(sequence) - 4
    rows = []
    for count in (12, 13, 14):
        budget = Budget(oligo=54 + 24 + -(-span // count))
        row: list[object] = [count]
        for width in (16, 32, 64):
            answer = greedy_set_split(
                sequence,
                budget=budget,
                profile=ligase,
                reserved=RESERVED,
                count=count,
                window=width,
            )
            row += [show(answer), f"{answer.seconds * 1000:.0f}"]
        rows.append(row)
    table(
        f"{number}. The set-fidelity greedy's window widened — 2,276 bp cargo, oligo widened "
        "per row as in M1c",
        [
            "fragments",
            "window 16",
            "ms",
            "window 32",
            "ms",
            "window 64",
            "ms",
        ],
        rows,
    )


def heading(ligase: LigaseProfile | None) -> None:
    """Say what scored the tables below."""
    held = "FileS03_T4_18h_25C.xlsx (user-held T4 18 h 25 C)"
    print(f"- Scored by: {held if ligase else 'shipped BsaI matrix, the T4 file is ABSENT'}")
    print(f"- Reserved overhangs held out by name: {', '.join(RESERVED)}")
    print(f"- Overhead per oligo: {Budget().overhead} nt; shortest fragment {Budget().floor} nt")


def main() -> None:
    """Print the tables the argument names, or every table."""
    wanted = sys.argv[1:] or ["old", "m1", "m2", "m3", "m4", "m5"]
    ligase = profile()
    print("# Prototype measurements for issue #261\n")
    heading(ligase)
    large = product()
    print(f"- Largest shipped cargo: product.dna, {len(large)} bp")
    if "old" in wanted:
        sweep_fragments(large, ligase)
        sweep_caps(large, ligase)
        sweep_parts(ligase)
        sweep_oligo(ligase)
        sweep_pressure(ligase)
        sweep_reserved(ligase)
        sweep_scoring(large)
        oligo_count()
    if "m1" in wanted:
        sweep_three_parts(ligase, 300, "M1a")
        sweep_three_parts(ligase, 350, "M1b")
        sweep_three_counts(ligase, "M1c")
    if "m2" in wanted:
        sweep_budgets(ligase, "M2")
    if "m3" in wanted:
        sweep_scale(ligase, "M3")
    if "m4" in wanted or "m4a" in wanted:
        sweep_four_parts(ligase, 300, "M4a")
    if "m4" in wanted or "m4b" in wanted:
        sweep_four_parts(ligase, 350, "M4b")
    if "m5" in wanted:
        sweep_window(ligase, "M5")


if __name__ == "__main__":
    main()
