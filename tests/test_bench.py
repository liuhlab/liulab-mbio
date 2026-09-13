"""The shared bench package as a whole."""

import subprocess
import sys


def test_the_bench_imports_nothing_from_the_golden_gate_pipeline() -> None:
    # A fresh interpreter: this one has already imported every module under `src/`.
    loaded = subprocess.run(
        [sys.executable, "-c", "import sys, liulab_mbio.bench; print(*sys.modules)"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    assert [name for name in loaded if name.startswith("liulab_mbio.goldengate")] == []
