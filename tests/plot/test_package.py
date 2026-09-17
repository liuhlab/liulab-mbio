"""The plot package as a whole."""

import subprocess
import sys


def test_importing_the_package_and_its_command_loads_no_converter() -> None:
    # A fresh interpreter: this one has already imported every module under `src/`.
    loaded = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, liulab_mbio.plot.cli, liulab_mbio.plot.convert; print(*sys.modules)",
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    assert not {"vl_convert", "pypdf"} & set(loaded)
