"""The layer rule: what a module may import is what sits below it."""

import subprocess
import sys

import pytest

#: Every module under this prefix is one cloning method's pipeline. Nothing below the pipelines
#: may load one, and nor may a pipeline that is not a cloning method.
CLONING = "liulab_mbio.goldengate"


@pytest.mark.parametrize("package", ["liulab_mbio.bench", "liulab_mbio.library"])
def test_no_cloning_pipeline_is_loaded(package: str) -> None:
    # A fresh interpreter: this one has already imported every module under `src/`.
    loaded = subprocess.run(
        [sys.executable, "-c", f"import sys, {package}; print(*sys.modules)"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()

    assert [name for name in loaded if name.startswith(CLONING)] == []
