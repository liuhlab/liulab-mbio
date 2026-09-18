"""The layer rule: what a module may import is what sits below it."""

import subprocess
import sys

import pytest

#: Every module under this prefix is one cloning method's pipeline, save `SHARED`. Nothing below
#: the pipelines may load a method, and nor may a pipeline that is not a cloning method.
CLONING = "liulab_mbio.cloning"

#: What is not a method: the plan layer every plan is built on, and the package holding it.
SHARED = ("liulab_mbio.cloning", "liulab_mbio.cloning.plan")


@pytest.mark.parametrize("package", ["liulab_mbio.bench", "liulab_mbio.library"])
def test_no_cloning_pipeline_is_loaded(package: str) -> None:
    # A fresh interpreter: this one has already imported every module under `src/`.
    loaded = subprocess.run(
        [sys.executable, "-c", f"import sys, {package}; print(*sys.modules)"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()

    assert [name for name in loaded if name.startswith(CLONING) and name not in SHARED] == []
