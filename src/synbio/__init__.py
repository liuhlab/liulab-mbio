"""Pipelines for one named method, built on `mbio`.

This package imports `mbio`; nothing in `mbio` imports this one. `AGENTS.md`
states the boundary and the test for which side a module is on.

Only the version is re-exported, as `mbio` does. Import a name from the module that
owns it.
"""

from importlib.metadata import version

#: One distribution ships both import packages, so this is `mbio.__version__`. No
#: version string is written by hand anywhere in this repo.
__version__ = version("liulab-mbio")

__all__ = ["__version__"]
