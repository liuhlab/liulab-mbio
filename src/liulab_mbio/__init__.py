"""Molecular biology design tools: sequences, enzymes, primers and Golden Gate cloning.

Only the version is re-exported. Import a name from the module that owns it --
``from liulab_mbio.io import read_record`` -- because names are spelled twice across modules
(`Check`, `Fragment`, `Junction`) and a flat re-export would import every dependency here.
"""

from importlib.metadata import version

#: Read from the installed distribution's metadata, which hatch-vcs fills from the newest
#: git tag. No version string is written by hand anywhere in this repo.
__version__ = version("liulab-mbio")

__all__ = ["__version__"]
