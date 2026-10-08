"""Golden Gate cloning: enzyme and overhang choice, assembly, and the reaction that joins it.

`plan_assembly` is the way in: it takes a vector and any number of inserts and runs the whole
design, and `Plan.write` puts the product, the primer sheet and the bench protocol in one
directory.
The submodules are the steps it is made of -- `design`, `assembly`, `oligos` and `steps`, where
`steps` is the protocol's own steps and their order. What any pipeline shares is
`liulab_mbio.bench`: NEB's enzymes, assembly reaction and cycling are
`liulab_mbio.bench.goldengate`, and the rules every overhang set is held to are
`liulab_mbio.overhangs`. `design` is not re-exported here: import it by module.
"""

from liulab_mbio.cloning.goldengate.plan import DEFAULT_HOST, Files, Plan, Site, plan_assembly
from liulab_mbio.cloning.plan import Orientation

__all__ = [
    "DEFAULT_HOST",
    "Files",
    "Orientation",
    "Plan",
    "Site",
    "plan_assembly",
]
