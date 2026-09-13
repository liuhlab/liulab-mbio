"""Golden Gate cloning: enzyme and overhang choice, assembly, and the reaction that joins it.

`plan_assembly` is the way in: it takes a vector and any number of inserts and runs the whole
design, and `Plan.write` puts the product, the primer sheet and the bench protocol in one
directory.
The submodules are the steps it is made of -- `design`, `assembly`, `bench` and `steps` -- plus
`ligase`, which reads a ligase fidelity matrix the user holds. `bench` is the assembly reaction
and its cycling; what any pipeline shares is `liulab_mbio.bench`. `design` and `ligase` are not
re-exported here: import them by module.
"""

from liulab_mbio.goldengate.plan import (
    DEFAULT_HOST,
    Files,
    Orientation,
    Phenotype,
    Plan,
    Site,
    plan_assembly,
    primer_sheet,
)

__all__ = [
    "DEFAULT_HOST",
    "Files",
    "Orientation",
    "Phenotype",
    "Plan",
    "Site",
    "plan_assembly",
    "primer_sheet",
]
