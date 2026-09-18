"""Gibson assembly: the overlap at each junction, the product it makes, and its one reaction.

`plan_gibson` is the way in: it takes a vector and its insert and runs the whole design, and
`Plan.write` puts the plasmid, the oligo sheet and the bench protocol in one directory.
The submodules are the steps it is made of -- `design`, `assembly`, `bench`, `oligos` and
`steps`. `bench` is every assembly product's own documented numbers and the reaction and
incubation built from them, and `steps` the protocol's own steps and their order; what any
pipeline shares is `liulab_mbio.bench`. `design` is not re-exported here: import it by module.
"""

from liulab_mbio.cloning.gibson.plan import DEFAULT_HOST, Files, Plan, Site, plan_gibson

__all__ = [
    "DEFAULT_HOST",
    "Files",
    "Plan",
    "Site",
    "plan_gibson",
]
