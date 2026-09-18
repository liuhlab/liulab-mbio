"""Gateway cloning: two att sites recombine, and the reaction rewrites the sites themselves.

`plan_gateway` is the way in: it takes an entry clone and a destination vector, simulates the
LR reaction between them, and `Plan.write` puts the expression clone and the bench protocol in
one directory.

The submodules are the steps it is made of. `att` owns the att site sequences, which one pairs
with which, and the arithmetic a junction follows; `recombination` simulates one reaction on
two records; `bench` is each reaction's own components, incubation and stop; `steps` is the
protocol's steps and their order. What any pipeline shares is `liulab_mbio.bench`, and
`liulab_mbio.sites` is enzyme cut sites, which have nothing to do with an att site.
"""

from liulab_mbio.cloning.gateway.plan import DEFAULT_HOST, Files, Plan, plan_gateway

__all__ = [
    "DEFAULT_HOST",
    "Files",
    "Plan",
    "plan_gateway",
]
