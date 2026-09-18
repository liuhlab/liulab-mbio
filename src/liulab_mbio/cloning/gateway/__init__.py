"""Gateway cloning: two att sites recombine, and the reaction rewrites the sites themselves.

`plan_gateway` is the way in: it takes a destination vector with either an entry clone or an
attB-flanked insert and its donor vector, simulates the recombinations, and `Plan.write` puts
the clones and the bench protocol in one directory.

The submodules are the steps it is made of. `att` owns the att site sequences, which one pairs
with which, and the arithmetic a junction follows; `recombination` simulates one reaction on
two records; `bench` is each reaction's own components, incubation and stop; `steps` is the
protocol's steps and their order. What any pipeline shares is `liulab_mbio.bench`, and
`liulab_mbio.sites` is enzyme cut sites, which have nothing to do with an att site.
"""

from liulab_mbio.cloning.gateway.bench import DEFAULT_HOST
from liulab_mbio.cloning.gateway.plan import ENTRY_FILE, Files, Plan, plan_gateway

__all__ = [
    "DEFAULT_HOST",
    "ENTRY_FILE",
    "Files",
    "Plan",
    "plan_gateway",
]
