"""Gateway cloning: two att sites recombine, and the reaction rewrites the sites themselves.

`plan_gateway` is the way in: it takes a destination vector with either an entry clone or an
insert and its donor vector, designs the attB primers an insert carrying no att site needs,
simulates the recombinations, and `Plan.write` puts the clones, the oligo sheet and the bench
protocol in one directory.

The submodules are the steps it is made of. `att` owns the att site sequences, which one pairs
with which, and the arithmetic a junction follows; `design` owns the attB primer tail and the
PCR that puts it on an insert; `recombination` simulates one reaction on two records; `checks`
carries the verdicts spanning the experiment rather than one reaction; `oligos` says what each
ordered oligo is for; `bench` is each reaction's own components, incubation and stop; `steps` is
the protocol's steps and their order. What any pipeline shares is `liulab_mbio.bench`, and
`liulab_mbio.sites` is enzyme cut sites, which have nothing to do with an att site.
"""

from liulab_mbio.cloning.gateway.bench import DEFAULT_HOST
from liulab_mbio.cloning.gateway.design import FUSIONS, Amplicon, Fusion
from liulab_mbio.cloning.gateway.plan import ENTRY_FILE, Files, Plan, plan_gateway

__all__ = [
    "DEFAULT_HOST",
    "ENTRY_FILE",
    "FUSIONS",
    "Amplicon",
    "Files",
    "Fusion",
    "Plan",
    "plan_gateway",
]
