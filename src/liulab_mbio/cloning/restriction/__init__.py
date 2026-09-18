"""Restriction and ligation cloning: two digests, a gel, and a ligase.

`plan_restriction` is the way in: it takes the vector, the plasmid the insert is cut out of and
the enzymes to cut both with, and runs the whole design, and `Plan.write` puts the product, the
oligo sheet and the bench protocol in one directory.
The submodules are the steps it is made of -- `digest`, `ligation`, `bench`, `verdicts` and
`steps`.
`digest` cuts a plasmid and says which enzyme left each end, `ligation` joins two pieces into a
circle and reads off what each junction spells, `bench` is the method's own documented numbers
and the reactions built from them, `verdicts` the checks a plan carries and the one nothing
sourced can judge, and `steps` the protocol's own steps and their order. What
any pipeline shares is `liulab_mbio.bench`, and whether two ends anneal at all is
`liulab_mbio.overhangs`.
"""

from liulab_mbio.cloning.restriction.plan import DEFAULT_HOST, Files, Plan, plan_restriction

__all__ = [
    "DEFAULT_HOST",
    "Files",
    "Plan",
    "plan_restriction",
]
