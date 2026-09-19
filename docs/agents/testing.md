---
search:
  exclude: true
---

# Testing

Where a test goes, and what it may cost. The model is
[grugbrain.dev on testing](https://grugbrain.dev/#grug-on-testing): tests at the seams the
code already has, a few end-to-end tests per method, and none that repeats another.

## Keep rule

A test checks a promise, from a docstring, an ADR or `CONTEXT.md`, through a public function.
It uses real inputs and mocks none of the package's own parts. It goes when another test
checks the same fact, or when it checks incidental behaviour that a refactor keeping
behaviour would break. Never weaken an assertion to make a test faster.

## Seams

Test where the code already divides. Add no seam, and change no API, to move a check.

| Seam | Tested for |
| --- | --- |
| Each pipeline's way in: `plan_assembly`, `plan_gibson`, `plan_restriction`, `plan_gateway`, `plan_library` | the plan whole, and the files it writes |
| The CLI | its wiring |
| Below a plan: each method's `design`, `assembly`, `digest`, `ligation`, `recombination`, `checks`, `bench`, `oligos` and `steps`; `bench.validation`; `primers.design` and `primers.placement` | each detail they decide |

## End to end, per method

A method gets these end-to-end tests and no more:

- **The common path**: one shared plan, built once in a session- or module-scoped fixture
  that every plan test reads, as `tests/conftest.py` and `tests/cloning/gateway/conftest.py`
  do.
- **Critical edge cases**: one plan for each route the method takes differently, such as
  several inserts, a reversed insert, a linear vector, a site across the origin or blunt ends.
- **Same inputs, same bytes**: one test that plans the shared plan's inputs again and
  compares the files byte for byte.
- **The CLI**: one real run whose single call exercises every option's wiring. Refusal tests
  stay, since they fail before any planning.

## Detail checks

A detail one function decides, such as which overhang, which overlap or which band, is tested
on that function, not read off a whole plan. A failure then points at the function that broke,
and a refactor elsewhere leaves the test alone. The check stays at plan level only when no
function below exposes it. `tests/bench/test_validation.py` shows the shape: it calls the
colony PCR check directly, with no plan.

## Regression tests

Every fixed bug keeps its test, at the lowest level that still reproduces the bug.

## Time

pytest on CI has a budget of 30 s. No gate enforces it: runner speed varies, and a time gate
would fail a correct pull request. The suite runs in parallel workers, a file at a time each,
set in pytest's own options so the laptop and CI run alike. Over budget, trim by the rules
above, or make a plan's first design cheaper: every worker rebuilds the caches
`primers.placement` fills, so adding workers buys less than it looks.

A plan test's cost is its first primer design on a sequence the run has not seen:
`primers.placement` remembers what each sequence primes on a template, so a later design on
that template is cheap. Build an edge-case plan from sequences the suite already designs on,
such as the pUC19 and GFP records in `tests/conftest.py`, wherever the case allows.
