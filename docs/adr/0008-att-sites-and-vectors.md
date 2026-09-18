---
search:
  exclude: true
---

# The att sequences ship, no vector does, and a vector is recognised by searching for its sites

Gateway needs the eight att sequences in three places: the tail a primer carries, the arithmetic
a reaction runs, and the records the tests are built from. It also has to recognise a donor or a
destination vector in a file someone hands it. `docs/research/gateway-cloning.md` settles the
licence question for both — the sequences may ship, on two independent footings, and a vendor
vector catalogue could legally ship too.

So: the eight 25 bp regions ship as constants `liulab_mbio.cloning.gateway.att` owns. No vector
catalogue ships. And a vector is recognised by **finding** its sites rather than by matching a
shipped constant — the 7 bp overlap exactly, the rest of the region within one base.

The note records two places where the vendor's own sequences drift between products. Matching a
constant against a drifted vector refuses a vector that works, and the refusal reads like a
wrong file. A search that tolerates drift outside the overlap accepts it, and still names what
it looked for when it finds nothing. The tolerance is measured rather than chosen: two sites
sharing an overlap differ by as few as two bases, so one mismatch is as far as a window may
drift and stay nearer one site than any other.

## Considered options

- **Ship a vendor vector catalogue**, so a user could name pDONR221 instead of supplying a file.
  It is legally possible and still wrong: the catalogue is a copy of sequences that drift, so it
  goes stale silently, and the plan would answer for maps nobody in this lab has read. The
  vectors on a bench are the only ones that matter, and the user has those files.
- **Match each site exactly.** Cheaper, and it fails on the one input the note says to expect.
  A drifted flank is a working vector, so an exact match turns correct work away.

## Consequences

Both vectors are the user's own files, always. A record carrying no site is refused with a
message naming every site looked for and the tolerance applied. The tolerance is derived from
the shipped regions, so it cannot drift away from them. Nothing in the package compares a
record's bases to an att constant, which is what keeps the two decisions one.
