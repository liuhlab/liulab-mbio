---
name: sequence-verify
description: >-
  Check one clone's sequencing results against the plasmid it should be with `mbio`: Sanger
  `.ab1` reads or a whole-plasmid consensus are placed on the planned product, each junction and
  each insert gets a verdict, and every disagreement is placed on the record with the feature it
  falls in. Use whenever someone's sequencing has come back and they want to know whether a
  clone is right: Sanger reads of a miniprep, a whole-plasmid consensus from a vendor, which
  colony to keep, whether a mismatch sits in an insert or in the backbone, or a cloning protocol
  from this package reaching its sequencing step.
---

# Sequence verification

`mbio.verification` judges the reads. Your job is to hand it the record the clone should be and
every result the sample gave, then pass its verdicts on as they are. Never call a base, a
verdict or a position from your own reading of a trace or an alignment: the package works each
one out, and nothing checks one you invented (`docs/adr/0002-editable-protocols.md`).

## Run it

```bash
pixi run mbio sequence-verify PRODUCT RESULT... [--feature NAME]
```

`PRODUCT` is the record the clone should be, such as the `product.dna` a cloning plan wrote.
Each `RESULT` is read by its extension: `.ab1` is a Sanger trace, and FASTA, GenBank or `.dna`
is a whole-plasmid consensus, trusted whole. Give every result one sample gave in one call, such
as a read from each side of the insert: a base counts as read when any trusted result covers it.

What is judged comes from the product. A plan from this package tags each junction it writes,
and the inserts are the stretches between them. A record no plan wrote carries no tag, so name
its regions with `--feature`, once per feature. With neither, every disagreement is listed and
nothing is judged.

It prints one line per region, each result's own check, and whether the clone is verified, and
exits non-zero when it is not.

## From Python

```python
from mbio.verification import regions, verify
from mbio.verification.trace import read_trace
```

`verify(expected, results, regions)` returns a `Verification`: `checks`, one per region in
record order; `disagreements`, each a span on `expected` naming the results that show it and the
region or feature it falls in; `placements`, where each result landed and on which strand;
`status`; and `verified`. Spans are 0-based and half-open (`docs/adr/0001-coordinates.md`).
Read the docstrings rather than reconstructing a call:

```bash
pixi run python -c "from mbio.verification import verify; help(verify)"
```

## What a verdict means

| Verdict | On a region |
| --- | --- |
| pass | a trusted base of some result covers every base, and nothing disagrees |
| fail | every trusted result covering a base shares a disagreement there, or the base or the sample is mixed |
| warn | two trusted results contradict each other at a base, or more than 10 reads but fewer than 20 stand behind one |
| none | part of it went unread and nothing failed |

A clone is verified only when every region passes. Nothing failing is not enough: a junction no
read reached carries no verdict, and the clone is not verified until a read covers it.

A disagreement outside every region carries no verdict. It is named with the feature it falls
in, because whether it matters depends on what that feature does, and that call is the user's.
Pass it on rather than calling it harmless.

A consensus alone cannot show a mixed sample, since the commonest plasmid becomes the consensus.
The command says so for each consensus it reads, and so should you.

`molecular-cloning` plans the clone; come here once its reads are back.
