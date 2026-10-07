---
search:
  exclude: true
---

# A coding sequence is what the bases spell, not what the feature type says

`liulab_mbio.sites.domesticate` treats a feature as a coding sequence when the record types it
`CDS`, **or** when the feature's own bases spell a whole protein: `ATG`, whole codons of
definite bases, one stop as the last codon and none before it. Nothing else changes. A site in
no such feature is still reported rather than edited.

The type field is the weakest part of an annotation. Addgene's feature list, a SnapGene export
and most files a user brings in type every element `misc_feature`, so a designer that trusts the
type refuses a free synonymous change on a real open reading frame. pLVX-TetOne-Puro-GFP is the
case: `docs/research/working-vector-plvx-tetone.md` translated PuroR and AmpR and found both
clean, yet BsmBI at 5636 inside PuroR came back as outside a coding sequence.

The feature is still what supplies the frame and the strand. Only its claim to code is checked,
against the bases the record already holds, so the rule adds no number and no search.

## Considered options

- **Type the two open reading frames `CDS` in our copy of the GenBank.** It fixes one fixture
  and leaves every other loosely typed vector refusing changes it could make. The divergence
  from the source file is then ours to keep correct for ever.
- **Find an open reading frame where no feature covers the site.** Six frames over a circular
  plasmid offer several frames across one six-base site, and separating a real frame from a
  chance one needs a minimum length nothing here measures. That is a guessed number, so it is
  not taken. A record that annotates nothing over a site keeps the old answer.
- **Accept any feature that is a whole number of codons.** Too loose: a promoter divisible by
  three would be recoded.

## Consequences

Every vector whose open reading frames are typed loosely is domesticated like one that types
them `CDS`. A partial coding sequence — no start, or the stop trimmed off, as one of the two
EGFP annotations on pLVX has it — is not read as coding unless the record types it, which is the
conservative direction: the site is reported, not edited in a frame nobody confirmed.
