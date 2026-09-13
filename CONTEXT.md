# Context

## Glossary

### Sequence record

One DNA sequence with its topology, features, primers and notes: what a sequence file holds,
whatever its format.
_Avoid_: SeqRecord, file, map

### Feature

A named, typed annotation of a sequence record, such as a CDS or a promoter, lying on a strand
over one or more segments.
_Avoid_: annotation, region

### Segment

One contiguous span of a feature. Most features have one; a promoter drawn with its -35 and -10
boxes has several.
_Avoid_: range, interval, part

### Across the origin

Of a span on a circular sequence record: running through the last base and on into the first.
_Avoid_: wrapping, wraparound
