"""PROTOTYPE, throwaway: an mbio record as SeqViz props, plain and with mbio's 28 enzymes."""

import json
import sys

from mbio.io import read_record

ENZYMES = json.load(open("src/mbio/data/enzymes.json"))["enzymes"]


def props(path):
    record = read_record(path)
    n = len(record.sequence)
    annotations, translations = [], []
    for f in record.features:
        if f.type == "source":
            continue
        start, end = f.segments[0].start, f.segments[-1].end
        item = {"name": f.name, "start": start, "end": end if end <= n else end - n,
                "direction": int(f.strand)}
        if f.color:
            item["color"] = f.color
        annotations.append(item)
        if f.type == "CDS":
            translations.append({"start": start, "end": end, "direction": int(f.strand) or 1})
    primers = [
        {"name": p.name, "start": b.start, "end": b.end, "direction": int(b.strand), "color": "#a389d4"}
        for p in record.primers for b in p.binding_sites
    ]
    base = {"name": record.name or path.split("/")[-2], "seq": record.sequence, "annotations": annotations,
            "primers": primers, "translations": translations}
    print(path, n, record.topology, len(annotations), "features", len(primers), "primer sites")
    return base


# SeqViz knows enzymes by its own lowercase keys and lacks some of mbio's, so each goes in as a
# custom enzyme: the site padded with N out to the furthest cut, cut offsets as mbio counts them.
custom = [
    {"name": e["name"], "rseq": e["site"] + "N" * max(0, max(e["top_cut"], e["bottom_cut"]) - len(e["site"])),
     "fcut": e["top_cut"], "rcut": e["bottom_cut"]}
    for e in ENZYMES if e["top_cut"] >= 0 and e["bottom_cut"] >= 0
]
print("custom enzymes", len(custom), "of", len(ENZYMES))
json.dump({"plain": props(sys.argv[1]), "crowded": {**props(sys.argv[2]), "enzymes": custom}},
          open(sys.argv[3], "w"))

