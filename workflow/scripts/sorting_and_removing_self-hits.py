#!/usr/bin/env python3

import argparse
import gzip
import sys

parser = argparse.ArgumentParser(
    description="Remove self-hits, canonicalize BLAST pairs, and sort for ANI/AF calculation"
)
parser.add_argument("--input", required=True, help="Input BLAST TSV")
parser.add_argument("--output", required=True, help="Output sorted TSV.gz")
args = parser.parse_args()

records = []

with open(args.input, "r") as fin:
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        q, s = parts[0], parts[1]

        # remove self-hits
        if q == s:
            continue

        # canonicalize pair
        a, b = sorted([q, s])

        # store: canonical_pair + original line
        records.append((a, b, parts))

# sort by canonical pair
records.sort(key=lambda x: (x[0], x[1]))

# write gzipped output
with gzip.open(args.output, "wt") as fout:
    for a, b, parts in records:
        fout.write("\t".join([a, b] + parts) + "\n")
