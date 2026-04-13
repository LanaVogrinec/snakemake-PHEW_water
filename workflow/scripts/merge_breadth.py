#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import re


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--breadth-files", nargs="+", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def get_sample_name(f):
    fname = os.path.basename(f)
    m = re.match(r"09_(.+)_breadth\.tsv", fname)
    if not m:
        raise ValueError(f"Unexpected filename format: {fname}")
    return m.group(1)


def main():
    args = parse_args()

    dfs = []

    for f in args.breadth_files:
        sample = get_sample_name(f)

        df = pd.read_csv(f, sep="\t", header=None, names=["contig", sample])
        df = df.set_index("contig")

        dfs.append(df)

    merged = pd.concat(dfs, axis=1).fillna(0)
    merged.to_csv(args.out, sep="\t")


if __name__ == "__main__":
    main()