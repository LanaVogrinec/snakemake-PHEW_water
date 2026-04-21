#!/usr/bin/env python3

import argparse
import pandas as pd
import os


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--breadth-files", nargs="+", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def get_sample_name(f):
    fname = os.path.basename(f)
    sample = fname.split("_breadth")[0].replace("09_", "")
    return sample


def main():
    args = parse_args()

    dfs = []

    for f in args.breadth_files:
        sample = get_sample_name(f)

        df = pd.read_csv(f, sep="\t", header=None, names=["rep_contig", sample])
        df = df.set_index("rep_contig")

        dfs.append(df)

    merged = pd.concat(dfs, axis=1).fillna(0)
    merged.to_csv(args.out, sep="\t")


if __name__ == "__main__":
    main()